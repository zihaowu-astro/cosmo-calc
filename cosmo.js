"use strict";

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
const C_KMS = 299792.458;          // speed of light [km/s]
const MPC_PER_GLY = 306.601;       // 1 Gly = 306.601 Mpc
const SEC_PER_GYR = 3.15576e16;    // seconds in a Gyr (Julian)
const KM_PER_MPC = 3.0856775815e19;

// ---------------------------------------------------------------------------
// Fixed cosmologies (flat LambdaCDM)
//   Planck 2018: TT,TE,EE+lowE+lensing+BAO (Planck 2018 VI, Table 2)
//   WMAP9:       Hinshaw et al. 2013, nine-year WMAP+eCMB+BAO+H0
// ---------------------------------------------------------------------------
const COSMOLOGIES = {
  planck18: { name: "Planck 2018", H0: 67.66, OmegaM: 0.3111, OmegaL: 0.6889 },
  wmap9:    { name: "WMAP9",       H0: 69.32, OmegaM: 0.2865, OmegaL: 0.7135 },
};

// Build the full parameter set for a given cosmology, including radiation.
function buildModel(key) {
  const base = COSMOLOGIES[key];
  const h = base.H0 / 100;
  // Radiation density (photons + 3 massless neutrino species), per Wright.
  const OmegaR = 4.165e-5 / (h * h);
  const OmegaK = 1 - base.OmegaM - base.OmegaL - OmegaR;
  return {
    ...base,
    h,
    OmegaR,
    OmegaK,
    DH: C_KMS / base.H0,             // Hubble distance [Mpc]
    tH: KM_PER_MPC / base.H0 / SEC_PER_GYR, // Hubble time [Gyr]
  };
}

// Dimensionless Hubble parameter E(z) = H(z)/H0.
function Ez(z, m) {
  const zp1 = 1 + z;
  return Math.sqrt(
    m.OmegaR * zp1 ** 4 +
    m.OmegaM * zp1 ** 3 +
    m.OmegaK * zp1 ** 2 +
    m.OmegaL
  );
}

// ---------------------------------------------------------------------------
// Numerical integration: adaptive Simpson's rule
// ---------------------------------------------------------------------------
function adaptiveSimpson(f, a, b, tol = 1e-9, maxDepth = 50) {
  function simpson(a, b, fa, fb, fm) {
    return ((b - a) / 6) * (fa + 4 * fm + fb);
  }
  function recurse(a, b, fa, fb, fm, whole, tol, depth) {
    const m = (a + b) / 2;
    const lm = (a + m) / 2;
    const rm = (m + b) / 2;
    const flm = f(lm);
    const frm = f(rm);
    const left = simpson(a, m, fa, fm, flm);
    const right = simpson(m, b, fm, fb, frm);
    if (depth <= 0 || Math.abs(left + right - whole) <= 15 * tol) {
      return left + right + (left + right - whole) / 15;
    }
    return (
      recurse(a, m, fa, fm, flm, left, tol / 2, depth - 1) +
      recurse(m, b, fm, fb, frm, right, tol / 2, depth - 1)
    );
  }
  if (b <= a) return 0;
  const fa = f(a);
  const fb = f(b);
  const m = (a + b) / 2;
  const fm = f(m);
  const whole = simpson(a, b, fa, fb, fm);
  return recurse(a, b, fa, fb, fm, whole, tol, maxDepth);
}

// ---------------------------------------------------------------------------
// Cosmological quantities at redshift z
// ---------------------------------------------------------------------------
function compute(z, m) {
  // Comoving radial distance: D_C = D_H * integral_0^z dz'/E(z')
  const DC = m.DH * adaptiveSimpson((zp) => 1 / Ez(zp, m), 0, z);

  // Transverse comoving distance D_M accounts for curvature.
  let DM;
  const sqrtOk = Math.sqrt(Math.abs(m.OmegaK));
  if (m.OmegaK > 1e-9) {
    DM = (m.DH / sqrtOk) * Math.sinh((sqrtOk * DC) / m.DH);
  } else if (m.OmegaK < -1e-9) {
    DM = (m.DH / sqrtOk) * Math.sin((sqrtOk * DC) / m.DH);
  } else {
    DM = DC; // flat
  }

  const DA = DM / (1 + z);          // angular diameter distance
  const DL = DM * (1 + z);          // luminosity distance

  // Comoving volume enclosed (general curvature form).
  let VC;
  const x = DM / m.DH;
  if (m.OmegaK > 1e-9) {
    VC = ((4 * Math.PI * m.DH ** 3) / (2 * m.OmegaK)) *
      (x * Math.sqrt(1 + m.OmegaK * x * x) -
        (1 / sqrtOk) * Math.asinh(sqrtOk * x));
  } else if (m.OmegaK < -1e-9) {
    VC = ((4 * Math.PI * m.DH ** 3) / (2 * m.OmegaK)) *
      (x * Math.sqrt(1 + m.OmegaK * x * x) -
        (1 / sqrtOk) * Math.asin(sqrtOk * x));
  } else {
    VC = (4 / 3) * Math.PI * DM ** 3;
  }

  // Ages via integral over scale factor a = 1/(1+z):
  //   t(a) = t_H * integral_0^a da' / (a' * E(a'))
  const Ea = (a) => {
    return Math.sqrt(
      m.OmegaR / a ** 2 +
      m.OmegaM / a +
      m.OmegaK +
      m.OmegaL * a ** 2
    );
    // note: a'*E(a') = sqrt(OmegaR/a^2 + OmegaM/a + OmegaK + OmegaL*a^2)
  };
  const ageIntegrand = (a) => 1 / Ea(a);
  const aTarget = 1 / (1 + z);
  const ageNow = m.tH * adaptiveSimpson(ageIntegrand, 0, 1);
  const ageAtZ = m.tH * adaptiveSimpson(ageIntegrand, 0, aTarget);
  const lookback = ageNow - ageAtZ;

  // Angular size scale: physical kpc per arcsecond.
  const kpcPerArcsec = (DA * 1000 * Math.PI) / (180 * 3600);

  // Radial scale: proper (physical) line-of-sight distance per dz = 0.01.
  //   dl_proper/dz = D_H / (E(z) (1+z))   [Mpc] -> pkpc, times dz
  const pkpcPerDz01 = (m.DH / (Ez(z, m) * (1 + z))) * 1000 * 0.01;

  return {
    z,
    Ez: Ez(z, m),
    ageNow,
    ageAtZ,
    lookback,
    comovingRadialMpc: DC,
    comovingRadialGly: DC / MPC_PER_GLY,
    DA_Mpc: DA,
    DA_Gly: DA / MPC_PER_GLY,
    DL_Mpc: DL,
    DL_Gly: DL / MPC_PER_GLY,
    kpcPerArcsec,
    pkpcPerDz01,
    comovingVolumeGpc3: VC / 1e9, // Mpc^3 -> Gpc^3
    // Distance modulus for per-frequency flux density F_nu, including the
    // bandwidth-compression K-correction term for a flat-F_nu spectrum:
    //   mu = 5 log10(D_L/10pc) - 2.5 log10(1+z)
    distmod: 5 * Math.log10((DL * 1e6) / 10) - 2.5 * Math.log10(1 + z),
  };
}

// ---------------------------------------------------------------------------
// Formatting helpers
// ---------------------------------------------------------------------------
function fmt(x, digits = 4) {
  if (!isFinite(x)) return "—";
  if (x === 0) return "0";
  const abs = Math.abs(x);
  if (abs >= 1e5 || abs < 1e-3) return x.toExponential(3);
  return x.toLocaleString(undefined, {
    maximumSignificantDigits: digits + 1,
    minimumSignificantDigits: Math.min(digits, 4),
  });
}

// ---------------------------------------------------------------------------
// UI
// ---------------------------------------------------------------------------
const form = document.getElementById("calc-form");
const cosmoSel = document.getElementById("cosmology");
const zInput = document.getElementById("redshift");
const paramsEl = document.getElementById("params");
const resultsEl = document.getElementById("results");
const primaryBody = document.querySelector("#primary-table tbody");
const secondaryBody = document.querySelector("#secondary-table tbody");
const plotsEl = document.getElementById("plots");

// ---------------------------------------------------------------------------
// SVG line-plot helpers (dependency-free)
// ---------------------------------------------------------------------------
function niceTicks(min, max, n) {
  const span = max - min || 1;
  const step0 = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const ticks = [];
  for (let t = Math.ceil(min / step) * step; t <= max + step * 1e-6; t += step) {
    ticks.push(t);
  }
  return ticks;
}

function tickLabel(v) {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-2) return v.toExponential(0);
  return Number(v.toPrecision(3)).toLocaleString();
}

function linePlot(subtitle, pts, curZ, curVal) {
  const W = 260, H = 196, padL = 42, padR = 12, padT = 26, padB = 38;
  const xmin = 0, xmax = 20;
  const ys = pts.map((p) => p[1]).filter(isFinite);
  let ymin = Math.min(...ys), ymax = Math.max(...ys);
  if (ymin === ymax) { ymin -= 1; ymax += 1; }
  const padY = (ymax - ymin) * 0.08;
  ymin -= padY; ymax += padY;
  const sx = (z) => padL + ((z - xmin) / (xmax - xmin)) * (W - padL - padR);
  const sy = (v) => H - padB - ((v - ymin) / (ymax - ymin)) * (H - padT - padB);

  let d = "", started = false;
  for (const [z, v] of pts) {
    if (!isFinite(v)) { started = false; continue; }
    d += (started ? "L" : "M") + sx(z).toFixed(1) + " " + sy(v).toFixed(1) + " ";
    started = true;
  }

  let g = "";
  for (const t of niceTicks(ymin, ymax, 4)) {
    if (t < ymin || t > ymax) continue;
    const y = sy(t).toFixed(1);
    g += `<line class="grid" x1="${padL}" y1="${y}" x2="${W - padR}" y2="${y}"/>`;
    g += `<text class="tick" x="${padL - 4}" y="${(+y + 3).toFixed(1)}" text-anchor="end">${tickLabel(t)}</text>`;
  }
  for (const t of [0, 5, 10, 15, 20]) {
    const x = sx(t).toFixed(1);
    g += `<text class="tick" x="${x}" y="${H - padB + 12}" text-anchor="middle">${t}</text>`;
  }

  let marker = "";
  if (curZ >= xmin && curZ <= xmax) {
    const x = sx(curZ).toFixed(1);
    const labelLeft = curZ > (xmin + xmax) * 0.7;
    marker += `<line x1="${x}" y1="${padT}" x2="${x}" y2="${H - padB}" stroke="var(--marker)" stroke-width="1.1" stroke-dasharray="4 3"/>`;
    marker += `<text class="ztag" x="${labelLeft ? +x - 3 : +x + 3}" y="${padT + 9}" text-anchor="${labelLeft ? "end" : "start"}">z = ${(+curZ).toPrecision(3) / 1}</text>`;
    if (isFinite(curVal) && curVal >= ymin && curVal <= ymax) {
      marker += `<circle cx="${x}" cy="${sy(curVal).toFixed(1)}" r="3.2" fill="var(--marker)"/>`;
    }
  }

  const xMid = (padL + W - padR) / 2;
  const html = `<figure class="panel"><svg viewBox="0 0 ${W} ${H}">
    <text class="psubtitle" x="${xMid}" y="14" text-anchor="middle">${subtitle}</text>
    ${g}
    <line class="axis" x1="${padL}" y1="${padT}" x2="${padL}" y2="${H - padB}"/>
    <line class="axis" x1="${padL}" y1="${H - padB}" x2="${W - padR}" y2="${H - padB}"/>
    <path d="${d}" fill="none" stroke="var(--accent)" stroke-width="1.6"/>
    ${marker}
    <g class="hover" style="display:none">
      <line class="hline" y1="${padT}" y2="${H - padB}"/>
      <circle class="hdot" r="3"/>
      <text class="htext htext1" y="${padT + 12}"></text>
      <text class="htext htext2" y="${padT + 24}"></text>
    </g>
    <rect class="capture" x="${padL}" y="${padT}" width="${W - padL - padR}" height="${H - padT - padB}" fill="transparent"/>
    <text class="axtitle" x="${xMid}" y="${H - 5}" text-anchor="middle" font-style="italic">z</text>
  </svg></figure>`;
  return { html, meta: { W, H, padL, padR, padT, padB, xmin, xmax, ymin, ymax, pts } };
}

function attachHover(svg, meta, unit, dec) {
  const { W, H, padL, padR, padT, padB, xmin, xmax, ymin, ymax, pts } = meta;
  const sx = (z) => padL + ((z - xmin) / (xmax - xmin)) * (W - padL - padR);
  const sy = (v) => H - padB - ((v - ymin) / (ymax - ymin)) * (H - padT - padB);
  const hg = svg.querySelector(".hover");
  const hline = svg.querySelector(".hline");
  const hdot = svg.querySelector(".hdot");
  const ht1 = svg.querySelector(".htext1");
  const ht2 = svg.querySelector(".htext2");
  const z0 = pts[0][0];
  const zN = pts[pts.length - 1][0];
  const step = pts.length > 1 ? pts[1][0] - pts[0][0] : 0.05;

  function zAt(evt) {
    const rect = svg.getBoundingClientRect();
    const vbx = ((evt.clientX - rect.left) / rect.width) * W;
    const z = xmin + ((vbx - padL) / (W - padL - padR)) * (xmax - xmin);
    return Math.max(z0, Math.min(zN, z));
  }

  function move(evt) {
    const z = zAt(evt);
    const idx = Math.max(0, Math.min(pts.length - 1, Math.round((z - z0) / step)));
    const zv = pts[idx][0];
    const val = pts[idx][1];
    if (!isFinite(val)) { hg.style.display = "none"; return; }
    const x = sx(zv);
    hline.setAttribute("x1", x);
    hline.setAttribute("x2", x);
    hdot.setAttribute("cx", x);
    hdot.setAttribute("cy", sy(val));
    const left = zv > (xmin + xmax) * 0.5;
    const tx = left ? x - 5 : x + 5;
    const anchor = left ? "end" : "start";
    ht1.setAttribute("x", tx);
    ht1.setAttribute("text-anchor", anchor);
    ht1.textContent = `z = ${zv.toFixed(1)}`;
    ht2.setAttribute("x", tx);
    ht2.setAttribute("text-anchor", anchor);
    ht2.textContent = `${val.toLocaleString(undefined, { minimumFractionDigits: dec, maximumFractionDigits: dec })} ${unit}`;
    hg.style.display = "";
  }

  svg.addEventListener("mousemove", move);
  svg.addEventListener("mouseleave", () => { hg.style.display = "none"; });
  svg.addEventListener("click", (evt) => {
    zInput.value = Number(zAt(evt).toFixed(2));
    render();
  });
}

function renderPlots(m, z, r) {
  const curve = [];
  for (let zz = 0.05; zz <= 20.0001; zz += 0.05) {
    const c = compute(zz, m);
    curve.push({ z: zz, scale: c.kpcPerArcsec, dm: c.distmod, age: c.ageAtZ });
  }
  const it = (s) => `<tspan font-style="italic">${s}</tspan>`;
  const panels = [
    { sub: `Angular scale (kpc/″)`, key: "scale", unit: "kpc/″", cur: r.kpcPerArcsec, dec: 1 },
    { sub: `Distance modulus, ${it("μ")} (mag)`, key: "dm", unit: "mag", cur: r.distmod, dec: 1 },
    { sub: `Age, ${it("t")}(${it("z")}) (Gyr)`, key: "age", unit: "Gyr", cur: r.ageAtZ, dec: 2 },
  ];

  const metas = [];
  plotsEl.innerHTML = panels.map((p) => {
    const { html, meta } = linePlot(p.sub, curve.map((c) => [c.z, c[p.key]]), z, p.cur);
    metas.push(meta);
    return html;
  }).join("");

  plotsEl.querySelectorAll(".panel svg").forEach((svg, i) => {
    attachHover(svg, metas[i], panels[i].unit, panels[i].dec);
  });
}

function row(label, value, unit) {
  const tr = document.createElement("tr");
  const u = unit ? ` <span class="unit">${unit}</span>` : "";
  tr.innerHTML = `<td class="label">${label}</td><td class="value">${value}${u}</td>`;
  return tr;
}

function showParams(m) {
  const sep = `<span class="sep">|</span>`;
  const OmegaR = m.OmegaR.toExponential(2).replace("e-", "&times;10<sup>&minus;") + "</sup>";
  paramsEl.innerHTML =
    `<span class="model-params">` +
    `<i>H</i><sub>0</sub> = ${m.H0} km s<sup>&minus;1</sup> Mpc<sup>&minus;1</sup>${sep}` +
    `&Omega;<sub>m</sub> = ${m.OmegaM}${sep}` +
    `&Omega;<sub>&Lambda;</sub> = ${m.OmegaL}${sep}` +
    `&Omega;<sub>r</sub> = ${OmegaR}${sep}` +
    `flat (&Omega;<sub>k</sub> = 0)` +
    `</span>`;
}

function render() {
  const m = buildModel(cosmoSel.value);
  const z = parseFloat(zInput.value);
  showParams(m);

  if (!isFinite(z) || z < 0) {
    resultsEl.hidden = true;
    return;
  }

  const r = compute(z, m);

  const sig4 = (x) => Number(x.toPrecision(4)).toLocaleString();
  const sig3 = (x) => Number(x.toPrecision(3)).toLocaleString();
  const int = (x) => Math.round(x).toLocaleString();
  primaryBody.innerHTML = "";
  primaryBody.append(
    row("Angular size scale", sig3(r.kpcPerArcsec), "kpc/&Prime;"),
    row(`Distance modulus, <i>&mu;</i> <a class="info" href="distance-modulus.html" target="_blank" rel="noopener" title="Includes the &minus;2.5 log(1+z) band-shift term. Why? Click to read.">&#9432;</a>`, sig4(r.distmod), "mag"),
    row("Age at redshift <i>z</i>", sig4(r.ageAtZ), "Gyr"),
  );

  secondaryBody.innerHTML = "";
  secondaryBody.append(
    row("Age of the Universe (<i>z</i> = 0)", sig4(r.ageNow), "Gyr"),
    row("Lookback time", sig4(r.lookback), "Gyr"),
    row("Comoving radial distance, <i>D</i><sub>C</sub>", int(r.comovingRadialMpc), "Mpc"),
    row("Angular-diameter distance, <i>D</i><sub>A</sub>", int(r.DA_Mpc), "Mpc"),
    row("Luminosity distance, <i>D</i><sub>L</sub>", int(r.DL_Mpc), "Mpc"),
    row("Comoving volume (&lt; <i>z</i>)", sig4(r.comovingVolumeGpc3), "Gpc<sup>3</sup>"),
    row("Hubble parameter, <i>H</i>(<i>z</i>)", sig4(m.H0 * r.Ez), "km s<sup>&minus;1</sup> Mpc<sup>&minus;1</sup>"),
    row("Expansion rate, <i>E</i>(<i>z</i>) = <i>H</i>(<i>z</i>)/<i>H</i><sub>0</sub>", r.Ez.toFixed(2)),
    row("Radial scale (&Delta;<i>z</i> = 0.01)", int(r.pkpcPerDz01), "pkpc"),
  );

  renderPlots(m, z, r);
  resultsEl.hidden = false;
}

form.addEventListener("submit", (e) => {
  e.preventDefault();
  render();
});
cosmoSel.addEventListener("change", render);

// Initial render on load.
render();
