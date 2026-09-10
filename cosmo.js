"use strict";

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
const C_KMS = 299792.458;          // speed of light [km/s]
const MPC_PER_GLY = 306.601;       // 1 Gly = 306.601 Mpc
const SEC_PER_GYR = 3.15576e16;    // seconds in a Gyr (Julian)
const KM_PER_MPC = 3.0856775815e19;
const AB_ZP_NJY = 3.631e12;        // AB zero point, 3631 Jy, expressed in nJy
const G_KPC = 4.30091e-6;          // G [kpc (km/s)^2 / Msun]

// ---------------------------------------------------------------------------
// Fixed cosmologies (flat LambdaCDM). OmegaM is the published value, which counts
// the massive neutrino as matter; neutrinos are then modelled as massless below.
//   Planck 2018: TT,TE,EE+lowE+lensing+BAO (Planck 2018 VI, Table 2)
//   Planck 2015: TT,TE,EE+lowP+lensing+ext (Planck 2015 XIII, Table 4)
//   Planck 2013: Planck+WP+highL+BAO       (Planck 2013 XVI, Table 5)
//   WMAP9:       Hinshaw et al. 2013, nine-year WMAP+eCMB+BAO+H0 (Table 4)
//   WMAP7:       Komatsu et al. 2011, seven-year WMAP+BAO+H0 ML  (Table 1)
//   h70:         the values conventionally adopted in the literature, no reference
// ---------------------------------------------------------------------------
const COSMOLOGIES = {
  planck18: { name: "Planck 2018",   H0: 67.66, OmegaM: 0.3111, OmegaL: 0.6889 },
  planck15: { name: "Planck 2015",   H0: 67.74, OmegaM: 0.3089, OmegaL: 0.6911 },
  planck13: { name: "Planck 2013",   H0: 67.77, OmegaM: 0.3086, OmegaL: 0.6914 },
  wmap9:    { name: "WMAP9",         H0: 69.32, OmegaM: 0.2865, OmegaL: 0.7135 },
  wmap7:    { name: "WMAP7",         H0: 70.40, OmegaM: 0.2720, OmegaL: 0.7280 },
  h70:      { name: "H₀ = 70, Ωₘ = 0.3", H0: 70.00, OmegaM: 0.3000, OmegaL: 0.7000 },
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

  // Same scale in comoving units: comoving Mpc per arcminute.
  const cMpcPerArcmin = (DM * Math.PI) / (180 * 60);

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
    cMpcPerArcmin,
    pkpcPerDz01,
    comovingVolumeGpc3: VC / 1e9, // Mpc^3 -> Gpc^3
    // Distance modulus for per-frequency flux density F_nu, including the
    // bandwidth-compression K-correction term for a flat-F_nu spectrum:
    //   mu = 5 log10(D_L/10pc) - 2.5 log10(1+z)
    distmod: 5 * Math.log10((DL * 1e6) / 10) - 2.5 * Math.log10(1 + z),
  };
}

// ---------------------------------------------------------------------------
// Halo scaling relations at redshift z
// ---------------------------------------------------------------------------
// M200 is the mass inside R200, the radius where the mean enclosed density is
// 200 rho_crit(z), with rho_crit = 3 H(z)^2 / (8 pi G):
//   M200 = (4 pi / 3) 200 rho_crit R200^3 = 100 H(z)^2 R200^3 / G
// The circular velocity there is V200 = sqrt(G M200 / R200) = 10 H(z) R200,
// so equivalently M200 = V200^3 / (10 G H(z)).
//   M200 [Msun], Hz [km/s/Mpc]  ->  R200 [kpc], V200 [km/s]
function halo(M200, Hz) {
  const H = Hz / 1000;                                    // km/s/kpc
  const R200 = Math.cbrt((G_KPC * M200) / (100 * H * H));
  return { R200, V200: 10 * H * R200 };
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

// Flux densities span decades; keep ~3 significant digits without exponents.
function fmtFlux(x) {
  if (!isFinite(x)) return "—";
  const d = x >= 100 ? 0 : x >= 10 ? 1 : x >= 1 ? 2 : 3;
  return x.toLocaleString(undefined, { minimumFractionDigits: d, maximumFractionDigits: d });
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
const haloPlotsEl = document.getElementById("halo-plots");

// ---------------------------------------------------------------------------
// SVG line-plot helpers (dependency-free)
// ---------------------------------------------------------------------------
// One panel geometry, shared by every plot (viewBox units).
// The panel title lives in HTML above the SVG, so that formulas can be set in
// MathML; padT only has to clear the topmost gridline label and the z tag.
const PLOT = { W: 380, H: 226, padL: 57, padR: 16, padT: 14, padB: 44 };
const MUV_MIN = -24, MUV_MAX = -16, MUV_STEP = 0.1;   // M_UV grid of the flux panel

function niceTicks(min, max, n) {
  if (!isFinite(min) || !isFinite(max)) return [];
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
  const s = a >= 1e4 || a < 1e-2 ? v.toExponential(0) : Number(v.toPrecision(3)).toLocaleString();
  return s.replace("-", "−");
}

// Ticks for a log10 axis, given the axis range in log10 units.
function logTicks(lmin, lmax) {
  if (!isFinite(lmin) || !isFinite(lmax)) return [];
  const mant = lmax - lmin >= 3 ? [1] : lmax - lmin >= 1.5 ? [1, 3] : [1, 2, 5];
  const ticks = [];
  for (let k = Math.floor(lmin); k <= Math.ceil(lmax); k++) {
    for (const s of mant) ticks.push(k + Math.log10(s));
  }
  return ticks;
}

const SUPERSCRIPT = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };

function supNum(k) {
  return String(k).split("").map((ch) => SUPERSCRIPT[ch]).join("");
}

function logTickLabel(l) {
  const k = Math.round(l);
  if (Math.abs(l - k) > 1e-9) return tickLabel(Math.pow(10, l));  // 2, 5, ... ticks
  return "10" + supNum(k);
}

// ---------------------------------------------------------------------------
// Native MathML for the formulas in the panel titles. No dependency, and
// msubsup stacks the subscript under the exponent the way TeX's R_{200}^{3}
// does, which a single SVG <text> run cannot do.
// ---------------------------------------------------------------------------
const THIN = `<mspace width="0.19em"/>`;
const eq = (body) => `<math display="block">${body}</math>`;
const mi = (s) => `<mi>${s}</mi>`;
const msub = (base, sub) => `<msub>${mi(base)}<mn>${sub}</mn></msub>`;
const msup = (base, sup) => `<msup>${base}<mn>${sup}</mn></msup>`;
const msubsup = (base, sub, sup) => `<msubsup>${mi(base)}<mn>${sub}</mn><mn>${sup}</mn></msubsup>`;
const frac = (numer, denom) => `<mfrac><mrow>${numer}</mrow><mrow>${denom}</mrow></mfrac>`;
const paren = (head, inner) =>
  `<mrow>${head}<mo stretchy="false">(</mo>${mi(inner)}<mo stretchy="false">)</mo></mrow>`;

// A log axis carrying 2/3/5 mantissa ticks labels every tick as a plain number,
// so that it reads "3, 10, 30, 100" rather than mixing "3" with "10²".
function plainLogLabel(l) {
  const v = Math.pow(10, l);
  if (v >= 1e5 || v < 1e-3) return logTickLabel(l);
  return String(Number(v.toPrecision(3)));
}

// Typeset a positive number in scientific notation, e.g. 3.2×10¹¹.
function sci(v, d = 2) {
  const e = Math.floor(Math.log10(v));
  const a = Number((v / Math.pow(10, e)).toFixed(d));
  return `${a}×10${supNum(e)}`;
}

// Typeset a number with a proper minus sign.
function num(v, d) {
  return (d === undefined ? String(v) : v.toFixed(d)).replace("-", "−");
}
// Turn an axis spec into a concrete axis. `min`/`max`/`ticks` are optional:
// given explicitly they are used as-is, otherwise they follow the data.
function makeAxis(spec, vals) {
  const log = !!spec.log;
  const to = (v) => (log ? Math.log10(v) : v);
  const auto = spec.min === undefined || spec.max === undefined;
  const seen = vals.map(to);
  let min = spec.min !== undefined ? to(spec.min) : Math.min(...seen);
  let max = spec.max !== undefined ? to(spec.max) : Math.max(...seen);
  if (!isFinite(min) || !isFinite(max)) { min = 0; max = log ? 3 : 1; }  // degenerate input (z = 0)
  else if (min === max) { min -= 1; max += 1; }
  else if (auto) { const pad = (max - min) * 0.08; min -= pad; max += pad; }
  const ticks = (spec.ticks ? spec.ticks.map(to) : log ? logTicks(min, max) : niceTicks(min, max, 4))
    .filter((t) => t >= min && t <= max);
  const subDecade = log && ticks.some((t) => Math.abs(t - Math.round(t)) > 1e-9);
  const fmt = !log ? tickLabel : subDecade ? plainLogLabel : logTickLabel;
  return { log, to, min, max, ticks, label: spec.label, fmt };
}

// Is a value plottable on an axis of this kind?
function onAxis(v, log) {
  return isFinite(v) && (!log || v > 0);
}

// A single panel: title on top, both axes labelled, optionally a second y axis
// on the right that is the left one shifted by a constant. `opts.index` names
// the axis whose values form the uniform grid the hover snaps to.
function linePlot(title, pts, opts) {
  const { W, H, padL, padT, padB } = PLOT;
  const padR = opts.padR ?? PLOT.padR;
  const X = makeAxis(opts.x, pts.map((p) => p[0]).filter((v) => onAxis(v, opts.x.log)));
  const Y = makeAxis(opts.y, pts.map((p) => p[1]).filter((v) => onAxis(v, opts.y.log)));
  const sxT = (t) => padL + ((t - X.min) / (X.max - X.min)) * (W - padL - padR);
  const syT = (t) => H - padB - ((t - Y.min) / (Y.max - Y.min)) * (H - padT - padB);
  const sx = (v) => sxT(X.to(v));
  const sy = (v) => syT(Y.to(v));

  let d = "", started = false;
  for (const [x, y] of pts) {
    if (!onAxis(x, X.log) || !onAxis(y, Y.log)) { started = false; continue; }
    d += (started ? "L" : "M") + sx(x).toFixed(1) + " " + sy(y).toFixed(1) + " ";
    started = true;
  }

  let g = "";
  for (const t of Y.ticks) {
    const y = syT(t).toFixed(1);
    g += `<line class="grid" x1="${padL}" y1="${y}" x2="${W - padR}" y2="${y}"/>`;
    g += `<text class="tick" x="${padL - 5}" y="${(+y + 3.5).toFixed(1)}" text-anchor="end">${Y.fmt(t)}</text>`;
  }
  for (const t of X.ticks) {
    g += `<text class="tick" x="${sxT(t).toFixed(1)}" y="${H - padB + 16}" text-anchor="middle">${X.fmt(t)}</text>`;
  }

  let mk = "";
  if (opts.marker && opts.marker.x >= X.min && opts.marker.x <= X.max) {
    const x = sx(opts.marker.x).toFixed(1);
    const flip = opts.marker.x > X.min + 0.7 * (X.max - X.min);
    mk += `<line x1="${x}" y1="${padT}" x2="${x}" y2="${H - padB}" stroke="var(--marker)" stroke-width="1.1" stroke-dasharray="4 3"/>`;
    mk += `<text class="ztag" x="${flip ? +x - 4 : +x + 4}" y="${padT + 10}" text-anchor="${flip ? "end" : "start"}">${opts.marker.label}</text>`;
    const v = opts.marker.value;
    if (onAxis(v, Y.log) && Y.to(v) >= Y.min && Y.to(v) <= Y.max) {
      mk += `<circle cx="${x}" cy="${sy(v).toFixed(1)}" r="3.2" fill="var(--marker)"/>`;
    }
  }
  if (opts.tag) {
    mk += `<text class="ztag" x="${W - padR - 4}" y="${padT + 10}" text-anchor="end">${opts.tag}</text>`;
  }

  const xMid = (padL + W - padR) / 2;
  const yMid = (padT + H - padB) / 2;

  // Second y axis: the same quantity shifted by a constant (M_UV -> m_UV).
  let right = "";
  if (opts.yRight && isFinite(opts.yRight.offset)) {
    const { offset, label } = opts.yRight;
    right += `<line class="axis" x1="${W - padR}" y1="${padT}" x2="${W - padR}" y2="${H - padB}"/>`;
    for (const t of niceTicks(Y.min + offset, Y.max + offset, 4)) {
      if (t - offset < Y.min || t - offset > Y.max) continue;
      const y = syT(t - offset).toFixed(1);
      right += `<line class="axis" x1="${W - padR}" y1="${y}" x2="${W - padR + 4}" y2="${y}"/>`;
      right += `<text class="tick" x="${W - padR + 7}" y="${(+y + 3.5).toFixed(1)}" text-anchor="start">${tickLabel(t)}</text>`;
    }
    right += `<text class="axtitle" transform="rotate(90 ${W - 9} ${yMid})" x="${W - 9}" y="${yMid}" text-anchor="middle">${label}</text>`;
  }

  const html = `<figure class="panel"><div class="ptitle"><span>${title}</span></div><svg viewBox="0 0 ${W} ${H}">
    ${g}
    <line class="axis" x1="${padL}" y1="${padT}" x2="${padL}" y2="${H - padB}"/>
    <line class="axis" x1="${padL}" y1="${H - padB}" x2="${W - padR}" y2="${H - padB}"/>
    <path d="${d}" fill="none" stroke="var(--accent)" stroke-width="1.6"/>
    ${mk}
    ${right}
    <g class="hover" style="display:none">
      <line class="hline"/>
      <circle class="hdot" r="3"/>
      <text class="htext"></text>
      <text class="htext"></text>
      <text class="htext"></text>
    </g>
    <rect class="capture" x="${padL}" y="${padT}" width="${W - padL - padR}" height="${H - padT - padB}" fill="transparent"/>
    <text class="axtitle" x="${xMid}" y="${H - 6}" text-anchor="middle">${X.label}</text>
    <text class="axtitle" transform="rotate(-90 15 ${yMid})" x="15" y="${yMid}" text-anchor="middle">${Y.label}</text>
  </svg></figure>`;
  return { html, meta: { X, Y, padR, pts, index: opts.index === "y" ? 1 : 0 } };
}

function attachHover(svg, meta, opts) {
  const { W, H, padL, padT, padB } = PLOT;
  const { X, Y, padR, pts, index } = meta;
  const sx = (v) => padL + ((X.to(v) - X.min) / (X.max - X.min)) * (W - padL - padR);
  const sy = (v) => H - padB - ((Y.to(v) - Y.min) / (Y.max - Y.min)) * (H - padT - padB);
  const hg = svg.querySelector(".hover");
  const hline = svg.querySelector(".hline");
  const hdot = svg.querySelector(".hdot");
  const texts = [...svg.querySelectorAll(".htext")];
  // The points are a uniform grid along the indexed axis in that axis's own
  // units (log10 of the value for a log axis), so snapping is one division.
  const A = index === 0 ? X : Y;
  const v0 = A.to(pts[0][index]);
  const vN = A.to(pts[pts.length - 1][index]);
  const step = pts.length > 1 ? A.to(pts[1][index]) - v0 : 0.05;
  const lo = Math.min(v0, vN), hi = Math.max(v0, vN);

  // Position of the pointer along the indexed axis, in axis units.
  function valueAt(evt) {
    const rect = svg.getBoundingClientRect();
    const u = index === 0
      ? A.min + ((((evt.clientX - rect.left) / rect.width) * W - padL) / (W - padL - padR)) * (A.max - A.min)
      : A.min + ((H - padB - ((evt.clientY - rect.top) / rect.height) * H) / (H - padT - padB)) * (A.max - A.min);
    return Math.max(lo, Math.min(hi, u));
  }

  function move(evt) {
    const i = Math.max(0, Math.min(pts.length - 1, Math.round((valueAt(evt) - v0) / step)));
    const [xv, yv] = pts[i];
    if (!onAxis(xv, X.log) || !onAxis(yv, Y.log)) { hg.style.display = "none"; return; }
    const x = sx(xv), y = sy(yv);
    hdot.setAttribute("cx", x);
    hdot.setAttribute("cy", y);
    if (index === 0) {                       // vertical guide, labels pinned to the top
      hline.setAttribute("x1", x); hline.setAttribute("x2", x);
      hline.setAttribute("y1", padT); hline.setAttribute("y2", H - padB);
    } else {                                 // horizontal guide, labels ride with it
      hline.setAttribute("x1", padL); hline.setAttribute("x2", W - padR);
      hline.setAttribute("y1", y); hline.setAttribute("y2", y);
    }
    const lines = opts.texts(xv, yv);
    const flip = x > (padL + W - padR) / 2;
    const top = index === 0 ? padT + 13 : Math.min(y + 15, H - padB - 4 - 13 * (lines.length - 1));
    texts.forEach((el, k) => {
      const html = lines[k] || "";
      el.innerHTML = html;
      el.setAttribute("x", flip ? x - 5 : x + 5);
      el.setAttribute("y", top + 13 * k);
      el.setAttribute("text-anchor", flip ? "end" : "start");
    });
    hg.style.display = "";
  }

  svg.addEventListener("mousemove", move);
  svg.addEventListener("mouseleave", () => { hg.style.display = "none"; });
  if (opts.onClick) {
    svg.addEventListener("click", (evt) => {
      const u = valueAt(evt);
      opts.onClick(A.log ? Math.pow(10, u) : u);
    });
  }
}

function renderPlots(m, z, r) {
  const it = (s) => `<tspan font-style="italic">${s}</tspan>`;
  // subscript, followed by text that returns to the baseline
  const sb = (s, after) => `<tspan dy="0.28em" font-size="0.72em">${s}</tspan><tspan dy="-0.28em">${after}</tspan>`;

  const zCurve = [];
  for (let zz = 0.05; zz <= 20.0001; zz += 0.05) {
    const c = compute(zz, m);
    zCurve.push({ z: zz, scale: c.kpcPerArcsec, dm: c.distmod, age: c.ageAtZ });
  }
  // Apparent AB magnitude at the current redshift: m_UV = M_UV + (m_UV - M_UV).
  const fluxCurve = [];
  for (let i = 0; i <= Math.round((MUV_MAX - MUV_MIN) / MUV_STEP); i++) {
    const M = MUV_MIN + i * MUV_STEP;
    fluxCurve.push([AB_ZP_NJY * Math.pow(10, -0.4 * (M + r.distmod)), M]);
  }

  const zAxis = { min: 0, max: 20, ticks: [0, 5, 10, 15, 20], label: it("z") };
  const zOpts = {
    x: zAxis, index: "x",
    onClick: (x) => { zInput.value = Number(x.toFixed(2)); render(); },
  };
  const zMark = (value) => ({ x: z, value, label: `${it("z")} = ${Number(z.toPrecision(3))}` });
  const zText = (x) => `${it("z")} = ${x.toFixed(1)}`;
  const zTag = `${it("z")} = ${Number(z.toPrecision(3))}`;

  // Row 1: the two z panels; row 2: magnitude and the flux it implies.
  const panels = [
    {
      title: "Angular scale",
      pts: zCurve.map((c) => [c.z, c.scale]),
      opts: { ...zOpts, y: { label: "Angular scale (kpc/″)" }, marker: zMark(r.kpcPerArcsec) },
      texts: (x, y) => [zText(x), `${y.toFixed(1)} kpc/″`],
    },
    {
      title: "Age",
      pts: zCurve.map((c) => [c.z, c.age]),
      opts: { ...zOpts, y: { label: "Age (Gyr)" }, marker: zMark(r.ageAtZ) },
      texts: (x, y) => [zText(x), `${y.toFixed(2)} Gyr`],
    },
    {
      title: "<i>m</i><sub>UV</sub> − <i>M</i><sub>UV</sub>",
      pts: zCurve.map((c) => [c.z, c.dm]),
      opts: { ...zOpts, y: { label: `${it("m")}${sb("UV", " − ")}${it("M")}${sb("UV", " (mag)")}` }, marker: zMark(r.distmod) },
      texts: (x, y) => [zText(x), `${y.toFixed(2)} mag`],
    },
    {
      title: "Flux density, <i>f</i><sub>ν</sub>",
      pts: fluxCurve,
      opts: {
        x: { log: true, label: `${it("f")}${sb("ν", " (nJy)")}` },
        y: { min: MUV_MIN, max: MUV_MAX, ticks: [-24, -22, -20, -18, -16],
             label: `${it("M")}${sb("UV", " (mag)")}` },
        yRight: { offset: r.distmod, label: `${it("m")}${sb("UV", " (mag)")}` },
        padR: 46, index: "y", tag: zTag,
      },
      texts: (x, y) => [
        `${it("M")}${sb("UV", " = " + num(y, 1))}`,
        `${it("m")}${sb("UV", " = " + num(y + r.distmod, 2))}`,
        `${fmtFlux(x)} nJy`,
      ],
    },
  ];

  const metas = [];
  plotsEl.innerHTML = panels.map((p) => {
    const { html, meta } = linePlot(p.title, p.pts, p.opts);
    metas.push(meta);
    return html;
  }).join("");

  plotsEl.querySelectorAll(".panel svg").forEach((svg, i) => {
    attachHover(svg, metas[i], { texts: panels[i].texts, onClick: panels[i].opts.onClick });
  });
}

// R200 and V200 against M200 at the current redshift; both axes logarithmic.
const LOG_M_MIN = 9, LOG_M_MAX = 14, LOG_M_STEP = 0.01;

function renderHaloPlots(m, z, r) {
  const it = (s) => `<tspan font-style="italic">${s}</tspan>`;
  const sb = (s, after) => `<tspan dy="0.28em" font-size="0.72em">${s}</tspan><tspan dy="-0.28em">${after}</tspan>`;
  const M200 = (after) => `${it("M")}${sb("200", after)}`;
  const R200 = (after) => `${it("R")}${sb("200", after)}`;
  const V200 = (after) => `${it("V")}${sb("200", after)}`;

  const Hz = m.H0 * r.Ez;
  const curve = [];
  for (let i = 0; i <= Math.round((LOG_M_MAX - LOG_M_MIN) / LOG_M_STEP); i++) {
    const M = Math.pow(10, LOG_M_MIN + i * LOG_M_STEP);
    curve.push({ M, ...halo(M, Hz) });
  }

  const mAxis = {
    log: true, min: Math.pow(10, LOG_M_MIN), max: Math.pow(10, LOG_M_MAX),
    label: `${M200(" (")}${it("M")}${sb("⊙", ")")}`,
  };
  const zTag = `${it("z")} = ${Number(z.toPrecision(3))}`;
  const sig3 = (x) => Number(x.toPrecision(3)).toLocaleString();
  const texts = (M) => {
    const c = curve[Math.round((Math.log10(M) - LOG_M_MIN) / LOG_M_STEP)];
    // R200 also as an angle on the sky, via the angular scale at this redshift.
    // At z = 0 the angular scale vanishes, so there is no angle to quote.
    const arcsec = c.R200 / r.kpcPerArcsec;
    const angle = isFinite(arcsec) ? ` (${sig3(arcsec)}″)` : "";
    return [
      `${M200(" = " + sci(c.M) + " ")}${it("M")}${sb("⊙", "")}`,
      R200(" = " + sig3(c.R200) + " kpc" + angle),
      V200(" = " + sig3(c.V200) + " km/s"),
    ];
  };

  const panels = [
    {
      // M200 = 100 H(z)^2 R200^3 / G
      title: eq(`${msub("M", 200)}<mo>=</mo><mn>100</mn>${THIN}` +
        `${msup(paren(mi("H"), "z"), 2)}${THIN}${msubsup("R", 200, 3)}` +
        `<mo>/</mo>${mi("G")}`),
      pts: curve.map((c) => [c.M, c.R200]),
      opts: { x: mAxis, y: { log: true, label: R200(" (kpc)") }, index: "x", tag: zTag },
    },
    {
      // M200 = V200^3 / (10 G H(z))
      title: eq(`${msub("M", 200)}<mo>=</mo>` +
        frac(msubsup("V", 200, 3), `<mn>10</mn>${THIN}${mi("G")}${THIN}${paren(mi("H"), "z")}`)),
      pts: curve.map((c) => [c.M, c.V200]),
      opts: { x: mAxis, y: { log: true, label: V200(" (km s⁻¹)") }, index: "x", tag: zTag },
    },
  ];

  const metas = [];
  haloPlotsEl.innerHTML = panels.map((p) => {
    const { html, meta } = linePlot(p.title, p.pts, p.opts);
    metas.push(meta);
    return html;
  }).join("");

  haloPlotsEl.querySelectorAll(".panel svg").forEach((svg, i) => {
    attachHover(svg, metas[i], { texts });
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
    row("Angular size scale",
      `${sig3(r.kpcPerArcsec)} <span class="unit">pkpc/arcsec</span>` +
      `<span class="alt">${sig3(r.cMpcPerArcmin)} <span class="unit">cMpc/arcmin</span></span>`),
    row(`<i>m</i><sub>UV</sub> &minus; <i>M</i><sub>UV</sub> <a class="info" href="distance-modulus.html" target="_blank" rel="noopener" title="Includes the &minus;2.5 log(1+z) band-shift term. Why? Click to read.">&#9432;</a>`, sig4(r.distmod), "mag"),
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
  renderHaloPlots(m, z, r);
  resultsEl.hidden = false;
}

form.addEventListener("submit", (e) => {
  e.preventDefault();
  render();
});
cosmoSel.addEventListener("change", render);

// Initial render on load.
render();
