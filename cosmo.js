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
  planck18: { name: "Planck 2018", H0: 67.36, OmegaM: 0.3153, OmegaL: 0.6847 },
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
    comovingVolumeGpc3: VC / 1e9, // Mpc^3 -> Gpc^3
    distmod: 5 * Math.log10((DL * 1e6) / 10), // distance modulus
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
const tbody = document.querySelector("#results-table tbody");

function row(label, value, unit) {
  const tr = document.createElement("tr");
  const u = unit ? ` <span class="unit">${unit}</span>` : "";
  tr.innerHTML = `<td class="label">${label}</td><td class="value">${value}${u}</td>`;
  return tr;
}

function showParams(m) {
  paramsEl.innerHTML =
    `<code>H&#8320; = ${m.H0}</code>` +
    `<code>&Omega;<sub>M</sub> = ${m.OmegaM}</code>` +
    `<code>&Omega;<sub>&Lambda;</sub> = ${m.OmegaL}</code>` +
    `<code>&Omega;<sub>R</sub> = ${m.OmegaR.toExponential(3)}</code>` +
    `<code>&Omega;<sub>k</sub> = ${m.OmegaK.toExponential(2)}</code>` +
    `<code>flat &Lambda;CDM</code>`;
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
  tbody.innerHTML = "";
  tbody.append(
    row("Age of the universe now", fmt(r.ageNow), "Gyr"),
    row("Age at redshift z", fmt(r.ageAtZ), "Gyr"),
    row("Light travel time (lookback)", fmt(r.lookback), "Gyr"),
    row("Comoving radial distance", `${fmt(r.comovingRadialMpc, 5)}`, "Mpc"),
    row("&nbsp;&nbsp;&hookrightarrow; equivalently", `${fmt(r.comovingRadialGly)}`, "Gly"),
    row("Angular diameter distance D<sub>A</sub>", `${fmt(r.DA_Mpc, 5)}`, "Mpc"),
    row("&nbsp;&nbsp;&hookrightarrow; equivalently", `${fmt(r.DA_Gly)}`, "Gly"),
    row("Angular size scale", fmt(r.kpcPerArcsec), "kpc/&Prime;"),
    row("Luminosity distance D<sub>L</sub>", `${fmt(r.DL_Mpc, 5)}`, "Mpc"),
    row("&nbsp;&nbsp;&hookrightarrow; equivalently", `${fmt(r.DL_Gly)}`, "Gly"),
    row("Distance modulus", fmt(r.distmod), "mag"),
    row("Comoving volume within z", fmt(r.comovingVolumeGpc3), "Gpc&sup3;"),
    row("E(z) = H(z)/H&#8320;", fmt(r.Ez)),
  );
  resultsEl.hidden = false;
}

form.addEventListener("submit", (e) => {
  e.preventDefault();
  render();
});
cosmoSel.addEventListener("change", render);

// Initial render on load.
render();
