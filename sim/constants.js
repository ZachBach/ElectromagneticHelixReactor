/* Physical constants (SI) and argon collision data for the Phase 1 pitch-physics
   simulator. Unlike the browser engine (normalized units, tuned coefficients),
   everything here is dimensional so confinement times come out in seconds. */
'use strict';

const QE = 1.602176634e-19;    // elementary charge, C
const EPS0 = 8.8541878128e-12; // vacuum permittivity, F/m
const ME = 9.1093837015e-31;   // electron mass, kg
const AMU = 1.66053906660e-27; // atomic mass unit, kg
const MAR = 39.948 * AMU;      // argon mass, kg
const KB = 1.380649e-23;       // Boltzmann constant, J/K
const EV = QE;                 // 1 eV in J

// Argon ionization threshold, eV.
const EION = 15.76;

/* Electron-argon momentum-transfer cross section, [eV, m²].
   Approximate log-log table capturing the Ramsauer-Townsend minimum near
   0.25 eV and the broad peak near 15 eV, anchored to the representative value
   used in the paper (2e-20 m² at a few eV). Adequate for the relative pitch
   scan; swap in LXCat (Phelps/Biagi) tables when absolute rates matter. */
const SIGMA_EL_TABLE = [
  [0.01, 4.0e-20], [0.1, 6.0e-21], [0.25, 1.5e-21], [1.0, 8.0e-21],
  [3.0, 2.0e-20], [10.0, 1.0e-19], [15.0, 1.5e-19], [50.0, 7.0e-20],
  [100.0, 4.0e-20], [1000.0, 1.0e-20],
];

/* Electron-impact ionization cross section of argon, [eV, m²].
   Approximate: rises from threshold to ~2.9e-20 m² near 100 eV. */
const SIGMA_ION_TABLE = [
  [16.0, 2.0e-22], [20.0, 3.0e-21], [30.0, 1.4e-20], [50.0, 2.3e-20],
  [100.0, 2.9e-20], [200.0, 2.4e-20], [1000.0, 1.2e-20],
];

/* Ar+ + Ar: charge exchange and elastic, treated as energy-independent at
   low energy. Combined magnitude consistent with the paper's σ_in ≈ 1e-18 m². */
const SIGMA_CX = 5.0e-19;
const SIGMA_I_EL = 5.0e-19;

// Log-log interpolation with flat extrapolation at both ends.
function lookupLogLog(table, E) {
  if (E <= table[0][0]) return table[0][1];
  const last = table.length - 1;
  if (E >= table[last][0]) return table[last][1];
  let i = 1;
  while (table[i][0] < E) i++;
  const [e0, s0] = table[i - 1], [e1, s1] = table[i];
  const f = (Math.log(E) - Math.log(e0)) / (Math.log(e1) - Math.log(e0));
  return Math.exp(Math.log(s0) + f * (Math.log(s1) - Math.log(s0)));
}

function sigmaElectronElastic(E) { return lookupLogLog(SIGMA_EL_TABLE, E); }
function sigmaElectronIonization(E) {
  return E > EION ? lookupLogLog(SIGMA_ION_TABLE, E) : 0;
}

// Neutral density from fill pressure (mTorr) and gas temperature (K).
function neutralDensity(pmTorr, Tgas) {
  return (pmTorr * 0.133322) / (KB * Tgas); // 1 mTorr = 0.133322 Pa
}

module.exports = {
  QE, EPS0, ME, AMU, MAR, KB, EV, EION,
  SIGMA_CX, SIGMA_I_EL,
  sigmaElectronElastic, sigmaElectronIonization, neutralDensity,
};
