/* Module 3 — Monte Carlo collisions with the neutral argon background
   (test-particle MCC). The runner picks dt so that ν_max·dt ≪ 1, which lets us
   sample directly with P = ν(E)·dt per step instead of the null-collision
   machinery.

   Electrons: elastic (isotropic redirect using the momentum-transfer cross
   section, exact 2m/M·(1−cosθ) fractional energy loss) and ionization
   (threshold subtracted, random energy split with the untracked secondary,
   isotropic redirect; events are counted for later phases).

   Ions (Ar+): charge exchange (ion velocity replaced by a neutral Maxwellian
   sample — the fast neutral leaves, a slow ion remains) and equal-mass
   hard-sphere elastic scattering against a sampled Maxwellian partner. */
'use strict';

const C = require('./constants');

// Deterministic RNG (mulberry32) + cached Box-Muller gaussian.
function makeRng(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function makeGauss(rng) {
  let spare = null;
  return function () {
    if (spare !== null) { const v = spare; spare = null; return v; }
    let u, v, s;
    do { u = rng() * 2 - 1; v = rng() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt(-2 * Math.log(s) / s);
    spare = v * m;
    return u * m;
  };
}

const LUT_N = 512, E_MIN = 0.005, E_MAX = 2000; // eV, log grid

// Largest electron collision rate ν = n_n·σ_tot(E)·v(E) over the LUT range,
// used by the runner to bound ν·dt before building the collider.
function electronMaxRate(nn) {
  const lnMin = Math.log(E_MIN), span = Math.log(E_MAX) - lnMin;
  let max = 0;
  for (let i = 0; i < LUT_N; i++) {
    const E = Math.exp(lnMin + (span * i) / (LUT_N - 1));
    const v = Math.sqrt((2 * E * C.EV) / C.ME);
    const nu = nn * (C.sigmaElectronElastic(E) + C.sigmaElectronIonization(E)) * v;
    if (nu > max) max = nu;
  }
  return max;
}

/* Returns { collide(vel, i3) → 0 none | 1 elastic | 2 ionization }. */
function makeElectronCollider(nn, dt, rng) {
  const lnMin = Math.log(E_MIN), span = Math.log(E_MAX) - lnMin;
  const pEl = new Float64Array(LUT_N), pIon = new Float64Array(LUT_N);
  for (let i = 0; i < LUT_N; i++) {
    const E = Math.exp(lnMin + (span * i) / (LUT_N - 1));
    const v = Math.sqrt((2 * E * C.EV) / C.ME);
    pEl[i] = nn * C.sigmaElectronElastic(E) * v * dt;
    pIon[i] = nn * C.sigmaElectronIonization(E) * v * dt;
  }
  const idxScale = (LUT_N - 1) / span;
  const dm = (2 * C.ME) / C.MAR; // elastic energy-transfer coefficient

  function isotropic(vel, i3, speed) {
    const cth = rng() * 2 - 1, sth = Math.sqrt(1 - cth * cth), ph = rng() * 2 * Math.PI;
    vel[i3] = speed * sth * Math.cos(ph);
    vel[i3 + 1] = speed * sth * Math.sin(ph);
    vel[i3 + 2] = speed * cth;
    return cth; // caller uses dot with old direction, not this
  }

  return {
    collide(vel, i3) {
      const vx = vel[i3], vy = vel[i3 + 1], vz = vel[i3 + 2];
      const v2 = vx * vx + vy * vy + vz * vz;
      const E = (0.5 * C.ME * v2) / C.EV;
      let k = ((Math.log(E) - lnMin) * idxScale) | 0;
      if (k < 0) k = 0; else if (k > LUT_N - 1) k = LUT_N - 1;
      const pe = pEl[k], pi = pIon[k];
      const u = rng();
      if (u >= pe + pi) return 0;
      const vmag = Math.sqrt(v2);
      if (u < pe) {
        // elastic: isotropic new direction; energy loss uses the actual
        // scattering angle between old and new directions
        const ox = vx / vmag, oy = vy / vmag, oz = vz / vmag;
        isotropic(vel, i3, 1);
        const cosChi = ox * vel[i3] + oy * vel[i3 + 1] + oz * vel[i3 + 2];
        const speed = vmag * Math.sqrt(Math.max(1 - dm * (1 - cosChi), 0));
        vel[i3] *= speed; vel[i3 + 1] *= speed; vel[i3 + 2] *= speed;
        return 1;
      }
      // ionization: pay the threshold, split the remainder with the secondary
      const Eres = Math.max(E - C.EION, 0.01) * (0.5 + 0.5 * rng());
      isotropic(vel, i3, Math.sqrt((2 * Eres * C.EV) / C.ME));
      return 2;
    },
  };
}

// Upper bound on the ion collision rate for timestep choice: fastest ions the
// runner will see are a few × thermal; vRef supplies that scale.
function ionMaxRate(nn, vRef) {
  return nn * (C.SIGMA_CX + C.SIGMA_I_EL) * vRef;
}

/* Returns { collide(vel, i3) → 0 none | 1 elastic | 2 charge exchange }.
   mIon (kg) sets the ion AND neutral mass — CX partners are the same species
   (defaults to argon; reduced-mass model gases pass a lighter value). */
function makeIonCollider(nn, Tgas, dt, rng, gauss, mIon) {
  const m = mIon || C.MAR;
  const vthn = Math.sqrt((C.KB * Tgas) / m); // per-component neutral thermal
  const sigTot = C.SIGMA_CX + C.SIGMA_I_EL;
  const pcx = C.SIGMA_CX / sigTot;
  return {
    collide(vel, i3) {
      const vx = vel[i3], vy = vel[i3 + 1], vz = vel[i3 + 2];
      // effective relative speed: ion speed combined with the neutral thermal
      // spread (partner is only sampled when a collision actually fires)
      const g = Math.sqrt(vx * vx + vy * vy + vz * vz + 3 * vthn * vthn);
      if (rng() >= nn * sigTot * g * dt) return 0;
      const wx = gauss() * vthn, wy = gauss() * vthn, wz = gauss() * vthn;
      if (rng() < pcx) {
        // charge exchange: the ion identity jumps to the neutral
        vel[i3] = wx; vel[i3 + 1] = wy; vel[i3 + 2] = wz;
        return 2;
      }
      // equal-mass hard-sphere elastic: isotropize the relative velocity in COM
      const cmx = (vx + wx) / 2, cmy = (vy + wy) / 2, cmz = (vz + wz) / 2;
      const gx = vx - wx, gy = vy - wy, gz = vz - wz;
      const gm = Math.sqrt(gx * gx + gy * gy + gz * gz);
      const cth = rng() * 2 - 1, sth = Math.sqrt(1 - cth * cth), ph = rng() * 2 * Math.PI;
      vel[i3] = cmx + (gm / 2) * sth * Math.cos(ph);
      vel[i3 + 1] = cmy + (gm / 2) * sth * Math.sin(ph);
      vel[i3 + 2] = cmz + (gm / 2) * cth;
      return 1;
    },
  };
}

module.exports = { makeRng, makeGauss, electronMaxRate, makeElectronCollider, ionMaxRate, makeIonCollider };
