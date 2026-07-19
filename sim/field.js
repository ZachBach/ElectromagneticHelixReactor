/* Module 1 — static helical magnetic field: B = Bz(r) ẑ + Bθ(r) θ̂, Br = 0
   (identically divergence-free). Two models, matching the web app:

   'screw'    — Bz constant, Bθ ∝ r (uniform axial current density).
   'beltrami' — force-free constant-α field, Bz = B0·J0(αr), Bθ = B0·J1(αr),
                which satisfies ∇×B = αB.

   The scan variable is ratio = Bθ(R)/Bz at the wall. The field is normalized so
   that |B(R)| = Bwall for every ratio: a pitch scan changes geometry at constant
   field magnitude at the wall radius, which is the closest divergence-free
   realization of the paper's "constant |B|" condition in Prediction 4 (a field
   with |B| literally constant everywhere and tunable pitch does not exist for
   these profiles). */
'use strict';

// Bessel J0/J1 by series — accurate for x ≲ 8, and here x < 2.405 always.
function besselJ(n, x) {
  let term = 1;
  for (let k = 1; k <= n; k++) term *= (x / 2) / k;
  let sum = term;
  const x2 = -(x * x) / 4;
  for (let k = 1; k < 40; k++) {
    term *= x2 / (k * (k + n));
    sum += term;
    if (Math.abs(term) < 1e-14) break;
  }
  return sum;
}

/* makeField({ model, Bwall, ratio, R }) → {
     at(x, y, z, out)  — fills out[0..2] with (Bx, By, Bz), T
     pitchAt(r)        — local field-line pitch L_p = 2πr·Bz/Bθ, m
     Bmax              — upper bound on |B| in the volume (for timestep choice)
     meta              — resolved parameters
   } */
function makeField(opts) {
  const model = opts.model || 'screw';
  const Bwall = opts.Bwall;      // |B| at r = R, tesla
  const ratio = opts.ratio;      // Bθ(R)/Bz(R), dimensionless (≥ 0)
  const R = opts.R;              // wall radius, m

  if (model === 'screw') {
    const Bz = Bwall / Math.sqrt(1 + ratio * ratio);
    const BthWall = ratio * Bz;
    const k = BthWall / R;       // Bθ(r)/r, constant
    return {
      at(x, y, z, out) {
        out[0] = -y * k;
        out[1] = x * k;
        out[2] = Bz;
        return out;
      },
      pitchAt(r) {
        const bth = k * r;
        return bth > 0 ? (2 * Math.PI * r * Bz) / bth : Infinity;
      },
      Bmax: Bwall, // |B|² = Bz² + (kr)² is maximal at r = R
      meta: { model, Bwall, ratio, R, Bz, BthWall },
    };
  }

  if (model === 'beltrami') {
    // Solve J1(αR)/J0(αR) = ratio on (0, j_{0,1}) — monotone increasing.
    let lo = 1e-9, hi = 2.4048;
    for (let i = 0; i < 60; i++) {
      const mid = (lo + hi) / 2;
      if (besselJ(1, mid) / besselJ(0, mid) < ratio) lo = mid; else hi = mid;
    }
    const aR = (lo + hi) / 2;
    const alpha = aR / R;
    const j0R = besselJ(0, aR), j1R = besselJ(1, aR);
    const B0 = Bwall / Math.sqrt(j0R * j0R + j1R * j1R);
    return {
      at(x, y, z, out) {
        const r = Math.sqrt(x * x + y * y);
        const ar = alpha * r;
        // Bθ/r → B0·α/2 as r → 0
        const bthOverR = r < 1e-9 ? B0 * alpha / 2 : (B0 * besselJ(1, ar)) / r;
        out[0] = -y * bthOverR;
        out[1] = x * bthOverR;
        out[2] = B0 * besselJ(0, ar);
        return out;
      },
      pitchAt(r) {
        const bth = besselJ(1, alpha * r), bz = besselJ(0, alpha * r);
        return bth > 0 ? (2 * Math.PI * r * bz) / bth : Infinity;
      },
      Bmax: B0, // J0² + J1² decreases from 1 at the axis
      meta: { model, Bwall, ratio, R, alpha, B0 },
    };
  }

  throw new Error('unknown field model: ' + model);
}

module.exports = { makeField, besselJ };
