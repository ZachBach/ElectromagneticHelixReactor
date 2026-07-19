/* Phase 3A — quasi-static axisymmetric Poisson solver on a cell-centered
   (r,z) grid:

     (1/r) d/dr ( r dphi/dr ) + d²phi/dz² = -rho/eps0

   Finite volumes with SOR. Boundaries are grounded conductors (Dirichlet
   phi = 0 at r = R and both end plates), applied at cell faces via ghost
   folding. The r = 0 axis is regular by construction: the r_{i-1/2} face
   coefficient vanishes for the first cell, so no axis condition is needed.
   Validated in selftest.js against a manufactured solution.               */
'use strict';

class PoissonRZ {
  constructor(Nr, Nz, R, L) {
    this.Nr = Nr; this.Nz = Nz; this.R = R; this.L = L;
    this.dr = R / Nr; this.dz = L / Nz;
    this.phi = new Float64Array(Nr * Nz);
    this.Er = new Float64Array(Nr * Nz);
    this.Ez = new Float64Array(Nr * Nz);
    // face-coefficient arrays: aW[i]·phiW + aE[i]·phiE + az·(phiS+phiN)
    this.aW = new Float64Array(Nr);
    this.aE = new Float64Array(Nr);
    for (let i = 0; i < Nr; i++) {
      const rc = (i + 0.5) * this.dr;
      this.aW[i] = i / (rc * this.dr);           // r_{i-1/2}/(rc·dr²)·dr = i/(rc·dr)
      this.aE[i] = (i + 1) / (rc * this.dr);
    }
    this.az = 1 / (this.dz * this.dz);
  }

  /* SOR solve of lap(phi) = -S with S = rho/eps0, warm-starting from the
     previous phi. Returns sweeps used. */
  solve(S, opts) {
    const { Nr, Nz, aW, aE, az, phi } = this;
    const omega = (opts && opts.omega) || 1.9;
    const maxIter = (opts && opts.maxIter) || 4000;
    const tol = (opts && opts.tol) || 1e-4;
    let sweep = 0;
    for (; sweep < maxIter; sweep++) {
      let maxD = 0, maxP = 1e-9;
      for (let j = 0; j < Nz; j++) {
        const jOff = j * Nr;
        for (let i = 0; i < Nr; i++) {
          const c = jOff + i;
          // Dirichlet phi=0 at faces: ghost = -phi_C folds into the diagonal
          let diag = aW[i] + aE[i] + 2 * az;
          let sum = S[c];
          if (i > 0) sum += aW[i] * phi[c - 1];
          if (i + 1 < Nr) sum += aE[i] * phi[c + 1]; else diag += aE[i];
          if (j > 0) sum += az * phi[c - Nr]; else diag += az;
          if (j + 1 < Nz) sum += az * phi[c + Nr]; else diag += az;
          const nu = (1 - omega) * phi[c] + (omega * sum) / diag;
          const d = Math.abs(nu - phi[c]);
          if (d > maxD) maxD = d;
          const a = Math.abs(nu);
          if (a > maxP) maxP = a;
          phi[c] = nu;
        }
      }
      if (maxD < tol * maxP) { sweep++; break; }
    }
    return sweep;
  }

  /* E = -grad phi, cell-centered, using the same ghost conventions
     (mirror at the axis, Dirichlet at wall and end plates). */
  gradients() {
    const { Nr, Nz, dr, dz, phi, Er, Ez } = this;
    for (let j = 0; j < Nz; j++) {
      const jOff = j * Nr;
      for (let i = 0; i < Nr; i++) {
        const c = jOff + i;
        const pW = i > 0 ? phi[c - 1] : phi[c];           // mirror across axis
        const pE = i + 1 < Nr ? phi[c + 1] : -phi[c];     // wall ghost
        Er[c] = -(pE - pW) / (2 * dr);
        const pS = j > 0 ? phi[c - Nr] : -phi[c];         // end-plate ghosts
        const pN = j + 1 < Nz ? phi[c + Nr] : -phi[c];
        Ez[c] = -(pN - pS) / (2 * dz);
      }
    }
  }
}

module.exports = { PoissonRZ };
