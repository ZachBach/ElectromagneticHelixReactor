/* Phase 2 seed — loss bookkeeping and (r,z) density deposition. */
'use strict';

const KIND = { RADIAL: 1, END_LOW: 2, END_HIGH: 3 };

/* Records each particle's loss time and channel; survivors are censored at
   Tmax, so with survivors present the mean is a lower bound (flagged). */
class LossStats {
  constructor(n) {
    this.times = new Float64Array(n);
    this.kinds = new Uint8Array(n);
    this.count = 0;
    this.cap = n;
  }
  record(t, kind) {
    this.times[this.count] = t;
    this.kinds[this.count] = kind;
    this.count++;
  }
  summary(nTotal, Tmax) {
    const n = this.count;
    const survivors = nTotal - n;
    let sum = 0, radial = 0, endLow = 0, endHigh = 0;
    for (let i = 0; i < n; i++) {
      sum += this.times[i];
      if (this.kinds[i] === KIND.RADIAL) radial++;
      else if (this.kinds[i] === KIND.END_LOW) endLow++;
      else endHigh++;
    }
    // restricted mean: survivors contribute Tmax (lower bound on true mean)
    const mean = (sum + survivors * Tmax) / nTotal;
    let var_ = 0;
    for (let i = 0; i < n; i++) { const d = this.times[i] - mean; var_ += d * d; }
    for (let i = 0; i < survivors; i++) { const d = Tmax - mean; var_ += d * d; }
    const se = Math.sqrt(var_ / nTotal / Math.max(nTotal - 1, 1));
    const sorted = this.times.slice(0, n).sort();
    // median over all particles incl. censored survivors at Tmax
    const mid = (nTotal / 2) | 0;
    const median = mid < n ? sorted[mid] : Tmax;
    return {
      nTotal, nLost: n, survivors,
      surviveFrac: survivors / nTotal,
      tauMean: mean, tauSE: se, tauMedian: median,
      fracRadial: n ? radial / n : 0,
      fracEnd: n ? (endLow + endHigh) / n : 0,
      fracEndLow: n ? endLow / n : 0,
      fracEndHigh: n ? endHigh / n : 0,
      censored: survivors > 0,
    };
  }
}

/* Axisymmetric density accumulator: deposit(x, y, z, w) adds weight w to the
   (r,z) cell; profile rows come back volume-normalized so uniform density
   reads flat in r. */
class DensityRZ {
  constructor(nr, nz, R, L) {
    this.nr = nr; this.nz = nz; this.R = R; this.L = L;
    this.grid = new Float64Array(nr * nz);
  }
  deposit(x, y, z, w) {
    const r = Math.sqrt(x * x + y * y);
    let ir = ((r / this.R) * this.nr) | 0;
    let iz = ((z / this.L) * this.nz) | 0;
    if (ir >= this.nr) ir = this.nr - 1;
    if (iz < 0) iz = 0; else if (iz >= this.nz) iz = this.nz - 1;
    this.grid[iz * this.nr + ir] += w;
  }
  // rows of [rCenter, zCenter, densityPerVolume]
  rows() {
    const out = [];
    const dz = this.L / this.nz, dr = this.R / this.nr;
    for (let iz = 0; iz < this.nz; iz++) {
      for (let ir = 0; ir < this.nr; ir++) {
        const r0 = ir * dr, r1 = r0 + dr;
        const vol = Math.PI * (r1 * r1 - r0 * r0) * dz;
        out.push([r0 + dr / 2, (iz + 0.5) * dz, this.grid[iz * this.nr + ir] / vol]);
      }
    }
    return out;
  }
}

module.exports = { KIND, LossStats, DensityRZ };
