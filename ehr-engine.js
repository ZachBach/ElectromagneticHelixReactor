/* EHR Simulator — reduced-physics engine.
   Normalized units: chamber radius R=1, length z∈[-2,2] (periodic), q/m=-1 for electrons,
   B~O(1) so gyro period ~2π. Energy bookkeeping: E[eV] = 5·|v|².
   Reduced model notes:
   - B = Bz ẑ + Bθ(r) θ̂ with Bθ ∝ r (uniform axial current density).
   - RF heating: stochastic kicks along B̂, spatially patterned by the helical antenna
     wave cos(mθ − k_z z − ω t). Coupling envelope peaks at M=1 (the core hypothesis,
     encoded as antenna/field-line alignment) and shifts deposition edge→core with coupling.
   - Monte Carlo collisions vs neutrals; nn(z,t) carries the acoustic modulation.
   - Ionization events deposit into a 64×64×128 density grid with volumetric loss (decay),
     so the steady-state grid maps the ionization source distribution. */
(function () {
  'use strict';
  const TWO_PI = Math.PI * 2;
  const NX = 64, NY = 64, NZ = 128;
  const EVK = 5.0;          // eV per v²
  const EION = 15.8;        // argon ionization threshold, eV
  const DT = 0.085;         // substep (normalized time)
  const SUBSTEPS = 3;
  const SIGMA_EL = 0.12, SIGMA_ION = 0.42;
  const KICK = 0.55;        // heating kick scale
  const DECAY_TAU = 1.1;    // density grid loss time, real seconds

  function clamp(v, a, b) { return v < a ? a : v > b ? b : v; }

  let _g = null;
  function gauss() {
    if (_g !== null) { const v = _g; _g = null; return v; }
    let u = 0, v = 0, s = 0;
    do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
    const m = Math.sqrt(-2 * Math.log(s) / s);
    _g = v * m; return u * m;
  }

  class EHREngine {
    constructor(N) {
      this.params = {
        N: N || 18000,
        rfPower: 42, rfFreq: 13.56, mMode: 1, kzA: 1.6,
        Bz: 1.0, btRatio: 0.62,
        pressure: 3.0, acA: 0.3, acLambda: 1.3, acFreq: 0.4,
        couplingMode: 'kernel', fieldModel: 'screw', alpha: 1.2, acMode: 'traveling',
        bias: 80, dustN: 300, dustQ: 4, dustLam: 0.06,
      };
      this.geomDirty = true;
      this.grid = new Float32Array(NX * NY * NZ);
      this.gridMax = 0.001;
      // precomputed radius per (ix,iy) voxel column
      this.rTab = new Float32Array(NX * NY);
      let maskCount = 0;
      for (let iy = 0; iy < NY; iy++) for (let ix = 0; ix < NX; ix++) {
        const x = (ix + 0.5) / NX * 2 - 1, y = (iy + 0.5) / NY * 2 - 1;
        const r = Math.sqrt(x * x + y * y);
        this.rTab[iy * NX + ix] = r;
        if (r <= 1) maskCount++;
      }
      this.maskCount = maskCount;
      // diagnostics buffers
      this.radial = new Float32Array(28); this._radCnt = new Float32Array(28);
      this.axial = new Float32Array(64);
      this.gpi = new Float32Array(NX * NY);
      this.diag = { nCore: 0, nEdge: 0, Rec: 0, nMean: 0, dnn: 0, ionRate: 0, meanE: 0, resTime: 0, refl: 0, waferUnif: 0, ionFlux: 0, ionEw: 0, mFit: 0, mSnr: 0, acCorr: 0 };
      // history rings (10 Hz)
      const CAP = 1200;
      this.hist = { cap: CAP, len: 0, head: 0, t: new Float32Array(CAP), M: new Float32Array(CAP), n: new Float32Array(CAP), rec: new Float32Array(CAP), dnn: new Float32Array(CAP), ion: new Float32Array(CAP) };
      this.trail = { M: [], n: [], rec: [], dnn: [], ion: [], mode: [] }; // logged operating points
      this.sweep = { active: false, points: null, idx: 0, phase: 0, timer: 0, save: 0 };
      this.sweepRuns = []; // completed sweeps, tagged with coupling mode (max 4)
      this.tReal = 0; this.tSim = 0;
      this._probeAcc = 0; this._trailAcc = 0;
      this._ionCount = 0; this._wallCount = 0; this._collCount = 0;
      this.pushRate = 0;
      this.bdotLen = 160;
      this.bdot = new Float32Array(8 * this.bdotLen);
      this.bdotHead = 0;
      this._cq = new Float32Array(4);
      this._nCal = 0; this._nSmooth = 0; this._chuckCount = 0;
      this.nRes = 0; this.fRes = 0;
      this.dust = window.EHRDust ? new window.EHRDust(this) : null;
      this.alloc();
      this.reset();
    }

    alloc() {
      const N = this.params.N;
      this.pos = new Float32Array(3 * N);
      this.vel = new Float32Array(3 * N);
    }

    setN(N) { this.params.N = N; this.alloc(); this.reset(); }

    reset() {
      const N = this.params.N, p = this.pos, v = this.vel;
      for (let i = 0; i < N; i++) this.spawn(i, p, v);
      this.grid.fill(0); this.gridMax = 0.001;
      this.hist.len = 0; this.hist.head = 0;
      this.tReal = 0; this.tSim = 0;
      this._ionCount = 0; this._wallCount = 0; this._collCount = 0;
      this.diag.ionRate = 0; this.diag.resTime = 0;
      this._nCal = 0; this._nSmooth = 0; this._chuckCount = 0;
      if (this._cq) this._cq.fill(0);
      if (this.bdot) this.bdot.fill(0);
      if (this.dust) this.dust.reseed();
    }

    spawn(i, p, v) {
      const r = Math.sqrt(Math.random()) * 0.75, th = Math.random() * TWO_PI;
      p[3 * i] = r * Math.cos(th); p[3 * i + 1] = r * Math.sin(th);
      p[3 * i + 2] = Math.random() * 4 - 2;
      const vth = 0.5; // ~3.75 eV thermal-ish
      v[3 * i] = gauss() * vth; v[3 * i + 1] = gauss() * vth; v[3 * i + 2] = gauss() * vth;
    }

    setParam(k, val) {
      this.params[k] = val;
      if (k === 'Bz' || k === 'btRatio' || k === 'kzA' || k === 'mMode' || k === 'alpha' || k === 'fieldModel') this.geomDirty = true;
    }

    get M() {
      const p = this.params;
      const ktheta = p.mMode / 1.0; // at wall radius R=1
      return (ktheta / p.kzA) * (1 / Math.max(this._bthOverBz(1), 1e-4));
    }
    // local Bθ/Bz ratio for either field model
    _bthOverBz(r) {
      const p = this.params;
      if (p.fieldModel === 'beltrami') {
        const j0 = this._besselJ(0, p.alpha * r);
        return this._besselJ(1, p.alpha * r) / Math.max(j0, 1e-3);
      }
      return p.btRatio * r;
    }
    pitchAt(r) { return TWO_PI * r / Math.max(this._bthOverBz(r), 1e-4); } // L_p = 2πr·Bz/Bθ
    _getFieldLUT() {
      const p = this.params;
      const key = p.fieldModel + '|' + p.Bz + '|' + p.btRatio + '|' + p.alpha;
      if (this._fKey === key) return this._fLUT;
      const bthOverR = new Float32Array(130), bzArr = new Float32Array(130);
      for (let i = 0; i <= 129; i++) {
        const r = i / 128;
        if (p.fieldModel === 'beltrami') {
          bthOverR[i] = p.Bz * (r < 1e-4 ? p.alpha / 2 : this._besselJ(1, p.alpha * r) / r);
          bzArr[i] = p.Bz * this._besselJ(0, p.alpha * r);
        } else {
          bthOverR[i] = p.Bz * p.btRatio;
          bzArr[i] = p.Bz;
        }
      }
      this._fKey = key; this._fLUT = { bthOverR, bz: bzArr };
      return this._fLUT;
    }
    // coupling: 'kernel' = hypothesis prior peaked at M=1 (+ non-resonant floor);
    // 'emergent' = volume-averaged (k̂·b̂)²·F from antenna/field geometry — no M prior.
    get coupling() {
      if (this.params.couplingMode === 'emergent') return this._emergentCoupling();
      if (this.params.couplingMode === 'wave') return this._waveCoupling() * this._dispFactor();
      const lm = Math.log(this.M);
      return 0.13 + 0.87 * Math.exp(-(lm * lm) / (2 * 0.42 * 0.42));
    }
    _emergentCoupling() {
      const p = this.params;
      const key = p.kzA + '|' + p.mMode + '|' + p.btRatio + '|' + p.fieldModel + '|' + p.alpha;
      if (this._cplKey === key) return this._cplVal;
      const LC2 = 16; // coherence length² (device scale)
      let s = 0, wsum = 0;
      for (let i = 0; i < 32; i++) {
        const r = (i + 0.5) / 32;
        const kth = p.mMode / Math.max(r, 0.04);
        const bth = this._bthOverBz(r); // Bθ/Bz
        const bnorm = Math.sqrt(bth * bth + 1);
        const knorm = Math.sqrt(kth * kth + p.kzA * p.kzA);
        const align = (kth * bth + p.kzA) / (knorm * bnorm);
        const kmis = (kth - p.kzA * bth) / bnorm; // k⊥ within flux surface (Bz-normalized)
        const F = 1 / (1 + LC2 * kmis * kmis);
        const wgt = 2 * r * (0.35 + 0.65 * Math.exp(-(1 - r) / 0.5));
        s += align * align * F * wgt; wsum += wgt;
      }
      this._cplKey = key; this._cplVal = s / wsum;
      return this._cplVal;
    }
    // WAVE mode: emergent coupling evaluated at the antenna (wall), deposition via Bessel eigenmode
    _waveCoupling() {
      const p = this.params;
      const key = 'w' + p.kzA + '|' + p.mMode + '|' + p.btRatio + '|' + p.fieldModel + '|' + p.alpha;
      if (this._cplKey === key) return this._cplVal;
      const bt = this._bthOverBz(1), m = p.mMode, kz = p.kzA;
      const bnorm = Math.sqrt(bt * bt + 1), knorm = Math.sqrt(m * m + kz * kz);
      const align = (m * bt + kz) / (knorm * bnorm); // k̂·b̂ at r=R
      const kmis = (m - kz * bt) / bnorm;            // k⊥ at wall ∝ (M−1)
      const F = 1 / (1 + 16 * kmis * kmis);
      this._cplKey = key; this._cplVal = 0.1 + 0.9 * align * align * F;
      return this._cplVal;
    }
    // Helicon dispersion (WAVE mode): ω ∝ B·k_z·k_tot/n_e ⇒ resonant density n_res(f,B).
    // Heating feedback locks n̄ onto the dispersion curve — f-scans move the locked density.
    _dispFactor() {
      const p = this.params;
      if (!this._nCal) { this.nRes = 0; this.fRes = 0; return 1; }
      const ktot = Math.sqrt(p.mMode * p.mMode + p.kzA * p.kzA);
      const K0 = 3.019; // Bz·kz·k_tot at defaults (1 · 1.6 · √(1+1.6²))
      const w = Math.max(p.rfFreq / 13.56, 0.05);
      this.nRes = this._nCal * (p.Bz * p.kzA * ktot) / (K0 * w);
      const nS = Math.max(this._nSmooth, 1e-3);
      this.fRes = 13.56 * (p.Bz * p.kzA * ktot) * this._nCal / (K0 * nS);
      const x = Math.log(nS / this.nRes) / 0.45;
      return 0.2 + 0.8 / (1 + x * x);
    }
    _besselJ(n, x) { // series, fine for x ≤ 8
      let term = 1;
      for (let k = 1; k <= n; k++) term *= (x / 2) / k;
      let sum = term;
      const x2 = -(x * x) / 4;
      for (let k = 1; k < 30; k++) {
        term *= x2 / (k * (k + n));
        sum += term;
        if (Math.abs(term) < 1e-12) break;
      }
      return sum;
    }
    _getWaveLUT() { // deposition ∝ J_{m-1}(j_{m,1}·r)², volume-normalized
      const m = this.params.mMode;
      if (this._waveLUTm === m && this._waveLUT) return this._waveLUT;
      const jz = { 1: 3.8317, 2: 5.1356, 3: 6.3802 }[m] || 3.8317;
      const lut = new Float32Array(64);
      let mean = 0;
      for (let i = 0; i < 64; i++) {
        const r = (i + 0.5) / 64;
        const v = this._besselJ(m - 1, jz * r);
        lut[i] = v * v;
        mean += lut[i] * 2 * r / 64;
      }
      const inv = 1 / Math.max(mean, 1e-6);
      for (let i = 0; i < 64; i++) lut[i] *= inv;
      this._waveLUTm = m; this._waveLUT = lut;
      return lut;
    }
    get resonant() { // 0..1 resonant fraction (controls edge→core deposition shift)
      const lm = Math.log(this.M);
      return Math.exp(-(lm * lm) / (2 * 0.42 * 0.42));
    }

    fieldAt(x, y, z, out) {
      const lut = this._getFieldLUT();
      const r = Math.sqrt(x * x + y * y);
      const i = Math.min(129, (r * 128) | 0);
      out[0] = -y * lut.bthOverR[i];
      out[1] = x * lut.bthOverR[i];
      out[2] = lut.bz[i];
      return out;
    }

    nn(z) { // neutral density incl. acoustic modulation (uses current tReal)
      const p = this.params;
      const n0 = p.pressure * 0.5;
      const k = TWO_PI / p.acLambda;
      const ph = TWO_PI * p.acFreq * this.tReal;
      const fac = p.acMode === 'standing'
        ? 1 + p.acA * Math.sin(k * z) * Math.cos(ph)
        : 1 + p.acA * Math.sin(k * z - ph);
      return n0 * Math.max(0.05, fac);
    }

    step(dtMs) {
      const dtReal = clamp(dtMs, 0, 50) / 1000;
      this.tReal += dtReal;
      if (this.sweep.active) this._sweepUpdate(dtReal);

      const p = this.params, N = p.N, pos = this.pos, vel = this.vel, grid = this.grid;
      const qmh = -0.5 * DT; // (q/m)·dt/2 with q/m=-1
      const flut = this._getFieldLUT();
      const fBthR = flut.bthOverR, fBzArr = flut.bz;
      const m = p.mMode, kz = p.kzA;
      const wRF = (p.rfFreq / 13.56) * 2.4; // visualization-scale wave rotation
      const sqrtP = Math.sqrt(p.rfPower / 40);
      const cpl = this.coupling, res = this.resonant;
      const kickBase = KICK * sqrtP * cpl * DT;
      const emergent = p.couplingMode === 'emergent';
      const wave = p.couplingMode === 'wave';
      const waveLUT = wave ? this._getWaveLUT() : null;
      const kickWave = KICK * 1.15 * sqrtP * cpl * DT;
      const kickEm = KICK * 1.9 * sqrtP * DT;
      const LC2 = 16;
      const nn0 = p.pressure * 0.5;
      const kac = TWO_PI / p.acLambda;
      const acPhaseT = TWO_PI * p.acFreq * this.tReal;
      const acA = p.acA;
      const standing = p.acMode === 'standing';
      const cosAc = Math.cos(acPhaseT);
      let ion = 0, wall = 0, coll = 0, chuck = 0;

      for (let s = 0; s < SUBSTEPS; s++) {
        this.tSim += DT;
        const wPhaseT = wRF * this.tReal;
        for (let i = 0; i < N; i++) {
          const i3 = 3 * i;
          let x = pos[i3], y = pos[i3 + 1], z = pos[i3 + 2];
          let vx = vel[i3], vy = vel[i3 + 1], vz = vel[i3 + 2];

          const r2 = x * x + y * y;
          const r = Math.sqrt(r2);
          // --- Boris rotation (E=0) ---
          const fi = r >= 1.0078 ? 129 : (r * 128) | 0;
          const bor = fBthR[fi];
          const Bx = -bor * y, By = bor * x, Bzv = fBzArr[fi];
          const tx = qmh * Bx, ty = qmh * By, tz = qmh * Bzv;
          const t2 = tx * tx + ty * ty + tz * tz;
          const sf = 2 / (1 + t2);
          const sx = tx * sf, sy = ty * sf, sz = tz * sf;
          let ux = vx + (vy * tz - vz * ty);
          let uy = vy + (vz * tx - vx * tz);
          let uz = vz + (vx * ty - vy * tx);
          vx += uy * sz - uz * sy;
          vy += uz * sx - ux * sz;
          vz += ux * sy - uy * sx;

          // --- RF heating kick along B̂, helical antenna pattern ---
          const phase = m * Math.atan2(y, x) - kz * z - wPhaseT;
          const w = Math.cos(phase);
          const bmag = Math.sqrt(Bx * Bx + By * By + Bzv * Bzv) + 1e-9;
          let kick;
          if (wave) {
            const zEnv = 0.55 + 0.45 * Math.exp(-z * z / 1.96);
            kick = kickWave * w * waveLUT[Math.min(63, (r * 64) | 0)] * zEnv;
          } else if (emergent) {
            // E∥ projection (k̂·b̂) × transit coherence F(k⊥) — M-dependence emerges from geometry
            const kth = m / Math.max(r, 0.04);
            const knorm = Math.sqrt(kth * kth + kz * kz);
            const align = (kth * bor * r + kz * Bzv) / (knorm * bmag);
            const kmis = (kth * Bzv - kz * bor * r) / bmag;
            const F = 1 / (1 + LC2 * kmis * kmis);
            const env = 0.35 + 0.65 * Math.exp(-(1 - Math.min(r, 1)) / 0.5);
            kick = kickEm * w * align * align * F * env;
          } else {
            const edgeProf = Math.exp(-(1 - Math.min(r, 1)) / 0.22);
            const coreProf = Math.exp(-r2 / 0.32);
            kick = kickBase * w * (edgeProf * (1 - res) + coreProf * res);
          }
          vx += kick * Bx / bmag; vy += kick * By / bmag; vz += kick * Bzv / bmag;

          // --- energy cap (proxy for radiative/inelastic losses) ---
          let E = EVK * (vx * vx + vy * vy + vz * vz);
          if (E > 75) { const f = Math.sqrt(75 / E); vx *= f; vy *= f; vz *= f; E = 75; }

          // --- Monte Carlo collision vs neutrals ---
          const acFac = standing ? 1 + acA * Math.sin(kac * z) * cosAc : 1 + acA * Math.sin(kac * z - acPhaseT);
          const nnLoc = nn0 * (acFac < 0.05 ? 0.05 : acFac);
          const vmag = Math.sqrt(vx * vx + vy * vy + vz * vz) + 1e-9;
          const nuEl = nnLoc * SIGMA_EL * vmag;
          const nuIon = E > EION ? nnLoc * SIGMA_ION * vmag * (1 - EION / E) : 0;
          const Pc = (nuEl + nuIon) * DT;
          if (Math.random() < Pc) {
            coll++;
            if (Math.random() < nuIon / (nuEl + nuIon)) {
              // ionization: pay threshold, deposit ion in voxel
              E = Math.max(0.4, E - EION);
              const vnew = Math.sqrt(E / EVK);
              const cth = Math.random() * 2 - 1, sth = Math.sqrt(1 - cth * cth), ph2 = Math.random() * TWO_PI;
              vx = vnew * sth * Math.cos(ph2); vy = vnew * sth * Math.sin(ph2); vz = vnew * cth;
              const ix = clamp((x + 1) * 0.5 * NX | 0, 0, NX - 1);
              const iy = clamp((y + 1) * 0.5 * NY | 0, 0, NY - 1);
              const iz = clamp((z + 2) * 0.25 * NZ | 0, 0, NZ - 1);
              grid[(iz * NY + iy) * NX + ix] += 1;
              ion++;
            } else {
              // elastic: isotropic scatter, small energy transfer to neutral
              const sp = vmag * 0.99;
              const cth = Math.random() * 2 - 1, sth = Math.sqrt(1 - cth * cth), ph2 = Math.random() * TWO_PI;
              vx = sp * sth * Math.cos(ph2); vy = sp * sth * Math.sin(ph2); vz = sp * cth;
            }
          }

          // --- advance ---
          x += vx * DT; y += vy * DT; z += vz * DT;
          if (z > 2) {
            if (x * x + y * y < 0.6084) { // r<0.78 — wafer chuck absorbs
              chuck++;
              const rr = Math.sqrt(Math.random()) * 0.75, th = Math.random() * TWO_PI;
              x = rr * Math.cos(th); y = rr * Math.sin(th); z = -1.9 + Math.random() * 1.6;
              const vth = 0.5;
              vx = gauss() * vth; vy = gauss() * vth; vz = gauss() * vth;
            } else { z = 4 - z; vz = -vz; }
          } else if (z < -2) { z = -4 - z; vz = -vz; }
          if (x * x + y * y > 1) {
            // wall absorption + cold recycling
            wall++;
            const rr = Math.sqrt(Math.random()) * 0.75, th = Math.random() * TWO_PI;
            x = rr * Math.cos(th); y = rr * Math.sin(th); z = Math.random() * 4 - 2;
            const vth = 0.5;
            vx = gauss() * vth; vy = gauss() * vth; vz = gauss() * vth;
          }
          pos[i3] = x; pos[i3 + 1] = y; pos[i3 + 2] = z;
          vel[i3] = vx; vel[i3 + 1] = vy; vel[i3 + 2] = vz;
        }
      }
      this._ionCount += ion; this._wallCount += wall; this._collCount += coll; this._chuckCount += chuck;
      this.pushRate = this.pushRate * 0.9 + 0.1 * (N * SUBSTEPS / Math.max(dtReal, 1e-3) / 1e6);

      // density loss (volumetric + wall) → steady state maps ionization source
      const dec = Math.exp(-dtReal / DECAY_TAU);
      for (let i = 0; i < grid.length; i++) grid[i] *= dec;

      // dust dynamics ride on the live ion grid
      if (this.dust) this.dust.update(dtReal);

      // probes at 10 Hz
      this._probeAcc += dtReal;
      if (this._probeAcc >= 0.1) { this._probePass(this._probeAcc); this._probeAcc = 0; }
      // operating-point trail at ~1.5 s cadence (skip during sweep; sweep logs its own)
      this._trailAcc += dtReal;
      if (this._trailAcc >= 1.5 && !this.sweep.active) {
        this._trailAcc = 0;
        if (this.tReal > 3) this._logTrail();
      }
    }

    _logTrail() {
      const t = this.trail, d = this.diag;
      t.M.push(this.M); t.n.push(d.nMean); t.rec.push(d.Rec); t.dnn.push(d.dnn); t.ion.push(d.ionRate); t.mode.push(this.params.couplingMode);
      if (t.M.length > 400) { t.M.shift(); t.n.shift(); t.rec.shift(); t.dnn.shift(); t.ion.shift(); t.mode.shift(); }
    }
    clearTrail() { const t = this.trail; t.M.length = t.n.length = t.rec.length = t.dnn.length = t.ion.length = t.mode.length = 0; if (this.sweep.points) this.sweep.points.length = 0; this.sweepRuns.length = 0; }

    _probePass(dt) {
      const g = this.grid, rT = this.rTab;
      const radial = this.radial, radCnt = this._radCnt, axial = this.axial, gpi = this.gpi;
      radial.fill(0); radCnt.fill(0); axial.fill(0); gpi.fill(0);
      let sC = 0, cC = 0, sE = 0, cE = 0, tot = 0, cnt = 0, max = 0;
      for (let iz = 0; iz < NZ; iz++) {
        const zOff = iz * NY * NX, azBin = iz >> 1;
        for (let iy = 0; iy < NY; iy++) {
          const yOff = zOff + iy * NX, rOff = iy * NX;
          for (let ix = 0; ix < NX; ix++) {
            const r = rT[rOff + ix];
            if (r > 1) continue;
            const d = g[yOff + ix];
            tot += d; cnt++;
            if (d > max) max = d;
            gpi[rOff + ix] += d;
            axial[azBin] += d;
            const rb = Math.min(27, (r * 28) | 0);
            radial[rb] += d; radCnt[rb]++;
            if (r < 0.3) { sC += d; cC++; } else if (r > 0.8) { sE += d; cE++; }
          }
        }
      }
      for (let b = 0; b < 28; b++) radial[b] = radCnt[b] ? radial[b] / radCnt[b] : 0;
      for (let b = 0; b < 64; b++) axial[b] /= (this.maskCount * 2);
      const d = this.diag;
      d.nCore = cC ? sC / cC : 0;
      d.nEdge = cE ? sE / cE : 0;
      d.Rec = d.nCore > 1e-6 ? d.nEdge / d.nCore : 0;
      d.nMean = cnt ? tot / cnt : 0;
      this.gridMax = Math.max(0.001, this.gridMax * 0.7 + max * 0.3);
      d.ionRate = d.ionRate * 0.7 + 0.3 * (this._ionCount / Math.max(dt, 1e-3) / 1000); // k events/s
      d.resTime = this._wallCount > 0 ? this.params.N / (this._wallCount / Math.max(dt, 1e-3)) : 0;
      d.refl = (1 - this.coupling) * this.params.rfPower;
      this._ionCount = 0; this._wallCount = 0; this._collCount = 0;
      // mean energy (sampled)
      let se = 0; const S = Math.min(2000, this.params.N), v = this.vel, stride = Math.max(1, (this.params.N / S) | 0);
      let ns = 0;
      for (let i = 0; i < this.params.N; i += stride) { const i3 = 3 * i; se += v[i3] * v[i3] + v[i3 + 1] * v[i3 + 1] + v[i3 + 2] * v[i3 + 2]; ns++; }
      d.meanE = EVK * se / Math.max(ns, 1);
      // ---- wafer: uniformity across chuck footprint + ion flux/energy ----
      const pp = this.params;
      let wS = 0, wS2 = 0, wC = 0;
      for (let iy = 0; iy < NY; iy++) for (let ix = 0; ix < NX; ix++) {
        if (rT[iy * NX + ix] > 0.7) continue;
        let col = 0;
        for (let iz = 120; iz < 128; iz++) col += g[(iz * NY + iy) * NX + ix];
        wS += col; wS2 += col * col; wC++;
      }
      const wMu = wC ? wS / wC : 0;
      d.waferUnif = wMu > 1e-6 ? Math.sqrt(Math.max(wS2 / wC - wMu * wMu, 0)) / wMu : 0;
      d.ionFlux = d.ionFlux * 0.7 + 0.3 * (this._chuckCount / Math.max(dt, 1e-3) / 1000);
      this._chuckCount = 0;
      d.ionEw = pp.bias + 0.5 * d.meanE;
      // ---- dispersion calibration + smoothed density ----
      this._nSmooth = this._nSmooth * 0.8 + d.nMean * 0.2;
      if (!this._nCal && this.tReal > 5) this._nCal = Math.max(d.nMean, 1e-3);
      // ---- Prediction 3: correlation of n_e(z) with acoustic drive n_n(z) ----
      const nnA = this.axialNeutral();
      let axm = 0, nnm = 0;
      for (let b = 0; b < 64; b++) { axm += axial[b]; nnm += nnA[b]; }
      axm /= 64; nnm /= 64;
      let c01 = 0, c00 = 0, c11 = 0;
      for (let b = 0; b < 64; b++) {
        const a0 = axial[b] - axm, a1 = nnA[b] - nnm;
        c01 += a0 * a1; c00 += a0 * a0; c11 += a1 * a1;
      }
      d.acCorr = (c00 > 1e-12 && c11 > 1e-12) ? c01 / Math.sqrt(c00 * c11) : 0;
      // ---- Prediction 1: synthetic B-dot probe array (8 probes @ r=0.75R, z=0) ----
      const wRFp = (pp.rfFreq / 13.56) * 2.4;
      const ampB = Math.sqrt(Math.max(this.coupling, 0)) * Math.min(this._nCal ? d.nMean / this._nCal : 0.5, 1.5);
      const sigN = 0.22 / (0.12 + ampB * ampB * 2);
      const hb = this.bdotHead;
      let re1 = 0, im1 = 0, re2 = 0, im2 = 0, re3 = 0, im3 = 0;
      for (let i = 0; i < 8; i++) {
        const th = i / 8 * TWO_PI;
        const s = ampB * Math.cos(pp.mMode * th - wRFp * this.tReal) + gauss() * sigN;
        this.bdot[hb * 8 + i] = s;
        re1 += s * Math.cos(th); im1 += s * Math.sin(th);
        re2 += s * Math.cos(2 * th); im2 += s * Math.sin(2 * th);
        re3 += s * Math.cos(3 * th); im3 += s * Math.sin(3 * th);
      }
      this.bdotHead = (hb + 1) % this.bdotLen;
      const cq = this._cq;
      cq[1] = cq[1] * 0.9 + 0.1 * (re1 * re1 + im1 * im1);
      cq[2] = cq[2] * 0.9 + 0.1 * (re2 * re2 + im2 * im2);
      cq[3] = cq[3] * 0.9 + 0.1 * (re3 * re3 + im3 * im3);
      let best = 1;
      if (cq[2] > cq[best]) best = 2;
      if (cq[3] > cq[best]) best = 3;
      d.mFit = best;
      d.mSnr = cq[best] / ((cq[1] + cq[2] + cq[3] - cq[best]) / 2 + 1e-9);
      // history push
      const h = this.hist, idx = h.head;
      h.t[idx] = this.tReal; h.M[idx] = this.M; h.n[idx] = d.nMean; h.rec[idx] = d.Rec; h.ion[idx] = d.ionRate;
      // δn/n over last ~100 samples of nMean
      const K = Math.min(100, h.len + 1);
      let mu = 0;
      for (let k = 0; k < K; k++) mu += h.n[(idx - k + h.cap) % h.cap];
      mu /= K;
      let va = 0;
      for (let k = 0; k < K; k++) { const dv = h.n[(idx - k + h.cap) % h.cap] - mu; va += dv * dv; }
      d.dnn = mu > 1e-6 ? Math.sqrt(va / K) / mu : 0;
      h.dnn[idx] = d.dnn;
      h.head = (h.head + 1) % h.cap;
      if (h.len < h.cap) h.len++;
    }

    // Recent history as plain arrays (newest last), n samples
    histSlice(key, n) {
      const h = this.hist, len = Math.min(n, h.len);
      const out = new Float32Array(len), tArr = new Float32Array(len);
      for (let k = 0; k < len; k++) {
        const idx = (h.head - len + k + h.cap) % h.cap;
        out[k] = h[key][idx]; tArr[k] = h.t[idx];
      }
      return { t: tArr, v: out };
    }

    computePSD(nBins) {
      const { v } = this.histSlice('n', 256);
      const n = v.length;
      if (n < 32) return null;
      let mu = 0; for (let i = 0; i < n; i++) mu += v[i]; mu /= n;
      const bins = Math.min(nBins || 96, n >> 1);
      const psd = new Float32Array(bins), freq = new Float32Array(bins);
      const fs = 10; // Hz sampling
      for (let k = 1; k <= bins; k++) {
        let re = 0, im = 0;
        const w = TWO_PI * k / n;
        for (let i = 0; i < n; i++) {
          const s = (v[i] - mu) * (0.5 - 0.5 * Math.cos(TWO_PI * i / (n - 1))); // Hann
          re += s * Math.cos(w * i); im -= s * Math.sin(w * i);
        }
        psd[k - 1] = (re * re + im * im) / n;
        freq[k - 1] = k * fs / n;
      }
      return { freq, psd };
    }

    axialNeutral() { // analytic nn(z)/nn0 overlay, 64 pts
      const out = new Float32Array(64), p = this.params;
      const k = TWO_PI / p.acLambda, ph = TWO_PI * p.acFreq * this.tReal;
      for (let b = 0; b < 64; b++) {
        const z = (b + 0.5) / 64 * 4 - 2;
        const fac = p.acMode === 'standing'
          ? 1 + p.acA * Math.sin(k * z) * Math.cos(ph)
          : 1 + p.acA * Math.sin(k * z - ph);
        out[b] = Math.max(0.05, fac);
      }
      return out;
    }

    // ---- field lines: integrate dx/ds = B̂ ----
    fieldLines(nLines) {
      const lines = [], p = this.params;
      const radii = [0.28, 0.52, 0.76, 0.96];
      const perRing = Math.max(4, Math.round((nLines || 72) / radii.length));
      const B = [0, 0, 0];
      for (const r0 of radii) {
        for (let j = 0; j < perRing; j++) {
          const th0 = j / perRing * TWO_PI;
          let x = r0 * Math.cos(th0), y = r0 * Math.sin(th0), z = -2;
          const pts = [], mags = [];
          const h = 0.055;
          for (let s = 0; s < 900; s++) {
            this.fieldAt(x, y, z, B);
            const bm = Math.sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]) + 1e-9;
            pts.push(x, y, z); mags.push(bm);
            x += B[0] / bm * h; y += B[1] / bm * h; z += B[2] / bm * h;
            if (z > 2 || x * x + y * y > 1.1) break;
          }
          if (pts.length > 6) lines.push({ pts: new Float32Array(pts), mag: new Float32Array(mags) });
        }
      }
      return lines;
    }

    // ---- parameter sweeps: M (via field twist) or f (RF frequency) ----
    startSweep(variable) {
      if (this.sweep.active) return;
      const v = variable === 'freq' ? 'freq' : 'M';
      const pts = [];
      const n = 21;
      for (let i = 0; i < n; i++) {
        const t = i / (n - 1);
        const x = v === 'freq' ? 2 + t * 25 : Math.exp(Math.log(0.25) + (Math.log(4) - Math.log(0.25)) * t);
        pts.push({ x, n: 0, rec: 0, dnn: 0, ion: 0, samples: 0 });
      }
      this.sweep = { active: true, variable: v, points: pts, idx: 0, phase: 0, timer: 0,
        save: { bt: this.params.btRatio, al: this.params.alpha, f: this.params.rfFreq }, done: false };
      this._applySweepX(pts[0].x);
    }
    cancelSweep() {
      if (!this.sweep.active) return;
      this._restoreSweep();
      this.sweep.active = false;
    }
    _restoreSweep() {
      const s = this.sweep.save;
      this.params.btRatio = s.bt; this.params.alpha = s.al; this.params.rfFreq = s.f;
      this.geomDirty = true;
    }
    _applySweepX(x) {
      const p = this.params;
      if (this.sweep.variable === 'freq') { p.rfFreq = x; return; }
      if (p.fieldModel === 'beltrami') {
        // M(α) = (m/kz)·J0(α)/J1(α) is monotone decreasing — bisect
        let lo = 0.05, hi = 2.3;
        for (let i = 0; i < 40; i++) {
          p.alpha = (lo + hi) / 2;
          if (this.M > x) lo = p.alpha; else hi = p.alpha;
        }
        p.alpha = (lo + hi) / 2;
      } else {
        // M = (m/kz)·(1/btRatio) ⇒ btRatio = m/(kz·M)
        p.btRatio = clamp(p.mMode / (p.kzA * x), 0.02, 4);
      }
      this.geomDirty = true;
    }
    _sweepUpdate(dt) {
      const sw = this.sweep;
      sw.timer += dt;
      const SETTLE = sw.variable === 'freq' ? 1.6 : 1.0, MEASURE = 1.0;
      const pt = sw.points[sw.idx];
      if (sw.phase === 0) {
        if (sw.timer >= SETTLE) { sw.phase = 1; sw.timer = 0; }
      } else {
        const d = this.diag;
        pt.n += d.nMean; pt.rec += d.Rec; pt.dnn += d.dnn; pt.ion += d.ionRate; pt.samples++;
        if (sw.timer >= MEASURE) {
          pt.n /= pt.samples; pt.rec /= pt.samples; pt.dnn /= pt.samples; pt.ion /= pt.samples;
          sw.idx++;
          if (sw.idx >= sw.points.length) {
            this.sweepRuns.push({
              variable: sw.variable, mode: this.params.couplingMode, field: this.params.fieldModel,
              points: sw.points.map(p => ({ x: p.x, n: p.n, rec: p.rec, dnn: p.dnn, ion: p.ion })),
            });
            if (this.sweepRuns.length > 4) this.sweepRuns.shift();
            this._restoreSweep();
            sw.active = false; sw.done = true;
          } else {
            this._applySweepX(sw.points[sw.idx].x);
            sw.phase = 0; sw.timer = 0;
          }
        }
      }
    }
    get sweepProgress() {
      const sw = this.sweep;
      if (!sw.points || !sw.points.length) return 0;
      return sw.active ? sw.idx / sw.points.length : (sw.done ? 1 : 0);
    }
  }

  EHREngine.GRID = { NX, NY, NZ };
  window.EHREngine = EHREngine;
})();
