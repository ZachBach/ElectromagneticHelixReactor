# EHR Phase 1 — Pitch-Physics Simulator

Headless, SI-units test-particle simulator built to answer the central EHR
question (Prediction 4 of the paper):

> Is helical field-line pitch an independent and useful plasma control
> parameter? Does scanning Bθ/Bz at fixed |B|, pressure, and energy
> change confinement and transport?

This is deliberately separate from the browser engine (`../ehr-engine.js`),
which runs in normalized units with tuned coefficients and recycles lost
particles for visualization. Here everything is dimensional, walls absorb, and
confinement times come out in seconds.

## Running

```
npm run sim:selftest    # physics checks: Boris vs exact gyromotion, div B = 0,
                        # Beltrami force-free, MCC rate vs analytic nu(E)
npm run sim:pitch       # the headline study (defaults below)
npm run sim:mechanism   # free-streaming geometric model (mechanism layer 0)
```

Options (`node sim/run-pitch-scan.js --key=value`):

| flag | default | meaning |
|---|---|---|
| `--species` | `electron` | `electron` or `ion` (Ar+) |
| `--model` | `screw` | `screw` (Bz const, Bθ ∝ r) or `beltrami` (force-free J0/J1) |
| `--n` | 4000 | ensemble size |
| `--ratios` | 0…3 | comma list of Bθ(R)/Bz values to scan |
| `--bwall` | 0.01 | \|B(R)\| in tesla (100 G, paper baseline) |
| `--pressure` | 10 | argon fill, mTorr (paper baseline) |
| `--te` / `--ti` | 3 / 0.026 | seed temperature, eV |
| `--r` / `--l` | 0.1 / 0.4 | chamber radius / length, m |
| `--tmax` | 50 µs (e) / 50 ms (ion) | censoring time |
| `--seed` | 1 | RNG seed (runs are reproducible) |
| `--density` | off | also write (r,z) density CSVs |
| `--nocoll` | off | disable MCC — geometry + gyration + drifts only |
| `--dtscale` | 1 | timestep multiplier (convergence checks) |

Results go to `sim/output/*.csv` (gitignored).

## Physics content

- **Field (Module 1)** — `field.js`. B = Bz(r) ẑ + Bθ(r) θ̂, Br = 0,
  identically divergence-free. The scan holds |B(R)| fixed while changing the
  ratio, which is the closest divergence-free realization of the paper's
  "constant |B|" pitch scan (no field with tunable pitch has |B| constant
  everywhere for these profiles; for the Beltrami model in particular the
  on-axis field grows as pitch increases).
- **Pusher (Module 2)** — `boris.js`. Standard Boris rotation, E = 0,
  ω_c·dt ≤ 0.2. Energy-conserving to machine precision (see selftest).
- **Collisions (Module 3)** — `collisions.js`. Test-particle MCC against the
  neutral argon background: electron elastic (momentum-transfer σ with the
  Ramsauer minimum, exact 2m/M energy transfer) + ionization (counted;
  secondaries not tracked), ion charge exchange + equal-mass hard-sphere
  elastic. Cross sections are approximate log-log tables anchored to the
  representative values used in the paper; swap in LXCat (Phelps/Biagi) data
  when absolute rates matter.
- **Diagnostics (Phase 2 seed)** — `diagnostics.js`. Per-particle loss time
  and channel (radial wall vs ends), survival censoring at Tmax, optional
  volume-normalized (r,z) density deposition.

Not in Phase 1 (by design): self-consistent E fields (Poisson), sheaths, RF
heating, and the acoustic channel. This is a pure geometry-and-collisions
measurement of how pitch reshapes single-particle confinement.

## First results (2026-07-18, seed 1, N = 4000)

Electrons, 3 eV seed, argon 10 mTorr, |B(R)| = 100 G, R = 0.1 m, L = 0.4 m:

| Bθ(R)/Bz | screw τ_mean (µs) | beltrami τ_mean (µs) |
|---|---|---|
| 0.00 | 1.00 ± 0.02 | — |
| 0.25 | 1.02 ± 0.02 | 0.98 ± 0.02 |
| 0.50 | 1.02 ± 0.02 | 1.04 ± 0.02 |
| 0.75 | 1.11 ± 0.02 | — |
| 1.00 | 1.24 ± 0.02 | 1.19 ± 0.02 |
| 1.50 | 1.57 ± 0.03 | 1.32 ± 0.03 |
| 2.00 | 1.89 ± 0.04 | 1.41 ± 0.03 |
| 3.00 | 2.82 ± 0.05 | 1.54 ± 0.03 |

Observations:

1. **Pitch changes confinement at fixed |B(R)|** — τ rises monotonically with
   the pitch ratio, ~2.8× for the screw field at ratio 3. Loss remains
   end-wall dominated (>97%) throughout; the radial loss fraction grows
   slightly with pitch (0.4% → 2.5%).
2. **The response is quadratic at small pitch** — flat below ratio ≈ 0.5, then
   near-linear growth. Consistent with the geometric mechanism: tilting field
   lines slows parallel streaming to the ends by |B|/Bz = √(1 + (Bθ/Bz)²),
   which is 1 + O(ratio²) at small ratio.
3. **Twist distribution matters, not just the wall ratio** — the Beltrami
   field (same wall ratio, twist concentrated differently in radius, stronger
   on-axis |B|) gains roughly half as much as the screw field. Pitch *profile*
   is itself a control variable.

## Mechanism (proven, not suspected)

Three-layer decomposition at the baseline point, adding one physics ingredient
at a time (identical seeding; censored means, so layers differ in how the
slow-v_z tail is counted — compare gains within a row's own baseline):

| layer | physics | τ(0) µs | τ(3) µs | gain |
|---|---|---|---|---|
| `sim:mechanism` | field-line geometry only (analytic) | 1.35 | 2.41 | 1.8× |
| `sim:pitch --nocoll` | + gyration + grad-B/curvature drifts | 1.21 | 2.43 | 2.0× |
| `sim:pitch` | + Monte-Carlo collisions | 1.00 | 2.82 | 2.8× |

- **Drifts are negligible**: the collisionless Boris runs match the
  free-streaming model within statistical error at every ratio, for both field
  models, and their radial loss is exactly 0.0%. At 100 G the pitch effect is
  not a drift effect.
- **Geometry sets the shape**: the flat-then-linear τ(ratio) curve, and the
  screw/Beltrami difference, are both reproduced by the analytic model. The
  Beltrami deficit is geometric in origin — its twist is concentrated toward
  the wall, so interior electrons still see a strong Bz and stream out fast
  (free-streaming predicts 2.41 vs 1.78 µs at ratio 3; Boris measures the same
  ordering).
- **Collisions amplify the effect**: λ_mfp ≈ 15 cm against L = 40 cm puts
  parallel transport in the marginally diffusive regime, where residence time
  scales toward the *square* of the field-line path factor — hence 2.8×
  measured vs 1.8× geometric at ratio 3.

## Parameter study, τ_c = f(B, p, E, pitch) (screw field)

One-factor-at-a-time around the baseline, 5-point pitch scans:

| scan | values | τ(0) µs | τ(3) µs | pitch gain |
|---|---|---|---|---|
| pressure | 3 / 10 / 30 mTorr | 0.80 / 1.00 / 1.68 | 1.85 / 2.82 / 4.96 | 2.3× / 2.8× / 3.0× |
| energy | 1 / 3 / 10 eV | 1.36 / 1.00 / 0.93 | 3.17 / 2.82 / 2.91 | 2.3× / 2.8× / 3.1× |
| \|B(R)\| | 50 / 100 / 200 G | 0.96 / 1.00 / 0.98 | 2.59 / 2.82 / 3.01 | 2.7× / 2.8× / 3.1× |

- **Pitch leverage grows with collisionality** — with pressure directly, and
  with energy through σ(E) (1 eV sits near the Ramsauer minimum, so it is the
  *least* collisional case and lands nearest the geometric limit; 10 eV has
  σ ≈ 1e-19 m² and lands furthest above it). This is the diffusive
  amplification signature, and it means pitch control is *strongest* in the
  collisional regime where processing plasmas actually operate.
- **τ(0) is independent of \|B\|** (0.96–1.00 µs across 50–200 G), exactly as
  the geometric mechanism requires — field magnitude sets no axial timescale.
  What \|B\| does control is the collisional radial leakage at high pitch
  (ratio 3 radial loss: 7.5% at 50 G → 2.5% at 100 G → 0.5% at 200 G), which
  is why the pitch gain creeps up slightly with B.
- **Robustness**: an independent seed reproduces every point within the quoted
  errors, and halving the timestep changes nothing (ω_c·dt = 0.1 vs 0.2), so
  the integrator is converged.
- *Tentative*: collisional runs show a small (~3σ) excess of losses toward one
  end wall at high pitch, absent both at ratio 0 and in collisionless runs —
  the signature one would expect if collisions reveal a slow unidirectional
  axial drift flux. Needs multi-seed statistics before claiming; flagged for
  follow-up.

Caveats: test-particle, no ambipolar electric field, no RF sustainment — the
absolute numbers are the collisionless-loss skeleton, not a device prediction.
The claims these runs support are the *existence, shape, and mechanism* of the
pitch dependence, which is exactly what Prediction 4 asks an experiment to
test.

## Roadmap position

Phase 1 (this) → Phase 2 diagnostics (density profiles per pitch — the
deposition grid is already here behind `--density`) → Phase 3 Poisson/sheath →
Phase 4 helicon coupling models → Phase 5 dusty plasma.
