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
npm run sim:alpha -- sim/output/<scan>.csv   # fit the transport-law exponent
npm run sim:ambipolar   # Phase 3: self-consistent E-field feedback study
```

Options (`node sim/run-pitch-scan.js --key=value`):

| flag | default | meaning |
|---|---|---|
| `--species` | `electron` | `electron` or `ion` (Ar+) |
| `--model` | `screw` | `screw` (Bz const, Bθ ∝ r), `beltrami` (force-free J0/J1), or `powerlaw` (Bθ ∝ r^n) |
| `--nexp` | 1 | powerlaw twist exponent n |
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

## Study A — the confinement exponent

Define the seeding-averaged field-line path factor
⟨F⟩ = ⟨|B|/Bz⟩ (computed numerically from the field model by `fit-alpha.js`)
and fit τ ∝ ⟨F⟩^α. Median-based fits (robust to the censored slow-v_z tail):

| run | α (median fit) |
|---|---|
| free-streaming, all four profiles | 0.95 – 0.99 |
| Boris, collisions off (screw) | 0.96 ± 0.02 |
| full, 3 mTorr | 1.19 ± 0.03 |
| full, 10 mTorr | 1.33 ± 0.03 |
| full, 30 mTorr | 1.39 ± 0.05 |
| full, 10 mTorr, powerlaw n = 0.5 / 2 / 4 | 1.41 / 1.16 / 1.15 (± 0.04–0.07) |

The collisionless exponent is **α = 1**, measured to a few percent across
every profile — residence time proportional to path length, as free streaming
requires. Collisions push α up toward the diffusive limit of 2, monotonically
with L/λ_mfp (λ ≈ 45 / 15 / 5 cm at 3 / 10 / 30 mTorr against L = 40 cm).
Mean-based fits give the same picture shifted slightly (collisionless means
are biased low by censoring). At fixed pressure the exponent is consistent
across profiles at the ±0.15 level — the residual spread is the current
precision limit of the law below.

## Study B — the pitch-profile map

`powerlaw` fields Bθ ∝ (r/R)^n at fixed wall ratio and fixed |B(R)| separate
"how much twist" from "where the twist lives". Full-physics τ_mean (µs),
10 mTorr, τ(0) = 1.00 for every profile:

| Bθ(R)/Bz | n=0.5 | n=1 (screw) | n=2 | n=4 | beltrami |
|---|---|---|---|---|---|
| 0.5 | 1.10 | 1.02 | 1.03 | 0.99 | 1.04 |
| 1.0 | 1.32 | 1.24 | 1.14 | 1.05 | 1.19 |
| 2.0 | 2.26 | 1.89 | 1.49 | 1.18 | 1.41 |
| 3.0 | **3.66** | 2.82 | 1.96 | **1.38** | 1.54 |

At the *same wall pitch ratio*, moving twist inward (n = 0.5) versus outward
(n = 4) changes confinement by **2.7×**. The profile is a control dimension as
strong as the pitch magnitude itself. The free-streaming model reproduces the
full ordering (2.92 / 2.41 / 2.07 / 1.61 µs at ratio 3), so this too is
geometric: interior twist slows the interior electrons that dominate the
ensemble; wall-concentrated twist (n = 4, and Beltrami, whose J1/J0 ratio
profile behaves like n ≈ 2–4) leaves the core streaming at full speed.

## Study C — loss topology

Per-channel loss rates Γ_ch ≈ (channel fraction)/τ_mean, ratio 0 → 3 (µs⁻¹):

| condition | Γ_end | Γ_radial | Γ_rad/Γ_end at ratio 3 |
|---|---|---|---|
| 3 mTorr | 1.26 → 0.54 | 0.000 → 0.001 | 0.2% |
| 10 mTorr (baseline) | 0.99 → 0.35 | 0.004 → 0.009 | 2.5% |
| 30 mTorr | 0.58 → 0.18 | 0.014 → 0.020 | 11% |
| 50 G | 1.02 → 0.36 | 0.016 → 0.029 | 8% |
| 200 G | 1.02 → 0.33 | 0.000 → 0.002 | 0.5% |
| 10 eV | 1.06 → 0.31 | 0.018 → 0.032 | 10% |

The entire confinement gain is **end-channel suppression** — Γ_end falls
~3× at every condition. The radial channel is parasitic: it *rises* with
pitch and with collisionality, and falls steeply with B (≈ r_L²ν collisional
cross-field steps accumulated over the longer dwell time). It stays 1–2
orders of magnitude below the end channel everywhere tested, but its growth
sets the eventual limit: the pitch gain must saturate when Γ_rad approaches
Γ_end, and at 30 mTorr or 10 eV the parasitic channel is already at ~10%.

## A predictive transport law

Combining the studies, electron confinement in these geometries follows

    τ(profile, ratio) ≈ τ₀ · ⟨|B|/Bz⟩^α,   α = α(L/λ_mfp)

with α = 1 collisionless, rising toward 2 with parallel collisionality, and
the entire field geometry entering through the single scalar ⟨|B|/Bz⟩
(seeding-averaged). Out-of-sample test: with α = 1.33 fitted on the screw
scan alone, the law predicts the *Beltrami* medians — a profile never used in
any fit — to within 1–2% at ratios 1.5, 2, and 3, and all powerlaw profiles
to within ±8%:

| profile at ratio 3 | ⟨F⟩ | predicted τ_med | measured | error |
|---|---|---|---|---|
| screw (in-sample) | 2.09 | 1.76 | 1.75 | +1% |
| powerlaw n=0.5 | 2.50 | 2.23 | 2.43 | −8% |
| powerlaw n=2 | 1.65 | 1.28 | 1.21 | +6% |
| powerlaw n=4 | 1.29 | 0.92 | 0.89 | +4% |
| **beltrami (out-of-sample)** | 1.39 | 1.02 | 1.01 | +1% |

Caveats: test-particle, no ambipolar electric field, no RF sustainment — the
absolute numbers are the collisionless-loss skeleton, not a device prediction.
The claims these runs support are the *existence, shape, mechanism, and
scaling law* of the pitch dependence, which is exactly what Prediction 4 asks
an experiment to test.

## Phase 3 — ambipolar electric fields

`poisson.js` solves the axisymmetric Poisson equation on a cell-centered
(r,z) grid (finite volumes, SOR, grounded Dirichlet walls; validated against
a manufactured solution to 0.01% — see selftest). `run-ambipolar.js` closes
the loop: CIC deposit → ρ = e(n_i − n_e) → solve → E = −∇φ → Boris with
half-E kicks, every 5 steps, over a **frozen ion background** equal to the
initial electron cloud (immobile-ion approximation on electron timescales).

Numerical validity (printed at startup, learned the hard way): the
full-depletion potential e·n0·R²/4ε₀ must sit at a few Te — at helicon
densities the frozen-ion quasi-static model runs away to kV; the default
n0 = 3e11 m⁻³ caps it at 4.5 Te and resolves λ_D on the default grid.
Feedback-off runs reproduce the Phase 1 runner within error bars from an
independent ensemble-stepping implementation (cross-validation of both).

**Does the pitch effect survive self-generated fields?** Screw field,
100 G, 10 mTorr, ratios 0 / 1 / 3:

| | τ_mean (µs) | survival at 20 µs | φ_well |
|---|---|---|---|
| feedback OFF | 1.02 / 1.19 / 2.76 | ~0% | — |
| feedback ON | >6.9 / >7.1 / >8.1 | 22.3 / 23.2 / 24.2% | 11.5 V ≈ 3.8 Te |

(Control at doubled n0: well deepens to ~6 Te and the same +10–18% pitch
ordering holds, so the conclusion is not an artifact of the charge cap.)

The ambipolar well **dominates and compresses the pitch effect** — from
2.7× to ~1.1–1.2× under these conditions. Once the escape of the fast tail
charges the plasma to a Te-scale potential, the barrier — whose depth is set
by charge balance, not field geometry — controls electron loss, and it is
nearly the same at every pitch. The pitch dependence survives only as a
consistent ordering (survival, τ, and leak-phase transport all still rank
with pitch) plus the familiar growth of the radial channel.

**Interpretation, and the sharpened hypothesis.** This model freezes the
ions, which at 100 G are unmagnetized (H_i ≈ 0.2 — the paper's own regime
analysis) and therefore pitch-blind. In a real discharge the steady state is
ambipolar: the potential adjusts until ion and electron losses balance, so
bulk particle confinement tracks the *ion* channel — which ignores pitch
below ion magnetization. The Phase 3 result therefore refines Prediction 4
into a regime statement:

> Below the ion-magnetization boundary, pitch should primarily shape
> electron energy transport, profiles, and wall heat flux rather than bulk
> particle confinement; strong pitch control of bulk confinement should
> switch on across H_i ≈ 1 (≈530 G at 10 mTorr — point B on the paper's
> regime map).

That is directly testable by the paper's proposed B-field scan, and it is a
falsifiable refinement this simulator produced rather than assumed.

### Regime I structure: the well is flat, the profiles move

Quantitative follow-up (`analyze-fields.js` on the per-ratio field
snapshots; 6-point pitch scan, frozen ions, 10 µs):

| ratio | φ_max V | surv% | τ_leak µs | radial% | ne peak radius m |
|---|---|---|---|---|---|
| 0 | 9.31 | 32.0 | 34.8 | 0.3 | 0.028 |
| 0.5 | 9.53 | 31.8 | 33.9 | 0.6 | 0.032 |
| 1 | 9.13 | 33.0 | 28.6 | 1.0 | 0.047 |
| 1.5 | 9.34 | 30.3 | 22.3 | 1.3 | 0.045 |
| 2 | 8.73 | 33.6 | 20.6 | 2.9 | 0.047 |
| 3 | 10.01 | 34.9 | 17.6 | 3.8 | 0.045 |

1. **The well is rigorously pitch-flat** — depth scatter ±7% with no trend,
   and the *shape* is fixed too (radial half-max 0.061–0.068 m, axial FWHM
   0.29–0.30 m at every ratio). This is the quantitative reason the
   confinement gain collapses: the barrier that now controls loss simply
   does not respond to field geometry.
2. **Pitch opens the slow radial drain on the trapped population** — τ_leak
   falls monotonically 34.8 → 17.6 µs as the parasitic radial channel
   (Study C) grows from 0.3% to 3.8%. In Regime I, pitch trades a small
   early-time confinement gain for a faster late-time leak: it redistributes
   loss between channels rather than reducing it.
3. **Pitch reshapes where the trapped electrons live** even though it cannot
   deepen the well: the density peak moves off-axis by ~60% (0.028 →
   0.045 m), the on-axis hollow partially fills only at ratio 3
   (core/edge 0 → 0.41), and comparing the 10 µs and 20 µs snapshots the
   high-pitch profile is *stationary in time* while the ratio-0 profile
   degrades as the population decays (σ/μ 25.7% → 45.5% at ratio 0 vs
   ~29–31% at both times for ratio 3). In an ambipolarly confined plasma,
   pitch acts on the density *distribution* — the processing-relevant
   quantity — not on the barrier.

Proving the ion side numerically needs mobile ions, which the runner now
supports (`--ions=mobile`, with `--mion` setting a reduced ion/neutral mass
so the H_i = 1 boundary falls in an affordable field range; Γ_i/Γ_e → 1 is
reported as the ambipolar-equilibrium check). The B-scan crossing H_i = 1 —
pitch gain vs ion Hall parameter — is the current experiment.

## Roadmap position

The project's driving question has evolved with the results: from
"pitch → confinement" (Phase 1, answered: yes, via a validated geometric
transport law) to "pitch → density distribution → potential structure →
transport" (Phase 3: the plasma's self-organization decides which pitch
effects survive). Current sequence: mobile-ion B-scan across H_i = 1
(in progress — the direct test of the regime-gated hypothesis) → floating
walls / sheath model → biased wafer chuck → helicon coupling models →
dusty plasma. RF wave solvers and full PIC remain deliberately out of scope
until the electrostatic story is complete.
