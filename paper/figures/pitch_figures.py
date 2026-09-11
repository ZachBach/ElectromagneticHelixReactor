"""Figures and audited numbers for paper #2 (paper/pitch-confinement.tex).

"Pitch-Dependent Electron Confinement in Helical Magnetic Geometries" takes every
figure and every quoted number from simulator runs archived in ../data/ -- copies
of sim/output/ CSVs, which are gitignored. Each file is checked against the run
parameters the paper attributes to it before it is used (a mismatch stops the
script), and every derived number the text quotes is printed, so the draft can
be audited against its data.

Outputs, next to this script:
  pitch-transport-law.pdf/.png          Fig. 1: field geometry, pitch scans, collapse on <|B|/Bz>
  pitch-gain-hierarchy.pdf/.png         Fig. 2: pitch gain from test particles to steady state
  pitch-equilibrium-structure.pdf/.png  Fig. 3: sustained-discharge radial structure

Run: python pitch_figures.py
"""

import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# ---- constants and baseline, matching sim/constants.js and the runners ----
QE = 1.602176634e-19
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
AMU = 1.66053906660e-27
M_AR = 39.948 * AMU
SIGMA_IN = 5.0e-19 + 5.0e-19          # SIGMA_CX + SIGMA_I_EL, m^2
R_WALL, L_CHAMBER = 0.1, 0.4          # m
TE_EV, T_GAS = 3.0, 300.0
N0_PHASE3 = 3e11                      # m^-3, space-charge scale of the Phase 3 runs
R_CORE = 0.4                          # core radius (units of R) for the inventory fraction

# ---- ink and marks (series colours checked with the dataviz palette validator) ----
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"

# Twist profiles are an ordered family (twist moves outward as n grows; darkest =
# most twist in the core) plus the Beltrami comparator in orange. Markers differ
# as well, so identity never rides on colour alone -- the figures must survive
# greyscale print.
PROFILES = [  # key, legend label, field model, power-law exponent, colour, marker
    ("n05", r"power law $n=0.5$", "powerlaw", 0.5, "#104281", "v"),
    ("screw", r"screw ($n=1$)", "screw", 1.0, "#256abf", "o"),
    ("n2", r"power law $n=2$", "powerlaw", 2.0, "#5598e7", "s"),
    ("n4", r"power law $n=4$", "powerlaw", 4.0, "#86b6ef", "D"),
    ("beltrami", "Beltrami", "beltrami", None, ORANGE, "^"),
]
RATIO_COLOURS = {0: "#86b6ef", 1: "#2a78d6", 3: "#104281"}   # one-hue ramp, darker = more pitch

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "grid.linestyle": "-",
    "axes.axisbelow": True, "legend.frameon": False,
    "lines.linewidth": 1.4, "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
    "savefig.dpi": 300, "pdf.fonttype": 42,
})


# ======================================================================
# data access
# ======================================================================

def load(name, **expect):
    """Parse a simulator CSV ('#' metadata lines, one header row, numeric rows).
    Stops if any expected run parameter differs from what the file records."""
    meta, header, rows = {}, None, []
    with open(os.path.join(DATA, name)) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                meta.update(re.findall(r"([A-Za-z_][A-Za-z0-9_]*)=(\S+)", line))
            elif header is None:
                header = line.split(",")
            else:
                rows.append([float(v) for v in line.split(",")])
    for key, want in expect.items():
        got = meta.get(key)
        ok = got is not None and (
            np.isclose(float(got), float(want)) if isinstance(want, (int, float)) else got == want)
        if not ok:
            raise SystemExit(f"{name}: expected {key}={want}, file records {got}")
    table = np.array(rows)
    return meta, {h: table[:, i] for i, h in enumerate(header)}


def at(c, col, ratio):
    hit = np.isclose(c["ratio"], ratio)
    if not hit.any():
        raise SystemExit(f"ratio {ratio} not in run")
    return float(c[col][hit][0])


def load_fields(name):
    """Time-averaged (r, z) snapshot -> r, z, and phi / ne / ni as [z, r] grids."""
    table = np.loadtxt(os.path.join(DATA, name), delimiter=",", skiprows=1)
    r, z = np.unique(table[:, 0]), np.unique(table[:, 1])
    grid = table[:, 2:].reshape(len(z), len(r), 3)   # rows run z-outer, r-inner
    return r, z, grid[..., 0], grid[..., 1], grid[..., 2]


def structure(r, z, phi, ne):
    """Python replica of sim/analyze-fields.js (same definitions, same bands),
    plus the core inventory fraction used in the paper."""
    nz, nr = phi.shape
    band = slice(int(nz * 0.25), int(nz * 0.75))
    phi_max = phi.max()
    plateau = phi[band, 0].mean()
    jm = nz // 2
    below = np.nonzero(phi[jm, :] < phi_max / 2)[0]
    r_half = r[below[0]] if below.size else r[-1]
    above = np.nonzero(phi[:, 0] >= phi_max / 2)[0]
    z_fwhm = z[above[-1]] - z[above[0]]
    prof = ne[band, :].mean(axis=0)
    r_peak = r[np.argmax(prof)]
    dr = r[1] - r[0]
    i07 = min(nr - 1, int(np.floor(0.7 * (r[-1] + dr / 2) / dr)))
    core_edge = prof[0] / prof[i07] if prof[i07] > 0 else np.inf
    w = r[: i07 + 1]
    mu = np.sum(prof[: i07 + 1] * w) / np.sum(w)
    unif = np.sqrt(np.sum((prof[: i07 + 1] - mu) ** 2 * w) / np.sum(w)) / mu
    return dict(phi_max=phi_max, plateau=plateau, r_half=r_half, z_fwhm=z_fwhm,
                r_peak=r_peak, core_edge=core_edge, unif=unif, prof=prof,
                phi_prof=phi[band, :].mean(axis=0),
                f_core=core_fraction(r, ne, 0.25, 0.75),
                f_core_halves=(core_fraction(r, ne, 0.25, 0.5), core_fraction(r, ne, 0.5, 0.75)))


def core_fraction(r, ne, z_lo, z_hi):
    """Share of the electron inventory inside r < R_CORE*R, over a z-band. An
    integral, so far less sensitive to per-cell sampling noise than a peak radius;
    the two half-bands give a sampling-noise estimate."""
    nz = ne.shape[0]
    w = ne[int(nz * z_lo):int(nz * z_hi), :].mean(axis=0) * r
    return float(w[r < R_CORE * R_WALL].sum() / w.sum())


# ======================================================================
# physics helpers (field geometry, Hall parameters, Debye length)
# ======================================================================

def bessel_j(n, x):
    """J_n by series, as sim/field.js; accurate for the x < 2.405 used here."""
    x = np.asarray(x, dtype=float)
    term = np.ones_like(x)
    for k in range(1, n + 1):
        term = term * (x / 2) / k
    total = term.copy()
    x2 = -(x * x) / 4
    for k in range(1, 40):
        term = term * x2 / (k * (k + n))
        total = total + term
    return total


def beltrami_alpha_r(ratio):
    lo, hi = 1e-9, 2.4048
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if float(bessel_j(1, mid) / bessel_j(0, mid)) < ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def path_factor_profile(model, ratio, r, nexp=1.0):
    """Local field-line path factor F(r) = |B|/Bz; the field normalization cancels."""
    if model in ("screw", "powerlaw"):
        n = 1.0 if model == "screw" else nexp
        t = ratio * (r / R_WALL) ** n
    elif model == "beltrami":
        ar = beltrami_alpha_r(ratio) * r / R_WALL
        t = bessel_j(1, ar) / bessel_j(0, ar)
    else:
        raise ValueError(model)
    return np.sqrt(1 + t * t)


def mean_path_factor(model, ratio, nexp=1.0):
    """<|B|/Bz> over the seeding region r < 0.9R, weight r dr -- as sim/fit-alpha.js."""
    r = 0.9 * R_WALL * (np.arange(400) + 0.5) / 400
    return float(np.sum(path_factor_profile(model, ratio, r, nexp) * r) / np.sum(r))


def fit_slope(x, y):
    """Least squares y = a + slope*x with the slope standard error of sim/fit-alpha.js."""
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    den = n * np.sum(x * x) - np.sum(x) ** 2
    slope = (n * np.sum(x * y) - np.sum(x) * np.sum(y)) / den
    a = (np.sum(y) - slope * np.sum(x)) / n
    sse = np.sum((y - a - slope * x) ** 2)
    se = np.sqrt((sse / max(n - 2, 1)) * (n / den))
    return a, slope, se


def neutral_density(p_mtorr):
    return p_mtorr * 0.133322 / (KB * T_GAS)


def hall_ion(b_gauss, p_mtorr, m_ion, speed):
    return (QE * b_gauss * 1e-4 / m_ion) / (neutral_density(p_mtorr) * SIGMA_IN * speed)


def c_s(m_ion):
    return np.sqrt(TE_EV * QE / m_ion)


def vbar_thermal(m_ion):
    return np.sqrt(8 * KB * T_GAS / (np.pi * m_ion))


def debye_length(n):
    return np.sqrt(EPS0 * TE_EV / (n * QE))


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2.5, width=0.6)


def panel_label(ax, text):
    ax.text(-0.02, 1.04, text, transform=ax.transAxes, fontsize=9, fontweight="bold",
            ha="right", va="bottom", color=INK)


# ======================================================================
# manifest -- which archived run backs which claim
# ======================================================================

P1 = dict(species="electron", Bwall_T=0.01, Tseed_eV=3, R_m=0.1, L_m=0.4)
FULL = {  # Boris + MCC at the baseline (10 mTorr, 100 G), N = 4000, seed 1
    "screw": ("pitch-scan-2026-07-19-05-10-21.csv", dict(model="screw")),
    "beltrami": ("pitch-scan-2026-07-19-05-10-47.csv", dict(model="beltrami")),
    "n05": ("pitch-scan-2026-07-19-05-51-50.csv", dict(model="powerlaw", nexp=0.5)),
    "n2": ("pitch-scan-2026-07-19-05-52-16.csv", dict(model="powerlaw", nexp=2)),
    "n4": ("pitch-scan-2026-07-19-05-52-34.csv", dict(model="powerlaw", nexp=4)),
}
FREE = {  # free-streaming geometric model, N = 20000
    "screw": ("mechanism-screw-2026-07-19-05-53-21.csv", dict(model="screw")),
    "beltrami": ("mechanism-beltrami-2026-09-10-21-43-38.csv", dict(model="beltrami")),
    "n05": ("mechanism-powerlaw-n0.5-2026-07-19-05-53-21.csv", dict(model="powerlaw", nexp=0.5)),
    "n2": ("mechanism-powerlaw-n2-2026-07-19-05-53-21.csv", dict(model="powerlaw", nexp=2)),
    "n4": ("mechanism-powerlaw-n4-2026-07-19-05-53-21.csv", dict(model="powerlaw", nexp=4)),
}
NOCOLL = {  # Boris with collisions off
    "screw": ("pitch-scan-2026-07-19-05-31-05.csv", dict(model="screw", nocoll="true")),
    "beltrami": ("pitch-scan-2026-07-19-05-31-15.csv", dict(model="beltrami", nocoll="true")),
}
SCANS = {  # one-factor-at-a-time around the baseline, screw field, N = 4000
    "baseline": ("pitch-scan-2026-07-19-05-10-21.csv", dict()),
    "3 mTorr": ("pitch-scan-2026-07-19-05-33-18.csv", dict(p_mTorr=3)),
    "30 mTorr": ("pitch-scan-2026-07-19-05-33-30.csv", dict(p_mTorr=30)),
    "1 eV": ("pitch-scan-2026-07-19-05-32-38.csv", dict(Tseed_eV=1)),
    "10 eV": ("pitch-scan-2026-07-19-05-33-00.csv", dict(Tseed_eV=10)),
    "50 G": ("pitch-scan-2026-07-19-05-34-24.csv", dict(Bwall_T=0.005)),
    "200 G": ("pitch-scan-2026-07-19-05-34-33.csv", dict(Bwall_T=0.02)),
    "seed 2": ("pitch-scan-2026-07-19-05-35-08.csv", dict(seed=2)),
    "dt/2": ("pitch-scan-2026-07-19-05-31-45.csv", dict(dtscale=0.5)),
}
FROZEN = {  # frozen-ion electrostatic feedback, screw, 100 G, 10 mTorr, N = 2000
    "off": ("ambipolar-screw-fboff-2026-07-19-08-57-46.csv", dict(feedback="false", n0_m3=3e11, Tmax_s=2e-5)),
    "on": ("ambipolar-screw-fbon-2026-07-19-08-57-51.csv", dict(feedback="true", n0_m3=3e11, Tmax_s=2e-5)),
    "on, 2 n0": ("ambipolar-screw-fbon-2026-07-19-09-04-33.csv", dict(feedback="true", n0_m3=6e11, Tmax_s=1e-5)),
    "6-point": ("ambipolar-screw-fbon-2026-07-19-09-20-38.csv", dict(feedback="true", n0_m3=3e11, Tmax_s=1e-5)),
}
FROZEN_FIELDS = {rt: f"ambipolar-fields-screw-fbon-ratio{rt:g}-2026-07-19-09-20-38.csv"
                 for rt in (0, 0.5, 1, 1.5, 2, 3)}
MOBILE = {  # decaying ensembles, mobile 1-amu ions, N = 1000, 25 us; (mTorr, gauss) -> file
    (10, 25): "ambipolar-screw-fbon-mob25G-2026-07-19-09-24-21.csv",
    (10, 50): "ambipolar-screw-fbon-mob50G-2026-07-19-09-30-12.csv",
    (10, 100): "ambipolar-screw-fbon-mob100G-2026-07-19-09-34-55.csv",
    (10, 200): "ambipolar-screw-fbon-mob200G-2026-07-19-09-40-40.csv",
    (10, 400): "ambipolar-screw-fbon-mob400G-2026-07-19-09-47-01.csv",
    (3, 50): "ambipolar-screw-fbon-mob50G-2026-07-19-09-56-42.csv",
    (3, 100): "ambipolar-screw-fbon-mob100G-2026-07-19-10-01-08.csv",
    (3, 200): "ambipolar-screw-fbon-mob200G-2026-07-19-10-06-41.csv",
    (3, 400): "ambipolar-screw-fbon-mob400G-2026-07-19-10-12-56.csv",
}
SUSTAINED = {  # source-sustained steady state, mobile 1-amu ions, 200 pairs/us
    "low": ("sustained-screw-100G-10mT-2026-07-19-10-56-23.csv", 100, 10, 5e-5, (0, 1, 3)),
    "high": ("sustained-screw-400G-3mT-2026-07-19-10-56-25.csv", 400, 3, 4e-5, (0, 3)),
}
SUSTAINED_FIELDS = {
    ("low", 0): "sustained-fields-screw-100G-10mT-ratio0-2026-07-19-10-56-23.csv",
    ("low", 1): "sustained-fields-screw-100G-10mT-ratio1-2026-07-19-10-56-23.csv",
    ("low", 3): "sustained-fields-screw-100G-10mT-ratio3-2026-07-19-10-56-23.csv",
    ("high", 0): "sustained-fields-screw-400G-3mT-ratio0-2026-07-19-10-56-25.csv",
    ("high", 3): "sustained-fields-screw-400G-3mT-ratio3-2026-07-19-10-56-25.csv",
}

full = {k: load(f, **P1, p_mTorr=10, n=4000, seed=1, **kw)[1] for k, (f, kw) in FULL.items()}
free = {k: load(f, species="electron", Bwall_T=0.01, Tseed_eV=3, n=20000, seed=1, **kw)[1]
        for k, (f, kw) in FREE.items()}
nocoll = {k: load(f, **P1, p_mTorr=10, n=4000, seed=1, **kw)[1] for k, (f, kw) in NOCOLL.items()}
scans = {}
for label, (f, kw) in SCANS.items():
    want = dict(P1, model="screw", p_mTorr=10, n=4000, seed=1)
    want.update(kw)
    scans[label] = load(f, **want)[1]
frozen = {k: load(f, model="screw", Bwall_T=0.01, p_mTorr=10, n=2000, Tseed_eV=3, **kw)[1]
          for k, (f, kw) in FROZEN.items()}
mobile = {}
for (p, b), f in MOBILE.items():
    meta, c = load(f, ions="mobile", mion_amu=1, feedback="true", model="screw", n=1000,
                   Bwall_T=b * 1e-4, p_mTorr=p, Tmax_s=2.5e-5, n0_m3=3e11)
    his = hall_ion(b, p, AMU, c_s(AMU))
    if "HiStar" in meta and not np.isclose(float(meta["HiStar"]), his, rtol=1e-6):
        raise SystemExit(f"{f}: H_i* recomputes to {his}, file records {meta['HiStar']}")
    mobile[(p, b)] = (his, c)
sustained = {}
for key, (f, b, p, ttotal, ratios) in SUSTAINED.items():
    meta, c = load(f, model="screw", mion_amu=1, Bwall_T=b * 1e-4, p_mTorr=p,
                   rate_per_us=200, Ttotal_s=ttotal)
    his = hall_ion(b, p, AMU, c_s(AMU))
    if not np.isclose(float(meta["HiStar"]), his, rtol=1e-6):
        raise SystemExit(f"{f}: H_i* recomputes to {his}, file records {meta['HiStar']}")
    sustained[key] = (b, p, his, ratios, c)
sustained_structure = {k: structure(*load_fields(f)[:4]) for k, f in SUSTAINED_FIELDS.items()}
frozen_structure = {rt: structure(*load_fields(f)[:4]) for rt, f in FROZEN_FIELDS.items()}


# ======================================================================
# numbers quoted in the text
# ======================================================================

def rule(title):
    print(f"\n== {title} ==")


rule("Path factor <|B|/Bz> over the seeding region")
for key, label, model, nexp, _, _ in PROFILES:
    print(f"  {key:9s}", "  ".join(f"r={rt:g}:{mean_path_factor(model, rt, nexp or 1):.3f}"
                                    for rt in (0.5, 1, 1.5, 2, 3)))

rule("Phase 1 pitch scans, baseline (tau in us)")
for key, *_ in PROFILES:
    c = full[key]
    rows = "  ".join(f"{rt:g}:{m * 1e6:.2f}+-{s * 1e6:.2f}"
                     for rt, m, s in zip(c["ratio"], c["tauMean_s"], c["tauSE_s"]))
    print(f"  {key:9s} {rows}")
TAU0_MEAN, TAU0_MED = at(full["screw"], "tauMean_s", 0), at(full["screw"], "tauMedian_s", 0)
G_TEST = at(full["screw"], "tauMean_s", 3) / TAU0_MEAN
print(f"  screw gain tau(3)/tau(0): mean {G_TEST:.2f}, median "
      f"{at(full['screw'], 'tauMedian_s', 3) / TAU0_MED:.2f}")
for key in ("n05", "n2", "n4", "beltrami"):
    print(f"  {key} gain (mean): {at(full[key], 'tauMean_s', 3) / TAU0_MEAN:.2f}")
print(f"  n=0.5 vs n=4 at ratio 3 (mean): "
      f"{at(full['n05'], 'tauMean_s', 3) / at(full['n4'], 'tauMean_s', 3):.2f}x")

rule("Mechanism layers, screw: tau(0) -> tau(3) [us], gain")
FS0_MEAN, FS0_MED = at(free["screw"], "tauMean_s", 0), at(free["screw"], "tauMedian_s", 0)
for label, c in (("free streaming", free["screw"]), ("Boris, no collisions", nocoll["screw"]),
                 ("Boris + MCC", full["screw"])):
    for col in ("tauMean_s", "tauMedian_s"):
        t0, t3 = at(c, col, 0), at(c, col, 3)
        print(f"  {label:22s} {col:12s} {t0 * 1e6:.3f} -> {t3 * 1e6:.3f}  gain {t3 / t0:.2f}")
    print(f"  {'':22s} survivors at ratio 3: {100 * at(c, 'surviveFrac', 3):.2f}%")
print("  drift check, median(Boris no-coll)/median(free streaming) - 1:")
for key in ("screw", "beltrami"):
    diffs = [f"{rt:g}:{100 * (at(nocoll[key], 'tauMedian_s', rt) / at(free[key], 'tauMedian_s', rt) - 1):+.1f}%"
             for rt in nocoll[key]["ratio"]]
    print(f"    {key:9s}", "  ".join(diffs))
    print(f"    {key:9s} radial loss fraction, no collisions: max {nocoll[key]['fracRadialLoss'].max():.4f}")
print(f"  free-streaming Beltrami tau_mean(3) = {at(free['beltrami'], 'tauMean_s', 3) * 1e6:.2f} us "
      f"vs screw {at(free['screw'], 'tauMean_s', 3) * 1e6:.2f} us")

rule("Transport-law exponent alpha (median fit, mean fit)")


def alpha_of(c, model, nexp):
    x = np.log([mean_path_factor(model, rt, nexp) for rt in c["ratio"]])
    return fit_slope(x, np.log(c["tauMedian_s"])), fit_slope(x, np.log(c["tauMean_s"])), len(x)


ALPHA = {}
for group, runs in (("Boris + MCC", full), ("free streaming", free)):
    for key, label, model, nexp, _, _ in PROFILES:
        (a_md, s_md, e_md), (_, s_mn, e_mn), npts = alpha_of(runs[key], model, nexp or 1)
        ALPHA[(group, key)] = (a_md, s_md)
        print(f"  {group:15s} {key:9s} alpha_med {s_md:.2f}+-{e_md:.2f}  alpha_mean {s_mn:.2f}+-{e_mn:.2f}  ({npts} pts)")
for key in ("screw", "beltrami"):
    (_, s_md, e_md), _, npts = alpha_of(nocoll[key], key, 1)
    print(f"  {'no collisions':15s} {key:9s} alpha_med {s_md:.2f}+-{e_md:.2f}  ({npts} pts)")
for label, c in scans.items():
    (_, s_md, e_md), (_, s_mn, e_mn), npts = alpha_of(c, "screw", 1)
    print(f"  {'screw scan':15s} {label:9s} alpha_med {s_md:.2f}+-{e_md:.2f}  alpha_mean {s_mn:.2f}+-{e_mn:.2f}  ({npts} pts)")

rule("Out-of-sample: screw median fit predicts the other profiles")
A_SCREW, ALPHA_SCREW = ALPHA[("Boris + MCC", "screw")]
for key, label, model, nexp, _, _ in PROFILES:
    c = full[key]
    cells = []
    for rt, med in zip(c["ratio"], c["tauMedian_s"]):
        if rt < 0.5:
            continue
        F = mean_path_factor(model, rt, nexp or 1)
        pred = np.exp(A_SCREW) * F ** ALPHA_SCREW
        cells.append(f"{rt:g}: F={F:.2f} pred {pred * 1e6:.2f} meas {med * 1e6:.2f} ({100 * (pred / med - 1):+.0f}%)")
    print(f"  {key:9s}", "  ".join(cells))

rule("Loss topology: channel rates (1/us) ratio 0 -> 3")
for label, c in scans.items():
    ge = [at(c, "fracEndLoss", rt) / at(c, "tauMean_s", rt) * 1e-6 for rt in (0, 3)]
    gr = [at(c, "fracRadialLoss", rt) / at(c, "tauMean_s", rt) * 1e-6 for rt in (0, 3)]
    print(f"  {label:9s} tau {at(c, 'tauMean_s', 0) * 1e6:.2f} -> {at(c, 'tauMean_s', 3) * 1e6:.2f}  "
          f"G_end {ge[0]:.2f} -> {ge[1]:.2f}   G_rad {1e3 * gr[0]:.2f} -> {1e3 * gr[1]:.2f} e-3   "
          f"G_rad/G_end(3) {100 * gr[1] / ge[1]:.2f}%   gain {at(c, 'tauMean_s', 3) / at(c, 'tauMean_s', 0):.2f}")

rule("End-to-end loss asymmetry (fracEndLow - fracEndHigh), significance")
for label, c in list(scans.items()) + [("no-coll", nocoll["screw"]), ("n=0.5", full["n05"])]:
    if "fracEndLow" not in c:
        continue
    cells = []
    for rt in (0, 3):
        lo, hi, nl = at(c, "fracEndLow", rt), at(c, "fracEndHigh", rt), at(c, "nLost", rt)
        d = lo - hi
        cells.append(f"r={rt}: {d:+.3f} ({d / np.sqrt((lo + hi - d * d) / nl):+.1f} sigma)")
    print(f"  {label:9s}", "   ".join(cells), f"  radial(3) {at(c, 'fracRadialLoss', 3):.3f}")

rule("Frozen-ion feedback")
for label, c in frozen.items():
    g = at(c, "tauMean_s", 3) / at(c, "tauMean_s", 0)
    print(f"  {label:9s}", "  ".join(
        f"{rt:g}: tau {m * 1e6:.2f} surv {100 * s:.1f}% leak {lk * 1e6:.1f} phi {ph:.2f} rad {100 * rd:.1f}%"
        for rt, m, s, lk, ph, rd in zip(c["ratio"], c["tauMean_s"], c["surviveFrac"], c["tauLeak_s"],
                                        c["phiMax_V"], c["fracRadialLoss"])), f"| gain {g:.2f}")
G_FROZEN = at(frozen["on"], "tauMean_s", 3) / at(frozen["on"], "tauMean_s", 0)
for n0 in (3e11, 6e11):
    lam = debye_length(n0)
    cap = QE * n0 * R_WALL ** 2 / (4 * EPS0)
    print(f"  n0 = {n0:.0e}: lambda_D = {100 * lam:.2f} cm, R/lambda_D = {R_WALL / lam:.1f}, "
          f"lambda_D/dr = {lam / (R_WALL / 48):.1f}, full-depletion phi = {cap:.1f} V = {cap / TE_EV:.1f} Te")
print("  structure (6-point scan snapshots): phi_max, plateau, r_half, z_fwhm, r_peak, core/edge, sigma/mu, f_core (halves)")
for rt, s in frozen_structure.items():
    print(f"    {rt:>3g}: {s['phi_max']:.2f} {s['plateau']:.2f} {s['r_half']:.3f} {s['z_fwhm']:.3f} "
          f"{s['r_peak']:.3f} {s['core_edge']:.2f} {100 * s['unif']:.1f}%  "
          f"f_core {s['f_core']:.3f} ({s['f_core_halves'][0]:.3f}, {s['f_core_halves'][1]:.3f})")

rule("Mobile ions, decaying ensembles")
for (p, b), (his, c) in sorted(mobile.items()):
    hi_th = hall_ion(b, p, AMU, vbar_thermal(AMU))
    g = at(c, "tauMean_s", 3) / at(c, "tauMean_s", 0)
    print(f"  {p:>2} mTorr {b:>3} G  H_i {hi_th:.2f}  H_i* {his:.2f}  gain {g:.2f}  "
          f"survE {100 * at(c, 'surviveFrac', 0):.1f}->{100 * at(c, 'surviveFrac', 3):.1f}%  "
          f"survI {100 * at(c, 'surviveFracIon', 0):.1f}->{100 * at(c, 'surviveFracIon', 3):.1f}%  "
          f"Gi/Ge {at(c, 'GiOverGe', 0):.2f}  phi {at(c, 'phiMax_V', 0):.2f}")

rule("Sustained steady state")
W_SUPER = N0_PHASE3 * np.pi * (0.9 * R_WALL) ** 2 * (0.8 * L_CHAMBER) / 2000
V_CHAMBER = np.pi * R_WALL ** 2 * L_CHAMBER
for key, (b, p, his, ratios, c) in sustained.items():
    t0 = at(c, "tauEff_s", 0)
    base = sustained_structure[(key, 0)]
    for rt in ratios:
        ne_mean = at(c, "NeSS", rt) * W_SUPER / V_CHAMBER
        ni_mean = at(c, "NiSS", rt) * W_SUPER / V_CHAMBER
        s = sustained_structure[(key, rt)]
        print(f"  {b} G {p} mTorr (H_i*={his:.2f}) ratio {rt}: tau_eff {at(c, 'tauEff_s', rt) * 1e6:.2f} us "
              f"({100 * (at(c, 'tauEff_s', rt) / t0 - 1):+.1f}%), Gi/Ge {at(c, 'GiOverGe', rt):.2f}, "
              f"e-radial {100 * at(c, 'fracRadialE', rt):.1f}%, Ni/Ne {ni_mean / ne_mean:.2f}, "
              f"<ne> {ne_mean:.2e}, R/lambda_D {R_WALL / debye_length(ne_mean):.1f}")
        print(f"      phi_max {s['phi_max']:.2f} plateau {s['plateau']:.2f} ({100 * (s['plateau'] / base['plateau'] - 1):+.0f}%) "
              f"r_half {s['r_half']:.3f} z_fwhm {s['z_fwhm']:.3f} r_peak {s['r_peak']:.3f} "
              f"core/edge {s['core_edge']:.2f} sigma/mu {100 * s['unif']:.1f}%")
        print(f"      f_core(r<{R_CORE}R) {s['f_core']:.3f} (half-bands {s['f_core_halves'][0]:.3f}, "
              f"{s['f_core_halves'][1]:.3f}; change vs ratio 0 {100 * (s['f_core'] / base['f_core'] - 1):+.0f}%)")

rule("Argon ion-magnetization field, sim cross sections")
for p in (10, 3):
    b_th = 100 / hall_ion(100, p, M_AR, vbar_thermal(M_AR))
    b_star = 100 / hall_ion(100, p, M_AR, c_s(M_AR))
    print(f"  {p:>2} mTorr: H_i = 1 at {b_th:.0f} G (thermal), H_i* = 1 at {b_star:.0f} G (exit speed), "
          f"ratio {b_star / b_th:.2f}")
print(f"  1 amu model gas at 10 mTorr: H_i* = 1 at {100 / hall_ion(100, 10, AMU, c_s(AMU)):.0f} G; "
      f"at 3 mTorr: {100 / hall_ion(100, 3, AMU, c_s(AMU)):.0f} G")


# ======================================================================
# Fig. 1 -- field geometry, pitch scans, collapse onto <|B|/Bz>
# ======================================================================

fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(6.7, 2.55))
rr = np.linspace(0, 1, 300)
handles = []
for key, label, model, nexp, colour, marker in PROFILES:
    ax_a.plot(rr, path_factor_profile(model, 3.0, rr * R_WALL, nexp or 1), color=colour, lw=1.4)
    handles.append(Line2D([], [], color=colour, marker=marker, lw=1.4, ms=4.5,
                          mec="white", mew=0.6, label=label))
ax_a.axvline(0.9, color=MUTED, lw=0.6)
ax_a.text(0.915, 1.06, "seeding edge", rotation=90, fontsize=6.5, color=INK2, ha="left", va="bottom")
ax_a.text(0.03, 3.2, r"$B_\theta(R)/B_z = 3$", fontsize=7.5, color=INK, ha="left", va="top")
ax_a.set_xlim(0, 1)
ax_a.set_ylim(1, 3.3)
ax_a.set_xlabel(r"$r/R$")
ax_a.set_ylabel(r"path factor $|B|/B_z$")

for key, label, model, nexp, colour, marker in PROFILES:
    c = full[key]
    ratios, mean, se = c["ratio"], c["tauMean_s"] * 1e6, c["tauSE_s"] * 1e6
    if not np.isclose(ratios[0], 0):   # every profile is the same uniform field at ratio 0
        ratios = np.r_[0, ratios]
        mean = np.r_[TAU0_MEAN * 1e6, mean]
        se = np.r_[at(full["screw"], "tauSE_s", 0) * 1e6, se]
    ax_b.errorbar(ratios, mean, yerr=se, color=colour, marker=marker, ms=4.2, mec="white", mew=0.6,
                  lw=1.3, elinewidth=0.8, capsize=0)
ax_b.set_xlabel(r"pitch ratio $B_\theta(R)/B_z$")
ax_b.set_ylabel(r"$\tau_{\mathrm{mean}}$  [$\mu$s]")
ax_b.set_xlim(-0.1, 3.15)
ax_b.set_ylim(0.8, 3.9)

f_line = np.linspace(1.0, 2.6, 60)
for key, label, model, nexp, colour, marker in PROFILES:
    for runs, tau0, filled in ((full, TAU0_MED, True), (free, FS0_MED, False)):
        c = runs[key]
        F = [mean_path_factor(model, rt, nexp or 1) for rt in c["ratio"]]
        ax_c.plot(F, c["tauMedian_s"] / tau0, ls="none", marker=marker, ms=4.4,
                  color=colour, mfc=colour if filled else "white", mec="white" if filled else colour,
                  mew=0.6 if filled else 0.9, zorder=4 if filled else 3)
ax_c.plot(f_line, np.exp(A_SCREW) / TAU0_MED * f_line ** ALPHA_SCREW, color=INK, lw=0.9, zorder=2)
ax_c.plot(f_line, f_line, color=MUTED, lw=0.9, zorder=2)
ax_c.text(1.9, np.exp(A_SCREW) / TAU0_MED * 1.9 ** ALPHA_SCREW * 1.13,
          rf"$\alpha={ALPHA_SCREW:.2f}$ (screw fit)", fontsize=6.8, color=INK, ha="right")
ax_c.text(2.05, 2.05 * 0.84, r"$\alpha=1$", fontsize=6.8, color=INK2, ha="left")
ax_c.set_xscale("log")
ax_c.set_yscale("log")
for axis, ticks in ((ax_c.xaxis, [1, 1.5, 2, 2.5]), (ax_c.yaxis, [1, 1.5, 2, 3, 4])):
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_locator(NullLocator())
ax_c.set_xlim(0.97, 2.65)
ax_c.set_ylim(0.9, 4.2)
ax_c.set_xlabel(r"mean path factor $\langle |B|/B_z \rangle$")
ax_c.set_ylabel(r"$\tau_{\mathrm{med}} \,/\, \tau_{\mathrm{med}}(0)$")
ax_c.legend(handles=[Line2D([], [], ls="none", marker="o", color=INK2, mec="white", mew=0.6, ms=4.4,
                            label="Boris + MCC, 10 mTorr"),
                     Line2D([], [], ls="none", marker="o", mfc="white", mec=INK2, mew=0.9, ms=4.4,
                            label="free streaming")],
            loc="upper left", fontsize=6.8, handletextpad=0.2, borderaxespad=0.2)

for ax, tag in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
    style(ax)
    panel_label(ax, tag)
fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0),
           handletextpad=0.4, columnspacing=1.2)
fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=1.2)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"pitch-transport-law.{ext}"))
plt.close(fig)


# ======================================================================
# Fig. 2 -- pitch gain from test particles to sustained steady state
# ======================================================================

fig, ax = plt.subplots(figsize=(5.3, 3.1))
ax.axhline(1.0, color=MUTED, lw=0.6, zorder=1)
ax.axvline(1.0, color=MUTED, lw=0.6, zorder=1)
ax.text(1.06, 2.55, r"$H_i^*=1$", fontsize=7, color=INK2, ha="left")
# reference levels carry direct labels in the right margin, clear of the data
for level, text in ((G_TEST, "test particles,\nno electric field"),
                    (G_FROZEN, "frozen ions,\nself-consistent field")):
    ax.axhline(level, color=INK2, lw=0.9, zorder=1)
    ax.annotate(f"{text}\n({level:.2f})", xy=(1.0, level), xycoords=("axes fraction", "data"),
                xytext=(5, 0), textcoords="offset points", ha="left", va="center",
                fontsize=7, color=INK2, annotation_clip=False)
for p, colour, marker in ((10, BLUE, "o"), (3, ORANGE, "s")):
    pts = sorted((his, at(c, "tauMean_s", 3) / at(c, "tauMean_s", 0))
                 for (pp, _), (his, c) in mobile.items() if pp == p)
    ax.plot(*zip(*pts), color=colour, marker=marker, ms=5.5, lw=1.4, mec="white", mew=0.8, zorder=3,
            label=f"mobile ions, decaying ensemble, {p} mTorr")
sus = [(his, at(c, "tauEff_s", 3) / at(c, "tauEff_s", 0)) for (_, _, his, _, c) in sustained.values()]
ax.plot(*zip(*sus), ls="none", color=AQUA, marker="D", ms=6.5, mec="white", mew=0.9, zorder=4,
        label="mobile ions, sustained steady state")
for his, g in sus:
    ax.annotate(f"{100 * (g - 1):+.0f}%", (his, g), xytext=(0, -11), textcoords="offset points",
                fontsize=7, color=INK, ha="center", va="top")
ax.set_xscale("log")
ax.xaxis.set_major_locator(FixedLocator([0.03, 0.1, 0.3, 1, 3]))
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlim(0.025, 4)
ax.set_ylim(0.75, 3.05)
ax.set_xlabel(r"exit-speed ion Hall parameter $H_i^* = \omega_{ci}/\nu_{in}(c_s)$")
ax.set_ylabel(r"pitch gain  $\tau(B_\theta/B_z=3)\,/\,\tau(0)$")
ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), fontsize=7, handletextpad=0.4)
style(ax)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"pitch-gain-hierarchy.{ext}"))
plt.close(fig)


# ======================================================================
# Fig. 3 -- sustained-discharge radial structure
# ======================================================================

fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.5), sharex=True, sharey="row")
ratio_handles = {}
for col, key in enumerate(("low", "high")):
    b, p, his, ratios, _ = sustained[key]
    for rt in ratios:
        r = load_fields(SUSTAINED_FIELDS[(key, rt)])[0]
        s = sustained_structure[(key, rt)]
        ratio_handles[rt], = axes[0, col].plot(r / R_WALL, s["prof"] / 1e11, color=RATIO_COLOURS[rt],
                                               lw=1.2, label=rf"$B_\theta(R)/B_z = {rt}$")
        axes[1, col].plot(r / R_WALL, s["phi_prof"], color=RATIO_COLOURS[rt], lw=1.2)
    axes[0, col].set_title(rf"{b} G, {p} mTorr  ($H_i^* = {his:.2f}$)", color=INK)
    for row in range(2):
        axes[row, col].axvline(0.7, color=MUTED, lw=0.6)
        style(axes[row, col])
    axes[1, col].set_xlabel(r"$r/R$")
axes[0, 0].text(0.69, 0.02, r"$r=0.7R$", transform=axes[0, 0].get_xaxis_transform(),
                fontsize=6.5, color=INK2, ha="right", va="bottom")
axes[0, 0].set_ylabel(r"$n_e$  [$10^{11}$ m$^{-3}$]")
axes[1, 0].set_ylabel(r"$\phi$  [V]")
axes[1, 0].set_xlim(0, 1)
axes[0, 0].set_ylim(bottom=0)
axes[1, 0].set_ylim(bottom=0)
for ax, tag in zip(axes.flat, "abcd"):
    panel_label(ax, tag)
fig.legend(handles=[ratio_handles[k] for k in (0, 1, 3)], loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 1.0), columnspacing=1.6)
fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=1.0, w_pad=1.2)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"pitch-equilibrium-structure.{ext}"))
plt.close(fig)

print("\nsaved pitch-transport-law, pitch-gain-hierarchy, pitch-equilibrium-structure (.pdf/.png)")
