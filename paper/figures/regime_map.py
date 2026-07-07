"""EHR regime map: magnetic field vs. neutral pressure with Hall-parameter contours.

Generates the B-p operating-regime figure for the EHR paper (Section: Magnetization
Regime). A species is magnetized when its Hall parameter H_s = omega_cs / nu_sn > 1,
with nu_sn = n_n * sigma_sn * vbar_s and n_n = p / (k_B T_gas). Both H_e and H_i
scale as B/p, so regime boundaries are straight lines on log-log axes.

Representative argon parameters (matching the paper):
  T_gas = 300 K, T_e = 3 eV, sigma_en = 2e-20 m^2, sigma_in = 1e-18 m^2.

Outputs ehr-regime-map.pdf (for LaTeX) and ehr-regime-map.png next to this script.
Run: python regime_map.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- physical constants (SI) ----
E = 1.602176634e-19
KB = 1.380649e-23
ME = 9.1093837015e-31
MI = 39.948 * 1.66053906660e-27   # argon

# ---- gas / plasma parameters ----
T_GAS = 300.0                      # K
TE_EV = 3.0                        # electron temperature, eV
SIGMA_EN = 2e-20                   # e-Ar momentum-transfer cross section, m^2 (few eV)
SIGMA_IN = 1e-18                   # Ar+-Ar (elastic + charge exchange), m^2

VBAR_E = np.sqrt(8 * TE_EV * E / (np.pi * ME))
VBAR_I = np.sqrt(8 * KB * T_GAS / (np.pi * MI))


def hall(B_gauss, p_mtorr, mass, sigma, vbar):
    """Hall parameter H = (qB/m) / (n_n sigma vbar)."""
    B = B_gauss * 1e-4                                # G -> T
    n_n = (p_mtorr * 1e-3 * 133.322) / (KB * T_GAS)   # mTorr -> m^-3
    return (E * B / mass) / (n_n * sigma * vbar)


# ---- grid ----
B = np.logspace(-1, np.log10(5000), 400)      # 0.1 G .. 5 kG
p = np.logspace(0, 3, 400)                    # 1 mTorr .. 1 Torr
BB, PP = np.meshgrid(B, p)
He = hall(BB, PP, ME, SIGMA_EN, VBAR_E)
Hi = hall(BB, PP, MI, SIGMA_IN, VBAR_I)

# ---- figure ----
INK = "#1f2937"; MUTED = "#6b7280"
FILLS = ["#eff6ff", "#bfdbfe", "#60a5fa"]     # sequential: degree of magnetization

fig, ax = plt.subplots(figsize=(6.0, 4.3), dpi=300)

# regimes: 0 unmagnetized, 1 electron-magnetized (helicon), 2 fully magnetized
regime = np.zeros_like(He)
regime[(He >= 1) & (Hi < 1)] = 1
regime[Hi >= 1] = 2
ax.contourf(BB, PP, regime, levels=[-0.5, 0.5, 1.5, 2.5], colors=FILLS)

# boundaries
ax.contour(BB, PP, He, levels=[1], colors=[INK], linewidths=1.4)
ax.contour(BB, PP, Hi, levels=[1], colors=[INK], linewidths=1.4)
ci = ax.contour(BB, PP, Hi, levels=[0.01, 0.1, 10], colors=[MUTED],
                linewidths=0.7, linestyles="dotted")
ax.clabel(ci, fmt=lambda v: rf"$H_i={v:g}$", fontsize=7, colors=MUTED,
          manual=[(53, 100), (212, 40), (797, 1.5)])

# region labels (direct labels: identity never rides on fill color alone)
ax.text(0.35, 500, "UNMAGNETIZED\n$H_e<1$", fontsize=8, color=INK,
        ha="center", va="center")
ax.text(20, 30, "ELECTRON-MAGNETIZED\n(helicon regime)  $H_e>1,\\ H_i<1$",
        fontsize=8.5, color=INK, ha="center", va="center")
ax.text(1800, 2.6, "FULLY\nMAGNETIZED\n$H_i>1$", fontsize=8, color=INK,
        ha="center", va="center")
ax.text(19, 300, "$H_e=1$", fontsize=8, color=INK, rotation=45)
ax.text(2400, 30, "$H_i=1$", fontsize=8, color=INK, rotation=45)

# operating points: A = baseline, B = ion-magnetization onset (H_i = 1), same pressure
B_ONSET = 100.0 / hall(100, 10, MI, SIGMA_IN, VBAR_I)
ax.annotate("", xy=(1000, 10), xytext=(115, 10),
            arrowprops=dict(arrowstyle="-|>", color="#b91c1c", lw=1.3))
ax.plot([100], [10], marker="o", ms=7, mfc="#b91c1c", mec="white", mew=1.2, zorder=5)
ax.plot([B_ONSET], [10], marker="s", ms=6.5, mfc="white", mec="#b91c1c", mew=1.4, zorder=5)
ax.text(100, 13.5, "A", fontsize=9, fontweight="bold", color="#b91c1c", ha="center", va="bottom")
ax.text(B_ONSET, 13.5, "B", fontsize=9, fontweight="bold", color="#b91c1c", ha="center", va="bottom")
ax.text(12, 2.05,
        "A: baseline (100 G, 10 mTorr)\n"
        rf"B: ion-magnetization onset ($H_i=1$, $\approx${B_ONSET:.0f} G)"
        "\nred arrow: proposed $B$ scan",
        fontsize=7.5, color="#b91c1c", ha="center", va="center")

# the map explains itself: definition + scaling of the contoured quantity
ax.text(0.85, 4.2,
        r"$H_s=\dfrac{\omega_{c,s}}{\nu_{s,n}} \propto \dfrac{B}{p}$",
        fontsize=9, color=INK, ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.45", fc="white", ec=MUTED, lw=0.8, alpha=0.9))

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(B[0], B[-1]); ax.set_ylim(p[0], p[-1])
ax.set_xlabel("Magnetic field  $B$  [G]", fontsize=10)
ax.set_ylabel("Neutral pressure  $p$  [mTorr]", fontsize=10)
ax.tick_params(labelsize=8.5, colors=INK)
for s in ax.spines.values():
    s.set_color(MUTED); s.set_linewidth(0.8)
ax.set_title("EHR operating regimes in argon ($T_e=3$ eV, $T_{gas}=300$ K)",
             fontsize=10, color=INK)

fig.tight_layout()
out = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(out, "ehr-regime-map.pdf"))
fig.savefig(os.path.join(out, "ehr-regime-map.png"))
print("H_i at (100 G, 10 mTorr):", round(float(hall(100, 10, MI, SIGMA_IN, VBAR_I)), 3))
print("H_e at (100 G, 10 mTorr):", round(float(hall(100, 10, ME, SIGMA_EN, VBAR_E)), 1))
print("B where H_i=1 at 10 mTorr [G]:",
      round(float(100 / hall(100, 10, MI, SIGMA_IN, VBAR_I)), 0))
print("saved ehr-regime-map.pdf / .png")
