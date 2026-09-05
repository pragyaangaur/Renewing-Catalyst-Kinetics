"""Generate every figure for the study from the saved result archives."""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

import thermo
from reactor import required_tof
from kmc_geometry import make_random, make_islands, pair_stats

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "..", "results")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 160, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.bbox": "tight",
})

C_WIN = "#f2c14e"
C_A, C_B, C_C, C_D = "#1b6ca8", "#c1443c", "#2e8b57", "#7d5ba6"
WINDOW = (673.0, 823.0)


def _window(ax, label=True):
    ax.axvspan(*WINDOW, color=C_WIN, alpha=0.22, lw=0, zorder=0)
    ax.axvline(745, color="#8a6d1f", ls=":", lw=1.2, zorder=1)
    if label:
        ax.text(745, ax.get_ylim()[1], " 745 K", ha="left", va="top",
                fontsize=7.5, color="#8a6d1f")


# ------------------------------------------------------------------- figure 1

def fig1_thermo():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 3.9))

    Ts = np.linspace(600, 1000, 500)
    ax1.semilogy(Ts, thermo.equilibrium_so3_ppm(Ts), color=C_A, lw=2,
                 label="equilibrium SO$_3$ over CuSO$_4$/CuO (bulk)")
    for so3, lab, col in [(900, "900 ppm SO$_3$  (1000 ppm SO$_2$, 90% conv)", C_B),
                          (89, "89 ppm SO$_3$  (100 ppm SO$_2$ duty)", C_C)]:
        ax1.axhline(so3, color=col, ls="--", lw=1.3)
        T_on = thermo.sulfation_onset_temperature(so3)
        ax1.plot([T_on], [so3], "o", color=col, ms=6, zorder=5)
        ax1.annotate(f"{T_on:.0f} K", (T_on, so3), textcoords="offset points",
                     xytext=(6, -12), fontsize=8, color=col, fontweight="bold")
    ax1.set_ylim(1e-2, 1e6)
    _window(ax1, label=False)
    ax1.axvline(745, color="#8a6d1f", ls=":", lw=1.2)
    ax1.set_xlabel("temperature (K)")
    ax1.set_ylabel("SO$_3$ partial pressure (ppmv)")
    ax1.set_title("Sulphation is thermodynamically favoured\nacross most of the range of interest",
                  fontsize=9.5)
    ax1.legend(fontsize=7, loc="upper left")
    ax1.text(690, 3e-2, "673 to 823 K", fontsize=7.5, color="#8a6d1f", ha="center")

    ds = np.linspace(0, 90, 200)
    onset = np.array([thermo.surface_corrected_onset(900.0, d) for d in ds])
    ax2.plot(ds, onset, color=C_A, lw=2)
    ax2.axhspan(*WINDOW, color=C_WIN, alpha=0.22, lw=0)
    ax2.axhline(745, color="#8a6d1f", ls=":", lw=1.2)
    d745 = np.interp(-745.0, -onset, ds)
    ax2.plot([d745], [745], "o", color=C_B, ms=7, zorder=5)
    ax2.annotate(f"{d745:.0f} kJ/mol puts the\nboundary at 745 K",
                 (d745, 745), textcoords="offset points", xytext=(10, 26),
                 fontsize=8, color=C_B, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=C_B, lw=1))
    ax2.set_xlabel("surface sulphate destabilisation vs bulk CuSO$_4$ (kJ/mol)")
    ax2.set_ylabel("sulphation onset temperature (K)")
    ax2.set_title("A surface sulphate is less stable than bulk CuSO$_4$,\n"
                  "which moves the boundary into the range of interest", fontsize=9.5)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_thermodynamic_boundary.png"))
    plt.close(fig)
    print("  fig1_thermodynamic_boundary.png")


# ------------------------------------------------------------------- figure 2

def fig2_headroom():
    d = np.load(os.path.join(RES, "s2_headroom.npz"))
    f, tof, req = d["f_CuSO4"], d["tof"], float(d["req_tof"][0])

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogy(f * 100, np.maximum(tof, 1e-6), color=C_A, lw=2,
                label="mean field turnover frequency")
    ax.axhline(req, color=C_B, lw=1.8, ls="--",
               label=f"reactor requirement, {req:.2f} s$^{{-1}}$")
    ax.fill_between(f * 100, req, 1e5, color=C_C, alpha=0.10, lw=0)
    ax.text(8, 2.2e3, "surplus activity", color=C_C, fontsize=9, fontweight="bold")

    ok = tof >= req
    fc = f[ok].max() * 100
    ax.axvline(fc, color="#555", ls=":", lw=1.2)
    ax.annotate(f"mean field says {fc:.1f}%\nsulphation is tolerable",
                (fc, req * 30), textcoords="offset points", xytext=(-118, 10),
                fontsize=8, color="#333",
                arrowprops=dict(arrowstyle="->", color="#555", lw=1))
    ax.set_xlabel("steady state sulphated fraction of the surface (%)")
    ax.set_ylabel("turnover frequency (per site per second)")
    ax.set_ylim(1e-3, 3e4)
    ax.set_title("The reactor needs about one turnover per site per second,\n"
                 "and a clean oxide surface offers roughly nine thousand", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_activity_headroom.png"))
    plt.close(fig)
    print("  fig2_activity_headroom.png")


# ------------------------------------------------------------------- figure 3

def fig3_kmc_vs_mf():
    d = np.load(os.path.join(RES, "kmc_critical.npz"))
    mf_f, kmc_f = d["mf_f"], d["kmc_f"]
    mf_t, kmc_t = d["mf_tof"], d["kmc_tof"]
    iso, req = d["frac_isolated"], float(d["req"][0])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    ax.semilogy(mf_f * 100, np.maximum(mf_t, 1e-5), "o-", color=C_A, lw=2, ms=5,
                label="mean field")
    ax.semilogy(kmc_f * 100, np.maximum(kmc_t, 1e-5), "s-", color=C_B, lw=2, ms=5,
                label="exact lattice KMC")
    ax.axhline(req, color="#444", ls="--", lw=1.5,
               label=f"reactor requirement, {req:.2f} s$^{{-1}}$")
    ax.set_xlabel("sulphated fraction of the surface (%)")
    ax.set_ylabel("turnover frequency (per site per second)")
    ax.set_ylim(1e-4, 5e4)
    ax.set_title("Mean field overstates surviving activity\nby up to two orders of magnitude",
                 fontsize=10)
    ax.legend(fontsize=8, loc="lower left")

    for xf, yk, ym in zip(kmc_f * 100, kmc_t, mf_t):
        if ym > 0 and yk > 0 and ym / yk > 5:
            ax.annotate("", xy=(xf, yk), xytext=(xf, ym),
                        arrowprops=dict(arrowstyle="<->", color="#999", lw=0.9))

    ax2.plot(kmc_f * 100, iso * 100, "o-", color=C_D, lw=2, ms=5)
    ax2.set_xlabel("sulphated fraction of the surface (%)")
    ax2.set_ylabel("surviving oxide sites with no oxide neighbour (%)")
    ax2.set_title("Because the surviving oxide sites\nstop touching one another", fontsize=10)
    ax2.axhline(50, color="#999", ls=":", lw=1)
    ax2.text(30, 53, "half the remaining sites are stranded", fontsize=7.5, color="#666")
    ax2.set_ylim(-3, 105)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_kmc_vs_meanfield.png"))
    plt.close(fig)
    print("  fig3_kmc_vs_meanfield.png")


# ------------------------------------------------------------------- figure 4

def fig4_geometry():
    d = np.load(os.path.join(RES, "geometry.npz"))
    cov, tof, conn = d["cov"], d["tof"], d["connected"]
    labels = [str(x) for x in d["labels"]]
    req = float(d["req"][0])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2))
    cols = [C_B, C_A, C_C, C_D]
    marks = ["o", "s", "^", "D"]
    for j, lab in enumerate(labels):
        ax.semilogy(cov * 100, np.maximum(tof[:, j], 1e-4), marks[j] + "-",
                    color=cols[j], lw=2, ms=5, label=lab)
    ax.axhline(req, color="#444", ls="--", lw=1.5,
               label=f"reactor requirement, {req:.2f} s$^{{-1}}$")
    ax.set_xlabel("sulphated fraction of the surface (%)")
    ax.set_ylabel("turnover frequency (per site per second)")
    ax.set_ylim(1e-3, 3e4)
    ax.set_title("At the same amount of sulphate, the pattern it forms\n"
                 "decides whether the catalyst still works", fontsize=10)
    ax.legend(fontsize=7.5, loc="lower left")
    ax.annotate("random sulphation is dead here",
                (90, 1e-3), textcoords="offset points", xytext=(-40, 34),
                fontsize=8, color=C_B, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=C_B, lw=1))

    for j, lab in enumerate(labels):
        ax2.plot(cov * 100, conn[:, j] * 100, marks[j] + "-", color=cols[j],
                 lw=2, ms=5, label=lab)
    ax2.set_xlabel("sulphated fraction of the surface (%)")
    ax2.set_ylabel("oxide sites retaining an oxide neighbour (%)")
    ax2.set_title("Islands keep the surviving oxide connected", fontsize=10)
    ax2.legend(fontsize=7.5, loc="lower left")
    ax2.set_ylim(0, 105)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_sulphation_geometry.png"))
    plt.close(fig)
    print("  fig4_sulphation_geometry.png")


# ------------------------------------------------------------------- figure 5

def fig5_snapshots():
    L, cov = 64, 0.90
    rng = np.random.default_rng(4)
    confs = [("random", make_random(L, cov, rng)),
             ("islands, 64 seeds", make_islands(L, cov, 64, rng)),
             ("islands, 16 seeds", make_islands(L, cov, 16, rng)),
             ("islands, 4 seeds", make_islands(L, cov, 4, rng))]
    d = np.load(os.path.join(RES, "geometry.npz"))
    ci = int(np.argmin(np.abs(d["cov"] - cov)))
    tofs = d["tof"][ci]

    cmap = ListedColormap(["#d8d8d8", "#1b6ca8", "#c9a227"])
    fig, axes = plt.subplots(1, 4, figsize=(11.4, 3.5))
    for ax, (lab, ph), tf in zip(axes, confs, tofs):
        ax.imshow(ph.reshape(L, L), cmap=cmap, norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3),
                  interpolation="nearest")
        c, _ = pair_stats(ph, L)
        ax.set_title(f"{lab}\nTOF {tf:,.0f} s$^{{-1}}$   connected {c*100:.0f}%",
                     fontsize=8.5)
        ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    fig.suptitle("All four surfaces are 90 percent sulphated. "
                 "Blue is the surviving copper oxide, gold is sulphate.",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_lattice_snapshots.png"))
    plt.close(fig)
    print("  fig5_lattice_snapshots.png")


# ------------------------------------------------------------------- figure 6

def fig6_context():
    s1 = np.load(os.path.join(RES, "s1_temperature.npz"))
    s7 = np.load(os.path.join(RES, "s7_axial.npz"))

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.0))
    ax.semilogy(s1["T"], np.maximum(s1["tof"], 1e-4), color=C_A, lw=2,
                label="available turnover frequency")
    ax.semilogy(s1["T"], s1["req_tof"], color=C_B, ls="--", lw=1.8,
                label="reactor requirement")
    ax.set_ylim(1e-2, 1e5)
    _window(ax)
    ax.set_xlabel("temperature (K)")
    ax.set_ylabel("turnover frequency (per site per second)")
    ax.set_title("Available activity against what the duty needs", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")

    ax3 = ax2.twinx()
    for T, col in [(700, C_C), (745, C_A), (800, C_D)]:
        ax2.plot(s7[f"z_{T}"], s7[f"fSO4_{T}"] * 100, color=col, lw=2, label=f"{T} K")
        ax3.plot(s7[f"z_{T}"], s7[f"ySO2_{T}"] * 1e6, color=col, lw=1, ls=":")
    ax2.set_xlabel("residence time along the bed (s)")
    ax2.set_ylabel("sulphated fraction (%)  [solid]")
    ax3.set_ylabel("SO$_2$ remaining (ppm)  [dotted]")
    ax3.grid(False)
    ax2.set_title("The bed sulphates at its front, not its back,\n"
                  "because sulphate forms from adsorbed SO$_3$", fontsize=10)
    ax2.legend(fontsize=8, title="temperature", title_fontsize=8)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_temperature_and_axial.png"))
    plt.close(fig)
    print("  fig6_temperature_and_axial.png")


if __name__ == "__main__":
    print("figures:")
    fig1_thermo()
    fig2_headroom()
    fig3_kmc_vs_mf()
    fig4_geometry()
    fig5_snapshots()
    fig6_context()
    print("done")
