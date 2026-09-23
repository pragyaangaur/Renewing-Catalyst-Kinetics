"""Figures for the paper, drawn only from saved results.

Nothing is recomputed here. Every panel reads a file under results/, so the figures
always match the data that is archived alongside them. Output goes to
paper/figures/ as vector PDF, with PNG copies for quick viewing.
"""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "paper" / "figures"

plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.spines.top": False,
    "axes.spines.right": False, "savefig.bbox": "tight", "savefig.dpi": 300,
    "font.family": "serif", "mathtext.fontset": "dejavuserif", "lines.linewidth": 1.3,
})

# Okabe and Ito colours, distinguishable under the common colour vision deficiencies
BLUE, ORANGE, GREEN, RED, PURPLE, SKY, GREY = ("#0072B2", "#E69F00", "#009E73", "#D55E00",
                                               "#CC79A7", "#56B4E9", "#555555")
LATTICE_COLOURS = {"honeycomb": ORANGE, "square": BLUE, "triangular": GREEN}


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", metadata={"CreationDate": None})
    fig.savefig(FIG / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"  wrote paper/figures/{name}.pdf")


def label(ax, text):
    ax.text(0.0, 1.06, text, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="right")


def fig_result_a():
    rows = json.loads((RES / "structure" / "result_a.json").read_text())
    reach = json.loads((RES / "structure" / "result_a_reachability.json").read_text())
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw=dict(width_ratios=[1.35, 1]))
    floor = 1e-25
    rng = np.random.default_rng(0)
    for n in sorted({r["n"] for r in rows}):
        vals = np.array([r["tof"] for r in rows if r["n"] == n])
        x = n + rng.uniform(-0.18, 0.18, len(vals))
        live = n >= 3
        a.scatter(x, np.maximum(vals, floor), s=9, lw=0, alpha=0.8,
                  color=BLUE if live else RED,
                  label=("strictly positive" if live else "exactly zero") if n in (1, 3) else None)
    a.axhline(floor, color=RED, lw=0.8, ls="--")
    a.text(1.5, floor * 30, "exact zero", ha="center", color=RED, fontsize=7)
    a.axvspan(0.5, 2.5, color=RED, alpha=0.07, lw=0)
    a.set_yscale("log")
    a.set_ylim(floor / 10, 1e5)
    a.set_xticks(range(1, 6))
    a.set_xlabel("sites in the isolated active patch")
    a.set_ylabel(r"stationary product rate per site (s$^{-1}$)")
    a.legend(loc="lower right", frameon=False)
    label(a, "a")

    sizes = sorted({r["n"] for r in reach})
    graphs = [sum(r["n"] == n for r in reach) for n in sizes]
    states = [max(r["states_with_diffusion"] for r in reach if r["n"] == n) for n in sizes]
    b.bar(sizes, graphs, color=[RED if n < 3 else BLUE for n in sizes], width=0.7)
    for n, g in zip(sizes, graphs):
        b.text(n, g * 1.15, str(g), ha="center", fontsize=7)
    b.set_yscale("log")
    b.set_ylim(0.7, 4000)
    b.set_xlabel("sites in the patch")
    b.set_ylabel("connected graphs checked")
    b2 = b.twinx()
    b2.plot(sizes, states, "o-", color=GREY, ms=3, lw=1)
    b2.set_yscale("log")
    b2.set_ylabel("reachable configurations", color=GREY)
    b2.tick_params(axis="y", colors=GREY)
    b2.spines["right"].set_visible(True)
    label(b, "b")
    fig.tight_layout()
    save(fig, "fig1_result_a")


def fig_dead_fraction():
    pop = json.loads((RES / "cluster" / "population.json").read_text())
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for lattice, d in pop.items():
        c = LATTICE_COLOURS[lattice]
        p = np.array(d["p"])
        a.plot(1 - p, d["dead_closed_form"], color=c, label=f"{lattice}, z = {d['z']}")
        a.plot(1 - p, d["dead_lattice_oxygen"], color=c, ls=":", lw=1.0)
        mean = np.array([s["mean"][0] + s["mean"][1] for s in d["sampled"]])
        a.plot(1 - p[::3], mean[::3], "o", color=c, ms=3.2, mfc="white", mew=0.9)
        # residual panel, in units of the sampling error with the counting floor
        err = np.array([np.hypot(s["sem"][0], s["sem"][1]) for s in d["sampled"]])
        floor = np.sqrt(2.0 * np.array(d["dead_closed_form"]) / (8 * 512 * 512 * p))
        b.plot(1 - p, (mean - np.array(d["dead_closed_form"])) / np.maximum(err, floor), "o-",
               color=c, ms=2.5, lw=0.7, label=lattice)
    a.plot([], [], color=GREY, ls=":", label="lattice oxygen control")
    a.plot([], [], "o", color=GREY, mfc="white", ms=3.2, label="sampled, 8 lattices of 512$^2$")
    a.set_xlabel(r"deactivated fraction $1-p$")
    a.set_ylabel("dead share of surviving sites")
    a.set_xlim(0, 1)
    a.set_ylim(0, 1.02)
    a.legend(frameon=False, loc="upper left")
    label(a, "a")
    b.axhspan(-2, 2, color=GREY, alpha=0.12, lw=0)
    b.axhline(0, color=GREY, lw=0.8)
    b.set_xlabel(r"deactivated fraction $1-p$")
    b.set_ylabel("(sampled $-$ exact) / s.e.")
    b.set_ylim(-4.5, 4.5)
    b.legend(frameon=False, loc="lower left", ncol=3)
    label(b, "b")
    fig.tight_layout()
    save(fig, "fig2_dead_fraction")


def fig_expansion():
    exp = json.loads((RES / "cluster" / "expansion.json").read_text())
    post = json.loads((RES / "firstpass_postmortem.json").read_text())
    lat = json.loads((RES / "lattice_check.json").read_text())
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.3))

    d = exp["illustrative_745K"]
    p = np.array(d["p"])
    sizes = sorted(d["contribution_by_size"], key=int)
    bottom = np.zeros_like(p)
    shades = ["#c6dbef", "#6baed6", "#2171b5", "#08306b"]
    for s, col in zip(sizes, shades):
        v = np.array(d["contribution_by_size"][s])
        a.fill_between(1 - p, bottom, bottom + v, color=col, lw=0, label=f"{s}-site components")
        bottom += v
    for row in post:
        x = row["coverage"]
        a.plot(x, row["exact_stationary_tof"], "D", color=ORANGE, ms=5, mec="black", mew=0.5,
               zorder=5)
        a.errorbar(x, row["long_run"]["tof"], yerr=2 * row["long_run"]["sem"], fmt="s",
                   color=GREEN, ms=4, mec="black", mew=0.5, capsize=2, zorder=6)
        a.plot(x, max(row["first_pass_tof"], 3e-4), "v" if row["first_pass_tof"] == 0 else "x",
               color=RED, ms=6, zorder=7)
    a.plot([], [], "D", color=ORANGE, mec="black", mew=0.5, ms=5,
           label="exact, first pass lattice")
    a.plot([], [], "s", color=GREEN, mec="black", mew=0.5, ms=4, label="long KMC run, same lattice")
    a.plot([], [], "x", color=RED, ms=6, label="first pass KMC")
    a.plot([], [], "v", color=RED, ms=6, label="first pass KMC, zero (at floor)")
    a.set_yscale("log")
    a.set_ylim(3e-4, 1.5)
    a.set_xlabel(r"randomly deactivated fraction $1-p$")
    a.set_ylabel(r"product rate per lattice site (s$^{-1}$)")
    handles, labels = a.get_legend_handles_labels()
    label(a, "a")

    rows = lat["rows"]
    colours = {0.1: ORANGE, 0.2: BLUE, 0.3: GREEN}
    for r in rows:
        x = r["p"] + 0.008 * (r["realisation"] - 3.5)
        b.errorbar(x, 100 * (r["kmc"] / r["exact"] - 1), yerr=200 * r["sem"] / r["exact"],
                   fmt="o", ms=3, color=colours[r["p"]], capsize=1.5, lw=0.8)
    b.axhline(0, color=GREY, lw=0.8)
    b.set_xticks(list(colours))
    b.set_xlim(0.06, 0.34)
    b.set_xlabel(r"active fraction $p$, eight lattices each")
    b.set_ylabel("KMC minus exact (%, bars 2 s.e.)")
    s = lat["summary"]
    b.text(0.98, 0.97, f"{s['runs']} lattices\n$\\chi^2$ = {s['chi2']:.1f} on {s['dof']} dof",
           transform=b.transAxes, ha="right", va="top", fontsize=7)
    label(b, "b")
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=6.5)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    save(fig, "fig3_expansion_and_validation")


def fig_result_b():
    rows = json.loads((RES / "structure" / "result_b.json").read_text())
    with (RES / "audit" / "renewal.csv").open() as fh:
        renewal = list(csv.DictReader(fh))
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.2, 2.5))
    on = [r for r in rows if r["regen_on"]]
    off = [r for r in rows if not r["regen_on"]]
    for group, col, mk, lab in ((off, GREEN, "o", "renewal only"),
                                (on, PURPLE, "s", "with regeneration")):
        a.loglog([r["predicted"] for r in group], [r["actual"] for r in group], mk, color=col,
                 ms=3, alpha=0.8, lw=0, label=lab)
    lo = min(r["actual"] for r in rows if r["actual"] > 0) * 0.3
    hi = max(r["actual"] for r in rows) * 3
    a.plot([lo, hi], [lo, hi], color=GREY, lw=0.8, ls="--")
    a.set_xlabel(r"$Y_\mathrm{Cu}$ from the identity")
    a.set_ylabel(r"$Y_\mathrm{Cu}$ from the stationary solve")
    a.legend(frameon=False, loc="upper left")
    label(a, "a")

    for group, col, mk in ((off, GREEN, "o"), (on, PURPLE, "s")):
        b.semilogy([r["f_sulfate"] for r in group], [r["actual"] / r["ratio"] for r in group], mk,
                   color=col, ms=3, alpha=0.8, lw=0)
    b.axhline(1.0, color=RED, lw=1.1, ls="--", label=r"ceiling $Y_\mathrm{Cu}=k_d/k_s$")
    b.legend(frameon=False, loc="upper left", fontsize=6.5)
    b.set_xlabel(r"mean sulphated fraction $f_\mathrm{sulphate}$")
    b.set_ylabel(r"$Y_\mathrm{Cu}\,/\,(k_d/k_s)$")
    label(b, "b")

    for name, col in (("path3", ORANGE), ("path4", SKY), ("star4", PURPLE), ("square", BLUE)):
        sel = [r for r in renewal if r["motif"] == name and float(r["barrier"]) == 1.45
               and float(r["oxidation"]) == 39120]
        c.loglog([float(r["tof"]) for r in sel], [float(r["copper_yield"]) for r in sel], "o-",
                 color=col, ms=2.5, lw=0.9, label=name)
    c.axhline(107.0, color=RED, lw=1.0, ls="--", label="ceiling")
    c.set_xlabel(r"product per exposed site (s$^{-1}$)")
    c.set_ylabel("product per Cu atom exposed")
    c.legend(frameon=False, loc="lower left", ncol=2, fontsize=6, columnspacing=0.8)
    label(c, "c")
    fig.tight_layout()
    save(fig, "fig4_result_b")


def fig_controls():
    controls = json.loads((RES / "audit" / "controls.json").read_text())
    with (RES / "audit" / "motifs.csv").open() as fh:
        motifs = list(csv.DictReader(fh))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.6))
    ratios = np.array([r["ratio"] for r in controls["joint_uncertainty"]])
    a.scatter(np.arange(1, len(ratios) + 1), ratios, s=8, color=BLUE, lw=0,
              label="dual-site mechanism, 100 joint rate draws")
    redox = [r["ratio"] for r in controls["redox_shape_sweep"]]
    a.axhspan(min(redox), max(redox), color=ORANGE, alpha=0.35, lw=0,
              label="lattice oxygen control, all refill ratios")
    a.axhline(1, color=GREY, lw=0.8)
    a.set_yscale("log")
    a.set_xlabel("parameter draw")
    a.set_ylabel("star / chain rate, four sites")
    a.legend(frameon=False, loc="upper left")
    label(a, "a")
    for name, col in (("path3", ORANGE), ("path4", SKY), ("star4", PURPLE), ("square", BLUE)):
        sel = [r for r in motifs if float(r["temperature"]) == 745 and r["motif"] == name]
        b.plot([max(float(r["diffusion"]), 1) for r in sel], [float(r["tof"]) for r in sel], "o-",
               color=col, ms=3, label=name)
    b.set_xscale("log")
    b.set_xlabel(r"hop rate per eligible bond (s$^{-1}$; 0 drawn at 1)")
    b.set_ylabel(r"product per site (s$^{-1}$)")
    b.legend(frameon=False)
    label(b, "b")
    fig.tight_layout()
    save(fig, "fig5_controls")


def main():
    fig_result_a()
    fig_dead_fraction()
    fig_expansion()
    fig_result_b()
    fig_controls()


if __name__ == "__main__":
    main()
