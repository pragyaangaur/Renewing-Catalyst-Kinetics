"""Figures for the two rate independent results.

Both panels are drawn from results/structure/, which structure_theory.py writes.
Nothing is recomputed here, so the figure always matches the saved run.
"""

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "structure"
FIG = ROOT / "figures"

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 160, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.bbox": "tight",
})

DEAD = "#c1443c"
LIVE = "#1b6ca8"
ACCENT = "#2e8b57"
MUTED = "#7d5ba6"


def panel_a(ax):
    rows = json.loads((RES / "result_a.json").read_text())
    sizes = sorted({r["n"] for r in rows})

    for n in sizes:
        vals = [r["tof"] for r in rows if r["n"] == n]
        live = n >= 3
        colour = LIVE if live else DEAD
        # exact zeros cannot be drawn on a log axis, so they sit on a floor line
        floor = 1e-14
        y = [max(v, floor) for v in vals]
        x = np.full(len(y), n) + np.random.default_rng(n).uniform(-0.13, 0.13, len(y))
        ax.scatter(x, y, s=16, color=colour, alpha=0.75, lw=0,
                   label=("productive" if live else "exactly zero") if n in (1, 3) else None)

    ax.axhline(1e-14, color=DEAD, lw=1.0, ls="--")
    ax.axvspan(0.5, 2.5, color=DEAD, alpha=0.09, lw=0)
    ax.set_yscale("log")
    ax.set_ylim(3e-15, 1e4)
    ax.set_xticks(sizes)
    ax.set_xlabel("sites in the connected patch")
    ax.set_ylabel("stationary product rate per site (s$^{-1}$)")
    ax.set_title("Result A: one and two site patches are exactly dead\n"
                 "for every rate assignment", fontsize=10)
    ax.text(1.5, 3e-13, "no reachable state\nholds O and SO$_2$ together",
            ha="center", fontsize=7.5, color=DEAD)
    ax.legend(fontsize=8, loc="lower right")


def panel_b(ax):
    rows = json.loads((RES / "result_b.json").read_text())
    on = [r for r in rows if r["regen_on"]]
    off = [r for r in rows if not r["regen_on"]]

    for group, colour, marker, label in [
        (off, ACCENT, "o", "renewal only"),
        (on, MUTED, "s", "with chemical regeneration"),
    ]:
        pred = [r["predicted"] for r in group]
        act = [r["actual"] for r in group]
        ax.loglog(pred, act, marker, color=colour, ms=4.5, alpha=0.7, lw=0, label=label)

    lo = min(r["actual"] for r in rows if r["actual"] > 0) * 0.4
    hi = max(r["actual"] for r in rows) * 2.5
    ax.plot([lo, hi], [lo, hi], color="#444", lw=1.2, ls="--", label="identity")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("predicted $Y_{Cu}$ from the identity")
    ax.set_ylabel("measured $Y_{Cu}$ from the stationary solve")
    ax.set_title("Result B: the copper efficiency identity holds\n"
                 "to 6e-6 relative over 200 randomised cases", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")


def panel_c(ax):
    rows = json.loads((RES / "result_b.json").read_text())
    off = [r for r in rows if not r["regen_on"]]
    on = [r for r in rows if r["regen_on"]]

    for group, colour, marker, label in [
        (off, ACCENT, "o", "renewal only"),
        (on, MUTED, "s", "with chemical regeneration"),
    ]:
        x = [r["f_sulfate"] for r in group]
        y = [r["actual"] / r["ratio"] for r in group]
        ax.semilogy(x, y, marker, color=colour, ms=4.5, alpha=0.7, lw=0, label=label)

    ax.axhline(1.0, color=DEAD, lw=1.6, ls="--")
    ax.text(0.02, 1.35, "ceiling for renewal only: $Y_{Cu} \\leq k_d/k_s$",
            fontsize=8, color=DEAD, fontweight="bold")
    ax.set_xlabel("mean sulphated fraction of the patch")
    ax.set_ylabel("$Y_{Cu}$ divided by the branching ratio $k_d/k_s$")
    ax.set_title("Renewal alone cannot beat the branching ratio.\n"
                 "Chemical regeneration is the only way past it", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")


def main():
    FIG.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.3))
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    fig.tight_layout()
    out = FIG / "structure_results.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
