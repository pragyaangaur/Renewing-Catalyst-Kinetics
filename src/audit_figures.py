"""Regenerate follow-up research figures from saved numerical results."""

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/renewal-kinetics-matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "audit"
FIGURES = ROOT / "figures"


def main():
    controls = json.loads((RESULTS / "controls.json").read_text())
    with (RESULTS / "renewal.csv").open() as handle:
        renewal = list(csv.DictReader(handle))
    with (RESULTS / "motifs.csv").open() as handle:
        motifs = list(csv.DictReader(handle))
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150})

    figure, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout="constrained")
    selected = [row for row in motifs if float(row["temperature"]) == 745
                and float(row["diffusion"]) == 0 and row["motif"] in
                ("pair", "path3", "path4", "star4", "square")]
    axes[0].bar([row["motif"] for row in selected], [float(row["tof"]) for row in selected], color="#237c91")
    axes[0].set(title="Exact rates on isolated oxide patches", ylabel="SO3 / exposed site / second")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].text(0, 0.12, "exact zero", ha="center", fontsize=8)
    ratios = [row["ratio"] for row in controls["joint_uncertainty"]]
    axes[1].scatter(range(1, 101), ratios, s=15, color="#237c91", label="Original mechanism: 100 rate sets")
    axes[1].axhline(1, color="black", lw=1)
    redox = [row["ratio"] for row in controls["redox_shape_sweep"]]
    axes[1].axhspan(min(redox), max(redox), color="#e29942", alpha=0.4, label="Redox control range")
    axes[1].set(yscale="log", xlabel="Seeded parameter draw", ylabel="Branched / chain turnover",
                title="Shape ranking is not a unique diagnostic")
    axes[1].legend(fontsize=7)
    for name in ("path3", "path4", "star4", "square"):
        selected = [row for row in motifs if float(row["temperature"]) == 745 and row["motif"] == name]
        axes[2].plot([float(row["diffusion"]) for row in selected],
                     [float(row["tof"]) for row in selected], "o-", label=name)
    axes[2].set(xscale="symlog", xlabel="Hop rate per eligible bond (1/s)",
                ylabel="SO3 / exposed site / second", title="Diffusion changes the shape advantage")
    axes[2].legend(fontsize=8)
    figure.suptitle("Mechanism audit | illustrative kinetics, not calibrated predictions")
    figure.savefig(FIGURES / "audit_mechanism.png")
    plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout="constrained")
    for name in ("path3", "path4", "star4", "square"):
        selected = [row for row in renewal if row["motif"] == name
                    and float(row["barrier"]) == 1.45 and float(row["oxidation"]) == 39120]
        axes[0].loglog([float(row["renewal"]) for row in selected],
                       [float(row["tof"]) for row in selected], "o-", ms=3, label=name)
        axes[1].loglog([float(row["tof"]) for row in selected],
                       [float(row["copper_yield"]) for row in selected], "o-", ms=3, label=name)
    axes[0].set(xlabel="Patch reset frequency (1/s)", ylabel="SO3 / exposed site / second",
                title="Fast reoxidation: renewal raises throughput")
    axes[1].set(xlabel="SO3 / exposed site / second", ylabel="SO3 molecules / Cu atom exposed",
                title="Higher throughput can spend more copper")
    axes[0].legend(fontsize=8)
    ceilings = {entry["barrier"]: entry["ceiling"] for entry in controls["renewal_identity"]}
    for barrier, color in ((0.95, "#e29942"), (1.45, "#237c91")):
        selected = [row for row in renewal if float(row["barrier"]) == barrier]
        axes[2].scatter([float(row["sulfate_coverage"]) for row in selected],
                        [float(row["copper_yield"]) / ceilings[barrier] for row in selected],
                        s=10, alpha=0.6, color=color, label=f"Sulphation barrier {barrier} eV")
    axes[2].plot([0, 1], [0, 1], "k--", lw=1, label="Exact balance identity")
    axes[2].set(xlabel="Mean sulphated fraction of patch", ylabel="Yield / branching ceiling",
                title="364 cases obey one material-efficiency law", xlim=(-0.03, 1.03), ylim=(-0.03, 1.03))
    axes[2].legend(fontsize=7)
    figure.suptitle("Renewal model | irreversible sulphation; each reset exposes one Cu layer")
    figure.savefig(FIGURES / "audit_renewal.png")
    plt.close(figure)
    print("Saved audit_mechanism.png and audit_renewal.png")


if __name__ == "__main__":
    main()
