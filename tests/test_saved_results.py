"""The numbers quoted in the paper and docs, read back from the saved results."""

import csv
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "results" / "audit"


def test_exact_motif_rates_quoted_in_the_paper():
    rows = {r["motif"]: float(r["tof"]) for r in csv.DictReader((AUDIT / "motifs.csv").open())
            if float(r["temperature"]) == 745 and float(r["diffusion"]) == 0}
    assert rows["pair"] == 0.0 and rows["single"] == 0.0
    assert rows["path3"] == pytest.approx(1.229, abs=5e-4)
    assert rows["path4"] == pytest.approx(1.855, abs=5e-4)
    assert rows["star4"] == pytest.approx(5.448, abs=5e-4)
    assert rows["square"] == pytest.approx(3.500, abs=5e-4)
    assert rows["star4"] / rows["path4"] == pytest.approx(2.94, abs=5e-3)


def test_stochastic_checks_all_passed():
    checks = json.loads((AUDIT / "checks.json").read_text())
    assert all(c["passed"] for c in checks)
    controls = json.loads((AUDIT / "controls.json").read_text())
    assert controls["renewal_ssa_check"]["exact"] == pytest.approx(0.05733, abs=5e-5)


def test_uncertainty_sweep_range_quoted_in_the_paper():
    controls = json.loads((AUDIT / "controls.json").read_text())
    ratios = [r["ratio"] for r in controls["joint_uncertainty"]]
    assert min(ratios) == pytest.approx(0.783, abs=5e-4)
    assert max(ratios) == pytest.approx(272.2, abs=0.05)
    assert sum(r > 1 for r in ratios) == 88
    redox = [r["ratio"] for r in controls["redox_shape_sweep"]]
    assert min(redox) == pytest.approx(0.712, abs=5e-4)
    assert max(redox) == pytest.approx(0.994, abs=5e-4)


def test_renewal_ceilings_and_identity():
    controls = json.loads((AUDIT / "controls.json").read_text())
    ceilings = {e["barrier"]: e for e in controls["renewal_identity"]}
    assert ceilings[0.95]["ceiling"] == pytest.approx(0.04436, abs=5e-6)
    assert ceilings[1.45]["ceiling"] == pytest.approx(107.0, abs=0.05)
    assert sum(e["points"] for e in ceilings.values()) == 364


def test_geometry_sampling_quoted_in_the_paper():
    controls = json.loads((AUDIT / "controls.json").read_text())
    big = {g["coordination"]: g for g in controls["geometry"] if g["side"] == 128}
    assert big[4]["mean"][2] == pytest.approx(0.13213, abs=5e-5)
    assert big[6]["mean"][2] == pytest.approx(0.21253, abs=5e-5)


def test_lattice_kmc_validation_passed():
    data = json.loads((ROOT / "results" / "lattice_check.json").read_text())
    assert data["summary"]["all_passed"]
    assert data["summary"]["runs"] == 24


def test_first_pass_postmortem():
    rows = {round(r["coverage"], 2): r for r in
            json.loads((ROOT / "results" / "firstpass_postmortem.json").read_text())}
    ninety, eighty_five = rows[0.9], rows[0.85]
    # the replayed first pass reproduces the numbers it originally reported
    assert ninety["observed_events"] == 0
    assert eighty_five["first_pass_tof"] == pytest.approx(0.1567, abs=5e-5)
    # the surface at 90 percent was never dead
    assert ninety["exact_stationary_tof"] > 0.02
    assert ninety["largest_component"] <= 6 and eighty_five["largest_component"] <= 6
    for row in (ninety, eighty_five):
        assert abs(row["long_run"]["z"]) < 3


def test_paper_numbers_file_is_current():
    import subprocess
    import sys
    target = ROOT / "paper" / "numbers.tex"
    if not target.exists():
        pytest.skip("the paper source is not part of this checkout")
    before = target.read_text()
    subprocess.run([sys.executable, str(ROOT / "src" / "paper_numbers.py")], check=True,
                   capture_output=True)
    assert target.read_text() == before
