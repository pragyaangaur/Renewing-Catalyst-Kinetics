"""Regenerate every result and figure in this repository, then write a manifest.

Usage, from the repository root:

    python3 run_all.py            the current study, about an hour on a laptop
    python3 run_all.py --all      also the superseded first pass, about 20 minutes more
    python3 run_all.py --manifest-only

The steps run in dependency order. Each is a separate process, so a failure stops
the run at the step that failed. The manifest, results/MANIFEST.json, records the
SHA-256 checksum of every public file under results/, figures/ and paper/figures/, the
Python and library versions, the platform, and the wall time of each step. Two runs
on different machines should give identical checksums for every deterministic
file. The stochastic files are seeded, so they also reproduce exactly with the same
library versions, and within their stated errors otherwise.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable

CURRENT = [
    ["src/structure_theory.py"],
    ["src/site_loss_check.py"],
    ["src/ensemble_audit.py", "--checks"],
    ["src/ensemble_audit.py", "--renewal"],
    ["src/research_controls.py"],
    ["src/cluster_expansion.py"],
    ["src/lattice_validation.py"],
    ["src/firstpass_postmortem.py"],
    ["src/structure_figures.py"],
    ["src/audit_figures.py"],
    ["src/paper_figures.py"],
    ["src/paper_numbers.py"],
]

# The superseded first pass, kept so that its withdrawn claims can be checked.
# These scripts write relative to src/, so they run from there.
FIRST_PASS = [
    ["thermo.py"], ["validate.py"], ["kmc_critical.py"], ["kmc_geometry.py"],
    ["run_scans.py", "all"], ["make_figures.py"],
]


def run(step, cwd):
    start = time.time()
    print(f"\n=== {' '.join(step)}", flush=True)
    subprocess.run([PY, *step], cwd=cwd, check=True)
    return round(time.time() - start, 1)


def versions():
    out = {"python": sys.version.split()[0]}
    for mod in ("numpy", "scipy", "matplotlib", "networkx", "mpmath", "pytest"):
        try:
            out[mod] = __import__(mod).__version__
        except ImportError:
            out[mod] = None
    return out


def ignored(paths):
    """The subset of paths that .gitignore excludes from the public record."""
    try:
        out = subprocess.run(["git", "check-ignore", "--stdin"], cwd=ROOT, text=True,
                             input="\n".join(paths), capture_output=True).stdout
    except FileNotFoundError:
        return set()
    return set(out.split())


def manifest(timings):
    files = {}
    candidates = [str(p.relative_to(ROOT)) for folder in ("results", "figures", "paper/figures")
                  for p in sorted((ROOT / folder).rglob("*")) if p.is_file()]
    skip = ignored(candidates)
    for rel in candidates:
        path = ROOT / rel
        if rel not in skip and path.name not in ("MANIFEST.json", ".DS_Store"):
            files[rel] = dict(
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                bytes=path.stat().st_size)
    data = dict(generated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                platform=platform.platform(), machine=platform.machine(),
                versions=versions(), step_seconds=timings, files=files)
    (ROOT / "results" / "MANIFEST.json").write_text(json.dumps(data, indent=1))
    print(f"\nwrote results/MANIFEST.json covering {len(files)} files")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="also rerun the superseded first pass")
    ap.add_argument("--manifest-only", action="store_true")
    args = ap.parse_args()
    timings = {}
    if not args.manifest_only:
        for step in CURRENT:
            timings[" ".join(step)] = run(step, ROOT)
        if args.all:
            for step in FIRST_PASS:
                timings["src/" + " ".join(step)] = run(step, ROOT / "src")
        timings["pytest"] = run(["-m", "pytest", "-q", "tests"], ROOT)
    manifest(timings)


if __name__ == "__main__":
    main()
