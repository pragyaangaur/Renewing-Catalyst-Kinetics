"""Check the lattice KMC engine against an exact answer on a whole lattice.

The earlier engine checks used a single three-site patch. This one runs the full
rejection free engine in kmc_core on a 64 by 64 periodic square lattice with a
random fraction of sites deactivated and frozen, which is the setting the first
pass geometry study used.

On such a lattice the exact stationary rate is known. Deactivated sites take part
in no process, so the active components evolve independently, and the lattice rate
is the sum of the exact rates of its components. Every component on the square lattice with six or
fewer sites is isomorphic to one of the eighteen shapes that cluster_expansion.py
has already solved with structure_theory.solve, so the rates are looked up from
results/cluster/expansion.json rather than solved again. To keep every component
exactly solvable, components with more than six sites are also deactivated before
the run. What is compared is
therefore one KMC trajectory against one exact number for the same frozen lattice.

The rates are synthetic and of order one, as in the renewal check of
research_controls. The illustrative 745 K rates are too stiff for a pure Python
engine to sample a small rate on a whole lattice, which is the failure that
docs/CORRECTIONS.md records. The engine is the same in both cases, so a check at
tractable rates tests the event selection, the pool bookkeeping and the clock.

Errors are estimated by batch means over twenty equal batches of the measurement
window. The acceptance test is five standard errors plus one percent, declared
before the runs, and a pooled chi-square over all runs is reported as well.
"""

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import stats
from scipy.ndimage import label

from kmc_core import (KMC, Lattice, CUO, CUSO4, EMPTY, N_SITE_PROC, N_BOND_PROC,
                      P_ADS_SO2, P_DES_SO2, P_DES_SO3, B_ADS_O2, B_DES_O2, B_LH)
from cluster_expansion import SYNTHETIC
from structure_theory import solve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "lattice_check.json"

SIDE = 64
S_MAX = 6
BURN = 500.0
MEASURE = 10_000.0
BATCHES = 20
COVERAGES = (0.10, 0.20, 0.30)
REALISATIONS = 8

KMC_RATES = {("site", p): 0.0 for p in range(N_SITE_PROC)}
KMC_RATES.update({("bond", p): 0.0 for p in range(N_BOND_PROC)})
KMC_RATES.update({("site", P_ADS_SO2): 1.0, ("site", P_DES_SO2): 1.0,
                  ("site", P_DES_SO3): 1.0, ("bond", B_ADS_O2): 1.0,
                  ("bond", B_DES_O2): 0.1, ("bond", B_LH): 1.0})
EXACT_RATES = SYNTHETIC
EXPANSION = ROOT / "results" / "cluster" / "expansion.json"


def frozen_lattice(p, seed):
    """Random active sites, with components larger than S_MAX deactivated."""
    rng = np.random.default_rng(seed)
    active = rng.random((SIDE, SIDE)) < p
    labels, n = label(active, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
    # join labels that touch across the periodic boundary
    parent = list(range(n + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in list(zip(labels[0], labels[-1])) + list(zip(labels[:, 0], labels[:, -1])):
        if a and b:
            parent[find(a)] = find(b)
    roots = np.array([find(i) for i in range(n + 1)])
    comp = roots[labels]
    sizes = np.bincount(comp.ravel(), minlength=n + 1)
    big = (sizes > S_MAX)
    big[0] = False
    active &= ~big[comp]
    return active


def components(active):
    """Edge lists of the active components on the periodic lattice."""
    idx = {int(i): k for k, i in enumerate(np.flatnonzero(active.ravel()))}
    seen, comps = set(), []
    for start in idx:
        if start in seen:
            continue
        stack, members = [start], []
        seen.add(start)
        while stack:
            v = stack.pop()
            members.append(v)
            r, c = divmod(v, SIDE)
            for nb in (r * SIDE + (c + 1) % SIDE, r * SIDE + (c - 1) % SIDE,
                       ((r + 1) % SIDE) * SIDE + c, ((r - 1) % SIDE) * SIDE + c):
                if nb in idx and nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        local = {v: k for k, v in enumerate(sorted(members))}
        edges = set()
        for v in members:
            r, c = divmod(v, SIDE)
            for nb in (r * SIDE + (c + 1) % SIDE, ((r + 1) % SIDE) * SIDE + c):
                if nb in local:
                    edges.add(tuple(sorted((local[v], local[nb]))))
        comps.append((len(members), sorted(edges)))
    return comps


def shape_table():
    """The exact rate of every component shape, from the cluster expansion run."""
    data = json.loads(EXPANSION.read_text())["synthetic"]
    if any(abs(data["rates"][k] - v) > 0 for k, v in EXACT_RATES.items()):
        raise ValueError("expansion.json was solved with different rates")
    return [(row["size"], nx.Graph([tuple(e) for e in row["edges"]]), row["R"])
            for row in data["shapes"]]


def exact_rate(active, table):
    """Exact stationary product rate per lattice site, summed over components."""
    total = 0.0
    for size, edges in components(active):
        if size < 3:
            continue  # exactly zero by Result A
        g = nx.Graph(edges)
        for s, shape, rate in table:
            if s == size and nx.is_isomorphic(shape, g):
                total += rate
                break
        else:
            total += solve(size, edges, EXACT_RATES)["tof"] * size
    return total / active.size


def run(job):
    p, rep = job
    seed = int(1000 * p) * 100 + rep
    active = frozen_lattice(p, seed)
    expected = exact_rate(active, shape_table())

    sim = KMC(Lattice(SIDE), KMC_RATES, seed=seed + 7, init_phase=CUO)
    sim.phase = np.where(active.ravel(), CUO, CUSO4).astype(np.int8)
    sim.ads[:] = EMPTY
    for i in range(sim.lat.N):
        sim._refresh_site(i)
    for b in range(sim.lat.n_bonds):
        sim._refresh_bond(b)

    while sim.t < BURN:
        sim.step()
    edges = BURN + MEASURE * np.arange(1, BATCHES + 1) / BATCHES
    counts, start, k = [], sim.n_so3_net, 0
    t_start = sim.t
    while k < BATCHES:
        before_t, before_n = sim.t, sim.n_so3_net
        sim.step()
        # an event that crosses a batch edge is assigned to the batch it lands in
        while k < BATCHES and sim.t >= edges[k]:
            counts.append(before_n - start)
            start = before_n
            k += 1
    widths = np.diff(np.concatenate([[t_start], edges]))
    rates = np.array(counts) / (widths * active.size)
    mean = float(rates.mean())
    sem = float(rates.std(ddof=1) / np.sqrt(BATCHES))
    return dict(p=p, realisation=rep, seed=seed, active_sites=int(active.sum()),
                components=len(components(active)), exact=expected, kmc=mean, sem=sem,
                events=int(sim.n_events), z=(mean - expected) / sem,
                passed=bool(abs(mean - expected) <= 5 * sem + 0.01 * expected))


def main():
    jobs = [(p, rep) for p in COVERAGES for rep in range(REALISATIONS)]
    with ProcessPoolExecutor() as pool:
        rows = list(pool.map(run, jobs))
    print(f"Lattice KMC against the exact component sum, {SIDE}x{SIDE}, components <= {S_MAX}")
    print(f"  {'p':>5} {'rep':>4} {'active':>7} {'exact':>11} {'KMC':>11} {'sem':>9} {'z':>6}")
    for r in rows:
        print(f"  {r['p']:5.2f} {r['realisation']:4d} {r['active_sites']:7d} {r['exact']:11.5f} "
              f"{r['kmc']:11.5f} {r['sem']:9.5f} {r['z']:6.2f}")
    z = np.array([r["z"] for r in rows])
    chi2 = float(np.sum(z ** 2))
    summary = dict(runs=len(rows), all_passed=all(r["passed"] for r in rows),
                   chi2=chi2, dof=len(rows), chi2_p_value=float(stats.chi2.sf(chi2, len(rows))),
                   max_abs_z=float(np.max(np.abs(z))),
                   max_rel_dev=float(max(abs(r["kmc"] - r["exact"]) / r["exact"] for r in rows)),
                   total_events=int(sum(r["events"] for r in rows)))
    print(f"\n  runs {summary['runs']}, all passed {summary['all_passed']}, "
          f"chi2 {chi2:.1f} on {len(rows)} dof (p = {summary['chi2_p_value']:.3f}), "
          f"max |z| {summary['max_abs_z']:.2f}, max relative deviation {summary['max_rel_dev']:.4f}")
    OUT.write_text(json.dumps(dict(summary=summary, rows=rows, side=SIDE, s_max=S_MAX,
                                   burn=BURN, measure=MEASURE, batches=BATCHES,
                                   rates={str(k): v for k, v in KMC_RATES.items()}), indent=1))
    print(f"  wrote {OUT.relative_to(ROOT)}")
    if not summary["all_passed"]:
        raise AssertionError("lattice KMC disagrees with the exact component sum")


if __name__ == "__main__":
    main()
