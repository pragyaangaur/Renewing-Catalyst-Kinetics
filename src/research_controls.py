"""Mechanism, rate uncertainty, lattice-size and renewal validation controls."""

import csv
import itertools
import json
from pathlib import Path

import numpy as np

from ensemble_audit import MOTIFS, exact, simulate, transitions
from kmc_core import (
    KMC, CUO, N_SITE_PROC, N_BOND_PROC, P_ADS_SO2, P_DES_SO2,
    P_DES_SO3, B_ADS_O2, B_DES_O2, B_LH, P_SULFATE, P_OXIDISE_CU,
)
from params import build_rates

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "audit"


def redox(size, edges, consumption=1.0, refill=1.0):
    states = list(itertools.product((0, 1), repeat=size))
    indices = {state: index for index, state in enumerate(states)}
    generator = np.zeros((len(states), len(states)))
    reward = np.zeros(len(states))
    for index, state in enumerate(states):
        for site, oxygen in enumerate(state):
            if oxygen:
                target = list(state)
                target[site] = 0
                generator[index, indices[tuple(target)]] += consumption
                reward[index] += consumption
        for left, right in edges:
            if state[left] == state[right] == 0:
                target = list(state)
                target[left] = target[right] = 1
                generator[index, indices[tuple(target)]] += refill
        generator[index, index] = -generator[index].sum()
    balance = generator.T.copy()
    balance[-1] = 1
    rhs = np.zeros(len(states))
    rhs[-1] = 1
    probability = np.linalg.solve(balance, rhs)
    assert np.min(probability) >= -1e-10
    assert np.max(np.abs(generator.T @ probability)) < 1e-9
    return float(probability @ reward / size)


def component_counts(side, coordination, seed):
    rng = np.random.default_rng(seed)
    total = side * side
    active = np.zeros(total, dtype=bool)
    active[rng.choice(total, round(0.1 * total), replace=False)] = True
    offsets = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    if coordination == 6:
        offsets += [(1, 1), (-1, -1)]
    seen = set()
    counts = [0, 0, 0]
    for site in np.flatnonzero(active):
        if site in seen:
            continue
        pending = [int(site)]
        seen.add(int(site))
        size = 0
        while pending:
            current = pending.pop()
            size += 1
            row, column = divmod(current, side)
            for delta_row, delta_column in offsets:
                neighbor = ((row + delta_row) % side) * side + (column + delta_column) % side
                if active[neighbor] and neighbor not in seen:
                    seen.add(neighbor)
                    pending.append(neighbor)
        counts[min(size - 1, 2)] += size
    return [count / int(active.sum()) for count in counts]


class Graph:
    def __init__(self, size, edges):
        self.N = size
        self.n_bonds = len(edges)
        self.bond_ends = np.array(edges, dtype=int).reshape(-1, 2)
        self.site_bonds = [[index for index, edge in enumerate(edges) if site in edge]
                           for site in range(size)]


def core_check(rates):
    frozen = {(kind, process): 0.0 for kind, count in (("site", N_SITE_PROC), ("bond", N_BOND_PROC))
              for process in range(count)}
    for key in (("site", P_ADS_SO2), ("site", P_DES_SO2), ("site", P_DES_SO3),
                ("bond", B_ADS_O2), ("bond", B_DES_O2), ("bond", B_LH)):
        frozen[key] = rates[key]
    size, edges = MOTIFS["path3"]
    samples = []
    for seed in range(12):
        simulation = KMC(Graph(size, edges), frozen, seed=seed, init_phase=CUO)
        count = 0
        while simulation.t < 40:
            before = simulation.n_so3_net
            simulation.step()
            if 20 <= simulation.t <= 40:
                count += simulation.n_so3_net - before
        samples.append(count / (20 * size))
        for pool in simulation.site_pools + simulation.bond_pools:
            members = pool.items[:pool.n]
            assert len(set(members)) == pool.n
            assert all(pool.pos[member] == position for position, member in enumerate(members))
    expected = exact(size, edges, rates)["tof"]
    sem = float(np.std(samples, ddof=1) / np.sqrt(len(samples)))
    assert abs(np.mean(samples) - expected) < 5 * sem + 0.03 * expected
    return dict(exact=expected, mean=float(np.mean(samples)), sem=sem, samples=samples)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rates = build_rates(745, 0.001, 0.05, 0.0)
    results = {"original_engine_check": core_check(rates)}
    print("Original KMC engine agrees with exact path3 solution", flush=True)
    assert abs(redox(2, [(0, 1)]) - 0.4) < 1e-12
    results["mechanisms"] = [{"motif": name, "lh": exact(size, edges, rates)["tof"],
                               "redox": redox(size, edges)}
                              for name, (size, edges) in list(MOTIFS.items())[:7]]
    results["redox_shape_sweep"] = []
    for refill in np.logspace(-2, 2, 41):
        results["redox_shape_sweep"].append(dict(refill=float(refill),
            ratio=redox(*MOTIFS["star4"], refill=refill) / redox(*MOTIFS["path4"], refill=refill)))
    rng = np.random.default_rng(20260905)
    uncertain_keys = [("site", P_ADS_SO2), ("site", P_DES_SO2), ("site", P_DES_SO3),
                      ("bond", B_ADS_O2), ("bond", B_DES_O2), ("bond", B_LH)]
    results["joint_uncertainty"] = []
    for sample in range(100):
        factors = 10 ** rng.uniform(-2, 2, len(uncertain_keys))
        varied = dict(rates)
        for key, factor in zip(uncertain_keys, factors):
            varied[key] *= factor
        diffusion = float(10 ** rng.uniform(0, 6))
        path = exact(*MOTIFS["path4"], varied, diffusion)["tof"]
        star = exact(*MOTIFS["star4"], varied, diffusion)["tof"]
        results["joint_uncertainty"].append(dict(sample=sample, factors=factors.tolist(),
            diffusion=diffusion, path=path, star=star, ratio=star / path))
    print("Joint rate uncertainty and redox shape sweeps completed", flush=True)
    results["rate_sensitivity"] = []
    for process in (B_ADS_O2, B_DES_O2, B_LH):
        for factor in (0.01, 1.0, 100.0):
            varied = dict(rates)
            varied[("bond", process)] *= factor
            entry = {"process": process, "factor": factor}
            for name in ("pair", "path3", "star4", "path4"):
                entry[name] = exact(*MOTIFS[name], varied)["tof"]
            assert abs(entry["pair"]) < 1e-8
            assert entry["path3"] > 0
            results["rate_sensitivity"].append(entry)
    geometry = []
    for side in (32, 64, 128):
        for coordination in (4, 6):
            samples = np.array([component_counts(side, coordination, seed) for seed in range(40)])
            geometry.append(dict(side=side, coordination=coordination,
                                 mean=samples.mean(axis=0).tolist(),
                                 sem=(samples.std(axis=0, ddof=1) / np.sqrt(40)).tolist()))
    results["geometry"] = geometry
    results["square_bernoulli_limit"] = {"single": 0.9 ** 4, "pair": 4 * 0.1 * 0.9 ** 6,
                                         "three_plus": 1 - 0.9 ** 4 - 4 * 0.1 * 0.9 ** 6}
    results["renewal_structural_checks"] = []
    renewing = dict(rates)
    renewing[("site", P_OXIDISE_CU)] = 1
    renewing[("site", P_SULFATE)] = 1
    for renewal in (0.01, 1, 10000):
        result = exact(*MOTIFS["pair"], renewing, renewal=renewal)
        assert abs(result["tof"]) < 1e-8
        results["renewal_structural_checks"].append(result)
    synthetic = {key: 0.0 for key in rates}
    for key in (("site", P_ADS_SO2), ("site", P_DES_SO2), ("site", P_DES_SO3),
                ("site", P_OXIDISE_CU), ("bond", B_ADS_O2), ("bond", B_LH)):
        synthetic[key] = 1.0
    synthetic[("bond", B_DES_O2)] = 0.1
    synthetic[("site", P_SULFATE)] = 0.1
    expected = exact(*MOTIFS["square"], synthetic, renewal=0.1)["tof"]
    samples = [simulate(*MOTIFS["square"], synthetic, seed, 1000, 100,
                        renewal=0.1)["tof"] for seed in range(8)]
    sem = float(np.std(samples, ddof=1) / np.sqrt(len(samples)))
    assert abs(np.mean(samples) - expected) <= 5 * sem + 0.03 * expected
    results["renewal_ssa_check"] = dict(exact=expected, mean=float(np.mean(samples)),
        sem=sem, samples=samples, duration=1000, burn=100, renewal=0.1,
        rates={str(key): value for key, value in synthetic.items()})
    results["renewal_identity"] = []
    with (OUTPUT / "renewal.csv").open() as handle:
        renewal_rows = list(csv.DictReader(handle))
    for barrier in (0.95, 1.45):
        branch_rates = build_rates(745, 0.001, 0.05, 0,
                                   params={"Ea_sulfate": barrier})
        ceiling = branch_rates[("site", P_DES_SO3)] / branch_rates[("site", P_SULFATE)]
        selected = [row for row in renewal_rows if float(row["barrier"]) == barrier]
        errors = [abs(float(row["copper_yield"]) - ceiling * float(row["sulfate_coverage"]))
                  for row in selected]
        assert max(errors) < 1e-5 * max(ceiling, 1)
        assert all(float(row["copper_yield"]) <= ceiling * (1 + 1e-6) for row in selected)
        results["renewal_identity"].append(dict(barrier=barrier, ceiling=ceiling,
            points=len(selected), max_absolute_error=max(errors)))
    (OUTPUT / "controls.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
