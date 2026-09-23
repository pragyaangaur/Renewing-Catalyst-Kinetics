"""Why the first pass saw a dead surface, computed exactly.

The first pass geometry study (kmc_geometry.py) ran lattice KMC on a 64 by 64
square lattice with 90 percent of the sites sulphated at random and reported a
turnover of exactly zero. docs/CORRECTIONS.md withdrew that conclusion on general
grounds, because a finite run with no events does not show a zero rate. This script
replaces the general argument with numbers.

The frozen lattice of that run is rebuilt from its seed. Every active component in
it has five sites or fewer, and one at 85 percent has six or fewer, so the lattice
can be solved exactly component by component. Three things are computed.

1.  The exact stationary turnover of that lattice, from structure_theory.solve.
2.  The exact expected number of product events inside the window the first pass
    measured, starting from the empty surface it started from. This uses the
    transient solution of the master equation, p(t) = p(0) exp(Qt), with the
    cumulative product obtained from an augmented generator. The exponential is
    taken densely by scaling and squaring, since the rates span seven decades and
    Krylov or uniformisation methods would need millions of steps. It shows whether the
    first pass was unlucky or simply not yet in steady state.
3.  A much longer run of the same KMC engine on the same lattice, compared with
    the exact stationary turnover.

The rates are the first pass rates at 745 K with the phase frozen, which are
illustrative and not fitted to any measurement.
"""

import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

import kmc_geometry as kg
from cluster_expansion import illustrative_rates
from kmc_core import (KMC, Lattice, CUO, EMPTY, P_SULFATE, P_DESULFATE, P_SPALL_SULF,
                      P_SPALL_OX, P_OXIDISE_CU)
from lattice_validation import components
from params import build_rates
from structure_theory import build_chain, solve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "firstpass_postmortem.json"
SIDE = 64
COVERAGES = [0.0, 0.30, 0.50, 0.65, 0.75, 0.85, 0.90, 0.94, 0.97]
EVENTS = 700_000
LONG_BURN = 1.0      # seconds of simulated time
LONG_MEASURE = 4.0
BATCHES = 20


def first_pass_lattice(i):
    """The random geometry that kmc_geometry.main built for coverage index i."""
    rng = np.random.default_rng(100 + i)
    return kg.make_random(SIDE, COVERAGES[i], rng)


def first_pass_rates():
    k = dict(build_rates(kg.T, kg.DUTY["y_SO2"], kg.DUTY["y_O2"], kg.DUTY["y_SO3"],
                         params=kg.PARAMS))
    for pid in (P_SULFATE, P_DESULFATE, P_SPALL_SULF, P_SPALL_OX, P_OXIDISE_CU):
        k[("site", pid)] = 0.0
    return k


def engine(phase, seed):
    sim = KMC(Lattice(SIDE), first_pass_rates(), seed=seed, init_phase=CUO)
    sim.phase = phase.copy()
    sim.ads[:] = EMPTY
    for j in range(sim.lat.N):
        sim._refresh_site(j)
    for b in range(sim.lat.n_bonds):
        sim._refresh_bond(b)
    return sim


def replay_first_pass(phase, seed):
    """The first pass measurement, repeated event for event from its seed."""
    sim = engine(phase, seed)
    n = 0
    while n < EVENTS // 3 and sim.step():
        n += 1
    t0, d0, a0 = sim.t, sim.n_so3_des, sim.n_so3_readsorbed
    while n < EVENTS and sim.step():
        n += 1
    net = (sim.n_so3_des - d0) - (sim.n_so3_readsorbed - a0)
    return t0, sim.t - t0, int(net)


def exact_component_terms(phase, rates, t0, dt):
    """Stationary rate and expected product in [t0, t0 + dt] for each component shape."""
    active = (phase == CUO).reshape(SIDE, SIDE)
    shapes = []
    for size, edges in components(active):
        g = nx.Graph(edges)
        g.add_nodes_from(range(size))
        for entry in shapes:
            if entry["size"] == size and nx.is_isomorphic(entry["graph"], g):
                entry["count"] += 1
                break
        else:
            shapes.append(dict(size=size, graph=g, edges=edges, count=1))
    rows = []
    for entry in shapes:
        size, edges = entry["size"], entry["edges"]
        stationary = solve(size, edges, rates)["tof"] * size
        chain = build_chain(size, edges, rates)
        n = len(chain["states"])
        # augmented generator: the last coordinate integrates the product rate
        aug = np.zeros((n + 1, n + 1))
        aug[:n, :n] = chain["gen"].toarray()
        aug[:n, n] = chain["prod"]
        x0 = np.zeros(n + 1)
        x0[0] = 1.0  # the empty patch, which is where the first pass started
        x_t0 = x0 @ expm(aug * t0)
        x_t1 = x0 @ expm(aug * (t0 + dt))
        expected = float(x_t1[-1] - x_t0[-1])
        rows.append(dict(size=size, edges=edges, count=entry["count"],
                         stationary_rate=stationary, expected_in_window=expected,
                         states=n))
    return rows


def long_run(phase, seed):
    sim = engine(phase, seed)
    while sim.t < LONG_BURN:
        sim.step()
    edges = LONG_BURN + LONG_MEASURE * np.arange(1, BATCHES + 1) / BATCHES
    counts, start, k, t_start = [], sim.n_so3_net, 0, sim.t
    while k < BATCHES:
        before = sim.n_so3_net
        sim.step()
        while k < BATCHES and sim.t >= edges[k]:
            counts.append(before - start)
            start = before
            k += 1
    widths = np.diff(np.concatenate([[t_start], edges]))
    rates = np.array(counts) / (widths * SIDE * SIDE)
    return float(rates.mean()), float(rates.std(ddof=1) / np.sqrt(BATCHES)), int(sim.n_events)


def main():
    rates = illustrative_rates()
    out = []
    for i in (5, 6):
        phase = first_pass_lattice(i)
        seed = 200 + i * 10
        t0, dt, observed = replay_first_pass(phase, seed)
        terms = exact_component_terms(phase, rates, t0, dt)
        stationary = sum(r["stationary_rate"] * r["count"] for r in terms) / (SIDE * SIDE)
        expected = sum(r["expected_in_window"] * r["count"] for r in terms)
        mean, sem, events = long_run(phase, seed + 1)
        row = dict(coverage=COVERAGES[i], active_sites=int((phase == CUO).sum()),
                   largest_component=max(r["size"] for r in terms),
                   burn_in_end=t0, window=dt, observed_events=observed,
                   first_pass_tof=observed / (SIDE * SIDE * dt),
                   exact_stationary_tof=stationary,
                   expected_events_stationary=stationary * SIDE * SIDE * dt,
                   expected_events_transient=expected,
                   poisson_p_zero_transient=float(np.exp(-expected)),
                   long_run=dict(burn=LONG_BURN, measure=LONG_MEASURE, tof=mean, sem=sem,
                                 events=events, z=(mean - stationary) / sem),
                   shapes=[{k: v for k, v in r.items() if k != "graph"} for r in terms])
        out.append(row)
        print(f"sulphated {row['coverage']:.2f}: first pass window {dt:.4g} s after "
              f"{t0:.4g} s, observed {observed} events, TOF {row['first_pass_tof']:.4g}")
        print(f"  exact stationary TOF {stationary:.5g}; expected events in that window "
              f"{row['expected_events_stationary']:.3g} at steady state, "
              f"{expected:.3g} from the empty start (P(0) = {row['poisson_p_zero_transient']:.3g})")
        print(f"  long run of the same engine: TOF {mean:.5g} +/- {sem:.2g} "
              f"(z = {row['long_run']['z']:.2f}, {events} events)")
    OUT.write_text(json.dumps(out, indent=1))
    print(f"  wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
