"""
Does the spatial pattern of sulphation matter, independently of how much of it
there is?

The previous result showed that mean field badly overestimates surviving activity
at high sulphate coverage, because the productive chemistry needs two adjacent
oxide sites and dilute oxide sites stop touching one another. If that reading is
right then coverage alone does not determine activity. The same sulphate coverage
arranged as a few compact islands leaves the remaining oxide connected, while the
same coverage scattered at random strands it.

This script tests that directly. Lattices are built at matched sulphate coverage
with controlled geometry, the phase field is then frozen, and only the adsorbate
chemistry is allowed to run. Freezing the phase removes the confound that
different geometries would otherwise relax to different coverages, so what is
measured is the geometric effect and nothing else.

Geometries
    random   each site sulphated independently, the implicit assumption of mean field
    islands  sulphate grown by Eden cluster growth from a controlled number of
             nucleation seeds, so fewer seeds gives larger and more compact islands
"""

import numpy as np

from kmc_core import Lattice, KMC, CU, CUO, CUSO4, EMPTY
from params import build_rates
from reactor import required_tof, ILLUSTRATIVE_CONVERSION

T = 745.0
DUTY = dict(y_SO2=1000e-6, y_O2=0.05, y_SO3=1e-9, tau=1.0)
PARAMS = dict(Ea_sulfate=1.45, Ea_desulfate=1.55)


def make_random(L, cov, rng):
    n = L * L
    phase = np.full(n, CUO, dtype=np.int8)
    n_s = int(round(cov * n))
    idx = rng.choice(n, size=n_s, replace=False)
    phase[idx] = CUSO4
    return phase


def make_islands(L, cov, n_seeds, rng):
    """Eden growth of n_seeds sulphate islands up to the target coverage."""
    n = L * L
    phase = np.full(n, CUO, dtype=np.int8)
    n_s = int(round(cov * n))
    if n_s == 0:
        return phase
    n_seeds = max(1, min(n_seeds, n_s))

    seeds = rng.choice(n, size=n_seeds, replace=False)
    for s in seeds:
        phase[s] = CUSO4
    count = n_seeds

    def nbrs(i):
        r, c = divmod(i, L)
        return (r * L + (c + 1) % L, r * L + (c - 1) % L,
                ((r + 1) % L) * L + c, ((r - 1) % L) * L + c)

    # perimeter set of oxide sites adjacent to sulphate
    perim = set()
    for s in seeds:
        for j in nbrs(int(s)):
            if phase[j] == CUO:
                perim.add(j)

    while count < n_s and perim:
        j = int(rng.choice(list(perim)))
        perim.discard(j)
        if phase[j] != CUO:
            continue
        phase[j] = CUSO4
        count += 1
        for m in nbrs(j):
            if phase[m] == CUO:
                perim.add(m)

    # if growth stalled, top up at random
    while count < n_s:
        j = int(rng.integers(n))
        if phase[j] == CUO:
            phase[j] = CUSO4
            count += 1
    return phase


def pair_stats(phase, L):
    grid = (phase.reshape(L, L) == CUO)
    q = grid.mean()
    if q == 0:
        return 0.0, 0.0
    nb = (np.roll(grid, 1, 0) | np.roll(grid, -1, 0)
          | np.roll(grid, 1, 1) | np.roll(grid, -1, 1))
    connected = float((grid & nb).sum()) / float(grid.sum())
    # fraction of lattice bonds that join two oxide sites
    bonds = (grid & np.roll(grid, -1, 0)).sum() + (grid & np.roll(grid, -1, 1)).sum()
    return connected, float(bonds) / (2 * L * L)


def measure_tof(phase, L, n_events=700_000, seed=0):
    """
    Freeze the phase field and run only the adsorbate chemistry, then report the
    steady state turnover frequency per total site.
    """
    k = build_rates(T, DUTY["y_SO2"], DUTY["y_O2"], DUTY["y_SO3"], params=PARAMS)
    k = dict(k)
    # freeze every process that would change the phase field
    from kmc_core import (P_SULFATE, P_DESULFATE, P_SPALL_SULF,
                          P_SPALL_OX, P_OXIDISE_CU)
    for pid in (P_SULFATE, P_DESULFATE, P_SPALL_SULF, P_SPALL_OX, P_OXIDISE_CU):
        k[("site", pid)] = 0.0

    lat = Lattice(L)
    s = KMC(lat, k, seed=seed, init_phase=CUO)
    s.phase = phase.copy()
    s.ads[:] = EMPTY
    for i in range(lat.N):
        s._refresh_site(i)
    for b in range(lat.n_bonds):
        s._refresh_bond(b)

    n_eq = n_events // 3
    n = 0
    while n < n_eq and s.step():
        n += 1
    t0, d0, a0 = s.t, s.n_so3_des, s.n_so3_readsorbed
    while n < n_events and s.step():
        n += 1
    dt = s.t - t0
    if dt <= 0:
        return 0.0
    net = (s.n_so3_des - d0) - (s.n_so3_readsorbed - a0)
    return net / (lat.N * dt)


def main():
    L = 64
    req = required_tof(ILLUSTRATIVE_CONVERSION, DUTY["y_SO2"], DUTY["tau"], T)
    covs = [0.0, 0.30, 0.50, 0.65, 0.75, 0.85, 0.90, 0.94, 0.97]
    geoms = [("random", None), ("islands, 64 seeds", 64),
             ("islands, 16 seeds", 16), ("islands, 4 seeds", 4)]

    print(f"Turnover at matched sulphate coverage, phase field frozen, L = {L}")
    print(f"Reactor requirement: TOF >= {req:.3f} per site per second\n")
    header = f"{'coverage':>9}" + "".join(f"{g[0]:>20}" for g in geoms)
    print(header)
    print("-" * len(header))

    table = np.zeros((len(covs), len(geoms)))
    conn = np.zeros_like(table)
    for i, cov in enumerate(covs):
        rng = np.random.default_rng(100 + i)
        cells = []
        for j, (name, seeds) in enumerate(geoms):
            if seeds is None:
                ph = make_random(L, cov, rng)
            else:
                ph = make_islands(L, cov, seeds, rng)
            c, _ = pair_stats(ph, L)
            tof = measure_tof(ph, L, seed=200 + i * 10 + j)
            table[i, j] = tof
            conn[i, j] = c
            flag = "" if tof >= req else " *"
            cells.append(f"{tof:15.4g}{flag:>5}")
        print(f"{cov:9.2f}" + "".join(cells))

    print("\n  * marks a point that fails the reactor requirement\n")
    print("Fraction of surviving oxide sites having at least one oxide neighbour:")
    print(header)
    print("-" * len(header))
    for i, cov in enumerate(covs):
        print(f"{cov:9.2f}" + "".join(f"{conn[i, j]:20.3f}" for j in range(len(geoms))))

    np.savez_compressed("../results/geometry.npz",
                        cov=np.array(covs), tof=table, connected=conn,
                        labels=np.array([g[0] for g in geoms]),
                        req=np.array([req]))
    print("\n  wrote ../results/geometry.npz")


if __name__ == "__main__":
    main()
