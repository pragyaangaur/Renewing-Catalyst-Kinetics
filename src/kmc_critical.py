"""
Exact KMC check of the central result, at the point where mean field is most
likely to be wrong.

The headroom result says the surface may be almost entirely sulphated and still
carry the required duty, because the few surviving oxide sites turn over fast.
Mean field reaches that conclusion by treating the surviving sites as randomly
scattered and independent. That assumption is exactly what fails if sulphate
grows in patches, because the productive step is a Langmuir-Hinshelwood reaction
between adsorbates on neighbouring sites. An isolated oxide site with no oxide
neighbour cannot run that step at all.

So the question this script answers is: at high sulphate coverage, does the
lattice keep enough connected oxide to sustain turnover, or does mean field
overstate the surviving activity?

The comparison is run along a sweep of sulphate stability, and additionally the
spatial structure of the surviving oxide is measured by cluster analysis.
"""

import sys
import numpy as np

from kmc_core import Lattice, KMC, CUO, CU, CUSO4
from params import build_rates
from meanfield import steady_state
from reactor import required_tof, ILLUSTRATIVE_CONVERSION

DUTY = dict(y_SO2=1000e-6, y_O2=0.05, y_SO3=1e-9, tau=1.0)
T = 745.0
LIVE = dict(Ea_sulfate=1.45)


def oxide_clusters(phase, L):
    """
    Label connected clusters of oxide sites on the periodic lattice and return
    the size distribution and the fraction of oxide sites that have at least one
    oxide neighbour. Isolated oxide sites cannot run the bimolecular step.
    """
    grid = (phase.reshape(L, L) == CUO)
    seen = np.zeros_like(grid, dtype=bool)
    sizes = []
    n_isolated = 0
    for r in range(L):
        for c in range(L):
            if not grid[r, c]:
                continue
            nb = (grid[(r + 1) % L, c] or grid[(r - 1) % L, c]
                  or grid[r, (c + 1) % L] or grid[r, (c - 1) % L])
            if not nb:
                n_isolated += 1
            if seen[r, c]:
                continue
            stack = [(r, c)]
            seen[r, c] = True
            n = 0
            while stack:
                a, b = stack.pop()
                n += 1
                for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    x, y = (a + da) % L, (b + db) % L
                    if grid[x, y] and not seen[x, y]:
                        seen[x, y] = True
                        stack.append((x, y))
            sizes.append(n)
    n_ox = int(grid.sum())
    return {
        "n_oxide": n_ox,
        "n_clusters": len(sizes),
        "max_cluster": max(sizes) if sizes else 0,
        "mean_cluster": float(np.mean(sizes)) if sizes else 0.0,
        "frac_isolated": (n_isolated / n_ox) if n_ox else 0.0,
    }


def seed_from_meanfield(sim, mf, rng):
    """
    Initialise the lattice phase field by drawing each site independently from
    the mean field steady state composition.

    This is what makes the comparison meaningful. Sulphate decomposition can be
    as slow as 1e-2 per second while the surface chemistry runs at 1e6 per
    second, so a KMC started from a clean oxide spends its whole event budget on
    fast processes and never relaxes the slow phase variables. Starting at the
    mean field composition removes that problem and turns the run into a direct
    test of two things: whether the composition is stationary there, and whether
    the spatial arrangement of the surviving oxide changes the turnover that mean
    field predicts from the same composition.
    """
    n = sim.lat.N
    u = rng.random(n)
    f_so4 = mf["f_CuSO4"]
    f_cu = mf["f_Cu"]
    for i in range(n):
        if u[i] < f_so4:
            sim.phase[i] = CUSO4
        elif u[i] < f_so4 + f_cu:
            sim.phase[i] = CU
        else:
            sim.phase[i] = CUO
        sim.ads[i] = 0
    for i in range(n):
        sim._refresh_site(i)
    for b in range(sim.lat.n_bonds):
        sim._refresh_bond(b)


def run_point(Ea_desulf, L=48, max_events=3_000_000, seed=3, seeded=True):
    p = dict(LIVE, Ea_desulfate=float(Ea_desulf))
    k = build_rates(T, DUTY["y_SO2"], DUTY["y_O2"], DUTY["y_SO3"], params=p)

    mf0 = steady_state(k)
    lat = Lattice(L)
    s = KMC(lat, k, seed=seed, init_phase=CUO)
    if seeded:
        seed_from_meanfield(s, mf0, np.random.default_rng(seed + 991))
    f_init = float(np.count_nonzero(s.phase == CUSO4)) / lat.N

    n_eq = max_events // 2
    n = 0
    while n < n_eq and s.step():
        n += 1

    t0, d0, a0, g0 = s.t, s.n_so3_des, s.n_so3_readsorbed, s.n_so3_from_desulf
    acc = None
    n_acc = 0
    clus = []
    while n < max_events and s.step():
        n += 1
        if n % 10000 == 0:
            c = s.coverages()
            if acc is None:
                acc = {kk: 0.0 for kk in c}
            for kk, vv in c.items():
                acc[kk] += vv
            n_acc += 1
        if n % 250000 == 0:
            clus.append(oxide_clusters(s.phase, L))

    dt = s.t - t0
    net = (s.n_so3_des - d0) - (s.n_so3_readsorbed - a0) + (s.n_so3_from_desulf - g0)
    cov = {kk: vv / n_acc for kk, vv in acc.items()} if n_acc else s.coverages()
    tof = net / (lat.N * dt) if dt > 0 else 0.0

    cl = {}
    if clus:
        for key in clus[0]:
            cl[key] = float(np.mean([c[key] for c in clus]))

    mf = steady_state(k)
    return {
        "Ea_desulf": Ea_desulf,
        "kmc_f_CuSO4": cov["f_CUSO4"], "kmc_f_CuO": cov["f_CUO"], "kmc_tof": tof,
        "mf_f_CuSO4": mf["f_CuSO4"], "mf_f_CuO": mf["f_CuO"], "mf_tof": mf["tof_net"],
        "sim_time": dt, "clusters": cl, "f_init": f_init,
        "drift": cov["f_CUSO4"] - f_init,
    }


def main():
    req = required_tof(ILLUSTRATIVE_CONVERSION, DUTY["y_SO2"], DUTY["tau"], T)
    print(f"Reactor requirement at {T:.0f} K: TOF >= {req:.4f} per site per second\n")
    print("Lattice seeded at the mean field composition. 'drift' is how far the KMC")
    print("composition moved from that seed, so a small drift means the mean field")
    print("steady state is genuinely stationary under the exact dynamics.\n")
    print(f"{'Ea_desulf':>10} {'f_SO4 MF':>9} {'f_SO4 KMC':>10} {'drift':>8} "
          f"{'TOF MF':>11} {'TOF KMC':>11} {'KMC/MF':>8} "
          f"{'isolated':>9} {'maxclus':>8} {'t sim':>9} {'meets?':>7}")
    print("-" * 118)

    rows = []
    for ed in [1.30, 1.55, 1.75, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40]:
        r = run_point(ed)
        ratio = r["kmc_tof"] / r["mf_tof"] if r["mf_tof"] > 0 else float("nan")
        cl = r["clusters"]
        meets = "yes" if r["kmc_tof"] >= req else "NO"
        print(f"{ed:10.2f} {r['mf_f_CuSO4']:9.5f} {r['kmc_f_CuSO4']:10.5f} "
              f"{r['drift']:+8.4f} "
              f"{r['mf_tof']:11.4g} {r['kmc_tof']:11.4g} {ratio:8.3f} "
              f"{cl.get('frac_isolated', 0):9.3f} {cl.get('max_cluster', 0):8.0f} "
              f"{r['sim_time']:9.2e} {meets:>7}")
        rows.append(r)

    np.savez_compressed("../results/kmc_critical.npz",
                        Ea_desulf=np.array([r["Ea_desulf"] for r in rows]),
                        kmc_f=np.array([r["kmc_f_CuSO4"] for r in rows]),
                        mf_f=np.array([r["mf_f_CuSO4"] for r in rows]),
                        kmc_tof=np.array([r["kmc_tof"] for r in rows]),
                        mf_tof=np.array([r["mf_tof"] for r in rows]),
                        frac_isolated=np.array([r["clusters"].get("frac_isolated", 0)
                                                for r in rows]),
                        max_cluster=np.array([r["clusters"].get("max_cluster", 0)
                                              for r in rows]),
                        drift=np.array([r["drift"] for r in rows]),
                        sim_time=np.array([r["sim_time"] for r in rows]),
                        req=np.array([req]))
    print("\n  wrote ../results/kmc_critical.npz")


if __name__ == "__main__":
    main()
