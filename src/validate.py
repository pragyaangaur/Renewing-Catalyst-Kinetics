"""
Validate the mean field closure against exact lattice KMC.

The mean field model is only worth using for the dense parameter maps if it
agrees with KMC where KMC can be afforded. Points are chosen to span the live
catalytic regime, the poisoned regime, and the boundary between them.
"""

import time
import numpy as np

from kmc_core import Lattice, KMC, CUO
from params import build_rates
from meanfield import steady_state

CASES = [
    # label,               T,     Ea_sulfate, Ea_desulfate, k_spall_sulf
    ("live, low sulfation",  745.0, 1.45, 1.55, 1.0e-3),
    ("live, hotter",         800.0, 1.45, 1.55, 1.0e-3),
    ("boundary",             745.0, 1.30, 1.60, 1.0e-3),
    ("boundary, cooler",     700.0, 1.35, 1.55, 1.0e-3),
    ("poisoned",             745.0, 0.95, 1.80, 1.0e-3),
    ("renewal rescued",      745.0, 1.10, 1.70, 1.0e2),
]


def run_kmc(T, extra, L=32, t_eq=None, t_sample=None, seed=0, max_events=1_200_000):
    k = build_rates(T, 1000e-6, 0.05, 500e-6, params=extra)
    lat = Lattice(L)
    s = KMC(lat, k, seed=seed, init_phase=CUO)

    # equilibrate for a fixed fraction of the event budget, then sample
    n_eq = max_events // 3
    n = 0
    while n < n_eq and s.step():
        n += 1
    t0 = s.t
    d0 = s.n_so3_des
    a0 = s.n_so3_readsorbed
    g0 = s.n_so3_from_desulf
    cov_acc = None
    n_cov = 0
    while n < max_events and s.step():
        n += 1
        if n % 5000 == 0:
            c = s.coverages()
            if cov_acc is None:
                cov_acc = {kk: 0.0 for kk in c}
            for kk, vv in c.items():
                cov_acc[kk] += vv
            n_cov += 1
    dt = s.t - t0
    N = lat.N
    net = (s.n_so3_des - d0) - (s.n_so3_readsorbed - a0) + (s.n_so3_from_desulf - g0)
    gross = (s.n_so3_des - d0)
    cov = {kk: vv / n_cov for kk, vv in cov_acc.items()} if n_cov else s.coverages()
    return {
        "tof_gross": gross / (N * dt) if dt > 0 else 0.0,
        "tof_net": net / (N * dt) if dt > 0 else 0.0,
        "sim_time": dt,
        **cov,
    }


def main():
    print(f"{'case':22s} {'src':5s} {'f_CuSO4':>9s} {'f_CuO':>8s} {'f_Cu':>8s} "
          f"{'th_SO3':>9s} {'TOF net (1/s)':>14s}")
    print("-" * 82)
    rows = []
    for label, T, ea_s, ea_d, ksp in CASES:
        extra = dict(Ea_sulfate=ea_s, Ea_desulfate=ea_d, k_spall_sulf=ksp)
        k = build_rates(T, 1000e-6, 0.05, 500e-6, params=extra)

        mf = steady_state(k)
        t0 = time.time()
        km = run_kmc(T, extra, seed=11)
        el = time.time() - t0

        print(f"{label:22s} {'MF':5s} {mf['f_CuSO4']:9.4f} {mf['f_CuO']:8.4f} "
              f"{mf['f_Cu']:8.4f} {mf['th_SO3']:9.3e} {mf['tof_net']:14.4g}")
        print(f"{'':22s} {'KMC':5s} {km['f_CUSO4']:9.4f} {km['f_CUO']:8.4f} "
              f"{km['f_CU']:8.4f} {km['th_SO3']:9.3e} {km['tof_net']:14.4g}"
              f"   [{el:.0f}s, t={km['sim_time']:.2e}s]")
        rat = (km['tof_net'] / mf['tof_net']) if mf['tof_net'] > 0 else float('nan')
        print(f"{'':22s} {'ratio KMC/MF on TOF':>40s} = {rat:8.3f}")
        print()
        rows.append((label, mf, km))
    return rows


if __name__ == "__main__":
    main()
