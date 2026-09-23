"""The stochastic engines against exact answers, kept short enough to run often."""

import numpy as np

import ensemble_audit as ea
import research_controls as rc
from kmc_core import (KMC, CUO, N_SITE_PROC, N_BOND_PROC, P_ADS_SO2, P_DES_SO2,
                      P_DES_SO3, B_ADS_O2, B_DES_O2, B_LH)
from params import build_rates


def synthetic():
    rates = {("site", p): 0.0 for p in range(N_SITE_PROC)}
    rates.update({("bond", p): 0.0 for p in range(N_BOND_PROC)})
    rates.update({("site", P_ADS_SO2): 1.0, ("site", P_DES_SO2): 1.0, ("site", P_DES_SO3): 1.0,
                  ("bond", B_ADS_O2): 1.0, ("bond", B_DES_O2): 0.1, ("bond", B_LH): 1.0})
    return rates


def test_redox_pair_rate_matches_its_closed_form():
    for consume, refill in ((1.0, 1.0), (0.3, 7.0), (5.0, 0.2)):
        expected = 2 * consume * refill / (3 * refill + 2 * consume)
        assert abs(rc.redox(2, [(0, 1)], consume, refill) - expected) < 1e-12


def test_pool_based_kmc_matches_exact_three_site_path():
    rates = synthetic()
    size, edges = ea.MOTIFS["path3"]
    samples = []
    for seed in range(6):
        sim = KMC(rc.Graph(size, edges), rates, seed=seed, init_phase=CUO)
        count = 0
        while sim.t < 2200:
            before = sim.n_so3_net
            sim.step()
            if sim.t > 200:
                count += sim.n_so3_net - before
        samples.append(count / (2000 * size))
        for pool in sim.site_pools + sim.bond_pools:
            members = pool.items[:pool.n]
            assert len(set(members)) == pool.n
            assert all(pool.pos[m] == i for i, m in enumerate(members))
    expected = ea.exact(size, edges, rates)["tof"]
    sem = np.std(samples, ddof=1) / np.sqrt(len(samples))
    assert abs(np.mean(samples) - expected) < 5 * sem + 0.02 * expected


def test_gillespie_reference_matches_exact_square():
    rates = synthetic()
    size, edges = ea.MOTIFS["square"]
    expected = ea.exact(size, edges, rates)["tof"]
    samples = [ea.simulate(size, edges, rates, seed, 400, 50)["tof"] for seed in range(6)]
    sem = np.std(samples, ddof=1) / np.sqrt(len(samples))
    assert abs(np.mean(samples) - expected) < 5 * sem + 0.02 * expected


def test_ensemble_audit_exact_matches_structure_theory_solver():
    """Two independently written exact solvers agree on the illustrative motifs."""
    import structure_theory as st
    from kmc_core import P_ADS_SO3
    rates = build_rates(745, 0.001, 0.05, 0.0)
    r = dict(st.BASE_RATES, k_ads_so2=rates[("site", P_ADS_SO2)],
             k_des_so2=rates[("site", P_DES_SO2)], k_des_so3=rates[("site", P_DES_SO3)],
             k_ads_o2=rates[("bond", B_ADS_O2)], k_des_o2=rates[("bond", B_DES_O2)],
             k_lh=rates[("bond", B_LH)], k_ads_so3=rates[("site", P_ADS_SO3)])
    for name in ("path3", "triangle", "star4", "square", "path5"):
        size, edges = ea.MOTIFS[name]
        a = ea.exact(size, edges, rates)["tof"]
        b = st.solve(size, edges, r)["tof"]
        assert abs(a - b) <= 1e-9 * b
