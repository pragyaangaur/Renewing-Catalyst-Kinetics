"""Result A: one and two site patches are dead, three sites are enough."""

import numpy as np
import pytest

import structure_theory as st


def test_atlas_matches_bruteforce_enumeration_up_to_five_sites():
    import networkx as nx
    for n in range(2, 6):
        atlas = [nx.Graph(e) for e in st.connected_graphs(n)]
        brute = [nx.Graph(e) for e in st.connected_graphs_bruteforce(n)]
        assert len(atlas) == len(brute)
        for g in brute:
            assert sum(nx.is_isomorphic(g, h) for h in atlas) == 1


def test_connected_graph_counts_are_oeis_a001349():
    assert [len(st.connected_graphs(n)) for n in range(1, 8)] == [1, 1, 2, 6, 21, 112, 853]


@pytest.mark.parametrize("diffusion", [False, True])
def test_pair_is_irreducible_and_never_productive(diffusion):
    res = st.reachability(2, [(0, 1)], diffusion)
    assert res["irreducible"]
    assert not res["productive"]
    assert res["states"] == 5  # empty, SO2 on either site, SO2 on both, O on both


@pytest.mark.parametrize("n", [3, 4, 5])
def test_every_graph_of_three_or_more_sites_is_productive_and_irreducible(n):
    for edges in st.connected_graphs(n):
        for diffusion in (False, True):
            res = st.reachability(n, edges, diffusion)
            assert res["productive"] and res["irreducible"]


def test_every_connected_graph_on_three_vertices_or_more_has_a_two_edge_path():
    for n in range(3, 7):
        for edges in st.connected_graphs(n):
            degree = np.bincount(np.array(edges).ravel(), minlength=n)
            assert degree.max() >= 2


def test_pair_rate_is_exactly_zero_for_random_rates_and_diffusion():
    rng = np.random.default_rng(11)
    for _ in range(20):
        r = dict(st.BASE_RATES)
        for key in ("k_ads_so2", "k_des_so2", "k_des_so3", "k_ads_o2", "k_des_o2", "k_lh",
                    "k_diff"):
            r[key] = float(10.0 ** rng.uniform(-3, 7))
        assert st.solve(2, [(0, 1)], r)["tof"] == 0.0


def test_three_site_path_rate_is_positive_for_random_rates():
    rng = np.random.default_rng(12)
    for _ in range(20):
        r = dict(st.BASE_RATES)
        for key in ("k_ads_so2", "k_des_so2", "k_des_so3", "k_ads_o2", "k_des_o2", "k_lh"):
            r[key] = float(10.0 ** rng.uniform(-2, 7))
        assert st.solve(3, [(0, 1), (1, 2)], r)["tof"] > 0.0


def test_saved_result_a_is_consistent():
    import json
    rows = json.loads((st.ROOT / "results" / "structure" / "result_a.json").read_text())
    reach = json.loads((st.ROOT / "results" / "structure" /
                        "result_a_reachability.json").read_text())
    assert len(reach) == 996
    assert all(r["irreducible"] for r in reach)
    assert all(r["productive"] == (r["n"] >= 3) for r in reach)
    assert all(r["tof"] == 0.0 for r in rows if r["n"] < 3)
    assert all(r["tof"] > 0.0 for r in rows if r["n"] >= 3)
