"""Result B: product per copper atom, its ceiling, and the site-loss correspondence."""

import numpy as np
import pytest

import structure_theory as st

MOTIFS = [(3, [(0, 1), (1, 2)]), (4, [(0, 1), (0, 2), (0, 3)]),
          (4, [(0, 1), (1, 2), (2, 3), (0, 3)])]


def draw(rng, regen):
    r = dict(st.BASE_RATES)
    r["k_sulf"] = float(10.0 ** rng.uniform(1, 6))
    r["k_des_so3"] = float(10.0 ** rng.uniform(1, 6))
    r["k_desulf"] = float(10.0 ** rng.uniform(-3, 3)) if regen else 0.0
    r["gamma"] = float(10.0 ** rng.uniform(-2, 3))
    r["k_ox"] = float(10.0 ** rng.uniform(2, 5))
    r["k_diff"] = float(10.0 ** rng.uniform(-1, 5))
    return r


@pytest.mark.parametrize("regen", [False, True])
def test_identity_holds_to_machine_precision(regen):
    rng = np.random.default_rng(5 + regen)
    for size, edges in MOTIFS:
        for _ in range(4):
            r = draw(rng, regen)
            out = st.solve(size, edges, r)
            spend = r["gamma"] * size
            ratio = r["k_des_so3"] / r["k_sulf"]
            rho = out["R_regen"] / spend
            y = out["R_product"] / spend
            assert y == pytest.approx((1 + ratio) * rho + ratio * out["f_sulfate"], rel=1e-12)
            if not regen:
                assert y <= ratio * (1 + 1e-12)


def test_sulphate_balance_closes():
    rng = np.random.default_rng(9)
    r = draw(rng, True)
    out = st.solve(4, MOTIFS[2][1], r)
    removed = out["R_regen"] + r["gamma"] * 4 * out["f_sulfate"]
    assert out["R_sulfation"] == pytest.approx(removed, rel=1e-12)


def test_pair_makes_nothing_even_with_renewal_and_regeneration():
    rng = np.random.default_rng(10)
    r = draw(rng, True)
    out = st.solve(2, [(0, 1)], r)
    assert out["R_product"] == 0.0 and out["R_sulfation"] == 0.0


def test_saved_site_loss_relations():
    import json
    rows = json.loads((st.ROOT / "results" / "structure" / "site_loss_check.json").read_text())
    live = [r for r in rows if not r["dead"]]
    assert len(live) == 150
    assert max(r["err_inverse_S_cu"] for r in live) < 1e-14
    assert max(r["err_inverse_S_sulf"] for r in live) < 1e-12
    assert max(r["err_consumed"] for r in live) < 1e-14
    assert max(r["err_ceiling"] for r in live if r["err_ceiling"] is not None) < 1e-14
