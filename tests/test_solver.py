"""The stationary solver is the foundation of every exact number in the repository."""

import mpmath
import numpy as np
import pytest
from scipy.sparse import coo_matrix

import structure_theory as st


def random_generator(n, spread, rng):
    """A dense irreducible generator with rates spread over the given decades."""
    q = np.zeros((n, n))
    for i in range(n):
        q[i, (i + 1) % n] = 10.0 ** rng.uniform(-spread, spread)
        for j in rng.choice(n, 2, replace=False):
            if j != i:
                q[i, j] += 10.0 ** rng.uniform(-spread, spread)
    np.fill_diagonal(q, -q.sum(axis=1))
    return q


def mp_stationary(q):
    n = q.shape[0]
    a = mpmath.matrix(q.T.tolist())
    for j in range(n):
        a[n - 1, j] = 1
    rhs = mpmath.zeros(n, 1)
    rhs[n - 1] = 1
    return np.array([float(x) for x in mpmath.lu_solve(a, rhs)])


@pytest.mark.parametrize("seed", range(6))
def test_gth_is_entrywise_accurate_against_60_digit_arithmetic(seed):
    mpmath.mp.dps = 60
    rng = np.random.default_rng(seed)
    q = random_generator(12, 6, rng)
    exact = mp_stationary(q)
    gth = st.gth_solve(q)
    assert np.all(gth > 0)
    assert np.max(np.abs(gth - exact) / exact) < 1e-11


def test_solve_rejects_nothing_on_a_live_patch_and_residual_is_tiny():
    out = st.solve(3, [(0, 1), (1, 2)], dict(st.BASE_RATES))
    assert out["tof"] > 0
    assert out["residual"] < 1e-14


def test_tiny_rates_match_high_precision():
    """The case that exposed the old least squares solve: a rate near 1e-21."""
    mpmath.mp.dps = 60
    r = dict(st.BASE_RATES)
    rng = np.random.default_rng(3)
    for key in ("k_ads_so2", "k_des_so2", "k_des_so3", "k_ads_o2", "k_des_o2", "k_lh"):
        r[key] = float(10.0 ** rng.uniform(-2, 7))
    edges = [(0, 1), (1, 2)]
    out = st.solve(3, edges, r)
    # rebuild the same chain in 60 digit arithmetic
    start = (st.EMPTY,) * 3
    states, index, rows, cols, vals, prod = [start], {start: 0}, [], [], [], []
    i = 0
    while i < len(states):
        p = 0.0
        for nxt, rate, pr, _, _ in st.transitions(states[i], edges, r):
            if nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
            rows.append(i); cols.append(index[nxt]); vals.append(rate)
            p += rate * pr
        prod.append(p)
        i += 1
    n = len(states)
    q = coo_matrix((vals, (rows, cols)), shape=(n, n)).toarray()
    np.fill_diagonal(q, -q.sum(axis=1))
    pi = mp_stationary(q)
    exact = float(pi @ np.array(prod)) / 3
    assert out["tof"] == pytest.approx(exact, rel=1e-10)
