"""Lattice animals, the dead fraction law and the cluster rate expansion."""

import math

import numpy as np
import pytest

import cluster_expansion as ce


@pytest.mark.parametrize("lattice", list(ce.LATTICES))
def test_animal_counts_match_oeis(lattice):
    animals = ce.enumerate_animals(lattice, 7)
    assert [len(animals[s]) for s in range(1, 8)] == ce.OEIS[lattice][:7]


@pytest.mark.parametrize("lattice", list(ce.LATTICES))
def test_sum_rule_is_exact_through_the_enumerated_order(lattice):
    _, poly = ce.perimeter_polynomials(lattice, 7)
    coeffs = ce.series_sum_rule(poly)
    assert np.max(np.abs(coeffs[:8])) < 1e-9
    assert abs(coeffs[8]) > 1e-9  # the first order that the truncation misses


@pytest.mark.parametrize("lattice,t2", [("honeycomb", 4), ("square", 6), ("triangular", 8)])
def test_pair_perimeter(lattice, t2):
    animals = ce.enumerate_animals(lattice, 2)
    assert {ce.perimeter(lattice, a) for a in animals[2]} == {t2}


@pytest.mark.parametrize("lattice", list(ce.LATTICES))
def test_closed_form_dead_fraction_equals_expansion(lattice):
    _, poly = ce.perimeter_polynomials(lattice, 2)
    p = np.linspace(0.01, 0.99, 50)
    n = ce.cluster_numbers(poly, p)
    assert np.allclose((n[1] + 2 * n[2]) / p, ce.dead_fraction(lattice, p), rtol=1e-12)


def test_square_lattice_at_ninety_percent_blocking():
    assert ce.dead_fraction("square", 0.1) == pytest.approx(0.9 ** 4 + 0.4 * 0.9 ** 6)
    assert ce.dead_fraction("square", 0.1) == pytest.approx(0.868676, abs=1e-6)


@pytest.mark.parametrize("lattice", list(ce.LATTICES))
def test_sampled_lattices_agree_with_the_closed_form(lattice):
    rows = ce.sampled_population(lattice, [0.1, 0.3], side=256, reps=6, seed=3)
    for row in rows:
        dead = row["mean"][0] + row["mean"][1]
        err = math.hypot(row["sem"][0], row["sem"][1])
        assert abs(dead - float(ce.dead_fraction(lattice, row["p"]))) < 5 * err + 1e-3


def test_saved_expansion_is_internally_consistent():
    import json
    data = json.loads((ce.OUT / "expansion.json").read_text())
    for name in ("illustrative_745K", "synthetic"):
        d = data[name]
        assert len(d["shapes"]) == 18
        assert d["max_residual"] < 1e-12
        assert all(row["R"] > 0 for row in d["shapes"])
        total = np.sum([v for v in d["contribution_by_size"].values()], axis=0)
        assert np.allclose(total, d["tof_truncated"], rtol=1e-12)
