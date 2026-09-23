"""Exact activity of a randomly deactivated surface, built from its small components.

When a fraction 1 - p of the sites on a lattice is deactivated at random and then
frozen, the surviving active sites fall into connected components. Deactivated
sites take part in no process, so each component evolves independently of every
other, and the stationary state of the whole surface is the product of the
stationary states of its components. The product rate per lattice site is then a
sum over component shapes

    TOF(p) = sum over lattice animals A of  w_A * p^|A| * (1 - p)^t(A) * R(A)

where the sum runs over fixed lattice animals counted once per translation class,
w_A is one over the number of sites in the unit cell, t(A) is the number of
perimeter sites of A, and R(A) is the total stationary product rate of the
isolated component, which structure_theory.solve computes exactly. The factor
p^|A| (1 - p)^t(A) is the probability that a given placement of A is occupied and
completely surrounded by deactivated sites. This is the classical cluster number
expansion of site percolation (Sykes and Glen 1976, Stauffer and Aharony 1994),
paired here with exact master equation rates.

Two consequences are used in the paper.

The dead fraction. By Result A, a component of one or two sites never makes
product under the dual-site mechanism. The fraction of surviving sites that sit in
such components is exact and closed form on any lattice

    D(p) = (1 - p)^z + z p (1 - p)^t2

where z is the coordination number and t2 = 2z - 2 - c is the perimeter of a
nearest neighbour pair, with c the number of common neighbours of the pair. That
gives t2 = 4 on the honeycomb lattice, 6 on the square lattice and 8 on the
triangular lattice. Under the lattice oxygen control, where a pair is productive,
only isolated sites are dead and D(p) = (1 - p)^z.

The low coverage rate. Truncating the sum at components of s_max sites gives the
exact contribution of every component up to that size. What is left out is the
contribution of larger components, and the mass fraction of active sites in them
is known exactly from the same expansion.

The enumeration is checked against the published animal counts (OEIS A001168 for
the square lattice, A001207 for the triangular lattice and A001420 for the
honeycomb lattice) and against the sum rule that the cluster numbers must
reproduce p term by term. The cluster size statistics are checked against direct
sampling of random lattices.
"""

import argparse
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from structure_theory import BASE_RATES, solve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "cluster"

# Coordination number, common neighbours of an adjacent pair, sites per unit cell
LATTICES = {
    "honeycomb": dict(z=3, common=0, cell=2),
    "square": dict(z=4, common=0, cell=1),
    "triangular": dict(z=6, common=2, cell=1),
}

# Published fixed animal counts for sizes 1, 2, 3 and so on. The honeycomb entries
# are fixed polyiamonds, which count animals per unit cell of two sites.
OEIS = {
    "square": [1, 2, 6, 19, 63, 216, 760, 2725, 9910],          # A001168
    "triangular": [1, 3, 11, 44, 186, 814, 3652, 16689],        # A001207
    "honeycomb": [2, 3, 6, 14, 36, 94, 250, 675, 1838],         # A001420
}


# ------------------------------------------------------------------- geometry

def neighbours(lattice, site):
    """Nearest neighbours of a site. The honeycomb lattice uses brick wall coordinates."""
    x, y = site
    if lattice == "square":
        return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]
    if lattice == "triangular":
        return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1), (x + 1, y + 1), (x - 1, y - 1)]
    if lattice == "honeycomb":
        return [(x + 1, y), (x - 1, y), (x, y + 1) if (x + y) % 2 == 0 else (x, y - 1)]
    raise ValueError(lattice)


def normalise(lattice, animal):
    """Translate an animal to a canonical position using lattice translations only."""
    mx, my = min(animal)
    if lattice == "honeycomb":
        # translations must preserve the sublattice, so the shift has an even sum
        dy = -my + ((mx + my) % 2)
        return frozenset((x - mx, y + dy) for x, y in animal)
    return frozenset((x - mx, y - my) for x, y in animal)


def enumerate_animals(lattice, s_max):
    """
    Every fixed lattice animal up to s_max sites, one per translation class.

    Grown one site at a time from the single site. Each class is stored in its
    normalised form, so a class reached along several growth paths is kept once.
    On the honeycomb lattice both sublattice origins are seeded, because the two
    sites of the unit cell are not related by a translation.
    """
    if lattice == "honeycomb":
        level = {normalise(lattice, {(0, 0)}), normalise(lattice, {(0, 1)})}
    else:
        level = {frozenset({(0, 0)})}
    out = {1: sorted(level, key=sorted)}
    for s in range(2, s_max + 1):
        nxt = set()
        for animal in level:
            for site in animal:
                for nb in neighbours(lattice, site):
                    if nb not in animal:
                        nxt.add(normalise(lattice, animal | {nb}))
        level = nxt
        out[s] = sorted(level, key=sorted)
    return out


def perimeter(lattice, animal):
    return len({nb for site in animal for nb in neighbours(lattice, site)} - set(animal))


def animal_graph(lattice, animal):
    """The adjacency graph of an animal, relabelled onto 0 .. s-1."""
    sites = sorted(animal)
    index = {site: i for i, site in enumerate(sites)}
    edges = set()
    for site in sites:
        for nb in neighbours(lattice, site):
            if nb in index:
                edges.add(tuple(sorted((index[site], index[nb]))))
    return len(sites), sorted(edges)


def perimeter_polynomials(lattice, s_max):
    """
    For each size s, the list of (count per site, perimeter) pairs, so that the
    number of components of size s per lattice site is sum count * p^s (1-p)^t.
    """
    cell = LATTICES[lattice]["cell"]
    animals = enumerate_animals(lattice, s_max)
    poly = {}
    for s, group in animals.items():
        tally = {}
        for a in group:
            t = perimeter(lattice, a)
            tally[t] = tally.get(t, 0) + 1
        poly[s] = sorted((c / cell, t) for t, c in tally.items())
    return animals, poly


def cluster_numbers(poly, p):
    """n_s(p), the number of components of s sites per lattice site."""
    p = np.asarray(p, dtype=float)
    return {s: sum(c * p ** s * (1 - p) ** t for c, t in terms) for s, terms in poly.items()}


def dead_fraction(lattice, p, mechanism="dual-site"):
    """Closed form share of surviving sites in components that can never make product."""
    z = LATTICES[lattice]["z"]
    t2 = 2 * z - 2 - LATTICES[lattice]["common"]
    p = np.asarray(p, dtype=float)
    single = (1 - p) ** z
    if mechanism == "lattice-oxygen":
        return single
    return single + z * p * (1 - p) ** t2


def series_sum_rule(poly):
    """
    Coefficients of p - sum_s s n_s(p) as a polynomial in p.

    Below the percolation threshold every active site belongs to a finite
    component, so the full sum reproduces p exactly and the truncated sum must
    agree with p through order p^s_max. Nonzero low order coefficients would mean
    the enumeration or the perimeters are wrong.
    """
    s_max = max(poly)
    t_max = max(t for terms in poly.values() for _, t in terms)
    coeffs = np.zeros(s_max + t_max + 1)
    coeffs[1] = 1.0
    for s, terms in poly.items():
        for c, t in terms:
            # s * c * p^s * (1 - p)^t expanded binomially
            for k in range(t + 1):
                coeffs[s + k] -= s * c * math.comb(t, k) * (-1) ** k
    return coeffs


# ------------------------------------------------------------------ sampling

def sample_lattice(lattice, side, p, rng):
    """Component sizes of the active sites on a periodic random lattice."""
    n = side * side
    active = rng.random(n) < p
    x, y = np.divmod(np.arange(n), side)
    rows, cols = [], []
    if lattice == "square":
        shifts = [(1, 0), (0, 1)]
    elif lattice == "triangular":
        shifts = [(1, 0), (0, 1), (1, 1)]
    else:
        shifts = [(1, 0)]
    for dx, dy in shifts:
        nb = ((x + dx) % side) * side + (y + dy) % side
        keep = active & active[nb]
        rows.append(np.flatnonzero(keep)); cols.append(nb[keep])
    if lattice == "honeycomb":
        # the vertical bond goes up from even sites; side must be even for periodicity
        up = (x + y) % 2 == 0
        nb = x * side + (y + 1) % side
        keep = active & up & active[nb]
        rows.append(np.flatnonzero(keep)); cols.append(nb[keep])
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    idx = np.flatnonzero(active)
    local = -np.ones(n, dtype=int); local[idx] = np.arange(len(idx))
    graph = coo_matrix((np.ones(len(rows)), (local[rows], local[cols])),
                       shape=(len(idx), len(idx)))
    _, labels = connected_components(graph, directed=False)
    return np.bincount(labels), len(idx)


def sampled_population(lattice, p_values, side=512, reps=8, seed=7):
    """Mass fractions of active sites in components of size 1, 2 and 3 or more."""
    rng = np.random.default_rng(seed)
    rows = []
    for p in p_values:
        fr = []
        for _ in range(reps):
            sizes, n_active = sample_lattice(lattice, side, p, rng)
            mass = np.bincount(sizes, weights=sizes, minlength=4)
            fr.append([mass[1] / n_active, mass[2] / n_active,
                       1 - (mass[1] + mass[2]) / n_active])
        fr = np.array(fr)
        rows.append(dict(p=float(p), mean=fr.mean(0).tolist(),
                         sem=(fr.std(0, ddof=1) / np.sqrt(reps)).tolist()))
    return rows


# --------------------------------------------------------------------- rates

def isomorphism_classes(lattice, animals, s_min=3):
    """Group the animals of each size by graph isomorphism and keep one representative."""
    classes = []
    for s, group in animals.items():
        if s < s_min:
            continue
        buckets = {}
        for a in group:
            size, edges = animal_graph(lattice, a)
            g = nx.Graph(edges)
            key = (s, len(edges), tuple(sorted(d for _, d in g.degree())))
            for cls in buckets.setdefault(key, []):
                if nx.is_isomorphic(cls["graph"], g):
                    cls["animals"].append(a)
                    break
            else:
                buckets[key].append(dict(size=s, edges=edges, graph=g, animals=[a]))
        for bucket in buckets.values():
            classes.extend(bucket)
    return classes


def rate_table(lattice, classes, rates):
    """Exact stationary product rate of each component shape, one solve per shape."""
    cell = LATTICES[lattice]["cell"]
    rows = []
    for c in classes:
        out = solve(c["size"], c["edges"], rates)
        perims = {}
        for a in c["animals"]:
            t = perimeter(lattice, a)
            perims[t] = perims.get(t, 0) + 1
        rows.append(dict(size=c["size"], edges=c["edges"], R=out["tof"] * c["size"],
                         states=out["states"], residual=out["residual"],
                         poly=[[n / cell, t] for t, n in sorted(perims.items())]))
    return rows


def tof_series(table, p_values):
    """Contribution of each component size to the rate per lattice site."""
    p = np.asarray(p_values, dtype=float)
    by_size = {}
    for row in table:
        w = sum(c * p ** row["size"] * (1 - p) ** t for c, t in row["poly"])
        by_size[row["size"]] = by_size.get(row["size"], 0.0) + w * row["R"]
    return by_size


def illustrative_rates():
    """The 745 K rate set of the first pass geometry study, with the phase frozen."""
    from kmc_core import (P_ADS_SO2, P_DES_SO2, P_DES_SO3, P_ADS_SO3,
                          B_ADS_O2, B_DES_O2, B_LH)
    from params import build_rates
    k = build_rates(745.0, 1000e-6, 0.05, 1e-9, params=dict(Ea_sulfate=1.45, Ea_desulfate=1.55))
    r = dict(BASE_RATES)
    r.update(k_ads_so2=k[("site", P_ADS_SO2)], k_des_so2=k[("site", P_DES_SO2)],
             k_des_so3=k[("site", P_DES_SO3)], k_ads_so3=k[("site", P_ADS_SO3)],
             k_ads_o2=k[("bond", B_ADS_O2)], k_des_o2=k[("bond", B_DES_O2)],
             k_lh=k[("bond", B_LH)])
    return r


SYNTHETIC = dict(BASE_RATES, k_ads_so2=1.0, k_des_so2=1.0, k_des_so3=1.0,
                 k_ads_o2=1.0, k_des_o2=0.1, k_lh=1.0)


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--count-max", type=int, default=8,
                    help="largest animal enumerated for counts and the sum rule")
    ap.add_argument("--rate-max", type=int, default=6,
                    help="largest component solved exactly for the rate expansion")
    ap.add_argument("--part", choices=("all", "population", "expansion"), default="all")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.part in ("all", "population"):
        population_part(args)
    if args.part in ("all", "expansion"):
        expansion_part(args)


def population_part(args):
    print("Lattice animals and the dead fraction law")
    grid = np.round(np.arange(0.02, 0.99, 0.02), 2)
    population = {}
    for lattice in LATTICES:
        s_max = min(args.count_max, len(OEIS[lattice]))
        animals, poly = perimeter_polynomials(lattice, s_max)
        counts = [len(animals[s]) for s in range(1, s_max + 1)]
        if counts != OEIS[lattice][:s_max]:
            raise AssertionError(f"{lattice} animal counts {counts} disagree with OEIS")
        coeffs = series_sum_rule(poly)
        worst_low = float(np.max(np.abs(coeffs[:s_max + 1])))
        if worst_low > 1e-9:
            raise AssertionError(f"{lattice} sum rule fails below order {s_max + 1}")
        sampled = sampled_population(lattice, grid)
        n = cluster_numbers(poly, grid)
        exact_single = n[1] / grid
        exact_pair = 2 * n[2] / grid
        closed = dead_fraction(lattice, grid)
        dev = max(abs(exact_single[i] + exact_pair[i] - closed[i]) for i in range(len(grid)))
        # Standard errors from eight replicas can come out as zero when the dead
        # fraction is tiny, so each one is floored at the counting error expected
        # for that many dead sites, with pairs counted as correlated doublets.
        n_total = 8 * 512 * 512 * grid
        floor = np.sqrt(2.0 * closed / n_total)
        z_scores = [abs((row["mean"][0] + row["mean"][1]) - closed[i])
                    / max(np.hypot(row["sem"][0], row["sem"][1]), floor[i])
                    for i, row in enumerate(sampled)]
        population[lattice] = dict(
            z=LATTICES[lattice]["z"], counts=counts, oeis=OEIS[lattice][:s_max],
            sum_rule_low_order_max=worst_low, sum_rule_first_nonzero_order=int(
                next(k for k in range(len(coeffs)) if abs(coeffs[k]) > 1e-9)),
            closed_form_vs_expansion_max=float(dev),
            p=grid.tolist(), dead_closed_form=closed.tolist(),
            dead_lattice_oxygen=dead_fraction(lattice, grid, "lattice-oxygen").tolist(),
            sampled=sampled, max_abs_z=float(max(z_scores)),
            dead_at_p_0_1=float(dead_fraction(lattice, 0.1)))
        print(f"  {lattice:10s} z={LATTICES[lattice]['z']}  counts {counts} match OEIS; "
              f"sum rule exact to order p^{population[lattice]['sum_rule_first_nonzero_order'] - 1}; "
              f"dead at p=0.1 {population[lattice]['dead_at_p_0_1']:.5f}; "
              f"sampling max |z| {population[lattice]['max_abs_z']:.2f}")
    (OUT / "population.json").write_text(json.dumps(population, indent=1))
    print(f"  wrote {(OUT / 'population.json').relative_to(ROOT)}")


def expansion_part(args):

    print(f"\nExact rate expansion on the square lattice, components up to {args.rate_max} sites")
    animals, poly = perimeter_polynomials("square", args.rate_max)
    classes = isomorphism_classes("square", animals)
    print(f"  {sum(len(animals[s]) for s in animals if s >= 3)} animals of 3 or more sites "
          f"fall into {len(classes)} graph isomorphism classes")
    expansion = {}
    for name, rates in (("illustrative_745K", illustrative_rates()), ("synthetic", SYNTHETIC)):
        table = rate_table("square", classes, rates)
        p_values = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
        contrib = tof_series(table, p_values)
        n = cluster_numbers(poly, p_values)
        covered = sum(s * n[s] for s in n) / p_values
        expansion[name] = dict(
            rates=rates, shapes=table, p=p_values.tolist(),
            contribution_by_size={str(s): v.tolist() for s, v in contrib.items()},
            tof_truncated=sum(contrib.values()).tolist(),
            mass_fraction_beyond=(1 - covered).tolist(),
            max_residual=max(r["residual"] for r in table))
        print(f"  {name}: {len(table)} shapes solved, largest residual "
              f"{expansion[name]['max_residual']:.1e}")
        for i, p in enumerate(p_values):
            parts = "  ".join(f"s={s}: {contrib[s][i]:.3e}" for s in sorted(contrib))
            print(f"    p={p:.2f}  TOF per site {sum(contrib.values())[i]:.4e}   {parts}"
                  f"   mass beyond {1 - covered[i]:.2e}")
    (OUT / "expansion.json").write_text(json.dumps(expansion, indent=1))
    print(f"\n  wrote {(OUT / 'expansion.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
