"""Rate-independent structure results for the dual-site oxidation model.

Two statements are established here. Both hold for every strictly positive choice
of rate constants, so neither depends on the illustrative barriers used elsewhere
in this repository. That is the point of the module. Results that survive any
rate assignment are not weakened by the fact that the barriers are unsourced.

Result A, the three-site theorem. On a connected patch of active sites, the
stationary rate of product formation is exactly zero when the patch has one or two
sites, and strictly positive when it has three or more. This is proved below and
checked by exhaustive enumeration. The reachability and irreducibility checks run
on all 996 connected graphs with one to seven sites, with and without diffusion.
Stationary rates are solved numerically, for randomised rate constants, on every
connected graph with up to five sites.

Result B, the copper efficiency identity with regeneration. When sulphated sites
can be recovered chemically as well as by discarding the layer, the product
obtained per copper atom exposed is

    Y_Cu = (1 + k_d/k_s) * rho + (k_d/k_s) * f_sulphate

where k_d and k_s are the product-release and sulphation rates from the shared
adsorbed intermediate, f_sulphate is the mean sulphated fraction, and rho is the
number of chemical regeneration events per copper atom exposed. Setting rho to
zero recovers the renewal-only ceiling reported in the follow-up study. A nonzero
rho lifts the ceiling without bound, which identifies chemical regeneration as the
only route past it.

The species codes are shared with kmc_core so that the two models stay comparable.
"""

import argparse
import itertools
import json
from pathlib import Path

import networkx as nx
import numpy as np
from networkx.generators.atlas import graph_atlas_g
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

EMPTY, O, SO2, SO3, SO4, BARE = 0, 1, 2, 3, 4, 5
ROOT = Path(__file__).resolve().parents[1]

# Oxygen atoms carried by a site in each species, used for the balance check.
OXYGEN_CONTENT = {EMPTY: 0, O: 1, SO2: 2, SO3: 3, SO4: 4, BARE: 0}
SULFUR_CONTENT = {EMPTY: 0, O: 0, SO2: 1, SO3: 1, SO4: 1, BARE: 0}

BASE_RATES = dict(
    k_ads_so2=3.7e4, k_des_so2=1.7e6, k_des_so3=1.7e5, k_ads_o2=4.1e3,
    k_des_o2=0.3, k_lh=7.9e5, k_sulf=0.0, k_desulf=0.0, k_ox=3.9e4, gamma=0.0,
    k_diff=0.0, k_ads_so3=0.0,
)


# --------------------------------------------------------------------- dynamics

def transitions(state, edges, r):
    """Every elementary move out of one configuration, with its accounting."""
    out = []

    def add(changes, rate, product=0.0, regen=0.0, sulfation=0.0):
        if rate <= 0.0:
            return
        nxt = list(state)
        for site, species in changes:
            nxt[site] = species
        out.append((tuple(nxt), rate, product, regen, sulfation))

    for i, sp in enumerate(state):
        if sp == EMPTY:
            add([(i, SO2)], r["k_ads_so2"])
            # readsorption of product from the gas, counted as negative product so
            # that the product rate is the net rate the lattice KMC reports
            add([(i, SO3)], r.get("k_ads_so3", 0.0), product=-1.0)
        elif sp == SO2:
            add([(i, EMPTY)], r["k_des_so2"])
        elif sp == SO3:
            add([(i, EMPTY)], r["k_des_so3"], product=1.0)
            add([(i, SO4)], r["k_sulf"], sulfation=1.0)
        elif sp == SO4:
            # chemical regeneration: the sulphate decomposes and the SO3 leaves,
            # so the site returns to service without spending a copper atom
            add([(i, EMPTY)], r["k_desulf"], product=1.0, regen=1.0)
        elif sp == BARE:
            add([(i, EMPTY)], r["k_ox"])

    for a, b in edges:
        pair = (state[a], state[b])
        if pair == (EMPTY, EMPTY):
            add([(a, O), (b, O)], r["k_ads_o2"])
        elif pair == (O, O):
            add([(a, EMPTY), (b, EMPTY)], r["k_des_o2"])
        elif pair in ((SO2, O), (O, SO2)):
            s, o = (a, b) if pair[0] == SO2 else (b, a)
            add([(s, SO3), (o, EMPTY)], r["k_lh"])
        if r["k_diff"] > 0 and max(pair) < SO4 and ((pair[0] == EMPTY) != (pair[1] == EMPTY)):
            add([(a, pair[1]), (b, pair[0])], r["k_diff"])

    if r["gamma"] > 0.0:
        # the whole layer is discarded and replaced by fresh metallic copper
        add([(i, BARE) for i in range(len(state))], r["gamma"])
    return out


def gth_solve(generator):
    """
    Stationary distribution by the Grassmann, Taksar and Heyman reduction.

    The algorithm never subtracts, so it cannot lose a small probability through
    cancellation, and each entry of the result carries a small relative error no
    matter how widely the rate constants are spread (O'Cinneide, Numerische
    Mathematik 65 (1993) 109). That is the property the random parameter draws
    need. A least squares solve was the default here until September 2026 and was
    found to be wrong by up to twelve orders of magnitude on the smallest rates,
    which docs/CORRECTIONS.md records. The generator must be irreducible.
    """
    a = np.array(generator.todense() if hasattr(generator, "todense") else generator,
                 dtype=float)
    n = a.shape[0]
    for k in range(n - 1):
        scale = a[k, k + 1:].sum()
        if scale <= 0.0:
            raise ArithmeticError("GTH reduction met a non positive pivot")
        a[k + 1:, k] /= scale
        a[k + 1:, k + 1:] += np.outer(a[k + 1:, k], a[k, k + 1:])
    x = np.zeros(n)
    x[n - 1] = 1.0
    for k in range(n - 2, -1, -1):
        x[k] = x[k + 1:] @ a[k + 1:, k]
    return x / x.sum()


def build_chain(size, edges, r):
    """
    The master equation on one patch, restricted to the configurations reachable
    from its starting state. The start is the empty patch, or bare copper when the
    patch renews. Returns the states, the generator with rows as the source state,
    and the per-state rates of net product release, regeneration and sulphation.
    """
    start = (BARE,) * size if r["gamma"] > 0.0 else (EMPTY,) * size
    states = [start]
    index = {start: 0}
    rows, cols, vals = [], [], []
    prod, regen, sulf = [], [], []

    i = 0
    while i < len(states):
        state = states[i]
        total = p_i = g_i = s_i = 0.0
        for nxt, rate, p, g, s in transitions(state, edges, r):
            if nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
            rows.append(i); cols.append(index[nxt]); vals.append(rate)
            total += rate
            p_i += rate * p
            g_i += rate * g
            s_i += rate * s
        rows.append(i); cols.append(i); vals.append(-total)
        prod.append(p_i); regen.append(g_i); sulf.append(s_i)
        i += 1

    n = len(states)
    gen = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return dict(states=states, gen=gen, prod=np.array(prod), regen=np.array(regen),
                sulf=np.array(sulf))


def solve(size, edges, r):
    """Stationary solution of the master equation on one patch."""
    chain = build_chain(size, edges, r)
    states, gen = chain["states"], chain["gen"]
    prod, regen, sulf = chain["prod"], chain["regen"], chain["sulf"]
    n = len(states)

    _, labels = connected_components(gen, directed=True, connection="strong")
    closed = set(labels)
    src, dst = gen.nonzero()
    for a, b in zip(src, dst):
        if labels[a] != labels[b]:
            closed.discard(labels[a])
    if len(closed) != 1:
        raise ValueError("expected exactly one reachable closed class")

    # States outside the closed class are transient and carry zero stationary
    # probability, so the stationary solve runs on the closed class alone. That
    # keeps the generator irreducible, which the GTH reduction requires.
    keep = np.flatnonzero(labels == closed.pop())
    sub = gen[keep][:, keep].tolil()
    sub.setdiag(0.0)
    sub = sub.tocsr()
    sub = sub - coo_matrix((np.asarray(sub.sum(axis=1)).ravel(),
                            (np.arange(len(keep)), np.arange(len(keep)))),
                           shape=sub.shape).tocsr()
    p = np.zeros(n)
    p[keep] = gth_solve(sub)

    scale = max(float(np.max(-gen.diagonal())), 1.0)
    residual = float(np.max(np.abs(gen.T @ p))) / scale
    if np.min(p) < 0.0 or residual > 1e-12:
        raise ArithmeticError("stationary solution failed its checks")

    f_sulf = float(p @ [s.count(SO4) / size for s in states])
    return {
        "R_product": float(p @ prod),
        "R_regen": float(p @ regen),
        "R_sulfation": float(p @ sulf),
        "f_sulfate": f_sulf,
        "tof": float(p @ prod) / size,
        "states": n,
        "residual": residual,
    }


# ------------------------------------------------------- graph enumeration

def connected_graphs(n):
    """
    Every connected graph on n vertices up to isomorphism, as an edge list.

    Taken from the Atlas of Graphs of Read and Wilson, which networkx ships as
    graph_atlas_g and which is complete up to seven vertices. The counts for one to
    seven vertices are 1, 1, 2, 6, 21, 112 and 853, which is OEIS A001349. The
    brute force enumerator below reproduces the atlas up to five vertices, and the
    test suite checks that the two agree.
    """
    if not 1 <= n <= 7:
        raise ValueError("the graph atlas is complete only up to seven vertices")
    return [sorted(tuple(sorted(e)) for e in g.edges())
            for g in graph_atlas_g() if g.number_of_nodes() == n and nx.is_connected(g)]


def connected_graphs_bruteforce(n):
    """Every connected graph on n labelled vertices, deduplicated by isomorphism."""
    verts = list(range(n))
    all_edges = list(itertools.combinations(verts, 2))
    seen = {}
    for mask in range(1 << len(all_edges)):
        edges = [all_edges[i] for i in range(len(all_edges)) if mask >> i & 1]
        if len(edges) < n - 1:
            continue
        adj = {v: set() for v in verts}
        for a, b in edges:
            adj[a].add(b); adj[b].add(a)
        stack, seen_v = [0], {0}
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if w not in seen_v:
                    seen_v.add(w); stack.append(w)
        if len(seen_v) != n:
            continue
        canon = min(
            tuple(sorted(tuple(sorted((perm[a], perm[b]))) for a, b in edges))
            for perm in map(lambda p: dict(zip(verts, p)), itertools.permutations(verts))
        )
        seen.setdefault(canon, edges)
    return list(seen.values())


def adsorbate_moves(state, edges, diffusion=False):
    """
    Every configuration one elementary step away, with rates treated as merely
    possible. This is the adsorbate network of Result A: SO2 adsorption and
    desorption, dissociative O2 adsorption onto an empty pair, recombinative O2
    desorption, the surface reaction and SO3 desorption, and optionally the hop of
    an adsorbate onto an empty neighbour.
    """
    out = []
    for i, sp in enumerate(state):
        if sp == EMPTY:
            out.append(state[:i] + (SO2,) + state[i + 1:])
        elif sp in (SO2, SO3):
            out.append(state[:i] + (EMPTY,) + state[i + 1:])
    for a, b in edges:
        pair = (state[a], state[b])
        changes = []
        if pair == (EMPTY, EMPTY):
            changes.append(((a, O), (b, O)))
        elif pair == (O, O):
            changes.append(((a, EMPTY), (b, EMPTY)))
        elif pair in ((SO2, O), (O, SO2)):
            sx, ox = (a, b) if pair[0] == SO2 else (b, a)
            changes.append(((sx, SO3), (ox, EMPTY)))
        if diffusion and (pair[0] == EMPTY) != (pair[1] == EMPTY):
            changes.append(((a, pair[1]), (b, pair[0])))
        for mv in changes:
            nxt = list(state)
            for site, species in mv:
                nxt[site] = species
            out.append(tuple(nxt))
    return out


def reachability(size, edges, diffusion=False):
    """
    Exact integer test of the two facts the proof of Result A rests on.

    productive   some configuration holding SO3 is reachable from the empty patch
    irreducible  the empty patch can be reached again from every reachable
                 configuration, so the reachable set is one communicating class

    When both hold, the chain restricted to the reachable set is finite and
    irreducible, so every reachable configuration has strictly positive stationary
    probability, and the product rate is strictly positive for every positive rate
    assignment. When productive is false the product rate is exactly zero. Neither
    test uses a rate constant, so the answer depends on the patch shape alone.
    Floating point rates can fall to solver precision when rate constants span many
    orders of magnitude, which says nothing about whether the true rate is zero.
    """
    start = (EMPTY,) * size
    index = {start: 0}
    states = [start]
    succ = []
    i = 0
    while i < len(states):
        nexts = []
        for nxt in adsorbate_moves(states[i], edges, diffusion):
            if nxt not in index:
                index[nxt] = len(states)
                states.append(nxt)
            nexts.append(index[nxt])
        succ.append(nexts)
        i += 1
    pred = [[] for _ in states]
    for u, nexts in enumerate(succ):
        for v in nexts:
            pred[v].append(u)
    back = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u in pred[v]:
            if u not in back:
                back.add(u)
                stack.append(u)
    return {
        "productive": any(SO3 in st for st in states),
        "irreducible": len(back) == len(states),
        "states": len(states),
    }


def productive_state_reachable(size, edges, diffusion=False):
    """Whether any configuration holding SO3 is reachable from the empty patch."""
    return reachability(size, edges, diffusion)["productive"]


# ------------------------------------------------------------------ result A

def result_a(max_reach=7, max_numeric=5, draws=8, seed=0):
    """
    The three-site theorem, checked exhaustively.

    Scope. The adsorbate network of adsorbate_moves on a connected patch whose
    neighbours are permanently inactive, with every rate strictly positive and
    diffusion optional. Sulphation is switched off, since an irreversible sink with
    no recovery would make every patch dead in the long run for a trivial reason.

    Proof. Product can only leave a site that holds SO3, and the only way to make
    SO3 is the reaction between an adsorbed SO2 and an adsorbed O on neighbouring
    sites. So a productive configuration needs SO2 and O present at the same time.

    On a single site there is no neighbour, so the reaction can never fire.

    On two sites the only source of adsorbed oxygen is O2, which lands on a pair of
    empty neighbours and fills both. After it lands the patch holds O on both sites,
    so there is no empty site left for SO2 to occupy, and a hop needs an empty site
    too. Removing an O requires the reverse step, which takes both O atoms away
    together. Starting from the empty patch, no reachable configuration holds an O
    and an SO2 at the same time, so the stationary product rate is exactly zero.

    On three or more connected sites there is always a path a-b-c. Oxygen adsorbs on
    the pair (a, b), SO2 then adsorbs on the still empty c, and the reaction fires
    across (b, c). That sequence is reachable from the empty patch. The empty patch
    is also reachable from every configuration. SO2 and SO3 can always leave. Two
    neighbouring O atoms can leave together. An O atom with no O neighbour has some
    neighbour, which can be emptied and then filled with SO2 so the reaction removes
    the O. The reachable set is therefore a single communicating class, every state
    in it has positive stationary probability, and the product rate is strictly
    positive.

    Every connected graph on three or more vertices contains such a path, because a
    connected graph with no vertex of degree two or more is a single edge at most.
    The result therefore depends on the patch having three sites and not on any rate.

    Returns the numerical rows, for graphs up to max_numeric sites, and one
    reachability row per graph up to max_reach sites.
    """
    rng = np.random.default_rng(seed)
    rows, reach_rows = [], []
    for n in range(1, max_reach + 1):
        graphs = connected_graphs(n) if n > 1 else [[]]
        for g, edges in enumerate(graphs):
            plain = reachability(n, edges, diffusion=False)
            hop = reachability(n, edges, diffusion=True)
            for tag, res in (("no diffusion", plain), ("diffusion", hop)):
                if res["productive"] != (n >= 3) or not res["irreducible"]:
                    raise AssertionError(
                        f"three-site theorem failed on n={n} edges={edges} ({tag})")
            degree = np.bincount(np.array(edges, dtype=int).ravel(), minlength=n)
            reach_rows.append({"n": n, "graph": g, "edges": len(edges),
                               "max_degree": int(degree.max()) if n > 1 else 0,
                               "productive": plain["productive"],
                               "irreducible": plain["irreducible"] and hop["irreducible"],
                               "productive_with_diffusion": hop["productive"],
                               "states": plain["states"],
                               "states_with_diffusion": hop["states"]})
            if n > max_numeric:
                continue
            for d in range(draws):
                r = dict(BASE_RATES)
                if d:  # first draw uses the illustrative rates, the rest are random
                    for key in ("k_ads_so2", "k_des_so2", "k_des_so3",
                                "k_ads_o2", "k_des_o2", "k_lh"):
                        r[key] = float(10.0 ** rng.uniform(-2, 7))
                    r["k_diff"] = float(10.0 ** rng.uniform(-2, 6)) if d % 2 else 0.0
                out = solve(n, edges, r)
                rows.append({"n": n, "graph": g, "edges": len(edges), "draw": d,
                             "tof": out["tof"], "residual": out["residual"],
                             "states": out["states"], "predicted_zero": n < 3,
                             "reachable": plain["productive"],
                             "rates": {k: v for k, v in r.items()
                                       if k not in ("k_sulf", "k_desulf", "k_ox", "gamma")}})
    return rows, reach_rows


# ------------------------------------------------------------------ result B

def result_b(draws=200, seed=1):
    """
    The copper efficiency identity, including a chemical regeneration channel.

    Derivation. Product leaves either by desorbing from the adsorbed SO3 or by the
    decomposition of a sulphated site. Both sulphation and product desorption draw
    on the same adsorbed SO3 population, so their rates stand in the fixed ratio

        R_product_direct / R_sulfation = k_d / k_s

    At stationary state the sulphate created must equal the sulphate destroyed, and
    sulphate is destroyed either by chemical regeneration or by discarding the layer

        R_sulfation = R_regen + gamma * n * f_sulfate

    Writing rho for the regenerations per copper atom exposed, R_regen / (gamma * n),
    and dividing the total product by the copper spend gamma * n gives

        Y_Cu = (1 + k_d/k_s) * rho + (k_d/k_s) * f_sulfate

    With rho set to zero this is the renewal-only ceiling, which cannot exceed the
    branching ratio k_d/k_s no matter how the patch is shaped or how often it resets.
    With regeneration switched on the first term grows without limit as the reset
    frequency falls, so chemical regeneration is the only way past that ceiling.
    """
    rng = np.random.default_rng(seed)
    motifs = {
        "pair": (2, [(0, 1)]),
        "path3": (3, [(0, 1), (1, 2)]),
        "star4": (4, [(0, 1), (0, 2), (0, 3)]),
        "square": (4, [(0, 1), (1, 2), (2, 3), (0, 3)]),
    }
    rows = []
    for d in range(draws):
        name = list(motifs)[d % len(motifs)]
        size, edges = motifs[name]
        r = dict(BASE_RATES)
        r["k_sulf"] = float(10.0 ** rng.uniform(1, 6))
        r["k_des_so3"] = float(10.0 ** rng.uniform(1, 6))
        r["k_desulf"] = 0.0 if d % 3 == 0 else float(10.0 ** rng.uniform(-3, 3))
        r["gamma"] = float(10.0 ** rng.uniform(-2, 3))
        r["k_ox"] = float(10.0 ** rng.uniform(2, 5))
        r["k_diff"] = float(10.0 ** rng.uniform(-1, 5)) if d % 2 else 0.0
        out = solve(size, edges, r)

        spend = r["gamma"] * size
        ratio = r["k_des_so3"] / r["k_sulf"]
        rho = out["R_regen"] / spend
        predicted = (1.0 + ratio) * rho + ratio * out["f_sulfate"]
        actual = out["R_product"] / spend
        rows.append({
            "motif": name, "ratio": ratio, "rho": rho,
            "f_sulfate": out["f_sulfate"], "predicted": predicted, "actual": actual,
            "abs_error": abs(predicted - actual),
            "rel_error": abs(predicted - actual) / max(abs(actual), 1e-30),
            "regen_on": r["k_desulf"] > 0,
        })
    return rows


# ---------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-reach", type=int, default=7,
                    help="largest patch for the exhaustive reachability check")
    ap.add_argument("--max-sites", type=int, default=5,
                    help="largest patch for the numerical stationary solves")
    ap.add_argument("--draws", type=int, default=8)
    args = ap.parse_args()

    print("Result A: the three-site theorem")
    rows_a, reach = result_a(args.max_reach, args.max_sites, args.draws)
    print(f"  {'sites':>6} {'graphs':>7} {'productive':>11} {'irreducible':>12} {'max states':>11}")
    for n in sorted({r["n"] for r in reach}):
        rs = [r for r in reach if r["n"] == n]
        print(f"  {n:6d} {len(rs):7d} {sum(r['productive'] for r in rs):11d} "
              f"{sum(r['irreducible'] for r in rs):12d} "
              f"{max(r['states_with_diffusion'] for r in rs):11d}")
    print(f"  graphs checked exhaustively, with and without diffusion: {len(reach)}")

    # a patch predicted to be dead must give exactly zero in the numerics too
    hard_zero = [r for r in rows_a if r["predicted_zero"] and r["tof"] != 0.0]
    live = [r for r in rows_a if not r["predicted_zero"]]
    print(f"\n  numerical stationary solves: {len(rows_a)}  (live {len(live)})")
    print(f"  patches predicted dead with any nonzero rate: {len(hard_zero)}")
    print(f"  live patches with a rate that is not strictly positive: "
          f"{sum(r['tof'] <= 0.0 for r in live)}")
    print(f"  smallest live rate {min(r['tof'] for r in live):.4g}, largest "
          f"{max(r['tof'] for r in live):.4g} per site per second")
    print(f"  largest relative stationary residual {max(r['residual'] for r in rows_a):.2e}")

    print("\nResult B: copper efficiency identity with regeneration")
    rows_b = result_b()
    err = max(r["rel_error"] for r in rows_b)
    on = [r for r in rows_b if r["regen_on"]]
    off = [r for r in rows_b if not r["regen_on"]]
    print(f"  cases: {len(rows_b)}  (regeneration on {len(on)}, off {len(off)})")
    print(f"  maximum relative error in the identity: {err:.3e}")
    if off:
        worst = max(r["actual"] / r["ratio"] for r in off)
        print(f"  renewal only, max of Y_Cu / (k_d/k_s): {worst:.6f}  (ceiling is 1)")
    if on:
        best = max(r["actual"] / r["ratio"] for r in on)
        print(f"  with regeneration, max of Y_Cu / (k_d/k_s): {best:.4g}  (no ceiling)")

    out = ROOT / "results" / "structure"
    out.mkdir(parents=True, exist_ok=True)
    (out / "result_a.json").write_text(json.dumps(rows_a, indent=1))
    (out / "result_a_reachability.json").write_text(json.dumps(reach, indent=1))
    (out / "result_b.json").write_text(json.dumps(rows_b, indent=1))
    print(f"\n  wrote {out.relative_to(ROOT)}/result_a.json, result_a_reachability.json"
          " and result_b.json")


if __name__ == "__main__":
    main()
