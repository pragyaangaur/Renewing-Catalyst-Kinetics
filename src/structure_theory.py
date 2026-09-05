"""Rate-independent structure results for the dual-site oxidation model.

Two statements are established here. Both hold for every strictly positive choice
of rate constants, so neither depends on the illustrative barriers used elsewhere
in this repository. That is the point of the module. Results that survive any
rate assignment are not weakened by the fact that the barriers are unsourced.

Result A, the three-site theorem. On a connected patch of active sites, the
stationary rate of product formation is exactly zero when the patch has one or two
sites, and strictly positive when it has three or more. This is proved below and
checked by exhaustive enumeration of every connected graph up to five sites, and up
to six sites when the module is run with --max-sites 6.

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

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

EMPTY, O, SO2, SO3, SO4, BARE = 0, 1, 2, 3, 4, 5
ROOT = Path(__file__).resolve().parents[1]

# Oxygen atoms carried by a site in each species, used for the balance check.
OXYGEN_CONTENT = {EMPTY: 0, O: 1, SO2: 2, SO3: 3, SO4: 4, BARE: 0}
SULFUR_CONTENT = {EMPTY: 0, O: 0, SO2: 1, SO3: 1, SO4: 1, BARE: 0}

BASE_RATES = dict(
    k_ads_so2=3.7e4, k_des_so2=1.7e6, k_des_so3=1.7e5, k_ads_o2=4.1e3,
    k_des_o2=0.3, k_lh=7.9e5, k_sulf=0.0, k_desulf=0.0, k_ox=3.9e4, gamma=0.0,
    k_diff=0.0,
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

    The algorithm only ever adds and multiplies positive quantities, so it cannot
    produce a negative probability through cancellation. That makes it the reliable
    choice when rate constants span many orders of magnitude, which is exactly the
    regime the random parameter draws explore. It costs more than a direct solve,
    so it is used as a fallback rather than as the default.
    """
    a = np.array(generator.todense(), dtype=float)
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


def solve(size, edges, r):
    """Stationary solution of the master equation on one patch."""
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

    _, labels = connected_components(gen, directed=True, connection="strong")
    closed = set(labels)
    src, dst = gen.nonzero()
    for a, b in zip(src, dst):
        if labels[a] != labels[b]:
            closed.discard(labels[a])
    if len(closed) != 1:
        raise ValueError("expected exactly one reachable closed class")

    # Rate constants here span many orders of magnitude, which makes the balance
    # matrix badly conditioned. These state spaces are small, so a dense least
    # squares solve is used. It is far more robust than a sparse direct solve and
    # the cost is irrelevant at this size.
    scale = max(float(np.max(-gen.diagonal())), 1.0)
    bal = np.asarray((gen.T / scale).todense())
    bal = np.vstack([bal, np.ones((1, n))])
    rhs = np.zeros(n + 1)
    rhs[-1] = 1.0
    p, *_ = np.linalg.lstsq(bal, rhs, rcond=None)

    # Stationary distributions here can be extremely skewed, so rounding noise on
    # the near zero states is judged relative to the largest probability. When the
    # direct solve does worse than that, fall back to the stable reduction.
    def acceptable(vec):
        if vec is None or vec.sum() <= 0:
            return None
        peak = float(np.max(np.abs(vec)))
        if np.min(vec) < -1e-6 * peak:
            return None
        out = np.clip(vec, 0.0, None)
        out = out / out.sum()
        if float(np.max(np.abs(gen.T @ out))) / scale > 1e-8:
            return None
        return out

    good = acceptable(p)
    if good is None:
        good = acceptable(gth_solve(gen))
    if good is None:
        raise ArithmeticError("both stationary solvers failed their checks")
    p = good
    residual = float(np.max(np.abs(gen.T @ p))) / scale

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


def productive_state_reachable(size, edges):
    """
    Exact test of whether any configuration holding SO3 can be reached from the
    empty patch, using integer logic only.

    Every rate is treated as merely possible rather than as a number, so the answer
    depends on the patch shape alone. This is the honest form of the theorem check.
    Floating point turnover values can fall to solver precision when rate constants
    span many orders of magnitude, which says nothing about whether the rate is
    truly zero.
    """
    start = (EMPTY,) * size
    seen = {start}
    stack = [start]
    while stack:
        state = stack.pop()
        if SO3 in state:
            return True
        for i, sp in enumerate(state):
            if sp == EMPTY:
                nxt = list(state); nxt[i] = SO2
                nxt = tuple(nxt)
                if nxt not in seen:
                    seen.add(nxt); stack.append(nxt)
            elif sp == SO2:
                nxt = list(state); nxt[i] = EMPTY
                nxt = tuple(nxt)
                if nxt not in seen:
                    seen.add(nxt); stack.append(nxt)
        for a, b in edges:
            pair = (state[a], state[b])
            moves = []
            if pair == (EMPTY, EMPTY):
                moves.append(((a, O), (b, O)))
            elif pair == (O, O):
                moves.append(((a, EMPTY), (b, EMPTY)))
            elif pair in ((SO2, O), (O, SO2)):
                sx, ox = (a, b) if pair[0] == SO2 else (b, a)
                moves.append(((sx, SO3), (ox, EMPTY)))
            for mv in moves:
                nxt = list(state)
                for site, species in mv:
                    nxt[site] = species
                nxt = tuple(nxt)
                if nxt not in seen:
                    seen.add(nxt); stack.append(nxt)
    return False


# ------------------------------------------------------------------ result A

def result_a(max_n=6, draws=8, seed=0):
    """
    The three-site theorem, checked exhaustively.

    Proof. Product can only leave a site that holds SO3, and the only way to make
    SO3 is the reaction between an adsorbed SO2 and an adsorbed O on neighbouring
    sites. So a productive configuration needs SO2 and O present at the same time.

    On a single site there is no neighbour, so the reaction can never fire.

    On two sites the only source of adsorbed oxygen is O2, which lands on a pair of
    empty neighbours and fills both. After it lands the patch holds O on both sites,
    so there is no empty site left for SO2 to occupy. Removing an O requires the
    reverse step, which takes both O atoms away together. Starting from the empty
    patch, no reachable configuration holds an O and an SO2 at the same time, so the
    stationary product rate is exactly zero.

    On three or more connected sites there is always a path a-b-c. Oxygen adsorbs on
    the pair (a, b), SO2 then adsorbs on the still empty c, and the reaction fires
    across (b, c). That sequence is reachable from the empty patch and has positive
    probability, so the stationary rate is strictly positive.

    Every connected graph on three or more vertices contains such a path, because a
    connected graph with no vertex of degree two or more is a single edge at most.
    The result therefore depends on the patch having three sites and not on any rate.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for n in range(1, max_n + 1):
        graphs = connected_graphs(n) if n > 1 else [[]]
        for edges in graphs:
            reach = productive_state_reachable(n, edges)
            if reach != (n >= 3):
                raise AssertionError(f"three-site theorem failed on n={n} edges={edges}")
            for d in range(draws):
                r = dict(BASE_RATES)
                if d:  # first draw uses the illustrative rates, the rest are random
                    for key in ("k_ads_so2", "k_des_so2", "k_des_so3",
                                "k_ads_o2", "k_des_o2", "k_lh"):
                        r[key] = float(10.0 ** rng.uniform(-2, 7))
                    r["k_diff"] = float(10.0 ** rng.uniform(-2, 6)) if d % 2 else 0.0
                out = solve(n, edges, r)
                rows.append({"n": n, "edges": len(edges), "draw": d,
                             "tof": out["tof"], "predicted_zero": n < 3,
                             "reachable": reach})
    return rows


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
    # six sites is exhaustively checkable but the isomorphism reduction is slow in
    # pure Python, so the default stops at five and six is opt in
    ap.add_argument("--max-sites", type=int, default=5)
    ap.add_argument("--draws", type=int, default=8)
    args = ap.parse_args()

    print("Result A: the three-site theorem")
    rows_a = result_a(args.max_sites, args.draws)
    # the theorem is a statement about reachability, so that is what is checked
    bad = [r for r in rows_a if r["reachable"] == r["predicted_zero"]]
    # a patch predicted to be dead must give exactly zero in the numerics too
    hard_zero = [r for r in rows_a if r["predicted_zero"] and abs(r["tof"]) > 1e-14]
    by_n = {}
    for r in rows_a:
        by_n.setdefault(r["n"], []).append(r)
    print(f"  {'sites':>6} {'graphs':>7} {'cases':>7} {'min TOF':>12} {'max TOF':>12}")
    for n in sorted(by_n):
        rs = by_n[n]
        g = len({r["edges"] for r in rs}) if n < 3 else len(connected_graphs(n))
        tofs = [r["tof"] for r in rs]
        print(f"  {n:6d} {g:7d} {len(rs):7d} {min(tofs):12.4g} {max(tofs):12.4g}")
    print(f"  exact reachability cases: {len(rows_a)}, violations: {len(bad)}")
    print(f"  patches predicted dead that produced anything: {len(hard_zero)}")
    tiny = [r for r in rows_a if not r["predicted_zero"] and abs(r["tof"]) < 1e-9]
    print(f"  live patches whose rate fell to solver precision: {len(tiny)}")
    print("  those are extreme random rate draws, and the reachability test above")
    print("  confirms their true rate is positive.")

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
    (out / "result_b.json").write_text(json.dumps(rows_b, indent=1))
    print(f"\n  wrote {out.relative_to(ROOT)}/result_a.json and result_b.json")


if __name__ == "__main__":
    main()
