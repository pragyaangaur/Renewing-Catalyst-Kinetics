"""Check the copper efficiency identity against the cumulative site-loss selectivity.

This module answers a question put by Aditya Bhan, Distinguished McKnight
University Professor at the University of Minnesota, who suggested that the
product obtained per copper atom in this repository should be the inverse of the
cumulative site-loss selectivity defined by Foley, Johnson and Bhan, ACS Catalysis
9 (2019) 7065, doi:10.1021/acscatal.9b01106. That paper treats active sites as a
consumable reactant and defines the cumulative site-loss selectivity as the total
moles of sites lost divided by the total moles of reactant converted to effluent
product, so its inverse is the product obtained per site lost.

The question only has an answer once "a site lost" is pinned to an event in this
model, and the model has two candidates.

    S_Cu   = (copper atoms discarded) / (product molecules)
    S_sulf = (sulphation events)      / (product molecules)

The two differ because a reset discards the whole layer, including copper that had
never sulphated. Three relations are checked numerically here.

    1.  Y_Cu = 1 / S_Cu

        Exact and definitional. Counting a discarded copper atom as the lost site
        makes the copper efficiency the inverse cumulative site-loss selectivity,
        which is the correspondence suggested.

    2.  1 / S_sulf = k_d/k_s + rho / (rho + f_sulphate)  * ... see below

        Counting a sulphation event as the lost site instead gives
        1/S_sulf = R_product / R_sulfation, and the identity rearranges to

            Y_Cu = (rho + f_sulphate) / S_sulf

        so the two inverse selectivities differ by the number of sulphation events
        per copper atom discarded. With renewal only, rho is zero, that factor is
        f_sulphate, and 1/S_sulf collapses to the ceiling k_d/k_s exactly.

    3.  1 / S_consumed = Y_Cu + f_sulphate

        The variant of the 2019 definition that puts all reactant consumed in the
        denominator rather than reactant converted to effluent product. Sulphur is
        also consumed by the sulphate that leaves with the discarded layer.

The outcome is that the ceiling k_d/k_s reported in this repository is exactly the
inverse site-loss selectivity of the sulphation channel, and f_sulphate is the
fraction of the discarded copper that had actually deactivated. The gap between
the copper efficiency and the ceiling is therefore not a kinetic effect. It is the
cost of throwing away copper that was still in service.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from structure_theory import BASE_RATES, ROOT, solve

MOTIFS = {
    "pair": (2, [(0, 1)]),
    "path3": (3, [(0, 1), (1, 2)]),
    "star4": (4, [(0, 1), (0, 2), (0, 3)]),
    "square": (4, [(0, 1), (1, 2), (2, 3), (0, 3)]),
}


def rel(a, b):
    """Relative difference, guarded against a zero denominator."""
    return abs(a - b) / max(abs(b), 1e-30)


def run(draws=200, seed=1):
    """
    Draw the same ensemble as result_b in structure_theory and compare the copper
    efficiency with both readings of the cumulative site-loss selectivity.

    The seed and the draw order match result_b on purpose, so the two tables cover
    the same cases and can be read side by side.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(draws):
        name = list(MOTIFS)[d % len(MOTIFS)]
        size, edges = MOTIFS[name]
        r = dict(BASE_RATES)
        r["k_sulf"] = float(10.0 ** rng.uniform(1, 6))
        r["k_des_so3"] = float(10.0 ** rng.uniform(1, 6))
        r["k_desulf"] = 0.0 if d % 3 == 0 else float(10.0 ** rng.uniform(-3, 3))
        r["gamma"] = float(10.0 ** rng.uniform(-2, 3))
        r["k_ox"] = float(10.0 ** rng.uniform(2, 5))
        r["k_diff"] = float(10.0 ** rng.uniform(-1, 5)) if d % 2 else 0.0
        out = solve(size, edges, r)

        spend = r["gamma"] * size          # copper atoms discarded per unit time
        ratio = r["k_des_so3"] / r["k_sulf"]
        rho = out["R_regen"] / spend
        y_cu = out["R_product"] / spend

        # sulphur leaving with the discarded layer, so the reactant consumed but
        # not converted to effluent product
        s_buried = spend * out["f_sulfate"]

        if out["R_product"] <= 0.0:
            # a one or two site patch makes no product at all, by result A, so every
            # site-loss selectivity here is a division by zero and is undefined
            rows.append({"motif": name, "regen_on": r["k_desulf"] > 0, "dead": True})
            continue

        s_cu = spend / out["R_product"]
        s_sulf = out["R_sulfation"] / out["R_product"]
        s_consumed = spend / (out["R_product"] + s_buried)

        rows.append({
            "motif": name,
            "regen_on": r["k_desulf"] > 0,
            "dead": False,
            "ratio": ratio,
            "rho": rho,
            "f_sulfate": out["f_sulfate"],
            "Y_Cu": y_cu,
            "S_cu": s_cu,
            "S_sulf": s_sulf,
            "S_consumed": s_consumed,
            # relation 1, the correspondence as suggested
            "err_inverse_S_cu": rel(y_cu, 1.0 / s_cu),
            # relation 2, the sulphation reading of a lost site
            "err_inverse_S_sulf": rel(y_cu, (rho + out["f_sulfate"]) / s_sulf),
            # with renewal only the sulphation reading is the ceiling itself
            "err_ceiling": rel(1.0 / s_sulf, ratio) if r["k_desulf"] == 0.0 else None,
            # relation 3, the reactant consumed variant of the denominator
            "err_consumed": rel(1.0 / s_consumed, y_cu + out["f_sulfate"]),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--draws", type=int, default=200)
    args = ap.parse_args()

    rows = run(args.draws)
    dead = [r for r in rows if r["dead"]]
    live = [r for r in rows if not r["dead"]]
    off = [r for r in live if not r["regen_on"]]
    on = [r for r in live if r["regen_on"]]

    def worst(key, subset=live):
        vals = [r[key] for r in subset if r[key] is not None]
        return max(vals) if vals else float("nan")

    print("Cumulative site-loss selectivity, Foley, Johnson and Bhan 2019")
    print(f"  cases: {len(live)}  (regeneration on {len(on)}, off {len(off)})")
    print(f"  dead patches excluded, where no product forms and the ratio is "
          f"undefined: {len(dead)}")
    print()
    print("  1. copper atoms discarded read as the sites lost")
    print(f"     Y_Cu = 1 / S_Cu, maximum relative error {worst('err_inverse_S_cu'):.3e}")
    print()
    print("  2. sulphation events read as the sites lost")
    print(f"     Y_Cu = (rho + f_sulphate) / S_sulf, maximum relative error "
          f"{worst('err_inverse_S_sulf'):.3e}")
    print(f"     renewal only, 1 / S_sulf = k_d/k_s, maximum relative error "
          f"{worst('err_ceiling', off):.3e}")
    if off:
        gap = [r["Y_Cu"] * r["S_sulf"] for r in off]
        print(f"     renewal only, Y_Cu / (1 / S_sulf) ranges over "
              f"[{min(gap):.4g}, {max(gap):.4g}], and it equals f_sulphate")
    print()
    print("  3. reactant consumed rather than converted to effluent product")
    print(f"     1 / S_consumed = Y_Cu + f_sulphate, maximum relative error "
          f"{worst('err_consumed'):.3e}")

    out = ROOT / "results" / "structure"
    out.mkdir(parents=True, exist_ok=True)
    (out / "site_loss_check.json").write_text(json.dumps(rows, indent=1))
    print(f"\n  wrote {(out / 'site_loss_check.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
