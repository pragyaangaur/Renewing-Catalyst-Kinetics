"""Exact finite-component master equations and independent Gillespie checks.

Frozen oxide graphs inherit the original adsorbate rates. Only states reachable
from an empty surface are included. Diffusion exchanges an adsorbate and a vacancy
across an oxide edge. This is an idealised model audit, not fitted CuO chemistry.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

from kmc_core import (
    EMPTY, O, SO2, SO3, P_ADS_SO2, P_DES_SO2, P_ADS_SO3, P_DES_SO3,
    B_ADS_O2, B_DES_O2, B_LH, P_SULFATE, P_OXIDISE_CU,
)
from params import build_rates

ROOT = Path(__file__).resolve().parents[1]
MOTIFS = {
    "single": (1, []),
    "pair": (2, [(0, 1)]),
    "path3": (3, [(0, 1), (1, 2)]),
    "triangle": (3, [(0, 1), (1, 2), (0, 2)]),
    "path4": (4, [(0, 1), (1, 2), (2, 3)]),
    "star4": (4, [(0, 1), (0, 2), (0, 3)]),
    "square": (4, [(0, 1), (1, 2), (2, 3), (0, 3)]),
    "path5": (5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
    "square_tail": (5, [(0, 1), (1, 2), (2, 3), (0, 3), (3, 4)]),
    "rectangle6": (6, [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]),
}


def transitions(state, edges, rates, diffusion=0.0, renewal=0.0):
    events = []

    def add(changes, rate, product=0, oxygen=0, sulfur=0):
        if rate <= 0:
            return
        target = list(state)
        for site, species in changes:
            target[site] = species
        events.append((tuple(target), rate, product, oxygen, sulfur))

    for site, species in enumerate(state):
        if species == EMPTY:
            add([(site, SO2)], rates[("site", P_ADS_SO2)], sulfur=1)
            add([(site, SO3)], rates[("site", P_ADS_SO3)], product=-1)
        elif species == SO2:
            add([(site, EMPTY)], rates[("site", P_DES_SO2)], sulfur=-1)
        elif species == SO3:
            add([(site, EMPTY)], rates[("site", P_DES_SO3)], product=1)
            if renewal:
                add([(site, 4)], rates[("site", P_SULFATE)])
        elif species == 5:
            add([(site, EMPTY)], rates[("site", P_OXIDISE_CU)], oxygen=1)
    for left, right in edges:
        pair = state[left], state[right]
        if pair == (EMPTY, EMPTY):
            add([(left, O), (right, O)], rates[("bond", B_ADS_O2)], oxygen=2)
        elif pair == (O, O):
            add([(left, EMPTY), (right, EMPTY)], rates[("bond", B_DES_O2)], oxygen=-2)
        elif pair in ((SO2, O), (O, SO2)):
            sulfur_site, oxygen_site = (left, right) if pair[0] == SO2 else (right, left)
            add([(sulfur_site, SO3), (oxygen_site, EMPTY)], rates[("bond", B_LH)])
        if max(pair) < 4 and ((pair[0] == EMPTY) != (pair[1] == EMPTY)):
            add([(left, pair[1]), (right, pair[0])], diffusion)
    if renewal:
        discarded = sum(species in (SO2, SO3, 4) for species in state)
        extra_oxygen = sum((1, 2, 1, 2, 2, 0)[species] for species in state)
        add([(site, 5) for site in range(len(state))], renewal,
            sulfur=-discarded, oxygen=-extra_oxygen)
    return events


def exact(size, edges, rates, diffusion=0.0, renewal=0.0):
    states = [((5 if renewal else EMPTY),) * size]
    indices = {states[0]: 0}
    rows, columns, values, rewards, oxygen_fluxes, sulfur_fluxes = [], [], [], [], [], []
    for index, state in enumerate(states):
        total = reward = oxygen_flux = sulfur_flux = 0.0
        for target, rate, product, oxygen, sulfur in transitions(state, edges, rates, diffusion, renewal):
            if target not in indices:
                indices[target] = len(states)
                states.append(target)
            rows.append(index)
            columns.append(indices[target])
            values.append(rate)
            total += rate
            reward += rate * product
            oxygen_flux += rate * oxygen
            sulfur_flux += rate * sulfur
        rows.append(index)
        columns.append(index)
        values.append(-total)
        rewards.append(reward)
        oxygen_fluxes.append(oxygen_flux)
        sulfur_fluxes.append(sulfur_flux)
    generator = coo_matrix((values, (rows, columns)), shape=(len(states), len(states))).tocsr()
    _, labels = connected_components(generator, directed=True, connection="strong")
    closed = set(labels)
    source, target = generator.nonzero()
    for start, finish in zip(source, target):
        if labels[start] != labels[finish]:
            closed.discard(labels[start])
    if len(closed) != 1:
        raise ValueError("Expected one reachable closed communicating class")
    balance = generator.T.tolil()
    scale = max(float(np.max(-generator.diagonal())), 1.0)
    balance /= scale
    balance[-1, :] = np.ones(len(states))
    rhs = np.zeros(len(states))
    rhs[-1] = 1.0
    probability = spsolve(balance.tocsc(), rhs)
    residual = float(np.max(np.abs(generator.T @ probability)))
    if np.min(probability) < -1e-8 or residual > scale * 1e-9:
        raise ArithmeticError("Stationary distribution failed numerical checks")
    reward = float(probability @ rewards)
    oxygen_flux = float(probability @ oxygen_fluxes)
    sulfur_flux = float(probability @ sulfur_fluxes)
    if max(abs(oxygen_flux - reward), abs(sulfur_flux - reward)) > 1e-5 * max(abs(reward), 1):
        raise ArithmeticError("Steady-state elemental flux balance failed")
    if abs(sulfur_flux - reward) > 1e-5 * max(abs(reward), 1):
        raise ArithmeticError("Sulfur input must equal gas product plus discarded sulfur")
    return {
        "tof": reward / size, "size": size, "states": len(states),
        "residual": residual, "oxygen_balance_error": oxygen_flux - reward,
        "sulfur_balance_error": sulfur_flux - reward,
        "oxygen_coverage": float(probability @ [state.count(O) / size for state in states]),
        "sulfate_coverage": float(probability @ [state.count(4) / size for state in states]),
        "bare_coverage": float(probability @ [state.count(5) / size for state in states]),
    }


def simulate(size, edges, rates, seed, duration, burn, diffusion=0.0, renewal=0.0):
    rng = np.random.default_rng(seed)
    state = ((5 if renewal else EMPTY),) * size
    time = 0.0
    production = 0
    events = 0
    end = burn + duration
    while time < end:
        options = transitions(state, edges, rates, diffusion, renewal)
        hazards = np.array([event[1] for event in options])
        total = hazards.sum()
        if total <= 0:
            break
        event_time = time + rng.exponential(1.0 / total)
        if event_time > end:
            break
        chosen = int(np.searchsorted(np.cumsum(hazards), rng.random() * total))
        state, _, product, _, _ = options[chosen]
        if event_time >= burn:
            production += product
        time = event_time
        events += 1
    return {"tof": production / (size * duration), "events": events, "seed": seed}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checks", action="store_true")
    parser.add_argument("--renewal", action="store_true")
    parser.add_argument("--checks-only", action="store_true")
    args = parser.parse_args()
    output = ROOT / "results" / "audit"
    output.mkdir(parents=True, exist_ok=True)
    rates = build_rates(745, 0.001, 0.05, 0.0)
    if args.renewal:
        rows = []
        for oxidation in (0.16, 39120.0):
            for barrier in (0.95, 1.45):
                rates = build_rates(745, 0.001, 0.05, 0.0, params={"Ea_sulfate": barrier})
                rates[("site", P_OXIDISE_CU)] = oxidation
                for renewal in np.logspace(-2, 4, 13):
                    for name in ("single", "pair", "path3", "triangle", "path4", "star4", "square"):
                        size, edges = MOTIFS[name]
                        result = exact(size, edges, rates, renewal=float(renewal))
                        rows.append(dict(motif=name, oxidation=oxidation, barrier=barrier,
                                         renewal=float(renewal), copper_yield=result["tof"] / renewal,
                                         **result))
                print(f"Renewal sweep completed: oxidation={oxidation}, barrier={barrier}", flush=True)
        with (output / "renewal.csv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
        return
    rows = []
    for temperature in (() if args.checks_only else (673, 745, 823)):
        for diffusion in (0.0, 1e2, 1e4, 1e6):
            for name, (size, edges) in MOTIFS.items():
                result = exact(size, edges, build_rates(temperature, 0.001, 0.05, 0.0), diffusion)
                rows.append(dict(motif=name, temperature=temperature, diffusion=diffusion, **result))
        print(f"Exact motif sweep completed at {temperature} K", flush=True)
    if rows:
        with (output / "motifs.csv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    if args.checks or args.checks_only:
        checks = []
        for name in ("pair", "path3", "triangle", "square"):
            size, edges = MOTIFS[name]
            expected = exact(size, edges, rates)["tof"]
            duration = 20.0
            samples = [simulate(size, edges, rates, seed, duration, duration) for seed in range(12)]
            observed = np.array([sample["tof"] for sample in samples])
            sem = observed.std(ddof=1) / np.sqrt(len(observed))
            passed = abs(observed.mean() - expected) <= 5 * sem + 0.03 * expected + 1e-8
            check = dict(motif=name, exact=expected, mean=float(observed.mean()), sem=float(sem),
                         passed=bool(passed), duration=duration, samples=samples)
            checks.append(check)
            print(json.dumps({key: value for key, value in check.items() if key != "samples"}), flush=True)
        (output / "checks.json").write_text(json.dumps(checks, indent=2))
        if not all(check["passed"] for check in checks):
            raise AssertionError("Independent SSA disagrees with stationary solution")
    for row in rows:
        if row["temperature"] == 745 and row["diffusion"] == 0:
            print(row["motif"], row["tof"], flush=True)


if __name__ == "__main__":
    main()
