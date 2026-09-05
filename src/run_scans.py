"""
Parameter scans for the copper oxide surface and reactor model.

Each scan is written to results/ as a .npz archive so that figures and the report
can be regenerated without rerunning the sweeps.

S1  temperature sweep, surface state and reactor conversion
S2  activity headroom: conversion against sulphated fraction
S3  two dimensional map, temperature against sulphate stability
S4  two dimensional map, temperature against sulphation barrier
S5  self renewal: spallation rate against copper reoxidation rate
S6  SO2 inlet concentration sweep, testing the low partial pressure argument
S7  axial profiles through the reactor
"""

import os
import json
import numpy as np

from params import build_rates, BASE
from meanfield import steady_state
from reactor import march, required_tof, DEFAULT_AV, ILLUSTRATIVE_CONVERSION

RESULTS = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(RESULTS, exist_ok=True)

# Nominal illustrative duty used throughout unless a scan varies it.
DUTY = dict(y_SO2_in=1000e-6, y_O2=0.05, tau=1.0)
TARGET_CONV = ILLUSTRATIVE_CONVERSION
T_TEST = 745.0

# A parameter set in the live catalytic regime, used as the reference case.
# Ea_sulfate and Ea_desulfate here are not measured; they are the values the
# inversion in S2 and S3 shows are required to reproduce the measured datum.
LIVE = dict(Ea_sulfate=1.45, Ea_desulfate=1.55)


def save(name, **arrays):
    path = os.path.join(RESULTS, name + ".npz")
    np.savez_compressed(path, **arrays)
    print(f"  wrote {os.path.relpath(path)}")


# --------------------------------------------------------------------------- S1

def scan_temperature():
    print("S1 temperature sweep")
    Ts = np.linspace(620.0, 900.0, 57)
    out = {k: [] for k in ["T", "tof", "f_CuSO4", "f_CuO", "f_Cu", "conv", "req_tof"]}
    for T in Ts:
        k = build_rates(T, DUTY["y_SO2_in"], DUTY["y_O2"], 1e-9, params=LIVE)
        ss = steady_state(k)
        r = march(T, DUTY["y_SO2_in"], DUTY["y_O2"], DUTY["tau"], params=LIVE)
        out["T"].append(T)
        out["tof"].append(ss["tof_net"])
        out["f_CuSO4"].append(ss["f_CuSO4"])
        out["f_CuO"].append(ss["f_CuO"])
        out["f_Cu"].append(ss["f_Cu"])
        out["conv"].append(r["conversion"])
        out["req_tof"].append(required_tof(TARGET_CONV, DUTY["y_SO2_in"], DUTY["tau"], T))
    save("s1_temperature", **{k: np.array(v) for k, v in out.items()})
    return out


# --------------------------------------------------------------------------- S2

def scan_headroom():
    """
    How much of the surface may be sulphated before the reactor stops meeting the
    measured conversion?

    Sulphate stability is walked from labile to refractory. That drives the steady
    state sulphated fraction across its whole range, from a nearly clean oxide to a
    nearly fully sulphated one. For each point the surface turnover frequency is
    compared with the turnover frequency the reactor actually needs, which the
    sizing calculation puts at about 0.89 per site per second.

    The crossing point is the central quantitative result of the study.
    """
    print("S2 activity headroom")
    eds = np.linspace(1.05, 2.60, 260)
    fs, tofs, convs = [], [], []
    for ed in eds:
        p = dict(LIVE, Ea_desulfate=float(ed))
        k = build_rates(T_TEST, DUTY["y_SO2_in"], DUTY["y_O2"], 1e-9, params=p)
        ss = steady_state(k)
        fs.append(ss["f_CuSO4"])
        tofs.append(max(ss["tof_net"], 0.0))
        r = march(T_TEST, DUTY["y_SO2_in"], DUTY["y_O2"], DUTY["tau"],
                  n_slices=25, params=p)
        convs.append(r["conversion"])
    fs, tofs, convs = np.array(fs), np.array(tofs), np.array(convs)

    req = required_tof(TARGET_CONV, DUTY["y_SO2_in"], DUTY["tau"], T_TEST)

    # crossing: largest sulphated fraction whose turnover still clears the requirement
    ok = tofs >= req
    f_crit = float(fs[ok].max()) if ok.any() else float("nan")
    active_crit = 1.0 - f_crit
    tof_clean = float(tofs.min()), float(tofs.max())

    save("s2_headroom", Ea_desulfate=eds, f_CuSO4=fs, tof=tofs, conv=convs,
         req_tof=np.array([req]), f_crit=np.array([f_crit]))
    print(f"  turnover frequency required by the reactor : {req:.4f} 1/s")
    print(f"  turnover frequency of a clean oxide surface: {tof_clean[1]:.4g} 1/s")
    print(f"  headroom factor                            : {tof_clean[1]/req:.4g}")
    print(f"  critical sulphated fraction                : {f_crit:.6f}")
    print(f"  surviving active fraction at that point    : {active_crit:.3e}")
    return eds, fs, tofs, f_crit


# --------------------------------------------------------------------------- S3

def map_T_vs_desulfation():
    print("S3 map, temperature against sulphate stability")
    Ts = np.linspace(620.0, 900.0, 40)
    eds = np.linspace(1.20, 2.20, 40)
    conv = np.zeros((len(eds), len(Ts)))
    fso4 = np.zeros_like(conv)
    for i, ed in enumerate(eds):
        for j, T in enumerate(Ts):
            p = dict(LIVE, Ea_desulfate=float(ed))
            r = march(T, DUTY["y_SO2_in"], DUTY["y_O2"], DUTY["tau"],
                      n_slices=25, params=p)
            conv[i, j] = r["conversion"]
            fso4[i, j] = r["f_CuSO4"].mean()
    save("s3_map_desulf", T=Ts, Ea_desulfate=eds, conv=conv, f_CuSO4=fso4)
    return Ts, eds, conv


# --------------------------------------------------------------------------- S4

def map_T_vs_sulfation():
    print("S4 map, temperature against sulphation barrier")
    Ts = np.linspace(620.0, 900.0, 40)
    eas = np.linspace(0.85, 1.70, 40)
    conv = np.zeros((len(eas), len(Ts)))
    fso4 = np.zeros_like(conv)
    for i, ea in enumerate(eas):
        for j, T in enumerate(Ts):
            p = dict(LIVE, Ea_sulfate=float(ea))
            r = march(T, DUTY["y_SO2_in"], DUTY["y_O2"], DUTY["tau"],
                      n_slices=25, params=p)
            conv[i, j] = r["conversion"]
            fso4[i, j] = r["f_CuSO4"].mean()
    save("s4_map_sulf", T=Ts, Ea_sulfate=eas, conv=conv, f_CuSO4=fso4)
    return Ts, eas, conv


# --------------------------------------------------------------------------- S5

def scan_renewal():
    """
    The self renewal claim. Spallation of sulphated material exposes fresh
    metallic copper, which reoxidises and returns to service. This scans the
    spallation rate against the rate at which fresh copper reoxidises.

    A slow reoxidation sticking coefficient represents diffusion limited regrowth
    of a thick oxide scale, a fast one represents first monolayer oxidation.
    """
    print("S5 self renewal, spallation against reoxidation")
    ksp = np.geomspace(1e-6, 1e6, 49)
    s0s = [1e-6, 1e-4, 1e-2, 1e-1]
    conv = np.zeros((len(s0s), len(ksp)))
    fso4 = np.zeros_like(conv)
    fcu = np.zeros_like(conv)
    # deliberately start from a poisoning parameter set, so that renewal has
    # something to rescue
    poison = dict(Ea_sulfate=1.15, Ea_desulfate=1.75)
    for i, s0 in enumerate(s0s):
        for j, ks in enumerate(ksp):
            p = dict(poison, k_spall_sulf=float(ks), S0_ox_Cu=float(s0))
            r = march(T_TEST, DUTY["y_SO2_in"], DUTY["y_O2"], DUTY["tau"],
                      n_slices=25, params=p)
            conv[i, j] = r["conversion"]
            fso4[i, j] = r["f_CuSO4"].mean()
            fcu[i, j] = r["f_Cu"].mean()
    save("s5_renewal", k_spall=ksp, S0_ox=np.array(s0s), conv=conv,
         f_CuSO4=fso4, f_Cu=fcu)
    for i, s0 in enumerate(s0s):
        jb = int(np.argmax(conv[i]))
        print(f"  S0_ox = {s0:.0e}: best conversion {conv[i, jb]*100:6.2f}% "
              f"at k_spall = {ksp[jb]:.2e} 1/s")
    return ksp, s0s, conv


# --------------------------------------------------------------------------- S6

def scan_so2_concentration():
    """
    The specification argues that the low SO3 partial pressures characteristic of
    dilute flue gas are what keep the copper oxide from sulphating. That is a
    testable statement: raising the inlet SO2, and therefore the SO3 produced,
    should push the surface toward sulphate and cost conversion.
    """
    print("S6 inlet SO2 sweep")
    ys = np.geomspace(50e-6, 50000e-6, 45)
    conv, fso4, tof = [], [], []
    for y in ys:
        r = march(T_TEST, float(y), DUTY["y_O2"], DUTY["tau"], params=LIVE)
        conv.append(r["conversion"])
        fso4.append(r["f_CuSO4"].mean())
        tof.append(r["tof"].mean())
    save("s6_so2", y_SO2_in=ys, conv=np.array(conv), f_CuSO4=np.array(fso4),
         tof=np.array(tof))
    return ys, np.array(conv), np.array(fso4)


# --------------------------------------------------------------------------- S7

def axial_profiles():
    print("S7 axial profiles")
    prof = {}
    for T in [700.0, 745.0, 800.0]:
        r = march(T, DUTY["y_SO2_in"], DUTY["y_O2"], DUTY["tau"],
                  n_slices=120, params=LIVE)
        prof[f"z_{int(T)}"] = r["z"]
        prof[f"ySO2_{int(T)}"] = r["y_SO2"]
        prof[f"ySO3_{int(T)}"] = r["y_SO3"]
        prof[f"fSO4_{int(T)}"] = r["f_CuSO4"]
        prof[f"tof_{int(T)}"] = r["tof"]
        print(f"  T = {T:.0f} K: conversion {r['conversion']*100:5.2f}%, "
              f"inlet sulphate {r['f_CuSO4'][0]:.4f}, outlet sulphate "
              f"{r['f_CuSO4'][-1]:.4f}")
    save("s7_axial", **prof)
    return prof


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] if len(sys.argv) > 1 else ["all"]
    run_all = "all" in which

    summary = {}
    if run_all or "s1" in which:
        scan_temperature()
    if run_all or "s2" in which:
        _, _, _, fcrit = scan_headroom()
        summary["f_crit"] = fcrit
    if run_all or "s3" in which:
        map_T_vs_desulfation()
    if run_all or "s4" in which:
        map_T_vs_sulfation()
    if run_all or "s5" in which:
        scan_renewal()
    if run_all or "s6" in which:
        scan_so2_concentration()
    if run_all or "s7" in which:
        axial_profiles()

    if summary:
        with open(os.path.join(RESULTS, "summary.json"), "w") as fh:
            json.dump(summary, fh, indent=2)
    print("done")
