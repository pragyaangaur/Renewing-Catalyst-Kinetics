"""
Thermodynamic baseline for the CuO/CuSO4 sulfation boundary.

This frames the kinetic Monte Carlo work against a thermodynamic boundary rather
than an assertion. See docs/CORRECTIONS.md: the sharply worded onset temperature
from the first pass was withdrawn, and this calculation should be read as an
order of magnitude guide rather than a firm operating limit.

Reaction studied (decomposition direction):

    CuSO4(s)  ->  CuO(s) + SO3(g)        Kp = p_SO3 / p_ref

If the actual SO3 partial pressure in the gas exceeds Kp, bulk copper sulphate is
the thermodynamically stable phase and the oxide is predicted to sulfate.

Thermochemical data are standard 298 K values from the NIST-JANAF / CODATA
compilations. Heat capacities are included via a Shomate-style constant-Cp
correction over the range of interest, which is adequate for a boundary estimate
at the few-kelvin level.
"""

import numpy as np

R = 8.314462618  # J/mol/K

# Standard formation enthalpy (J/mol) and absolute entropy (J/mol/K) at 298.15 K
SPECIES = {
    "CuSO4_s": {"Hf": -771_400.0, "S": 109.2, "Cp": 100.0},
    "CuO_s":   {"Hf": -157_300.0, "S":  42.6, "Cp":  47.0},
    "SO3_g":   {"Hf": -395_700.0, "S": 256.8, "Cp":  62.0},
}

T_REF = 298.15


def delta_props():
    """Reaction deltas for CuSO4(s) -> CuO(s) + SO3(g) at 298.15 K."""
    dH = SPECIES["CuO_s"]["Hf"] + SPECIES["SO3_g"]["Hf"] - SPECIES["CuSO4_s"]["Hf"]
    dS = SPECIES["CuO_s"]["S"] + SPECIES["SO3_g"]["S"] - SPECIES["CuSO4_s"]["S"]
    dCp = SPECIES["CuO_s"]["Cp"] + SPECIES["SO3_g"]["Cp"] - SPECIES["CuSO4_s"]["Cp"]
    return dH, dS, dCp


def delta_G(T):
    """Gibbs energy of decomposition at temperature T (K), with constant-Cp correction."""
    dH0, dS0, dCp = delta_props()
    T = np.asarray(T, dtype=float)
    dH = dH0 + dCp * (T - T_REF)
    dS = dS0 + dCp * np.log(T / T_REF)
    return dH - T * dS


def Kp_atm(T):
    """Equilibrium SO3 partial pressure over the CuSO4/CuO couple, in atm."""
    return np.exp(-delta_G(T) / (R * np.asarray(T, dtype=float)))


def equilibrium_so3_ppm(T):
    """Same quantity expressed as ppmv at 1 atm total pressure."""
    return Kp_atm(T) * 1e6


def sulfation_onset_temperature(p_so3_ppm, lo=500.0, hi=1400.0):
    """
    Temperature above which CuO is the stable phase at the given SO3 level.
    Below this temperature bulk CuSO4 is thermodynamically favoured.
    """
    target = p_so3_ppm / 1e6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if Kp_atm(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def surface_corrected_onset(p_so3_ppm, destabilisation_kJ):
    """
    Repeat the onset calculation for a *surface* sulphate rather than bulk
    crystalline CuSO4.

    A chemisorbed sulphate on an oxide surface does not gain the full lattice
    energy of the three dimensional CuSO4 crystal, so it is less stable by some
    destabilisation energy. This shifts the onset temperature down. The magnitude
    is the single most important unknown in the sulphation argument.
    """
    dH0, dS0, dCp = delta_props()
    dH0_surf = dH0 - destabilisation_kJ * 1000.0  # easier to decompose

    def kp(T):
        dH = dH0_surf + dCp * (T - T_REF)
        dS = dS0 + dCp * np.log(T / T_REF)
        return np.exp(-(dH - T * dS) / (R * T))

    target = p_so3_ppm / 1e6
    lo, hi = 300.0, 1400.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if kp(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    dH0, dS0, dCp = delta_props()
    print("CuSO4(s) -> CuO(s) + SO3(g)")
    print(f"  dH(298) = {dH0/1000:8.1f} kJ/mol")
    print(f"  dS(298) = {dS0:8.1f} J/mol/K")
    print(f"  dCp     = {dCp:8.1f} J/mol/K")
    print()

    print("Equilibrium SO3 over the CuSO4/CuO couple (bulk phases):")
    print(f"  {'T (K)':>8} {'p_SO3 eq (ppm)':>18}")
    for T in [673, 700, 745, 760, 790, 823, 850, 900, 950]:
        print(f"  {T:8d} {equilibrium_so3_ppm(T):18.3g}")
    print()

    # Illustrative duty: 1000 ppm SO2 inlet at a range of conversions.
    for so2_in, conv in [(500, 0.90), (1000, 0.90), (2000, 0.90), (1000, 0.5)]:
        so3 = so2_in * conv
        T_on = sulfation_onset_temperature(so3)
        print(f"  {so2_in:5d} ppm SO2 in, {conv*100:4.1f}% conv -> {so3:6.0f} ppm SO3 "
              f"| bulk sulfation onset {T_on:6.1f} K")
    print()

    print("Surface sulphate destabilisation sensitivity (900 ppm SO3):")
    print(f"  {'destab (kJ/mol)':>18} {'onset T (K)':>14}")
    for d in [0, 10, 20, 30, 40, 50, 60, 70, 80]:
        print(f"  {d:18d} {surface_corrected_onset(900.0, d):14.1f}")
