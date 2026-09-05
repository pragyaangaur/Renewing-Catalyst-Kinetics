"""
Rate constant construction for the copper oxide surface model.

Provenance and honesty note
---------------------------
Relevant CuO/SO2 DFT studies exist, including Liu et al. (2017),
doi:10.1016/j.apsusc.2017.01.088, and Barua and Padak (2025),
doi:10.1021/acs.energyfuels.5c03924. No barrier table from those studies has been
used here. The values below are illustrative assumptions, not measured or
verified literature-derived parameters. See docs/LIMITATIONS.md.

No experimental data is used to fit these values. A single conversion measurement
could not identify them in any case, especially without a measured feed composition
and site density. Where a scan needs a concrete conversion target it uses a round
illustrative figure, stated as such at the point of use.

Adsorption rate constants are built from kinetic theory impingement flux, so they
are multiplied by assumed sticking factors. Desorption and surface reaction constants use a
transition state prefactor of 1e13 per second.
"""

import math
from kmc_core import (
    P_OXIDISE_CU, P_ADS_SO2, P_DES_SO2, P_DES_SO3, P_SULFATE,
    P_DESULFATE, P_SPALL_SULF, P_SPALL_OX, P_ADS_SO3,
    B_ADS_O2, B_DES_O2, B_LH,
)

KB_EV = 8.617333262e-5      # eV/K
KB_J = 1.380649e-23         # J/K
NA = 6.02214076e23
NU = 1.0e13                 # s^-1, transition state prefactor
A_SITE = 1.02e-19           # m^2, area per surface site (~3.2 angstrom square)
P_ATM = 101325.0            # Pa

MASS = {"SO2": 64.066e-3 / NA, "SO3": 80.066e-3 / NA, "O2": 31.998e-3 / NA}

# ------------------------------------------------------------------ base barriers

BASE = dict(
    Ea_des_SO2    = 1.00,   # eV, molecular SO2 desorption from oxide
    Ea_ads_O2     = 0.30,   # eV, activated dissociative adsorption
    S0_O2         = 0.05,
    S0_SO2        = 0.30,
    S0_SO3        = 0.30,
    Ea_des_O2     = 2.00,   # eV, recombinative oxygen desorption
    Ea_LH         = 1.05,   # eV, SO2* + O* -> SO3*, rate determining step
    Ea_des_SO3    = 1.15,   # eV, SO3 desorption, the productive branch
    Ea_sulfate    = 0.95,   # eV, SO3* + CuO -> CuSO4, the poisoning branch  [SCANNED]
    Ea_desulfate  = 1.80,   # eV, CuSO4 -> CuO + SO3(g)                      [SCANNED]
    # Reoxidation of freshly exposed metallic copper. The lattice model resolves a
    # single surface layer, and first monolayer oxidation of clean copper in an
    # oxygen bearing gas at these temperatures is fast and close to non activated.
    # Thick scale regrowth is diffusion limited and much slower; that regime is
    # reached by scanning S0_ox_Cu downward and is treated explicitly in scan S5.
    Ea_ox_Cu      = 0.20,   # eV, metallic Cu -> CuO in situ
    S0_ox_Cu      = 0.10,   #                                                 [SCANNED]
    k_spall_sulf  = 1.0e-3, # s^-1, sulphate layer spallation exposing fresh Cu [SCANNED]
    k_spall_ox    = 1.0e-6, # s^-1, oxide erosion exposing fresh Cu
)


def impingement(p_pa, T, species):
    """Molecular impingement rate per surface site, s^-1."""
    m = MASS[species]
    return p_pa * A_SITE / math.sqrt(2.0 * math.pi * m * KB_J * T)


def arr(Ea_eV, T, nu=NU):
    return nu * math.exp(-Ea_eV / (KB_EV * T))


def build_rates(T, y_SO2, y_O2, y_SO3, p_total=P_ATM, params=None):
    """
    Assemble the full rate constant dictionary.

    T       temperature, K
    y_SO2   SO2 mole fraction in the gas contacting the surface
    y_O2    O2 mole fraction
    y_SO3   SO3 mole fraction (sets the readsorption / back pressure term)
    """
    q = dict(BASE)
    if params:
        q.update(params)

    p_SO2 = y_SO2 * p_total
    p_O2 = y_O2 * p_total
    p_SO3 = y_SO3 * p_total

    k = {}
    k[("site", P_ADS_SO2)]    = q["S0_SO2"] * impingement(p_SO2, T, "SO2")
    k[("site", P_ADS_SO3)]    = q["S0_SO3"] * impingement(p_SO3, T, "SO3")
    k[("site", P_DES_SO2)]    = arr(q["Ea_des_SO2"], T)
    k[("site", P_DES_SO3)]    = arr(q["Ea_des_SO3"], T)
    k[("site", P_SULFATE)]    = arr(q["Ea_sulfate"], T)
    k[("site", P_DESULFATE)]  = arr(q["Ea_desulfate"], T)
    k[("site", P_SPALL_SULF)] = q["k_spall_sulf"]
    k[("site", P_SPALL_OX)]   = q["k_spall_ox"]
    k[("site", P_OXIDISE_CU)] = q["S0_ox_Cu"] * impingement(p_O2, T, "O2") \
                                * math.exp(-q["Ea_ox_Cu"] / (KB_EV * T))

    k[("bond", B_ADS_O2)] = q["S0_O2"] * impingement(p_O2, T, "O2") \
                            * math.exp(-q["Ea_ads_O2"] / (KB_EV * T))
    k[("bond", B_DES_O2)] = arr(q["Ea_des_O2"], T)
    k[("bond", B_LH)]     = arr(q["Ea_LH"], T)
    return k


if __name__ == "__main__":
    T = 745.0
    k = build_rates(T, 1000e-6, 0.05, 500e-6)
    names = {
        ("site", P_OXIDISE_CU): "Cu -> CuO (in situ)",
        ("site", P_ADS_SO2): "SO2 adsorption",
        ("site", P_DES_SO2): "SO2 desorption",
        ("site", P_DES_SO3): "SO3 desorption  [PRODUCT]",
        ("site", P_ADS_SO3): "SO3 readsorption",
        ("site", P_SULFATE): "sulfation       [POISON]",
        ("site", P_DESULFATE): "desulfation",
        ("site", P_SPALL_SULF): "sulphate spallation",
        ("site", P_SPALL_OX): "oxide spallation",
        ("bond", B_ADS_O2): "O2 dissociative adsorption",
        ("bond", B_DES_O2): "O2 recombinative desorption",
        ("bond", B_LH): "SO2* + O* -> SO3*  [RDS]",
    }
    print(f"Rate constants at T = {T:.0f} K, 1000 ppm SO2, 5% O2, 500 ppm SO3\n")
    for key, nm in names.items():
        print(f"  {nm:32s} {k[key]:12.4g} s^-1")
    print()
    br = k[("site", P_SULFATE)] / k[("site", P_DES_SO3)]
    print(f"  Branching ratio sulfation / SO3 desorption = {br:.1f}")
    print("  Values above 1 mean an adsorbed SO3 is more likely to poison its")
    print("  own site than to leave as product.")
