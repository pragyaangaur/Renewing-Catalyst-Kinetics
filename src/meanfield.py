"""
Mean field (site approximation) companion to the lattice KMC.

Built from exactly the same elementary process list as kmc_core, with pair
processes closed at the site level, so that a bond between neighbours is free
free with probability theta_free squared and so on. Coordination number z = 4
matches the square lattice used by the KMC.

Purpose. The KMC is exact but costs roughly 6e4 events per second in pure
Python, which is too slow for the dense two dimensional parameter maps this
study needs. The mean field model runs a full steady state in milliseconds. The
workflow is the standard multiscale one: validate mean field against KMC on a
spread of points, map parameter space with mean field, then return to KMC where
the two disagree, because disagreement is the signature of spatial correlation
and is itself a result.

State vector
    f_Cu     fraction of sites that are bare metallic copper
    f_SO4    fraction of sites that are copper sulphate (dead)
    th_O     atomic oxygen coverage, referred to total sites
    th_SO2   adsorbed SO2 coverage
    th_SO3   adsorbed SO3 coverage
with f_CuO = 1 - f_Cu - f_SO4 and th_free = f_CuO - th_O - th_SO2 - th_SO3.
"""

import numpy as np
from scipy.integrate import solve_ivp

from kmc_core import (
    P_OXIDISE_CU, P_ADS_SO2, P_DES_SO2, P_DES_SO3, P_SULFATE,
    P_DESULFATE, P_SPALL_SULF, P_SPALL_OX, P_ADS_SO3,
    B_ADS_O2, B_DES_O2, B_LH,
)

Z = 4.0  # square lattice coordination


def unpack(k):
    return (
        k[("site", P_OXIDISE_CU)], k[("site", P_ADS_SO2)], k[("site", P_DES_SO2)],
        k[("site", P_DES_SO3)], k[("site", P_ADS_SO3)], k[("site", P_SULFATE)],
        k[("site", P_DESULFATE)], k[("site", P_SPALL_SULF)], k[("site", P_SPALL_OX)],
        k[("bond", B_ADS_O2)], k[("bond", B_DES_O2)], k[("bond", B_LH)],
    )


def rhs(t, y, k):
    f_Cu, f_SO4, th_O, th_SO2, th_SO3 = y
    (k_ox, k_aSO2, k_dSO2, k_dSO3, k_aSO3, k_sulf,
     k_desulf, k_spS, k_spO, k_aO2, k_dO2, k_lh) = k

    f_Cu = max(f_Cu, 0.0)
    f_SO4 = max(f_SO4, 0.0)
    f_CuO = max(1.0 - f_Cu - f_SO4, 0.0)
    th_O = max(th_O, 0.0)
    th_SO2 = max(th_SO2, 0.0)
    th_SO3 = max(th_SO3, 0.0)
    th_free = max(f_CuO - th_O - th_SO2 - th_SO3, 0.0)

    # elementary fluxes, all per total site per second
    r_ox      = k_ox * f_Cu                      # Cu -> CuO
    r_aSO2    = k_aSO2 * th_free
    r_dSO2    = k_dSO2 * th_SO2
    r_aSO3    = k_aSO3 * th_free
    r_dSO3    = k_dSO3 * th_SO3
    r_sulf    = k_sulf * th_SO3                  # CuO + SO3* -> CuSO4
    r_desulf  = k_desulf * f_SO4
    r_spS     = k_spS * f_SO4                    # sulphate spallation -> fresh Cu
    r_spO     = k_spO * f_CuO                    # oxide erosion -> fresh Cu
    r_aO2     = Z * k_aO2 * th_free * th_free    # produces 2 O per event, folded in
    r_dO2     = Z * k_dO2 * th_O * th_O
    r_lh      = Z * k_lh * th_SO2 * th_O         # SO2* + O* -> SO3* + free

    # adsorbates are lost when their host site changes phase
    frac = (th_O / f_CuO, th_SO2 / f_CuO, th_SO3 / f_CuO) if f_CuO > 1e-12 else (0.0, 0.0, 0.0)

    df_Cu = -r_ox + r_spS + r_spO
    df_SO4 = r_sulf - r_desulf - r_spS
    dth_O = r_aO2 - r_dO2 - r_lh - r_spO * frac[0]
    dth_SO2 = r_aSO2 - r_dSO2 - r_lh - r_spO * frac[1]
    dth_SO3 = r_lh + r_aSO3 - r_dSO3 - r_sulf - r_spO * frac[2]

    return [df_Cu, df_SO4, dth_O, dth_SO2, dth_SO3]


def steady_state(k_dict, t_end=1.0e5, y0=None, rtol=1e-8, atol=1e-12):
    """
    Integrate to steady state. Returns a dict of state and observable rates.

    t_end is generous because sulphate spallation can be as slow as 1e-6 per
    second and the slow phase variables must be allowed to converge.
    """
    k = unpack(k_dict)
    if y0 is None:
        y0 = [0.0, 0.0, 0.0, 0.0, 0.0]  # start fully CuO, clean
    sol = solve_ivp(rhs, (0.0, t_end), y0, args=(k,), method="LSODA",
                    rtol=rtol, atol=atol, dense_output=False)
    y = sol.y[:, -1]
    f_Cu, f_SO4, th_O, th_SO2, th_SO3 = [max(v, 0.0) for v in y]
    f_CuO = max(1.0 - f_Cu - f_SO4, 0.0)
    th_free = max(f_CuO - th_O - th_SO2 - th_SO3, 0.0)

    (k_ox, k_aSO2, k_dSO2, k_dSO3, k_aSO3, k_sulf,
     k_desulf, k_spS, k_spO, k_aO2, k_dO2, k_lh) = k

    tof_gross = k_dSO3 * th_SO3                 # SO3 desorption events per site per s
    tof_net = tof_gross - k_aSO3 * th_free + k_desulf * f_SO4
    r_lh = Z * k_lh * th_SO2 * th_O

    return {
        "ok": sol.success,
        "f_Cu": f_Cu, "f_CuO": f_CuO, "f_CuSO4": f_SO4,
        "th_O": th_O, "th_SO2": th_SO2, "th_SO3": th_SO3, "th_free": th_free,
        "tof_gross": tof_gross,       # per site per second
        "tof_net": tof_net,
        "r_lh": r_lh,
        "r_sulfation": k_sulf * th_SO3,
        "r_desulfation": k_desulf * f_SO4,
        "r_spall": k_spS * f_SO4 + k_spO * f_CuO,
        "cu_consumption": k_ox * f_Cu,   # sites of metallic Cu oxidised per site per s
    }


def trajectory(k_dict, t_end, n_points=400, y0=None):
    k = unpack(k_dict)
    if y0 is None:
        y0 = [0.0, 0.0, 0.0, 0.0, 0.0]
    ts = np.geomspace(max(t_end * 1e-12, 1e-12), t_end, n_points)
    sol = solve_ivp(rhs, (0.0, t_end), y0, args=(k,), method="LSODA",
                    t_eval=ts, rtol=1e-8, atol=1e-12)
    return sol.t, sol.y
