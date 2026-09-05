"""
Coupling the surface model to a one dimensional plug flow reactor.

This connects site level turnover to a reactor duty. The conversion target used
below is an illustrative sizing scenario and not a calibration against measured
data. See docs/CORRECTIONS.md, which withdraws the turnover requirement that an
earlier version of this work derived here.

Marching scheme. The reactor is divided into axial slices. In each slice the gas
composition is treated as locally uniform, the surface is relaxed to its steady
state under that local composition, and the resulting turnover frequency is used
to advance the gas composition into the next slice. Because SO3 accumulates along
the bed, the sulphation driving force rises with axial position, so the surface
state is not uniform. The model resolves that gradient rather than assuming it
away.

Site density. The substrate is a gas permeable copper structure. Its geometric
specific surface area per unit reactor volume is the bridge between turnover
frequency and volumetric rate.
"""

import numpy as np
from params import build_rates, A_SITE
from meanfield import steady_state

R_GAS = 8.314462618
P_ATM = 101325.0

# Specific geometric surface area of the substrate, m^2 per m^3 of reactor volume.
# A honeycomb or wire mesh in this duty typically falls between 500 and 3000.
DEFAULT_AV = 1000.0
SITES_PER_M2 = 1.0 / A_SITE          # about 9.8e18

# Conversion target used when sizing. This is a round illustrative figure chosen
# to make the scenario concrete. It is not fitted to any measurement.
ILLUSTRATIVE_CONVERSION = 0.90


def site_density(a_v=DEFAULT_AV):
    """Active surface sites per cubic metre of reactor volume."""
    return a_v * SITES_PER_M2


def gas_number_density(T, p=P_ATM):
    """Molecules per cubic metre."""
    return p / (1.380649e-23 * T)


def required_tof(conversion, y_SO2_in, tau, T, a_v=DEFAULT_AV, p=P_ATM):
    """
    Turnover frequency, per site per second, that a plug flow reactor needs on
    average to reach the stated conversion in the stated residence time.
    """
    n_gas = gas_number_density(T, p)
    c_in = y_SO2_in * n_gas                      # SO2 molecules per m^3
    rate_needed = conversion * c_in / tau        # molecules per m^3 per s
    return rate_needed / site_density(a_v)


def march(T, y_SO2_in, y_O2, tau, n_slices=60, a_v=DEFAULT_AV,
          params=None, y_SO3_in=0.0, p=P_ATM, verbose=False):
    """
    March the plug flow reactor and return axial profiles.

    Returns a dict of arrays indexed by slice, plus outlet conversion.
    """
    n_gas = gas_number_density(T, p)
    S_v = site_density(a_v)
    dt = tau / n_slices

    y_SO2 = y_SO2_in
    y_SO3 = y_SO3_in

    z, ySO2, ySO3, tof, fSO4, fCuO, fCu = [], [], [], [], [], [], []

    for i in range(n_slices):
        k = build_rates(T, max(y_SO2, 1e-12), y_O2, max(y_SO3, 1e-12), p_total=p,
                        params=params)
        ss = steady_state(k)

        rate = ss["tof_net"] * S_v          # molecules per m^3 per s produced as SO3
        d_c = rate * dt                     # molecules per m^3 converted in this slice
        d_y = d_c / n_gas

        # cannot convert more SO2 than is present
        d_y = min(d_y, y_SO2)

        z.append((i + 0.5) * tau / n_slices)
        ySO2.append(y_SO2)
        ySO3.append(y_SO3)
        tof.append(ss["tof_net"])
        fSO4.append(ss["f_CuSO4"])
        fCuO.append(ss["f_CuO"])
        fCu.append(ss["f_Cu"])

        y_SO2 -= d_y
        y_SO3 += d_y
        y_SO2 = max(y_SO2, 0.0)

    conv = (y_SO2_in - y_SO2) / y_SO2_in if y_SO2_in > 0 else 0.0
    return {
        "z": np.array(z), "y_SO2": np.array(ySO2), "y_SO3": np.array(ySO3),
        "tof": np.array(tof), "f_CuSO4": np.array(fSO4),
        "f_CuO": np.array(fCuO), "f_Cu": np.array(fCu),
        "conversion": conv, "y_SO2_out": y_SO2, "y_SO3_out": y_SO3,
    }


if __name__ == "__main__":
    T = 745.0
    y_in = 1000e-6
    tau = 1.0
    print("Illustrative reactor sizing scenario, not a calibration")
    print(f"  T = {T:.0f} K, inlet SO2 = {y_in*1e6:.0f} ppm, residence time = {tau:.1f} s")
    print(f"  gas number density        = {gas_number_density(T):.3e} molecules/m3")
    print()
    print(f"  {'a_v (m2/m3)':>12} {'sites/m3':>12} {'required TOF (1/s)':>20}")
    for a_v in [200, 500, 1000, 2000, 3000]:
        print(f"  {a_v:12.0f} {site_density(a_v):12.3e} "
              f"{required_tof(ILLUSTRATIVE_CONVERSION, y_in, tau, T, a_v):20.4f}")
    print()
    print("  The required turnover frequency is of order one per site per second.")
    print("  Compare that with the turnover frequencies the surface model produces")
    print("  in its live regime, which are of order 1e3 to 1e4 per site per second.")
