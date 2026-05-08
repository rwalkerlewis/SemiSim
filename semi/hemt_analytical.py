"""
Closed-form HEMT analytical references (M17, ADR 0016).

Pure-Python helpers for the `benchmarks/hemt_2d/` 2DEG sheet-density
verifier. The module is import-clean (no dolfinx) so it can be exercised
by `tests/test_hemt_analytical.py` on the pure-Python CI matrix.

The classical-electrostatic relation between the 2DEG sheet density
`n_s` and the gate bias `V_GS` in the linear regime (above the HEMT
threshold) follows Pozela & Reklaitis 1980 and Sze 3rd ed Section 7.4:

    n_s(V_GS) = (eps_AlGaAs * eps_0 / (q * d_barrier))
                * (V_GS - V_T_HEMT)

where `d_barrier` is the AlGaAs barrier thickness (m), `eps_AlGaAs` is
the static relative permittivity of AlGaAs at the chosen Al fraction,
and `V_T_HEMT` is the HEMT threshold voltage. The threshold is the
Schottky-gate-induced flat-band shift minus the conduction-band
discontinuity:

    V_T_HEMT = phi_B - delta_E_C - (q * N_D_barrier * d_barrier^2)
                                   / (2 * eps_AlGaAs * eps_0)

with `phi_B` the metal-semiconductor barrier height (eV), `delta_E_C =
chi_GaAs - chi_AlGaAs` the Anderson-rule conduction-band offset (eV),
and `N_D_barrier` the modulation-doping density (m^-3).

ADR 0016 documents that this classical reference is the V&V gate for
the discontinuous-coefficient case in M17; quantum confinement
corrections (which would close the ~15 % gap to a fully-quantum
Poisson-Schrodinger reference) are a future milestone.
"""
from __future__ import annotations

from .constants import EPS0, Q


def hemt_threshold_voltage(
    *,
    barrier_height_eV: float,
    delta_Ec_eV: float,
    N_D_barrier_m3: float,
    d_barrier_m: float,
    eps_r_barrier: float,
) -> float:
    """Classical HEMT threshold voltage `V_T_HEMT` in volts.

    The threshold is the gate bias at which the 2DEG charge sheet
    just begins to populate. Below threshold the channel is
    depleted; above threshold the linear-regime formula in
    `hemt_2deg_classical` applies.

    Parameters
    ----------
    barrier_height_eV
        Metal-semiconductor (Schottky) barrier height phi_B in eV.
    delta_Ec_eV
        Anderson-rule conduction-band offset chi_GaAs - chi_AlGaAs
        in eV.
    N_D_barrier_m3
        Donor concentration in the AlGaAs barrier (m^-3).
    d_barrier_m
        AlGaAs barrier thickness (m).
    eps_r_barrier
        Relative permittivity of AlGaAs.
    """
    if d_barrier_m <= 0.0:
        raise ValueError("d_barrier_m must be positive")
    if eps_r_barrier <= 0.0:
        raise ValueError("eps_r_barrier must be positive")
    space_charge_term = (
        Q * N_D_barrier_m3 * d_barrier_m ** 2
        / (2.0 * eps_r_barrier * EPS0)
    )
    return float(barrier_height_eV - delta_Ec_eV - space_charge_term)


def hemt_2deg_classical(
    V_GS: float,
    *,
    barrier_height_eV: float,
    delta_Ec_eV: float,
    N_D_barrier_m3: float,
    d_barrier_m: float,
    eps_r_barrier: float,
) -> float:
    """Classical 2DEG sheet density `n_s(V_GS)` in cm^-2.

    Returns zero below the HEMT threshold V_T_HEMT (the channel is
    depleted). Above threshold the linear-regime relation gives
    `n_s = (eps_AlGaAs eps_0 / (q d_barrier)) * (V_GS - V_T_HEMT)`,
    converted to cm^-2 (per-area-unit, the standard HEMT reporting
    convention).
    """
    V_T = hemt_threshold_voltage(
        barrier_height_eV=barrier_height_eV,
        delta_Ec_eV=delta_Ec_eV,
        N_D_barrier_m3=N_D_barrier_m3,
        d_barrier_m=d_barrier_m,
        eps_r_barrier=eps_r_barrier,
    )
    delta_V = float(V_GS) - V_T
    if delta_V <= 0.0:
        return 0.0
    n_s_per_m2 = (eps_r_barrier * EPS0 / (Q * d_barrier_m)) * delta_V
    # Convert m^-2 to cm^-2 for the device-physics reporting unit
    # (m^-2 -> cm^-2 = m^-2 * 1e-4).
    return float(n_s_per_m2 * 1.0e-4)


def hemt_2deg_classical_si_units(
    V_GS: float,
    *,
    barrier_height_eV: float,
    delta_Ec_eV: float,
    N_D_barrier_m3: float,
    d_barrier_m: float,
    eps_r_barrier: float,
) -> float:
    """Same as `hemt_2deg_classical` but returns m^-2 (SI) directly."""
    return hemt_2deg_classical(
        V_GS,
        barrier_height_eV=barrier_height_eV,
        delta_Ec_eV=delta_Ec_eV,
        N_D_barrier_m3=N_D_barrier_m3,
        d_barrier_m=d_barrier_m,
        eps_r_barrier=eps_r_barrier,
    ) * 1.0e4
