"""
Pure-Python tests for the M17 HEMT classical-electrostatic analytical
reference (semi/hemt_analytical.py). These exercise the closed-form
2DEG sheet-density relation used by `benchmarks/hemt_2d/` for the M17
acceptance gate.
"""
from __future__ import annotations

import pytest

from semi.constants import EPS0, Q
from semi.hemt_analytical import (
    hemt_2deg_classical,
    hemt_2deg_classical_si_units,
    hemt_threshold_voltage,
)
from semi.materials import MATERIALS


def _algaas_gaas_params():
    """Standard AlGaAs_0p3 / GaAs HEMT parameters (M17, ADR 0016)."""
    al = MATERIALS["AlGaAs_0p3"]
    gaas = MATERIALS["GaAs"]
    return {
        "barrier_height_eV": 0.95,
        "delta_Ec_eV": float(gaas.chi - al.chi),  # ~0.33 eV
        "N_D_barrier_m3": 5.0e24,                 # 5e18 cm^-3
        "d_barrier_m": 30.0e-9,                   # 30 nm
        "eps_r_barrier": float(al.epsilon_r),     # 12.0
    }


def test_threshold_voltage_in_normally_on_range():
    """The benchmark parameters (N_D = 5e18 cm^-3 in 30 nm AlGaAs)
    over-deplete the modulation-doped barrier so V_T_HEMT falls in
    the deeply-negative range [-3, -2] V; the device is strongly
    normally-on across the full [0, 1] V sweep, which is the simplest
    operating regime for a 2DEG-vs-V_GS V&V comparison (the channel
    is in linear regime throughout, no sub-threshold roll-off to
    model). A realistic depletion-mode HEMT lives in this part of
    the parameter space; an enhancement-mode HEMT would lower N_D
    (e.g. 1e18 cm^-3) to push V_T closer to zero."""
    p = _algaas_gaas_params()
    V_T = hemt_threshold_voltage(**p)
    assert -3.0 < V_T < -2.0, (
        f"V_T_HEMT {V_T:.3f} V outside expected [-3.0, -2.0] V band "
        f"for the M17 benchmark parameters"
    )


def test_2deg_zero_below_threshold():
    p = _algaas_gaas_params()
    V_T = hemt_threshold_voltage(**p)
    n_s = hemt_2deg_classical(V_T - 0.5, **p)
    assert n_s == pytest.approx(0.0)


def test_2deg_linear_regime_in_textbook_range():
    """At V_GS = 1.0 V (well above threshold) the classical 2DEG
    sheet density is in the textbook AlGaAs/GaAs HEMT range, which
    spans roughly 1e12 - 1e13 cm^-2 depending on modulation-doping
    density and barrier thickness."""
    p = _algaas_gaas_params()
    n_s = hemt_2deg_classical(1.0, **p)
    assert 1.0e12 < n_s < 1.0e13, (
        f"n_s(V_GS=1.0) = {n_s:.3e} cm^-2 outside the textbook "
        f"AlGaAs/GaAs HEMT range (1e12 - 1e13)"
    )


def test_2deg_si_units_consistency():
    """The cm^-2 and m^-2 return values agree to the conversion factor."""
    p = _algaas_gaas_params()
    V_GS = 0.7
    n_s_cm2 = hemt_2deg_classical(V_GS, **p)
    n_s_m2 = hemt_2deg_classical_si_units(V_GS, **p)
    assert n_s_m2 == pytest.approx(n_s_cm2 * 1.0e4, rel=1.0e-12)


def test_2deg_linear_in_v_gs_above_threshold():
    """The classical formula is linear in V_GS - V_T. Test by sampling
    two points and checking the ratio."""
    p = _algaas_gaas_params()
    V_T = hemt_threshold_voltage(**p)
    n_s_a = hemt_2deg_classical(V_T + 0.3, **p)
    n_s_b = hemt_2deg_classical(V_T + 0.6, **p)
    assert n_s_b == pytest.approx(2.0 * n_s_a, rel=1.0e-9)


def test_2deg_slope_matches_capacitance_per_unit_area():
    """The slope d(n_s) / d(V_GS) = eps_AlGaAs eps_0 / (q d_barrier),
    which is the parallel-plate capacitance per unit area divided by
    q. For 12.0 * 8.854e-12 / (1.6e-19 * 30e-9) ~= 2.21e21 m^-2/V or
    2.21e17 cm^-2/V."""
    p = _algaas_gaas_params()
    V_T = hemt_threshold_voltage(**p)
    # Sample slope between V_T+0.4 and V_T+0.5 (linear regime).
    n_s_lo = hemt_2deg_classical(V_T + 0.4, **p)
    n_s_hi = hemt_2deg_classical(V_T + 0.5, **p)
    measured_slope_per_cm2 = (n_s_hi - n_s_lo) / 0.1  # cm^-2 / V
    expected_slope_per_m2 = (
        p["eps_r_barrier"] * EPS0 / (Q * p["d_barrier_m"])
    )
    expected_slope_per_cm2 = expected_slope_per_m2 * 1.0e-4
    assert measured_slope_per_cm2 == pytest.approx(
        expected_slope_per_cm2, rel=1.0e-9
    )


def test_invalid_d_barrier_raises():
    p = _algaas_gaas_params()
    p["d_barrier_m"] = 0.0
    with pytest.raises(ValueError):
        hemt_threshold_voltage(**p)


def test_invalid_eps_r_raises():
    p = _algaas_gaas_params()
    p["eps_r_barrier"] = 0.0
    with pytest.raises(ValueError):
        hemt_threshold_voltage(**p)
