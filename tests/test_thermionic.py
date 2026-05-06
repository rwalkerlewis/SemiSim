"""
Pure-Python tests for the M16.5 thermionic-emission helpers exposed on
`Scaling` and `Material`. The Schottky surface form itself lives in
`semi/physics/drift_diffusion.py` (UFL) and is exercised by the FEM
smoke test in `tests/fem/test_schottky_surface.py` and end-to-end by
the `benchmarks/schottky_1d` verifier.
"""
from __future__ import annotations

import math

import pytest

from semi.constants import KB, M0, Q
from semi.materials import get_material
from semi.scaling import Scaling


def _v_th_textbook(m_star_rel: float, T: float = 300.0) -> float:
    return math.sqrt(KB * T / (2.0 * math.pi * m_star_rel * M0))


def test_v_n_thermal_si_at_300K_matches_richardson_form():
    """`sc.v_n_thermal` returns the effective Richardson velocity
    ``v_R = sqrt(kT / (2 pi m_n*))`` (Sze 3rd ed Section 3.4 derivation
    of the thermionic-emission boundary condition; equivalent to
    ``A* T^2 / (q N_C)``). Distinct from the Maxwell-Boltzmann
    average thermal speed (``sqrt(8 kT / (pi m*))``) which is also
    sometimes called the "thermal velocity" in textbooks."""
    sc = Scaling(L0=1.0e-6, C0=1.0e22, T=300.0, mu0=1.0e-1, n_i=1.0e16,
                 m_n_star=0.26, m_p_star=0.39)
    expected = _v_th_textbook(0.26, T=300.0)
    assert sc.v_n_thermal == pytest.approx(expected, rel=1e-12)
    # Si Richardson velocity at 300 K with m_n* = 0.26 m_0:
    # ~5.3e4 m/s = ~5.3e6 cm/s.
    assert sc.v_n_thermal == pytest.approx(5.27e4, rel=0.02)


def test_v_p_thermal_si_at_300K_matches_richardson_form():
    sc = Scaling(L0=1.0e-6, C0=1.0e22, T=300.0, mu0=1.0e-1, n_i=1.0e16,
                 m_n_star=0.26, m_p_star=0.39)
    expected = _v_th_textbook(0.39, T=300.0)
    assert sc.v_p_thermal == pytest.approx(expected, rel=1e-12)
    # Si hole Richardson velocity at 300 K with m_p* = 0.39 m_0:
    # ~4.3e4 m/s = ~4.3e6 cm/s.
    assert sc.v_p_thermal == pytest.approx(4.31e4, rel=0.02)


def test_v_n_thermal_raises_when_m_star_unset():
    sc = Scaling(L0=1.0e-6, C0=1.0e22, T=300.0, mu0=1.0e-1, n_i=1.0e16)
    with pytest.raises(ValueError, match="m_n_star"):
        _ = sc.v_n_thermal


def test_v_p_thermal_raises_when_m_star_unset():
    sc = Scaling(L0=1.0e-6, C0=1.0e22, T=300.0, mu0=1.0e-1, n_i=1.0e16,
                 m_n_star=0.26)
    with pytest.raises(ValueError, match="m_p_star"):
        _ = sc.v_p_thermal


def test_n_eq_metal_supply_closed_form():
    """
    Closed-form metal-side equilibrium electron density at three barrier
    heights, computed directly from N_C and phi_B per Sze eq 3.4.5
    (n_eq = N_C exp(-phi_B / V_t)). Verifies the relation that the
    thermionic-emission Robin form embeds.
    """
    Si = get_material("Si")
    V_t = KB * 300.0 / Q
    cases = [
        # (phi_B in eV, expected n_eq in m^-3, comment)
        (0.85, Si.Nc * math.exp(-0.85 / V_t), "Pt-on-n-Si"),
        (0.70, Si.Nc * math.exp(-0.70 / V_t), "moderate barrier"),
        (0.10, Si.Nc * math.exp(-0.10 / V_t), "low barrier"),
    ]
    for phi_B, expected, _ in cases:
        observed = Si.Nc * math.exp(-phi_B / V_t)
        assert observed == pytest.approx(expected, rel=0.0, abs=1e-12)


def test_n_eq_reduces_to_nc_at_zero_barrier():
    """In the phi_B -> 0 limit (perfect ohmic), n_eq -> N_C. Beyond N_D,
    the Schottky boundary effectively floods with electrons; the
    thermionic-emission rate becomes much faster than diffusion and the
    contact behaves ohmically."""
    Si = get_material("Si")
    V_t = KB * 300.0 / Q
    n_eq_phi_zero = Si.Nc * math.exp(-0.0 / V_t)
    assert n_eq_phi_zero == pytest.approx(Si.Nc, rel=1e-12)


def test_si_material_carries_thermionic_masses():
    Si = get_material("Si")
    assert Si.m_n_star == pytest.approx(0.26, rel=1e-9)
    assert Si.m_p_star == pytest.approx(0.39, rel=1e-9)


def test_other_materials_default_thermionic_masses_to_none():
    # Materials that have not been characterized for Schottky contacts
    # carry None placeholders so the failure surfaces clearly when a
    # benchmark requests a thermionic-emission BC.
    for name in ("SiO2", "HfO2", "Si3N4"):
        mat = get_material(name)
        assert mat.m_n_star is None
        assert mat.m_p_star is None


def test_scaling_threads_thermionic_masses_from_material():
    """`make_scaling_from_config` populates `m_n_star` / `m_p_star` from
    the reference material; this is what `_build_schottky_surface_forms`
    consumes via `sc.v_n_thermal`."""
    from semi.scaling import make_scaling_from_config

    cfg = {
        "schema_version": "2.5.0",
        "name": "test",
        "dimension": 1,
        "mesh": {"source": "builtin",
                 "extents": [[0.0, 1.0e-6]],
                 "resolution": [10]},
        "regions": {"si": {"material": "Si"}},
        "doping": [
            {"region": "si",
             "profile": {"type": "uniform", "N_D": 1.0e16, "N_A": 0.0}},
        ],
        "contacts": [],
        "physics": {"temperature": 300.0},
    }
    Si = get_material("Si")
    sc = make_scaling_from_config(cfg, Si)
    assert sc.m_n_star == pytest.approx(0.26)
    assert sc.m_p_star == pytest.approx(0.39)
    # And the derived thermal velocity is positive.
    assert sc.v_n_thermal > 0.0
    assert sc.v_p_thermal > 0.0
