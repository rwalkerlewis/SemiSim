"""
Pure-Python unit tests for the M17 heterojunction-field machinery.

The dolfinx-dependent FEM-side coverage lives in
`tests/fem/test_heterojunction_assembly.py`; this file exercises the
material-resolution logic and the `Material.n_i_at_T` helper, neither
of which touch dolfinx.
"""
from __future__ import annotations

import math

import pytest

from semi.constants import cm3_to_m3, thermal_voltage
from semi.materials import MATERIALS, Material, get_material
from semi.physics.heterojunction import (
    _resolve_reference_chi_eV,
    _resolve_region_material,
    cfg_uses_heterojunction,
)
from semi.scaling import Scaling


def test_resolve_no_overrides_returns_database_material():
    """Without `material_overrides`, the resolver returns the database
    material instance directly. v0.23.0 byte-identity hinges on this."""
    region = {"material": "Si", "tag": 1, "role": "semiconductor"}
    resolved = _resolve_region_material(region)
    assert resolved is MATERIALS["Si"]


def test_resolve_with_chi_override_returns_new_material():
    """A chi_eV override produces a new Material instance with the
    chi field replaced; the database material is unchanged."""
    region = {
        "material": "GaAs",
        "tag": 1,
        "role": "semiconductor",
        "material_overrides": {"chi_eV": 3.5},
    }
    resolved = _resolve_region_material(region)
    assert resolved is not MATERIALS["GaAs"]
    assert resolved.chi == pytest.approx(3.5)
    assert resolved.Eg == pytest.approx(MATERIALS["GaAs"].Eg)
    assert resolved.Nc == pytest.approx(MATERIALS["GaAs"].Nc)
    assert resolved.Nv == pytest.approx(MATERIALS["GaAs"].Nv)
    # The override on chi alone does not touch n_i (chi is not in the
    # n_i closed-form), so n_i stays at the database value.
    assert resolved.n_i == pytest.approx(MATERIALS["GaAs"].n_i)
    assert MATERIALS["GaAs"].chi == pytest.approx(4.07)


def test_resolve_with_eg_override_recomputes_n_i():
    """An Eg_eV override forces an n_i recompute via the closed-form
    sqrt(Nc Nv) exp(-Eg / (2 V_t)) so the override is self-consistent."""
    region = {
        "material": "Si",
        "tag": 1,
        "role": "semiconductor",
        "material_overrides": {"Eg_eV": 1.50},
    }
    resolved = _resolve_region_material(region)
    assert resolved.Eg == pytest.approx(1.50)
    V_t_300 = thermal_voltage(300.0)
    expected_n_i = math.sqrt(MATERIALS["Si"].Nc * MATERIALS["Si"].Nv) * math.exp(
        -1.50 / (2.0 * V_t_300)
    )
    assert resolved.n_i == pytest.approx(expected_n_i, rel=1e-12)


def test_resolve_with_nc_override_recomputes_n_i():
    """An Nc_per_cm3 override recomputes n_i as well."""
    region = {
        "material": "GaAs",
        "tag": 1,
        "role": "semiconductor",
        "material_overrides": {"Nc_per_cm3": 1.0e18},
    }
    resolved = _resolve_region_material(region)
    assert resolved.Nc == pytest.approx(cm3_to_m3(1.0e18))
    V_t_300 = thermal_voltage(300.0)
    expected_n_i = math.sqrt(resolved.Nc * MATERIALS["GaAs"].Nv) * math.exp(
        -MATERIALS["GaAs"].Eg / (2.0 * V_t_300)
    )
    assert resolved.n_i == pytest.approx(expected_n_i, rel=1e-12)


def test_resolve_full_override_set_propagates():
    """All four override fields can be set together and every one
    propagates into the resolved material."""
    region = {
        "material": "GaAs",
        "tag": 1,
        "role": "semiconductor",
        "material_overrides": {
            "chi_eV": 3.74,
            "Eg_eV": 1.798,
            "Nc_per_cm3": 5.5e17,
            "Nv_per_cm3": 1.0e19,
        },
    }
    resolved = _resolve_region_material(region)
    assert resolved.chi == pytest.approx(3.74)
    assert resolved.Eg == pytest.approx(1.798)
    assert resolved.Nc == pytest.approx(cm3_to_m3(5.5e17))
    assert resolved.Nv == pytest.approx(cm3_to_m3(1.0e19))


def test_resolve_does_not_mutate_database():
    """Repeated overrides do not accumulate on the global MATERIALS
    dict. After resolving with overrides, the database material is
    untouched (verified by reading `chi` after the fact)."""
    original_chi = MATERIALS["Si"].chi
    region = {
        "material": "Si",
        "tag": 1,
        "role": "semiconductor",
        "material_overrides": {"chi_eV": 9.99},
    }
    _resolve_region_material(region)
    assert MATERIALS["Si"].chi == pytest.approx(original_chi)
    assert MATERIALS["Si"].chi != pytest.approx(9.99)


def test_resolve_unknown_material_raises():
    region = {
        "material": "Unobtainium",
        "tag": 1,
        "role": "semiconductor",
    }
    with pytest.raises(KeyError):
        _resolve_region_material(region)


def test_n_i_at_T_returns_zero_for_insulator():
    """Insulators leave Nc / Nv / Eg at zero; n_i_at_T returns 0.0
    rather than raising or producing NaN."""
    sio2 = MATERIALS["SiO2"]
    assert sio2.n_i_at_T(300.0) == 0.0


def test_n_i_at_T_si_300k_within_tolerance():
    """For Si at 300 K the closed-form formula returns ~1.16e10 cm^-3
    while the empirical `Material.n_i` is the Altermatt 2003 consensus
    1.0e10 cm^-3. The two agree to within ~25 %; consumers that need
    the empirical value read `n_i` directly."""
    si = MATERIALS["Si"]
    formula = si.n_i_at_T(300.0)
    empirical = si.n_i
    ratio = formula / empirical
    assert 0.8 < ratio < 1.3, (
        f"n_i_at_T(300 K) for Si returned {formula:.3e} m^-3 "
        f"but empirical is {empirical:.3e} m^-3 (ratio {ratio:.2f})"
    )


def test_n_i_at_T_temperature_scaling():
    """n_i_at_T at T = 400 K is much larger than at T = 300 K because
    the Boltzmann factor relaxes (smaller Eg / V_t)."""
    si = MATERIALS["Si"]
    n_i_300 = si.n_i_at_T(300.0)
    n_i_400 = si.n_i_at_T(400.0)
    assert n_i_400 > n_i_300
    # Roughly two orders of magnitude rise from 300 K to 400 K for Si.
    assert n_i_400 / n_i_300 > 50.0


def test_n_i_at_T_algaas_smaller_than_gaas():
    """AlGaAs_0p3 has a wider gap (1.798 eV) than GaAs (1.424 eV), so
    its n_i at 300 K is exponentially smaller. This is the load-bearing
    physics of the HEMT 2DEG: the AlGaAs barrier is effectively
    intrinsic-free at room temperature relative to the GaAs channel."""
    n_i_gaas = MATERIALS["GaAs"].n_i_at_T(300.0)
    n_i_algaas = MATERIALS["AlGaAs_0p3"].n_i_at_T(300.0)
    assert n_i_algaas < n_i_gaas
    # Wider gap by ~0.37 eV -> exp(0.37 / 0.052) ~= 1100x suppression.
    assert n_i_gaas / n_i_algaas > 100.0


def test_algaas_0p3_in_database():
    """The AlGaAs_0p3 entry is registered and looks up cleanly via the
    standard `get_material` helper used by the rest of the engine."""
    mat = get_material("AlGaAs_0p3")
    assert mat.role == "semiconductor"
    assert mat.Eg == pytest.approx(1.798)
    assert mat.chi == pytest.approx(3.74)
    assert mat.epsilon_r == pytest.approx(12.0)
    assert mat.Nc == pytest.approx(cm3_to_m3(5.5e17))
    assert mat.Nv == pytest.approx(cm3_to_m3(1.0e19))
    assert mat.m_n_star is not None


def test_algaas_0p3_n_i_consistent_with_band_parameters():
    """The stored AlGaAs_0p3 n_i (1.85e3 cm^-3) is close to the
    closed-form value sqrt(Nc Nv) exp(-Eg / (2 V_t)) at 300 K. Unlike
    Si (where n_i is the empirical 1.0e10 cm^-3), AlGaAs_0p3 has no
    entrenched empirical value, so we keep n_i and the formula
    consistent at the database level."""
    mat = MATERIALS["AlGaAs_0p3"]
    formula = mat.n_i_at_T(300.0)
    ratio = formula / mat.n_i
    assert 0.8 < ratio < 1.3, (
        f"AlGaAs_0p3 n_i_at_T(300 K) {formula:.3e} differs from stored "
        f"n_i {mat.n_i:.3e} (ratio {ratio:.2f})"
    )


def test_material_n_i_at_T_method_signature():
    """The Material dataclass exposes n_i_at_T as an instance method
    accepting a single positional T argument. Smoke check: calling
    on a fresh in-test Material instance does not require a global
    state."""
    fake = Material(
        name="Fake",
        role="semiconductor",
        epsilon_r=10.0,
        Eg=1.0,
        chi=4.0,
        Nc=cm3_to_m3(1.0e19),
        Nv=cm3_to_m3(1.0e19),
        n_i=cm3_to_m3(1.0e10),
        mu_n=0.1,
        mu_p=0.05,
    )
    V_t = thermal_voltage(300.0)
    expected = math.sqrt(fake.Nc * fake.Nv) * math.exp(-fake.Eg / (2.0 * V_t))
    assert fake.n_i_at_T(300.0) == pytest.approx(expected, rel=1e-12)


def test_unknown_override_field_silently_ignored():
    """Forward-compat: an unknown override key (rejected at validate
    time by the schema additionalProperties:false gate) is silently
    ignored at the resolver level. The schema-test suite exercises
    the rejection path; this test guards the resolver against
    crashing on a hand-built cfg that bypassed validate()."""
    region = {
        "material": "Si",
        "tag": 1,
        "role": "semiconductor",
        "material_overrides": {"chi_eV": 4.05, "BogusUnused": 1.0},
    }
    resolved = _resolve_region_material(region)
    assert resolved.chi == pytest.approx(4.05)


def test_cfg_uses_heterojunction_rejects_non_dict_input():
    """Hand-built cfgs that bypass schema.validate may pass a list or
    None to the heterojunction-detection helper; the gate returns False
    rather than raising so the runners' fast-path branch stays selected.
    """
    assert cfg_uses_heterojunction(None) is False
    assert cfg_uses_heterojunction([]) is False
    assert cfg_uses_heterojunction("oops") is False


def test_cfg_uses_heterojunction_skips_non_dict_region_entries():
    """Iteration over a regions map that mixes dict and non-dict values
    skips the non-dict entries instead of raising. Single-material
    config with a stray non-dict entry still reports False."""
    cfg = {"left": {"material": "Si"}, "noise": "ignored"}
    assert cfg_uses_heterojunction(cfg) is False


def test_cfg_uses_heterojunction_returns_true_on_overrides():
    cfg = {
        "left": {"material": "Si"},
        "right": {"material": "GaAs", "material_overrides": {"chi_eV": 3.5}},
    }
    assert cfg_uses_heterojunction(cfg) is True


def test_cfg_uses_heterojunction_returns_true_on_explicit_flag():
    cfg = {
        "left": {"material": "Si", "heterojunction": True},
        "right": {"material": "GaAs"},
    }
    assert cfg_uses_heterojunction(cfg) is True


def test_cfg_uses_heterojunction_returns_false_on_empty_overrides():
    """An explicit `material_overrides: {}` is functionally a no-op and
    keeps the byte-identity branch active."""
    cfg = {"left": {"material": "Si", "material_overrides": {}}}
    assert cfg_uses_heterojunction(cfg) is False


def _scaling_for_material(mat: Material) -> Scaling:
    return Scaling(
        L0=1.0e-6,
        C0=1.0e23,
        T=300.0,
        mu0=mat.mu_n,
        n_i=mat.n_i,
        N_C=mat.Nc,
        N_V=mat.Nv,
        m_n_star=mat.m_n_star,
        m_p_star=mat.m_p_star,
        E_g=mat.Eg,
    )


def test_resolve_reference_chi_picks_first_semiconductor_region():
    """_resolve_reference_chi_eV walks the regions map in insertion order
    and returns chi (eV) for the first entry whose role is semiconductor
    with a positive chi. Insulators are skipped."""
    sc = _scaling_for_material(MATERIALS["Si"])
    cfg = {
        "oxide": {"material": "SiO2", "role": "insulator"},
        "channel": {"material": "GaAs", "role": "semiconductor"},
    }
    chi = _resolve_reference_chi_eV(cfg, sc)
    assert chi == pytest.approx(MATERIALS["GaAs"].chi)


def test_resolve_reference_chi_skips_non_dict_region_entries():
    """Non-dict entries in the region map (hand-built cfgs) are skipped
    rather than raising."""
    sc = _scaling_for_material(MATERIALS["Si"])
    cfg = {
        "garbage": "not-a-dict",
        "channel": {"material": "Si", "role": "semiconductor"},
    }
    chi = _resolve_reference_chi_eV(cfg, sc)
    assert chi == pytest.approx(MATERIALS["Si"].chi)


def test_resolve_reference_chi_falls_back_to_silicon_when_all_insulator():
    """When no semiconductor region exists in the cfg, the helper falls
    back to Si.chi so downstream BC math still has a sane reference."""
    sc = _scaling_for_material(MATERIALS["Si"])
    cfg = {
        "oxide_top": {"material": "SiO2", "role": "insulator"},
        "oxide_bot": {"material": "SiO2", "role": "insulator"},
    }
    chi = _resolve_reference_chi_eV(cfg, sc)
    assert chi == pytest.approx(MATERIALS["Si"].chi)
