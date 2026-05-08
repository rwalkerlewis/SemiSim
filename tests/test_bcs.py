"""
Pure-Python tests for `semi.bcs.resolve_contacts` and `ContactBC`.

These exercise only the JSON-config parsing, kind validation, facet
resolution, and the bias-sweep voltages override. The actual dolfinx
BC construction is covered in `tests/fem/test_bcs.py`.
"""
from __future__ import annotations

import pytest

from semi.bcs import ContactBC, resolve_contacts


def _cfg_one_ohmic():
    return {
        "mesh": {
            "facets_by_plane": [
                {"name": "left",  "tag": 1, "axis": 0, "value": 0.0},
                {"name": "right", "tag": 2, "axis": 0, "value": 1.0},
            ],
        },
        "contacts": [
            {"name": "anode", "type": "ohmic", "facet": "left", "voltage": 0.3},
        ],
    }


def _cfg_two_ohmic_with_sweep():
    return {
        "mesh": {
            "facets_by_plane": [
                {"name": "left",  "tag": 1, "axis": 0, "value": 0.0},
                {"name": "right", "tag": 2, "axis": 0, "value": 1.0},
            ],
        },
        "contacts": [
            {"name": "anode",   "type": "ohmic", "facet": "left",  "voltage": 0.0},
            {
                "name": "cathode", "type": "ohmic", "facet": "right",
                "voltage": 0.0,
                "voltage_sweep": {"start": 0.0, "stop": 0.6, "step": 0.05},
            },
        ],
    }


class _MockFacetTags:
    """Minimal stand-in for dolfinx.mesh.MeshTags exposing `.find(tag)`."""

    def __init__(self, tag_to_facets: dict[int, list[int]]):
        self._t2f = tag_to_facets

    def find(self, tag: int):
        return self._t2f.get(int(tag), [])


def test_resolve_contacts_single_ohmic():
    cfg = _cfg_one_ohmic()
    out = resolve_contacts(cfg)
    assert len(out) == 1
    bc = out[0]
    assert isinstance(bc, ContactBC)
    assert bc.name == "anode"
    assert bc.kind == "ohmic"
    assert bc.facet_tag == 1
    assert bc.V_applied == pytest.approx(0.3)
    assert bc.work_function is None


def test_resolve_contacts_two_contacts_default_voltages():
    cfg = _cfg_two_ohmic_with_sweep()
    out = resolve_contacts(cfg)
    assert [c.name for c in out] == ["anode", "cathode"]
    assert [c.facet_tag for c in out] == [1, 2]
    # Default voltage on both is 0.0 (no sweep override applied yet).
    assert out[0].V_applied == pytest.approx(0.0)
    assert out[1].V_applied == pytest.approx(0.0)


def test_resolve_contacts_unknown_kind_raises():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["type"] = "magic_unicorn"
    with pytest.raises(ValueError, match="Unknown contact kind"):
        resolve_contacts(cfg)


def test_resolve_contacts_unknown_facet_name_raises():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["facet"] = "ceiling"
    with pytest.raises(RuntimeError, match="not declared under mesh.facets_by_plane"):
        resolve_contacts(cfg)


def test_resolve_contacts_missing_mesh_facet_raises():
    cfg = _cfg_one_ohmic()
    facet_tags = _MockFacetTags({2: [10, 11]})  # tag 1 has no facets
    with pytest.raises(RuntimeError, match="no .*facets carry that tag"):
        resolve_contacts(cfg, facet_tags=facet_tags)


def test_resolve_contacts_facet_tags_pass_through():
    cfg = _cfg_two_ohmic_with_sweep()
    facet_tags = _MockFacetTags({1: [0], 2: [99]})
    out = resolve_contacts(cfg, facet_tags=facet_tags)
    # Successful path: facet_tags satisfied, no exception, both contacts back.
    assert [c.facet_tag for c in out] == [1, 2]


def test_resolve_contacts_voltages_override_updates_V_applied():
    cfg = _cfg_two_ohmic_with_sweep()
    # Simulate the bias sweep at three ramp points.
    for V_try in (0.10, 0.25, 0.55):
        out = resolve_contacts(cfg, voltages={"cathode": V_try})
        cathode = next(c for c in out if c.name == "cathode")
        anode = next(c for c in out if c.name == "anode")
        assert cathode.V_applied == pytest.approx(V_try), (
            f"voltages override should update V_applied to {V_try}, "
            f"got {cathode.V_applied}"
        )
        # Anode is unaffected when not present in the override dict.
        assert anode.V_applied == pytest.approx(0.0)


def test_resolve_contacts_skips_insulating():
    cfg = _cfg_one_ohmic()
    cfg["contacts"].append(
        {"name": "side", "type": "insulating", "facet": "right"}
    )
    out = resolve_contacts(cfg)
    # Insulating contact does not produce a Dirichlet ContactBC.
    assert [c.name for c in out] == ["anode"]


def test_resolve_contacts_integer_facet_ref():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["facet"] = 7  # raw integer tag, no name lookup
    out = resolve_contacts(cfg)
    assert out[0].facet_tag == 7


def test_resolve_contacts_gate_kind_carries_workfunction():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["type"] = "gate"
    cfg["contacts"][0]["workfunction"] = 4.65
    out = resolve_contacts(cfg)
    assert out[0].kind == "gate"
    assert out[0].work_function == pytest.approx(4.65)


def test_resolve_contacts_schottky_carries_barrier_height():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["type"] = "schottky"
    cfg["contacts"][0]["barrier_height_eV"] = 0.85
    out = resolve_contacts(cfg)
    assert out[0].kind == "schottky"
    assert out[0].barrier_height_eV == pytest.approx(0.85)


def test_resolve_contacts_non_schottky_leaves_barrier_none():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["barrier_height_eV"] = 0.85  # ignored on ohmic
    out = resolve_contacts(cfg)
    assert out[0].kind == "ohmic"
    assert out[0].barrier_height_eV is None


def test_resolve_contacts_gate_leaves_barrier_none():
    cfg = _cfg_one_ohmic()
    cfg["contacts"][0]["type"] = "gate"
    cfg["contacts"][0]["workfunction"] = 4.05
    out = resolve_contacts(cfg)
    assert out[0].kind == "gate"
    assert out[0].barrier_height_eV is None


def test_schottky_psi_eq_closed_form():
    """`_schottky_psi_eq` reduces to ln(N_C / n_i) - phi_B / V_t in
    scaled units; manual calculation matched to within 1e-12.
    """
    import math

    from semi.bcs import ContactBC, _schottky_psi_eq

    class _RefMat:
        def __init__(self, name, n_i, Nc):
            self.name = name
            self.n_i = n_i
            self.Nc = Nc

    class _Scaling:
        def __init__(self, V0):
            self.V0 = V0

    sc = _Scaling(V0=0.025852)  # kT/q at 300 K
    cases = [
        ("Si_Pt", 1.0e10 * 1.0e6, 2.86e19 * 1.0e6, 0.85),
        ("Si_Au", 1.0e10 * 1.0e6, 2.86e19 * 1.0e6, 0.80),
        ("GaAs", 2.1e6 * 1.0e6, 4.7e17 * 1.0e6, 0.70),
    ]
    for name, n_i, Nc, phi_b in cases:
        ref_mat = _RefMat(name=name, n_i=n_i, Nc=Nc)
        c = ContactBC(
            name="anode", kind="schottky", facet_tag=1, V_applied=0.0,
            barrier_height_eV=phi_b,
        )
        observed = _schottky_psi_eq(c, ref_mat, sc)
        expected = math.log(Nc / n_i) - phi_b / sc.V0
        assert observed == pytest.approx(expected, rel=0.0, abs=1e-12), (
            f"case {name!r}: expected {expected}, got {observed}"
        )


def test_schottky_psi_eq_missing_barrier_raises():
    from semi.bcs import ContactBC, _schottky_psi_eq

    class _RefMat:
        name = "Si"
        n_i = 1.0e16
        Nc = 2.86e25

    class _Scaling:
        V0 = 0.025852

    c = ContactBC(
        name="anode", kind="schottky", facet_tag=1, V_applied=0.0,
        barrier_height_eV=None,
    )
    with pytest.raises(ValueError, match="barrier_height_eV"):
        _schottky_psi_eq(c, _RefMat(), _Scaling())


def test_schottky_psi_eq_zero_nc_raises():
    from semi.bcs import ContactBC, _schottky_psi_eq

    class _RefMat:
        name = "Insulator"
        n_i = 1.0e16
        Nc = 0.0

    class _Scaling:
        V0 = 0.025852

    c = ContactBC(
        name="anode", kind="schottky", facet_tag=1, V_applied=0.0,
        barrier_height_eV=0.85,
    )
    with pytest.raises(ValueError, match="N_C"):
        _schottky_psi_eq(c, _RefMat(), _Scaling())


def test_contact_bc_dataclass_equality_and_repr():
    a = ContactBC(name="anode", kind="ohmic", facet_tag=1,
                  V_applied=0.3, work_function=None)
    b = ContactBC(name="anode", kind="ohmic", facet_tag=1,
                  V_applied=0.3, work_function=None)
    c = ContactBC(name="anode", kind="ohmic", facet_tag=1,
                  V_applied=0.4, work_function=None)
    assert a == b
    assert a != c
    r = repr(a)
    assert "ContactBC" in r
    assert "anode" in r
    assert "ohmic" in r


def test_contact_bc_carries_barrier_height_field():
    a = ContactBC(name="anode", kind="schottky", facet_tag=1,
                  V_applied=0.0, work_function=None,
                  barrier_height_eV=0.85)
    assert a.barrier_height_eV == pytest.approx(0.85)
    assert a.kind == "schottky"


# M17 Phase D: ohmic-equilibrium psi with local-material chi shift.

def _scaling_for_si_test():
    from semi.materials import get_material
    from semi.scaling import Scaling
    Si = get_material("Si")
    return Si, Scaling(
        L0=1.0e-6, C0=1.0e23, T=300.0,
        mu0=Si.mu_n, n_i=Si.n_i, N_C=Si.Nc, N_V=Si.Nv,
        m_n_star=Si.m_n_star, m_p_star=Si.m_p_star, E_g=Si.Eg,
    )


def test_ohmic_psi_eq_hat_no_overrides_matches_v0_23_0():
    """When `local_chi_eV` and `local_n_i` are both None, the helper
    returns the v0.23.0 single-material formula bit-identically."""
    import math

    import numpy as np

    from semi.bcs import _ohmic_psi_eq_hat

    Si, sc = _scaling_for_si_test()
    N_net = 1.0e23
    expected = float(np.arcsinh(N_net / (2.0 * Si.n_i)))
    got = _ohmic_psi_eq_hat(N_net, Si, sc)
    assert math.isclose(got, expected, rel_tol=1e-12, abs_tol=0.0)


def test_ohmic_psi_eq_hat_chi_shift_local_equals_ref():
    """Passing the reference material's own chi / n_i as the local
    values must collapse to the v0.23.0 formula. Regression gate for
    byte-identity on every existing benchmark (no benchmark uses
    overrides)."""
    import math

    import numpy as np

    from semi.bcs import _ohmic_psi_eq_hat

    Si, sc = _scaling_for_si_test()
    N_net = 5.0e22
    expected = float(np.arcsinh(N_net / (2.0 * Si.n_i)))
    got = _ohmic_psi_eq_hat(
        N_net, Si, sc,
        local_chi_eV=Si.chi, local_n_i=Si.n_i,
    )
    assert math.isclose(got, expected, rel_tol=1e-12, abs_tol=0.0)


def test_ohmic_psi_eq_hat_chi_shift_anderson_rule():
    """A 0.5 eV shift in the local material's chi translates into a
    `0.5 / V_t` additive shift in the dimensionless psi. This is the
    Anderson-rule band alignment (M17, ADR 0016)."""
    import math

    from semi.bcs import _ohmic_psi_eq_hat

    Si, sc = _scaling_for_si_test()
    N_net = 5.0e22
    baseline = _ohmic_psi_eq_hat(N_net, Si, sc)
    chi_local = Si.chi - 0.5
    shifted = _ohmic_psi_eq_hat(
        N_net, Si, sc,
        local_chi_eV=chi_local, local_n_i=Si.n_i,
    )
    expected_shift = (chi_local - Si.chi) / sc.V0
    assert math.isclose(shifted - baseline, expected_shift,
                        rel_tol=1e-12, abs_tol=0.0)
    # The chi-shift contribution must be many V_t units (a 0.5 eV shift
    # at 300 K is ~19 V_t units), so this is not a numerical trivium.
    assert abs(shifted - baseline) > 5.0


def test_ohmic_psi_eq_hat_local_n_i_change():
    """Lowering local n_i (wider gap, e.g. AlGaAs) increases the
    `asinh(N_net / (2 n_i))` term, biasing the contact more strongly
    in the same doping environment."""
    from semi.bcs import _ohmic_psi_eq_hat

    Si, sc = _scaling_for_si_test()
    N_net = 1.0e22
    baseline = _ohmic_psi_eq_hat(N_net, Si, sc)
    # AlGaAs_0p3-ish: ~7 orders of magnitude smaller n_i.
    shifted = _ohmic_psi_eq_hat(
        N_net, Si, sc,
        local_chi_eV=Si.chi, local_n_i=Si.n_i * 1.0e-7,
    )
    assert shifted > baseline
    # The shift is `asinh(N / (2 n_i_smaller)) - asinh(N / (2 n_i))`,
    # which grows logarithmically; for a 1e7 ratio in n_i the shift
    # is on the order of `ln(1e7)` ~= 16 (in V_t units).
    assert (shifted - baseline) > 10.0
