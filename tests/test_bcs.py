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
