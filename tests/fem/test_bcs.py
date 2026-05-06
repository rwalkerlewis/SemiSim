"""
FEM tests for `semi.bcs.build_psi_dirichlet_bcs` and
`semi.bcs.build_dd_dirichlet_bcs`.

These compare the extracted helpers against a copy of the legacy
inline implementation that previously lived in `semi/run.py`. Keeping
the reference embedded here means the assertion stays meaningful even
after the inline copies are removed in the follow-up commit.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.bcs import (
    build_dd_dirichlet_bcs,
    build_psi_dirichlet_bcs,
    resolve_contacts,
)
from semi.doping import build_profile
from semi.materials import get_material
from semi.mesh import build_mesh
from semi.scaling import make_scaling_from_config


def _legacy_evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn) -> float:
    """Verbatim copy of the pre-extraction helper from semi/run.py."""
    tdim = msh.topology.dim
    msh.topology.create_connectivity(fdim, 0)
    f2v = msh.topology.connectivity(fdim, 0)
    coords = msh.geometry.x
    verts = f2v.links(int(facets[0]))
    centroid = coords[verts, :tdim].mean(axis=0)
    pt = centroid.reshape(tdim, 1)
    return float(N_raw_fn(pt)[0])


def _legacy_build_ohmic_bcs_psi(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn):
    """Verbatim copy of the pre-extraction inline helper from semi/run.py."""
    from dolfinx import fem
    from petsc4py import PETSc

    fdim = msh.topology.dim - 1
    tag_by_name = {p["name"]: int(p["tag"])
                   for p in cfg["mesh"].get("facets_by_plane", [])}
    bcs = []
    for contact in cfg["contacts"]:
        if contact["type"] != "ohmic":
            continue
        facet_ref = contact["facet"]
        tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)
        facets = facet_tags.find(tag)
        N_net = _legacy_evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
        psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
        V_applied = contact.get("voltage", 0.0)
        psi_bc = psi_eq_hat + V_applied / sc.V0
        dofs = fem.locate_dofs_topological(V, fdim, facets)
        bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs, V))
    return bcs


def _legacy_build_dd_ohmic_bcs(cfg, spaces, msh, facet_tags, sc, ref_mat,
                               N_raw_fn, voltages):
    """Verbatim copy of the pre-extraction inline DD helper from semi/run.py."""
    from dolfinx import fem
    from petsc4py import PETSc

    fdim = msh.topology.dim - 1
    tag_by_name = {p["name"]: int(p["tag"])
                   for p in cfg["mesh"].get("facets_by_plane", [])}
    bcs = []
    for contact in cfg["contacts"]:
        if contact["type"] != "ohmic":
            continue
        name = contact["name"]
        facet_ref = contact["facet"]
        tag = tag_by_name[facet_ref] if isinstance(facet_ref, str) else int(facet_ref)
        facets = facet_tags.find(tag)
        N_net = _legacy_evaluate_doping_at_facet(msh, facets, fdim, N_raw_fn)
        psi_eq_hat = float(np.arcsinh(N_net / (2.0 * ref_mat.n_i)))
        V_app = float(voltages.get(name, contact.get("voltage", 0.0)))
        V_hat = V_app / sc.V0
        psi_bc = psi_eq_hat + V_hat
        phi_bc = V_hat
        dofs_psi = fem.locate_dofs_topological(spaces.V_psi, fdim, facets)
        dofs_n = fem.locate_dofs_topological(spaces.V_phi_n, fdim, facets)
        dofs_p = fem.locate_dofs_topological(spaces.V_phi_p, fdim, facets)
        bcs.append(fem.dirichletbc(PETSc.ScalarType(psi_bc), dofs_psi, spaces.V_psi))
        bcs.append(fem.dirichletbc(PETSc.ScalarType(phi_bc), dofs_n, spaces.V_phi_n))
        bcs.append(fem.dirichletbc(PETSc.ScalarType(phi_bc), dofs_p, spaces.V_phi_p))
    return bcs


def _pn_1d_cfg(*, anode_voltage: float = 0.0, cathode_voltage: float = 0.0):
    """A minimal symmetric 1D pn junction config sufficient for BC tests."""
    return {
        "name": "pn_1d_bc_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 2.0e-6]],
            "resolution": [40],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 2.0e-6]]}
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 2.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step",
                    "axis": 0,
                    "location": 1.0e-6,
                    "N_D_left": 0.0,
                    "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic",
             "voltage": float(anode_voltage)},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": float(cathode_voltage)},
        ],
        "physics": {"temperature": 300.0, "statistics": "boltzmann"},
        "solver": {"type": "equilibrium"},
    }


def _bc_value(bc):
    """Return the (scalar) Dirichlet value attached to a dolfinx bc."""
    arr = np.asarray(bc.g.value)
    return float(arr.flat[0])


def _bc_dofs(bc):
    """Return the sorted DOF index array attached to a dolfinx bc."""
    return np.sort(np.asarray(bc.dof_indices()[0], dtype=np.int64))


def test_build_psi_dirichlet_bcs_matches_legacy_at_equilibrium():
    from dolfinx import fem

    cfg = _pn_1d_cfg(anode_voltage=0.0, cathode_voltage=0.0)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    bcs_old = _legacy_build_ohmic_bcs_psi(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn)
    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs_new = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    assert len(bcs_old) == len(bcs_new) == 2
    for bc_old, bc_new in zip(bcs_old, bcs_new, strict=True):
        assert _bc_value(bc_new) == pytest.approx(_bc_value(bc_old), rel=0.0, abs=1e-14)
        assert np.array_equal(_bc_dofs(bc_new), _bc_dofs(bc_old))


def test_build_psi_dirichlet_bcs_matches_legacy_at_forward_bias():
    from dolfinx import fem

    # Anode at +0.4 V exercises the V_applied / V_t branch.
    cfg = _pn_1d_cfg(anode_voltage=0.4, cathode_voltage=0.0)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    bcs_old = _legacy_build_ohmic_bcs_psi(cfg, V, msh, facet_tags, sc, ref_mat, N_raw_fn)
    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs_new = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    for bc_old, bc_new in zip(bcs_old, bcs_new, strict=True):
        assert _bc_value(bc_new) == pytest.approx(_bc_value(bc_old), rel=0.0, abs=1e-14)
        assert np.array_equal(_bc_dofs(bc_new), _bc_dofs(bc_old))


def _mos_like_cfg(V_gate: float = 0.0, phi_ms: float = 0.0):
    """
    Minimal two-region 2D config with one ohmic body contact and one
    gate contact. Used to cover the gate branch in build_psi_dirichlet_bcs
    and to confirm build_dd_dirichlet_bcs skips the gate entirely.
    """
    return {
        "name": "mos_bc_test",
        "dimension": 2,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-7], [0.0, 1.0e-7]],
            "resolution": [8, 10],
            "regions_by_box": [
                {"name": "silicon", "tag": 1,
                 "bounds": [[0.0, 1.0e-7], [0.0, 0.7e-7]]},
                {"name": "oxide", "tag": 2,
                 "bounds": [[0.0, 1.0e-7], [0.7e-7, 1.0e-7]]},
            ],
            "facets_by_plane": [
                {"name": "body", "tag": 1, "axis": 1, "value": 0.0},
                {"name": "gate", "tag": 2, "axis": 1, "value": 1.0e-7},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
            "oxide":   {"material": "SiO2", "tag": 2, "role": "insulator"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17},
            },
        ],
        "contacts": [
            {"name": "body", "facet": "body", "type": "ohmic", "voltage": 0.0},
            {"name": "gate", "facet": "gate", "type": "gate",
             "voltage": float(V_gate), "workfunction": float(phi_ms)},
        ],
        "physics": {"temperature": 300.0, "statistics": "boltzmann"},
        "solver": {"type": "equilibrium"},
    }


def test_build_psi_dirichlet_bcs_gate_zero_workfunction():
    """V_gate = 0, phi_ms = 0 -> psi_hat = 0 at the gate facet."""
    from dolfinx import fem

    cfg = _mos_like_cfg(V_gate=0.0, phi_ms=0.0)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    # One ohmic body BC + one gate BC.
    assert len(bcs) == 2
    gate_bc = bcs[1]
    assert _bc_value(gate_bc) == pytest.approx(0.0, abs=1e-14)


def test_build_psi_dirichlet_bcs_gate_applied_voltage_scaled_by_Vt():
    """V_gate = 1.5 V, phi_ms = 0 -> psi_hat = 1.5 / V_t."""
    from dolfinx import fem

    cfg = _mos_like_cfg(V_gate=1.5, phi_ms=0.0)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    gate_bc = bcs[1]
    expected = 1.5 / sc.V0
    assert _bc_value(gate_bc) == pytest.approx(expected, rel=1e-12)


def test_build_psi_dirichlet_bcs_gate_honors_workfunction():
    """psi_hat(gate) = (V_gate - phi_ms) / V_t."""
    from dolfinx import fem

    V_gate = 1.0
    phi_ms = 0.25
    cfg = _mos_like_cfg(V_gate=V_gate, phi_ms=phi_ms)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    gate_bc = bcs[1]
    expected = (V_gate - phi_ms) / sc.V0
    assert _bc_value(gate_bc) == pytest.approx(expected, rel=1e-12)


def test_build_dd_dirichlet_bcs_includes_gate_psi_only():
    """Gate contact contributes a psi-only Dirichlet (phi_n / phi_p skipped).

    Pre-M14.3, build_dd_dirichlet_bcs skipped gate contacts entirely on
    the rationale that gates sit on the oxide side and the Slotboom
    continuity blocks live on the semiconductor submesh. That left
    bias_sweep with no way to apply a gate BC during a sweep, which
    broke the M14.3 mosfet_2d Pao-Sah verifier (every step reported
    SNES iterations=0 because the BCs were unchanged across V_gate
    values). Fix: gate contacts now contribute a single psi
    Dirichlet here (phi_n / phi_p remain skipped). The expected
    BC count for a body-ohmic + gate config is therefore 3 (body
    psi/phi_n/phi_p) + 1 (gate psi) = 4.
    """
    from semi.physics.drift_diffusion import make_dd_block_spaces

    V_gate = 0.5
    cfg = _mos_like_cfg(V_gate=V_gate, phi_ms=0.0)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    spaces = make_dd_block_spaces(msh)
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_dd_dirichlet_bcs(spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    assert len(bcs) == 4

    # build_dd_dirichlet_bcs processes contacts in config order:
    # body (ohmic): bcs[0]=psi, bcs[1]=phi_n, bcs[2]=phi_p
    # gate (gate):  bcs[3]=psi  ← the gate psi Dirichlet
    # Verify the gate BC has the expected scaled value without relying on
    # dolfinx FunctionSpace Python-object identity (bc.function_space is
    # not guaranteed to be the same object across dolfinx versions).
    gate_bc = bcs[3]
    expected = V_gate / sc.V0
    assert _bc_value(gate_bc) == pytest.approx(expected, rel=1e-12)


def test_build_dd_dirichlet_bcs_gate_voltage_overrides_via_resolve():
    """When `voltages={gate_name: V_new}` is passed to resolve_contacts,
    the gate psi Dirichlet must reflect V_new (regression test for the
    M14.3 mosfet_2d failure where gate BC was never updated mid-sweep).
    """
    from semi.physics.drift_diffusion import make_dd_block_spaces

    cfg = _mos_like_cfg(V_gate=0.0, phi_ms=0.0)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    spaces = make_dd_block_spaces(msh)
    N_raw_fn = build_profile(cfg["doping"])

    V_swept = 0.7
    contacts = resolve_contacts(
        cfg, facet_tags=facet_tags, voltages={"gate": V_swept},
    )
    bcs = build_dd_dirichlet_bcs(spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    # Gate psi is bcs[3] (contacts processed in config order: body ohmic
    # contributes bcs[0-2], gate contributes bcs[3]).
    assert len(bcs) == 4
    gate_bc = bcs[3]
    expected = V_swept / sc.V0
    assert _bc_value(gate_bc) == pytest.approx(expected, rel=1e-12)


def _schottky_1d_cfg(*, anode_voltage: float = 0.0,
                     barrier_height_eV: float = 0.85):
    """Trivial 1D mesh with a Schottky anode and an ohmic cathode."""
    return {
        "name": "schottky_1d_bc_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [40],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6]]}
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {"type": "uniform", "N_D": 1.0e16, "N_A": 0.0},
            }
        ],
        "contacts": [
            {"name": "anode", "facet": "anode", "type": "schottky",
             "barrier_height_eV": float(barrier_height_eV),
             "voltage": float(anode_voltage)},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "physics": {"temperature": 300.0, "statistics": "boltzmann"},
        "solver": {"type": "equilibrium"},
    }


def test_build_psi_dirichlet_bcs_schottky_value_matches_metal_fermi_level():
    """Psi at a Schottky facet equals
    ln(N_C / n_i) - phi_B / V_t + V_applied / V_t (scaled units).
    """
    import math

    from dolfinx import fem

    cfg = _schottky_1d_cfg(anode_voltage=0.0, barrier_height_eV=0.85)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    # Two contacts -> two psi BCs (Schottky anode + ohmic cathode).
    assert len(bcs) == 2
    schottky_bc = bcs[0]
    expected = math.log(ref_mat.Nc / ref_mat.n_i) - 0.85 / sc.V0
    assert _bc_value(schottky_bc) == pytest.approx(expected, rel=1e-12)


def test_build_psi_dirichlet_bcs_schottky_applied_bias_shifts_psi():
    """Forward bias V_a > 0 shifts psi at the Schottky facet by V_a / V_t."""
    import math

    from dolfinx import fem

    V_anode = 0.3
    cfg = _schottky_1d_cfg(anode_voltage=V_anode, barrier_height_eV=0.85)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn)

    schottky_bc = bcs[0]
    expected = math.log(ref_mat.Nc / ref_mat.n_i) - 0.85 / sc.V0 + V_anode / sc.V0
    assert _bc_value(schottky_bc) == pytest.approx(expected, rel=1e-12)


def test_build_dd_dirichlet_bcs_schottky_psi_only_no_phi_bcs():
    """A Schottky contact contributes a single psi Dirichlet; phi_n
    and phi_p rows are not pinned at the Schottky facet (they are
    handled by the Phase C Robin form on the continuity equations).
    The cathode-only ohmic contact contributes 3 BCs; the anode
    Schottky contact contributes 1 (psi only).
    """
    import math

    from semi.physics.drift_diffusion import make_dd_block_spaces

    cfg = _schottky_1d_cfg(anode_voltage=0.0, barrier_height_eV=0.85)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    spaces = make_dd_block_spaces(msh)
    N_raw_fn = build_profile(cfg["doping"])

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_dd_dirichlet_bcs(
        spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
    )

    # Anode (schottky): bcs[0] = psi.
    # Cathode (ohmic): bcs[1] = psi, bcs[2] = phi_n, bcs[3] = phi_p.
    assert len(bcs) == 4
    schottky_bc = bcs[0]
    expected = math.log(ref_mat.Nc / ref_mat.n_i) - 0.85 / sc.V0
    assert _bc_value(schottky_bc) == pytest.approx(expected, rel=1e-12)


def test_build_dd_dirichlet_bcs_schottky_voltage_override_via_resolve():
    """Mid-sweep override updates the Schottky psi BC."""
    import math

    from semi.physics.drift_diffusion import make_dd_block_spaces

    cfg = _schottky_1d_cfg(anode_voltage=0.0, barrier_height_eV=0.85)
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    spaces = make_dd_block_spaces(msh)
    N_raw_fn = build_profile(cfg["doping"])

    V_anode = 0.4
    contacts = resolve_contacts(
        cfg, facet_tags=facet_tags, voltages={"anode": V_anode},
    )
    bcs = build_dd_dirichlet_bcs(
        spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
    )

    assert len(bcs) == 4
    schottky_bc = bcs[0]
    expected = math.log(ref_mat.Nc / ref_mat.n_i) - 0.85 / sc.V0 + V_anode / sc.V0
    assert _bc_value(schottky_bc) == pytest.approx(expected, rel=1e-12)


def test_build_dd_dirichlet_bcs_matches_legacy_at_bias_step():
    from semi.physics.drift_diffusion import make_dd_block_spaces

    cfg = _pn_1d_cfg()
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    spaces = make_dd_block_spaces(msh)
    N_raw_fn = build_profile(cfg["doping"])

    # Mid-ramp bias step: cathode held at 0, anode pushed to +0.25 V via the
    # voltages override (the same code path the bias-sweep driver uses).
    voltages = {"anode": 0.25, "cathode": 0.0}

    bcs_old = _legacy_build_dd_ohmic_bcs(
        cfg, spaces, msh, facet_tags, sc, ref_mat, N_raw_fn, voltages,
    )
    contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=voltages)
    bcs_new = build_dd_dirichlet_bcs(
        spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
    )

    # Each ohmic contact contributes three Dirichlet entries
    # (psi, phi_n, phi_p), so two contacts -> six BCs total.
    assert len(bcs_old) == len(bcs_new) == 6
    for bc_old, bc_new in zip(bcs_old, bcs_new, strict=True):
        assert _bc_value(bc_new) == pytest.approx(_bc_value(bc_old), rel=0.0, abs=1e-14)
        assert np.array_equal(_bc_dofs(bc_new), _bc_dofs(bc_old))


def test_build_dd_dirichlet_bcs_skips_insulating_contacts():
    """Contacts with kind 'insulating' produce no Dirichlet entries.

    Exercises the ``continue`` branch inside ``build_dd_dirichlet_bcs``.
    Insulating contacts are dropped at ``resolve_contacts`` time, but the
    builder is also defensive against an explicit insulating ContactBC.
    """
    from semi.bcs import ContactBC
    from semi.physics.drift_diffusion import make_dd_block_spaces

    cfg = _pn_1d_cfg()
    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)
    spaces = make_dd_block_spaces(msh)
    N_raw_fn = build_profile(cfg["doping"])

    ohmic_contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    insulating = ContactBC(
        name="fake_insulating", kind="insulating", facet_tag=999,
        V_applied=0.0,
    )
    contacts = [insulating] + list(ohmic_contacts)

    bcs = build_dd_dirichlet_bcs(
        spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
    )
    # Only the two ohmic contacts (3 BCs each) should appear; the
    # insulating contact is silently skipped.
    assert len(bcs) == 6
