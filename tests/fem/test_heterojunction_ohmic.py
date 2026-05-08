"""
M17 Phase D: ohmic-contact equilibrium psi reads local chi.

These tests exercise the BC-builder integration with the heterojunction
field machinery on a tiny 1D mesh; they run only inside the
Dockerized FEM environment.
"""
from __future__ import annotations

import numpy as np
import pytest


def _build_two_region_1d(L: float = 1.0e-6, N: int = 20):
    """1D mesh tagged 1 on left half (Si), 2 on right half (GaAs)."""
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    msh = dmesh.create_interval(MPI.COMM_WORLD, N, [0.0, L])
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim, 0)
    msh.topology.create_entities(tdim)

    from dolfinx import mesh as _dmesh
    n_local = msh.topology.index_map(tdim).size_local
    cells = np.arange(n_local, dtype=np.int32)
    midpoints = _dmesh.compute_midpoints(msh, tdim, cells)[:, :1]
    tags = np.where(midpoints[:, 0] < 0.5 * L, 1, 2).astype(np.int32)
    cell_tags = _dmesh.meshtags(msh, tdim, cells, tags)
    return msh, cell_tags


def _build_left_right_facets(msh, L: float = 1.0e-6):
    """Tag facet 1 at x=0, facet 2 at x=L."""
    from dolfinx import mesh as dmesh
    fdim = msh.topology.dim - 1
    msh.topology.create_entities(fdim)

    def left(x):
        return np.isclose(x[0], 0.0)

    def right(x):
        return np.isclose(x[0], L)

    left_facets = dmesh.locate_entities_boundary(msh, fdim, left)
    right_facets = dmesh.locate_entities_boundary(msh, fdim, right)
    indices = np.concatenate([left_facets, right_facets]).astype(np.int32)
    values = np.concatenate([
        np.ones_like(left_facets, dtype=np.int32),
        2 * np.ones_like(right_facets, dtype=np.int32),
    ])
    facet_tags = dmesh.meshtags(msh, fdim, indices, values)
    return facet_tags


def test_psi_dirichlet_local_chi_no_op_for_single_material():
    """Two-region cfg with the same material on both regions and no
    overrides: the BC builder must produce psi values bit-identical
    to the v0.23.0 path (chi_local == chi_ref everywhere, n_i_local
    == n_i_ref, so the Anderson-rule shift is zero)."""
    from dolfinx import fem

    from semi.bcs import build_psi_dirichlet_bcs, resolve_contacts
    from semi.materials import get_material
    from semi.scaling import Scaling

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d(L=L)
    facet_tags = _build_left_right_facets(msh, L=L)
    Si = get_material("Si")
    sc = Scaling(L0=L, C0=1.0e23, T=300.0,
                 mu0=Si.mu_n, n_i=Si.n_i,
                 N_C=Si.Nc, N_V=Si.Nv,
                 m_n_star=Si.m_n_star, m_p_star=Si.m_p_star,
                 E_g=Si.Eg)
    cfg = {
        "mesh": {"facets_by_plane": [
            {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
            {"name": "cathode", "tag": 2, "axis": 0, "value": L},
        ]},
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "regions": {
            "left":  {"material": "Si", "tag": 1, "role": "semiconductor"},
            "right": {"material": "Si", "tag": 2, "role": "semiconductor"},
        },
    }
    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    bcs_baseline = build_psi_dirichlet_bcs(
        V_psi, msh, facet_tags, contacts, sc, Si, lambda x: np.full(x.shape[1], 1.0e23),
    )
    bcs_local = build_psi_dirichlet_bcs(
        V_psi, msh, facet_tags, contacts, sc, Si, lambda x: np.full(x.shape[1], 1.0e23),
        regions_cfg=cfg["regions"], cell_tags=cell_tags,
    )
    # Both BC lists carry the same Dirichlet values to byte precision.
    for bc_a, bc_b in zip(bcs_baseline, bcs_local, strict=True):
        assert np.allclose(bc_a.g.value, bc_b.g.value, atol=1.0e-14, rtol=0.0)


def test_psi_dirichlet_anderson_shift_on_heterojunction():
    """Two-region heterojunction (Si + AlGaAs_0p3) with ohmic contacts
    on both ends: the BC builder applies the Anderson-rule chi shift
    on the right contact (in the AlGaAs region). The shift magnitude
    matches `(chi_AlGaAs - chi_Si) / V_t`."""
    from dolfinx import fem

    from semi.bcs import build_psi_dirichlet_bcs, resolve_contacts
    from semi.materials import get_material
    from semi.scaling import Scaling

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d(L=L)
    facet_tags = _build_left_right_facets(msh, L=L)
    Si = get_material("Si")
    AlGaAs = get_material("AlGaAs_0p3")
    sc = Scaling(L0=L, C0=1.0e23, T=300.0,
                 mu0=Si.mu_n, n_i=Si.n_i,
                 N_C=Si.Nc, N_V=Si.Nv,
                 m_n_star=Si.m_n_star, m_p_star=Si.m_p_star,
                 E_g=Si.Eg)
    cfg = {
        "mesh": {"facets_by_plane": [
            {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
            {"name": "cathode", "tag": 2, "axis": 0, "value": L},
        ]},
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "regions": {
            "left":  {"material": "Si",         "tag": 1, "role": "semiconductor"},
            "right": {"material": "AlGaAs_0p3", "tag": 2, "role": "semiconductor",
                      "heterojunction": True},
        },
    }
    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    bcs_local = build_psi_dirichlet_bcs(
        V_psi, msh, facet_tags, contacts, sc, Si, lambda x: np.full(x.shape[1], 1.0e23),
        regions_cfg=cfg["regions"], cell_tags=cell_tags,
    )
    # First contact is on the Si region (left); BC value is the
    # baseline psi_eq_hat(Si). Second contact is in AlGaAs; BC carries
    # the chi shift plus the asinh term with the AlGaAs n_i.
    bc_anode_value = float(bcs_local[0].g.value)
    bc_cathode_value = float(bcs_local[1].g.value)
    # Anode (Si region) is exactly the v0.23.0 value.
    expected_anode = float(np.arcsinh(1.0e23 / (2.0 * Si.n_i)))
    assert bc_anode_value == pytest.approx(expected_anode, rel=1e-9)
    # Cathode (AlGaAs region) is shifted by chi_AlGaAs - chi_Si and
    # uses the AlGaAs n_i.
    expected_chi_shift = (AlGaAs.chi - Si.chi) / sc.V0
    expected_cathode = expected_chi_shift + float(
        np.arcsinh(1.0e23 / (2.0 * AlGaAs.n_i))
    )
    assert bc_cathode_value == pytest.approx(expected_cathode, rel=1e-6)
