"""
FEM smoke tests for the M16.6 BBT and TAT tunneling kernels in
``build_dd_block_residual``.

The acceptance gate proper is the analytical Kane reverse-breakdown
benchmark in `benchmarks/zener_1d/`; these tests just verify the
form assembles in every flag combination, the BBT and TAT branches
are exercised in the gated docker-fem-tests coverage job (matching
the M16.5 lesson learned that schottky pre-solve coverage required
a separate FEM unit test; the zener_1d benchmark CI job runs in a
separate matrix entry whose coverage is not merged into the gated
suite).
"""
from __future__ import annotations

import numpy as np

from semi.bcs import build_dd_dirichlet_bcs, resolve_contacts
from semi.doping import build_profile
from semi.materials import get_material
from semi.mesh import build_mesh
from semi.scaling import make_scaling_from_config


def _zener_1d_cfg(*, bbt: bool = True, tat: bool = True):
    """Trivial 1D heavy-doped abrupt junction cfg with tunneling on."""
    return {
        "name": "zener_1d_smoke",
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
                "profile": {
                    "type": "step", "axis": 0, "location": 5.0e-7,
                    "N_D_left": 0.0, "N_A_left": 1.0e19,
                    "N_D_right": 1.0e19, "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {"name": "anode", "facet": "anode", "type": "ohmic",
             "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "fermi_dirac",
            "tunneling": {
                "bbt": bool(bbt),
                "tat": bool(tat),
                "A_kane": 4.0e14,
                "B_kane": 1.9e7,
                "F_kT": 1.4e7,
                "alpha": 2.0,
            },
        },
        "solver": {"type": "bias_sweep"},
    }


def _build_residual(cfg, *, recomb_cfg):
    """Helper to assemble the DD residual."""
    from dolfinx import fem

    from semi.physics.drift_diffusion import (
        build_dd_block_residual,
        make_dd_block_spaces,
    )

    ref_mat = get_material("Si")
    sc = make_scaling_from_config(cfg, ref_mat)
    msh, _cell_tags, facet_tags = build_mesh(cfg)

    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi
    N_raw_fn = build_profile(cfg["doping"])
    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    spaces.psi.interpolate(
        lambda x: np.arcsinh(N_raw_fn(x) / (2.0 * ref_mat.n_i))
    )
    spaces.phi_n.x.array[:] = 0.0
    spaces.phi_p.x.array[:] = 0.0
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.scatter_forward()

    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_dd_dirichlet_bcs(
        spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
    )
    F_list = build_dd_block_residual(
        spaces, N_hat_fn, sc, ref_mat.epsilon_r,
        1.0, 1.0, 1.0e-1, 1.0e-1, 0.0,
        mobility_cfg={"model": "constant", "mu_n": 1400.0, "mu_p": 450.0},
        facet_tags=facet_tags,
        recomb_cfg=recomb_cfg,
        statistics_cfg={"statistics": "fermi_dirac"},
    )
    return F_list, spaces, bcs


def test_tunneling_off_path_assembles():
    """Both flags off: bit-identical pre-M16.6 form (Auger/SRH only)."""
    cfg = _zener_1d_cfg(bbt=False, tat=False)
    recomb_cfg = {
        "srh": True, "tau_n": 1.0e-7, "tau_p": 1.0e-7,
        "E_t": 0.0, "auger": False,
        "bbt": False, "tat": False,
    }
    F_list, _spaces, _bcs = _build_residual(cfg, recomb_cfg=recomb_cfg)
    assert len(F_list) == 3


def test_bbt_on_path_assembles():
    """bbt=True exercises the Kane band-to-band branch in the form."""
    cfg = _zener_1d_cfg(bbt=True, tat=False)
    recomb_cfg = {
        "srh": True, "tau_n": 1.0e-7, "tau_p": 1.0e-7,
        "E_t": 0.0, "auger": False,
        "bbt": True, "tat": False,
        "A_kane": 4.0e14, "B_kane": 1.9e7,
    }
    F_list, _spaces, _bcs = _build_residual(cfg, recomb_cfg=recomb_cfg)
    assert len(F_list) == 3


def test_tat_on_path_assembles():
    """tat=True exercises the Hurkx TAT enhancement branch."""
    cfg = _zener_1d_cfg(bbt=False, tat=True)
    recomb_cfg = {
        "srh": True, "tau_n": 1.0e-7, "tau_p": 1.0e-7,
        "E_t": 0.0, "auger": False,
        "bbt": False, "tat": True,
        "F_kT": 1.4e7, "alpha": 2.0,
    }
    F_list, _spaces, _bcs = _build_residual(cfg, recomb_cfg=recomb_cfg)
    assert len(F_list) == 3


def test_both_flags_on_path_assembles():
    """Both flags on: BBT generation subtracted from R_SRH * (1 + Gamma)."""
    cfg = _zener_1d_cfg(bbt=True, tat=True)
    recomb_cfg = {
        "srh": True, "tau_n": 1.0e-7, "tau_p": 1.0e-7,
        "E_t": 0.0, "auger": False,
        "bbt": True, "tat": True,
        "A_kane": 4.0e14, "B_kane": 1.9e7,
        "F_kT": 1.4e7, "alpha": 2.0,
    }
    F_list, _spaces, _bcs = _build_residual(cfg, recomb_cfg=recomb_cfg)
    assert len(F_list) == 3
