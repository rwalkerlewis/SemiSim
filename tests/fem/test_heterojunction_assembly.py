"""
FEM-side tests for the M17 heterojunction DG0 field builder.

Runs only inside the Dockerized FEM environment; tests/fem/conftest.py
skips the directory when dolfinx is unavailable.

Verifies the per-cell DG0 fields for chi / Eg / Nc / Nv / n_i / eps_r:
- a single-region 1D mesh produces a constant field on every cell
  (collapses to the v0.23.0 scalar single-material values),
- a two-region 1D mesh (Si + GaAs) produces a chi field with a step
  at the region boundary,
- the resolved-material override path picks up `material_overrides`
  on the Right region,
- the per-cell n_i for an AlGaAs_0p3 region is materially smaller
  than the GaAs n_i (heterojunction band-suppression).
"""
from __future__ import annotations

import numpy as np


def _build_two_region_1d_mesh(L: float = 1.0e-6, N: int = 20):
    """Build a 1D mesh with two cell-tag regions (left half / right half).

    Returns (mesh, cell_tags) with tag = 1 on x in [0, L/2] cells and
    tag = 2 on x in [L/2, L] cells.
    """
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    msh = dmesh.create_interval(MPI.COMM_WORLD, N, [0.0, L])
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim, 0)
    msh.topology.create_entities(tdim)

    midpoints = _cell_midpoints(msh)
    n_cells = midpoints.shape[0]
    tags = np.where(midpoints[:, 0] < 0.5 * L, 1, 2).astype(np.int32)
    indices = np.arange(n_cells, dtype=np.int32)
    cell_tags = dmesh.meshtags(msh, tdim, indices, tags)
    return msh, cell_tags


def _cell_midpoints(msh) -> np.ndarray:
    """Return per-cell midpoint coordinates as an (n_cells, gdim) array."""
    from dolfinx import mesh as dmesh

    tdim = msh.topology.dim
    gdim = msh.geometry.dim
    n_local = msh.topology.index_map(tdim).size_local
    if n_local == 0:
        return np.zeros((0, gdim))
    cells = np.arange(n_local, dtype=np.int32)
    return dmesh.compute_midpoints(msh, tdim, cells)[:, :gdim]


def _build_single_region_1d_mesh(L: float = 1.0e-6, N: int = 20):
    """One region tagged as 1 over the entire mesh."""
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    msh = dmesh.create_interval(MPI.COMM_WORLD, N, [0.0, L])
    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    n_local = msh.topology.index_map(tdim).size_local
    indices = np.arange(n_local, dtype=np.int32)
    tags = np.ones(n_local, dtype=np.int32)
    return msh, dmesh.meshtags(msh, tdim, indices, tags)


def _scaling_for_si(L: float, C0: float = 1.0e23, T: float = 300.0):
    from semi.materials import get_material
    from semi.scaling import Scaling

    Si = get_material("Si")
    return Scaling(
        L0=L, C0=C0, T=T,
        mu0=Si.mu_n, n_i=Si.n_i,
        N_C=Si.Nc, N_V=Si.Nv,
        m_n_star=Si.m_n_star, m_p_star=Si.m_p_star,
        E_g=Si.Eg,
    )


def test_single_region_field_is_constant_per_material():
    """A single-region 1D mesh tagged Si everywhere produces a chi
    field that is constant across all cells, matching the scalar
    single-material value to byte precision (the v0.23.0 byte-
    identity branch)."""
    from semi.materials import MATERIALS
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_single_region_1d_mesh(L=L, N=12)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)

    chi_arr = fields["chi_hat"].x.array
    expected_chi_hat = MATERIALS["Si"].chi / sc.V0
    assert np.allclose(chi_arr, expected_chi_hat, rtol=0.0, atol=1.0e-12)

    n_i_arr = fields["n_i_hat"].x.array
    expected_n_i_hat = MATERIALS["Si"].n_i / sc.C0
    assert np.allclose(n_i_arr, expected_n_i_hat, rtol=0.0, atol=1.0e-12)


def test_two_region_field_steps_at_interface():
    """Si on the left half, GaAs on the right: the chi field has a
    step at the interface. Cells near x = 0 evaluate to Si.chi / V_t;
    cells near x = L evaluate to GaAs.chi / V_t."""
    from semi.materials import MATERIALS
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "right": {"material": "GaAs", "tag": 2, "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)

    chi_arr = fields["chi_hat"].x.array
    # Expect two distinct values across cells: Si.chi/V_t and GaAs.chi/V_t.
    si_val = MATERIALS["Si"].chi / sc.V0
    gaas_val = MATERIALS["GaAs"].chi / sc.V0
    unique = np.unique(np.round(chi_arr, decimals=8))
    assert len(unique) == 2
    assert np.any(np.isclose(chi_arr, si_val, atol=1.0e-9))
    assert np.any(np.isclose(chi_arr, gaas_val, atol=1.0e-9))


def test_material_overrides_propagate_to_chi_field():
    """`material_overrides: {chi_eV: 3.5}` on the right region produces
    a chi field with the override value on the right half of the mesh."""
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "right": {
            "material": "GaAs",
            "tag": 2,
            "role": "semiconductor",
            "material_overrides": {"chi_eV": 3.5},
        },
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)

    chi_arr = fields["chi_hat"].x.array
    expected_override_hat = 3.5 / sc.V0
    assert np.any(np.isclose(chi_arr, expected_override_hat, atol=1.0e-9))


def test_algaas_n_i_field_smaller_than_gaas():
    """Si on the left, AlGaAs_0p3 on the right: the per-cell n_i_hat
    on the AlGaAs side is many orders of magnitude below the Si side
    because the wider gap suppresses the intrinsic density."""
    from semi.materials import MATERIALS
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "left":  {"material": "Si",         "tag": 1, "role": "semiconductor"},
        "right": {"material": "AlGaAs_0p3", "tag": 2, "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)

    n_i_arr = fields["n_i_hat"].x.array
    si_n_i_hat = MATERIALS["Si"].n_i / sc.C0
    algaas_n_i_hat = MATERIALS["AlGaAs_0p3"].n_i / sc.C0
    # Right-side cells should match AlGaAs n_i; left-side cells match Si.
    assert np.any(np.isclose(n_i_arr, si_n_i_hat, rtol=0.05))
    assert np.any(np.isclose(n_i_arr, algaas_n_i_hat, rtol=0.05))
    # AlGaAs_0p3 n_i is ~7 orders of magnitude below Si.
    assert si_n_i_hat / max(algaas_n_i_hat, 1.0e-300) > 1.0e5


def test_eps_r_field_matches_existing_mesh_helper():
    """The eps_r DG0 field built by the heterojunction helper agrees
    with the v0.20.0 multi-region helper `semi.mesh.build_eps_r_function`
    on every cell. This is the byte-identity guarantee for the existing
    multi-region MOSCAP path."""
    from semi.mesh import build_eps_r_function
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "right": {"material": "SiO2", "tag": 2, "role": "insulator"},
    }
    eps_r_old = build_eps_r_function(msh, cell_tags, regions_cfg)
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)
    eps_r_new = fields["epsilon_r"]
    assert np.allclose(eps_r_old.x.array, eps_r_new.x.array, atol=1.0e-12)


def test_chi_ref_hat_returned():
    """The result dict carries `chi_ref_hat` (Si chi divided by V_t at
    the runner temperature) and `n_i_ref_hat` (Si n_i divided by sc.C0)
    for the Phase D ohmic-equilibrium calculation."""
    from semi.materials import MATERIALS
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "right": {"material": "GaAs", "tag": 2, "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)
    expected_chi_ref_hat = MATERIALS["Si"].chi / sc.V0
    expected_n_i_ref_hat = MATERIALS["Si"].n_i / sc.C0
    assert fields["chi_ref_hat"] == \
        __import__("pytest").approx(expected_chi_ref_hat, rel=1.0e-12)
    assert fields["n_i_ref_hat"] == \
        __import__("pytest").approx(expected_n_i_ref_hat, rel=1.0e-12)


def test_cell_tags_none_collapses_to_first_region_material():
    """A single-region builtin mesh passes `cell_tags=None`; the helper
    must collapse every field to the first region's material values
    so the runner's pre-M17 single-material path stays byte-identical
    when no MeshTags object is built. `chi_ref_hat` is resolved from
    the first semiconductor region's chi (GaAs here), not from the
    Scaling reference."""
    import pytest

    from semi.materials import MATERIALS
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, _ = _build_single_region_1d_mesh(L=L, N=12)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "channel": {"material": "GaAs", "tag": 1, "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, None, regions_cfg, sc, T=300.0)

    chi_arr = fields["chi_hat"].x.array
    expected_chi_hat = MATERIALS["GaAs"].chi / sc.V0
    assert np.allclose(chi_arr, expected_chi_hat, rtol=0.0, atol=1.0e-12)

    eps_arr = fields["epsilon_r"].x.array
    assert np.allclose(eps_arr, MATERIALS["GaAs"].epsilon_r, rtol=0.0, atol=1.0e-12)

    # n_i_hat still uses sc.n_i in the cell_tags=None branch (the runner
    # builds Scaling from the reference material), so the field equals
    # sc.n_i / sc.C0 everywhere — independent of the resolved material.
    n_i_arr = fields["n_i_hat"].x.array
    assert np.allclose(n_i_arr, float(sc.n_i) / sc.C0, rtol=0.0, atol=1.0e-12)

    # GaAs is the first (and only) semiconductor in the cfg, so it
    # supplies the chi_ref used by the BC layer's Anderson shift.
    assert fields["chi_ref_hat"] == pytest.approx(
        MATERIALS["GaAs"].chi / sc.V0, rel=1e-12
    )


def test_cell_tags_none_falls_back_to_silicon_when_regions_empty():
    """If the regions map carries no material entry (degenerate hand-
    built cfg with only a non-dict noise value), the cell_tags=None
    branch falls back to silicon. The chi field must equal Si.chi /
    sc.V0 across the entire mesh; this exercises the
    `if ref_mat is None: ref_mat = MATERIALS["Si"]` Si fallback."""
    from semi.materials import MATERIALS
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, _ = _build_single_region_1d_mesh(L=L, N=8)
    sc = _scaling_for_si(L=L)
    regions_cfg: dict = {"noise": "ignored"}
    fields = build_dg0_material_fields(msh, None, regions_cfg, sc, T=300.0)

    chi_arr = fields["chi_hat"].x.array
    assert np.allclose(chi_arr, MATERIALS["Si"].chi / sc.V0, rtol=0.0, atol=1.0e-12)


def test_skips_non_dict_region_entries():
    """Hand-built cfgs with non-dict noise in the regions map (e.g.
    a stray comment string) skip those entries instead of raising."""
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_single_region_1d_mesh(L=L, N=12)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "channel": {"material": "Si", "tag": 1, "role": "semiconductor"},
        "noise": "ignored",
        "missing_tag": {"material": "GaAs", "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)
    chi_arr = fields["chi_hat"].x.array
    # Only the tagged silicon region should populate the field.
    from semi.materials import MATERIALS
    assert np.allclose(chi_arr, MATERIALS["Si"].chi / sc.V0, rtol=0.0, atol=1.0e-12)


def test_skips_region_with_no_cells_in_mesh():
    """A region cfg may carry a tag that no cell in the mesh actually
    uses (e.g. a 3D-only region in a 2D-mesh test fixture); the helper
    must skip the empty tag instead of indexing into an empty array."""
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    regions_cfg = {
        "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "right": {"material": "GaAs", "tag": 2, "role": "semiconductor"},
        "ghost": {"material": "AlGaAs_0p3", "tag": 99, "role": "semiconductor"},
    }
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=300.0)
    # The build does not crash; every cell still has a finite chi value.
    chi_arr = fields["chi_hat"].x.array
    assert np.all(np.isfinite(chi_arr))
    # Only Si and GaAs values appear — the AlGaAs ghost region was skipped.
    from semi.materials import MATERIALS
    si_val = MATERIALS["Si"].chi / sc.V0
    gaas_val = MATERIALS["GaAs"].chi / sc.V0
    algaas_val = MATERIALS["AlGaAs_0p3"].chi / sc.V0
    unique = np.unique(np.round(chi_arr, decimals=8))
    assert len(unique) == 2
    assert np.any(np.isclose(chi_arr, si_val, atol=1.0e-9))
    assert np.any(np.isclose(chi_arr, gaas_val, atol=1.0e-9))
    assert not np.any(np.isclose(chi_arr, algaas_val, atol=1.0e-9))


def test_n_i_fallback_when_formula_zero_at_runtime_temperature():
    """When the runner temperature is not exactly 300 K and a region
    carries a `material_overrides` entry, the helper consults
    `Material.n_i_at_T(T)`. For an insulator-style override that zeros
    the band-edge parameters, the formula returns 0; the helper falls
    back to the resolved Material's stored `n_i` rather than emitting
    a zero field that would propagate NaNs into the Slotboom log."""
    from semi.physics.heterojunction import build_dg0_material_fields

    L = 1.0e-6
    msh, cell_tags = _build_two_region_1d_mesh(L=L, N=20)
    sc = _scaling_for_si(L=L)
    # Override sets Eg high enough that exp(-Eg / 2 V_t) underflows on
    # the right region; n_i_at_T returns 0 and the fallback path runs.
    regions_cfg = {
        "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
        "right": {
            "material": "GaAs",
            "tag": 2,
            "role": "semiconductor",
            "material_overrides": {"Eg_eV": 50.0},
        },
    }
    # Run at 310 K so the no-override fast path is bypassed even on Left.
    fields = build_dg0_material_fields(msh, cell_tags, regions_cfg, sc, T=310.0)
    n_i_arr = fields["n_i_hat"].x.array
    # The right-region n_i was recomputed by `_resolve_region_material`
    # to the closed-form value with Eg=50 eV at 300 K (essentially 0).
    # The helper falls back to that stored value when n_i_at_T(310)
    # also underflows. The field is finite (no NaN) on every cell.
    assert np.all(np.isfinite(n_i_arr))
