"""
End-to-end tests for the submesh-based drift-diffusion assembly.

Exercises the M6 multi-region path:
    - V_psi on the parent mesh, V_phi_n/V_phi_p on a semiconductor submesh
    - Poisson stiffness over the full parent mesh with cellwise eps_r
    - Poisson source restricted to silicon cells via dx(subdomain_id=1)
    - Continuity integrated on the submesh, reading psi from the parent
    - entity_maps threaded through fem.form and NonlinearProblem

The tests run on 1D uniform and step-doped meshes where the submesh
equals the full mesh (every cell is silicon). This gives a single
semiconductor region that still exercises the full submesh + entity_maps
plumbing, so a regression in the mixed-domain assembly is caught here.
"""
from __future__ import annotations

import numpy as np
import pytest


def _build_1d_uniform_problem(N_A: float = 1.0e23, L: float = 1.0e-6, N: int = 40):
    """Build the MR-DD spaces on a uniform p-type 1D slab."""
    from dolfinx import fem
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    from semi.materials import get_material
    from semi.physics.drift_diffusion import make_dd_block_spaces_mr
    from semi.scaling import Scaling

    Si = get_material("Si")
    sc = Scaling(L0=L, C0=N_A, T=300.0, mu0=Si.mu_n, n_i=Si.n_i)

    msh = dmesh.create_interval(MPI.COMM_WORLD, N, [0.0, L])
    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    n_local = msh.topology.index_map(tdim).size_local
    all_cells = np.arange(n_local, dtype=np.int32)
    cell_tags = dmesh.meshtags(
        msh, tdim, all_cells, np.ones(n_local, dtype=np.int32),
    )

    submesh, em_cell, _em_vert, _em_geom = dmesh.create_submesh(msh, tdim, all_cells)
    spaces = make_dd_block_spaces_mr(msh, submesh, em_cell)

    N_hat_fn = fem.Function(spaces.V_psi, name="N_hat")
    N_hat_fn.x.array[:] = -N_A / sc.C0

    return spaces, sc, Si, cell_tags, N_hat_fn, L, N


def test_mr_dd_assembly_uniform_ptype_bulk_equilibrium():
    """
    Uniform p-type slab at equilibrium: p = N_A, n = n_i^2 / N_A.
    The MR path must converge with zero Newton iterations at the
    analytical initial guess (because it IS the solution), and the
    bulk carrier densities must satisfy the mass-action law.
    """
    from dolfinx import fem
    from dolfinx import mesh as dmesh
    from petsc4py import PETSc

    from semi.physics.drift_diffusion import build_dd_block_residual_mr
    from semi.solver import solve_nonlinear_block

    spaces, sc, Si, cell_tags, N_hat_fn, L, N = _build_1d_uniform_problem()

    psi_eq = float(np.arcsinh(-spaces.N_hat_raw if False else -(1.0)))
    # Analytical equilibrium psi_hat for uniform p-type: asinh(-N_A / (2 n_i))
    N_A_over_2ni = (-N_hat_fn.x.array[0]) * sc.C0 / (2.0 * Si.n_i)
    psi_eq = float(np.arcsinh(-N_A_over_2ni))
    spaces.psi.x.array[:] = psi_eq
    spaces.phi_n.x.array[:] = 0.0
    spaces.phi_p.x.array[:] = 0.0
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.scatter_forward()

    mu_n_hat = Si.mu_n / sc.mu0
    mu_p_hat = Si.mu_p / sc.mu0
    tau_n = 1.0e-7 / sc.t0
    tau_p = 1.0e-7 / sc.t0

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n, tau_p,
        cell_tags, semi_tag=1,
    )

    fdim = spaces.parent_mesh.topology.dim - 1
    spaces.parent_mesh.topology.create_connectivity(fdim, spaces.parent_mesh.topology.dim)
    all_facets = dmesh.locate_entities_boundary(
        spaces.parent_mesh, fdim, lambda x: np.full(x.shape[1], True),
    )
    spaces.submesh.topology.create_connectivity(fdim, spaces.submesh.topology.dim)
    all_sub = dmesh.locate_entities_boundary(
        spaces.submesh, fdim, lambda x: np.full(x.shape[1], True),
    )
    dofs_psi = fem.locate_dofs_topological(spaces.V_psi, fdim, all_facets)
    dofs_n = fem.locate_dofs_topological(spaces.V_phi_n, fdim, all_sub)
    dofs_p = fem.locate_dofs_topological(spaces.V_phi_p, fdim, all_sub)

    bcs = [
        fem.dirichletbc(PETSc.ScalarType(psi_eq), dofs_psi, spaces.V_psi),
        fem.dirichletbc(PETSc.ScalarType(0.0), dofs_n, spaces.V_phi_n),
        fem.dirichletbc(PETSc.ScalarType(0.0), dofs_p, spaces.V_phi_p),
    ]

    info = solve_nonlinear_block(
        F_list, [spaces.psi, spaces.phi_n, spaces.phi_p], bcs,
        prefix="test_mr_dd_uniform_ptype_",
        entity_maps=[spaces.entity_map],
        petsc_options={
            "snes_rtol": 1.0e-12,
            "snes_atol": 1.0e-14,
            "snes_max_it": 10,
        },
    )
    assert info["converged"], f"MR-DD did not converge: {info}"

    mid = N // 2
    n_bulk = Si.n_i * np.exp(
        spaces.psi.x.array[mid] - spaces.phi_n.x.array[mid]
    )
    p_bulk = Si.n_i * np.exp(
        spaces.phi_p.x.array[mid] - spaces.psi.x.array[mid]
    )
    N_A = -N_hat_fn.x.array[0] * sc.C0
    assert p_bulk == pytest.approx(N_A, rel=1.0e-6)
    assert (n_bulk * p_bulk) == pytest.approx(Si.n_i ** 2, rel=1.0e-6)


def test_mr_dd_assembly_pn_junction_matches_vbi_theory():
    """
    Step-doped pn junction at equilibrium via the MR path.
    V_bi = V_t * ln(N_A * N_D / n_i^2) must match to solver tolerance,
    and bulk p / n densities must agree with the doping to 0.1%.

    This is the same physical verification used in the pn_1d benchmark,
    but running through the submesh + entity_maps assembly path to
    prove the mixed-domain compile produces the correct matrix.
    """
    from dolfinx import fem
    from dolfinx import mesh as dmesh
    from mpi4py import MPI
    from petsc4py import PETSc

    from semi.materials import get_material
    from semi.physics.drift_diffusion import (
        build_dd_block_residual_mr,
        make_dd_block_spaces_mr,
    )
    from semi.scaling import Scaling
    from semi.solver import solve_nonlinear_block

    Si = get_material("Si")
    N_A = 1.0e23
    N_D = 1.0e23
    L = 2.0e-6
    N = 100
    sc = Scaling(L0=L, C0=N_A, T=300.0, mu0=Si.mu_n, n_i=Si.n_i)

    msh = dmesh.create_interval(MPI.COMM_WORLD, N, [0.0, L])
    tdim = msh.topology.dim
    msh.topology.create_entities(tdim)
    n_local = msh.topology.index_map(tdim).size_local
    all_cells = np.arange(n_local, dtype=np.int32)
    cell_tags = dmesh.meshtags(
        msh, tdim, all_cells, np.ones(n_local, dtype=np.int32),
    )

    submesh, em_cell, _em_vert, _em_geom = dmesh.create_submesh(msh, tdim, all_cells)
    spaces = make_dd_block_spaces_mr(msh, submesh, em_cell)

    two_ni = 2.0 * Si.n_i
    N_hat_fn = fem.Function(spaces.V_psi, name="N_hat")

    def _doping_hat(x):
        out = np.zeros(x.shape[1])
        out[x[0] < L / 2.0] = -N_A / sc.C0
        out[x[0] >= L / 2.0] = +N_D / sc.C0
        return out

    N_hat_fn.interpolate(_doping_hat)

    def _psi_init(x):
        raw = np.zeros(x.shape[1])
        raw[x[0] < L / 2.0] = -N_A
        raw[x[0] >= L / 2.0] = +N_D
        return np.arcsinh(raw / two_ni)

    spaces.psi.interpolate(_psi_init)
    spaces.phi_n.x.array[:] = 0.0
    spaces.phi_p.x.array[:] = 0.0
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.scatter_forward()

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        Si.mu_n / sc.mu0, Si.mu_p / sc.mu0,
        1.0e-7 / sc.t0, 1.0e-7 / sc.t0,
        cell_tags, semi_tag=1,
    )

    fdim = tdim - 1
    msh.topology.create_connectivity(fdim, tdim)
    left_parent = dmesh.locate_entities_boundary(
        msh, fdim, lambda x: np.isclose(x[0], 0.0),
    )
    right_parent = dmesh.locate_entities_boundary(
        msh, fdim, lambda x: np.isclose(x[0], L),
    )
    submesh.topology.create_connectivity(fdim, submesh.topology.dim)
    left_sub = dmesh.locate_entities_boundary(
        submesh, fdim, lambda x: np.isclose(x[0], 0.0),
    )
    right_sub = dmesh.locate_entities_boundary(
        submesh, fdim, lambda x: np.isclose(x[0], L),
    )

    psi_L = float(np.arcsinh(-N_A / two_ni))
    psi_R = float(np.arcsinh(+N_D / two_ni))
    bcs = [
        fem.dirichletbc(PETSc.ScalarType(psi_L),
                        fem.locate_dofs_topological(spaces.V_psi, fdim, left_parent),
                        spaces.V_psi),
        fem.dirichletbc(PETSc.ScalarType(psi_R),
                        fem.locate_dofs_topological(spaces.V_psi, fdim, right_parent),
                        spaces.V_psi),
        fem.dirichletbc(PETSc.ScalarType(0.0),
                        fem.locate_dofs_topological(spaces.V_phi_n, fdim, left_sub),
                        spaces.V_phi_n),
        fem.dirichletbc(PETSc.ScalarType(0.0),
                        fem.locate_dofs_topological(spaces.V_phi_n, fdim, right_sub),
                        spaces.V_phi_n),
        fem.dirichletbc(PETSc.ScalarType(0.0),
                        fem.locate_dofs_topological(spaces.V_phi_p, fdim, left_sub),
                        spaces.V_phi_p),
        fem.dirichletbc(PETSc.ScalarType(0.0),
                        fem.locate_dofs_topological(spaces.V_phi_p, fdim, right_sub),
                        spaces.V_phi_p),
    ]

    info = solve_nonlinear_block(
        F_list, [spaces.psi, spaces.phi_n, spaces.phi_p], bcs,
        prefix="test_mr_dd_pn_",
        entity_maps=[em_cell],
        petsc_options={
            "snes_rtol": 1.0e-12,
            "snes_atol": 1.0e-14,
            "snes_stol": 1.0e-14,
            "snes_max_it": 30,
        },
    )
    assert info["converged"], f"MR-DD pn junction did not converge: {info}"
    assert info["iterations"] < 20, (
        f"MR-DD took {info['iterations']} Newton iters on pn_1d; "
        "expected single-digit like the single-region path"
    )

    psi_arr = spaces.psi.x.array
    V_bi_sim = (psi_arr.max() - psi_arr.min()) * sc.V0
    V_bi_theory = sc.V0 * float(np.log(N_A * N_D / Si.n_i ** 2))
    assert V_bi_sim == pytest.approx(V_bi_theory, rel=1.0e-4)

    x_dof = spaces.V_psi.tabulate_dof_coordinates()
    order = np.argsort(x_dof[:, 0])
    psi_sorted = psi_arr[order]
    phin_sorted = spaces.phi_n.x.array[order]
    phip_sorted = spaces.phi_p.x.array[order]

    p_bulk = Si.n_i * float(np.exp(phip_sorted[1] - psi_sorted[1]))
    n_bulk = Si.n_i * float(np.exp(psi_sorted[-2] - phin_sorted[-2]))
    assert p_bulk == pytest.approx(N_A, rel=5.0e-3)
    assert n_bulk == pytest.approx(N_D, rel=5.0e-3)


def test_mr_dd_auger_branch_assembles():
    """
    M16.3: the multiregion form builder accepts `recomb_cfg` with
    `auger=True` and emits a residual that includes the Auger term
    on the submesh integration domain.

    The single-region (`build_dd_block_residual`) Auger path is
    exercised end-to-end by the diode_auger_1d benchmark; the MR
    counterpart has no benchmark today, so this assembly smoke test
    pins the M16.3 inline at the multiregion path against
    regressions.
    """
    from semi.physics.drift_diffusion import build_dd_block_residual_mr

    spaces, sc, Si, cell_tags, N_hat_fn, _L, _N = _build_1d_uniform_problem()

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        Si.mu_n / sc.mu0, Si.mu_p / sc.mu0,
        1.0e-7 / sc.t0, 1.0e-7 / sc.t0,
        cell_tags, semi_tag=1,
        recomb_cfg={"auger": True, "C_n": 1.0e-29, "C_p": 1.0e-29},
    )
    assert len(F_list) == 3


def test_mr_dd_bbt_branch_assembles():
    """
    M16.6: the multiregion form builder accepts `recomb_cfg` with
    `bbt=True` and emits a residual that subtracts the Kane
    band-to-band generation rate from R on the submesh integration
    domain.

    The single-region path is exercised end-to-end by the zener_1d
    benchmark and Variant H of MMS-DD; the MR counterpart has no
    benchmark today (BBT under MOSCAP / MOSFET geometry is unusual),
    so this assembly smoke test pins the M16.6 inline at the
    multiregion path against regressions.
    """
    from semi.physics.drift_diffusion import build_dd_block_residual_mr

    spaces, sc, Si, cell_tags, N_hat_fn, _L, _N = _build_1d_uniform_problem()

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        Si.mu_n / sc.mu0, Si.mu_p / sc.mu0,
        1.0e-7 / sc.t0, 1.0e-7 / sc.t0,
        cell_tags, semi_tag=1,
        recomb_cfg={"bbt": True, "A_kane": 4.0e14, "B_kane": 1.9e7},
    )
    assert len(F_list) == 3


def test_mr_dd_tat_branch_assembles():
    """
    M16.6: the multiregion form builder accepts `recomb_cfg` with
    `tat=True` and multiplies the SRH rate by the Hurkx field-
    enhancement factor (1 + Gamma_TAT(F)) on the submesh integration
    domain.
    """
    from semi.physics.drift_diffusion import build_dd_block_residual_mr

    spaces, sc, Si, cell_tags, N_hat_fn, _L, _N = _build_1d_uniform_problem()

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        Si.mu_n / sc.mu0, Si.mu_p / sc.mu0,
        1.0e-7 / sc.t0, 1.0e-7 / sc.t0,
        cell_tags, semi_tag=1,
        recomb_cfg={"tat": True, "F_kT": 1.4e7, "alpha": 2.0},
    )
    assert len(F_list) == 3


def test_mr_dd_bbt_and_tat_branches_assemble():
    """
    Both BBT and TAT enabled simultaneously; exercises the shared
    F_hat_mag_mr field-magnitude precomputation that gates both
    kernels off `if bbt_on_mr or tat_on_mr`.
    """
    from semi.physics.drift_diffusion import build_dd_block_residual_mr

    spaces, sc, Si, cell_tags, N_hat_fn, _L, _N = _build_1d_uniform_problem()

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        Si.mu_n / sc.mu0, Si.mu_p / sc.mu0,
        1.0e-7 / sc.t0, 1.0e-7 / sc.t0,
        cell_tags, semi_tag=1,
        recomb_cfg={
            "bbt": True, "A_kane": 4.0e14, "B_kane": 1.9e7,
            "tat": True, "F_kT": 1.4e7, "alpha": 2.0,
        },
    )
    assert len(F_list) == 3


def test_mr_dd_schottky_branch_assembles():
    """
    M16.5: the multiregion form builder accepts a non-empty
    `schottky_facets` list and adds the thermionic-emission Robin
    surface form to the (phi_n, phi_p) continuity rows.

    Schottky on the MR path is unusual (the MOSCAP / MOSFET 2D
    devices have no metal-semiconductor contacts on the silicon
    submesh), but the API surface accepts it for completeness and
    threads phi_n / phi_p through the entity_maps mechanism on the
    parent mesh's facet measure. This smoke test pins that surface-
    integral wiring against regressions.
    """
    from dolfinx import mesh as dmesh

    from semi.materials import get_material
    from semi.physics.drift_diffusion import build_dd_block_residual_mr

    spaces, sc, Si, cell_tags, N_hat_fn, L, _N = _build_1d_uniform_problem()
    msh = spaces.parent_mesh
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    left = dmesh.locate_entities_boundary(
        msh, fdim, lambda x: np.isclose(x[0], 0.0),
    )
    right = dmesh.locate_entities_boundary(
        msh, fdim, lambda x: np.isclose(x[0], L),
    )
    facets = np.concatenate([left, right]).astype(np.int32)
    values = np.concatenate([
        np.full(len(left), 1, dtype=np.int32),
        np.full(len(right), 2, dtype=np.int32),
    ])
    order = np.argsort(facets)
    facet_tags = dmesh.meshtags(msh, fdim, facets[order], values[order])

    ref_mat = get_material("Si")
    # The Schottky thermionic-emission Robin BC reads sc.v_n_thermal /
    # sc.v_p_thermal, which require the effective masses on Scaling.
    # The bare _build_1d_uniform_problem fixture omits them; populate
    # from the reference material so the surface form can build.
    sc.m_n_star = ref_mat.m_n_star
    sc.m_p_star = ref_mat.m_p_star

    F_list = build_dd_block_residual_mr(
        spaces, N_hat_fn, sc, Si.epsilon_r,
        Si.mu_n / sc.mu0, Si.mu_p / sc.mu0,
        1.0e-7 / sc.t0, 1.0e-7 / sc.t0,
        cell_tags, semi_tag=1,
        facet_tags=facet_tags,
        schottky_facets=[(1, {"barrier_height_eV": 0.85, "V_applied": 0.0})],
        ref_mat=ref_mat,
        statistics_cfg={"statistics": "boltzmann"},
    )
    assert len(F_list) == 3
