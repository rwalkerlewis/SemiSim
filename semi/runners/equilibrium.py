"""
Equilibrium Poisson runner (M1).

Solves -div(L_D^2 eps_r grad psi_hat) = p_hat - n_hat + N_hat under
Boltzmann statistics on the validated config and returns a populated
`SimulationResult`. No applied bias, no continuity equations.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_equilibrium(cfg: dict[str, Any], *, progress_callback=None):
    """M1 equilibrium Poisson solver (Boltzmann statistics)."""
    from dolfinx import fem

    from ..bcs import build_psi_dirichlet_bcs, resolve_contacts
    from ..doping import build_profile
    from ..mesh import build_mesh
    from ..physics.heterojunction import (
        build_dg0_material_fields,
        cfg_uses_heterojunction,
    )
    from ..physics.poisson import build_equilibrium_poisson_form
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear
    from ._common import reference_material

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, cell_tags, facet_tags = build_mesh(cfg)
    V = fem.functionspace(msh, ("Lagrange", 1))

    N_raw_fn = build_profile(cfg["doping"])

    def N_hat_expr(x):
        return N_raw_fn(x) / sc.C0

    N_hat_fn = fem.Function(V, name="N_net_hat")
    N_hat_fn.interpolate(N_hat_expr)

    psi = fem.Function(V, name="psi_hat")
    contacts = resolve_contacts(cfg, facet_tags=facet_tags)
    bcs = build_psi_dirichlet_bcs(
        V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
        regions_cfg=cfg.get("regions"), cell_tags=cell_tags,
    )

    two_ni = 2.0 * ref_mat.n_i

    def psi_init_expr(x):
        return np.arcsinh(N_raw_fn(x) / two_ni)

    psi.interpolate(psi_init_expr)
    for bc in bcs:
        bc.set(psi.x.array)
    psi.x.scatter_forward()

    phys = cfg.get("physics", {})
    stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}
    # M17: build the per-cell DG0 chi/Eg/Nc/Nv/n_i/eps_r fields when the
    # cfg opts into the heterojunction path; otherwise pass None and the
    # form builder runs the v0.23.0 scalar `ni_hat` Constant path.
    het_fields = None
    if cfg_uses_heterojunction(cfg.get("regions", {})):
        het_fields = build_dg0_material_fields(
            msh, cell_tags, cfg["regions"], sc, T=phys.get("temperature", 300.0),
        )
    F = build_equilibrium_poisson_form(
        V, psi, N_hat_fn, sc, ref_mat.epsilon_r,
        statistics_cfg=stat_cfg,
        heterojunction_fields=het_fields,
    )
    if progress_callback is not None:
        progress_callback({"type": "step_started", "bias_step": 0, "V_applied": 0.0})
    info = solve_nonlinear(F, psi, bcs, prefix=f"{cfg['name']}_", cfg=cfg)
    if progress_callback is not None:
        progress_callback({
            "type": "step_done", "bias_step": 0, "V_applied": 0.0,
            "iterations": int(info.get("iterations", 0)),
        })

    x_dof = V.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    if stat_cfg.get("statistics", "boltzmann") == "fermi_dirac":
        # Equilibrium FD: n = n_i * gamma_n * exp(psi), with phi_n = 0.
        from ..physics.statistics import gamma_n_blakemore, gamma_p_blakemore
        gamma_n = gamma_n_blakemore(psi.x.array, sc.eta_offset_n)
        gamma_p = gamma_p_blakemore(-psi.x.array, sc.eta_offset_p)
        n_phys = ref_mat.n_i * gamma_n * np.exp(psi.x.array)
        p_phys = ref_mat.n_i * gamma_p * np.exp(-psi.x.array)
    else:
        n_phys = ref_mat.n_i * np.exp(psi.x.array)
        p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc, solver_info=info,
    )
