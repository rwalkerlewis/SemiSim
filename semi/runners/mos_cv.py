"""
MOS capacitor C-V runner (M6).

.. deprecated:: M14.1
    Prefer :mod:`semi.runners.mos_cap_ac` (``solver.type == "mos_cap_ac"``)
    for new MOS C-V workloads. That runner solves the same multi-region
    equilibrium Poisson system but reports the differential capacitance
    C(V_gate) via analytic PDE sensitivity, which is exact to discretisation
    rather than approximated through `numpy.gradient` of Q_gate(V_gate).
    `mos_cv` is kept for backwards compatibility with consumers that already
    rely on the Q-then-finite-difference path; removal is a future cleanup.

Solves multi-region equilibrium Poisson at each gate bias in a sweep,
integrates the semiconductor space charge to recover Q_gate(V_gate),
and returns a `SimulationResult` whose `iv` rows carry `(V, Q_gate)`
pairs for the C-V verifier.

Why equilibrium Poisson (not block drift-diffusion). The gate contact
passes no current; with phi_n = phi_p = 0 pinned at the body, carriers
in the silicon are Boltzmann in psi and the full (psi, phi_n, phi_p)
block system collapses to the single-equation equilibrium Poisson. The
multi-region plumbing (cellwise eps_r on the parent mesh, space-charge
restricted to silicon cells via `dx(subdomain_id=semi_tag)`) still
applies; the oxide region carries only the Laplacian.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_mos_cv(cfg: dict[str, Any], *, progress_callback=None):
    """
    Run a MOS capacitor C-V sweep.

    The gate contact's `voltage_sweep` drives a sequence of V_gate
    values; at each V_gate the multi-region equilibrium Poisson system
    is solved and the gate charge per unit area is recorded in the
    SimulationResult's `iv` list as `{V, Q_gate, J=0}`. The body
    ohmic contact uses the standard equilibrium psi Dirichlet
    (asinh(N_net / 2 n_i)).
    """
    from dolfinx import fem
    from mpi4py import MPI

    from ..bcs import build_psi_dirichlet_bcs, resolve_contacts
    from ..constants import Q
    from ..doping import build_profile
    from ..mesh import build_eps_r_function, build_mesh
    from ..physics.heterojunction import (
        build_dg0_material_fields,
        cfg_uses_heterojunction,
    )
    from ..physics.poisson import build_equilibrium_poisson_form_mr
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear
    from ._common import reference_material

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, cell_tags, facet_tags = build_mesh(cfg)
    if cell_tags is None:
        raise RuntimeError(
            "mos_cv solver requires regions_by_box: the multi-region "
            "assembly cannot run without cell tags."
        )

    semi_tag = _resolve_semi_tag(cfg["regions"])
    eps_r_fn = build_eps_r_function(msh, cell_tags, cfg["regions"])

    V = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    N_hat_fn = fem.Function(V, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    psi = fem.Function(V, name="psi_hat")
    two_ni = 2.0 * ref_mat.n_i
    psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))

    phys = cfg.get("physics", {})
    stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}
    # M17: thread the heterojunction DG0 fields through when the cfg
    # opts in. MOSCAP runs on Si + SiO2 which is not a semiconductor
    # heterojunction; the threading is here so that future MOS-on-
    # heterojunction work has the seam in place.
    het_fields = None
    if cfg_uses_heterojunction(cfg.get("regions", {})):
        het_fields = build_dg0_material_fields(
            msh, cell_tags, cfg["regions"], sc,
            T=phys.get("temperature", 300.0),
        )
    F = build_equilibrium_poisson_form_mr(
        V, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
        statistics_cfg=stat_cfg,
        heterojunction_fields=het_fields,
    )

    sweep_contact, sweep_values = _resolve_gate_sweep(cfg)
    if sweep_contact is None:
        raise ValueError(
            "mos_cv requires a contact of type 'gate' with a voltage_sweep."
        )

    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["name"] == sweep_contact:
            continue
        if c["type"] not in ("ohmic", "gate"):
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    # Per-cycle charge form (rebuilt each iteration is overkill; rho_hat
    # is built off `psi` which is mutated in place, so one form suffices).
    # Under FD the same Blakemore prefactor that the residual sees must
    # be applied here, otherwise the gate-charge integral would not be
    # consistent with the solved psi profile.
    import ufl
    from petsc4py import PETSc

    from ..physics.poisson import _equilibrium_space_charge
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    dx_semi = ufl.Measure(
        "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
    )
    rho_hat_scaled = _equilibrium_space_charge(
        psi, ni_hat, N_hat_fn, sc, stat_cfg, msh,
    )
    charge_form = fem.form(rho_hat_scaled * dx_semi)

    extents = cfg["mesh"]["extents"]
    W_lat = float(extents[0][1] - extents[0][0])  # lateral extent, m

    iv_rows: list[dict[str, float]] = []
    last_info: dict[str, Any] = {}

    # Snapshot / restore lets us recover if a step fails, but the MOS
    # sweep is well-conditioned and we expect every step to converge.
    psi_backup = psi.x.array.copy()

    def solve_at(V_gate: float) -> dict[str, Any]:
        voltages = dict(static_voltages)
        voltages[sweep_contact] = float(V_gate)
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=voltages)
        bcs = build_psi_dirichlet_bcs(
            V, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
            regions_cfg=cfg.get("regions"), cell_tags=cell_tags,
        )
        for bc in bcs:
            bc.set(psi.x.array)
        psi.x.scatter_forward()
        tag = _fmt_gate_tag(V_gate)
        return solve_nonlinear(
            F, psi, bcs, prefix=f"{cfg['name']}_mos_cv_{tag}_",
            petsc_options={
                "snes_rtol": 1.0e-12,
                "snes_atol": 1.0e-14,
                "snes_max_it": 50,
            },
            cfg=cfg,
        )

    for V_gate in sweep_values:
        info = solve_at(V_gate)
        if not info["converged"]:
            psi.x.array[:] = psi_backup
            psi.x.scatter_forward()
            raise RuntimeError(
                f"mos_cv: SNES failed at V_gate={V_gate:+.4f} V "
                f"(reason={info.get('reason')})"
            )
        last_info = info
        psi_backup = psi.x.array.copy()

        rho_int_local = float(fem.assemble_scalar(charge_form))
        rho_int = msh.comm.allreduce(rho_int_local, op=MPI.SUM)
        # rho_int = integral of rho_hat [dimensionless] * dx [m^dim].
        # Actual semiconductor charge per unit depth (2D: C/m) is
        #     Q_semi_line = q * C_0 * rho_int
        # Per unit gate area (divide by lateral extent):
        #     Q_semi_area = Q_semi_line / W_lat        [C/m^2]
        # By charge conservation, Q_gate = -Q_semi (with the sign
        # convention that positive V_gate yields positive Q_gate).
        Q_semi = Q * sc.C0 * rho_int
        Q_gate = -Q_semi / W_lat

        iv_rows.append({
            "V": float(V_gate),
            "Q_gate": float(Q_gate),
            "J": 0.0,
        })
        if progress_callback is not None:
            progress_callback({
                "type": "step_done",
                "bias_step": len(iv_rows) - 1,
                "V_applied": float(V_gate),
                "iterations": int(info.get("iterations", 0)),
            })

    x_dof = V.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    n_phys = ref_mat.n_i * np.exp(psi.x.array)
    p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc,
        solver_info=last_info, iv=iv_rows, bias_contact=sweep_contact,
    )


def _resolve_semi_tag(regions_cfg: dict) -> int:
    """First region with role='semiconductor' wins. Raises if none."""
    for r in regions_cfg.values():
        if "tag" not in r:
            continue
        if r.get("role", "semiconductor") == "semiconductor":
            return int(r["tag"])
    raise ValueError("mos_cv: no region with role='semiconductor'.")


def _resolve_gate_sweep(cfg) -> tuple[str | None, list[float]]:
    """Find the gate contact's voltage_sweep and expand it to values."""
    for c in cfg["contacts"]:
        if c["type"] != "gate":
            continue
        sweep = c.get("voltage_sweep")
        if sweep is None:
            continue
        start = float(sweep["start"])
        stop = float(sweep["stop"])
        step = float(sweep["step"])
        if step <= 0.0:
            raise ValueError("voltage_sweep.step must be positive")
        direction = 1.0 if stop >= start else -1.0
        n = int(round((stop - start) / (direction * step))) + 1
        values = [start + i * direction * step for i in range(n)]
        return c["name"], values
    return None, []


def _fmt_gate_tag(v: float) -> str:
    """PETSc-prefix-safe encoding of a gate voltage."""
    return f"{v:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
