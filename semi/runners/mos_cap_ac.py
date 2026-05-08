"""
MOS capacitor differential C-V runner via analytic PDE sensitivity (M14.1).

Solves the same multi-region equilibrium Poisson system as
:mod:`semi.runners.mos_cv` at every gate bias, but reports the
small-signal differential capacitance C(V_gate) = dQ_gate/dV_gate by
linearising the PDE around the converged solution and solving one
additional linear system per bias point. This is the omega->0 limit of
the M14 AC sweep machinery applied to a multi-region semiconductor +
oxide stack with a gate contact, neither of which the general
:mod:`semi.runners.ac_sweep` runner supports today (it targets ohmic
two-terminal pn diodes).

Why analytic sensitivity (vs numerical dQ/dV)
---------------------------------------------
The legacy `mos_cv` runner records Q_gate(V_gate) at each bias and
the verifier extracts C_sim by `numpy.gradient(Q, V)`. That centred
finite difference is noisy for two reasons:

* It pollutes the C-V curve with O(h^2) discretisation noise from the
  bias step h and amplifies any SNES residual jitter into Q.
* The endpoints of the sweep degenerate to one-sided differences,
  losing two samples from the verifier window.

The analytic sensitivity solves K(psi_0) * delta_psi = -dF/dV at the
converged psi_0(V_gate). Because dF/dV is concentrated at the gate
Dirichlet row (only the gate BC depends on V_gate; psi_body is
V-independent) this reduces to a single linear solve with
delta_psi_gate = 1/V_t (the per-unit-V_gate BC perturbation in scaled
units) and homogeneous BCs elsewhere. C(V_gate) is then the integral

    C = (q C_0 / W_lat) integral_Si [ni_hat (exp(-psi_0) + exp(psi_0))]
                                   * delta_psi_hat dx

which is exact to discretisation (no FD noise, full sweep window).

Why a new runner instead of extending `ac_sweep` or `mos_cv`
------------------------------------------------------------
* `ac_sweep` is locked at M14 (ADR 0011, see CLAUDE.md). It assumes
  ohmic two-terminal devices in single-region (psi, n, p) primary-
  density form and would need substantial restructuring to handle the
  MOS submesh / gate-Dirichlet path.
* `mos_cv` is being deprecated by the same M14.1 task; extending it
  would block its removal.

The output is API-compatible with `mos_cv`: `iv` rows still carry
`{V, Q_gate, J: 0}` so the existing verifier and plotter keep working.
The differential capacitance is recorded as an extra `C_ac` field on
each row and exposed through `result.solver_info["C_ac"]` for
downstream consumers that want it without re-deriving via gradient.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_mos_cap_ac(cfg: dict[str, Any], *, progress_callback=None):
    """
    Run a MOS capacitor C(V_gate) sweep with analytic differential
    capacitance.

    The gate contact's `voltage_sweep` drives a sequence of V_gate
    values; at each V_gate the multi-region equilibrium Poisson system
    is solved (same as `run_mos_cv`), then the linearised system is
    solved for the per-unit-V_gate sensitivity delta_psi. iv rows
    carry `{V, Q_gate, J: 0, C_ac}` per bias point, where `C_ac` is
    in F/m^2 (2D) or F/m (1D), matching the units of `Q_gate`.
    """
    import ufl
    from dolfinx import fem
    from dolfinx.fem.petsc import LinearProblem
    from mpi4py import MPI
    from petsc4py import PETSc

    from ..bcs import build_psi_dirichlet_bcs, resolve_contacts
    from ..constants import Q
    from ..doping import build_profile
    from ..mesh import build_eps_r_function, build_mesh
    from ..physics.axisymmetric import build_equilibrium_poisson_form_axisym_mr
    from ..physics.heterojunction import (
        build_dg0_material_fields,
        cfg_uses_heterojunction,
    )
    from ..physics.poisson import build_equilibrium_poisson_form_mr
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear
    from ._common import reference_material

    coord_system = cfg.get("coordinate_system", "cartesian")
    is_axisym = coord_system == "axisymmetric"

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, cell_tags, facet_tags = build_mesh(cfg)
    if cell_tags is None:  # pragma: no cover - guarded by schema/build_mesh
        raise RuntimeError(
            "mos_cap_ac solver requires regions_by_box: the multi-region "
            "assembly cannot run without cell tags."
        )

    semi_tag = _resolve_semi_tag(cfg["regions"])
    eps_r_fn = build_eps_r_function(msh, cell_tags, cfg["regions"])

    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    N_raw_fn = build_profile(cfg["doping"])

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    psi = fem.Function(V_psi, name="psi_hat")
    two_ni = 2.0 * ref_mat.n_i
    psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))

    phys = cfg.get("physics", {})
    stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}
    # M17: build the per-cell DG0 fields when the cfg opts in. The
    # axisymmetric MOS form builder does not yet accept the
    # heterojunction_fields kwarg (axisym MOS-on-heterojunction is not
    # an existing benchmark and was not in M17 scope); the cartesian
    # MR form builder threads them through.
    het_fields = None
    if cfg_uses_heterojunction(cfg.get("regions", {})):
        het_fields = build_dg0_material_fields(
            msh, cell_tags, cfg["regions"], sc,
            T=phys.get("temperature", 300.0),
        )
    if is_axisym:
        F = build_equilibrium_poisson_form_axisym_mr(
            V_psi, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
            statistics_cfg=stat_cfg,
        )
    else:
        F = build_equilibrium_poisson_form_mr(
            V_psi, psi, N_hat_fn, sc, eps_r_fn, cell_tags, semi_tag,
            statistics_cfg=stat_cfg,
            heterojunction_fields=het_fields,
        )

    sweep_contact, sweep_values = _resolve_gate_sweep(cfg)
    if sweep_contact is None:
        raise ValueError(
            "mos_cap_ac requires a contact of type 'gate' with a voltage_sweep."
        )

    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["name"] == sweep_contact:
            continue
        if c["type"] not in ("ohmic", "gate"):
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    is_fd = stat_cfg.get("statistics", "boltzmann") == "fermi_dirac"
    if is_axisym:
        # r-weighted measures: dV_3D/(2*pi) = r dr dz on the meridian.
        r_coord = ufl.SpatialCoordinate(msh)[0]
        dx_semi = ufl.Measure(
            "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
        )
    else:
        dx_semi = ufl.Measure(
            "dx", domain=msh, subdomain_data=cell_tags, subdomain_id=int(semi_tag),
        )

    # Charge and sensitivity forms. The Boltzmann branch keeps the
    # explicit closed forms used pre-M16.4 so the assembled charge
    # integral is byte-identical to v0.19.0 on every existing
    # benchmark. The Fermi-Dirac branch builds rho_hat from the
    # generalized-Slotboom helpers and gets the sensitivity form via
    # UFL automatic differentiation (the Blakemore-corrected analogue
    # of `-d rho_hat / d psi`).
    delta_psi = fem.Function(V_psi, name="delta_psi_hat")
    if is_fd:
        from ..physics.poisson import _equilibrium_space_charge
        rho_hat_scaled = _equilibrium_space_charge(
            psi, ni_hat, N_hat_fn, sc, stat_cfg, msh,
        )
        drho_dpsi = ufl.diff(rho_hat_scaled, psi)
        if is_axisym:
            charge_form = fem.form(rho_hat_scaled * r_coord * dx_semi)
            sensitivity_form = fem.form(
                -drho_dpsi * delta_psi * r_coord * dx_semi
            )
        else:
            charge_form = fem.form(rho_hat_scaled * dx_semi)
            sensitivity_form = fem.form(
                -drho_dpsi * delta_psi * dx_semi
            )
    else:
        rho_hat_scaled = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn
        if is_axisym:
            charge_form = fem.form(rho_hat_scaled * r_coord * dx_semi)
            sensitivity_form = fem.form(
                ni_hat * (ufl.exp(-psi) + ufl.exp(psi)) * delta_psi * r_coord * dx_semi
            )
        else:
            charge_form = fem.form(rho_hat_scaled * dx_semi)
            sensitivity_form = fem.form(
                ni_hat * (ufl.exp(-psi) + ufl.exp(psi)) * delta_psi * dx_semi
            )

    # Bilinear form for the sensitivity solve: a(d, v) = d/dpsi F(psi, v)
    # at the converged psi. ufl.derivative gives this directly.
    trial_psi = ufl.TrialFunction(V_psi)
    a_form = ufl.derivative(F, psi, trial_psi)
    L_zero = fem.Constant(msh, PETSc.ScalarType(0.0)) * ufl.TestFunction(V_psi) * ufl.dx

    extents = cfg["mesh"].get("extents")
    if is_axisym:
        # In axisymmetric, integrate r ds along the gate facet to get
        # the gate's "characteristic length" used to convert the total
        # 3D-equivalent semiconductor charge into per-unit-gate-area.
        # Total gate area in 3D = 2*pi * integral_gate r ds_meridian,
        # and total semi charge in 3D = 2*pi * q * integral_Si rho r dx,
        # so Q_gate/A_gate = -q * integral_Si rho r dx / integral_gate r ds.
        gate_facet_tag = None
        for c in cfg["contacts"]:
            if c["type"] == "gate":
                gate_facet_tag = int(c["facet"])
                break
        if gate_facet_tag is None:  # pragma: no cover - schema-guarded
            raise ValueError("mos_cap_ac (axisym): no gate contact found.")
        ds_gate = ufl.Measure(
            "ds", domain=msh, subdomain_data=facet_tags, subdomain_id=gate_facet_tag,
        )
        L_gate_local = float(fem.assemble_scalar(fem.form(r_coord * ds_gate)))
        L_gate = msh.comm.allreduce(L_gate_local, op=MPI.SUM)
        if L_gate <= 0.0:  # pragma: no cover - mesh has no gate facets
            raise RuntimeError(
                f"mos_cap_ac (axisym): gate facet (tag={gate_facet_tag}) "
                f"has zero r-weighted length."
            )
        W_lat = L_gate
    else:
        if extents is None:  # pragma: no cover - cartesian path requires extents
            raise ValueError("mos_cap_ac (cartesian): cfg['mesh']['extents'] missing.")
        W_lat = float(extents[0][1] - extents[0][0])

    iv_rows: list[dict[str, float]] = []
    last_info: dict[str, Any] = {}
    psi_backup = psi.x.array.copy()

    def _solve_at(V_gate: float) -> dict[str, Any]:
        voltages = dict(static_voltages)
        voltages[sweep_contact] = float(V_gate)
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=voltages)
        bcs = build_psi_dirichlet_bcs(
            V_psi, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
            regions_cfg=cfg.get("regions"), cell_tags=cell_tags,
        )
        for bc in bcs:
            bc.set(psi.x.array)
        psi.x.scatter_forward()
        tag = _fmt_gate_tag(V_gate)
        return solve_nonlinear(
            F, psi, bcs, prefix=f"{cfg['name']}_mos_cap_ac_{tag}_",
            petsc_options={
                "snes_rtol": 1.0e-12,
                "snes_atol": 1.0e-14,
                "snes_max_it": 50,
            },
            cfg=cfg,
        )

    def _build_sensitivity_bcs(V_gate_dummy: float) -> list:
        """
        Dirichlet BCs for the sensitivity solve: delta_psi_hat = 1/V_t
        at the swept gate, 0 at every other Dirichlet boundary.

        psi_body_hat = arcsinh(N/(2 ni)) is V-independent, so its
        sensitivity is 0; for any extra ohmic / gate contact that is
        not the swept one, V is held fixed by the static_voltages map
        so the BC is also V-independent and its sensitivity is 0.

        We obtain the Dirichlet dof-set by re-using `build_psi_dirichlet_bcs`
        at a dummy V (the dofs are V-independent) and then overwriting
        the BC values: 1/V_t at the swept contact, 0 elsewhere.
        """
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages={})
        bcs: list = []
        fdim = msh.topology.dim - 1
        from petsc4py import PETSc as _PETSc
        for c in contacts:
            if c.kind not in ("ohmic", "gate"):
                continue
            facets = facet_tags.find(c.facet_tag)
            if len(facets) == 0:
                continue
            dofs = fem.locate_dofs_topological(V_psi, fdim, facets)
            value = 1.0 / sc.V0 if c.name == sweep_contact else 0.0
            bcs.append(fem.dirichletbc(_PETSc.ScalarType(value), dofs, V_psi))
        return bcs

    sensitivity_bcs = _build_sensitivity_bcs(0.0)

    sensitivity_problem = LinearProblem(
        a_form, L_zero, u=delta_psi, bcs=sensitivity_bcs,
        petsc_options_prefix=f"{cfg['name']}_mos_cap_ac_sens_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
        },
    )

    for V_gate in sweep_values:
        info = _solve_at(V_gate)
        if not info["converged"]:  # pragma: no cover - sweep is well-conditioned
            psi.x.array[:] = psi_backup
            psi.x.scatter_forward()
            raise RuntimeError(
                f"mos_cap_ac: SNES failed at V_gate={V_gate:+.4f} V "
                f"(reason={info.get('reason')})"
            )
        last_info = info
        psi_backup = psi.x.array.copy()

        rho_int_local = float(fem.assemble_scalar(charge_form))
        rho_int = msh.comm.allreduce(rho_int_local, op=MPI.SUM)
        Q_semi = Q * sc.C0 * rho_int
        Q_gate = -Q_semi / W_lat

        # Sensitivity solve.  K is rebuilt under the hood at each call
        # because `a_form` references `psi` which has just been mutated;
        # MUMPS reuses its symbolic factorisation. delta_psi is updated
        # in place with the per-unit-V_gate response.
        sensitivity_problem.solve()
        delta_psi.x.scatter_forward()

        # Numerator of C: integral_Si [ni_hat * (exp(-psi) + exp(psi)) * delta_psi] dx,
        # a dimensionless * length^dim quantity.
        sens_int_local = float(fem.assemble_scalar(sensitivity_form))
        sens_int = msh.comm.allreduce(sens_int_local, op=MPI.SUM)
        # C = q * C0 * (sens_int) / W_lat in F per unit transverse length
        # (2D mesh: F/m^2).  Sign matches dQ_gate/dV_gate = -dQ_semi/dV_gate.
        C_ac = Q * sc.C0 * sens_int / W_lat

        iv_rows.append({
            "V": float(V_gate),
            "Q_gate": float(Q_gate),
            "J": 0.0,
            "C_ac": float(C_ac),
        })
        if progress_callback is not None:
            progress_callback({
                "type": "step_done",
                "bias_step": len(iv_rows) - 1,
                "V_applied": float(V_gate),
                "iterations": int(info.get("iterations", 0)),
            })

    last_info = dict(last_info)
    last_info["C_ac"] = [float(r["C_ac"]) for r in iv_rows]

    x_dof = V_psi.tabulate_dof_coordinates()
    psi_phys = psi.x.array * sc.V0
    n_phys = ref_mat.n_i * np.exp(psi.x.array)
    p_phys = ref_mat.n_i * np.exp(-psi.x.array)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V_psi, psi=psi,
        psi_phys=psi_phys, n_phys=n_phys, p_phys=p_phys,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc,
        solver_info=last_info, iv=iv_rows, bias_contact=sweep_contact,
    )


def _resolve_semi_tag(regions_cfg: dict) -> int:
    for r in regions_cfg.values():
        if "tag" not in r:  # pragma: no cover - schema requires `tag`
            continue
        if r.get("role", "semiconductor") == "semiconductor":
            return int(r["tag"])
    raise ValueError("mos_cap_ac: no region with role='semiconductor'.")  # pragma: no cover


def _resolve_gate_sweep(cfg) -> tuple[str | None, list[float]]:
    for c in cfg["contacts"]:
        if c["type"] != "gate":
            continue
        sweep = c.get("voltage_sweep")
        if sweep is None:  # pragma: no cover - guarded above
            continue
        start = float(sweep["start"])
        stop = float(sweep["stop"])
        step = float(sweep["step"])
        if step <= 0.0:  # pragma: no cover - schema enforces step > 0
            raise ValueError("voltage_sweep.step must be positive")
        direction = 1.0 if stop >= start else -1.0
        n = int(round((stop - start) / (direction * step))) + 1
        values = [start + i * direction * step for i in range(n)]
        return c["name"], values
    return None, []


def _fmt_gate_tag(v: float) -> str:
    return f"{v:+.4f}".replace("+", "p").replace("-", "m").replace(".", "d")
