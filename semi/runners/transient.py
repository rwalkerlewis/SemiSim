"""
Transient drift-diffusion runner (BDF1 / BDF2) — Slotboom form.

Solves the coupled `(psi, phi_n, phi_p)` system in **Slotboom form**
with backward differentiation formula (BDF) time integration. Carrier
densities are recovered pointwise via the Slotboom relations:

    n = n_i exp(psi - phi_n)        [scaled, V_t = 1]
    p = n_i exp(phi_p - psi)

Because `n` and `p` are exponentials of the primary unknowns they are
strictly positive at every Newton iterate, eliminating the negative-
density failure modes of the (psi, n, p) primary-density formulation
documented in `docs/m13.1-followup-5-blocker.md`.

See `docs/adr/0014-slotboom-transient.md`. ADR 0014 supersedes ADR 0009.

Time integration
----------------
BDF order is configurable (1 = backward Euler, 2 = BDF2). For BDF2, a
single BDF1 step seeds the two-level history before switching. Fixed
timestep `dt` is used throughout (no adaptive time stepping in M13).
See `docs/adr/0010-bdf-time-integration.md`.

The time term applies BDF directly to the carrier *densities*
`n_ufl = n_i exp(psi - phi_n)` and `p_ufl = n_i exp(phi_p - psi)`.
UFL's automatic differentiation of `(alpha_0/dt) * n_ufl * v_n` against
`(psi, phi_n)` produces the chain-rule mass matrix

    M_n_diag(phi_n, phi_n) = -(alpha_0/dt) * n
    M_n_cross(psi, phi_n)  = +(alpha_0/dt) * n        [coupling to Poisson]

per ADR 0014 § Implementation. No hand-written Jacobian is needed.

Residual form (scaled units)
-----------------------------
Poisson (same as bias_sweep):
    L_D^2 eps_r grad(psi) . grad(v_psi)  -  (p - n + N) v_psi  =  0

Electron continuity (Slotboom):
    (alpha_0/dt) n v_n + (f_hist_n / dt) v_n
    +  L_0^2 mu_n n grad(phi_n) . grad(v_n)
    -  R v_n  =  0

Hole continuity (Slotboom):
    (alpha_0/dt) p v_p + (f_hist_p / dt) v_p
    +  L_0^2 mu_p p grad(phi_p) . grad(v_p)
    +  R v_p  =  0

with `f_hist_n`, `f_hist_p` per-DOF Functions storing
`sum_{k=1}^K alpha_k * n^{n+1-k}` and likewise for p, evaluated at
past converged Slotboom states.

Boundary conditions (ohmic contacts)
--------------------------------------
Reuses `build_dd_dirichlet_bcs` from `semi.bcs` exactly as
`run_bias_sweep` does. At an ohmic contact under applied bias V
(scaled `V_hat = V / V_t`):

    psi_bc   = arcsinh(N_net / (2 n_i)) + V_hat
    phi_n_bc = phi_p_bc = V_hat                         [Shockley boundary]

Time-varying contact voltage (M16.7)
------------------------------------
A contact may carry an optional `voltage_t` block instead of a fixed
`voltage` value. Two variants are supported:

  - `{"type": "table", "times": [...], "values": [...]}`: linear
    interpolation between `(times[i], values[i])` pairs; `t < times[0]`
    clamps to `values[0]` and `t > times[-1]` clamps to `values[-1]`.
  - `{"type": "step", "t0": ..., "v0": ..., "v1": ...}`: piecewise
    constant; returns `v0` for `t < t0` and `v1` for `t >= t0`.

The BC-ramp continuation (`bc_ramp_steps`) does not consume
`voltage_t`. It ramps to `voltage_t.values[0]` (table) or
`voltage_t.v0` (step). After the ramp, the time loop builds BCs from
`voltages_at_t(t_next)` at every step. Configs without `voltage_t`
on any contact are bit-identical to v0.22.0; the evaluator returns
the same dict at every t and the BC stack is built from the fixed
`voltage` values just as before.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from ..physics.drift_diffusion import DDBlockSpaces
    from ..results import TransientResult


def run_transient(
    cfg: dict[str, Any],
    *,
    progress_callback=None,
) -> TransientResult:
    """
    Coupled transient drift-diffusion solver (Slotboom + chain rule).

    Runs BDF1 or BDF2 time integration on the (psi, phi_n, phi_p)
    Slotboom block system, starting from a converged steady-state IC at
    `V_target` (BC-ramp continuation, ADR 0013) and stepping in time
    with the contact voltages held fixed.

    Parameters
    ----------
    cfg : dict
        Validated JSON config dict. ``solver.type`` must be
        ``"transient"``. Required solver sub-keys:

        - ``t_end`` (float): simulation end time in seconds.
        - ``dt`` (float): fixed time step in seconds.
        - ``order`` (int, default 2): BDF order (1 or 2).
        - ``max_steps`` (int, default 10000): safety cap on iterations.
        - ``bc_ramp_steps`` (int, default 10): number of steady-state
          sub-steps to ramp the bias from V=0 to its target before the
          time loop. ``0`` disables continuation (used by
          ``pn_1d_turnon``: the V=0 IC + step-bias-at-t=0 *is* the
          physical scenario being measured).

    progress_callback : callable or None
        If provided, called after each successful timestep with a dict::

            {"type": "step_done", "step": n, "t": t_current,
             "iterations": n_snes_iters}

    Returns
    -------
    TransientResult
    """
    import copy

    from dolfinx import fem

    from ..bcs import build_dd_dirichlet_bcs, resolve_contacts
    from ..doping import build_profile
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import make_dd_block_spaces
    from ..physics.heterojunction import (
        build_dg0_material_fields,
        cfg_uses_heterojunction,
    )
    from ..physics.slotboom import n_from_slotboom_np, p_from_slotboom_np
    from ..postprocess import (
        evaluate_partial_currents,
        fmt_tag,
        resolve_contact_facets,
    )
    from ..results import TransientResult
    from ..runners.equilibrium import run_equilibrium
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear_block
    from ..timestepping import BDFCoefficients
    from ._common import reference_material

    # ------------------------------------------------------------------
    # Parse solver parameters
    # ------------------------------------------------------------------
    solver_cfg = cfg.get("solver", {})
    t_end = float(solver_cfg["t_end"])
    dt_val = float(solver_cfg["dt"])
    order = int(solver_cfg.get("order", 2))
    max_steps = int(solver_cfg.get("max_steps", 10000))
    output_every = int(solver_cfg.get("output_every", 50))

    bdf = BDFCoefficients(order)
    alpha_0 = bdf.coeffs[0]

    # SNES tolerances. Defaults match `run_bias_sweep` (ADR 0008): the
    # Slotboom transient residual has the same scaling as the steady-
    # state residual it shares with bias_sweep, so the same tolerances
    # apply. See ADR 0014 § "Solver options".
    snes_opts = solver_cfg.get("snes", {}) or {}
    snes_petsc_options = {
        "snes_rtol": float(snes_opts.get("rtol", 1.0e-10)),
        "snes_atol": float(snes_opts.get("atol", 1.0e-7)),
        "snes_stol": float(snes_opts.get("stol", 1.0e-14)),
        "snes_max_it": int(snes_opts.get("max_it", 100)),
        # Bump MUMPS workspace so the (small) extra delayed pivots from
        # the rank-deficient Slotboom rows do not exhaust the default
        # 20% over-allocation. Without this MUMPS errors out with
        # INFO(1)=-9 / INFOG(2)=100 ("main internal real workspace too
        # small") on the first transient step that introduces non-zero
        # spatial residual on minority-side vertices.
        "mat_mumps_icntl_14": 200,
    }

    # Diagonal Jacobian regularisation for the time loop. The Slotboom
    # continuity equation drift-diffusion current J_n = -mu_n * n *
    # grad(phi_n) carries a factor of `n = ni * exp(psi - phi_n)` which
    # is below floating-point precision in the deep p-bulk; the entire
    # (phi_n) row of the assembled Jacobian is then numerically zero
    # at those vertices and MUMPS reports null pivots. The returned
    # Newton update has unbounded components in those rows, sending
    # exp(psi - phi_n) to inf and the next residual evaluation to NaN
    # (SNES_DIVERGED_FNORM_NAN). Adding ``epsilon * I`` to the
    # assembled block Jacobian removes the null pivots without
    # measurably perturbing the converged solution at the tolerances
    # used here. 1e-14 is small enough to leave SNES rtol=1e-10
    # unaffected and large enough to dominate the ~5e-30 minimum
    # pivot reported by MUMPS without the shift. ``bias_sweep`` does
    # not need this regularisation because its voltage-continuation
    # ramp keeps the residual near zero on the same minority-side
    # rows so MUMPS never needs to invert them.
    transient_jacobian_shift = float(
        solver_cfg.get("jacobian_shift", 1.0e-14)
    )

    # ------------------------------------------------------------------
    # Build mesh, scaling, doping
    # ------------------------------------------------------------------
    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])

    phys = cfg.get("physics", {})
    mob = phys.get("mobility", {})
    mu_n_SI = float(mob.get("mu_n", 1400.0)) * 1.0e-4
    mu_p_SI = float(mob.get("mu_p", 450.0)) * 1.0e-4
    mu_n_hat = mu_n_SI / sc.mu0
    mu_p_hat = mu_p_SI / sc.mu0

    rec = dict(phys.get("recombination", {}))
    tau_n_s = float(rec.get("tau_n", 1.0e-7))
    tau_p_s = float(rec.get("tau_p", 1.0e-7))
    E_t_eV = float(rec.get("E_t", 0.0))
    tau_n_hat = tau_n_s / sc.t0
    tau_p_hat = tau_p_s / sc.t0
    E_t_over_Vt = E_t_eV / sc.V0
    # M16.6: merge physics.tunneling into recomb_cfg; see bias_sweep
    # for the rationale.
    tun = phys.get("tunneling", {})
    for key, value in tun.items():
        rec[key] = value

    # ------------------------------------------------------------------
    # Slotboom block spaces and unknowns
    # ------------------------------------------------------------------
    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi
    V_phi_n = spaces.V_phi_n
    V_phi_p = spaces.V_phi_p
    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    ni_hat_val = sc.n_i / sc.C0

    # ------------------------------------------------------------------
    # Step 0: V=0 equilibrium for psi (sets the IC reference psi when
    # no BC ramp is requested, e.g. pn_1d_turnon).
    # ------------------------------------------------------------------
    eq_cfg = copy.deepcopy(cfg)
    eq_cfg["solver"] = {"type": "equilibrium"}
    for c in eq_cfg["contacts"]:
        c["voltage"] = 0.0
        c.pop("voltage_sweep", None)
    eq_result = run_equilibrium(eq_cfg)

    psi.x.array[:] = eq_result.psi.x.array
    psi.x.scatter_forward()
    # At V=0 thermal equilibrium phi_n = phi_p = 0 (the constant of
    # integration is fixed by ohmic-contact pinning).
    phi_n.x.array[:] = 0.0
    phi_p.x.array[:] = 0.0
    phi_n.x.scatter_forward()
    phi_p.x.scatter_forward()

    # ------------------------------------------------------------------
    # Step 0b: BC-ramp continuation to V_target steady state (ADR 0013).
    #
    # The Slotboom transient is positivity-preserving by construction
    # (n, p > 0 at every Newton iterate), so the V=0-to-V_target jump
    # the (n, p) form could not handle is no longer fatal. We keep the
    # BC-ramp continuation anyway because it walks the IC to a Newton-
    # warm state at V_target, which (a) lets the time loop start
    # satisfying the steady-state limit test (ADR 0014 § Validation)
    # and (b) keeps numerical kinetic-energy noise on the very first
    # step bounded.
    #
    # Configurable via `solver.bc_ramp_steps` (default 10). Setting it
    # to 0 disables continuation and uses the V=0 equilibrium as the
    # IC; that is the right choice for `pn_1d_turnon`, where the
    # step-bias-at-t=0 is the physical scenario being measured.
    # ------------------------------------------------------------------
    bc_ramp_steps = int(solver_cfg.get("bc_ramp_steps", 10))
    bc_ramp_v_factor = float(solver_cfg.get("bc_ramp_voltage_factor", 1.0))
    if not (0.0 < bc_ramp_v_factor <= 1.0):
        raise ValueError(
            "solver.bc_ramp_voltage_factor must lie in (0.0, 1.0]; "
            f"got {bc_ramp_v_factor}"
        )
    cont = _run_bc_continuation(cfg, bc_ramp_steps, bc_ramp_v_factor)
    if cont is not None:
        psi_arr_c, phi_n_arr_c, phi_p_arr_c = cont
        psi.x.array[:] = psi_arr_c
        phi_n.x.array[:] = phi_n_arr_c
        phi_p.x.array[:] = phi_p_arr_c
        psi.x.scatter_forward()
        phi_n.x.scatter_forward()
        phi_p.x.scatter_forward()

    # ------------------------------------------------------------------
    # Build contact facet info for IV recording
    # ------------------------------------------------------------------
    contact_facet_infos: dict[str, Any] = {}
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        info = resolve_contact_facets(cfg, msh, facet_tags, c["name"])
        contact_facet_infos[c["name"]] = info

    # ------------------------------------------------------------------
    # Static voltages from contacts (no sweep in transient).
    # `static_voltages` records the time-loop-target value for each ohmic
    # contact; for contacts with `voltage_t`, that target is the waveform
    # evaluated at t=0 and is only used for IV-row labeling. The actual
    # BC value at each step comes from `voltages_at_t(t_next)`.
    # ------------------------------------------------------------------
    voltages_at_t = _build_voltage_t_evaluator(cfg)
    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    # ------------------------------------------------------------------
    # Build the transient residual
    # ------------------------------------------------------------------
    f_hist_n = fem.Function(V_phi_n, name="hist_n")
    f_hist_p = fem.Function(V_phi_p, name="hist_p")
    f_hist_n.x.array[:] = 0.0
    f_hist_p.x.array[:] = 0.0

    # dt_hat = dt / t0 is the dimensionless time step consistent with the
    # scaled spatial coefficients (L0^2 * mu_hat, etc.) in the residual.
    dt_hat = dt_val / sc.t0

    from petsc4py import PETSc

    dt_const = fem.Constant(msh, PETSc.ScalarType(dt_hat))
    alpha0_const = fem.Constant(msh, PETSc.ScalarType(alpha_0))

    stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}

    # M17: build heterojunction DG0 fields once at runner setup so the
    # transient residual sees position-dependent n_i / eps_r when the
    # cfg opts in. Configs without material_overrides / heterojunction
    # are bit-identical to v0.23.0 because het_fields stays None.
    het_fields = None
    if cfg_uses_heterojunction(cfg.get("regions", {})):
        het_fields = build_dg0_material_fields(
            msh, cell_tags, cfg["regions"], sc,
            T=phys.get("temperature", 300.0),
        )

    F_list = _build_transient_residual(
        psi, phi_n, phi_p, N_hat_fn, f_hist_n, f_hist_p,
        spaces, sc, ref_mat.epsilon_r,
        mu_n_hat, mu_p_hat, tau_n_hat, tau_p_hat, E_t_over_Vt,
        dt_const, alpha0_const,
        recomb_cfg=rec,
        statistics_cfg=stat_cfg,
        heterojunction_fields=het_fields,
    )

    # ------------------------------------------------------------------
    # BCs at ohmic contacts. Identical to bias_sweep (Slotboom Shockley
    # boundary) -- delegate to build_dd_dirichlet_bcs.
    # ------------------------------------------------------------------
    def _build_transient_bcs(voltages: dict[str, float]) -> list:
        contacts = resolve_contacts(cfg, facet_tags=facet_tags,
                                    voltages=voltages)
        return build_dd_dirichlet_bcs(
            spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
            regions_cfg=cfg.get("regions"), cell_tags=cell_tags,
        )

    def _apply_bc_values(bcs: list):
        space_to_fn = {
            id(V_psi): psi,
            id(V_phi_n): phi_n,
            id(V_phi_p): phi_p,
        }
        for bc in bcs:
            fn = space_to_fn.get(id(bc.function_space))
            if fn is not None:
                bc.set(fn.x.array)
        for fn in (psi, phi_n, phi_p):
            fn.x.scatter_forward()

    # ------------------------------------------------------------------
    # History storage. We track carrier densities at past converged
    # Slotboom states; the BDF time term acts on n and p directly.
    # ------------------------------------------------------------------
    def _eval_n_p() -> tuple[np.ndarray, np.ndarray]:
        n_arr = n_from_slotboom_np(psi.x.array, phi_n.x.array, ni_hat_val)
        p_arr = p_from_slotboom_np(psi.x.array, phi_p.x.array, ni_hat_val)
        return np.asarray(n_arr).copy(), np.asarray(p_arr).copy()

    n0, p0 = _eval_n_p()
    n_hist: list[np.ndarray] = [n0]
    p_hist: list[np.ndarray] = [p0]

    def _trim_history():
        while len(n_hist) > order:
            n_hist.pop(0)
        while len(p_hist) > order:
            p_hist.pop(0)

    # ------------------------------------------------------------------
    # IV recording helper: evaluate_partial_currents operates directly
    # on the Slotboom (psi, phi_n, phi_p) Functions in `spaces` -- no
    # carrier-density conversion needed (unlike the (n, p) form runner).
    # ------------------------------------------------------------------
    def _record_all_iv(t_val: float, iv_rows: list):
        v_now = voltages_at_t(t_val)
        for c in cfg["contacts"]:
            if c["type"] != "ohmic":
                continue
            cname = c["name"]
            V_applied = float(v_now.get(cname, static_voltages.get(cname, 0.0)))
            finfo = contact_facet_infos.get(cname)
            if finfo is None:
                continue
            J_n, J_p = evaluate_partial_currents(
                spaces, sc, ref_mat, finfo, mu_n_SI, mu_p_SI,
            )
            J_total = J_n + J_p
            iv_rows.append({
                "t": float(t_val),
                "contact": cname,
                "V": float(V_applied),
                "J": float(J_total),
                "J_n": float(J_n),
                "J_p": float(J_p),
            })

    # ------------------------------------------------------------------
    # Time loop
    # ------------------------------------------------------------------
    iv_rows: list[dict[str, Any]] = []
    t_vals: list[float] = []
    fields_out: dict[str, list[np.ndarray]] = {
        "psi": [], "n": [], "p": [], "phi_n": [], "phi_p": [],
    }
    snap_t: list[float] = []
    n_steps_taken = 0
    step_count = 0

    t_current = 0.0
    _record_all_iv(t_current, iv_rows)
    t_vals.append(t_current)

    bdf1 = BDFCoefficients(1)

    while t_current < t_end - 0.5 * dt_val and step_count < max_steps:
        t_next = t_current + dt_val
        step_count += 1

        # Decide which BDF order to use this step (BDF1 seeds BDF2 history)
        if order == 2 and len(n_hist) >= 2:
            effective_order = 2
            use_bdf = bdf
        else:
            effective_order = 1
            use_bdf = bdf1

        alpha_k_use = use_bdf.coeffs[0]
        alpha0_const.value = PETSc.ScalarType(alpha_k_use)

        # Update history source functions:
        #   f_hist_n[i] = sum_{k=1}^K alpha_k * n^{n+1-k}[i]
        # Stored on the Function (no /dt here -- the residual already
        # divides by dt_const through the alpha0/dt prefactor; we keep
        # the symmetry by absorbing alpha_k coefficients here and the
        # /dt in the residual form).
        hist_n_arr = np.zeros(len(n_hist[-1]))
        hist_p_arr = np.zeros(len(p_hist[-1]))
        for k in range(1, effective_order + 1):
            idx = -k  # n_hist[-1] = n^n (most recent), n_hist[-2] = n^{n-1}
            hist_n_arr += use_bdf.coeffs[k] * n_hist[idx]
            hist_p_arr += use_bdf.coeffs[k] * p_hist[idx]
        f_hist_n.x.array[:] = hist_n_arr
        f_hist_p.x.array[:] = hist_p_arr
        f_hist_n.x.scatter_forward()
        f_hist_p.x.scatter_forward()

        # Build BCs (Slotboom Shockley boundary) and seed unknowns
        # with the BC values at the contacts (so SNES starts close to
        # the solution at each step). For contacts without
        # `voltage_t`, `voltages_at_t(t_next)` returns the fixed
        # `voltage` value at every step (bit-identical to v0.22.0).
        bcs = _build_transient_bcs(voltages_at_t(t_next))
        _apply_bc_values(bcs)

        tag = fmt_tag(t_next)
        # `fmt_tag` rounds to 4 decimal places (sub-microvolt
        # resolution); for time tags at sub-nanosecond `t_next` it
        # collapses adjacent steps to the same string. Disambiguate
        # by appending the step counter so the PETSc options database
        # gets a fresh prefix per solve and MUMPS does not reuse
        # cached symbolic-factor state across steps.
        info = solve_nonlinear_block(
            F_list, [psi, phi_n, phi_p], bcs,
            prefix=f"{cfg['name']}_tr_{tag}_s{step_count}_",
            petsc_options=snes_petsc_options,
            jacobian_shift=transient_jacobian_shift,
            cfg=cfg,
        )

        if not info["converged"]:
            raise RuntimeError(
                f"Transient SNES failed at t={t_next:.3e} s "
                f"(step {step_count}); reason={info['reason']}"
            )

        # Append converged carrier densities to history
        n_new, p_new = _eval_n_p()
        n_hist.append(n_new)
        p_hist.append(p_new)
        _trim_history()

        t_current = t_next
        n_steps_taken += 1

        # Record IV
        _record_all_iv(t_current, iv_rows)
        t_vals.append(t_current)

        # Snapshots
        if step_count % output_every == 0 or abs(t_current - t_end) < 0.5 * dt_val:
            n_phys = n_new * sc.C0
            p_phys = p_new * sc.C0
            fields_out["psi"].append(psi.x.array.copy() * sc.V0)
            fields_out["n"].append(n_phys)
            fields_out["p"].append(p_phys)
            # ADR 0014 (Limitations subsection): expose Slotboom
            # primary unknowns so MMS rate tests can compare primary
            # unknowns directly rather than the derived (n, p) pair.
            # The internal representation in `phi_n.x.array` is the
            # dimensionless quasi-Fermi potential in scaled units;
            # we multiply by sc.V0 here so the stored snapshot is
            # in volts (matching `psi` above).
            fields_out["phi_n"].append(phi_n.x.array.copy() * sc.V0)
            fields_out["phi_p"].append(phi_p.x.array.copy() * sc.V0)
            snap_t.append(t_current)

        if progress_callback is not None:
            progress_callback({
                "type": "step_done",
                "step": step_count,
                "t": float(t_current),
                "iterations": int(info.get("iterations", 0)),
            })

    x_dof = V_psi.tabulate_dof_coordinates()

    return TransientResult(
        t=t_vals,
        iv=iv_rows,
        fields=fields_out,
        meta={
            "order": order,
            "dt": dt_val,
            "n_steps_taken": n_steps_taken,
            "n_failed_steps": 0,
            "snap_t": snap_t,
        },
        x_dof=x_dof,
    )


def _build_transient_residual(
    psi, phi_n, phi_p, N_hat_fn,
    f_hist_n, f_hist_p,
    spaces: DDBlockSpaces,
    sc,
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float,
    dt_const,
    alpha0_const,
    *,
    recomb_cfg: dict | None = None,
    statistics_cfg: dict | None = None,
    heterojunction_fields: dict | None = None,
):
    """
    Build the three-block transient residual in Slotboom form.

    The Poisson and convection-diffusion blocks are identical in shape
    to `build_dd_block_residual` (steady-state Slotboom). The new piece
    is the BDF time term applied to the carrier density expressions
    `n_ufl = n_i exp(psi - phi_n)` and `p_ufl = n_i exp(phi_p - psi)`
    via UFL automatic differentiation.

    Parameters
    ----------
    psi, phi_n, phi_p : dolfinx.fem.Function
        Slotboom unknowns (will be solved for).
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping.
    f_hist_n, f_hist_p : dolfinx.fem.Function
        History sources (updated externally each timestep):

            f_hist_u[i] = sum_{k=1}^K alpha_k * u^{n+1-k}[i].
    spaces : DDBlockSpaces
    sc : Scaling
    eps_r : float | dolfinx.fem.Function
    mu_n_over_mu0, mu_p_over_mu0, tau_n_hat, tau_p_hat, E_t_over_Vt : float
    dt_const : dolfinx.fem.Constant
        Dimensionless `dt_hat = dt / t0`.
    alpha0_const : dolfinx.fem.Constant
        BDF `alpha_0` coefficient (updated each step when switching from
        BDF1 to BDF2).

    Returns
    -------
    list of ufl.Form
        ``[F_psi, F_phi_n, F_phi_p]``.
    """
    import math as _math

    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from ..physics.slotboom import n_from_slotboom, p_from_slotboom

    V_psi = spaces.V_psi
    V_phi_n = spaces.V_phi_n
    V_phi_p = spaces.V_phi_p
    msh = V_psi.mesh

    v_psi = ufl.TestFunction(V_psi)
    v_n = ufl.TestFunction(V_phi_n)
    v_p = ufl.TestFunction(V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(msh, PETSc.ScalarType(sc.L0 ** 2))
    # M17: heterojunction DG0 fields override the scalar eps_r and
    # n_i_hat when the cfg opts in. The Slotboom n / p substitutions
    # below pick up the per-cell n_i automatically (see ADR 0016).
    if heterojunction_fields is not None:
        eps_r_ufl = heterojunction_fields["epsilon_r"]
        ni_hat_c = heterojunction_fields["n_i_hat"]
    else:
        if isinstance(eps_r, (int, float)):
            eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
        else:
            eps_r_ufl = eps_r
        ni_hat_c = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n_c = fem.Constant(msh, PETSc.ScalarType(mu_n_over_mu0))
    mu_p_c = fem.Constant(msh, PETSc.ScalarType(mu_p_over_mu0))
    tau_n_c = fem.Constant(msh, PETSc.ScalarType(tau_n_hat))
    tau_p_c = fem.Constant(msh, PETSc.ScalarType(tau_p_hat))

    # Slotboom carrier densities (always positive). Under
    # statistics_cfg["statistics"] == "fermi_dirac" the helpers apply
    # the Blakemore prefactor; the Boltzmann default is bit-identical
    # to pre-M16.4. M16.4.
    eta_offset_n_val = (
        sc.eta_offset_n
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    eta_offset_p_val = (
        sc.eta_offset_p
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    n_ufl = n_from_slotboom(
        psi, phi_n, ni_hat_c,
        statistics_cfg=statistics_cfg, eta_offset_n=eta_offset_n_val,
    )
    p_ufl = p_from_slotboom(
        psi, phi_p, ni_hat_c,
        statistics_cfg=statistics_cfg, eta_offset_p=eta_offset_p_val,
    )

    # SRH (and optional Auger, M16.3) recombination, same shape as
    # bias_sweep so the steady-state limit of the transient is bit-
    # identical to the DD form builder's output.
    n1 = ni_hat_c * _math.exp(E_t_over_Vt)
    p1 = ni_hat_c * _math.exp(-E_t_over_Vt)
    np_minus_nieq = n_ufl * p_ufl - ni_hat_c * ni_hat_c
    R = np_minus_nieq / (
        tau_p_c * (n_ufl + n1) + tau_n_c * (p_ufl + p1)
    )
    if recomb_cfg is not None and recomb_cfg.get("auger", False):
        C_n_hat_val_t = float(recomb_cfg.get("C_n", 2.8e-31)) * 1.0e-12 \
            * sc.C0 ** 2 * sc.t0
        C_p_hat_val_t = float(recomb_cfg.get("C_p", 9.9e-32)) * 1.0e-12 \
            * sc.C0 ** 2 * sc.t0
        C_n_hat_t = fem.Constant(msh, PETSc.ScalarType(C_n_hat_val_t))
        C_p_hat_t = fem.Constant(msh, PETSc.ScalarType(C_p_hat_val_t))
        R = R + (C_n_hat_t * n_ufl + C_p_hat_t * p_ufl) * np_minus_nieq

    # Poisson block (identical to bias_sweep)
    rho_hat = p_ufl - n_ufl + N_hat_fn
    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
        - rho_hat * v_psi * ufl.dx
    )

    # Lumped mass measure for the BDF time term. P1 vertex quadrature
    # collapses the consistent (Galerkin) mass matrix
    #     M[i,j] = int phi_i phi_j dx
    # to its diagonal row-sum
    #     M_diag[i] = int phi_i dx,
    # which is the standard textbook lumped-mass approximation for
    # carrier continuity equations. In the Slotboom formulation lumping
    # is essential: the consistent mass `n_ufl * v_n * dx` couples DOFs
    # through the cell-support overlap of `v_n` with carrier values that
    # span ~25 OOM across the device, and at small dt the resulting
    # Jacobian block dominates the spatial Laplacian and pushes MUMPS
    # into pivot failure (`reason=-3`). Lumping localises the
    # (psi, phi_n) chain-rule coupling to a single DOF and makes the
    # Jacobian's mass block diagonal — exactly the conditioning regime
    # MUMPS handles cleanly. See ADR 0014 § Implementation.
    dx_lump = ufl.dx(metadata={
        "quadrature_rule": "vertex",
        "quadrature_degree": 1,
    })

    # Electron continuity (Slotboom + lumped BDF time term).
    #
    # The time term is BDF on `n_ufl(psi, phi_n)` evaluated with vertex
    # quadrature. UFL's auto-derivative of `(alpha_0/dt) n_ufl v_n`
    # against (psi, phi_n) yields the chain-rule mass diagonal (n on the
    # (phi_n, phi_n) diagonal with cross-term -n into psi). See ADR 0014
    # § Implementation.
    F_phi_n = (
        (alpha0_const / dt_const) * n_ufl * v_n * dx_lump
        + (f_hist_n / dt_const) * v_n * dx_lump
        + L0_sq * mu_n_c * n_ufl
        * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * ufl.dx
        - R * v_n * ufl.dx
    )

    # Hole continuity (Slotboom + lumped BDF time term).
    F_phi_p = (
        (alpha0_const / dt_const) * p_ufl * v_p * dx_lump
        + (f_hist_p / dt_const) * v_p * dx_lump
        + L0_sq * mu_p_c * p_ufl
        * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * ufl.dx
        + R * v_p * ufl.dx
    )

    return [F_psi, F_phi_n, F_phi_p]


def _run_bc_continuation(cfg: dict, n_steps: int, voltage_factor: float = 1.0):
    """Run BC-ramp continuation from V=0 to ``voltage_factor * V_target``.

    Used by `run_transient` to obtain a near-steady-state IC for the
    time loop. See ADR 0013. The strategy is formulation-agnostic --
    it walks `(psi, phi_n, phi_p)` from V=0 to the ramp target via
    `run_bias_sweep`'s adaptive halving controller, and the converged
    Slotboom triple becomes the time-loop IC. In the (psi, phi_n, phi_p)
    formulation (ADR 0014) no Slotboom-to-density conversion is needed.

    When ``voltage_factor < 1.0`` the ramp lands at a partial bias and
    the time loop steps from this partial-bias steady state to the full
    target voltage. This is required for the M13 BDF MMS rate test:
    ramping all the way to ``V_target`` puts the IC exactly at the
    fixed point and the time loop has nothing to integrate (pairwise
    diffs collapse to machine epsilon at every dt level), while
    ``bc_ramp_steps=0`` (V=0 IC + step bias) at the MMS device
    (1e15 doping, 20 um, V_F=0.1 V) drives Newton through carrier-
    density underflow at step 1 (ratio max/min n_ufl ~ 1e51) and
    leaves the post-step Jacobian too ill-conditioned for MUMPS at
    step 2. A factor of 0.5 hits the sweet spot: IC near the fixed
    point (well-conditioned) with a bias step large enough to drive
    a measurable transient.

    Parameters
    ----------
    cfg : dict
        Validated transient JSON config. Contact target voltages are
        read from ``cfg["contacts"][i]["voltage"]``.
    n_steps : int
        Number of steady-state sub-steps used to ramp the bias. The
        ramp is delegated to ``run_bias_sweep`` with a sweep step of
        ``|V_max| / n_steps`` on the contact with the largest
        ``|target voltage|``; ``run_bias_sweep``'s adaptive controller
        applies its own halving on top if any sub-step is hard. A
        value of ``0`` disables continuation entirely (caller's IC is
        kept untouched).
    voltage_factor : float, optional
        Multiplier in ``(0, 1]`` applied to each contact target voltage
        to obtain the ramp endpoint. Defaults to ``1.0`` (ramp to the
        full target).

    Returns
    -------
    None or tuple of three numpy.ndarray
        ``None`` if no ramp is required (n_steps <= 0, or every contact
        target is within 1e-12 V of zero). Otherwise the converged
        Slotboom triple ``(psi_hat, phi_n_hat, phi_p_hat)`` at the ramp
        endpoint, sampled on the same DOF ordering as the transient
        runner's mesh.
    """
    import copy

    if n_steps <= 0:
        return None

    # For contacts with `voltage_t` (M16.7), the ramp endpoint is the
    # waveform value at t=0: `values[0]` for tables, `v0` for steps.
    # The time loop then takes over with `voltages_at_t(t_next)`. The
    # plain `voltage` field is the fallback for contacts without
    # `voltage_t`.
    targets = {
        c["name"]: _ramp_target_voltage(c) * voltage_factor
        for c in cfg["contacts"]
        if c["type"] == "ohmic"
    }
    if not targets:
        return None

    ramp_contact, V_max = max(
        targets.items(),
        key=lambda kv: abs(kv[1]),
    )
    if abs(V_max) < 1.0e-12:
        return None

    cont_cfg = copy.deepcopy(cfg)
    for c in cont_cfg["contacts"]:
        c.pop("voltage_sweep", None)
        # M16.7: drop voltage_t before handing the cfg to bias_sweep;
        # the schema rejects voltage_t on non-transient solves. The
        # ramp endpoint is the waveform value at t=0, scaled by
        # voltage_factor and written into the contact's `voltage`.
        target_v = _ramp_target_voltage(c) * voltage_factor
        c.pop("voltage_t", None)
        if c["type"] == "ohmic":
            c["voltage"] = float(target_v)
    for c in cont_cfg["contacts"]:
        if c["name"] == ramp_contact:
            c["voltage"] = V_max
            c["voltage_sweep"] = {
                "start": 0.0,
                "stop": V_max,
                "step": abs(V_max) / n_steps,
            }
            break

    # Tight steady-state SNES tolerances so the IC handed to the
    # transient time loop is at the bias_sweep fixed point to machine
    # precision. The default bias_sweep atol of 1e-7 is too loose for
    # this purpose; the steady-state agreement gate is 1e-4 relative
    # in J, so the IC must be converged well below that.
    cont_cfg["solver"] = {
        "type": "bias_sweep",
        "snes": {
            "rtol": 1.0e-14,
            "atol": 1.0e-14,
            "stol": 1.0e-14,
            "max_it": 100,
        },
        "continuation": {
            "min_step": abs(V_max) / (n_steps * 8),
            "max_halvings": 3,
            "max_step": abs(V_max) / n_steps,
            "easy_iter_threshold": 4,
            "grow_factor": 1.5,
        },
    }

    from .bias_sweep import run_bias_sweep

    bs_result = run_bias_sweep(cont_cfg)
    return (
        bs_result.psi.x.array.copy(),
        bs_result.phi_n.x.array.copy(),
        bs_result.phi_p.x.array.copy(),
    )


def _ramp_target_voltage(contact: dict[str, Any]) -> float:
    """Return the BC-ramp endpoint voltage for a contact (M16.7).

    For a contact with a `voltage_t` block, the ramp lands at the
    waveform value at t=0: `values[0]` (table) or `v0` (step). For
    contacts without `voltage_t` the target is the plain `voltage`
    field (the v0.22.0 behavior).
    """
    vt = contact.get("voltage_t")
    if vt is None:
        return float(contact.get("voltage", 0.0))
    kind = vt.get("type")
    if kind == "table":
        return float(vt["values"][0])
    if kind == "step":
        t0 = float(vt["t0"])
        return float(vt["v1"]) if t0 <= 0.0 else float(vt["v0"])
    raise ValueError(
        f"voltage_t.type={kind!r} not recognized; expected 'table' or 'step' "
        f"(M16.7)"
    )


def _build_voltage_t_evaluator(cfg: dict[str, Any]):
    """Return a callable `voltages_at_t(t) -> dict[str, float]`.

    The dict maps each ohmic contact name to its applied voltage at
    time `t` (seconds, absolute simulation time). Contacts without a
    `voltage_t` block return their fixed `voltage` value at every t,
    bit-identical to the v0.22.0 behavior. (M16.7)
    """
    fixed: dict[str, float] = {}
    table_specs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    step_specs: dict[str, tuple[float, float, float]] = {}
    for c in cfg.get("contacts", []):
        if c.get("type") != "ohmic":
            continue
        name = c["name"]
        vt = c.get("voltage_t")
        if vt is None:
            fixed[name] = float(c.get("voltage", 0.0))
            continue
        kind = vt.get("type")
        if kind == "table":
            times = np.asarray(vt["times"], dtype=float)
            values = np.asarray(vt["values"], dtype=float)
            table_specs[name] = (times, values)
        elif kind == "step":
            step_specs[name] = (
                float(vt["t0"]), float(vt["v0"]), float(vt["v1"]),
            )
        else:
            raise ValueError(
                f"contact {name!r}: voltage_t.type={kind!r} not "
                f"recognized; expected 'table' or 'step' (M16.7)"
            )

    def voltages_at_t(t: float) -> dict[str, float]:
        out: dict[str, float] = dict(fixed)
        for name, (times, values) in table_specs.items():
            out[name] = float(np.interp(t, times, values))
        for name, (t0, v0, v1) in step_specs.items():
            out[name] = v1 if t >= t0 else v0
        return out

    return voltages_at_t
