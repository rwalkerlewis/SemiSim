"""
Coupled drift-diffusion runner with adaptive bias continuation.

Solves the three-block (psi, phi_n, phi_p) Slotboom residual under SRH
recombination, ramping the bias from V=0 to the final sweep endpoint
through `semi.continuation.AdaptiveStepController`. Per-step BC
construction is delegated to `semi.bcs`; per-step current evaluation
and IV recording are delegated to `semi.postprocess`.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def run_bias_sweep(
    cfg: dict[str, Any],
    *,
    post_step_hook=None,
    progress_callback=None,
):
    """
    Coupled drift-diffusion solver with optional bias ramping.

    `post_step_hook`, if provided, is called after every successful
    bias step with the signature

        post_step_hook(V_applied, spaces, iv_row)

    where `iv_row` is the freshly-appended dict in `result.iv`. The
    hook is expected to mutate `iv_row` in place (e.g. to record
    Phase 3 V&V continuity samples on top of the V / J the sweep
    already records). No return value is inspected.

    `progress_callback`, if provided, is called with a plain dict per
    successful bias step so the M10 HTTP server can surface SNES
    progress on its WebSocket channel.
    """
    from dolfinx import fem

    from ..bcs import build_dd_dirichlet_bcs, resolve_contacts
    from ..continuation import AdaptiveStepController, StepTooSmall
    from ..doping import build_profile
    from ..mesh import build_mesh
    from ..physics.drift_diffusion import build_dd_block_residual, make_dd_block_spaces
    from ..postprocess import (
        fmt_tag,
        record_iv,
        resolve_contact_facets,
    )
    from ..run import SimulationResult
    from ..scaling import make_scaling_from_config
    from ..solver import solve_nonlinear_block
    from ._common import reference_material

    ref_mat = reference_material(cfg)
    sc = make_scaling_from_config(cfg, ref_mat)

    msh, _cell_tags, facet_tags = build_mesh(cfg)

    N_raw_fn = build_profile(cfg["doping"])

    spaces = make_dd_block_spaces(msh)
    V_psi = spaces.V_psi

    N_hat_fn = fem.Function(V_psi, name="N_net_hat")
    N_hat_fn.interpolate(lambda x: N_raw_fn(x) / sc.C0)

    two_ni = 2.0 * ref_mat.n_i
    spaces.psi.interpolate(lambda x: np.arcsinh(N_raw_fn(x) / two_ni))
    spaces.phi_n.x.array[:] = 0.0
    spaces.phi_p.x.array[:] = 0.0
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.scatter_forward()

    phys = cfg.get("physics", {})

    # M16.5: configs with Schottky contacts need a proper equilibrium
    # Poisson pre-solve to set the psi profile near the metal-
    # semiconductor barrier; the doping-based asinh seed is correct in
    # the bulk but ignores the band-bending caused by the Schottky
    # Dirichlet, leaving the coupled DD seed solve far from a valid
    # Newton basin. The Poisson pre-solve uses the same psi Dirichlet
    # the DD solve does, populating spaces.psi with the depletion-
    # region profile so the DD seed solve can converge in a few
    # iterations. Devices with no Schottky contacts skip this branch
    # (bit-identical to v0.20.0).
    has_schottky_pre = any(
        c["type"] == "schottky" for c in cfg["contacts"]
    )
    if has_schottky_pre:
        from ..bcs import build_psi_dirichlet_bcs
        from ..physics.poisson import build_equilibrium_poisson_form
        from ..solver import solve_nonlinear

        contacts_eq = resolve_contacts(cfg, facet_tags=facet_tags)
        bcs_psi_only = build_psi_dirichlet_bcs(
            spaces.V_psi, msh, facet_tags, contacts_eq, sc, ref_mat,
            N_raw_fn,
        )
        for bc in bcs_psi_only:
            bc.set(spaces.psi.x.array)
        spaces.psi.x.scatter_forward()
        F_psi_eq = build_equilibrium_poisson_form(
            spaces.V_psi, spaces.psi, N_hat_fn, sc, ref_mat.epsilon_r,
            statistics_cfg={"statistics": phys.get("statistics", "boltzmann")},
        )
        solve_nonlinear(
            F_psi_eq, spaces.psi, bcs_psi_only,
            prefix=f"{cfg['name']}_schottky_psi_seed_", cfg=cfg,
        )

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
    # M16.6: merge the physics.tunneling sub-dict into the recomb_cfg
    # the form builder consumes, so the bbt/tat flags and Kane/Hurkx
    # parameters live alongside auger / C_n / C_p. The merge favours
    # tunneling-side keys on collision (none in practice; the schema
    # forbids overlapping property names via additionalProperties:
    # false on each sub-object).
    tun = phys.get("tunneling", {})
    for key, value in tun.items():
        rec[key] = value

    sweep_contact, sweep_values = _resolve_sweep(cfg)

    cont = cfg.get("solver", {}).get("continuation", {})
    max_halvings = int(cont.get("max_halvings", 6))
    min_step = float(cont.get("min_step", 1.0e-4))
    easy_iter_threshold = int(cont.get("easy_iter_threshold", 4))
    grow_factor = float(cont.get("grow_factor", 1.5))

    # SNES tolerances for the coupled DD block solve. Defaults match the
    # M12 close-out values documented in docs/adr/0008-snes-tolerances.md.
    # Per-config overrides live under solver.snes so legacy benchmarks
    # (pn_1d_bias, pn_1d_bias_reverse) can keep machine-tight convergence
    # required to pass <5% / <15% current-continuity gates.
    snes_opts = cfg.get("solver", {}).get("snes", {}) or {}
    snes_petsc_options = {
        "snes_rtol": float(snes_opts.get("rtol", 1.0e-10)),
        "snes_atol": float(snes_opts.get("atol", 1.0e-7)),
        "snes_stol": float(snes_opts.get("stol", 1.0e-14)),
        "snes_max_it": int(snes_opts.get("max_it", 100)),
    }

    stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}

    # M16.5: the Schottky thermionic-emission Robin form depends on
    # the per-contact applied voltage, so the residual is rebuilt per
    # bias step inside `solve_at`. Devices with no Schottky contacts
    # see an empty list and the Robin branch is skipped (bit-identical
    # to v0.20.0).
    has_schottky = any(c["type"] == "schottky" for c in cfg["contacts"])

    def _build_F_list(voltages_for_step: dict[str, float]):
        schottky_facets_step: list[tuple[int, dict]] = []
        if has_schottky:
            contacts_step = resolve_contacts(
                cfg, facet_tags=facet_tags, voltages=voltages_for_step,
            )
            for c in contacts_step:
                if c.kind != "schottky":
                    continue
                schottky_facets_step.append(
                    (int(c.facet_tag), {
                        "barrier_height_eV": float(c.barrier_height_eV),
                        "V_applied": float(c.V_applied),
                    })
                )
        return build_dd_block_residual(
            spaces, N_hat_fn, sc, ref_mat.epsilon_r,
            mu_n_hat, mu_p_hat, tau_n_hat, tau_p_hat, E_t_over_Vt,
            mobility_cfg=mob,
            facet_tags=facet_tags,
            recomb_cfg=rec,
            statistics_cfg=stat_cfg,
            schottky_facets=schottky_facets_step,
            ref_mat=ref_mat,
        )

    F_list = _build_F_list({})

    # Ohmic, gate, and (since M16.5) Schottky contacts can carry
    # static or swept voltages. The runner historically only swept
    # ohmic contacts; gate-sweep support was added so the mosfet_2d
    # benchmark can drive V_GS at fixed V_DS for the M14.3 Pao-Sah
    # verifier. M16.5 widened the dispatch to include `schottky` so the
    # schottky_1d benchmark can sweep its Schottky anode. Insulating
    # contacts contribute no Dirichlet row and stay untouched.
    static_voltages: dict[str, float] = {}
    for c in cfg["contacts"]:
        if c["type"] not in ("ohmic", "gate", "schottky"):
            continue
        if sweep_contact is not None and c["name"] == sweep_contact:
            continue
        static_voltages[c["name"]] = float(c.get("voltage", 0.0))

    if sweep_contact is None:
        v_sweep_list = [0.0]
        baked = float(next(
            (c.get("voltage", 0.0) for c in cfg["contacts"]
             if c["type"] == "ohmic" and c.get("voltage", 0.0) != 0.0),
            0.0,
        ))
        if baked != 0.0:
            v_sweep_list.append(baked)
    else:
        v_sweep_list = list(sweep_values)
        if not v_sweep_list or v_sweep_list[0] != 0.0:
            v_sweep_list = [0.0] + [v for v in v_sweep_list if v != 0.0]

    bipolar_legs = compute_bipolar_legs(v_sweep_list)

    # Nominal step size from the sweep spacing (e.g., voltage_sweep.step).
    # continuation.max_step defaults to this, which preserves the M2
    # halving-only behaviour when growth cannot exceed the initial step.
    if len(v_sweep_list) >= 2:
        nominal_step_abs = abs(v_sweep_list[1] - v_sweep_list[0])
    else:
        nominal_step_abs = 0.0
    max_step_abs = float(cont.get("max_step", nominal_step_abs))

    def snapshot():
        return (
            spaces.psi.x.array.copy(),
            spaces.phi_n.x.array.copy(),
            spaces.phi_p.x.array.copy(),
        )

    def restore(snap):
        spaces.psi.x.array[:] = snap[0]
        spaces.phi_n.x.array[:] = snap[1]
        spaces.phi_p.x.array[:] = snap[2]
        for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
            fn.x.scatter_forward()

    def solve_at(V_by_contact: dict[str, float], tag: str):
        contacts = resolve_contacts(cfg, facet_tags=facet_tags, voltages=V_by_contact)
        bcs = build_dd_dirichlet_bcs(
            spaces, msh, facet_tags, contacts, sc, ref_mat, N_raw_fn,
        )
        space_to_fn = {
            id(spaces.V_psi): spaces.psi,
            id(spaces.V_phi_n): spaces.phi_n,
            id(spaces.V_phi_p): spaces.phi_p,
        }
        for bc in bcs:
            fn = space_to_fn.get(id(bc.function_space))
            if fn is not None:
                bc.set(fn.x.array)
        for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
            fn.x.scatter_forward()
        F_list_step = _build_F_list(V_by_contact) if has_schottky else F_list
        return solve_nonlinear_block(
            F_list_step, [spaces.psi, spaces.phi_n, spaces.phi_p],
            bcs, prefix=f"{cfg['name']}_dd_{tag}_",
            petsc_options=snes_petsc_options,
            cfg=cfg,
        )

    sweep_facet_info = None
    if sweep_contact is not None:
        sweep_facet_info = resolve_contact_facets(
            cfg, msh, facet_tags, sweep_contact,
        )

    # Per-step current recording at all ohmic contacts (in addition to the
    # sweep contact). This lets the mosfet_2d Pao-Sah verifier read J_drain
    # while the gate is the swept contact (which has J ~ 0 in DC). The
    # legacy single-ohmic benchmarks (pn_1d_bias, etc.) still see the same
    # iv_row["V"] / iv_row["J"] values; the new entries are namespaced as
    # iv_row["J_<contact_name>"] and are additive.
    monitor_contacts: list[tuple[str, dict]] = []
    for c in cfg["contacts"]:
        if c["type"] != "ohmic":
            continue
        info_facet = resolve_contact_facets(cfg, msh, facet_tags, c["name"])
        monitor_contacts.append((c["name"], info_facet))

    def _record_per_contact_currents(iv_row: dict) -> None:
        from ..postprocess import evaluate_current_at_contact
        for name, finfo in monitor_contacts:
            J_name = evaluate_current_at_contact(
                spaces, sc, ref_mat, finfo, mu_n_SI, mu_p_SI,
            )
            iv_row[f"J_{name}"] = float(J_name)

    iv_rows: list[dict[str, float]] = []
    last_info: dict[str, Any] = {}

    # Seed solve at the first sweep entry (V=0 by construction above).
    V_seed = float(v_sweep_list[0])
    voltages = dict(static_voltages)
    if sweep_contact is not None:
        voltages[sweep_contact] = V_seed
    snap = snapshot()
    info = solve_at(voltages, fmt_tag(V_seed))
    if not info["converged"]:
        restore(snap)
        raise RuntimeError(f"SNES failed at seed bias V={V_seed:+.4f} V")
    last_info = info
    V_prev = V_seed
    record_iv(iv_rows, V_seed, spaces, sc, ref_mat,
              sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI)
    _record_per_contact_currents(iv_rows[-1])
    if post_step_hook is not None:
        post_step_hook(V_seed, spaces, iv_rows[-1])
    if progress_callback is not None:
        progress_callback({
            "type": "step_done", "bias_step": 0, "V_applied": float(V_seed),
            "iterations": int(info.get("iterations", 0)),
        })

    def ramp_leg(V_start: float, V_target: float):
        """Adaptive single-direction ramp from V_start to V_target.

        Mutates `spaces` in place, appends one iv row per successful step,
        and returns the SNES info from the last accepted solve. The
        AdaptiveStepController is local to the leg, so a bipolar sweep
        re-seeds growth/halving state at each turning point.
        """
        nonlocal last_info, V_prev
        if abs(V_target - V_start) <= 1.0e-12 or max_step_abs <= 0.0:
            return
        sign = 1.0 if V_target > V_start else -1.0
        initial_abs = nominal_step_abs if nominal_step_abs > 0.0 else max_step_abs
        initial_abs = min(initial_abs, max_step_abs)
        initial_step = sign * initial_abs

        controller = AdaptiveStepController(
            initial_step=initial_step,
            max_step_abs=max_step_abs,
            min_step_abs=min_step,
            easy_iter_threshold=easy_iter_threshold,
            grow_factor=grow_factor,
        )

        halvings = 0
        while abs(V_target - V_prev) > 1.0e-12:
            remaining = V_target - V_prev
            step = controller.clamp_to_endpoint(remaining)
            V_try = V_prev + step
            voltages = dict(static_voltages)
            if sweep_contact is not None:
                voltages[sweep_contact] = V_try
            snap = snapshot()
            try:
                info = solve_at(voltages, fmt_tag(V_try))
                converged = info["converged"]
            except Exception:
                converged = False
                info = {"converged": False, "iterations": -1, "reason": -99}

            if converged:
                last_info = info
                V_prev = V_try
                record_iv(iv_rows, V_try, spaces, sc, ref_mat,
                          sweep_contact, sweep_facet_info, mu_n_SI, mu_p_SI)
                _record_per_contact_currents(iv_rows[-1])
                if post_step_hook is not None:
                    post_step_hook(V_try, spaces, iv_rows[-1])
                if progress_callback is not None:
                    progress_callback({
                        "type": "step_done",
                        "bias_step": len(iv_rows) - 1,
                        "V_applied": float(V_try),
                        "iterations": int(info.get("iterations", 0)),
                    })
                halvings = 0
                controller.on_success(int(info.get("iterations", 0)))
                continue

            restore(snap)
            halvings += 1
            try:
                controller.on_failure()
            except StepTooSmall as exc:
                raise RuntimeError(
                    f"Bias ramp failed near V={V_try:+.4f} V: {exc}"
                ) from exc
            if halvings > max_halvings:
                raise RuntimeError(
                    f"Bias ramp failed near V={V_try:+.4f} V after "
                    f"{halvings} halvings (min_step={min_step})."
                )

    if bipolar_legs:
        # Walk: V_seed (=0) -> most-negative -> most-positive.
        for V_target in bipolar_legs:
            ramp_leg(V_prev, float(V_target))
    else:
        ramp_leg(V_prev, float(v_sweep_list[-1]))

    x_dof = V_psi.tabulate_dof_coordinates()
    psi_hat = spaces.psi.x.array
    phi_n_hat = spaces.phi_n.x.array
    phi_p_hat = spaces.phi_p.x.array
    n_hat = (ref_mat.n_i / sc.C0) * np.exp(psi_hat - phi_n_hat)
    p_hat = (ref_mat.n_i / sc.C0) * np.exp(phi_p_hat - psi_hat)

    return SimulationResult(
        cfg=cfg, mesh=msh, V=V_psi,
        psi=spaces.psi, phi_n=spaces.phi_n, phi_p=spaces.phi_p,
        psi_phys=psi_hat * sc.V0,
        phi_n_phys=phi_n_hat * sc.V0,
        phi_p_phys=phi_p_hat * sc.V0,
        n_phys=n_hat * sc.C0,
        p_phys=p_hat * sc.C0,
        x_dof=x_dof, N_hat=N_hat_fn, scaling=sc,
        solver_info=last_info, iv=iv_rows, bias_contact=sweep_contact,
    )


def compute_bipolar_legs(v_sweep_list: list[float]) -> list[float]:
    """Endpoints for a sign-spanning bias walk, empty for unipolar sweeps.

    A sweep list that contains both a negative entry and a positive entry
    (ignoring numerical-noise values within 1e-12 V of zero) is walked
    in two single-direction legs: V = 0 -> min(V) -> max(V). Unipolar
    sweeps return `[]`, preserving the original single-endpoint ramp
    path used by the pn-junction and MOS benchmarks.
    """
    has_neg = any(v < -1.0e-12 for v in v_sweep_list)
    has_pos = any(v > 1.0e-12 for v in v_sweep_list)
    if has_neg and has_pos:
        return [float(min(v_sweep_list)), float(max(v_sweep_list))]
    return []


def _resolve_sweep(cfg):
    """Return (contact_name, list_of_voltage_values) or (None, []) if none.

    Ohmic, gate, and (since M16.5) Schottky contacts are sweepable.
    Precedence follows the JSON contact order; in practice a
    benchmark declares at most one swept contact.
    """
    for c in cfg["contacts"]:
        if c["type"] not in ("ohmic", "gate", "schottky"):
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
