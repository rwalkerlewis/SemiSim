"""
Solver driver.

Wraps dolfinx.fem.petsc.NonlinearProblem (PETSc SNES) with default
options tuned for semiconductor problems:
    - Direct LU (MUMPS) linear solve per Newton iteration — robust for
      the block-ill-conditioned Jacobians that show up here
    - Backtracking line search — helps when initial guess is far from the
      root, especially with exponentials in the residual
    - SNES monitor on by default so we can diagnose convergence issues

For large 3D problems, LU is overkill and we'd switch to GAMG or
fieldsplit preconditioning, but that's a M6+ concern.
"""
from __future__ import annotations

import time
from typing import Any

DEFAULT_PETSC_OPTIONS = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "bt",
    "snes_rtol": 1.0e-10,
    "snes_atol": 1.0e-12,
    "snes_max_it": 50,
    "snes_monitor": None,
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}


# Option keys that dolfinx 0.10's NonlinearProblem stores in the PETSc
# options database under the SNES prefix but that PETSc never consumes
# from there. The MUMPS factor matrix is created lazily during PCSetUp,
# after SNESSetFromOptions has already run, so any
# `<prefix>mat_mumps_*` / `<prefix>pc_factor_*` entries are reported as
# "Option left" at finalize and have no effect on the solve. We extract
# these keys before passing the dict to NonlinearProblem and apply them
# via the direct petsc4py API after construction (see
# `_apply_factor_options`). The numeric suffix on the MUMPS keys is
# the (i)cntl index PETSc/MUMPS uses (e.g. mat_mumps_cntl_3 -> CNTL(3),
# the null-pivot threshold).
_MUMPS_CNTL_PREFIX = "mat_mumps_cntl_"
_MUMPS_ICNTL_PREFIX = "mat_mumps_icntl_"
_PC_FACTOR_DIRECT_KEYS = {
    "pc_factor_zeropivot",
    "pc_factor_shift_type",
    "pc_factor_shift_amount",
}


def _split_factor_options(opts: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Separate factor-matrix / MUMPS direct-API options from the rest.

    Returns ``(rest, factor_opts)`` where ``rest`` is the dict to hand
    to NonlinearProblem (consumed by SNES/KSP setFromOptions normally)
    and ``factor_opts`` collects the keys we will apply after the
    problem is constructed.
    """
    rest: dict[str, Any] = {}
    factor: dict[str, Any] = {}
    for k, v in opts.items():
        if k in _PC_FACTOR_DIRECT_KEYS:
            factor[k] = v
        elif k.startswith(_MUMPS_CNTL_PREFIX) or k.startswith(_MUMPS_ICNTL_PREFIX):
            factor[k] = v
        else:
            rest[k] = v
    return rest, factor


def _apply_factor_options(snes, factor_opts: dict[str, Any]) -> None:
    """
    Apply pc-factor and MUMPS controls via the direct petsc4py API.

    Forces the factor matrix to exist (PCFactorSetUpMatSolverType) and
    then sets the requested controls. Safe no-op when ``factor_opts``
    is empty.
    """
    if not factor_opts:
        return
    from petsc4py import PETSc

    ksp = snes.getKSP()
    pc = ksp.getPC()
    # Required so the factor matrix is materialised before we ask for it.
    pc.setFactorSetUpSolverType()

    shift_type = factor_opts.get("pc_factor_shift_type")
    shift_amount = factor_opts.get("pc_factor_shift_amount")
    if shift_type is not None or shift_amount is not None:
        st = None
        if shift_type is not None:
            st = getattr(PETSc.Mat.FactorShiftType, str(shift_type).upper())
        pc.setFactorShift(
            shift_type=st,
            amount=float(shift_amount) if shift_amount is not None else None,
        )
    if "pc_factor_zeropivot" in factor_opts:
        pc.setFactorPivot(float(factor_opts["pc_factor_zeropivot"]))

    # MUMPS controls require the factor matrix.
    needs_mumps = any(
        k.startswith(_MUMPS_CNTL_PREFIX) or k.startswith(_MUMPS_ICNTL_PREFIX)
        for k in factor_opts
    )
    if not needs_mumps:
        return
    Fmat = pc.getFactorMatrix()
    for k, v in factor_opts.items():
        if k.startswith(_MUMPS_CNTL_PREFIX):
            idx = int(k[len(_MUMPS_CNTL_PREFIX):])
            Fmat.setMumpsCntl(idx, float(v))
        elif k.startswith(_MUMPS_ICNTL_PREFIX):
            idx = int(k[len(_MUMPS_ICNTL_PREFIX):])
            Fmat.setMumpsIcntl(idx, int(v))


def apply_snes_line_search(snes, ls_type: str | None) -> None:
    """
    Set the SNES line-search type directly via petsc4py after the
    dolfinx NonlinearProblem has been constructed.

    The dolfinx 0.10 ``NonlinearProblem`` configures its SNES with
    ``SNESSetFromOptions`` before the SNES line-search context exists,
    so ``snes_linesearch_*`` options pushed through the
    ``petsc_options`` dict are silently dropped. Calling
    ``SNESLineSearch.setType`` directly on the line-search object is
    the working path (ADR 0018).

    ``ls_type`` is the PETSc line-search type name (``"bt"``,
    ``"nleqerr"``, ``"cp"``, ``"l2"``, ``"basic"``). ``None`` or
    ``""`` is a no-op so callers can pass the resolved cfg field
    unconditionally.
    """
    if not ls_type:
        return
    ls = snes.getLineSearch()
    ls.setType(str(ls_type))


def _install_jacobian_shift(snes, epsilon: float) -> None:
    """
    Wrap the SNES Jacobian callback to add ``epsilon * I`` to the
    assembled operator before the linear solve.

    Used by the transient runner to regularise rank-deficient (phi_n)
    / (phi_p) rows in the Slotboom continuity Jacobian where the
    minority carrier density is below floating-point precision (deep
    p-bulk / n-bulk respectively). Without the shift MUMPS reports
    null pivots and the returned Newton update has unbounded
    components in those rows, sending exp(psi - phi_n) to infinity
    and the next residual evaluation to NaN. The shift is small
    enough that it does not affect the converged solution at the
    tolerances used here.
    """
    if epsilon is None or epsilon == 0.0:
        return
    J_mat, P_mat, jac_callback = snes.getJacobian()
    orig_J_fn, orig_args, orig_kwargs = jac_callback
    args = orig_args or ()
    kwargs = orig_kwargs or {}

    def shifted_jacobian(snes_, X, J, P):
        orig_J_fn(snes_, X, J, P, *args, **kwargs)
        J.shift(epsilon)
        if P.handle != J.handle:
            P.shift(epsilon)

    snes.setJacobian(shifted_jacobian, J_mat, P_mat)


def _resolve_backend_options(
    cfg: dict[str, Any] | None,
    user_petsc_options: dict[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Combine backend-derived PETSc options with the runner's per-call
    overrides, and return ``(merged_options, backend_meta)``.

    ``backend_meta`` carries the resolved backend identity for the
    manifest. When ``cfg`` is None (legacy callers, in-process tests),
    the function falls back to the cpu-mumps defaults and an empty
    backend_meta dict so behavior is unchanged.
    """
    base = dict(DEFAULT_PETSC_OPTIONS)
    backend_meta: dict[str, Any] = {}
    if cfg is not None:
        from .compute import backend_settings_from_cfg
        info = backend_settings_from_cfg(cfg)
        backend_meta = {
            "backend_requested": info["requested"],
            "backend_resolved": info["resolved"],
            "device": info["device"],
            "preconditioner": info["preconditioner"],
            "linear_solver": info["linear_solver"],
        }
        for k, v in info["petsc_options"].items():
            if v is None:
                base.pop(k, None)
            else:
                base[k] = v
    if user_petsc_options:
        base.update(user_petsc_options)
    return base, backend_meta


def solve_nonlinear(F, u, bcs: list, prefix: str,
                    petsc_options: dict[str, Any] | None = None,
                    jacobian_shift: float = 0.0,
                    cfg: dict[str, Any] | None = None):
    """
    Solve a nonlinear variational problem F(u; v) = 0.

    Parameters
    ----------
    F : ufl.Form
        Residual form.
    u : dolfinx.fem.Function
        Unknown function (will be modified in place with the solution).
    bcs : list of dolfinx.fem.DirichletBC
        Dirichlet boundary conditions.
    prefix : str
        PETSc options prefix. Must be unique per problem (required in 0.10+).
    petsc_options : dict, optional
        Override DEFAULT_PETSC_OPTIONS.
    jacobian_shift : float, optional
        If non-zero, add ``jacobian_shift * I`` to the assembled
        Jacobian on every Newton step. Use this to regularise
        formulations with locally rank-deficient rows (the Slotboom
        transient time-loop is the canonical caller; see
        ``_install_jacobian_shift``). Default ``0.0`` (no shift).

    Returns
    -------
    dict
        {'iterations': int, 'reason': int, 'converged': bool}
    """
    from dolfinx.fem.petsc import NonlinearProblem

    opts, backend_meta = _resolve_backend_options(cfg, petsc_options)
    rest, factor_opts = _split_factor_options(opts)

    problem = NonlinearProblem(
        F, u,
        bcs=bcs,
        petsc_options_prefix=prefix,
        petsc_options=rest,
    )
    _apply_factor_options(problem.solver, factor_opts)
    _install_jacobian_shift(problem.solver, jacobian_shift)
    _t0 = time.monotonic()
    problem.solve()
    linear_solve_wall_s = time.monotonic() - _t0
    reason = problem.solver.getConvergedReason()
    n_iter = problem.solver.getIterationNumber()
    try:
        ksp_iters = int(problem.solver.getKSP().getIterationNumber())
    except Exception:
        ksp_iters = 0
    converged = reason > 0
    out = {
        "iterations": int(n_iter),
        "reason": int(reason),
        "converged": bool(converged),
        "ksp_iters": ksp_iters,
        "linear_solve_wall_s": float(linear_solve_wall_s),
        "problem": problem,  # kept so caller can inspect further
    }
    out.update(backend_meta)
    return out


def solve_nonlinear_block(
    F_list,
    u_list,
    bcs: list,
    prefix: str,
    petsc_options: dict[str, Any] | None = None,
    kind: str | None = None,
    entity_maps: list | None = None,
    jacobian_shift: float = 0.0,
    cfg: dict[str, Any] | None = None,
    line_search: str | None = None,
):
    """
    Solve a coupled block nonlinear problem via SNES.

    Wraps dolfinx 0.10's `NonlinearProblem` with its blocked-mode
    signature (F and u passed as sequences, plus optional `kind`).
    Dirichlet BCs are passed as a single flat list; each BC is tied to
    the subspace it was built against.

    Parameters
    ----------
    F_list : sequence of ufl.Form
        One residual form per block (for example [F_psi, F_phi_n, F_phi_p]).
    u_list : sequence of dolfinx.fem.Function
        Unknown Functions for each block (updated in place on return).
    bcs : list of dolfinx.fem.DirichletBC
        Flat list of Dirichlet BCs.
    prefix : str
        PETSc options prefix (must end with `_`).
    petsc_options : dict, optional
        Override DEFAULT_PETSC_OPTIONS.
    kind : str, optional
        PETSc matrix/vector kind. `None` for monolithic blocked (default,
        works with the MUMPS LU solver). Use `"nest"` if you supply a
        fieldsplit preconditioner.
    entity_maps : list of dolfinx.mesh.EntityMap, optional
        Passed through to `NonlinearProblem` for mixed-domain assembly
        (for example, when psi lives on the parent mesh and phi_n /
        phi_p live on a semiconductor submesh). `None` (the default)
        leaves the single-region path unchanged.
    jacobian_shift : float, optional
        Add ``jacobian_shift * I`` to the assembled block Jacobian on
        every Newton step. Used by the transient runner (~1e-12) to
        keep the LU factorisation well-defined where the Slotboom
        continuity equation has rank-deficient minority-side rows.
        Default ``0.0``.

    Returns
    -------
    dict
        Same keys as :func:`solve_nonlinear`.
    """
    from dolfinx.fem.petsc import NonlinearProblem

    opts, backend_meta = _resolve_backend_options(cfg, petsc_options)
    rest, factor_opts = _split_factor_options(opts)

    np_kwargs: dict[str, Any] = dict(
        bcs=list(bcs),
        petsc_options_prefix=prefix,
        petsc_options=rest,
        kind=kind,
    )
    if entity_maps is not None:
        np_kwargs["entity_maps"] = list(entity_maps)

    problem = NonlinearProblem(
        list(F_list), list(u_list),
        **np_kwargs,
    )
    _apply_factor_options(problem.solver, factor_opts)
    _install_jacobian_shift(problem.solver, jacobian_shift)
    apply_snes_line_search(problem.solver, line_search)
    _t0 = time.monotonic()
    problem.solve()
    linear_solve_wall_s = time.monotonic() - _t0
    reason = problem.solver.getConvergedReason()
    n_iter = problem.solver.getIterationNumber()
    try:
        ksp_iters = int(problem.solver.getKSP().getIterationNumber())
    except Exception:
        ksp_iters = 0
    converged = reason > 0
    out = {
        "iterations": int(n_iter),
        "reason": int(reason),
        "converged": bool(converged),
        "ksp_iters": ksp_iters,
        "linear_solve_wall_s": float(linear_solve_wall_s),
        "problem": problem,
    }
    out.update(backend_meta)
    return out
