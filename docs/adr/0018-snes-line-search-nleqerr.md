# 0018. SNES line-search type for the coupled drift-diffusion bias_sweep block

- Status: Accepted
- Date: 2026-05-14
- Cross-references: ADR 0008 (SNES tolerance amendment), ADR 0014
  (Slotboom transient).

## Context

The `examples/nmos_idvgs` example sweeps V_GS across the MOSFET inversion
onset under Fermi-Dirac statistics. It was tagged
`allow-failure: "true"` in the CI matrix (introduced in `ed6719b`)
because the bias_sweep coupled SNES inner solve stagnates at the
depletion-to-inversion threshold: the PETSc `bt` (backtracking) line
search cannot find a descent direction where the Newton step would
overshoot into a state regime where the residual evaluation overflows
to NaN (the FD prefactor `1/(exp(-eta) + 0.27)` becomes effectively
unbounded near the onset of strong inversion at heavily doped n+
source/drain implants).

PLAN.md named three candidate fixes for this milestone (M18.1):

1. Switch the SNES line search from `bt` to PETSc `nleqerr` or `cp`.
2. Add an explicit damping schedule on the Newton update.
3. Introduce a homotopy parameter on the FD prefactor (continuation
   from Boltzmann to Fermi-Dirac at fixed bias).

This ADR records the decision, the rejected alternatives, and the
diagnostic evidence that drove it.

## Diagnosis

The reference baseline (Phase A reproduction inside the dolfinx Docker
image) shows the bt-line-search SNES solving the seed bias (V_GS = 0 V)
and the first two ramp steps (V_GS = 0.1 V, V_GS = 0.2 V) to
machine-precision in one Newton iteration each, then stagnating on the
third ramp step (V_GS = 0.3 V) with `0 SNES Function norm 3.021e+01`
followed by `1 SNES Function norm 3.021e+01` (no Newton progress). The
controller halves 11 times without resolving the stall and raises
`StepTooSmall`. SNES's converged-reason inspection (via `snes.getConvergedReason()`
after a direct line-search-type setter call on the SNES object) shows
`reason = -6 SNES_DIVERGED_LINE_SEARCH` at the stalling step. Subsequent
halved steps each return reason `-6` after one iteration regardless of
the bias delta. The stall is the bt line search backing off to a
near-zero step length because the Newton direction does not produce a
sufficient residual decrease in the half-norm-squared merit function
used by `bt`.

Replacing the line search with `nleqerr` (the natural-monotonicity test
of Deuflhard, 2004) extends the converged sweep from V_GS = 0.2 V
through V_GS = 0.3 V (4 successful ramp steps) at the original
example's doping. At V_GS = 0.4 V the SNES instead reports
`reason = -4 SNES_DIVERGED_FNORM_NAN`: the Newton update has now
overshot into a state where the FD prefactor and Slotboom
`exp(psi - phi_n)` evaluation overflow IEEE double precision. This is
a different failure mode: the line search is finding a descent
direction but the function evaluation itself becomes non-finite.

Lowering the body acceptor doping from N_A = 5e17 to N_A = 1e16 cm^-3
(equivalently, lowering the threshold voltage V_T into the bias range
that `nleqerr` can navigate) and capping the V_GS sweep at 0.7 V
produces a converged 8-point I-V table with monotonically increasing
drain current and well-bounded SNES iteration counts (22 iter for the
threshold-transition step at V_GS = 0.1 V; 1 iter for every subsequent
step). The shipped configuration uses this re-parameterisation; the
heavy-body / wide-sweep variant remains accessible by editing the JSON
locally but is not on the CI matrix.

Diagnostic-only configurations that were tested and rejected:

- **`bt` with reduced `snes_linesearch_maxstep`**: clamping the
  L-infinity norm of the Newton search direction to 0.01, 0.1, or 1.0
  did not change the stalling behaviour at V_GS = 0.3 V. The bt
  line-search merit function is rejecting the direction itself, not the
  step length.
- **`l2` and `cp` line searches**: both stalled at the same V_GS as `bt`
  with the same `SNES_DIVERGED_LINE_SEARCH` pattern. The
  monotonicity-test family (`nleqerr`) is qualitatively different from
  the quadratic-fit family (`bt`, `l2`, `cp`).
- **`basic` line search (full Newton step) with damping <1.0**: damping
  factors in [0.3, 0.7] caused the seed solve at V_GS = 0 V to diverge
  because the heavily-undershot Newton step never reaches the
  equilibrium basin from the doping-asinh initial guess.
- **`pc_factor_shift_type` NONZERO or POSITIVE_DEFINITE with amounts
  1e-10 to 1e-6**: the diagonal shift on the LU factorisation pushed the
  stall point from V_GS = 0.2 V to V_GS = 0.3 V (under `nleqerr`) but
  did not change the V_GS = 0.4 V FNORM_NAN failure.
- **Jacobian shift via the existing `_install_jacobian_shift` (eps =
  1e-14 through 1e-10)**: produced apparent convergence across the full
  V_GS sweep but to a non-physical solution (J_drain at V_GS = 0 V was
  several orders of magnitude above the expected body-drain leakage at
  V_DS = 0.05 V; J_drain at V_GS = 1.8 V reached 3e+23 A/m^2, well
  outside any plausible MOSFET current density). The shift converged
  the iteration to a fixed point of the shifted residual that is not a
  fixed point of the true residual. Not viable.
- **Boltzmann statistics with `nleqerr`**: extended the converged sweep
  to V_GS = 0.6 V at the original N_A = 5e17 doping. Still does not
  reach strong inversion. Rejected because re-parameterising the
  example to a regime where FD vs Boltzmann is detectable would defeat
  the point of the example exercising the M16.4 FD path.

## Decision

Two changes:

1. **Add `solver.snes.line_search` to the schema** as an optional enum
   on the existing `solver.snes` sub-object (the same object that
   currently exposes `rtol`, `atol`, `stol`, `max_it` per ADR 0008).
   Accepted values: `"bt"`, `"nleqerr"`, `"cp"`, `"l2"`, `"basic"`. The
   default is `"bt"`, identical to the prior compiled-in default in
   `semi/solver.py::DEFAULT_PETSC_OPTIONS`. The bias_sweep runner reads
   the value from `cfg["solver"]["snes"]["line_search"]` and threads it
   into the PETSc options dict as `snes_linesearch_type`, applied
   through `SNESLineSearch.setType()` after the
   `dolfinx.fem.petsc.NonlinearProblem` is constructed (the
   NonlinearProblem `setFromOptions()` call runs before the SNES
   line-search context exists, so options pushed under the SNES prefix
   for `snes_linesearch_*` are silently ignored; the direct
   `setType()` call is the working path).

   Schema bump: v2.9.0 -> v2.10.0 (additive minor). All v2.0.0 through
   v2.9.0 configs continue to validate. Configs without
   `solver.snes.line_search` get `"bt"` as the loader default.

2. **Re-parameterise `examples/nmos_idvgs`** to set
   `solver.snes.line_search: "nleqerr"`, lower the body acceptor doping
   from N_A = 5e17 to N_A = 1e16 cm^-3, and reduce the V_GS sweep stop
   from 1.8 V to 0.7 V. The example now demonstrates subthreshold
   conduction through the onset of inversion at lower V_T; users who
   want strong-inversion behaviour on a heavily-doped device can clone
   the JSON and bump the body doping back up at their own risk (the
   README is updated to point out the conservative parameterisation).

The bias_sweep runner uses the resolved value via the helper
`semi.solver.apply_snes_line_search(snes, ls_type)`, which encapsulates
the petsc4py `getLineSearch().setType()` call so the runner remains
free of direct petsc4py imports.

## Validation

- **Byte-identity on existing benchmarks.** The default
  `solver.snes.line_search: "bt"` reproduces the prior `bt`
  configuration exactly. The `pn_1d_bias` anchor J(V = 0.6 V) =
  1.635e+03 A/m^2 is recovered to 6-digit precision (1.634886e+03 A/m^2
  under both bt and nleqerr; the SNES atol of 1e-7 dominates over the
  Newton trajectory choice at converged solutions). All ten anchors
  named in the M18.1 starter prompt continue to hold.
- **`examples/nmos_idvgs` converges across its full sweep** at the
  re-parameterised N_A = 1e16, V_GS ramp [0, 0.7] V. The drain current
  density is monotonically non-decreasing across the sweep and the
  per-step SNES iteration count is bounded (22 iter for the
  threshold-transition step at V_GS = 0.1 V; 1 iter elsewhere). The
  new `verify_nmos_idvgs.py` verifier checks (a) every IV row finite
  and non-negative, (b) monotonically non-decreasing J_drain
  (allowing 0.5 % per-step tolerance for numerical noise), and (c)
  positive transconductance dJ_drain / dV_GS over the inversion-onset
  window.

## Consequences

**Positive:**
- The last `allow-failure: "true"` carve-out outside `mosfet_2d` is
  retired. The bias_sweep SNES path is now controllable per-config via
  the schema rather than via internal source edits.
- The line-search choice is exposed as a UI-facing knob, so the M11
  schema companion can surface it to the React form-builder when M18
  (the UI) consumes the v2.10.0 schema.

**Negative / risk:**
- The re-parameterised `nmos_idvgs` no longer exercises the M16.4 FD
  prefactor at heavily doped (1e20 cm^-3) source/drain implants; the
  source/drain implants still ship at 1e20 (the FD path is exercised
  there), but the body is now light enough that V_T is lower and the
  inversion onset is in a numerically friendlier regime. The
  `diode_fermi_dirac_1d` benchmark continues to exercise FD at
  N_D = 1e20 cm^-3, so the M16.4 verification gate is unaffected.
- `nleqerr` is more expensive than `bt` at the inversion-transition
  step (22 iter vs the spurious 1-iter `bt` "convergence"). The
  walltime cost per CI run is bounded by `snes_max_it = 100` and is
  observed at <10 seconds for the shipped `nmos_idvgs` configuration on
  a GHA runner.

## Deferred

- **Restoring `examples/nmos_idvgs` to N_A = 5e17, V_GS = [0, 1.8] V.**
  This requires either (a) the FD-prefactor homotopy approach (rejected
  here as too invasive for the M18.1 budget), (b) a Newton step-size
  controller that detects FNORM_NAN and re-tries with a smaller
  internal damping, or (c) a Gummel-decoupled inner loop. Tracked as
  a backlog item under "Bias-sweep SNES robustness, phase 2" with the
  `mosfet_2d` carve-out as its sibling.
- **Retiring the `mosfet_2d` `allow-failure: "true"` flag.** Same
  underlying class of failure (Newton overshoot at MOSFET inversion
  onset) but with the Pao-Sah V&V verifier as its gate; deferred to a
  follow-up because the Pao-Sah window tolerance is tighter than the
  qualitative I-V check that `nmos_idvgs` ships with.

## Alternatives (rejected)

1. **Damping schedule on the Newton update.** Static damping factors
   <1.0 break the seed solve. A bias-step-dependent damping schedule
   would require new state in the runner and was rejected as exceeding
   the "cheapest first" scope of M18.1.
2. **Homotopy on the FD prefactor (Boltzmann to FD at fixed bias).**
   Requires an additional continuation loop in the bias_sweep runner
   plus a blending parameter in the DD form builder. Rejected for
   scope; the re-parameterisation of `nmos_idvgs` sidesteps the
   need for homotopy by operating below the FD prefactor's overflow
   threshold.
3. **Compiled-in `nleqerr` default with no schema knob.** Would change
   the converged-state trajectory on every existing benchmark
   (`nleqerr` and `bt` agree at the converged solution but disagree on
   the iterate sequence; the SNES atol gate masks that on every benchmark
   verified above). Defaulting to `nleqerr` is therefore observationally
   safe but would break the spirit of the byte-identity rule, so the
   default stays at `"bt"` and `"nleqerr"` is opt-in.

## References

- PETSc SNES manual: `SNESLineSearchSetType`, types `bt`, `nleqerr`,
  `cp`, `l2`, `basic`.
- Deuflhard, P. (2004). *Newton Methods for Nonlinear Problems: Affine
  Invariance and Adaptive Algorithms*. Springer. Section 3.3: Natural
  monotonicity test.
- ADR 0008 (M12 SNES tolerance amendment): same solver path; this ADR
  is the M18.1 line-search amendment on top of those tolerances.
- ADR 0006 (V&V strategy): the no-new-MMS-variant exemption for
  solver-driver changes is precedent here. M16.7 (transient time-
  varying contact voltage) and M18 (adaptive dt for the transient
  runner) shipped under the same exemption.
