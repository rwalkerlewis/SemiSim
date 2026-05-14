# 0017. Adaptive timestep for the transient runner

- Status: Accepted
- Date: 2026-05-09
- Milestone: M18

## Context

[ADR 0010](0010-bdf-time-integration.md) defined the BDF1 / BDF2
time integration in the transient runner shipped in M13 and stated
explicitly:

> Adaptive time stepping is not implemented in M13. Fixed dt is
> used throughout.

That choice carried a CI cost. The
`examples/power_diode_reverse_recovery/` example added in PR #85
exercises a `voltage_t.table` waveform with a +0.7 V to -2.0 V
transition at t = 50-60 ns. With fixed dt = 1 ns the BDF2 step at
the ramp boundary mixes a forward-conduction history with a sharp
reverse-bias BC update; SNES stagnates somewhere in the
transition or early-recovery window. The CI workflow comment
introduced in `ed6719b` recorded the carve-out:

> Both stagnate in regimes the current fixed-step SNES driver
> cannot handle: a V_GS bias-sweep through the MOSFET inversion
> onset under Fermi-Dirac statistics, and a forward-conduction-
> to-reverse-blocking transient on a long-base diode with stored
> minority carriers. ... the fix is runner work (adaptive dt for
> the transient runner; SNES line-search / stabilization upgrade
> for the bias_sweep continuation across the MOSFET threshold).
> Tag both as allow-failure: "true" so the catalogue matrix still
> runs end-to-end while the underlying solver follow-ups land.
> Drop the flag once the runner changes ... are in place.

The transient half of that carve-out is what this ADR resolves;
the bias_sweep half (`nmos_idvgs`) stays carved-out until a
separate stabilization follow-up lands.

[ADR 0014](0014-slotboom-transient.md) is the load-bearing
formulation for the transient block. It is unchanged by M18; the
adaptive controller drives the time loop, not the residual.

[`semi/continuation.py`](../../semi/continuation.py) already
provides an `AdaptiveStepController` with grow-on-easy-solves /
halve-on-failure semantics, used by the bias_sweep voltage-ramp.
The controller is signed-step but sign-agnostic in the inner
arithmetic; reusing it on the dt axis (always positive) is a
direct fit.

## Decision

1. **Reuse `AdaptiveStepController` on the dt axis.** No new
   continuation module, no parallel implementation. The
   bias_sweep ramp's controller is the load-bearing user; M18
   reuses it without modification, so any future tweak (line-
   search, larger grow factor, telemetry hooks) is shared.

2. **Add variable-step BDF2 in `semi/timestepping.py`.** A new
   classmethod `BDFCoefficients.variable_bdf2(omega)` returns the
   non-uniform-stencil triplet
   `((1+2*omega)/(1+omega), -(1+omega), omega**2/(1+omega))`
   for `omega = dt_n / dt_{n-1}`. At `omega == 1.0` the triplet
   collapses to `(1.5, -2.0, 0.5)` exactly, bit-identical to the
   uniform-step coefficients stored in
   `BDFCoefficients._SUPPORTED_COEFFS[2]`. The uniform-step
   table stays as the authoritative lookup; the variable-step
   classmethod is additive and does not change the existing
   constructor or `apply` method signatures. BDF1 is dt-agnostic
   and needs no variable form.

3. **Default off; opt-in via `solver.adaptive`.** Schema additive
   minor bump v2.8.0 -> v2.9.0 adds a `solver.adaptive` block
   with `enabled` (default `false`), `dt_min`, `dt_max`,
   `easy_iter_threshold` (default 4), `grow_factor` (default
   1.5), and `max_consecutive_failures` (default 6). Configs
   without the block (or with `enabled: false`) route to the
   existing fixed-dt code path, which is preserved verbatim and
   is the bit-identity branch. `dt_min` and `dt_max` are
   required when `enabled: true`; the schema cross-field
   validator rejects `dt_min > dt_max` and `dt < dt_min` or
   `dt > dt_max` so the initial dt sits inside the controller
   range. The block is rejected at validate time on non-transient
   solver types.

4. **Snapshot-and-retry on SNES failure.** The adaptive time
   loop in `semi/runners/transient.py` snapshots the Slotboom
   state (`psi.x.array`, `phi_n.x.array`, `phi_p.x.array`)
   before the BC seed and SNES solve. On failure: restore the
   snapshot, call `controller.on_failure()` (which halves dt or
   raises `StepTooSmall` at the floor), increment a contiguous-
   failure counter, and retry. On success: append converged
   carrier densities to the BDF history, advance `t_current`,
   and call `controller.on_success(n_iter)` if the step was not
   externally clamped. History is mutated only on success, so it
   does not need rollback.

5. **Clamp dt at endpoints and waveform breakpoints.** Two
   external clamp rules apply:

   - **Endpoint clamp 1 (`t_end`).** If `t_current + dt_current >
     t_end`, set `dt_current = t_end - t_current`. The clamp does
     not feed back into the controller (the next step would be
     past `t_end` anyway).
   - **Endpoint clamp 2 (waveform breakpoints).** For each
     ohmic contact's `voltage_t`, the `step` variant has a
     discontinuity at `t0` and the `table` variant has a slope
     change at every interior `times[i]`. The runner pre-
     computes the sorted list of breakpoints with `t > 0` once
     at setup, and on each step clamps `dt_current` to land
     exactly on the next breakpoint inside `(t_current,
     t_current + dt_current]`. The clamp does not feed back
     into the controller; `dt_prev` is still updated after a
     successful clamped step so the next variable-step BDF2
     `omega` is correct.

   Both clamps are required for `power_diode_reverse_recovery`:
   the `+0.7 V to -2.0 V` transition at t = 50-60 ns has 11
   interior `times[i]` samples that the integrator must land on
   exactly, and the t_end clamp ensures the run terminates at
   the documented 200 ns.

6. **First step is always BDF1.** BDF2 needs two history points;
   the seeding step is BDF1 regardless of `solver.order`. Once
   `len(n_hist) >= 2`, BDF2 takes over with
   `variable_bdf2(omega)` per step. `omega == 1.0` is the bit-
   identity hinge.

## Consequences

- **Bit-identity gate on every existing benchmark.** The
  fixed-dt branch is preserved verbatim under
  `if not adaptive_enabled:`; configs without `solver.adaptive`
  produce byte-identical results to v0.24.0 on
  `pn_1d_turnon`, `pn_1d_pulse`, `diode_sine_1d`,
  `rc_ac_sweep`, and every steady-state benchmark whose
  `solver.type` is not `"transient"`.
- **Audit case 07** is the V&V gate. M18 ships no MMS variant,
  on the same precedent as M16.7: per
  [ADR 0006](0006-verification-and-validation-strategy.md) the
  per-physics-module MMS rule applies to new physics kernels;
  M18 ships none. The transient runner's existing variant
  ladder (Variants A through H) gates the residual stack;
  audit case 07
  (`tests/audit/test_07_adaptive_dt_vs_fixed_dt.py`) gates the
  controller against the shipped fixed-dt reference within 1 %
  on `benchmarks/pn_1d_turnon`.
- **`power_diode_reverse_recovery` retires its
  `allow-failure: "true"` flag.** The example enables
  `solver.adaptive` (`dt_min = 1 ps`, `dt_max = 5 ns`); the
  controller halves dt across the +0.7 V to -2.0 V transition
  and grows back during the settled reverse-blocking tail. The
  `mosfet_2d` and `nmos_idvgs` `allow-failure: "true"` flags
  stay in place; both are bias_sweep SNES stabilization work
  and are not in M18's scope.
- **ADR 0010 amended.** The line "Adaptive time stepping is not
  implemented in M13. Fixed dt is used throughout." is now
  scoped: M13 shipped with fixed dt; M18 added the optional
  adaptive path on top. ADR 0010's BDF derivations and the
  uniform-step coefficients remain authoritative; the variable-
  step classmethod is additive and bit-identical at `omega = 1`.
- **Coverage gate at 95 holds without a follow-up commit.** The
  pure-Python tests in `tests/test_adaptive_dt_schema.py` and
  `tests/test_variable_bdf2.py` plus the gated FEM tests in
  `tests/fem/test_adaptive_dt_transient.py` cover every
  controller transition. M18 avoids the M16.5 / M16.6 follow-up-
  commit pattern in which new branches were exercised only by
  benchmarks (whose coverage is not merged into the gated job).
- **Bias-sweep follow-up.** The remaining half of the CI
  carve-out is the V_GS sweep across the MOSFET inversion onset
  under Fermi-Dirac statistics, which stagnates because the SNES
  default line search (`bt`, backtracking) cannot find a descent
  direction at the threshold. ADR-level decision; reuses the
  same `AdaptiveStepController` along the V_GS axis (already in
  place for forward-bias ramps).

## Alternatives considered

- **A new continuation module dedicated to time.** Rejected.
  `AdaptiveStepController` is sign-agnostic in the inner
  arithmetic; the dt axis just uses positive steps. Forking the
  controller would split future maintenance and risk drift
  between the bias-sweep and transient paths.
- **PETSc TSAdapt for adaptive dt.** Rejected. The transient
  runner is a hand-driven SNES loop, not a PETSc TS object, by
  design (ADR 0014). Wiring TSAdapt would require a TS
  reformulation that the M13 / M13.1 architecture choice
  explicitly avoided. The hand-driven loop is the right level
  for the controller to attach to.
- **Embedded error estimator (Milne / Petzold) instead of
  grow-on-easy-solves.** Rejected for now. The bias_sweep
  controller's heuristic (count consecutive easy SNES solves;
  grow by `grow_factor`) is shared between the two axes and
  makes both paths reason the same way. Local-truncation-error
  estimators are a candidate future enhancement when there is a
  benchmark that exercises a stiff-region adaptation pattern
  the heuristic cannot handle. None exists today.

## References

- [`semi/continuation.py`](../../semi/continuation.py)
  (`AdaptiveStepController`, `StepTooSmall`).
- [`semi/timestepping.py`](../../semi/timestepping.py)
  (`BDFCoefficients`, `variable_bdf2`).
- [`semi/runners/transient.py`](../../semi/runners/transient.py)
  (adaptive time loop; the bit-identity fixed-dt branch is
  preserved verbatim).
- [`schemas/input.v2.json`](../../schemas/input.v2.json)
  (`solver.adaptive` block; v2.9.0).
- [`tests/audit/test_07_adaptive_dt_vs_fixed_dt.py`](../../tests/audit/test_07_adaptive_dt_vs_fixed_dt.py)
  (V&V gate).
- [ADR 0010](0010-bdf-time-integration.md). M18 amends ADR 0010
  by lifting the "no adaptive dt" statement, scoped to the
  transient runner.
- [ADR 0014](0014-slotboom-transient.md). Unchanged by M18.
- [ADR 0006](0006-verification-and-validation-strategy.md).
  Same V&V scope reasoning as M16.7 (audit-only; no MMS
  variant).
