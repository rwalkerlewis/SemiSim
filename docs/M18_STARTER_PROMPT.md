## Context

You are working in `kronos-semi` at `v0.24.0` (post-M17, package
version `0.24.0`, schema `2.8.0`). The repo is at
`https://github.com/rwalkerlewis/kronos-semi`; `main` is at
`ed6719b` (`fix(ci): unblock main by reading cathode J in
schottky verifier; tag two known-broken examples allow-failure`).

Your assignment is **M18: Adaptive timestep for the transient
runner**. The M16 umbrella (M16.1 through M16.7) and M17
(heterojunctions) closed the position-dependent-physics expansion.
M18 is a runner-quality follow-up that is not on the original
roadmap as a numbered milestone; it is named in the CI carve-out
comment in `.github/workflows/ci.yml` as the exit criterion for
the `power_diode_reverse_recovery` example's `allow-failure: "true"`
tag.

The CI commit at HEAD names the carve-out and its exit criterion
verbatim:

> `nmos_idvgs` and `power_diode_reverse_recovery` (allow-failure
> holding pattern). Both stagnate in regimes the current
> fixed-step SNES driver cannot handle: a V_GS bias-sweep through
> the MOSFET inversion onset under Fermi-Dirac statistics, and a
> forward-conduction-to-reverse-blocking transient on a long-base
> diode with stored minority carriers. Config tuning (dt, BDF
> order, ramp duration, jacobian_shift, snes tolerances,
> voltage_t shape, bias swing) does not unstick either; the fix
> is runner work (adaptive dt for the transient runner; SNES
> line-search / stabilization upgrade for the bias_sweep
> continuation across the MOSFET threshold). Tag both as
> allow-failure: "true" so the catalogue matrix still runs
> end-to-end while the underlying solver follow-ups land. Drop
> the flag once the runner changes (adaptive dt for transient;
> SNES stabilization for bias_sweep across MOSFET threshold) are
> in place.

M18 closes the **transient half** of that carve-out and only that
half. The bias-sweep stabilization that unblocks `nmos_idvgs` is a
separate follow-up tracked as a candidate for the next milestone
after M18; this PR does not touch `nmos_idvgs` or `mosfet_2d`. The
`mosfet_2d` allow-failure flag stays in place.

This prompt **is** the starter. Save it as
`docs/M18_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below).

## Branch and PR rules

- Work on a fresh branch off `main`: `git checkout -b dev/m18-adaptive-dt-transient`.
- One milestone, one PR.
- Push every phase commit immediately. Open the PR after Phase 0.
- **Run all six phases consecutively without pausing for
  confirmation, status reports, or "should I continue" prompts.**
  Do not call `ask_user_input` or any equivalent between phases.
  The phases are an execution plan, not a checklist of
  separately-permissioned tasks. Run Phase 0 through Phase F end
  to end, push commits as they land, and open the PR after Phase 0
  so reviewers can watch the remaining phases land.

  The only legitimate reasons to stop and surface a question are:
  (a) a hard blocker that cannot be resolved by re-reading this
  prompt, the linked ADRs, or the codebase; (b) a `pytest` or
  `ruff` failure that persists after a reasonable repair attempt;
  or (c) the Stop conditions section reports green and the PR is
  ready for review. "Should I move on to the next phase" and "do
  you want me to also add X" are not legitimate reasons. Default
  behavior: pick the most defensible interpretation, document it
  inline, and continue.

- Title the PR `M18: Adaptive timestep for the transient runner`.
- **No AI-assistant credits in any shipped artifact.** Verify
  with `git log --format=%B dev/m18-adaptive-dt-transient | grep
  -i 'co-authored-by: claude'` before opening the PR. Empty
  output is the gate.

## Required reading (in order; ~45 minutes)

1. `PLAN.md` in full. Confirm `main` is at v0.24.0 (post-M17) and
   that no other PR is in flight on `dev/m18-*`. Read the
   "Next task" block: it currently names M19 (3D MOSFET capstone).
   M18 inserts ahead of M19. After M18 lands, the next reviewer
   either keeps M19 or promotes the bias-sweep stabilization
   follow-up that unblocks `nmos_idvgs`.
2. `docs/IMPROVEMENT_GUIDE.md` § 1 (Honest current state) and
   § 4 (the gap list near `examples/`). M18 is not yet enumerated
   in § 4; Phase F adds it as a Done entry.
3. `docs/ROADMAP.md` § Capability matrix.
4. `docs/adr/0010-bdf-time-integration.md`. Note the explicit
   line: "Adaptive time stepping is not implemented in M13. Fixed
   dt is used throughout." M18 retires that statement.
5. `docs/adr/0014-slotboom-transient.md`. The Slotboom transient
   formulation is unchanged by M18; only the time-loop driver
   changes.
6. **`semi/runners/transient.py` end to end.** The relevant
   locations:
   - L161-180 the solver-config parse block (`t_end`, `dt_val`,
     `order`, `max_steps`, `snes`).
   - L342-368 the dimensionless `dt_const` Constant and the
     residual builder call.
   - L409-416 history storage (`n_hist`, `p_hist`, `_trim_history`).
   - L464-540 the time loop body (BDF order selection, history
     update, BC build, SNES solve, history append).
   - L766-893 `_run_bc_continuation` (the BC-ramp; ADR 0013).
     M18 does not change this.
   - L918-960 `_build_voltage_t_evaluator` (M16.7). M18 reuses this
     but adds breakpoint-aware dt clamping.
7. **`semi/continuation.py` end to end.** `AdaptiveStepController`
   is the controller M18 reuses for dt scheduling. It is signed-
   step grow/halve with `easy_iter_threshold`, `grow_factor`,
   `min_step_abs`, `max_step_abs`, and `StepTooSmall` for the
   floor exhaustion case.
8. `semi/runners/bias_sweep.py` § "Snapshot/restore" block (around
   L260-273) and the `ramp_leg` adaptive controller usage (L352-422).
   This is the model. The transient version snapshots
   `(psi, phi_n, phi_p)` plus the carrier-density history
   `n_hist`, `p_hist` and the `t_current` cursor.
9. `semi/timestepping.py` end to end. The `BDFCoefficients` class
   currently hard-codes uniform-step BDF1 and BDF2 coefficients in
   `_SUPPORTED_COEFFS`. M18 adds a variable-step BDF2 path; uniform
   BDF2 must remain bit-identical when the dt ratio is exactly 1.0.
10. `examples/power_diode_reverse_recovery/diode_recovery.json`
    and `examples/power_diode_reverse_recovery/README.md`. The
    failing case M18 unblocks. The `voltage_t` table has a
    forward-conduction segment at +0.7 V (t in [0, 50] ns), a
    10 ns linear ramp from +0.7 V to -2.0 V (t in [50, 60] ns),
    and reverse blocking at -2.0 V (t in [60, 200] ns). The
    fixed-dt run stagnates somewhere in the transition or
    early-recovery window; the adaptive controller has to halve
    dt across the transition and grow it back during the settled
    reverse-blocking tail.
11. `tests/test_continuation.py` end to end. The test patterns
    M18's new dt-controller tests follow.
12. `.github/workflows/ci.yml` lines 220-243. The
    `allow-failure: "true"` carve-out for
    `power_diode_reverse_recovery` is what M18 retires. The
    `nmos_idvgs` and `mosfet_2d` flags stay.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** New `solver.transient.adaptive` block
  with subfields `enabled`, `dt_min`, `dt_max`, `easy_iter_threshold`,
  `grow_factor`, `max_consecutive_failures`. All optional under
  `additionalProperties: false`; the absence of the block is the
  default and is bit-identical to v0.24.0 fixed-dt.
- **Schema versioning is binding.** Additive minor bump v2.8.0 to
  v2.9.0. v2.0.0 through v2.8.0 inputs must continue to validate
  and produce bit-identical results to v0.24.0. Update
  `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Bit-identity gate.** Every existing benchmark runs bit-
  identical to v0.24.0. The default for `solver.transient.adaptive`
  is unset, which routes to the existing fixed-dt code path
  unchanged. Reuse all prior anchors:
  `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2`,
  `diode_velsat_1d` 56.27 % at 0.9 V and 0.19 % at 0.3 V,
  `diode_auger_1d` >20 % divergence at 0.9 V,
  `diode_fermi_dirac_1d` 7.37 % FD-vs-Boltzmann V_bi,
  `schottky_1d` worst-case <10 % thermionic match,
  `zener_1d` worst-case <20 % Kane match,
  `pn_1d_turnon` (M13 reference transient),
  `pn_1d_pulse` and `diode_sine_1d` (M16.7 demonstrators),
  `rc_ac_sweep` AC sweep,
  `algaas_gaas_heterojunction_1d` (M17 reference, if landed in
  benchmarks/; otherwise the tagged audit suffices).
- **Five layers, enforced.** The dt-controller adapter lives in
  `semi/runners/transient.py` (Layer 5). The schema block lives in
  `schemas/input.v2.json` (Layer 1) and the loader's conditional
  validation in `semi/schema.py` (Layer 2). The pure-Python
  `AdaptiveStepController` in `semi/continuation.py` (Layer 3) is
  reused without modification. No physics module is touched. No
  FEM types leak.
- **MMS scope.** Do **not** add an MMS Variant for M18. Same
  precedent as M16.7: the per-physics-module MMS rule from
  ADR 0006 applies to *new physics kernels*; M18 ships no new
  kernel. The transient runner already has MMS coverage at
  Variants A through H. M18 is a time-stepping driver change. A
  short note in `docs/PHYSICS.md` § Verification & Validation is
  sufficient (Phase F).
- **Per-runner threading.** Only the transient runner consumes
  `solver.transient.adaptive`. The schema declares the field at
  the `solver.transient` block; bias_sweep, ac_sweep, equilibrium,
  mos_cv, mos_cap_ac, and resistor_3d ignore it. The transient
  runner reads the block from `cfg["solver"]` (the existing parse
  block at L161-180); add the new keys there.
- **One milestone, one PR.** Do not bundle with M19 (3D MOSFET
  capstone), the bias-sweep SNES stabilization that unblocks
  `nmos_idvgs`, or any M14.2.x backlog item. The
  `mosfet_2d` `allow-failure: "true"` flag stays in place. The
  `nmos_idvgs` `allow-failure: "true"` flag stays in place.
  Only `power_diode_reverse_recovery`'s flag retires.
- **No em dashes in prose or code comments.**
- **No AI-assistant credits.**
- **Coverage gate is 95 in the gated `docker-fem-tests` job.**
  M16.5 and M16.6 each needed coverage-recovery follow-up commits
  because new branches were only exercised by benchmarks (whose
  coverage is not merged). Plan unit tests for the dt controller
  and the variable-step BDF2 path **in Phase E before** the
  example reactivation in Phase D. Specifically: every dt-
  controller transition (grow on N consecutive easy solves, halve
  on SNES non-convergence, floor exhaustion via `StepTooSmall`,
  endpoint clamp to `t_end`, breakpoint clamp at `voltage_t.step.t0`,
  variable-step BDF2 coefficient computation against a known dt
  ratio) gets at least one pure-Python test independent of the
  audit / benchmark matrix.

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/`
and `pytest tests/`. Do not advance with red tests. **Do not
pause between phases for confirmation or status reports.** Phase
gates are mechanical: green tests, push, advance to the next
phase. Run all six phases end to end.

---

### Phase 0: ship this starter prompt

Pure docs.

1. Save this entire file as `docs/M18_STARTER_PROMPT.md`. Strip
   nothing.
2. Append a one-line "Author M18 starter prompt" entry to
   `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under
   `[Unreleased]`.
3. Push the branch and open the PR. Verify the squash-body
   preview has no AI-credit trailer.

**Commit message:** `docs: ship M18 starter prompt (M18)`

**Acceptance:** the file exists at `docs/M18_STARTER_PROMPT.md`
on `dev/m18-adaptive-dt-transient`; CI green on docs-only diff.

---

### Phase A: schema surface for `solver.transient.adaptive`

Pure-Python. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.8.0 to 2.9.0. Bump
   `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`. Add a row to
   the `examples` array.
2. Add `solver.transient.adaptive` (sibling of `dt`, `t_end`,
   `order`, `max_steps`, `output_every`, `bc_ramp_steps`, `snes`):
   - Type `object`, `additionalProperties: false`.
   - `enabled`: boolean. Default `false`. When `false` (or the
     whole `adaptive` block is absent), the runner uses the
     existing fixed-dt path; the rest of the keys are ignored.
   - `dt_min`: number, units seconds. Required when
     `enabled: true`. Lower bound on the controller's dt; falling
     below this raises `StepTooSmall` and the runner aborts the
     step with `RuntimeError`.
   - `dt_max`: number, units seconds. Required when
     `enabled: true`. Upper bound on growth.
   - `easy_iter_threshold`: integer, default 4. SNES iteration
     count at or below which a solve is considered "easy" for the
     purpose of grow scheduling. Matches the bias-sweep convention
     in `semi/runners/bias_sweep.py` L368.
   - `grow_factor`: number, default 1.5. Multiplier applied to dt
     after `easy_iter_threshold` consecutive easy solves. Must
     be strictly greater than 1.0 (rejected by the
     `AdaptiveStepController` constructor otherwise).
   - `max_consecutive_failures`: integer, default 6. Maximum
     contiguous halvings before the runner raises `RuntimeError`.
     The floor-vs-`dt_min` check raises `StepTooSmall` first if
     dt would fall below `dt_min`; this counter catches the case
     where dt halves repeatedly while staying above `dt_min`.
   - Description: "Adaptive time-step controller for the transient
     runner. When `enabled: true`, dt is set by an
     `AdaptiveStepController` (`semi/continuation.py`) that grows
     on `easy_iter_threshold` consecutive easy SNES solves and
     halves on SNES non-convergence. The initial dt at the start
     of the time loop is `solver.transient.dt`. When `enabled:
     false` (or the block is absent), the runner uses the existing
     fixed-dt path. (M18)"
3. Conditional validation in `semi/schema.py`:
   - When `solver.transient.adaptive.enabled: true`, require
     `dt_min` and `dt_max`. Both must be positive; `dt_min <=
     dt <= dt_max` (the initial dt sits in the controller range).
   - `easy_iter_threshold >= 1`.
   - `grow_factor > 1.0` (the `AdaptiveStepController` constructor
     also enforces this; the schema rejection is a clearer error
     for users).
   - `max_consecutive_failures >= 1`.
   - For non-transient solver types (`solver.type` not equal to
     `"transient"`), reject configs where `solver.transient.adaptive`
     is set. Match the M16.5 conditional-validation pattern.
4. Default-fill: a config with no `solver.transient.adaptive` is
   bit-identical to v0.24.0 on every existing benchmark.
5. Pure-Python tests in a new
   `tests/test_adaptive_dt_schema.py`:
   - Every existing benchmark JSON validates unchanged against
     v2.9.0.
   - `solver.transient.adaptive` with all fields validates.
   - `enabled: true` without `dt_min` is rejected.
   - `enabled: true` without `dt_max` is rejected.
   - `dt_min > dt_max` is rejected.
   - `dt < dt_min` or `dt > dt_max` is rejected.
   - `grow_factor: 1.0` is rejected.
   - `solver.transient.adaptive` set on `solver.type:
     "bias_sweep"` is rejected.
6. Update `docs/schema/reference.md` versioning table with the
   v2.9.0 row; the change-log column reads "additive: adaptive
   time-step controller `solver.transient.adaptive` for the
   transient runner (M18)".

**Commit message:** `feat(schema): solver.transient.adaptive controller block (schema 2.9.0); existing solvers reject (M18)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change.

---

### Phase B: variable-step BDF2 in `semi/timestepping.py`

The current `BDFCoefficients` hard-codes uniform-step coefficients
`(3/2, -2, 1/2)` for BDF2. With variable dt, BDF2 needs the
non-uniform formula. When the ratio is exactly 1.0, the variable
formula must collapse to the uniform one (bit-identity).

1. Add a `BDFCoefficients.variable_bdf2(omega: float)` classmethod
   (or an instance method) that returns the three coefficients
   `(alpha_0, alpha_1, alpha_2)` for variable-step BDF2 at ratio
   `omega = dt_n / dt_{n-1}`:
   ```
   alpha_0 = (1 + 2*omega) / (1 + omega)
   alpha_1 = -(1 + omega)
   alpha_2 = omega**2 / (1 + omega)
   ```
   These reduce to `(3/2, -2, 1/2)` when `omega == 1.0`.
2. The runner consumes these coefficients per step; the
   `_SUPPORTED_COEFFS` table stays as the authoritative uniform-
   step lookup. The `variable_bdf2` method is additive and does
   not change the existing constructor or `apply` method
   signatures.
3. BDF1 is dt-agnostic (the formula `(u^{n+1} - u^n) / dt` has no
   ratio); no `variable_bdf1` is needed.
4. Pure-Python tests in `tests/test_variable_bdf2.py`:
   - `omega = 1.0` returns `(1.5, -2.0, 0.5)` exactly.
   - `omega = 2.0` returns `(5/3, -3, 4/3)` (a known
     hand-computed triplet).
   - `omega = 0.5` returns `(4/3, -1.5, 1/6)`.
   - Coefficient sum identity: `alpha_0 + alpha_1 + alpha_2 == 0`
     for all positive `omega` (consistency: the BDF2 stencil is
     exact for constants).
   - First-moment identity: `alpha_0 - alpha_2 * omega ==
     (1 + omega)` (consistency: exact for linear in time at the
     non-uniform stencil; this is the standard non-uniform-BDF2
     derivation check).
   - `omega <= 0.0` raises `ValueError`.
5. Update the `BDFCoefficients` module docstring to mention the
   variable-step path.

**Commit message:** `feat(timestepping): variable-step BDF2 coefficients (M18)`

**Acceptance:** new tests pass; existing tests pass; uniform-step
BDF2 unchanged in `_SUPPORTED_COEFFS`.

---

### Phase C: adaptive-dt controller in the transient runner

Layer 5. The change is local to `semi/runners/transient.py`. The
`AdaptiveStepController` from `semi/continuation.py` is reused
without modification.

1. At runner setup (after the existing `solver_cfg` parse at
   L161-180), parse the new `adaptive` block. If absent or
   `enabled: false`, set a sentinel `adaptive_enabled = False`
   and continue with the existing fixed-dt path. This branch is
   the bit-identity branch.
2. When `adaptive_enabled = True`:
   - Construct an `AdaptiveStepController(initial_step=dt_val,
     max_step_abs=dt_max, min_step_abs=dt_min,
     easy_iter_threshold=easy_iter_threshold,
     grow_factor=grow_factor)`. The "step" the controller manages
     is dt (always positive; sign is fixed because time advances
     forward).
   - Track a separate `dt_prev` cursor (initially `None`); the
     first step uses BDF1 regardless of `order` (BDF2 needs two
     history points; the seeding step is always BDF1). Once
     `len(n_hist) >= 2`, switch to BDF2 with `omega = dt_current /
     dt_prev` via `BDFCoefficients.variable_bdf2(omega)`. When
     `omega == 1.0` the coefficients are bit-identical to the
     uniform-BDF2 path used today.
3. Replace the existing time-loop body (L464-540) with an
   adaptive-aware version. The structure mirrors `bias_sweep`'s
   `ramp_leg`:
   - Snapshot `psi.x.array`, `phi_n.x.array`, `phi_p.x.array`,
     and the *length* of `n_hist` and `p_hist` (so on rejection,
     append-then-pop is well-defined). Also snapshot
     `t_current`, `dt_prev`, the controller state implicitly via
     `controller.step` and `controller._easy_count` (the
     controller is mutable; pair the snapshot with a paired
     `restore` that resets these).
   - Compute `dt_current = controller.step` and clamp via the
     two endpoint rules below.
   - **Endpoint clamp 1: t_end.** If `t_current + dt_current >
     t_end`, set `dt_current = t_end - t_current`. This clamp does
     not feed back into the controller (the next step would be
     past t_end anyway); record it as a special "endpoint" step.
   - **Endpoint clamp 2: voltage_t breakpoints.** For every
     contact with `voltage_t`, find the next discontinuity or
     slope change at `t > t_current`:
     - For `voltage_t.type == "step"`: the discontinuity is at
       `t0`. If `t_current < t0 <= t_current + dt_current`, clamp
       `dt_current = t0 - t_current` (the step lands exactly on
       `t0`; the BC is then `v1` at the new step per the M16.7
       evaluator).
     - For `voltage_t.type == "table"`: the slope changes at every
       interior `times[i]`. If any `times[i]` in
       `(t_current, t_current + dt_current]` exists, clamp
       `dt_current = min(times[i]) - t_current`. This prevents
       the integrator from straddling a slope change in the
       waveform, which is the principal cause of the
       `power_diode_reverse_recovery` stagnation.
     - When a clamp triggers, do not call
       `controller.on_success` after the solve (the dt was
       externally constrained, not the controller's choice). Do
       still update `dt_prev` after a successful clamped step.
4. Solve at the (possibly clamped) `dt_current`:
   - If success: append `n_new`, `p_new` to history,
     `_trim_history()`, advance `t_current`, set `dt_prev =
     dt_current`, call `controller.on_success(n_iter)` (unless
     the step was clamped per above). Record IV. Snapshot
     persistence: the snapshot goes out of scope, no restore
     needed.
   - If failure: restore from snapshot. Call
     `controller.on_failure()`; if it raises `StepTooSmall`,
     re-raise as `RuntimeError(f"Adaptive transient stalled at
     t={t_current:.3e} s (dt_min={dt_min:.3e} reached)")`.
     Increment a contiguous-failure counter; if it exceeds
     `max_consecutive_failures`, raise the same
     `RuntimeError`. Do not increment `step_count`. Restart the
     while-loop iteration with the new (smaller) dt.
   - Reset the contiguous-failure counter on every success.
5. The fixed-dt code path is preserved verbatim under
   `if not adaptive_enabled:`. No diff for the fixed-dt branch
   at the level of solved residuals or BC build order. The bit-
   identity gate is checked by re-running every existing
   transient benchmark in CI.
6. Two notes in the runner module docstring:
   - The `adaptive_enabled = False` branch is bit-identical to
     v0.24.0.
   - When `adaptive_enabled = True`, the first step is always
     BDF1 (history seeding); BDF2 takes over once two history
     points exist, using `variable_bdf2(omega)` per step.
7. New unit tests in `tests/fem/test_adaptive_dt_transient.py`
   (gated, requires dolfinx). Eight tiny configs (1D pn diode,
   8 cells, short t_end) covering the controller transitions:
   - `test_grows_after_easy_solves`: a benign config (steady-
     state-like, no voltage_t) where SNES converges in 2-3
     iterations every step. Asserts dt grows at least once over
     a 20-step window.
   - `test_halves_on_snes_failure`: a config tuned to fail the
     first solve at the initial dt (large initial dt, tight
     tolerances). Asserts dt halves at least once and the run
     completes.
   - `test_floor_exhaustion_raises`: a config tuned to fail
     even after halving to `dt_min` (e.g., a pathological
     reverse-bias jump that the engine cannot integrate).
     Asserts `RuntimeError` with the floor-exhaustion message.
   - `test_endpoint_clamps_to_t_end`: a config where the last
     step would overshoot `t_end`. Asserts the final `t_vals[-1]`
     equals `t_end` to within `1e-12 * t_end` and the controller
     is **not** notified (the endpoint clamp is external).
   - `test_step_breakpoint_clamp`: `voltage_t.type = "step"`
     with `t0 = 5e-9`. Initial dt = 1e-8 (would straddle t0).
     Asserts the first step lands exactly at t0 and the BC is
     `v1` at the post-step solve.
   - `test_table_breakpoint_clamp`: `voltage_t.type = "table"`
     with non-uniform `times`. Initial dt larger than the
     smallest interior `times[i+1] - times[i]`. Asserts the
     first step lands exactly at the next interior `times[i]`.
   - `test_bit_identity_when_disabled`: same config run with
     `adaptive: {enabled: false}` and with `adaptive` absent.
     Asserts byte-identity of the IV trace.
   - `test_variable_bdf2_at_dt_change`: a config where dt halves
     once and then grows back. Asserts the recorded coefficients
     at the post-grow step match `variable_bdf2(omega)` for the
     observed `omega`. (Add a debug hook on the runner exposing
     `(alpha_0, alpha_1, alpha_2)` per step under a config flag,
     or capture via `progress_callback`.)

**Commit message:** `feat(runners): adaptive-dt controller for transient.py with breakpoint-aware step clamping (M18)`

**Acceptance:** new tests pass; existing transient benchmarks
bit-identical to v0.24.0 (the adaptive block is opt-in).

---

### Phase D: turn on adaptive in `power_diode_reverse_recovery`

The headline deliverable. The example currently has
`allow-failure: "true"` in `.github/workflows/ci.yml` because the
fixed-dt run stagnates across the `voltage_t` ramp transition.
Phase D enables adaptive on the example and retires the flag.

1. Edit `examples/power_diode_reverse_recovery/diode_recovery.json`:
   - Add `solver.transient.adaptive`:
     ```
     "adaptive": {
       "enabled": true,
       "dt_min": 1.0e-12,
       "dt_max": 5.0e-9,
       "easy_iter_threshold": 4,
       "grow_factor": 1.5,
       "max_consecutive_failures": 6
     }
     ```
   - Bump the config's `schema_version` to `"2.9.0"`.
   - Optionally raise `solver.max_steps` from 250 to 2000 to
     accommodate the smaller dt during the transition (the
     adaptive controller will halve dt as small as 1 ps in the
     hardest neighborhood; the existing 250-step cap is sized for
     fixed dt = 1 ns).
   - Leave `solver.dt = 1.0e-9` as the initial dt (the controller
     starts there and adapts).
2. Update `examples/power_diode_reverse_recovery/README.md`:
   - Add a "How adaptive dt is used" subsection explaining the
     `adaptive` block, the breakpoint clamping at the two slope
     changes (t = 50 ns and t = 60 ns) plus every interior
     `voltage_t.times[i]` sample, and the expected halving behavior
     across the +0.7 V to -2.0 V transition.
   - Cite the M18 milestone and link `docs/M18_STARTER_PROMPT.md`.
3. Edit `.github/workflows/ci.yml`:
   - Remove the `allow-failure: "true"` flag from the
     `power_diode_reverse_recovery` matrix entry (lines 241-243
     in the current file).
   - Update the carve-out comment block (lines 220-235) to remove
     the `power_diode_reverse_recovery` reference and to leave the
     `nmos_idvgs` carve-out reasoning intact. The new comment says:
     ```
     # nmos_idvgs is the remaining allow-failure carve-out:
     # exercises a V_GS bias sweep across the MOSFET inversion
     # onset under FD statistics, where the current SNES
     # configuration in bias_sweep stagnates without a line-
     # search-strategy upgrade. This is bias_sweep work, not
     # transient work; M18 (adaptive dt for the transient
     # runner) does not address it. Tagged allow-failure: "true"
     # so the catalogue matrix is still exercised end-to-end
     # while the bias_sweep stabilization follow-up lands.
     ```
4. Verify the verifier (`scripts/run_benchmark.py
   verify_power_diode_reverse_recovery`, registered at L3643)
   passes with the new adaptive config:
   - `iv` has at least 20 anode samples.
   - `J(t)` finite everywhere (no NaN, no Inf).
   - Forward-conduction-window samples are positive.
   - Reverse-recovery-window samples have a negative minimum.
   - Final-time current is closer to zero than the recovery
     minimum (settling toward reverse blocking).
   The verifier's qualitative gates are unchanged from
   v0.23.x; only the upstream solve is now adaptive.
5. Document the observed dt trajectory in the example's
   README.md "Expected output" section (e.g., "the controller
   halves dt to ~50 ps across the t = 50 ns to t = 60 ns
   transition, then grows back toward `dt_max` during the settled
   reverse-blocking tail"). This is informational, not a gate.

**Commit message:** `feat(examples): power_diode_reverse_recovery uses adaptive dt; retire CI allow-failure (M18)`

**Acceptance:** `power_diode_reverse_recovery` runs end to end on
CI without `allow-failure: "true"`. The verifier exits 0.

---

### Phase E: regression bench and audit

The bit-identity gate is the headline V&V for M18. Add an
explicit regression bench that exercises the adaptive controller
on a benchmark with a known fixed-dt reference, and an audit case
that watches the dt trajectory shape.

1. New `tests/audit/test_07_adaptive_dt_vs_fixed_dt.py`. Audit
   case 07. Pattern matches `test_06_transient_fft_vs_ac_sweep.py`:
   - Load `benchmarks/pn_1d_turnon`. Run twice:
     (a) the shipped fixed-dt config (the v0.24.0 reference);
     (b) the same config with `solver.transient.adaptive.enabled =
     true, dt_min = 1e-13, dt_max = solver.dt` (i.e., the
     adaptive run is constrained to never exceed the fixed dt;
     this isolates the controller's halving response without
     comparing against a different-resolution trace).
   - Compare the IV traces at every recorded `t` value.
     Acceptance: the adaptive trace matches the fixed-dt trace
     within 1 % relative on the recorded I(t) at every
     comparison point. (The adaptive run may have more recorded
     points than the fixed-dt run if it halved during the
     transient; subsample to the fixed-dt grid via linear
     interpolation before comparing.)
   - Record CSV / markdown artifacts via `tests/audit/_helpers.py`.
2. Update `docs/PHYSICS_AUDIT.md` with audit case 07 entry.
3. Coverage: the unit tests in
   `tests/fem/test_adaptive_dt_transient.py` (Phase C step 7) and
   the schema tests in `tests/test_adaptive_dt_schema.py` (Phase
   A step 5) and the variable-step BDF2 tests in
   `tests/test_variable_bdf2.py` (Phase B step 4) cover the new
   branches in the gated `docker-fem-tests` job. Verify the
   coverage gate at 95 holds without a follow-up commit.

**Commit message:** `test(audit): activate case 07 (adaptive dt bit-bound vs fixed dt within 1 %) (M18)`

**Acceptance:** audit case 07 passes; coverage gate at 95 holds.

---

### Phase F: closeout

1. `PLAN.md`:
   - Move M18 from anywhere it currently appears (it does not, but
     once Phase 0 lands the starter prompt is in `docs/`) to
     "Completed work log" with PR number, deliverables, schema
     minor bump (2.8.0 to 2.9.0), and acceptance-test results
     (cite observed worst-case `|I_adaptive - I_fixed| / |I_fixed|`
     from audit case 07).
   - Update "Next task" to either **M19 (3D MOSFET capstone)**
     or **bias-sweep SNES line-search stabilization (unblocks
     `nmos_idvgs`)**. The maintainer chooses; do not author a
     starter prompt for either. Note that M18 unblocks the
     transient half of the CI carve-out and that the bias_sweep
     follow-up is the remaining half.
   - Refresh "Current state" with the new package version (bump
     0.24.0 to 0.25.0).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Add an M18 Done entry in § 4 with the deliverables (schema
     bump, runner extension, audit case 07, retirement of the
     `power_diode_reverse_recovery` allow-failure flag).
   - Update § 1 ("Honest current state") to mention adaptive dt
     as a shipped capability of the transient runner.
   - Append a § 9 changelog entry under `[0.25.0]`. Move the
     existing `[Unreleased]` block into `### Released`.
3. `docs/PHYSICS.md` § Verification & Validation: append a brief
   note that M18 is a runner-driver change with audit-only V&V
   (audit case 07), no MMS variant; same precedent as M16.7.
4. `docs/ROADMAP.md`:
   - Add an M18 row in the capability matrix as shipped.
   - Promote either M19 or the bias-sweep stabilization (whichever
     is next-task) to the top of the Planned column.
5. `CHANGELOG.md`:
   - New `[0.25.0]` entry with the M18 line items: schema 2.9.0
     additive (`solver.transient.adaptive`), variable-step BDF2
     in `semi/timestepping.py`, adaptive-dt controller in
     `semi/runners/transient.py`, breakpoint-aware step clamping
     for `voltage_t.step.t0` and `voltage_t.table.times[i]`,
     audit case 07 active, `power_diode_reverse_recovery`
     allow-failure retired.
   - Update the schema-banner comment at the top to mention
     v2.9.0.
6. `pyproject.toml` and `semi/__init__.py`: bump 0.24.0 to
   0.25.0.
7. ADR addition: `docs/adr/0017-adaptive-transient-dt.md`
   documenting the adaptive controller choice (reuse of
   `AdaptiveStepController`), the variable-step BDF2 path, the
   breakpoint clamping rule, and the relationship to ADR 0010
   ("adaptive time stepping is not implemented in M13"; this ADR
   amends that statement, scoped to the transient runner).
   Status: Accepted. Date: today. Milestone: M18.
8. Push every commit to `origin/dev/m18-adaptive-dt-transient`.
   Verify no AI-credit trailers anywhere
   (`git log --format=%B dev/m18-adaptive-dt-transient | grep -i
   'co-authored-by: claude'` returns empty).

**Commit message:** `docs: close out M18 (PLAN, IMPROVEMENT_GUIDE, PHYSICS, ROADMAP, CHANGELOG, ADR 0017)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by
      this PR.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact.
- [ ] `git log --format=%B dev/m18-adaptive-dt-transient | grep -i
      'co-authored-by: claude'` returns empty before opening and
      merging the PR.
- [ ] Pure-Python core remains dolfinx-free
      (`tests/test_lazy_imports.py` clean).
- [ ] Every existing benchmark is bit-identical to v0.24.0 when
      `solver.transient.adaptive` is unset (the default).
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004).
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `kronos_server` public API.
- [ ] Schema bumped per minor (2.8.0 to 2.9.0); v2.0.0 through
      v2.8.0 inputs still validate.
- [ ] MMS rate gate L2 >= 1.99 / H1 >= 0.99 active for Variants A
      through H (no new variant; M18 ships no new physics module).
- [ ] Coverage gate holds at 95 in the gated `docker-fem-tests`
      job. The adaptive-dt branches are covered by
      `tests/fem/test_adaptive_dt_transient.py`,
      `tests/test_adaptive_dt_schema.py`, and
      `tests/test_variable_bdf2.py`, not just by audit case 07.
- [ ] Audit case 07 active and passing within 1 %.
- [ ] `power_diode_reverse_recovery` runs without
      `allow-failure: "true"`.
- [ ] `mosfet_2d` and `nmos_idvgs` `allow-failure: "true"` flags
      retained, unchanged from v0.24.0.

## Anti-goals

- Do not add an MMS Variant for M18. M18 ships no new physics
  kernel; the per-physics-module MMS rule (ADR 0006) does not
  apply. Audit case 07 is the V&V gate.
- Do not modify `AdaptiveStepController` in `semi/continuation.py`.
  The bias-sweep ramp is a load-bearing user of that class; M18
  reuses it without modification.
- Do not generalize adaptive dt to non-transient runners.
  Bias-sweep already has its own adaptive controller for the
  voltage step. AC sweep, equilibrium, mos_cv, mos_cap_ac, and
  resistor_3d are not transient.
- Do not change anything in `semi/physics/` or `semi/bcs.py`.
  M18 is a runner-driver change plus a schema field; the
  residual builders are unchanged.
- Do not retire the `mosfet_2d` `allow-failure: "true"` flag.
- Do not retire the `nmos_idvgs` `allow-failure: "true"` flag.
  That carve-out is bias-sweep work and is out of scope for M18.
- Do not bundle with M19, the bias-sweep stabilization follow-up,
  or any M14.2.x backlog item.
- Do not add `Co-Authored-By` trailers, generated-by footers, or
  any other AI-credit marker.
- Do not introduce a new continuation module. Reuse the existing
  `semi.continuation.AdaptiveStepController` verbatim.

## Stop conditions

You are done when:

1. The M18 PR is opened on branch `dev/m18-adaptive-dt-transient`.
2. Both acceptance tests in this prompt pass in CI:
   - A1: audit case 07
     (`tests/audit/test_07_adaptive_dt_vs_fixed_dt.py`) active
     and passing within 1 %.
   - A2: every existing benchmark is bit-identical to v0.24.0
     when `solver.transient.adaptive` is unset.
3. Coverage tests
   (`tests/fem/test_adaptive_dt_transient.py` plus
   `tests/test_adaptive_dt_schema.py` plus
   `tests/test_variable_bdf2.py`) pass in the gated
   `docker-fem-tests` job; coverage gate holds at 95 without a
   follow-up commit.
4. `power_diode_reverse_recovery` runs end to end on CI without
   `allow-failure: "true"`. The verifier exits 0.
5. `mosfet_2d` and `nmos_idvgs` `allow-failure: "true"` flags
   are retained, unchanged.
6. PLAN.md, IMPROVEMENT_GUIDE.md, PHYSICS.md, ROADMAP.md,
   CHANGELOG.md, and ADR 0017 reflect the closeout.
7. Package version bumped 0.24.0 to 0.25.0 in `pyproject.toml`
   and `semi/__init__.py`.
8. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI
   assistant. Squash-merge body is also clean.
9. PR reviewed, CI green (modulo the documented `allow-failure`
   on `mosfet_2d` and `nmos_idvgs`, unchanged in scope from
   v0.24.0), merged.

## PR description template

```
## Summary

M18: Adaptive timestep for the transient runner. Closes the
transient half of the CI carve-out introduced in `ed6719b`. The
`power_diode_reverse_recovery` example exits CI with
`allow-failure: "true"` retired.

The adaptive controller reuses
`semi.continuation.AdaptiveStepController` (the same class that
drives the bias-sweep voltage ramp) on the dt axis. Variable-step
BDF2 coefficients are added to `semi/timestepping.py`; uniform-step
BDF2 (omega = 1.0) is bit-identical to v0.24.0. The time loop in
`semi/runners/transient.py` snapshots and restores the Slotboom
state plus the carrier-density history on SNES non-convergence,
halves dt, and retries. dt is clamped at `t_end` and at every
`voltage_t` waveform breakpoint (`step.t0` and every interior
`table.times[i]`) so the integrator never straddles a slope
change.

Schema additive minor bump v2.8.0 to v2.9.0
(`solver.transient.adaptive` block; absent or `enabled: false` is
the default and is bit-identical to v0.24.0). Audit case 07
records the adaptive-vs-fixed bit-bound match within 1 %.

The `mosfet_2d` and `nmos_idvgs` allow-failure flags stay; those
are bias-sweep SNES stabilization work and are the remaining
half of the CI carve-out.

## Acceptance tests

- [ ] A1 (bit-bound match): tests/audit/test_07_adaptive_dt_vs_fixed_dt.py
      adaptive-vs-fixed I(t) within 1 % at every recorded sample
      (observed: <fill in>)
- [ ] A2 (byte-identity): every existing benchmark bit-identical
      to v0.24.0 when `solver.transient.adaptive` is unset
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % at 0.9 V, 0.19 % at 0.3 V;
      diode_auger_1d >20 % divergence at 0.9 V;
      diode_fermi_dirac_1d 7.37 % FD-vs-Boltzmann V_bi at
      N_D = 1e20 cm^-3; schottky_1d worst-case <10 %;
      zener_1d worst-case <20 % Kane match;
      pn_1d_turnon, pn_1d_pulse, diode_sine_1d,
      rc_ac_sweep all bit-identical)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95 (in the gated
      docker-fem-tests job)
- [ ] pytest tests/audit/ (audit case 07 active)
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark power_diode_reverse_recovery
      (no allow-failure)
- [ ] docker compose run --rm benchmark pn_1d_turnon (bit-identical)
- [ ] All other benchmarks bit-identical to v0.24.0 via the
      existing CI matrix.

## Notes

- M18 ships no MMS variant; same precedent as M16.7. Audit case
  07 is the V&V gate. ADR 0006 V&V scope reasoning applies.
- ADR 0017 added documenting the adaptive controller choice and
  amending the "no adaptive dt" statement in ADR 0010.
- The `AdaptiveStepController` in `semi/continuation.py` is
  reused unchanged. Bias-sweep behavior is unaffected.
- Coverage gate at 95 held without a follow-up commit (the
  pattern from M16.5 and M16.6 that needed coverage-recovery
  follow-ups is avoided by Phase E unit tests).
```

## Hand-off

When M18 lands and is reviewed, the transient half of the CI
carve-out is closed. The two next-tier candidates are:

- **Bias-sweep SNES line-search stabilization (unblocks
  `nmos_idvgs`).** The remaining half of the CI carve-out. The
  V_GS sweep across the MOSFET inversion onset under Fermi-Dirac
  statistics stagnates because the SNES default line search
  (`bt`, backtracking) cannot find a descent direction at the
  threshold. Candidate fixes: switch to PETSc `nleqerr` or `cp`
  line search, add a damping schedule, or introduce a homotopy
  parameter on the FD prefactor. ADR-level decision.

- **M19: 3D MOSFET capstone.** Already documented in PLAN.md
  "Next task". Depends on M16.1 (mobility); unblocked.

Either is unblocked. The maintainer chooses; the next starter
prompt is authored against the chosen item in the same shape as
this prompt.
