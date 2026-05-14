## Context

You are working in `SemiSim` at `v0.22.0` (post-M16.6, package
version `0.22.0`, schema `2.6.0`). The repo is at
`https://github.com/rwalkerlewis/SemiSim`; `main` is at
`0db3ef6` (`M16.6: BBT and TAT tunneling (#83)`).

Your assignment is **M16.7: Time-varying transient contact voltage**,
the seventh and final physics-completeness slice of the M16 umbrella.
M16.1 (Caughey-Thomas mobility), M16.2 (Lombardi surface mobility),
M16.3 (Auger), M16.4 (Fermi-Dirac), M16.5 (Schottky), and M16.6 (BBT
and TAT tunneling) have all merged. M16.7 closes the umbrella; after
this PR, the gap list in PLAN.md, IMPROVEMENT_GUIDE.md, and
PHYSICS_INTRO.md reads "no remaining M16 gaps."

`PLAN.md` § "Next task" names M16.7 explicitly:

> **M16.7: transient FFT vs AC sweep validation** on a fresh branch
> `dev/m16.7-transient-ac`. Acceptance tests in
> `docs/IMPROVEMENT_GUIDE.md` § M16.7. M16.7 is the final M16 slice
> and ships no new physics: it adds a V&V cross-check between the
> M13.1 transient runner and the M14 small-signal AC sweep runner.
> The FFT of a transient response at a DC operating point must
> match the AC sweep at that operating point within an agreed
> tolerance.

The PLAN paraphrase under-sells the deliverable. The actual
IMPROVEMENT_GUIDE § M16.7 entry calls for *extending the transient
runner to accept a time-varying contact voltage* — currently the
runner takes a fixed bias for the entire time loop. Without that
extension, the FFT-vs-AC cross-check cannot run. So M16.7 is
*runner extension plus V&V*, not just V&V.

This prompt **is** the starter. Save it as
`docs/M16_7_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below).

## Branch and PR rules

- Work on a fresh branch off `main`: `git checkout -b dev/m16.7-transient-ac`.
- One milestone, one PR.
- Push every phase commit immediately. Open the PR after Phase 0.
- **Run all six phases consecutively without pausing for
  confirmation, status reports, or "should I continue" prompts.**
  Do not call `ask_user_input` or any equivalent between phases.
  Do not stop after Phase A and ask whether to proceed to Phase B.
  Do not stop after the audit reactivation passes and ask whether
  to ship the example benchmarks. The phases are an execution
  plan, not a checklist of separately-permissioned tasks. Run
  Phase 0 through Phase E end to end, push commits as they land,
  and open the PR after Phase 0 so reviewers can watch the
  remaining phases land.

  The only legitimate reasons to stop and surface a question are:
  (a) a hard blocker that cannot be resolved by re-reading this
  prompt, the linked ADRs, or the codebase (for example, a
  conflict between two of those sources where the conflict
  resolution is not obvious); (b) a `pytest` or `ruff` failure
  that persists after a reasonable repair attempt; or (c) the
  Stop conditions section reports green and the PR is ready for
  review. "Should I move on to the next phase" and "do you want
  me to also add X" are not legitimate reasons. Default behavior:
  pick the most defensible interpretation, document it inline,
  and continue.

  This applies to Claude Code, any other agent runner, and any
  human worker who reads this prompt and feels the urge to
  check in mid-flight. Do not check in mid-flight. Finish.

- Title the PR `M16.7: Time-varying transient contact voltage`.
- **No AI-assistant credits in any shipped artifact.** The
  `includeCoAuthoredBy: false` settings fix landed before M16.6
  and held; verify with
  `git log --format=%B dev/m16.7-transient-ac | grep -i
  'co-authored-by: claude'` before opening the PR. Empty output
  is the gate.

## Required reading (in order; ~50 minutes)

1. `PLAN.md` in full. Confirm `main` is at v0.22.0 (post-M16.6) and
   that no other PR is in flight on `dev/m16.7-*`. Read the M16.6
   entry in "Completed work log" carefully — the MMS Variant H
   experience (psi-only rate gate; phi blocks decoupled because
   Kane and Hurkx kernels have no carrier-density coupling on the
   continuity rows) is the relevant prior art if you find yourself
   tempted to add an MMS variant for M16.7. **Do not add one.**
   See "MMS scope" below.
2. `docs/IMPROVEMENT_GUIDE.md` § M16.7 (Why / Deliverable /
   Acceptance / Dependencies) and § 1 (Honest current state).
3. `docs/ROADMAP.md` § Capability matrix. M16.7 is the last row in
   the M16 umbrella; closing it lets the next reviewer compress
   the M16 rows into a single "shipped" entry and promote M17 or
   M19 as the new headline.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M16.7 touches
   the schema (Layer 2; new `contacts[].voltage_t` field), the
   transient runner (Layer 5; `semi/runners/transient.py`), and
   the audit suite (`tests/audit/test_06_transient_fft_vs_ac_sweep.py`).
   No new physics module.
5. `docs/PHYSICS.md` § 2.5 / § 3 (transient and AC sweep runners
   from M13.1 / M14). The cross-check that M16.7 enables is
   stated there as "audit case 06 is currently a `pytest.skip`."
6. **`tests/audit/test_06_transient_fft_vs_ac_sweep.py` end to
   end.** This is the file the worker has to make active. The
   placeholder records the AC reference value at V_DC = 0.4 V,
   F_HZ = 1.0e6, dV = 1.0e-3 V on the `rc_ac_sweep` benchmark
   geometry, then `pytest.skip`s with a tracking note: "comparing
   the FFT of I(t) under V(t) = V_DC + dV*sin(omega t) to
   ac_sweep Y(omega) requires either a `bc_voltage_callback`
   hook or a transient runner extension. Both are larger than
   the < 50 LOC bug-fix budget for this audit PR." M16.7 is
   that runner extension.
7. `tests/audit/_helpers.py` (`load_benchmark`, `require_dolfinx`,
   `write_csv`, `write_markdown`). The audit suite has its own
   recording conventions; the new active test must populate the
   same CSV / markdown artifacts as the skip-path does.
8. `semi/runners/transient.py` end to end. Two key locations:
   - L338-343 `_build_transient_bcs(voltages: dict)` builds the
     full DD Dirichlet stack from a per-contact voltage map.
     This already exists.
   - L422-458 the time loop body, where each step calls
     `_build_transient_bcs(static_voltages)` with a *fixed*
     dict. This is the line that has to become `voltages_at_t(
     t_next)` for time-varying contacts.
9. `schemas/input.v2.json` L482-517 (the existing
   `contacts[].voltage` and `contacts[].voltage_sweep` blocks).
   The new `voltage_t` field is a sibling.
10. `semi/runners/ac_sweep.py` `run_ac_sweep` and the `Y` admittance
    output. The audit test's reference value comes from this
    runner; do not modify it.
11. `benchmarks/rc_ac_sweep/rc_ac_sweep.json`. The geometry the
    audit test runs on (1D pn junction reverse-biased at -1.0 V,
    AC swept 1 Hz to 1 GHz). The audit test loads this benchmark
    and overrides the contact voltage to V_DC = 0.4 V.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** New `contacts[].voltage_t` field with
  two variants: `{type: "table", times: [...], values: [...]}` and
  `{type: "step", t0, v0, v1}`. Both must be expressible in
  `schemas/input.v2.json` (current strict v2.6.0 from M16.6),
  validated by `semi/schema.py`, and exercised by the audit case
  06 reactivation.
- **Schema versioning is binding.** Additive minor bump v2.6.0 to
  v2.7.0. v2.0.0 through v2.6.0 inputs must continue to validate
  and produce bit-identical results to v0.22.0 (no benchmark
  uses `voltage_t` today; every existing benchmark stays
  bit-identical by construction). Update `SCHEMA_SUPPORTED_MINOR`
  in `semi/schema.py`.
- **Mutual exclusion.** A contact cannot specify both `voltage_t`
  and `voltage_sweep`. The bias_sweep / equilibrium / ac_sweep
  runners ignore `voltage_t` (it is transient-only). The transient
  runner uses `voltage_t` if present; otherwise falls back to
  `voltage` (the existing fixed-bias behavior).
- **Five layers, enforced.** The `voltage_t` evaluation lives in
  `semi/runners/transient.py` (Layer 5). The schema block lives in
  `schemas/input.v2.json` (Layer 1) and the loader's conditional
  validation in `semi/schema.py` (Layer 2). No physics module is
  touched. No FEM-types leak.
- **MMS scope.** Do **not** add an MMS Variant I for M16.7. The
  per-physics-module MMS rule from ADR 0006 applies to *new
  physics kernels*; M16.7 ships no new kernels. Both the transient
  runner (M13.1) and the AC sweep runner (M14) already have MMS
  coverage. M16.7 is a *consistency check between two
  already-validated solvers*, not a new physics implementation.
  Audit case 06 is the V&V gate; that is what the IMPROVEMENT_GUIDE
  acceptance criteria specify. **If you find yourself reaching for
  an MMS variant, stop and re-read `docs/adr/0015-schottky-robin-bc.md`
  for the precedent on V&V departures.** A short note in
  `docs/PHYSICS.md` § Verification & Validation acknowledging
  M16.7's audit-only V&V is sufficient (Phase F step).
- **Per-runner threading.** Only the transient runner consumes
  `voltage_t`. The other six runners parse the field (so the
  schema validates) but raise a clear error if it is set on a
  contact for a non-transient solve. Match the M16.5 conditional-
  validation pattern: declare the schema field, do the
  conditional check in `semi/schema.py::validate`.
- **One milestone, one PR.** Do not bundle with M17 (heterojunctions),
  M19 (3D MOSFET), or M19.1 (MPI parallel). Do not retire the
  `mosfet_2d` `allow-failure: "true"` flag.
- **No em dashes in prose or code comments.**
- **No AI-assistant credits.** Settings layer is the durable fix.
- **Coverage gate is 95 in the gated `docker-fem-tests` job.**
  M16.5 needed a coverage-recovery follow-up commit; M16.6 needed
  *two* (one for `bias_sweep.py` schottky pre-solve, one for
  multi-region BBT/TAT/Schottky branches in `_mr` builders that
  no benchmark exercises). Plan unit tests **before** Phase C
  (the audit-test reactivation) and Phase D (the example
  benchmarks) so the new transient-runner branches are covered
  by the gated suite, not by the audit job or the benchmark
  matrix whose coverage is not merged. Specifically: every
  variant of `voltage_t` (`table` and `step`) gets at least
  one pytest in `tests/test_runners_transient.py` (or a new
  `tests/fem/test_transient_voltage_t.py`) that runs a tiny
  config and asserts the BC is set correctly per step. These
  tests are independent of audit case 06.
- **Bit-identity gate.** Every existing benchmark runs bit-
  identical to v0.22.0. Reuse all six prior anchors:
  `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2`,
  `diode_velsat_1d` 56.27 % @ 0.9 V, 0.19 % @ 0.3 V,
  `diode_auger_1d` >20 % divergence at 0.9 V,
  `diode_fermi_dirac_1d` 7.37 % FD-vs-Boltzmann V_bi,
  `schottky_1d` worst-case <10 % thermionic match,
  `zener_1d` worst-case <20 % Kane match.

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/`
and `pytest tests/`. Do not advance with red tests. **Do not
pause between phases for confirmation or status reports** (see
Branch and PR rules above for the full no-pause rule). Phase
gates are mechanical: green tests, push, advance to the next
phase. Run all six phases end to end.

---

### Phase 0: ship this starter prompt

Pure docs. Same shape as every M16.x Phase 0.

1. Save this entire file as `docs/M16_7_STARTER_PROMPT.md`.
   Strip nothing.
2. Append a one-line "Author M16.7 starter prompt" entry to
   `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under
   `[Unreleased]`.
3. Push the branch and open the PR. Verify the squash-body
   preview has no AI-credit trailer.

**Commit message:** `docs: ship M16.7 starter prompt (M16.7)`

**Acceptance:** the file exists at
`docs/M16_7_STARTER_PROMPT.md` on `dev/m16.7-transient-ac`; CI
green on docs-only diff.

---

### Phase A: schema surface for `voltage_t`

Pure-Python. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.6.0 to 2.7.0. Bump
   `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`. The `examples`
   array gains `"2.7.0"` ahead of existing entries.
2. Add `voltage_t` to the `contacts[]` block (sibling of
   `voltage` and `voltage_sweep`):
   - Type `object`, `additionalProperties: false`.
   - `type` field: enum `["table", "step"]`. Required.
   - For `type == "table"`:
     - `times`: array of numbers, monotonically increasing,
       length >= 2, units seconds.
     - `values`: array of numbers, same length as `times`,
       units volts.
     - The runner does linear interpolation for `t` between
       table points; clamps at the endpoints for `t < times[0]`
       or `t > times[-1]`.
   - For `type == "step"`:
     - `t0`: number, units seconds. The transition time.
     - `v0`: number, units volts. Voltage for `t < t0`.
     - `v1`: number, units volts. Voltage for `t >= t0`.
   - Description: "Time-varying applied voltage on this contact.
     Used only by the transient runner; ignored by bias_sweep,
     ac_sweep, equilibrium, mos_cv, mos_cap_ac, and resistor_3d.
     Mutually exclusive with `voltage_sweep`. The `table` variant
     interpolates linearly between (times[i], values[i]) pairs;
     the `step` variant returns `v0` for `t < t0` and `v1` for
     `t >= t0`. (M16.7)"
3. Conditional validation in `semi/schema.py`:
   - Reject `voltage_t` and `voltage_sweep` on the same contact.
   - For `type == "table"`: enforce `len(times) == len(values)
     >= 2` and strict monotonicity of `times`.
   - For non-transient solver types (`solver.type` is anything
     other than `"transient"`), reject configs where any contact
     has a `voltage_t` field set. Error message points the user
     at `voltage` or `voltage_sweep`.
4. Default-fill: a config with no `voltage_t` is bit-identical to
   v0.22.0 on every existing benchmark.
5. Pure-Python tests in a new `tests/test_voltage_t_schema.py`:
   - Every existing benchmark JSON validates unchanged against
     v2.7.0.
   - `voltage_t` `table` variant with valid times/values
     validates.
   - `voltage_t` `step` variant with t0/v0/v1 validates.
   - `times` non-monotonic is rejected.
   - `times` and `values` length mismatch is rejected.
   - `voltage_t` and `voltage_sweep` on the same contact is
     rejected.
   - `voltage_t` set on a `solver.type: "bias_sweep"` config is
     rejected.
   - Unknown `voltage_t.type` (`"sine"`) is rejected (note: the
     IMPROVEMENT_GUIDE only specifies `table` and `step`; the
     audit test uses the `table` variant with sampled-sinusoid
     values, so `sine` is not in scope here).
6. Update `docs/schema/reference.md` versioning table with the
   v2.7.0 row; the change-log column reads "additive: time-
   varying contact voltage `voltage_t` for the transient runner
   (M16.7)".

**Commit message:** `feat(schema): contacts[].voltage_t for transient runner (schema 2.7.0); existing solvers reject (M16.7)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change.

---

### Phase B: time-varying BC evaluation in the transient runner

Layer 5. The change is local to `semi/runners/transient.py`: a
`voltages_at_t(t)` callable replaces the fixed `static_voltages`
dict in the time-loop BC build.

1. Extend `semi/runners/transient.py`:
   - Add a private helper `_build_voltage_t_evaluator(cfg)` that
     parses the `voltage_t` block on each contact and returns a
     callable `voltages_at_t(t: float) -> dict[str, float]`.
     Contacts with no `voltage_t` return their fixed `voltage`
     value (the existing behavior).
   - For the `table` variant: cache the times and values arrays
     once at runner setup; the per-step lookup is
     `np.interp(t, times, values)`. Match endpoint-clamp behavior
     described in the schema docstring.
   - For the `step` variant: return `v0` if `t < t0` else `v1`.
   - At the time loop body (currently L457), replace
     `bcs = _build_transient_bcs(static_voltages)` with
     `bcs = _build_transient_bcs(voltages_at_t(t_next))`. Note
     that `t_next` is the *time the BC will hold over* (the
     step from `t_current` to `t_next`); using `t_next` matches
     the BDF-implicit convention since the residual is evaluated
     at the new time level.
   - For configs with no `voltage_t` on any contact, the
     evaluator returns the same dict at every t and the time-
     loop behavior is bit-identical to v0.22.0.
2. The BC ramp continuation (the pre-loop block where bias goes
   from 0 to its target) **does not** consume `voltage_t`. It
   ramps to the value of `voltage_t.values[0]` (table) or
   `voltage_t.v0` (step); after the ramp, the time loop takes
   over and `voltage_t` is in effect. Document this in the
   transient runner's module docstring.
3. Coverage tests in `tests/fem/test_transient_voltage_t.py`
   (new file). Eight tiny configs (1D pn diode, 4 cells, 5-10
   timesteps each) covering the full waveform / edge-case
   matrix:
   - `test_step_simple`: `voltage_t.type = "step"`, transition
     mid-loop, asserts the BC switches at exactly `t0`.
   - `test_step_at_t_zero`: `t0 = 0`, asserts the BC starts at
     `v1` (not `v0`).
   - `test_table_linear_ramp`: `times = [0, 1e-9, 2e-9, ...]`,
     `values` linear in time. Asserts mid-step interpolation is
     correct.
   - `test_table_sampled_sinusoid`: `values = V_DC + dV *
     sin(2*pi*f*times)`, the same waveform shape audit case 06
     uses. Asserts FFT bin alignment after run.
   - `test_table_nonuniform_times`: irregular `times` spacing.
     Asserts linear interpolation works across uneven samples.
   - `test_table_endpoint_clamp_low`: query at `t < times[0]`.
     Asserts the runner returns `values[0]` (the documented
     clamp behavior).
   - `test_table_endpoint_clamp_high`: query at `t > times[-1]`.
     Asserts the runner returns `values[-1]`.
   - `test_two_contacts_independent_voltage_t`: a config where
     two contacts each have their own `voltage_t` (one table,
     one step). Asserts the per-step BC dict has both values.
   Each test runs `run_transient` end-to-end (or instruments
   the runner via the existing `post_step_hook` to capture the
   per-step BC value) and asserts the per-step BC value matches
   the expected V(t) within 1e-12 relative.

   These are the gated-suite-coverage tests for the new
   branches. M16.5 and M16.6 both shipped with coverage gaps
   filled in follow-up commits because the new branches were
   only exercised by benchmarks (whose coverage isn't merged
   into the gate). **Do not repeat that pattern.** These
   tests run in the gated `docker-fem-tests` job.

**Commit message:** `feat(runners): transient.py voltage_t evaluator for time-varying contact voltage (M16.7)`

**Acceptance:** new tests pass; no FEM behavior change on
existing benchmarks (no benchmark uses `voltage_t` yet).

---

### Phase C: reactivate audit case 06

Replace the `pytest.skip` with the actual FFT-vs-AC comparison.
The reference (Y_ac) is already computed in the placeholder; you
add the transient-side FFT.

1. Edit `tests/audit/test_06_transient_fft_vs_ac_sweep.py`:
   - Keep the AC reference computation at the top.
   - Construct a transient config from the same benchmark with:
     - `solver.type: "transient"`.
     - `contacts[anode].voltage_t = {"type": "table", ...}`
       sampling `V(t) = V_DC + dV * sin(2*pi*F_HZ*t)` over
       N_periods >= 4 cycles at >= 32 samples per cycle (so the
       FFT bin at F_HZ has clean separation from neighbors and
       leakage is bounded).
     - `solver.transient.dt`: matches the sampling rate above.
     - `solver.transient.t_end`: matches the total sampled
       duration.
     - `solver.bc_ramp_steps`: ramp to `V_DC` (the table's
       starting value) before the time loop.
   - Run `run_transient` and capture I(t) at the same contact
     the AC sweep used for Y(omega).
   - Compute the FFT of I(t) at F_HZ. Convert to admittance:
     `Y_transient_fft = FFT(I)[bin_F] / FFT(V)[bin_F]`. Use a
     window function (Hann) to reduce leakage; document the
     window choice in the test docstring.
   - Assert `|Y_transient_fft - Y_ac| / |Y_ac| < 0.05` (5 %, per
     IMPROVEMENT_GUIDE § M16.7 acceptance test 1).
   - Replace the `pytest.skip` with the assertion.
   - Update the CSV / markdown writers to record both Y values
     and the relative error, with `status = "passed"` (or
     `"failed"` if the gate trips).
2. Tune the dt / t_end / N_periods knobs to make the gate clean
   on a typical CI runner. Suggested starting point:
   - F_HZ = 1.0e6 (per the existing placeholder).
   - dt = 1.0e-9 s (1000 samples per cycle).
   - N_periods = 8.
   - t_end = 8.0e-6 s (8000 timesteps).
   - Hann window applied to I(t) before FFT.
   This is heavy for a CI runner. If the audit job's wall-clock
   budget is tight, drop to N_periods = 4 and dt = 5e-9 s (200
   samples/cycle, 800 steps); document the choice in the test.
   The gate is 5 %; over-sampling beyond what the gate requires
   is wasted CI time.
3. Update `docs/PHYSICS_AUDIT.md` with the now-active status of
   case 06; remove the "skipped" note from any audit summary
   table.

**Commit message:** `test(audit): activate case 06 (transient FFT vs AC sweep) within 5 % (M16.7)`

**Acceptance test 1**: `pytest tests/audit/test_06_transient_fft_vs_ac_sweep.py`
exits 0 with the 5 % FFT-vs-AC match.

---

### Phase D: example benchmarks demonstrating `voltage_t`

The M16.x convention is "every milestone ships a benchmark."
M16.1 shipped `diode_velsat_1d`, M16.2 retuned `mosfet_2d`, M16.3
shipped `diode_auger_1d`, M16.4 shipped `diode_fermi_dirac_1d`,
M16.5 shipped `schottky_1d`, M16.6 shipped `zener_1d`. M16.7
breaks the pattern unless we ship at least one demonstration
benchmark per `voltage_t` variant.

The audit case is V&V; the benchmarks here are *demonstration*.
The verifiers can be lightweight (the audit case carries the
quantitative gate). Goal: a future user reading the repo learns
how to use `voltage_t` from a working config, not from the
audit-test code.

1. **`benchmarks/pn_1d_pulse/`** (the `step` variant
   demonstrator). Structure mirrors `benchmarks/pn_1d_turnon/`
   (which is the closest existing transient benchmark, but uses
   the legacy `bc_ramp_steps`-then-fixed-bias path):
   - `pn_pulse.json`: 1D pn junction (same geometry as
     `pn_1d_turnon`: 20 um, N_A = N_D = 1e17 cm^-3, tau = 1e-8
     s, 800 cells). Anode `voltage_t = {"type": "step", "t0":
     5e-9, "v0": 0.0, "v1": 0.6}`. `solver.transient.dt = 1e-10
     s, t_end = 5e-8 s` (250 timesteps; the step at t = 5 ns
     gives 50 pre-step and 200 post-step samples for clean
     turn-on dynamics on either side).
   - `README.md`: explain the difference from `pn_1d_turnon`
     (legacy fixed-bias path vs `voltage_t.step` path; the two
     should produce qualitatively the same I(t) for the post-
     step regime, and the worker can confirm this in the README).
   - Verifier in `scripts/run_benchmark.py`: `verify_pn_1d_pulse`.
     Lightweight: assert that I(t) is finite and non-NaN
     everywhere, that the post-step current at t = 50 ns
     matches the `pn_1d_turnon` final-time current within
     20 % (the two paths use slightly different solver
     warmup; tighter than 20 % over-claims the equivalence),
     and that the pre-step current is order-of-magnitude
     smaller than the post-step current. Mirror the
     `verify_pn_1d_turnon` structure.
   - Plotter `plot_pn_1d_pulse`: I(t) on linear axis, V(t)
     on the same plot (twin y-axis is fine).
   - CI matrix entry near `pn_1d_turnon`. Do **not** mark
     `allow-failure: true`.

2. **`benchmarks/diode_sine_1d/`** (the `table` variant
   demonstrator). Sinusoidal large-signal bias around a forward-
   bias DC operating point. Demonstrates the full audit-case-06
   shape but as a self-contained benchmark with its own README
   and verifier:
   - `diode_sine.json`: same 1D pn geometry. Anode `voltage_t
     = {"type": "table", "times": [...], "values": [...]}`
     where `values = 0.4 + 0.05 * sin(2 pi * 1e6 * times)`
     (V_DC = 0.4 V, dV = 50 mV, F = 1 MHz; large-signal, not
     small-signal, so the response has visible harmonic
     content). 4 periods, 200 samples per period, 800
     timesteps total.
   - `README.md`: device description; note that the small-
     signal limit of this benchmark recovers the AC sweep
     admittance Y(omega), and that audit case 06 carries the
     formal gate; this benchmark is illustrative, not a V&V
     gate.
   - Verifier `verify_diode_sine_1d`: lightweight. Assert I(t)
     finite and non-NaN; assert the FFT of I(t) has its
     dominant peak at F = 1 MHz (the input fundamental); assert
     the second-harmonic amplitude is at most 50 % of the
     fundamental (large-signal nonlinearity; if it dominates,
     the timestep is too coarse or the table is mis-built).
     Do **not** assert the fundamental amplitude against the
     small-signal Y; that is audit case 06's job and a tighter
     gate.
   - Plotter `plot_diode_sine_1d`: I(t) and V(t) on twin axes
     in the time domain; |I(omega)| on a log-y plot in the
     frequency domain showing the fundamental and the first
     few harmonics.
   - CI matrix entry near `pn_1d_turnon`. Do **not** mark
     `allow-failure: true`.

3. **`notebooks/10_voltage_t_variants.ipynb`** (optional; the
   pattern matches the optional notebooks shipped with M16.3
   `diode_auger_1d` and M16.4 `diode_fermi_dirac_1d`):
   - Cell 1: load `pn_1d_pulse`, run, plot I(t) and V(t).
     Annotate the step transition.
   - Cell 2: load `diode_sine_1d`, run, plot I(t) and the FFT.
     Show the fundamental and the first three harmonics.
   - Cell 3: explain the relationship to the M14 AC sweep and
     point at audit case 06 for the formal cross-check.
   - Number `10_*` (M16.6 shipped `09_zener_1d` per its
     starter prompt; if M16.6 shipped a different number,
     bump accordingly).

**Commit message:** `feat(benchmark): pn_1d_pulse and diode_sine_1d demonstrating voltage_t (M16.7)`

**Acceptance:** both new benchmark verifiers pass; CI matrix
entries green; existing benchmarks bit-identical to v0.22.0.

---

### Phase E: closeout

This is the final M16 closeout. Phase F in this PR also closes
the umbrella: the gap list in PLAN, IMPROVEMENT_GUIDE, and
PHYSICS_INTRO compress to "no remaining M16 gaps."

1. `PLAN.md`:
   - Move M16.7 from "Next task" to "Completed work log" with PR
     number, deliverables, schema minor bump (2.6.0 to 2.7.0),
     and acceptance-test results (cite observed worst-case
     `|Y_transient_fft - Y_ac| / |Y_ac|` from Phase C).
   - **Update the gap list at L276** from:
     ```
     **Physics gaps:** no transient FFT vs AC sweep validation.
     M16.7. (Field-dependent mobility shipped in M16.1; ...)
     ```
     to:
     ```
     **Physics gaps:** none in the M16 umbrella (all of M16.1
     through M16.7 shipped). Next-tier gaps in M17
     (heterojunctions, position-dependent chi and Eg) and M19
     (3D MOSFET capstone).
     ```
   - Set "Next task" to the maintainer's choice between
     **M17 (heterojunctions)** and **M19 (3D MOSFET capstone)**.
     Both are documented in IMPROVEMENT_GUIDE § 4 / § 5 with
     their own acceptance criteria. Do not author a starter
     prompt for either; that is the next reviewer's deliverable.
     Note that M19 depends on M16.1 (mobility) and is therefore
     unblocked; M17 depends on M16.4 (FD; needed at heterojunction
     barriers) and is also unblocked.
   - Refresh "Current state" with the new package version (bump
     0.22.0 to 0.23.0).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M16.7 Done in § 4 with a one-line summary and a
     CHANGELOG anchor.
   - **Update the gap line at L81-82** from:
     ```
     - **No transient FFT vs AC sweep validation.** M16.7.
       (M16.1 Caughey-Thomas mobility, ...)
     ```
     to:
     ```
     - **M16 umbrella complete.** All seven physics-completeness
       slices (M16.1 Caughey-Thomas mobility, M16.2 Lombardi
       surface mobility, M16.3 Auger, M16.4 Fermi-Dirac, M16.5
       Schottky contacts, M16.6 BBT and TAT tunneling, M16.7
       time-varying transient contact voltage with FFT-vs-AC
       audit) shipped; see § 4 for per-milestone Done entries.
       Next-tier gaps in M17 (heterojunctions) and M19 (3D
       MOSFET capstone).
     ```
   - Append a § 9 changelog entry under `[0.23.0]`. Move the
     existing `[Unreleased]` block into `### Released`.
3. `docs/PHYSICS_INTRO.md`:
   - § 6: append a transient-with-V(t) bullet:
     ```
     - **Time-varying transient.** Sinusoidal, table-interpolated,
       or step contact voltage in `run_transient` (M16.7,
       audit case 06). The audit suite cross-checks an FFT of
       I(t) against the M14 small-signal AC sweep at the same
       operating point.
     ```
     Insert after the BBT/TAT bullet from M16.6.
   - § 7: drop the "no transient FFT vs AC" bullet entirely. The
     remaining § 7 bullets are impact ionization, heterojunctions,
     and 3D MOSFET. Update the closing paragraph from "M16.6
     and M16.7 close these" (or whatever it now says) to "M17
     and M19 are the next-tier milestones; M16 is complete."
4. `docs/ROADMAP.md`:
   - Update the M16.7 row in the capability matrix from Planned
     to shipped.
   - Compress the seven M16 rows (M16.1 through M16.7) into a
     single shipped row labeled "M16: Physics completeness
     (Caughey-Thomas, Lombardi, Auger, FD, Schottky, BBT/TAT,
     transient FFT-vs-AC) — all shipped" if the table format
     supports it; otherwise leave them as seven shipped rows.
     The maintainer can adjust at review.
   - Promote either M17 or M19 (whichever is next-task) to the
     top of the Planned column.
5. `CHANGELOG.md`:
   - New `[0.23.0]` entry with the M16.7 line items.
   - Add a one-line "M16 umbrella complete" note at the top of
     the `[0.23.0]` entry.
   - Update the schema-banner comment at the top to mention
     v2.7.0.
6. `pyproject.toml` and `semi/__init__.py`: bump 0.22.0 to
   0.23.0.
7. Push every commit to `origin/dev/m16.7-transient-ac`. Verify
   no AI-credit trailers anywhere
   (`git log --format=%B dev/m16.7-transient-ac | grep -i
   'co-authored-by: claude'` returns empty).

**Commit message:** `docs: close out M16.7 and the M16 umbrella (PLAN, IMPROVEMENT_GUIDE, PHYSICS_INTRO, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by
      this PR.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact.
- [ ] `git log --format=%B dev/m16.7-transient-ac | grep -i
      'co-authored-by: claude'` returns empty before opening
      and merging the PR.
- [ ] Pure-Python core remains dolfinx-free
      (`tests/test_lazy_imports.py` clean).
- [ ] Every existing benchmark is bit-identical to v0.22.0
      (`voltage_t` is unset by default; the time-loop BC build
      is a no-op functional change in that path).
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004).
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `semisim_server` public
      API.
- [ ] Schema bumped per minor (2.6.0 to 2.7.0); v2.0.0 through
      v2.6.0 inputs still validate.
- [ ] MMS rate gate L2 >= 1.99 / H1 >= 0.99 active for Variants
      A through H (no new variant; M16.7 ships no new physics
      module).
- [ ] Coverage gate holds at 95 in the gated `docker-fem-tests`
      job. The `voltage_t` evaluator branches are covered by
      `tests/fem/test_transient_voltage_t.py`, not just by
      audit case 06.
- [ ] Audit case 06 is active and passing (no `pytest.skip`).

## Anti-goals

- Do not add an MMS Variant I. M16.7 ships no new physics
  kernel; the per-physics-module MMS rule (ADR 0006) does not
  apply. Audit case 06 is the V&V gate; it is what
  IMPROVEMENT_GUIDE § M16.7 acceptance criteria specify.
- Do not add a `sine` variant of `voltage_t`. The IMPROVEMENT_GUIDE
  specifies `table` and `step` only; sinusoidal V(t) is the
  audit test's input but it is materialized via the `table`
  variant with sampled values. Adding `sine` widens the schema
  beyond the milestone's scope.
- Do not generalize `voltage_t` to non-transient runners. The
  schema rejection on `solver.type != "transient"` is
  intentional; `voltage_sweep` is the right field for bias
  sweeps and `voltage` is right for AC. Future cross-cutting
  voltage models (e.g. a piecewise-linear ramp shared between
  bias_sweep and transient) are out of scope.
- Do not change anything in `semi/physics/` or `semi/bcs.py`.
  M16.7 is a runner extension plus a schema field; the BC
  builder API is unchanged because the only thing that varies
  is the per-step voltage *value* the existing builder consumes.
- Do not retire the `mosfet_2d` `allow-failure: "true"` flag.
- Do not bundle with M17, M19, M19.1, M20, or any M14.2.x
  backlog item.
- Do not add `Co-Authored-By` trailers, generated-by footers,
  or any other AI-credit marker.

## Stop conditions

You are done when:

1. The M16.7 PR is opened on branch `dev/m16.7-transient-ac`.
2. Both acceptance tests in `docs/IMPROVEMENT_GUIDE.md` § M16.7
   pass in CI:
   - A1: audit case 06 (`tests/audit/test_06_transient_fft_vs_ac_sweep.py`)
     active and passing within 5 %.
   - A2: every existing transient benchmark is bit-identical to
     v0.22.0 (`pn_1d_turnon` is the regression anchor).
3. Coverage tests `tests/fem/test_transient_voltage_t.py` (eight
   waveform / edge-case variants per Phase B step 3) pass in the
   gated `docker-fem-tests` job; coverage gate holds at 95
   without a follow-up commit.
4. Both example benchmarks pass: `pn_1d_pulse` (step variant,
   demonstrating switching transient) and `diode_sine_1d`
   (table variant, demonstrating large-signal sinusoidal drive)
   exit 0 with their lightweight verifiers passing.
5. PLAN.md, IMPROVEMENT_GUIDE.md, PHYSICS_INTRO.md, ROADMAP.md,
   CHANGELOG.md reflect the closeout *and* the M16 umbrella
   compression. The gap-list rewrites at PLAN.md L276 and
   IMPROVEMENT_GUIDE.md L81-82 land verbatim per Phase E.
6. Package version bumped 0.22.0 to 0.23.0 in `pyproject.toml`
   and `semi/__init__.py`.
7. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI
   assistant. Squash-merge body is also clean.
8. PR reviewed, CI green (modulo the documented `allow-failure`
   on `mosfet_2d`, unchanged in scope from v0.22.0), merged.

## PR description template

```
## Summary

M16.7: Time-varying transient contact voltage, the seventh and
final physics-completeness slice of the M16 umbrella. Ships no
new physics kernel; extends `semi/runners/transient.py` to
accept a `contacts[].voltage_t` block (variants `table` and
`step`) and reactivates audit case 06 (transient FFT of I(t)
under sinusoidal V(t) compared to the M14 small-signal AC sweep
admittance Y(omega)) within 5 %.

Two demonstration benchmarks ship the new feature: `pn_1d_pulse`
exercises the `step` variant with a switching transient at
t = 5 ns; `diode_sine_1d` exercises the `table` variant with a
1 MHz, 50 mV sinusoidal large-signal drive around V_DC = 0.4 V.
Both have lightweight verifiers (the formal V&V gate is audit
case 06).

Schema additive minor bump v2.6.0 to v2.7.0 (new
contacts[].voltage_t object). v2.0.0 through v2.6.0 inputs
continue to validate; no existing benchmark uses voltage_t, so
every prior benchmark is bit-identical to v0.22.0.

After this PR the M16 umbrella is complete. The gap list in
PLAN.md, IMPROVEMENT_GUIDE.md, and PHYSICS_INTRO.md compress to
"no remaining M16 gaps"; next-tier work is M17 (heterojunctions)
or M19 (3D MOSFET capstone).

## Acceptance tests

(Both numbered in docs/IMPROVEMENT_GUIDE.md § M16.7.)

- [ ] A1 (FFT match): tests/audit/test_06_transient_fft_vs_ac_sweep.py
      transient FFT vs ac_sweep Y(omega) at V_DC = 0.4 V,
      F = 1 MHz, dV = 1 mV agree within 5 %
      (observed: <fill in>)
- [ ] A2 (byte-identity): pn_1d_turnon bit-identical to
      v0.22.0; every other existing benchmark bit-identical
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
      diode_auger_1d >20 % divergence at 0.9 V;
      diode_fermi_dirac_1d 7.37 % FD-vs-Boltzmann V_bi at
      N_D=1e20 cm^-3; schottky_1d worst-case <10 %;
      zener_1d worst-case <20 % Kane match)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95 (in the gated
      docker-fem-tests job)
- [ ] pytest tests/audit/ (audit case 06 active)
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark pn_1d_pulse
      (step variant demonstrator)
- [ ] docker compose run --rm benchmark diode_sine_1d
      (table variant demonstrator)
- [ ] docker compose run --rm benchmark pn_1d_turnon
      (no voltage_t, bit-identical)
- [ ] All other benchmarks bit-identical to v0.22.0 via the
      existing CI matrix.

## Notes

- M16.7 ships no MMS variant; M16 umbrella tops out at Variant H
  (M16.6 BBT/TAT). The transient runner (M13.1) and AC sweep
  runner (M14) are the components being cross-checked, both
  already with MMS coverage.
- The mosfet_2d CI matrix entry remains `allow-failure: "true"`,
  unchanged in scope from M16.6.
- M16 umbrella complete after this PR. Per-milestone done
  entries in IMPROVEMENT_GUIDE § 4.
- Coverage gate held at 95 in the gated docker-fem-tests job
  via tests/fem/test_transient_voltage_t.py, avoiding the
  follow-up-commit pattern from M16.5 / M16.6.
```

## Hand-off

When M16.7 lands and is reviewed, the M16 umbrella closes. The
two next-tier candidates are:

- **M17: heterojunctions.** Position-dependent electron affinity
  chi(x) and band gap E_g(x) on multi-region meshes. Touches
  `semi/physics/poisson.py` and `semi/physics/slotboom.py`
  (the band-edge offsets are no longer constants). Depends on
  M16.4 (FD; the non-degenerate approximation breaks at
  heterojunction barriers).

- **M19: 3D MOSFET capstone.** A real 3D MOSFET on a gmsh-
  sourced unstructured mesh, exercising the M15 GPU linear-
  solver path and the M16.1 Caughey-Thomas mobility under
  non-trivial geometry. Largest single benchmark in the repo
  by problem size; expected to need MPI parallel orchestration
  (M19.1) for runtime to be tolerable.

Either is unblocked. The maintainer chooses; the next starter
prompt is authored against the chosen milestone in the same
shape as this prompt.
