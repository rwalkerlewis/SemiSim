## Context

You are working in `kronos-semi` at `v0.23.0` (post-M16.7, package
version `0.23.0`, schema `2.7.0`). The repo is at
`https://github.com/rwalkerlewis/kronos-semi`; `main` is at
`1a372b5` (`M16.7: Time-varying transient contact voltage (#84)`).

Your assignment is **the examples catalogue** — a new
`examples/` top-level directory containing three self-contained
practical device configs that exercise the M16 physics catalogue
in real engineering contexts. This is **not** an M-numbered
milestone. It is documentation / demonstration work that runs
parallel to the M17 / M19 next-task path.

The repo today has two layers of working configs:

- `benchmarks/` (13 entries) — each is a V&V gate against an
  analytical or textbook reference. Tight tolerances, narrow
  device parameters, intentionally minimal geometry. Not
  designed as starting points for user work.
- `notebooks/` — Colab walkthroughs of specific physics topics.
  Conceptual, not config-first.

Missing: the practical middle layer. A user who wants to simulate
"an NMOS like the ones in my fab" or "a Schottky diode at three
temperatures" has nowhere to clone-and-modify from. The
examples catalogue fills that gap with three self-contained
device configs, each exercising a different runner and a
different combination of M16 physics features.

Each example is illustrative, not a V&V gate. CI gates them on
"the run completes, the output JSON is well-formed, and the
recorded I-V values are finite and non-NaN." Nothing tighter.
The benchmarks directory remains the V&V layer; users who want
quantitative correctness gates look there. Examples are for
"how do I set this up in practice."

This prompt **is** the starter. Save it as
`docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md` as the first commit
of the PR (Phase 0 below). Note the file name is descriptive
rather than M-numbered; the convention is "M-numbered prompts
go in `docs/M*_STARTER_PROMPT.md`; non-milestone prompts get
descriptive names." Keeping these distinct preserves the
milestone-discipline of the M-numbered prompts.

## Branch and PR rules

- Work on a fresh branch off `main`:
  `git checkout -b dev/examples-catalogue`. Not an M-numbered
  branch; the parallel path is intentional.
- One scope, one PR. Do not bundle with M17 / M19 work or any
  M16.x maintenance.
- Push every phase commit immediately. Open the PR after
  Phase 0 so reviewers can watch examples land.
- **Run all five phases consecutively without pausing for
  confirmation, status reports, or "should I continue"
  prompts.** Do not call `ask_user_input` or any equivalent
  between phases. Run Phase 0 through Phase D end to end.
  The only legitimate reasons to stop are: (a) a hard blocker
  that re-reading this prompt and the codebase cannot
  resolve; (b) a `pytest` or `ruff` failure that persists
  after a reasonable repair attempt; (c) the Stop conditions
  section reports green and the PR is ready for review.
  "Should I move on to the next example" is not a legitimate
  reason. Default behavior: pick the most defensible
  interpretation, document it in a commit body, continue.
- Title the PR `examples: catalogue of practical device configs`.
  Use the PR description template at the bottom of this prompt.
- **No AI-assistant credits in any shipped artifact.** The
  `includeCoAuthoredBy: false` settings fix is presumed in
  effect (it held through M16.6 and M16.7); verify with
  `git log --format=%B dev/examples-catalogue | grep -i
  'co-authored-by: claude'` before opening the PR.

## Required reading (in order; ~25 minutes)

1. `PLAN.md` "Current state" and "Next task". Confirm `main` is
   at v0.23.0 (post-M16.7) and that no other PR is in flight on
   `dev/examples-*`. The "Next task" section names M17 or M19
   as the maintainer's choice; the examples PR is parallel and
   does not disrupt that.
2. `docs/IMPROVEMENT_GUIDE.md` § 4 / § 5. The M16 catalogue is
   the reference for which physics each example will exercise.
3. `benchmarks/` directory tour. Look at `pn_1d_bias`,
   `mosfet_2d`, `schottky_1d`, `pn_1d_turnon`, and the M16.7
   demonstrators (`pn_1d_pulse`, `diode_sine_1d`). The example
   configs share JSON structure but use different parameter
   regimes and different verifier shapes.
4. `schemas/input.v2.json` (current strict v2.7.0). The
   examples use only schema features that have already shipped;
   no schema bump.
5. `scripts/run_benchmark.py` end to end. The examples register
   with the same `@register("name")` mechanism the benchmarks
   use. Plotters are encouraged but optional; verifiers are
   minimal (smoke-test only).

## Conventions (project rules, not suggestions)

- **No schema changes.** v2.7.0 is the current strict version.
  The examples use only schema features shipped in M16.1-M16.7.
- **No new physics modules.** The examples consume the existing
  physics catalogue; they do not add to it.
- **Five layers, enforced.** The examples directory is
  pure-config (Layer 1, JSON) plus a small bit of Layer 5
  (verifier registration in `run_benchmark.py`). No changes to
  Layers 2-4.
- **Each example is self-contained.** `examples/foo/` contains
  `foo.json` plus `README.md`. The README is the user-facing
  documentation; it explains the device, the physics features
  exercised, the expected results, and the "what to change to
  adapt this to your problem" guide.
- **README structure (load-bearing).** Every example README has
  these sections, in this order:
  1. **What this is.** One paragraph: device type, geometry,
     intended use case.
  2. **Physics features exercised.** Bullet list pointing at
     the M-numbers and the relevant
     `physics.{mobility,recombination,statistics,...}` keys.
  3. **How to run.** Exact `docker compose run` command, then
     where to find the output.
  4. **Expected output.** Numerical values (not formal
     tolerances) the user should see; one or two
     "if you see something very different, here's what to
     check" lines.
  5. **How to adapt.** A "to model your own device" section
     that points at which parameters to change first.
- **Smoke-test verifiers, not V&V gates.** Each example
  registers a `verify_<name>` in `run_benchmark.py` that
  asserts: the run completed, the recorded I-V (or C-V, or
  transient I(t)) values are finite and non-NaN, and the
  result schema is well-formed. **Do not** assert numerical
  correctness against an analytical reference — that is what
  `benchmarks/` is for. Worth restating because the M16.x
  habit is to assert tight tolerances; the examples are
  illustrative and tight tolerances would create churn when
  parameters are tuned.
- **CI gating.** Each example gets a CI matrix entry in
  `.github/workflows/ci.yml`, near the existing benchmark
  entries. **Do not** mark `allow-failure: true`. The smoke
  test is the gate; if a run fails, the example is broken and
  the PR should not merge.
- **Naming.** `examples/<device_class>_<test_type>/` (snake
  case). Examples in this PR:
  - `examples/nmos_idvgs/`
  - `examples/schottky_iv_temperature/`
  - `examples/power_diode_reverse_recovery/`
- **One example, one phase, one commit.** Do not bundle.
- **No em dashes in prose or code comments.**
- **No AI-assistant credits.** Reiterated.
- **Coverage gate is 95 in `docker-fem-tests`.** Examples
  themselves do not touch `semi/`, so the gate is held by the
  existing tests. The verifier registrations in
  `run_benchmark.py` need a small `tests/test_examples_register.py`
  asserting each example name is registered (one-line
  pytest per example). This is the only addition to the
  gated suite.
- **Existing-benchmark byte-identity.** Every existing
  benchmark is bit-identical to v0.23.0 (the examples PR
  touches no existing config). The four prior anchors hold
  unchanged.

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check` (no-op for
JSON-only changes; still run it on any Python touched) and
`pytest tests/`. Do not advance with red tests. **Do not pause
between phases for confirmation.**

---

### Phase 0: ship this starter prompt and the examples scaffold

1. Save this entire file as
   `docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md`. Strip nothing.
2. Create `examples/README.md` (top-level entry point):
   - One paragraph explaining what `examples/` is for: "Each
     subdirectory is a self-contained practical device
     configuration. Examples are illustrative, not V&V gates;
     for analytical correctness gates, see `benchmarks/`."
   - A table of the three examples shipped in this PR, with
     columns: name, device class, runner, M16 features.
   - A "how to add an example" section pointing at this
     starter prompt for the conventions.
3. Append a one-line "Author examples catalogue starter prompt"
   entry to `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under
   `[Unreleased]`.
4. Push the branch and open the PR.

**Commit message:** `docs: ship examples catalogue starter prompt and scaffold`

**Acceptance:** `docs/EXAMPLES_CATALOGUE_STARTER_PROMPT.md` and
`examples/README.md` exist on `dev/examples-catalogue`; CI green
on docs-only diff; no AI-credit lines.

---

### Phase A: `examples/nmos_idvgs/` — NMOS Id-Vgs characterization

A practical NMOS transistor with realistic short-channel-ish
parameters and a four-curve Id-Vgs sweep at V_DS = 0.05, 0.5,
1.0, and 1.8 V. Exercises Caughey-Thomas mobility (M16.1) for
the channel and Lombardi surface mobility (M16.2) for the
inversion regime.

1. `examples/nmos_idvgs/nmos.json`. Geometry: 2D MOSFET, channel
   length 100 nm (longer than `mosfet_2d`'s toy geometry but
   short enough to show velocity saturation effects), gate oxide
   2 nm SiO2, p-type body N_A = 5e17 cm^-3, n+ source/drain
   N_D = 1e20 cm^-3 (degenerate; FD applies; set
   `physics.statistics: "fermi_dirac"`). Mesh: 200 x 80 cells
   (40 in oxide, 40 in silicon depth-wise; 200 in channel
   direction). `physics.mobility = {model: "lombardi",
   bulk_model: "caughey_thomas", interface_facet_tag: <oxide
   silicon facet tag>}`. Voltage sweep: parametric over
   V_DS in [0.05, 0.5, 1.0, 1.8] V with V_GS sweep [0, 1.8] V
   step 0.05 V at each V_DS.
2. `examples/nmos_idvgs/README.md` per the README structure
   convention. The "Expected output" section lists I_D at
   V_GS = 1.8 V, V_DS = 1.8 V (order-of-magnitude estimate based
   on textbook square-law I_D = (mu_eff C_ox W / 2 L) (V_GS -
   V_T)^2; the user is told to expect a value within a factor
   of 3 of that estimate, since real Lombardi mobility deviates
   from constant). Threshold voltage extraction guide using the
   linear-region tangent method; pointer at `mos_2d` for the
   companion C-V analysis.
3. Verifier `verify_nmos_idvgs` in `scripts/run_benchmark.py`:
   - Run completed (no exception, no NaN in output).
   - For every V_DS, the I-V curve is monotonically non-decreasing
     in V_GS (transistor is on the right side of subthreshold).
   - I_D at V_GS = 1.8 V, V_DS = 1.8 V is positive and finite
     (smoke check; no analytical bound).
4. Plotter `plot_nmos_idvgs`: four Id-Vgs curves overlaid on a
   semilog-y plot (subthreshold visible) plus the same curves
   on a linear-y plot (saturation visible). Two-panel layout.
5. CI matrix entry. Estimate the runtime: 200 x 80 mesh times
   four V_DS values times ~36 V_GS points per sweep is
   non-trivial; if the runtime exceeds 5 minutes on the CI
   runner, drop the V_GS step to 0.1 V or the V_DS sample
   count to two. Document the chosen knobs in the README.

**Commit message:** `examples: nmos_idvgs (Caughey-Thomas + Lombardi + FD)`

**Acceptance:** smoke verifier passes; CI green.

---

### Phase B: `examples/schottky_iv_temperature/` — Pt-on-n-Si Schottky at three temperatures

A practical Schottky diode characterization across temperature
(250 K, 300 K, 350 K). Exercises Schottky thermionic emission
(M16.5) under temperature variation. Runs three solves and
overlays the I-V curves.

1. `examples/schottky_iv_temperature/schottky_T.json`. Geometry
   matches `benchmarks/schottky_1d` structurally (1D Pt-on-n-Si,
   N_D = 1e16 cm^-3 uniform, 5 um total) but with three sibling
   configs in one directory or one parametric config that the
   verifier loops over temperature. The `physics.temperature`
   field already exists in the schema; the JSON uses the
   parametric mechanism (one config, the verifier overrides
   `physics.temperature` per run). If the parametric mechanism
   is not clean, ship three separate JSONs (`schottky_250K.json`,
   `schottky_300K.json`, `schottky_350K.json`) and the verifier
   loads all three. The README explains the choice.
2. Voltage sweep V_F in [0, 0.5] V step 0.025 V (matching the
   benchmark for parameter consistency).
3. `examples/schottky_iv_temperature/README.md` per the README
   structure convention. The "Expected output" section explains
   the temperature dependence: `J_sat = A* T^2 exp(-q phi_B /
   kT)`, so doubling T from 175 K to 350 K should *change* the
   saturation current by roughly 10x at this barrier height.
   The user can use this to verify their installation reproduces
   physically-reasonable temperature-dependent thermionic
   emission. Threshold-voltage tweak guide ("if you want to
   model a different barrier height, change `barrier_height_eV`
   in each config; barrier height is the dominant temperature-
   coefficient knob").
4. Verifier `verify_schottky_iv_temperature` in
   `scripts/run_benchmark.py`:
   - All three runs completed.
   - For each T, I-V curve is monotonically non-decreasing in
     V_F.
   - At V_F = 0.3 V (mid-sweep), I(350 K) > I(300 K) > I(250 K)
     monotonically (the temperature-dependent thermionic
     emission ordering must be physically correct; this is a
     smoke check, not a tight numerical gate).
5. Plotter `plot_schottky_iv_temperature`: three semilog-y
   I-V curves overlaid, color-coded by temperature.
6. CI matrix entry. Three small 1D runs are cheap; no runtime
   tuning needed.

**Commit message:** `examples: schottky_iv_temperature (Schottky + temperature dependence)`

**Acceptance:** smoke verifier passes; CI green.

---

### Phase C: `examples/power_diode_reverse_recovery/` — rectifier turn-off transient

A practical power-rectifier diode exercising the M16.7
voltage_t time-varying contact and M16.3 Auger recombination
at high injection. Demonstrates the standard reverse-recovery
waveform: forward conduction, abrupt voltage reversal, the
characteristic reverse-current spike during stored-charge
clearance, then steady reverse blocking.

1. `examples/power_diode_reverse_recovery/diode_recovery.json`.
   Geometry: 1D pn diode, N_A = 1e15 cm^-3, N_D = 5e15 cm^-3
   (light doping; allows high injection at modest forward
   bias and a meaningful stored charge), 100 um total length
   (long-base diode geometry; the diffusion length sets the
   recovery time scale). 800 cells. tau = 1e-7 s (carrier
   lifetime; controls Auger-vs-SRH balance and reverse-recovery
   shape). `physics.recombination = {srh: true, auger: true}`
   (Auger contributes meaningfully under forward conduction at
   high injection in a light-doped diode). `voltage_t` on the
   anode: piecewise-linear table from V = +0.7 V at t = 0 to
   V = -2.0 V at t = 100 ns, with a 10 ns transition at
   t = 50-60 ns. The materialized table has ~50 sample points
   over the full trajectory. Time loop: dt = 1e-9 s,
   t_end = 200 ns (200 timesteps).
2. `examples/power_diode_reverse_recovery/README.md` per the
   README structure convention. The "Expected output" section:
   peak forward current density during conduction
   (~1e6 A/m^2 order-of-magnitude); peak reverse current spike
   magnitude (~1e5 A/m^2 order-of-magnitude, smaller than
   forward by ~10x; this is the stored-charge clearance);
   recovery time scale (~50-100 ns; this is roughly tau = 100
   ns scaled by the depletion-region transit). The user is
   told to expect a "current dip and recovery to near-zero
   reverse current" qualitative shape; tightening the gate is
   benchmark-territory work.
3. Verifier `verify_power_diode_reverse_recovery` in
   `scripts/run_benchmark.py`:
   - Run completed.
   - I(t) is finite and non-NaN at every timestep.
   - max(I(t) for t < 50 ns) > 0 (forward conduction is
     positive).
   - min(I(t) for t > 60 ns) < 0 (reverse recovery dips
     negative).
   - I(t = 200 ns) > min(I(t)) (the diode is settling toward
     reverse-blocking; the final-time current is closer to
     zero than the recovery minimum).
4. Plotter `plot_power_diode_reverse_recovery`: I(t) and V(t)
   on twin y-axes, time on x-axis. Annotate the conduction,
   transition, recovery, and settled regions on the plot.
5. CI matrix entry. 1D 800-cell, 200-timestep transient is
   moderately expensive; estimate ~1-2 minutes on a typical
   runner. Acceptable.

**Commit message:** `examples: power_diode_reverse_recovery (voltage_t + Auger transient)`

**Acceptance:** smoke verifier passes; CI green.

---

### Phase D: closeout

1. `examples/README.md`: confirm the table-of-three-examples
   matches what shipped (Phases A-C). Add a "future examples"
   note pointing at the topics not yet covered (NPN BJT
   gummel plot once 3-terminal infrastructure is more mature;
   NMOS C-V at multiple body biases; tunnel diode forward I-V
   once forward-bias BBT is validated; SiC Schottky once SiC
   material parameters are added).
2. `tests/test_examples_register.py` (new file): three tiny
   pytests asserting each example's verifier name is
   registered in the `_VERIFIERS` dict in `run_benchmark.py`.
   This is the gated-suite-coverage gate for the new
   verifier registrations.
3. `PLAN.md` "Backlog": add a one-line entry under Backlog
   noting that the examples catalogue shipped (no milestone
   number; this is documentation work). The "Next task"
   field is unchanged (M17 or M19, maintainer's choice).
4. `docs/IMPROVEMENT_GUIDE.md` § 9 changelog: append an
   "examples catalogue: nmos_idvgs, schottky_iv_temperature,
   power_diode_reverse_recovery" entry under `[Unreleased]`.
   Move the existing `[Unreleased]` block (the Phase 0 entry
   from this PR) into `### Released` as part of the closeout.
5. `CHANGELOG.md`: this PR does not bump the package version
   (no API or schema change). Add an entry under `[Unreleased]`
   in CHANGELOG noting the examples directory ship; the next
   *physics* PR (M17 or M19) will roll this into its own
   version bump.
6. `docs/PHYSICS_INTRO.md` § 6: append a one-bullet entry
   pointing at the examples directory:
   ```
   - **Practical examples.** See `examples/` for self-contained
     device configs (NMOS Id-Vgs, Schottky temperature
     dependence, power diode reverse recovery) demonstrating
     M16 features in real engineering contexts. Distinct from
     `benchmarks/` (which gates V&V).
   ```
   Insert at the end of § 6.

**Commit message:** `docs: close out examples catalogue (PHYSICS_INTRO, CHANGELOG, IMPROVEMENT_GUIDE)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact.
- [ ] `git log --format=%B dev/examples-catalogue | grep -i
      'co-authored-by: claude'` returns empty.
- [ ] No schema change (v2.7.0 stays current).
- [ ] No `semi/` source change (Layer 2-4 untouched). Only
      `scripts/run_benchmark.py` (verifier registrations + plot
      helpers) and `tests/test_examples_register.py` are
      Python-touched.
- [ ] Every existing benchmark is bit-identical to v0.23.0
      (the examples PR touches no existing config).
- [ ] Coverage gate holds at 95 in the gated `docker-fem-tests`
      job. The new `tests/test_examples_register.py` adds three
      pytests; the coverage delta is bounded.
- [ ] Each example has a `README.md` matching the load-bearing
      structure convention.
- [ ] Each example has a CI matrix entry; none is marked
      `allow-failure: "true"`.

## Anti-goals

- Do not bump the schema. Examples use only v2.7.0 features.
- Do not start M17 (heterojunctions) or M19 (3D MOSFET) in
  this PR.
- Do not add a fourth example. Three is enough; the fourth
  candidates listed in `examples/README.md` "future examples"
  are explicit follow-up scope.
- Do not assert tight numerical correctness in the example
  verifiers. Examples are illustrative; tight tolerances make
  parameter tuning brittle. Smoke checks only.
- Do not refactor benchmarks/ to share infrastructure with
  examples/. The two directories serve different purposes;
  shared infrastructure would be premature unification.
- Do not add a tunnel-diode example. Forward-bias BBT under
  M16.6 has not been independently validated; ship a tunnel-
  diode example only after a future PR confirms the M16.6
  Kane formula handles forward-bias band-to-band tunneling
  cleanly.
- Do not ship optional notebooks for the examples. Notebooks
  are M16.x convention for benchmarks; examples have
  README.md as their narrative documentation.
- Do not retire the `mosfet_2d` `allow-failure: "true"` flag.
- Do not add `Co-Authored-By` trailers, generated-by footers,
  or any other AI-credit marker.

## Stop conditions

You are done when:

1. The PR is opened on branch `dev/examples-catalogue`.
2. All three example smoke verifiers pass in CI:
   - `examples/nmos_idvgs/`
   - `examples/schottky_iv_temperature/`
   - `examples/power_diode_reverse_recovery/`
3. `tests/test_examples_register.py` passes in the gated
   `docker-fem-tests` job; coverage gate holds at 95 without
   a follow-up commit.
4. `examples/README.md` table-of-examples matches what shipped.
5. PLAN.md, IMPROVEMENT_GUIDE.md (§ 9), PHYSICS_INTRO.md (§ 6),
   CHANGELOG.md reflect the examples catalogue ship.
6. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI
   assistant. Squash-merge body is also clean.
7. Every existing benchmark remains bit-identical to v0.23.0
   (cite the four anchors in the PR description: pn_1d_bias
   J(V=0.6 V) = 1.635e+03 A/m^2; diode_velsat_1d 56.27 % @
   0.9 V, 0.19 % @ 0.3 V; schottky_1d worst-case <10 %;
   zener_1d worst-case <20 % Kane match).
8. PR reviewed, CI green (modulo the documented `allow-failure`
   on `mosfet_2d`, unchanged in scope from v0.23.0), merged.

## PR description template

```
## Summary

Examples catalogue: a new `examples/` top-level directory with
three self-contained practical device configs that exercise the
M16 physics catalogue in real engineering contexts. Distinct
from `benchmarks/` (which gates V&V); examples are illustrative
clone-and-modify starting points for users.

Three examples ship in this PR:

- `examples/nmos_idvgs/` — NMOS Id-Vgs at V_DS = 0.05, 0.5, 1.0,
  1.8 V. Exercises Caughey-Thomas mobility (M16.1), Lombardi
  surface mobility (M16.2), and Fermi-Dirac statistics (M16.4)
  for the n+ source/drain.
- `examples/schottky_iv_temperature/` — Pt-on-n-Si Schottky
  diode I-V at T = 250 K, 300 K, 350 K. Exercises thermionic
  emission (M16.5) under temperature variation.
- `examples/power_diode_reverse_recovery/` — rectifier turn-off
  transient using a piecewise-linear voltage_t. Exercises
  voltage_t (M16.7) and Auger recombination (M16.3) at high
  injection.

No schema bump (v2.7.0 stays current). No `semi/` source change
beyond verifier registrations and plot helpers in
`scripts/run_benchmark.py`. Every existing benchmark is
bit-identical to v0.23.0.

## Acceptance

- [ ] examples/nmos_idvgs/ smoke verifier passes
- [ ] examples/schottky_iv_temperature/ smoke verifier passes
- [ ] examples/power_diode_reverse_recovery/ smoke verifier
      passes
- [ ] Existing benchmarks bit-identical to v0.23.0
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
      schottky_1d worst-case <10 %; zener_1d worst-case
      <20 % Kane match)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95 (in the gated
      docker-fem-tests job)
- [ ] docker compose run --rm benchmark nmos_idvgs
- [ ] docker compose run --rm benchmark schottky_iv_temperature
- [ ] docker compose run --rm benchmark power_diode_reverse_recovery
- [ ] docker compose run --rm benchmark pn_1d_bias
      (existing benchmark, bit-identical)

## Notes

- This PR does not advance the M-numbered roadmap. Next-task
  pointer (M17 heterojunctions or M19 3D MOSFET capstone) is
  unchanged.
- The `examples/` directory is documentation / demonstration
  layer; the V&V layer is `benchmarks/`. Future examples
  candidates (NPN BJT, NMOS C-V at multiple body biases, SiC
  Schottky, tunnel diode forward I-V) are listed in
  `examples/README.md` "future examples" section.
- No package version bump (no API or schema change). The next
  physics PR (M17 or M19) will roll this into its own version.
```

## Hand-off

When the examples catalogue lands, the next pickup is the
maintainer's choice between M17 (heterojunctions) and M19 (3D
MOSFET capstone), per the M16.7 PR's hand-off note.

The examples catalogue itself is extensible: future
contributors can add examples by following the README structure
convention and registering a smoke verifier in
`scripts/run_benchmark.py`. Each future example is its own PR
on `dev/examples-<name>` and does not require a milestone
number.
