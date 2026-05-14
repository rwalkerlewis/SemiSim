## Task: Bias-sweep SNES line-search stabilization (retire the `nmos_idvgs` `allow-failure` carve-out)

You are continuing development on SemiSim, a JSON-driven FEniCSx (dolfinx 0.10)
finite-element semiconductor device simulator. Read `PLAN.md` end to end before
writing any code. The "Current state", "Next task", and "Invariants" sections are
ground truth. This prompt corresponds to the first of the two unblocked Next-task
candidates: the bias-sweep SNES stabilization that unblocks `nmos_idvgs`. The
maintainer has selected this over M19 (3D MOSFET capstone) because it is the
smaller lift, it clears the last `allow-failure: "true"` CI carve-out outside
`mosfet_2d`, and a stable bias-sweep inner solve is a sane prerequisite for the
M19 capstone anyway.

### Problem statement

The `examples/nmos_idvgs` config sweeps V_GS across the MOSFET inversion onset
under Fermi-Dirac statistics. The run is tagged `allow-failure: "true"` in the CI
matrix (introduced in `ed6719b`) because the SNES inner solve stagnates at the
threshold: the default `bt` backtracking line search cannot find a descent
direction where the FD prefactor makes the Jacobian curvature steep at
depletion-to-inversion onset. The outer V_GS voltage step is already handled by
`semi.continuation.AdaptiveStepController` (the same class driving the
forward-bias ramp); the failure is in the SNES inner solve, not the outer step,
so do not look for the fix there.

### Scope and approach

This is an ADR-level decision. Evaluate the candidate fixes named in `PLAN.md`:

1. Switch the SNES line search from `bt` to PETSc `nleqerr` or `cp`.
2. Add an explicit damping schedule on the Newton update.
3. Introduce a homotopy parameter on the FD prefactor (continuation from
   Boltzmann to Fermi-Dirac at fixed bias).

Prototype the cheapest option first (line-search type is a one-line PETSc option
change) and only escalate to damping or homotopy if the line-search swap does not
converge `nmos_idvgs` across the full V_GS sweep. Whatever you land on, document
the decision and the rejected alternatives in a new ADR under `docs/adr/`
(next number in sequence; ADR 0017 was M18, so this is ADR 0018). Cross-reference
ADR 0008 (the M12 SNES tolerance amendment) since this touches the same solver
path.

This is a solver-stability change, not a new physics kernel, so it ships **no new
MMS variant** (ADR 0006's per-physics-module MMS rule applies to physics kernels;
the precedent is M18, which shipped audit-only V&V for the same reason). V&V for
this task is: (a) `nmos_idvgs` runs to completion across the full sweep with a
verifier that gates qualitative I-V shape, and (b) byte-identity on every existing
benchmark.

### Hard constraints (verify before delivering)

- **Byte-identity on every existing benchmark.** The line-search or damping change
  must not alter any converged result. Anchor: `pn_1d_bias` J(V=0.6 V) =
  1.635e+03 A/m^2. The full anchor set that must continue to hold: `diode_velsat_1d`,
  `diode_auger_1d`, `diode_fermi_dirac_1d`, `schottky_1d`, `zener_1d`, `pn_1d_turnon`,
  `pn_1d_pulse`, `diode_sine_1d`, `rc_ac_sweep`, `algaas_gaas_heterojunction`. If a
  new line-search type changes a converged Newton trajectory enough to perturb the
  last-digit result, that is a regression: gate the new behavior behind a config
  flag defaulting to the current `bt` behavior, or prove the perturbation is below
  the SNES atol floor and document it.
- **dolfinx 0.10 API only.** Use `dolfinx.fem.petsc.NonlinearProblem` with
  `petsc_options_prefix`. Do not introduce `dolfinx.nls.petsc.NewtonSolver`
  (deprecated, ADR 0003). SNES options go through the PETSc options database on
  the existing prefix.
- **Pure-Python core stays dolfinx-free** (ADR 0007). If you add a config knob,
  the schema-side validation lives in `semi/schema.py` and must import without
  dolfinx.
- **Slotboom primary variables** (ADR 0004). No SUPG or streamline-diffusion
  stabilization of raw continuity equations. The fix is in the nonlinear solver
  strategy, not the weak form.
- **No em dashes anywhere in repo prose** (Invariant 8).
- **Physics-style variable names are expected and Ruff-accepted** (Invariant 7).
  Do not PEP-8-rename `V_GS`, `N_A`, etc.

### Schema

If the fix exposes a user-facing knob (for example `solver.line_search` with an
enum, or `solver.damping` with a schedule), that is an additive minor schema bump
v2.9.0 -> v2.10.0: enum widening or a new optional sub-object, every new node
annotated with a UI-facing description, `additionalProperties: false` preserved,
default value chosen so configs without the field are bit-identical to v0.25.0.
Bump `SCHEMA_SUPPORTED_MINOR` accordingly and add schema-side assertions to a
`tests/test_*_schema.py` file covering default-fill, the new field's validation,
and v2.0.0-through-v2.9.0 forward compatibility.

If the fix is purely internal (a hard-coded line-search change with no user knob),
there is no schema change. Decide explicitly in the ADR and state which path you
took.

### Phase structure

All work happens on a new branch off `main`, never on `main` directly. Create
`dev/m18.1-bias-sweep-snes-stabilization` (or the maintainer's preferred number;
M18.1 is the natural slot) before Phase 0. Follow the repo's phase-letter commit
convention: one commit per phase, phase letter in the commit message. The branch
is merged via a single PR (see Phase F); do not push to `main`.

- **Phase 0:** Ship this starter prompt verbatim on the branch. Add an
  `[Unreleased]` entry to `docs/IMPROVEMENT_GUIDE.md` section 9.
- **Phase A (ADR + diagnosis):** Write ADR 0018. Reproduce the `nmos_idvgs`
  stagnation locally under Docker, capture the SNES convergence log at the failing
  V_GS step, and record in the ADR which line-search / damping / homotopy option
  resolves it and why the others were rejected.
- **Phase B (solver change):** Implement the chosen fix in the bias_sweep runner /
  SNES setup path. If a schema knob is added, wire it here and bump the schema
  (Phase B-schema folded in, or split into its own phase if cleaner).
- **Phase C (verifier + byte-identity proof):** Give `examples/nmos_idvgs` a real
  verifier gating qualitative I-V shape (monotonic increase past threshold,
  positive transconductance in the inversion window). Re-run every benchmark and
  confirm the anchor set above is byte-identical. Add the byte-identity check to
  the test suite if it is not already covered.
- **Phase D (CI):** Remove `allow-failure: "true"` from the `nmos_idvgs` matrix
  entry in `.github/workflows/ci.yml`. The `mosfet_2d` `allow-failure` flag stays
  (that is a separate, independently-tracked SNES depletion-onset issue; do not
  touch it). Confirm the coverage gate at 95 still holds; add unit tests for any
  new code path (line-search dispatch, damping schedule, homotopy loop) to keep
  gated coverage at 95 without a follow-up commit.
- **Phase E (closeout):** Update PLAN.md: move this task's summary to the
  append-only "Completed work log" (newest on top), rewrite "Current state" to the
  new reality, set "Next task" to M19 (3D MOSFET capstone, now the sole remaining
  unblocked candidate). Update `docs/ROADMAP.md`, `CHANGELOG.md` (`[0.26.0]`
  entry), and bump `pyproject.toml` and `semi/__init__.py` 0.25.0 -> 0.26.0. Mark
  the IMPROVEMENT_GUIDE entry Done with the observed converged V_GS-sweep result.
- **Phase F (PR):** Push the branch and open a pull request against `main`. Title
  it `M18.1: bias-sweep SNES line-search stabilization`. The PR description
  summarizes each phase, names ADR 0018 and the line-search decision, lists the
  schema bump (or states "no schema change"), states the new package version
  0.26.0, and confirms the acceptance test below is met. Do not merge; leave the
  PR open for maintainer review.

### Acceptance test

`nmos_idvgs` runs to completion across the full V_GS sweep in CI with no
`allow-failure` flag, its verifier passes, every benchmark in the anchor set is
byte-identical to v0.25.0, the coverage gate holds at 95, and ADR 0018 documents
the line-search decision with the rejected alternatives. Do not negotiate around
the byte-identity requirement: if the fix perturbs a converged result, either gate
it behind a default-off flag or prove the perturbation is sub-atol and document
the proof.
