# M16.1 starter prompt: Caughey-Thomas field-dependent mobility

Paste-ready prompt for Claude in VS Code. Mirrors the convention of
`docs/M9_STARTER_PROMPT.md`, `docs/M15_STARTER_PROMPT.md`, and
`docs/M14_3_STARTER_PROMPT.md`. Pick this up on a fresh branch
`dev/m16.1-caughey-thomas` after M14.3 has merged. Do not bundle
with any other milestone work.

---

You are working in `SemiSim` at v0.16.0 (post-M14.3). Your
assignment is **M16.1: Caughey-Thomas field-dependent mobility**, the
first physics-completeness slice of the M16 umbrella.

This prompt does not restate the M16.1 deliverable, the acceptance
tests, or the rationale; those live in `docs/IMPROVEMENT_GUIDE.md`
§ M16.1. Read the guide first; this prompt only tells you the order
in which to execute and which invariants must remain load-bearing.

## Required reading (do not skip; ~45 minutes)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.16.0 (post-M14.3) and
   that no other PR is in flight on `dev/m16.1-*`. The "Next task"
   section names M16.1; if it does not, stop and align with the
   maintainer before continuing.
2. `docs/IMPROVEMENT_GUIDE.md` § M16 (the umbrella context),
   § M16.1 (the milestone definition with Why / Deliverable /
   Acceptance / Dependencies), and § 1 (Honest current state).
3. `docs/ROADMAP.md` § Capability matrix, § Honest gap, and the
   M16.1 row.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M16.1 touches
   physics (Layer 4), the schema (Layer 2), and the benchmark
   (Layer 5). Layer 3 (pure-Python core) and the BC layer are not
   touched.
5. `docs/PHYSICS.md` for the scaled drift-diffusion derivation.
   The Caughey-Thomas mobility is a coefficient on the scaled
   continuity flux; understand which scale absorbs the field
   dependence.
6. `docs/adr/0001`, `0002`, `0004` (Slotboom variables for DD;
   M16.1 must remain in Slotboom form), `0006` (V&V strategy;
   M16.1 needs an MMS verifier), `0007` (BC interface; do not
   touch).
7. `docs/mms_dd_derivation.md` for the existing MMS-DD harness; the
   Caughey-Thomas variant extends it with a gradient-dependent
   mobility.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** The new schema entries
  `physics.mobility.model: "constant" | "caughey_thomas"` and the
  parameters `vsat_n`, `vsat_p`, `beta_n`, `beta_p` must be
  expressible in `schemas/input.v2.json` (M14.3 has bumped to v2),
  validated by `semi/schema.py`, and exercised by at least one
  benchmark JSON.
- **Schema versioning is binding.** This is an additive minor bump
  v2.0.0 to v2.1.0; v2.0.0 inputs (no `physics.mobility` block)
  must continue to validate and produce bit-identical results.
  Update `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Five layers, enforced.** The new `semi/physics/mobility.py` is
  Layer 4 (FEM); imports for `dolfinx`, `ufl`, `petsc4py` go
  inside function bodies, not at module scope. The mobility
  parameters live in the JSON schema (Layer 1) and the schema
  loader (Layer 2); pure Python.
- **Slotboom variables for DD.** ADR 0004 is locked. The
  Caughey-Thomas mobility multiplies the Slotboom flux; do not
  re-derive the continuity rows in primary-density form.
- **Every new physics module needs an MMS verifier.** ADR 0006.
  No exceptions. Acceptance L2 rate >= 1.99 and H1 rate >= 0.99 at
  the finest pair, gated in `scripts/run_verification.py mms_dd`.
- **Every new analytical model needs a benchmark.** The
  `diode_velsat_1d` benchmark is the analytical anchor.
- **One milestone, one PR.** Do not bundle M16.1 with M16.2 or any
  other physics work.
- **No em dashes in prose or code comments.** Use commas, periods,
  parentheses, or colons.
- **Physics-style variable names are allowed.** `mu_n0`, `vsat_n`,
  `beta_n`, `F_par`, etc.
- **Coverage gate is 95 on `semi/`.** The new mobility module must
  carry pure-Python unit tests for the closed-form math; the FEM
  wiring is covered by the MMS test and the benchmark.

## Five phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/` and
`pytest tests/`. Do not advance with red tests.

---

### Phase A: schema surface for the mobility dispatch

Pure-Python only. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.0.0 to 2.1.0. Confirm the
   major-version gate in `semi/schema.py` still accepts 2.1.0;
   bump `SCHEMA_SUPPORTED_MINOR` accordingly.
2. Add to `physics`:
   ```jsonc
   "mobility": {
     "type": "object",
     "additionalProperties": false,
     "properties": {
       "model": {
         "enum": ["constant", "caughey_thomas"],
         "default": "constant",
         "description": "Mobility model dispatch. constant is the V&V reference and the v0.15.0 baseline; caughey_thomas adds closed-form velocity saturation."
       },
       "vsat_n": {"type": "number", "default": 1.0e7,
                  "description": "Electron saturation velocity (cm/s). Default is the standard Si value."},
       "vsat_p": {"type": "number", "default": 8.0e6,
                  "description": "Hole saturation velocity (cm/s)."},
       "beta_n": {"type": "number", "default": 2.0,
                  "description": "Electron Caughey-Thomas exponent."},
       "beta_p": {"type": "number", "default": 1.0,
                  "description": "Hole Caughey-Thomas exponent."}
     }
   }
   ```
3. Default-fill: an input with no `physics.mobility` block produces
   an exact bit-equivalent solve to v0.16.0 (which is bit-equivalent
   to v0.15.0 plus the M14.3 housekeeping closures).
4. Pure-Python tests in `tests/test_mobility_schema.py`:
   - Every existing benchmark JSON validates unchanged against
     v2.1.0.
   - `physics.mobility.model: "caughey_thomas"` validates with
     defaults.
   - An unknown model name (`"foo"`) is rejected.
   - An extra property under `physics.mobility` (strict mode) is
     rejected.

**Commit message:** `feat(schema): physics.mobility dispatch (schema 2.1.0); constant default unchanged (M16.1)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change is observable on any benchmark.

---

### Phase B: the mobility builder

The closed-form Caughey-Thomas expression as a UFL builder. Layer 4;
imports inside function bodies.

1. New module `semi/physics/mobility.py`. Public API:
   ```python
   def caughey_thomas_mu(mu0, F_par, vsat, beta):
       """
       Caughey-Thomas field-dependent mobility:
           mu(F) = mu0 / (1 + (mu0 * F_par / vsat)**beta)**(1/beta)
       Inputs are UFL expressions or Constants in scaled units;
       returns a UFL expression in the same scaled units.
       """
   def constant_mu(mu0):
       """Identity wrapper; returns mu0 unchanged."""
   def build_mobility(physics_cfg, mu0_n, mu0_p, F_par_n, F_par_p):
       """Dispatch on physics.mobility.model and return (mu_n, mu_p) UFL expressions."""
   ```
2. `F_par` is the magnitude of the field component parallel to
   the current vector. For drift-diffusion in Slotboom form, the
   carrier-specific parallel field magnitude is what
   Caughey-Thomas wants; document the sign convention and the
   absolute-value choice in a docstring.
3. Pure-Python unit tests in
   `tests/test_mobility_closed_form.py`:
   - Limit `F_par -> 0`: `mu(F) -> mu0` to within 1e-12.
   - Limit `F_par -> inf`: `mu(F) * F_par -> vsat` to within 1%.
   - `beta = 2`, intermediate field: matches the closed form by
     hand calculation at three sample points.

**Commit message:** `feat(physics): Caughey-Thomas mobility UFL builder (M16.1)`

**Acceptance:** new tests pass; no FEM behavior change.

---

### Phase C: wire the builder into the runners

Add the dispatch hook in the drift-diffusion form builder and the
bias-sweep runner. Touch only what is needed.

1. `semi/physics/drift_diffusion.py`: thread the mobility dispatch
   through `build_dd_block_residual` and the multi-region variant.
   The `constant` branch is bit-identical; the
   `caughey_thomas` branch substitutes
   `caughey_thomas_mu(mu0, F_par, vsat, beta)` for `mu0` in the
   continuity flux.
2. `semi/runners/bias_sweep.py`: at runner setup, build
   `(mu_n, mu_p)` via `semi.physics.mobility.build_mobility` and
   pass them to the DD form builder.
3. The other runners (equilibrium, mos_cv, mos_cap_ac, transient,
   ac_sweep) do not need touching for this PR; the equilibrium and
   MOSCAP runners do not have continuity rows under bias, and the
   transient and AC runners can adopt mobility dispatch in a
   follow-up.
4. Existing benchmarks: with `physics.mobility.model: "constant"`
   (the default), every benchmark must produce results
   bit-identical to v0.16.0. Run the V&V suite and the benchmark
   matrix to confirm.

**Commit message:** `feat(runners): wire Caughey-Thomas mobility into bias_sweep and DD form (M16.1)`

**Acceptance:** every existing benchmark with the default
mobility produces bit-identical numerics to v0.16.0; CI green.

---

### Phase D: MMS verifier with a gradient-dependent mobility

ADR 0006 mandates an MMS verifier for every new physics module.
Extend the existing MMS-DD harness.

1. New file `tests/fem/test_mms_caughey_thomas.py`:
   - Manufactured solution: `psi_e = sin(k*x)`, `phi_n_e = a*sin(k*x)`,
     `phi_p_e = b*sin(k*x)` (or the existing MMS-DD reference if it
     fits). The mobility depends on `|grad(phi_n_e)|` and
     `|grad(phi_p_e)|` so the gradient is non-trivial.
   - Forcing terms in UFL (do not form `ufl.div(ufl.grad(...))`
     directly; the constant prefactor on a coarse mesh can collapse
     to numerical zero).
   - Mesh sweep with N in [50, 100, 200, 400].
   - Acceptance: finest-pair L2 rate >= 1.99 and H1 rate >= 0.99 on
     each block (psi, phi_n, phi_p).
2. Extend `semi/verification/mms_dd.py` to expose the new variant
   so `python scripts/run_verification.py mms_dd` includes it.
3. Update `docs/PHYSICS.md` § Verification & Validation with a
   one-paragraph note on the new MMS variant and the rate gate.
4. Update `docs/mms_dd_derivation.md` with the manufactured-source
   derivation for the Caughey-Thomas variant.

**Commit message:** `feat(verification): MMS-DD Caughey-Thomas variant (M16.1)`

**Acceptance test 1**: `python scripts/run_verification.py mms_dd`
includes the caughey_thomas variant and reports L2 rate >= 1.99 at
the finest pair.

---

### Phase E: the diode_velsat_1d benchmark

The analytical anchor: a 1D pn diode where Caughey-Thomas and
constant mobility predict measurably different I-V at high field
and converge at low field.

1. New directory `benchmarks/diode_velsat_1d/`:
   - `diode_velsat.json`: 1D pn diode at V_F in [0.5, 0.9] V, otherwise
     identical to `pn_1d_bias`. Two configs (or one config with a
     parametric switch): `physics.mobility.model: "constant"` and
     `physics.mobility.model: "caughey_thomas"`.
   - `README.md`: device description, analytical reasoning for the
     >5% divergence at V = 0.9 V and the <1% convergence at V = 0.5 V.
2. New verifier in `scripts/run_benchmark.py`:
   `verify_diode_velsat_1d`:
   - Run the constant-mu sweep, run the caughey-thomas sweep.
   - Assert `|I_CT(0.9) - I_const(0.9)| / I_const(0.9) > 0.05`.
   - Assert `|I_CT(0.5) - I_const(0.5)| / I_const(0.5) < 0.01`.
3. Wire the new benchmark into CI: add a step in
   `.github/workflows/ci.yml` (or the docker-fem job) that runs
   `python scripts/run_benchmark.py diode_velsat_1d`.
4. Notebook (optional but recommended for parity with shipped
   benchmarks): `notebooks/06_diode_velsat_1d.ipynb` showing the
   two I-V curves and the divergence plot.

**Commit message:** `feat(benchmark): diode_velsat_1d Caughey-Thomas vs constant mobility (M16.1)`

**Acceptance test 2**: `python scripts/run_benchmark.py
diode_velsat_1d` exits 0 with the divergence-vs-convergence
verifier passing.

---

## Phase F: close out

1. `PLAN.md`:
   - Move M16.1 from "Next task" to "Completed work log" with the
     PR number, deliverables, schema minor bump, and acceptance-test
     results.
   - Set "Next task" to **M16.2 Lombardi surface mobility** with a
     pointer to the (yet-to-be-authored) M16.2 starter prompt.
   - Refresh "Current state" with the new package version (bump
     0.16.0 to 0.17.0 for the schema minor and the new physics
     model).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M16.1 Done in § 4 with a one-line summary and a CHANGELOG
     anchor.
   - Append a § 9 changelog entry.
3. `docs/ROADMAP.md`:
   - Update the M16.1 row in the capability matrix from Planned to
     shipped.
4. `CHANGELOG.md`:
   - New `[0.17.0]` entry with the M16.1 line items.
5. Push every commit to origin immediately. Include the SHA, log
   line, and `git status` snapshot in each gate report.

**Commit message:** `docs: close out M16.1 (PLAN, IMPROVEMENT_GUIDE, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by this
      PR.
- [ ] Pure-Python core remains dolfinx-free.
- [ ] Constant-mobility path is bit-identical to v0.16.0 on every
      benchmark.
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004).
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `semisim_server` public API.
- [ ] Schema bumped per minor (2.0.0 to 2.1.0); v2.0.0 inputs still
      validate.
- [ ] MMS rate gate L2 >= 1.99 and H1 >= 0.99 active in
      `scripts/run_verification.py`.
- [ ] Coverage gate holds at 95.

## Anti-goals

- Do not start M16.2 (Lombardi) in this PR. Surface mobility is its
  own PR with its own MMS variant.
- Do not change anything in `semi/bcs.py`, `semi/scaling.py`, or
  `semi/runners/equilibrium.py`. The mobility dispatch goes through
  the bias_sweep runner and the DD form builder.
- Do not introduce a `mu_n` / `mu_p` field as a new primary unknown.
  Caughey-Thomas is closed-form; it does not add unknowns.
- Do not lower the MMS rate gate to "passes on the coarse mesh".
  The gate is L2 >= 1.99 and H1 >= 0.99 at the finest pair, every
  block. If a rate degrades, the implementation has a bug.
- Do not bundle this PR with the M14.2.x backlog or with M16.2
  through M16.7. Those are separate PRs.

## Stop conditions

You are done with this prompt when:

1. The M16.1 PR is opened on branch `dev/m16.1-caughey-thomas`.
2. Both acceptance tests in `docs/IMPROVEMENT_GUIDE.md` § M16.1
   are green in CI.
3. Acceptance test 3 (every existing benchmark with constant
   mobility produces bit-identical results to v0.16.0) is gated by
   the existing benchmark CI matrix; confirm in the PR description
   that the matrix passed.
4. PLAN.md, IMPROVEMENT_GUIDE.md, ROADMAP.md, and CHANGELOG.md
   reflect the close-out.
5. The package version has bumped 0.16.0 to 0.17.0 in
   `pyproject.toml` and `semi/__init__.py`.

## PR description template

```
## Summary

M16.1: Caughey-Thomas field-dependent mobility, the first
physics-completeness slice of the M16 umbrella. Closed-form,
no new unknowns, no change to Slotboom primary form.

Schema additive minor bump v2.0.0 to v2.1.0
(physics.mobility.model dispatch, default "constant"). Constant
branch is bit-identical to v0.16.0.

New benchmark: diode_velsat_1d demonstrates the >5% divergence
at V_F = 0.9 V and <1% convergence at V_F = 0.5 V between
constant and Caughey-Thomas mobility on a 1D pn diode.

New MMS variant: gradient-dependent mobility manufactured
solution. Finest-pair L2 rate >= 1.99 and H1 rate >= 0.99 on
each block.

## Acceptance tests

(Both are in docs/IMPROVEMENT_GUIDE.md § M16.1.)

- [x] A1: scripts/run_verification.py mms_dd caughey_thomas
      variant L2 >= 1.99 finest-pair
- [x] A2: scripts/run_benchmark.py diode_velsat_1d divergence
      and convergence verifier passes
- [x] A3: every existing benchmark with constant mobility is
      bit-identical to v0.16.0

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark diode_velsat_1d
- [ ] docker compose run --rm benchmark mosfet_2d
      (Pao-Sah, M14.3 verifier; constant mu still passes)
- [ ] docker compose run --rm benchmark pn_1d_bias
      (constant mu, bit-identical)
```

## Hand-off

When M16.1 lands and is reviewed, the next pickup is M16.2
(Lombardi surface mobility), authored at the time by the
contributor in the same shape as this prompt. Subsequent
milestones (M16.3 through M16.7, M19, M19.1, M20) each ship as
their own PR with their own starter prompt; the shape is
established by `docs/M14_3_STARTER_PROMPT.md` and this file.
