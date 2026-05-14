# M14.3 starter prompt: Housekeeping (cheap closes)

Paste-ready prompt for Claude in VS Code. Mirrors the convention of
`docs/M9_STARTER_PROMPT.md` and `docs/M15_STARTER_PROMPT.md`. Pick this
up on a fresh branch `dev/m14.3-housekeeping`. Do not bundle with any
other milestone work.

---

You are working in `SemiSim` at v0.15.0. Your assignment is
**M14.3: Housekeeping (cheap closes)**, the bridge PR between M15
(GPU linear-solver path) and M16.1 (Caughey-Thomas mobility, the
first physics-completeness slice).

This prompt does not restate the M14.3 deliverable, the acceptance
tests, or the rationale; those live in `docs/IMPROVEMENT_GUIDE.md`
§ M14.3. Read the guide first; this prompt only tells you the order
in which to execute and which invariants must remain load-bearing.

> **Note:** Phase A was struck after PR #68 landed (the original
> observation was a cached-render artifact in the external
> reviewer's tooling, not real staleness; see CHANGELOG `[Unreleased]`
> for the full retraction). The original phase lettering (B, C, D, E)
> is preserved below for traceability with PR #68's history.

## Required reading (do not skip; ~30 minutes)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.15.0 and that no other
   PR is in flight on `dev/m14.3-*`. The "Next task" section names
   M14.3; if it does not, stop and align with the maintainer before
   continuing.
2. `docs/IMPROVEMENT_GUIDE.md` § M14.3 (the milestone definition,
   including the Why / Deliverable / Acceptance / Dependencies block)
   and § 1 (Honest current state, for the v0.15.0 baseline).
3. `docs/ROADMAP.md` § Capability matrix and § Honest gap. M14.3 is
   listed as Planned.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M14.3 touches the
   schema (Layer 2), the FEM mesh ingest (Layer 4), and the
   benchmark / verifier (Layer 5). Layer 3 (pure-Python core) and
   Layer 4 physics kernels stay untouched.
5. `docs/adr/0001`, `0002`, `0006`, `0011`, `0012`. ADR 0006
   (V&V strategy) and ADR 0012 (SG primitives status) are directly
   relevant: the Pao-Sah verifier is a new analytical check, and the
   SG primitives go away in this PR.
6. `docs/PHYSICS.md` Section 6 (MOS reference). The new mosfet_2d
   verifier window choice belongs here.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** Every shipped feature must be expressible
  in `schemas/input.v1.json` (and v2.0.0 from this PR onward),
  validated by `semi/schema.py`, and exercised by at least one
  benchmark JSON.
- **Schema versioning is binding.** This PR ships a major bump
  (v1 to v2) because flipping `additionalProperties: false` is a
  schema-breaking change. Both schemas must be present in the tree,
  the engine must accept both for one minor cycle, and v1 inputs
  must log a deprecation warning.
- **Five layers, enforced.** The pure-Python core
  (`semi/constants.py`, `semi/materials.py`, `semi/scaling.py`,
  `semi/doping.py`, `semi/schema.py`, `semi/cv.py`,
  `semi/timestepping.py`, `semi/bcs.py`, `semi/results.py`,
  `semi/diode_analytical.py`, `semi/continuation.py`,
  `semi/compute.py`) does not import dolfinx, ufl, petsc4py, or
  mpi4py at module scope. The XDMF branch in `semi/mesh.py` is
  Layer 4; imports go inside the function body.
- **Every new analytical model needs a benchmark.** The Pao-Sah
  reference for mosfet_2d is the new analytical model; its place in
  the verifier is the benchmark.
- **One milestone, one PR.** Do not bundle M14.3 with M16.1 or any
  other physics work.
- **No em dashes in prose or code comments.** Use commas, periods,
  parentheses, or colons.
- **Physics-style variable names are allowed.** `N_A`, `V_t`,
  `psi_L`, `eps_Si`, etc. PEP 8 N rules are off globally in ruff.
- **Coverage gate ends this PR at 95.** It enters this PR at 92;
  the SG-removal step plus the new XDMF / Pao-Sah / strict-schema
  unit tests should restore it. Do not lower the gate; if the
  gate drops below 95, add unit tests rather than relax the
  threshold.

## Four phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/` and
`pytest tests/`. Do not advance with red tests.

---

### Phase B: tighten the mosfet_2d verifier with a Pao-Sah reference

Bring the `mosfet_2d` benchmark from a qualitative check to an
analytical-reference check. The Pao-Sah square-law model is the
correct linear-regime reference at low V_DS.

1. Add a Pao-Sah / square-law analytical reference to
   `scripts/run_benchmark.py`'s `verify_mosfet_2d`:
   - Linear-regime current
     `I_D = (W/L) * mu_n * C_ox * (V_GS - V_T) * V_DS` for
     `V_GS > V_T`, `V_DS << V_GS - V_T`.
   - Use the existing `semi/cv.py` MOSCAP analytical helpers for
     `V_T`; if they do not export the helper you need, extend
     `semi/cv.py` (Layer 3, pure-Python; no dolfinx imports).
2. Pick the verifier window: `V_DS = 0.05 V`, `V_GS in [V_T + 0.2, V_T + 0.6] V`.
   Tolerance: 20%.
3. Document the choice in `benchmarks/mosfet_2d/README.md` with a
   one-paragraph derivation. The window edges are: lower bound
   so we are above subthreshold (Boltzmann tail does not dominate);
   upper bound so velocity-saturation is still negligible
   (M16.1 will widen the window once Caughey-Thomas ships).
4. Document the verifier window choice in `docs/PHYSICS.md`
   Section 6 (the existing MOS reference). One short paragraph
   plus the equation; cross-reference the README.
5. Add a new pure-Python unit test in `tests/test_mosfet_2d_verifier.py`
   that mocks the IV result and asserts the verifier passes when
   the simulated I_D is within 20% and fails when outside.

**Acceptance test 2**: `python scripts/run_benchmark.py mosfet_2d`
exits 0 with the new analytical-reference verifier passing inside
the [V_T + 0.2, V_T + 0.6] V window at the documented 20%
tolerance.

**Commit message:** `feat(benchmark): mosfet_2d Pao-Sah analytical verifier (M14.3)`

---

### Phase C: implement XDMF mesh ingest

Wire the existing `NotImplementedError` branch in
`semi/mesh.py::_build_from_file`.

1. In `semi/mesh.py::_build_from_file`, add the XDMF branch:
   - Use `dolfinx.io.XDMFFile.read_mesh(comm, "r", path)` to load
     the mesh.
   - Use `XDMFFile.read_meshtags(...)` to load cell tags and facet
     tags. Propagate them as `cell_tags` and `facet_tags` matching
     the gmsh `.msh` branch.
   - Imports inside the function body, not at module scope.
2. Test: `tests/fem/test_mesh_xdmf.py` round-trips the resistor
   benchmark mesh through `box.msh` -> `box.xdmf` via meshio.
   Compute R from each load path; assert
   `abs(R_xdmf - R_msh) / R_msh < 1e-12`.
3. Update `docs/schema/reference.md` to note that
   `mesh.source: "file"` accepts both `.msh` and `.xdmf`.

**Acceptance test 3**: the resistor benchmark loaded from
`box.xdmf` produces R within 1e-12 relative of the same benchmark
loaded from `box.msh`.

**Commit message:** `feat(mesh): XDMF mesh ingest with cell-tag propagation (M14.3)`

---

### Phase D: strict-mode the input schema (v1 -> v2 major bump)

This is the largest single change in M14.3 and the only one that
breaks v1 compatibility. Do it with care; both v1 and v2 must coexist
for one minor cycle.

1. Copy `schemas/input.v1.json` to `schemas/input.v2.json`. In v2:
   - Set `additionalProperties: false` on every object node.
   - Bump `schema_version` constraint: `"const": "2.0.0"` (or a
     pattern that accepts only 2.x.y).
2. Keep `schemas/input.v1.json` unchanged for one minor cycle.
3. Update `semi/schema.py`:
   - `ENGINE_SUPPORTED_SCHEMA_MAJOR` now lists both 1 and 2 (or
     accepts a range).
   - Validation dispatches on `schema_version` major; v1 inputs log
     a `DeprecationWarning` ("Schema v1 is deprecated; migrate to
     v2.0.0") and continue.
   - `SCHEMA_SUPPORTED_MINOR` resets to 0 for the v2 line.
4. Migrate every benchmark JSON in `benchmarks/` to v2.0.0:
   - Bump `schema_version` to `"2.0.0"`.
   - Verify each one validates as strict v2 (no extra properties).
   - Fix any latent typos or unused fields the strict gate
     surfaces.
5. Add a unit test `tests/test_schema_strict.py`:
   - Every benchmark JSON validates as strict v2.
   - An intentional typo (`"voltag"` instead of `"voltage"`) in a
     contact entry is rejected with a clear error message that
     names the offending field.
   - A v1 input loads with a `DeprecationWarning`.
6. Update `docs/schema/reference.md` to document the v1 / v2
   coexistence and the deprecation timeline.

**Acceptance test 4**: every benchmark JSON in `benchmarks/`
validates as `schema_version: "2.0.0"` with `additionalProperties:
false` active; an intentional typo (`"voltag"` for `"voltage"`) is
rejected with a clear error message that names the offending field.

**Commit message:** `feat(schema): strict-mode v2.0.0 with additionalProperties false (M14.3)`

---

### Phase E: remove dead Scharfetter-Gummel primitives

ADR 0012 documented the SG primitives' status; M13.1 closed in
v0.14.1 with the Slotboom transient (ADR 0014) leaving the SG
primitives unreachable on the active code path. Delete them.

1. Delete `semi/fem/sg_assembly.py` (~792 LOC).
2. Delete `tests/fem/test_sg_assembly.py` if present (the
   dedicated SG test).
3. Confirm via grep that no other module imports
   `semi.fem.sg_assembly` or its public symbols. If any test or
   runner still references it, that is a sign the file is not
   actually dead; stop and reconcile with the maintainer before
   deletion.
4. Add a one-line note in ADR 0012's status section:
   "Reachable in git history at commit `<sha-of-deletion>` if
   M13.2 ever revives the SG flux path."
5. Raise the coverage gate in `pyproject.toml` from 92 to 95
   (back to the M5-era number).
6. Run `pytest --cov=semi --cov-fail-under=95` to confirm the
   gate holds. If it does not, the deletion exposed a real
   coverage gap; add unit tests rather than relax the gate.

**Acceptance test 5**: `semi/fem/sg_assembly.py` and its
dedicated test are gone, the coverage gate in `pyproject.toml`
is 95, and CI is green.

**Commit message:** `chore(fem): remove dead SG primitives; raise coverage gate to 95 (M14.3)`

---

## Phase F: close out

1. `PLAN.md`:
   - Move M14.3 from "Next task" to "Completed work log" with the
     PR number, deliverables, schema major bump, and acceptance-test
     results.
   - Set "Next task" to **M16.1 Caughey-Thomas field-dependent
     mobility** with a pointer to `docs/M16_1_STARTER_PROMPT.md`.
   - Refresh "Current state" with the new package version
     (bump 0.15.0 to 0.16.0; this PR ships the schema major bump,
     so the package version follows).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M14.3 Done in § 4 with a one-line summary and a CHANGELOG
     anchor.
   - Append a § 9 changelog entry.
3. `docs/ROADMAP.md`:
   - Update the M14.3 row in the capability matrix from Planned
     to shipped with the per-deliverable verifier summary.
4. `CHANGELOG.md`:
   - New `[0.16.0]` entry with the M14.3 line items.
5. Push every commit to origin immediately. The maintainer's
   convention is one push per commit (not at the end of the
   branch); include the SHA, log line, and `git status` snapshot in
   each gate report.

**Commit message:** `docs: close out M14.3 (PLAN, IMPROVEMENT_GUIDE, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by this
      PR. Use commas, parentheses, or colons.
- [ ] Pure-Python core remains dolfinx-free.
- [ ] CPU-MUMPS path (default) is bit-identical to v0.15.0 on every
      benchmark; new schema dispatch in v2.0.0 does not change the
      assembled forms.
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `semisim_server` public API.
- [ ] Schema major bumped (v1 to v2) per the schema-versioning
      invariant; both schemas coexist for one minor cycle; v1 inputs
      log a `DeprecationWarning` and still solve.
- [ ] No new ADR introduced unless an invariant is broken; the
      `additionalProperties: false` flip is a schema major bump
      already covered by the existing schema-versioning invariant
      and does not require a new ADR.
- [ ] Coverage gate holds at 95 after the SG removal.

## Anti-goals

- Do not start M16.1 in this PR. The Caughey-Thomas mobility is its
  own PR.
- Do not change anything in `semi/physics/*`, `semi/bcs.py`,
  `semi/scaling.py`, or `semi/runners/*` (other than the schema
  dispatch wiring).
- Do not delete the SG primitives' git history. Removing the file
  is enough; the history is the recovery path if M13.2 ever
  revives it.
- Do not raise the coverage gate above 95 in this PR. That is a
  separate decision.
- Do not bundle this PR with the M14.2.x backlog or with anything
  in M16. Those are separate PRs.

## Stop conditions

You are done with this prompt when:

1. The M14.3 PR is opened on branch `dev/m14.3-housekeeping`.
2. All four acceptance tests in `docs/IMPROVEMENT_GUIDE.md`
   § M14.3 are green in CI.
3. PLAN.md, IMPROVEMENT_GUIDE.md, ROADMAP.md, and CHANGELOG.md
   reflect the close-out.
4. The package version has bumped 0.15.0 to 0.16.0 in
   `pyproject.toml` and `semi/__init__.py`.
5. The schema files `schemas/input.v1.json` and
   `schemas/input.v2.json` both exist; v1 logs a deprecation
   warning; every benchmark JSON is at v2.0.0.
6. `semi/fem/sg_assembly.py` is gone, the coverage gate is 95,
   and CI is green.

## PR description template

```
## Summary

Housekeeping PR (M14.3) closing four small production-hardening
gaps before M16 physics work starts. Doc-only effects on PLAN /
IMPROVEMENT_GUIDE / ROADMAP / CHANGELOG; engine effects bounded
to the schema dispatch, the XDMF branch in semi/mesh.py, and the
mosfet_2d verifier in scripts/run_benchmark.py.

Schema major bump v1 to v2 (additionalProperties: false). Both
schemas coexist; v1 inputs log a DeprecationWarning.

Closes:
- mosfet_2d qualitative verifier (Phase B)
- semi/mesh.py XDMF NotImplementedError (Phase C)
- input schema additionalProperties: false enforcement (Phase D)
- semi/fem/sg_assembly.py dead-on-active-path (Phase E)

## Acceptance tests

(All four are in docs/IMPROVEMENT_GUIDE.md § M14.3.)

- [x] A1: scripts/run_benchmark.py mosfet_2d Pao-Sah verifier
- [x] A2: resistor benchmark from box.xdmf vs box.msh within 1e-12
- [x] A3: every benchmark JSON validates as v2.0.0; typo rejected
- [x] A4: SG primitives gone; coverage gate 95; CI green

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark mosfet_2d
- [ ] docker compose run --rm benchmark resistor_3d (msh and xdmf)
```

## Hand-off

When M14.3 lands and is reviewed, the next pickup is M16.1
(Caughey-Thomas field-dependent mobility). Paste
`docs/M16_1_STARTER_PROMPT.md` into a fresh Claude Code session on
branch `dev/m16.1-caughey-thomas`. Phase 2 of the post-M15 roadmap
refresh is a separate PR; it does not start until M14.3 is
merged.
