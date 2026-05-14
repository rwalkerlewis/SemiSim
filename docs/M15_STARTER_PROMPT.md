# M15 starter prompt: GPU linear solver path

Paste-ready prompt for Claude in VS Code. Mirrors the convention of
`docs/M9_STARTER_PROMPT.md`. If the prompt is useful, commit it as
`docs/M15_STARTER_PROMPT.md` after the first phase lands.

---

You are working in `SemiSim` at v0.14.1 (or v0.14.2 once PR #65 lands).
Your assignment is **M15: GPU linear solver path**.

This prompt does **not** restate the M15 deliverable, the acceptance tests,
or the GPU engineering plan. Those already live in
`docs/IMPROVEMENT_GUIDE.md` §4 M15 and §5. Read the guide; this prompt only
tells you the order in which to execute and which parts of the guide are
themselves out of date and need refreshing.

## Required reading (do not skip; ~45 minutes)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.14.1 and that PR #65 is
   merged. If PR #65 is still open, **stop and wait**; nothing in M15
   should land while a docs-cleanup PR is in flight.
2. `docs/IMPROVEMENT_GUIDE.md` sections §1 (Honest current state),
   §2 (Design principles, especially P5 and P6), §4 M15, §5 (GPU section,
   expanded), §6 (UI integration checklist), §8 (Anti-goals).
3. `docs/ARCHITECTURE.md` for the five-layer rule. The new compute-probe
   module lives in Layer 3 (pure Python core, no dolfinx). The GPU
   wiring lives in Layer 4.
4. `docs/adr/0001`, `0002`, `0003`, `0007`, `0011`. Skim the rest. ADR
   0011's errata on PETSc-real builds is directly relevant: complex
   PETSc on GPU is out of scope.
5. `docs/PHYSICS.md` only if you need to confirm a scaling convention.

## Triage

The M15 specification in `docs/IMPROVEMENT_GUIDE.md` is the work plan.
Follow it. Three things the guide does not state explicitly that you need
to know up front:

- §1 ("Honest current state") was written post-M8 and says "as of v0.8.0".
  M9 through M14.2 have shipped since. Triage decisions made against that
  stale section will be wrong. **Phase A below refreshes it before any
  GPU work begins.**
- M15's only listed dependency in the guide is M9 (artifact writer for
  timing instrumentation). M9 shipped in v0.9.0. M15 unblocks now.
- `PLAN.md` "Next task" lists "physics validation Phase 2" ahead of M15.
  The owner has explicitly re-prioritized M15 ahead of validation Phase
  2. Update `PLAN.md` in Phase A to reflect that; do not delete the
  validation line, move it to a backlog subsection.

Five phases. One commit per phase. Do not bundle.

---

## Phase A: refresh the development guide so triage is honest

No code. Documentation only. The point of doing this first is so that
every subsequent decision in this milestone is made against a current
view of the project.

1. `docs/IMPROVEMENT_GUIDE.md` §1 "Honest current state":
   - Replace the "as of v0.8.0" banner with the current `pyproject.toml`
     version.
   - Refresh "What exists and works" to include M9 through M14.2.
     Authoritative source: `PLAN.md` "Completed work log". Do not invent
     bullets; lift facts from there.
   - Refresh "What does not exist" by removing items that have shipped:
     result artifact serialization, server surface, schema versioning,
     transient solver, AC small-signal, axisymmetric coordinates.
   - The remaining honest gaps after the refresh are: GPU (M15), physics
     completeness (M16), heterojunctions (M17), Cartesian-2D MOSCAP
     variant and gate-driven HF C-V (both tracked as M14.2.x in
     `PLAN.md`).
2. `docs/IMPROVEMENT_GUIDE.md` §4 (Milestones):
   - Mark M9, M10, M11, M12, M13, M13.1, M14, M14.1, M14.2 as **Done**
     with a single-line summary and a date. Move the full deliverable
     blocks for those into a new appendix `§10 Shipped milestone detail`.
     Acceptance tests stay; they remain useful as regression criteria.
   - Leave M15, M16, M17, M18 sections intact.
3. `docs/IMPROVEMENT_GUIDE.md` §6 "UI integration checklist":
   - Tick off M9, M10, M11, `/materials`, `/schema`. They are done.
   - Leave `/capabilities` open and note that M15 will close it (GPU
     availability is a capability).
4. `docs/IMPROVEMENT_GUIDE.md` §9 "Change log":
   - Append a dated entry: "Refresh of §1, §4, §6 to reflect v0.14.x
     reality post-M14.2; introduce §10 shipped-milestone appendix."
5. `PLAN.md` "Next task":
   - Replace the v0.14.2 / validation Phase 2 lines with "M15: GPU
     linear solver path".
   - Move the "Physics validation suite, Phase 2" item into a new
     "Backlog" subsection if one does not already exist. Do not delete.
6. `README.md` cleanup. Even though `PLAN.md` and the dev guide use
   M-numbering, the README still has "Day 1 / Day 2 / Day 3 / Day 4 /
   Day 5 / Day 6 / Week 2 / Week 3" language. Bring the README into
   agreement with the rest of the project:
   - Remove every "Day N" and "Week N" reference from `README.md`.
   - In §Status, replace "Day 1 deliverable is the equilibrium Poisson
     solve..." with a milestone-anchored sentence that links to
     `PLAN.md` and `docs/ROADMAP.md`. The status content itself is
     short; the long status detail belongs in `PLAN.md` and
     `docs/ROADMAP.md` as it already does.
   - In §Design notes, drop the "(planned for Day 2+)" parenthetical
     from the Slotboom subheading. Slotboom shipped in M2.
   - In the inline Roadmap section of the README, replace the
     "Day 1 ✓ ... Day 6 ... Week 2 ... Week 3+" list with a short
     paragraph pointing at `docs/ROADMAP.md` and `PLAN.md` for the
     authoritative roadmap. The README is not the roadmap.
   - Rename `tests/check_day1_math.py` to
     `tests/check_analytical_math.py`. Update the reference in
     `CONTRIBUTING.md` "Before committing", and any CI workflow that
     calls it (`grep -rn check_day1_math .github/`).

Run `pytest tests/ -v` and `ruff check semi/ tests/`. No code paths
changed in this phase, so 238+ tests should still pass.

**Commit message:** `docs: refresh IMPROVEMENT_GUIDE for v0.14.x; remove residual day-language; rename check_day1_math`

---

## Phase B: schema surface for backend selection, no behavior change

Bump the input schema and add the user-facing knobs. CPU-MUMPS remains
the default; today's runs are bit-identical.

1. Bump `schemas/input.v1.json` minor: 1.3.0 → 1.4.0. Confirm the
   major-version gate in `semi/schema.py` (set in M11) still accepts
   1.4.0; bump `SCHEMA_SUPPORTED_MINOR` from 3 to 4.
2. Add to `solver`:
   ```jsonc
   "backend": {
     "enum": ["cpu-mumps", "gpu-amgx", "gpu-hypre", "auto"],
     "default": "cpu-mumps",
     "description": "Linear-solver backend. cpu-mumps is the V&V reference. auto picks the best available device at runtime."
   },
   "compute": {
     "type": "object",
     "additionalProperties": false,
     "properties": {
       "device":         {"enum": ["auto", "cpu", "cuda", "hip"], "default": "cpu"},
       "precision":      {"enum": ["float64"], "default": "float64",
                          "description": "Single precision is reserved for future M16+ work; only double is allowed today."},
       "preconditioner": {"enum": ["auto", "amgx", "hypre-boomeramg", "gamg"], "default": "auto"},
       "linear_solver":  {"enum": ["auto", "gmres", "bcgs", "cg"], "default": "auto"}
     }
   }
   ```
3. Cross-field validation in `semi/schema.py` (mirror the
   `coordinate_system` cross-field pattern from M14.2):
   - `backend == "cpu-mumps"` requires `compute.device in {"auto", "cpu"}`.
   - `backend in {"gpu-amgx", "gpu-hypre"}` requires
     `compute.device in {"auto", "cuda", "hip"}`.
   - `backend == "auto"` resolves at solve time, not validation time.
4. Default-fill: an input with no `backend` and no `compute` block
   produces an exact byte-equivalent solve to today.
5. Pure-Python tests in `tests/test_compute_schema.py`:
   - Every existing benchmark JSON validates unchanged against 1.4.0.
   - Cross-field rejections fire correctly.
   - Default-fill round-trips.
   - Schema-version major-gate still accepts 1.4.0 and rejects 2.0.0.

**Commit message:** `feat(schema): solver.backend and solver.compute (schema 1.4.0); CPU defaults unchanged`

**Acceptance:** all 238+ existing tests still pass; new tests pass; no
FEM behavior change is observable on any benchmark.

---

## Phase C: runtime backend probe, still no GPU code path

Pure-Python module that asks PETSc what it can do. Layer 3 only. Does
not import dolfinx, cupy, pycuda, or any vendor SDK directly.

1. New module `semi/compute.py`:
   - `available_backends() -> list[str]` returns the subset of
     `{"cpu-mumps", "gpu-amgx", "gpu-hypre"}` supported by the linked
     PETSc build. Probe via `petsc4py.PETSc.Mat.getTypes()` and
     `PETSc.Options()`. Detect `aijcusparse` / `aijhipsparse` for matrix
     types and `cuda` / `hip` for vector types. Detect `amgx` PC type.
     Detect `hypre` PC plus the `HYPRE_USE_GPU` build flag if
     readable; if not readable, treat hypre as CPU and return only
     `cpu-mumps` plus possibly `gpu-amgx`. Be conservative.
   - `device_info() -> dict`:
     ```jsonc
     {
       "engine_version": "...",
       "petsc_version": "...",
       "petsc_complex": false,
       "petsc_int64":   false,
       "backends_available": ["cpu-mumps", ...],
       "device_count":  null | int,
       "device_name":   null | str
     }
     ```
   - `resolve_backend(requested: str, available: list[str]) -> str`:
     - `"cpu-mumps"` always returns `"cpu-mumps"`.
     - `"gpu-amgx"` raises `ConfigError` with a clear message if AMGX
       is not in `available`.
     - `"gpu-hypre"` likewise.
     - `"auto"` prefers `gpu-amgx` over `gpu-hypre` over `cpu-mumps`.
   - Honor `SEMISIM_BACKEND` env var as the override default when
     `solver.backend == "auto"`.
2. `semisim_server`: add `GET /capabilities` returning `device_info()`
   plus the engine version and the schema version. This closes the
   open `/capabilities` line item in `docs/IMPROVEMENT_GUIDE.md` §6.
3. Pure-Python tests in `tests/test_compute_module.py`:
   - Mock `petsc4py.PETSc` for the `resolve_backend` logic tests.
   - Real probe for `available_backends()` so it adapts to whatever
     PETSc is in CI (CPU-only by default).
   - Round-trip a `device_info()` payload through the manifest schema.
4. Server tests in `tests/test_semisim_server.py` for the new endpoint.

**Commit message:** `feat(compute): runtime backend probe; GET /capabilities endpoint`

**Acceptance:** all tests pass; on the CI's CPU-only PETSc build,
`GET /capabilities` returns `backends_available: ["cpu-mumps"]` and
no GPU info.

---

## Phase D: wire the GPU path into the SNES solver

This is the core of M15. **Re-read `docs/IMPROVEMENT_GUIDE.md` §4 M15
and §5 before writing code.** The "What can actually run on GPU" table
in §5 is the implementation contract; do not deviate.

1. Touch only `semi/solver.py` and the runner files
   `semi/runners/{equilibrium,bias_sweep,mos_cv,transient,ac_sweep,mos_cap_ac}.py`.
   Do not touch the physics kernels (`semi/physics/*.py`), the BC layer
   (`semi/bcs.py`), or the scaling layer. If you find yourself
   tempted to, stop and re-read invariant 4 (pure-Python core) and
   ADR 0007 (BC interface).
2. Solve runs on device; assembly stays on host. Inject PETSc options
   via the `petsc_options_prefix` mechanism that ADR 0003 mandates:
   - `gpu-amgx`:
     `-mat_type aijcusparse -vec_type cuda -ksp_type gmres -pc_type amgx`
   - `gpu-hypre`:
     `-mat_type aijcusparse -vec_type cuda -ksp_type bcgs -pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_relax_type Chebyshev`
   - HIP variants: `aijhipsparse` and `vec_type hip`. Gate behind
     `compute.device == "hip"`.
   - `cpu-mumps` is unchanged from today's `DEFAULT_PETSC_OPTIONS`.
3. The Krylov method and preconditioner are user-overridable via
   `solver.compute.linear_solver` and `solver.compute.preconditioner`
   from Phase B. Resolve those at runtime in the runner before
   constructing `NonlinearProblem`. If the user picks a combination
   that PETSc rejects (e.g., `cg` on a non-symmetric AC system), let
   PETSc raise; do not pre-validate every combination.
4. Log the resolved backend, device, preconditioner, linear solver,
   KSP iteration count, and linear-solve wall time at INFO. These
   propagate into the manifest in step 6.
5. Per invariant P5: GPU only changes how fast, never what. Do not
   alter `result.psi_phys`, `result.iv`, `result.C`, or any other
   public field shape or units. Continue to round-trip everything
   through `make_scaling_from_config` (locked since M1).
6. Manifest extension. Bump `schemas/manifest.v1.json` minor and add
   under `solver`:
   ```jsonc
   "backend_requested":   "auto",
   "backend_resolved":    "gpu-amgx",
   "device":              "cuda",
   "device_name":         "NVIDIA A100",
   "linear_solver":       "gmres",
   "preconditioner":      "amgx",
   "ksp_iters":           <int>,
   "linear_solve_wall_s": <float>
   ```
   Update `semi/io/artifact.py` and `semi/io/reader.py` accordingly.

**Acceptance tests** (these are the gates from
`docs/IMPROVEMENT_GUIDE.md` §4 M15; do not invent new ones, do not
weaken the tolerances):

- **A1 — equivalence.** On every existing benchmark, the GPU backend
  and the CPU-MUMPS backend produce solutions within `1e-8` relative
  L2. Test file: `tests/fem/test_gpu_backend.py`. Skip cleanly when no
  GPU is available; the skip message must mention which backend was
  missing.
- **A2 — speedup.** New benchmark `benchmarks/poisson_3d_gpu/`:
  3D Poisson on a uniformly doped block at ≥500k DOFs. The wall-clock
  for the **linear-solve portion alone** (not assembly, not I/O) is
  ≥5x faster on GPU than on CPU-MUMPS on the same host. Capture the
  timing in the manifest. Do not measure SNES outer-loop time.
- **A3 — no silent fallback.** With `backend: "auto"` on a CPU-only
  host, the run completes on CPU and the manifest records
  `backend_resolved: "cpu-mumps"`. With `backend: "gpu-amgx"` on a
  CPU-only host, the run **fails fast** at validation time with a
  clear `ConfigError` from `semi.compute.resolve_backend`. This is
  invariant P5 read strictly: no silent CPU fallback when a GPU
  backend was explicitly requested.

**CI:** existing CI covers the CPU path only. Add a guarded job
`gpu-nightly` in `.github/workflows/` that runs on a `gpu` self-hosted
runner only, gated by `if: ${{ vars.GPU_RUNNER_AVAILABLE == 'true' }}`.
Do not attempt to build PETSc-CUDA inside GitHub-hosted runners; the
build is approximately 30 minutes and will time out the free runner.

**Commit message:** `feat(solver): GPU linear-solver path via PETSc CUDA/HIP (M15)`

---

## Phase E: close out

1. `PLAN.md`:
   - Move M15 in §Roadmap from `Planned` to `Done` with a date.
   - Append an entry to §"Completed work log" following the M14.2
     entry's format: PR number, deliverables, schema version bump,
     acceptance-test results, no-regressions confirmation,
     `pyproject.toml` version (bump minor: 0.14.2 → 0.15.0).
   - Refresh §"Current state" to reflect the new version.
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M15 Done in §4 with a one-line summary; move the M15
     deliverable detail block into the §10 appendix you created in
     Phase A.
   - Tick off `/capabilities` in §6 (closed in Phase C).
   - Append a §9 changelog entry.
3. `README.md` "What works today" feature list: add a one-liner
   about the optional GPU backend, linking to a new `docs/gpu.md`
   you write in this phase. `docs/gpu.md` covers:
   - Install: PETSc-with-CUDA / PETSc-with-HIP, AMGX, hypre with
     `--with-cuda`. Conda recipe, Docker recipe.
   - Known limits: assembly stays on host (P5 / §5 anti-goal); complex
     PETSc on GPU is out of scope (ADR 0011 errata); single precision
     is not yet enabled.
   - Acceptance numbers from A2.
4. `CHANGELOG.md`: bump entry to v0.15.0 with the M15 line.

**Commit message:** `docs: close out M15 (PLAN, IMPROVEMENT_GUIDE, README, CHANGELOG, gpu.md)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes anywhere in prose touched (PLAN.md invariant 8).
      Use commas, parentheses, or colons.
- [ ] Pure-Python core remains dolfinx-free. `semi.compute` imports
      only `petsc4py`, `os`, and stdlib (PLAN.md invariant 4 / ADR
      0007 / Layer 3).
- [ ] CPU-MUMPS path is bit-identical to pre-M15 on every benchmark
      (P5).
- [ ] `make_scaling_from_config` still on every solve path (P /
      M1-locked).
- [ ] No PETSc / UFL types leak into `semisim_server` public API (P6).
- [ ] Schema bumped per minor when fields are added; major-gate logic
      still passes (P4 / M11).
- [ ] No new ADR introduced unless an invariant is broken; if so,
      stop and draft an ADR in `docs/adr/` first.
- [ ] No `pip install` of `cupy` / `pycuda` / `pyamgx` as required
      dependencies; GPU support is detected via PETSc's compiled-in
      capabilities only.

## Anti-goals (from `docs/IMPROVEMENT_GUIDE.md` §8)

- Do not port assembly to GPU. The 80/20 in M15 is the linear solve.
  dolfinx-cuda is experimental; revisit post-M16.
- Do not rewrite any module in Rust or C++.
- Do not add a Python-side GPU library (cupy, JAX, torch) anywhere
  in `semi/`.
- Do not silently fall back from GPU to CPU on a request mismatch
  (A3).
- Do not introduce `mat_type` or `vec_type` fields into the JSON
  schema. Those are PETSc internals; the user-facing surface is
  `solver.backend` plus `solver.compute.device`.
- Do not change anything in `semi/physics/*`, `semi/bcs.py`, or
  `semi/scaling.py` for this milestone.
- Do not edit prior entries in `PLAN.md` "Completed work log". It
  is append-only.

## Working method

1. Read the required reading list. Do not skim.
2. Run the existing test suite once before touching anything; record
   the count and the coverage so you can compare at the end.
3. Phase A (docs only, no code).
4. Phase B (schema, no behavior change).
5. Phase C (probe module, no behavior change).
6. Phase D (the actual GPU work).
7. Phase E (close out).
8. After each phase: `pytest tests/ -v && ruff check semi/ tests/`.
   Do not advance to the next phase with red tests.

## Hand-off

When all five phases land, the obvious next picks are M14.2.x
(Cartesian-2D MOSCAP variant, ~2 days), the deferred validation
suite Phase 2, or M16.1 (Caughey-Thomas mobility). Pick exactly one
per the project rule, do not auto-pick; surface the choice to the
owner.
