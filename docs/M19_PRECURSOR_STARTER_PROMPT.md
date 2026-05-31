## M19 precursor starter prompt: 3D MOSFET geometry, mesh, equilibrium smoke

Authoritative worker prompt. Read this in full, then read `PLAN.md`,
`docs/IMPROVEMENT_GUIDE.md`, and `docs/ARCHITECTURE.md` before touching
code. Do not improvise on scope; if a phase deliverable disagrees with
what you find in the tree, stop and surface it in a commit message
note before proceeding.

## Mission

Stand up a 3D planar n-MOSFET as a gmsh-sourced unstructured device,
gate the mesh quality, and demonstrate that the existing pipeline
(schema validate, mesh ingest, equilibrium Poisson solve under
Slotboom DD primary unknowns at V_GS = V_DS = 0) runs to completion
on it. Ship one new benchmark, `mosfet_3d_eq`, with an equilibrium-
only verifier. No bias sweep, no Pao-Sah verifier, no GPU benchmark,
no MPI. Those are M19 proper and M19.1.

This is the narrower of the two M19 options PLAN.md allows ("the
gmsh-sourced 3D MOSFET geometry and mesh-quality gate"). Ship as
v0.26.0, branch `dev/m19-precursor-mosfet-3d-geometry`.

## Why this is the right next chunk

Full M19 (V_GS sweep, Pao-Sah verifier, M15 GPU comparison) is too
big for one PR and expected to need M19.1 MPI to be tolerable at
runtime. Splitting yields three concrete deliverables: (1) the
geometry and mesh-quality work in this PR, (2) the full sweep and
Pao-Sah gate in M19, (3) MPI orchestration in M19.1. This PR de-
risks (2) and (3) by surfacing any 3D gmsh / XDMF / runner issues
under a tractable equilibrium-only load.

The Copilot-assigned M18.1 (bias-sweep SNES line-search
stabilization) is orthogonal: it edits `semi/runners/bias_sweep.py`,
the schema, and CI YAML; this PR edits `benchmarks/`, `semi/mesh.py`
(if needed), `semi/verification/`, and adds a new gmsh `.geo`
template. The two can land in either order.

## Constraints (from PLAN.md, do not violate)

1. JSON only. No new config formats.
2. Nondimensionalization mandatory.
3. Mesh coordinates stay in meters. The `.geo` writes in meters and
   the ingest reads in meters; no centimeter conversions anywhere
   in the mesh layer.
4. Pure-Python core (`schema`, `materials`, `scaling`, `doping`,
   `constants`, `continuation`, `diode_analytical`, `bcs`) must not
   import dolfinx. The new verifier is FEM-side (Layer 4) and lives
   under `semi/verification/` alongside the existing Pao-Sah helper.
5. dolfinx 0.10 API only. `NonlinearProblem` with
   `petsc_options_prefix`.
6. Slotboom variables for DD. ADR 0004 stays.
7. Physics-style variable names allowed.
8. No em dashes anywhere in code or prose.
9. Docker for dev. The acceptance runs under
   `docker compose run --rm benchmark mosfet_3d_eq` and the CI matrix.

## Repo orientation

Read in this order before opening any file for edit:

1. `PLAN.md` "Current state" and "Next task". Source of truth.
2. `docs/IMPROVEMENT_GUIDE.md` § 1 (gap list) and the M19 entry.
3. `docs/ROADMAP.md` for the M19 acceptance test wording.
4. `docs/ARCHITECTURE.md` for the five-layer rules.
5. `docs/adr/0007-pure-python-core.md` (import rules).
6. `benchmarks/mosfet_2d/` for the 2D MOSFET precedent: schema
   shape, doping block conventions, Gaussian implant parameters,
   the existing Pao-Sah verifier wiring. The 3D config follows
   the same shape with an added third dimension.
7. `benchmarks/resistor_3d/` for the 3D mesh-ingest precedent:
   gmsh `.geo` template, physical-group tag conventions, builtin
   vs gmsh mesh switching in the schema.
8. `benchmarks/moscap_axisym_2d/` for the gmsh `.geo` style
   currently shipping in the repo.
9. `semi/mesh.py` `_build_from_file` for the gmsh / XDMF dispatch
   that the new benchmark consumes.
10. `semi/doping.py` for the Gaussian implant profile shape in 3D
    (the 2D MOSFET uses the same code path; confirm 3D works
    end-to-end on the precursor geometry before assuming it does).

If any of these have drifted from what this prompt describes,
surface it in the Phase 0 commit message.

## Geometry specification

Standard planar n-MOSFET, channel along x, width along y, depth
along z. All dimensions in meters in the `.geo`.

- Total domain: 1.5 um (x) by 1.0 um (y) by 2.0 um (z).
- Channel length L_g = 250 nm centered at x = 0.75 um.
- Source: x in [0, 0.5 um], y in [0, 1 um], z in [1.9, 2.0] um.
- Drain: x in [1.0, 1.5 um], y in [0, 1 um], z in [1.9, 2.0] um.
- Body: rest of the Si volume below z = 2.0 um.
- Gate oxide: 5 nm SiO2 over the channel only
  (x in [0.5, 1.0] um, y in [0, 1 um], z in [2.0, 2.005] um).
- Gate poly contact: top surface of the oxide
  (z = 2.005 um plane over the channel x and y window).
- Body contact: bottom face z = 0.
- Source contact: top face of the source n+ region.
- Drain contact: top face of the drain n+ region.

Doping (Gaussian implants on a p-type body):

- Body: N_A = 1e17 cm^-3 (uniform).
- Source / drain peaks: N_D = 1e20 cm^-3, sigma_x = 30 nm,
  sigma_z = 50 nm, sigma_y large (treat as uniform across width,
  or sigma_y = 200 nm if the 2D precedent's implant profile shape
  requires a finite value).

Physical groups in the `.geo` (mirror the M14.3 XDMF naming):

- Volume tags: `si_body = 1`, `si_source = 2` (optional region tag,
  may also be done via doping only), `si_drain = 3` (same caveat),
  `sio2_gate = 4`.
- Surface tags: `source = 10`, `drain = 11`, `gate = 12`,
  `body = 13`.

Decide once whether source and drain are separate region tags or
single Si region with doping-only differentiation, document the
choice in ADR 0019 (Phase 0), and stay consistent. The 2D MOSFET
precedent uses doping-only differentiation; follow that unless the
3D case forces otherwise.

## Mesh quality targets

The mesh must satisfy at ingest time:

- Min element edge length >= 1 nm.
- Max element edge length in the gate oxide volume <= 1 nm
  (the 5 nm oxide needs at least 5 cells across).
- Max element edge length in the channel Si (defined as the volume
  within 50 nm below the oxide and 100 nm laterally inside the
  channel edges) <= 10 nm.
- Max element edge length in the bulk Si (everywhere else) <= 100
  nm.
- Tetrahedral skewness max <= 0.85 (gmsh default Delaunay-based
  metric).
- Tetrahedral aspect ratio max <= 20.
- Total cell count target 200k to 500k.

The quality metrics live in a new pure-Python module
`semi/mesh/quality.py` (since this is post-ingest checking, it can
read dolfinx mesh objects but must not import dolfinx at module
scope; lazy imports in functions only).

## Phases

Phase 0, A, B, C, D, E, F, one commit each, one PR. Follow the
M18 / M16.7 commit-layout convention.

### Phase 0: starter prompt and ADR draft

Ship this file at `docs/M19_PRECURSOR_STARTER_PROMPT.md` verbatim
on the branch `dev/m19-precursor-mosfet-3d-geometry`. Open ADR 0019
(`docs/adr/0019-3d-mosfet-geometry-and-mesh-quality.md`) as a draft
(Status: Proposed). The ADR documents (a) source / drain region
strategy (separate tags vs doping-only), (b) the mesh quality
metrics and their thresholds, (c) the V&V scope departure (no MMS
variant, audit-only V&V via the equilibrium smoke verifier; same
precedent as ADR 0015, ADR 0017, ADR 0018). Move to Accepted at
the close of Phase F.

Merge Phase 0 via a small PR ahead of the rest of the milestone
(the M16.5 / M16.6 / M16.7 / M18 pattern).

### Phase A: gmsh geometry template

Create `benchmarks/mosfet_3d_eq/mosfet_3d.geo` per the Geometry
specification above. Hand-author the `.geo` (do not generate via
a Python gmsh API call at this stage; the precedent in
`benchmarks/moscap_axisym_2d/` is a hand-authored `.geo`).

The `.geo` defines points, lines, surfaces, volumes, physical
groups, and mesh-size fields. Use gmsh's `Field` mechanism for the
graded mesh: a `Distance` field anchored on the gate oxide and the
source/drain junction surfaces, and a `Threshold` field mapping
distance to target edge length per the Mesh quality targets.

Generate the mesh under docker-fem via:

```
docker compose run --rm dev gmsh \
  -3 -format msh4 \
  benchmarks/mosfet_3d_eq/mosfet_3d.geo \
  -o benchmarks/mosfet_3d_eq/mosfet_3d.msh
```

Commit both the `.geo` and the generated `.msh` so CI does not need
gmsh in the runtime image (the resistor_3d precedent commits the
`.msh`; do the same here). Document the gmsh version used (gmsh
4.x) in the benchmark README.

Tests for Phase A: a pure-Python `tests/test_mosfet_3d_geo.py`
asserts the `.msh` file exists, parses via meshio, has the expected
physical groups, and the expected element count range.

### Phase B: mesh quality module and ingest validation

Create `semi/mesh/quality.py` with:

- `check_mesh_quality(mesh: dolfinx.mesh.Mesh, meshtags, *, thresholds: dict) -> MeshQualityReport`
- `MeshQualityReport` dataclass with min/max edge length per
  physical group, skewness max, aspect ratio max, cell count, and
  a `passed: bool` summary.

The thresholds dict has the per-physical-group expected ranges
from the Mesh quality targets section above. The Si body / channel
/ oxide partitioning is done by looking up the cell physical
group tags.

Integrate into `semi/run.py` so that when a config sets
`mesh.quality_gate: true` (new schema field, see Phase C), the
quality check runs after `_build_from_file` and raises a
`MeshQualityError` on failure.

Tests for Phase B:
`tests/fem/test_mesh_quality.py` exercises the module on a known-
good mesh (the generated mosfet_3d.msh) and a deliberately-broken
mesh (a copy of the .geo with a too-large element size in the
oxide); the bad mesh must fail the gate.

### Phase C: schema additions

Add `mesh.quality_gate` boolean (default false) and
`mesh.quality_thresholds` sub-object (optional override of the
defaults in `semi/mesh/quality.py`) to `schemas/input.v2.json` and
`semi/schema.py`. Schema additive minor bump v2.9.0 -> v2.10.0;
`SCHEMA_SUPPORTED_MINOR` advances from 9 to 10.

If the Copilot M18.1 PR landed before this one and already bumped
the schema to v2.10.0, coordinate the bump (this PR ships v2.11.0;
both PRs additively contribute their own fields; the merge order
determines which one publishes v2.10.0 first). The bump number is
not load-bearing; the additive-only invariant is.

Default-fill: `quality_gate` defaults to false (existing benchmarks
do not need to be edited to opt in). The new `mosfet_3d_eq`
benchmark sets it to true.

Tests for Phase C:
`tests/test_mesh_quality_schema.py` exercises default-fill,
override threshold parsing, and every existing benchmark JSON
validates as v2.10.0 (or v2.11.0 if the Copilot bump landed first)
unchanged.

### Phase D: equilibrium benchmark config

Create `benchmarks/mosfet_3d_eq/mosfet_3d_eq.json`:

- `schema_version: "2.10.0"` (or v2.11.0 per Phase C note).
- `dimension: 3`.
- `mesh.source: "file"`, `mesh.path: "mosfet_3d.msh"`.
- `mesh.quality_gate: true`.
- `regions`: one Si region (tag 1, or three if separate source /
  drain regions per the Phase 0 ADR decision), one SiO2 region
  (tag 4).
- `doping`: body Gaussian implants per the Geometry specification.
- `contacts`: four ohmic / gate contacts at V = 0 each
  (V_GS = V_DS = V_BS = 0). The gate is type `"gate"` with
  `work_function_eV = 4.05` (n+ poly approximation) and the oxide
  region linked via `regions[].role: "insulator"`.
- `physics.statistics: "boltzmann"` (FD is M16.4; not needed for
  equilibrium at moderate doping; the source / drain peaks at 1e20
  cm^-3 are heavy but the equilibrium n in those regions is not
  the gate of this verifier; if equilibrium V_bi at the junctions
  needs FD, switch to `"fermi_dirac"` and document in the ADR).
- `physics.recombination`: SRH only, defaults.
- `solver.type: "equilibrium"`. Not `"bias_sweep"`; the equilibrium
  runner solves Poisson once and recovers n / p from Slotboom
  primary unknowns. This is the simplest valid runner for the
  precursor.

Tests for Phase D:
A `tests/test_mosfet_3d_eq_config.py` asserts the JSON validates
against the active schema minor.

### Phase E: equilibrium verifier

Create `semi/verification/mosfet_3d_eq.py` with
`verify_mosfet_3d_eq(result, cfg) -> VerificationReport`. The
verifier asserts:

1. n_max in the source and drain region volumes (cells whose centroid
   sits inside the implant boxes) is within a factor of 3 of N_D
   peak (1e20 cm^-3). The factor accounts for the Gaussian profile
   averaging across implant cells.
2. n_min in the channel (cells within 50 nm below the oxide,
   centered on the channel x-window) is below 1e12 cm^-3 at V_GS = 0.
3. p_max in the body bulk (cells deeper than 500 nm from the
   surface, away from the depletion regions) is within 10 % of
   N_A (1e17 cm^-3).
4. Charge neutrality holds in the body bulk: |n - p - N_D + N_A|
   per cell is below 1e15 cm^-3 in absolute terms (i.e. the
   numerical error in the equilibrium solve is small relative to
   the dopant density).
5. The depletion regions under the source / drain are present and
   geometrically reasonable: psi drops below psi_body by at least
   3 * V_t at the lateral edges of the source / drain n+ regions,
   indicating built-in junction potential is set up.

Register the verifier in `scripts/run_benchmark.py` under the
`verify_mosfet_3d_eq` name. The CI matrix gains a `mosfet_3d_eq`
entry without `allow-failure`.

Wall-clock budget for the CI run: target <= 5 minutes single-
threaded on the CI worker. If the equilibrium solve exceeds this,
coarsen the mesh in Phase A (raise the bulk-Si max edge to 200 nm)
rather than ship a slow CI gate. Document the trade in the ADR.

Tests for Phase E:
`tests/fem/test_mosfet_3d_eq_verifier.py` runs the full benchmark
under docker-fem and asserts each numerical gate holds. This is
also the integration test the CI matrix runs.

### Phase F: closeout and stale-doc sweep

The repo's README has been stale since v0.16.0 (the M14.4 README
rewrite removed milestone tags but newer milestones did not update
the Status section; the Out-of-scope list still mentions M16.x and
M17 even though both shipped). Closeout includes a doc sync from
v0.16.0 through v0.26.0, scoped to keep the diff focused on
factual accuracy rather than restructuring.

Update, in this order:

1. `PLAN.md`: move the "M19 precursor" line from "Next task" into
   Completed work log (append, never delete). Update "Current
   state" to v0.26.0 and to reflect M19 precursor shipped. Set the
   new "Next task" pointer to M19 proper (full V_GS sweep, Pao-Sah
   verifier, M15 GPU comparison) or, if M18.1 has not landed by
   the time this PR is ready, to that as the alternative. State
   both options; let the maintainer pick.
2. `README.md` § Status: bump the version line to v0.25.0 (the
   transient half) or v0.26.0 (this PR), refresh the capability
   bullets to include M16.1 through M16.7, M17, M18, and the M19
   precursor. Remove M16.1 through M16.7 and M17 from the
   "Out of scope today" list. Add the M19 precursor as in scope.
   Cite ROADMAP.md as the authoritative capability matrix.
3. `docs/IMPROVEMENT_GUIDE.md` § 1: prune the gap list to reflect
   shipped milestones. § 9 changelog under `[0.26.0]`.
4. `docs/ROADMAP.md`: add an M19 precursor row in the capability
   matrix, marked Done.
5. `CHANGELOG.md`: add a `[0.26.0]` entry under the existing
   Keep-a-Changelog structure. List the schema bump, the new
   `mesh.quality_gate` and `mesh.quality_thresholds` fields, the
   new `mosfet_3d_eq` benchmark, the new `semi/mesh/quality.py`
   module, the README factual-accuracy sweep.
6. `docs/PHYSICS_INTRO.md` §§ 6 and 7 if they reference the
   pre-M19 capability state.
7. `pyproject.toml` and `semi/__init__.py`: bump 0.25.0 -> 0.26.0.
8. Move ADR 0019 from Status: Proposed to Status: Accepted.

The doc sync is factual-accuracy only. Do not restructure, do not
rewrite sections for clarity, do not move content between files.
Touch the smallest set of lines needed to make the docs match the
state of `main` at this PR's merge.

## Acceptance tests (gate the PR)

1. `mosfet_3d_eq` benchmark runs to completion in CI without
   `allow-failure`. Wall clock <= 5 minutes single-threaded on the
   CI worker.
2. All five verifier gates in Phase E pass.
3. Mesh quality module reports all metrics within thresholds on
   the generated `mosfet_3d.msh`.
4. Existing-benchmark byte-identity. The schema bump is additive;
   every existing benchmark JSON validates as v2.10.0 (or v2.11.0
   if the Copilot bump landed first) and produces bit-identical
   results to v0.25.0. Canonical guards: `pn_1d_bias` J(V=0.6 V)
   = 1.635e+03 A/m^2, plus every anchor in PLAN.md "Current state".
5. Coverage gate at 95 on `semi/` holds. The new
   `semi/mesh/quality.py` and `semi/verification/mosfet_3d_eq.py`
   modules are covered by the new tests.
6. ADR 0019 is Accepted.
7. README and PLAN.md "Current state" are consistent with each
   other and with the actual main-branch capability set.

## Anti-goals

- **Do not add a V_GS or V_DS sweep.** Equilibrium only. Full
  sweep is M19 proper.
- **Do not add a Pao-Sah-style analytical verifier.** That gate
  needs bias-sweep data. Equilibrium verifier in this PR is
  geometric / shape-based only.
- **Do not run the M15 GPU backend in this benchmark.** Default
  CPU-MUMPS backend. GPU comparison is a separate benchmark
  setup in M19 proper.
- **Do not add MPI orchestration.** M19.1.
- **Do not add a new mesh format.** gmsh `.msh` (msh4 binary) is
  the locked input format for unstructured 3D meshes. XDMF is
  also supported per M14.3; the precursor ships gmsh because
  it has the most direct quality-gate metrics.
- **Do not refactor `semi/mesh.py` for clarity.** Targeted edits
  only.
- **Do not edit `semi/runners/bias_sweep.py` or `transient.py`.**
  This PR is equilibrium-only. The Copilot-assigned M18.1 owns
  bias_sweep edits.
- **Do not introduce a new dolfinx feature.** Use only the 0.10
  API surface already used elsewhere in the codebase.
- **Do not ship a notebook.** That's a follow-up after the
  benchmark runs reliably.
- **Do not ship an MMS variant.** No new physics kernel ships in
  this PR; per ADR 0006 amended, runner-driver / geometry PRs
  use audit-only V&V via the equilibrium verifier.

## Stop conditions (commit a note, do not push past)

- gmsh refuses to mesh the geometry (boolean intersection of the
  oxide and Si fails, physical group tagging breaks). Stop, write
  up the finding in the ADR draft, and consider a simpler
  geometry: drop the gate oxide as a separate volume and represent
  it as a thin lumped capacitance via a Robin BC on the channel
  surface. That is a deviation worth flagging before shipping.
- The equilibrium solve does not converge under default SNES
  settings. The 2D MOSFET equilibrium solve does converge today
  per the M16.2 acceptance test; the 3D analogue should also.
  If it does not, surface in the ADR; do not loosen SNES
  tolerances (ADR 0008 stays).
- Mesh cell count blows past 1M with the targets in this prompt.
  Coarsen the bulk-Si max edge first; if still over budget,
  shrink the device footprint (e.g. drop width to 500 nm). Do
  not ship a CI gate that takes longer than the 5-minute budget.
- The equilibrium verifier gate that needs the most attention is
  charge neutrality (gate 4). If the body bulk shows |n - p - N_D
  + N_A| above 1e15 cm^-3 at SNES convergence, the SNES tolerance
  may be too loose for this device size. Surface in the ADR,
  document the observed residual, and consider tightening on this
  benchmark only via a per-config solver override rather than a
  global tolerance change.

## File pointers (orientation, not blind-edit)

- `benchmarks/mosfet_2d/mosfet_2d.json`: 2D MOSFET schema shape,
  the doping block conventions, contacts list, physics block.
- `benchmarks/mosfet_2d/mosfet_2d.geo`: hand-authored 2D `.geo`
  (if shipped; if 2D MOSFET uses builtin mesh, fall back to the
  moscap_axisym_2d .geo for hand-authored style).
- `benchmarks/moscap_axisym_2d/moscap.geo`: confirmed hand-
  authored .geo shipping in the repo. Use as the style reference.
- `benchmarks/resistor_3d/`: 3D gmsh ingest precedent.
- `semi/mesh.py`: `_build_from_file` dispatch for gmsh / XDMF.
- `semi/doping.py`: Gaussian implant profile builder. Confirm
  the 3D code path on this profile exists; if it does not, add
  it (lazy import of dolfinx, per the pure-Python boundary).
- `semi/runners/equilibrium.py`: the runner this benchmark uses.
  Do not edit.
- `scripts/run_benchmark.py`: the CLI entry point that registers
  verifiers. Add `verify_mosfet_3d_eq` to the registry.
- `.github/workflows/ci.yml`: CI matrix to extend with the new
  benchmark entry.

## Branch and PR layout

- Branch: `dev/m19-precursor-mosfet-3d-geometry`.
- Commits: seven (Phase 0, A, B, C, D, E, F). One phase per commit.
- PR title: `M19 precursor: 3D MOSFET geometry, mesh quality gate,
  equilibrium smoke`.
- PR description: copy this prompt's "Mission" and "Acceptance
  tests" sections; add a one-line summary of mesh cell count, CI
  wall time, and any deviations from this prompt (sized to fit in
  the GitHub PR description without truncation).

End of prompt.
