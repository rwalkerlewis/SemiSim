# 0019. 3D MOSFET geometry, mesh quality gate, equilibrium smoke

- Status: Proposed
- Date: 2026-05-30
- Milestone: M19 precursor
- Cross-references: ADR 0006 (V&V strategy), ADR 0007 (pure-Python
  core), ADR 0015 (Schottky Robin BC V&V scope departure precedent),
  ADR 0017 (adaptive transient dt V&V scope departure precedent),
  ADR 0018 (SNES line-search type V&V scope departure precedent).

## Context

`PLAN.md` "Next task" identifies M19 (3D MOSFET capstone) as the
sole unblocked roadmap item after M18.1 closed the bias_sweep half
of the `ed6719b` CI carve-out. Full M19 as scoped in
`docs/IMPROVEMENT_GUIDE.md` and `docs/ROADMAP.md` carries three
deliverables that compound risk in a single PR:

1. A 3D unstructured planar n-MOSFET geometry, mesh, and material /
   contact wiring.
2. A V_GS sweep with Pao-Sah analytical-comparison verifier.
3. A GPU comparison against the M15 CPU-MUMPS path at the device's
   actual DOF count, with the V_GS sweep as the workload.

`docs/ROADMAP.md` notes M19.1 (MPI orchestration) as the natural
follow-up; the M19 wall-clock expectation under fixed-step bias
ramp at the device's DOF count makes single-rank execution
prohibitive for CI without first knowing whether the geometry,
mesh, and equilibrium pipeline even hold together at 3D.

This ADR records the decision to ship a precursor PR that scopes
down to the geometry, mesh quality gate, and an equilibrium-only
smoke benchmark, deferring (2) and (3) to M19 proper. The
precursor de-risks M19 by surfacing any 3D gmsh / XDMF / runner
issues under a tractable equilibrium load that runs in the CI
budget.

## Decisions

### 1. Source / drain region strategy: doping-only differentiation

The geometry defines two physical volumes: Si (tag 1) and SiO2
gate oxide (tag 4). The source and drain n+ regions are
differentiated by the Gaussian implant doping profile in the
JSON `doping` block, not by separate gmsh physical volume tags.

The 2D MOSFET precedent (`benchmarks/mosfet_2d/mosfet_2d.json`)
uses doping-only differentiation. The 3D analogue follows the
same convention because:

- The Slotboom DD residual is region-agnostic for n+ vs n
  junctions inside the same semiconductor material; only the
  effective doping `N_D(x) - N_A(x)` enters the Poisson coupling.
- A single Si region simplifies the gmsh `.geo` (one volume
  instead of three) and the schema `regions` array (one entry
  instead of three with identical material parameters).
- The mesh quality gate's per-region edge-length thresholds
  (Decision 2) are then keyed off geometric proximity to the
  oxide and to the implant centers rather than off region tags,
  which removes a class of failure mode where a cell at the
  source / channel boundary is assigned the wrong quality
  threshold.

Reserved tag assignments (`si_source = 2`, `si_drain = 3`) remain
unused in the `.geo` to preserve the option of switching to
separate-tag differentiation in M19 proper if the Pao-Sah verifier
needs per-region carrier averaging at the contact.

### 2. Mesh quality metrics and thresholds

A new pure-Python module `semi/mesh/quality.py` (lazy dolfinx
imports per ADR 0007) exposes `check_mesh_quality(mesh, meshtags,
*, thresholds)` returning a `MeshQualityReport` dataclass with:

- min and max element edge length per physical group
- tetrahedral skewness max (gmsh Delaunay metric)
- tetrahedral aspect ratio max
- total cell count
- a `passed: bool` summary

The default thresholds for the `mosfet_3d_eq` benchmark are:

| Region            | Max edge length |
|-------------------|-----------------|
| Gate oxide (4)    | <= 1 nm         |
| Channel Si        | <= 10 nm        |
| Bulk Si           | <= 100 nm       |

Plus mesh-wide gates: min edge >= 1 nm; tet skewness <= 0.85;
tet aspect ratio <= 20; cell count in [200k, 500k]. Channel Si
is defined by geometric proximity (within 50 nm below the oxide
and 100 nm laterally inside the channel x-window).

The gate is opt-in: schema additive field `mesh.quality_gate`
(boolean, default false). Existing benchmarks do not need to
opt in and are byte-identical at the no-quality-gate branch. An
optional `mesh.quality_thresholds` sub-object overrides the
defaults per config.

The thresholds are chosen so the 5 nm oxide carries at least 5
cells across (oxide max edge <= 1 nm), the inversion-layer
formation region in the channel resolves at <= 10 nm consistent
with the M16.2 Lombardi mobility's surface roughness lever, and
the bulk does not blow the cell count budget. The 200k-500k cell
band brackets the equilibrium solve's expected DOF count at ~3x
nodes (psi + phi_n + phi_p Slotboom primary unknowns) within
the 5-minute CI wall-clock budget on a single thread.

### 3. V&V scope: audit-only via the equilibrium smoke verifier

No MMS variant ships in this PR. The precursor introduces no new
physics kernel; the Poisson and DD residual builders are reused
verbatim from M14.3 onwards. Per ADR 0006 amended, runner-driver
and geometry PRs use audit-only V&V via a benchmark-specific
verifier. ADR 0015 (Schottky Robin BCs), ADR 0017 (adaptive
transient dt), and ADR 0018 (SNES line-search type) all shipped
under the same exemption.

The verifier in `semi/verification/mosfet_3d_eq.py` gates five
properties at V_GS = V_DS = V_BS = 0:

1. `n_max` in the source and drain implant volumes within a
   factor of 3 of the N_D peak (1e20 cm^-3); the factor accounts
   for Gaussian-profile averaging across implant cells.
2. `n_min` in the channel (cells within 50 nm below the oxide
   over the channel x-window) below 1e12 cm^-3.
3. `p_max` in the body bulk (cells deeper than 500 nm from the
   surface, away from depletion) within 10 % of N_A (1e17
   cm^-3).
4. Charge neutrality `|n - p - N_D + N_A|` per cell below 1e15
   cm^-3 absolute in the body bulk.
5. Built-in junction potential: psi drops below psi_body by at
   least 3 * V_t at the lateral edges of the source / drain n+
   regions.

These are geometric / shape-based gates; they do not compare to
a closed-form analytical solution. The full Pao-Sah comparison
against analytical I-V is the M19-proper deliverable and
requires a V_GS sweep.

The `mosfet_3d_eq` CI matrix entry ships without
`allow-failure`. The wall-clock budget is <= 5 minutes single-
threaded on the CI worker; if the equilibrium solve exceeds
that, the bulk-Si max edge is coarsened (raised from 100 nm
toward 200 nm) before shipping a slow CI gate.

## Validation

(Filled in at Phase F close; see the M19 precursor PR
description for the observed mesh cell count, equilibrium SNES
iteration count, and CI wall time at merge.)

## Consequences

**Positive:**
- M19 proper (V_GS sweep + Pao-Sah verifier + GPU comparison) is
  de-risked: the geometry, mesh, and ingest pipeline are
  validated under an equilibrium load before the bias-sweep
  controller is asked to converge through the inversion onset
  on a 3D device.
- The `semi/mesh/quality.py` module is reusable for any future
  3D unstructured-mesh benchmark (3D HEMT, 3D HBT) and for the
  M19 proper PR's bias-sweep mesh as well.
- The schema additive field `mesh.quality_gate` extends the
  M11 schema-companion surface; UI consumers see a documented
  knob for mesh-quality enforcement.

**Negative / risk:**
- The benchmark exercises the equilibrium runner only. The bias
  sweep, AC sweep, and transient runners are not exercised on
  3D unstructured meshes by this PR; surprises specific to those
  runners would land in M19 proper or follow-ups.
- The mesh quality gate's thresholds are calibrated for the
  precursor benchmark only. Other future 3D benchmarks may need
  bespoke threshold overrides; the `mesh.quality_thresholds`
  sub-object covers that case.

## Deferred

- **V_GS sweep + Pao-Sah analytical verifier.** M19 proper.
- **M15 GPU comparison on the 3D MOSFET DOF count.** M19 proper.
- **MPI parallel orchestration with `mpiexec -n {1,2,4}`.** M19.1.
- **Per-region carrier averaging at source / drain contacts.**
  If the Pao-Sah verifier in M19 proper benefits from separate
  `si_source = 2` and `si_drain = 3` region tags, the `.geo`
  promotes the reserved tag IDs into active physical volumes;
  the schema and verifier follow.

## Alternatives (rejected)

1. **Ship full M19 in one PR.** Rejected: bundles geometry risk,
   bias-sweep-convergence risk, Pao-Sah-window calibration risk,
   and GPU-comparison-tuning risk into a single review. The
   precursor split lets each risk land independently and lets
   M19.1 (MPI) trail M19 proper without blocking the geometry
   deliverable.
2. **Skip the mesh quality gate; rely on the runner's existing
   solver-failure signal.** Rejected: a too-coarse oxide mesh
   silently undersamples the inversion-layer field rather than
   producing a SNES failure. The quality gate catches mesh
   regressions at ingest time instead of leaving them to a
   downstream verifier-tolerance miss.
3. **Generate the `.geo` via a Python gmsh API call (gmsh-occ
   in-process) at benchmark runtime.** Rejected: introduces a
   runtime gmsh dependency, defeats the resistor_3d precedent
   of committing the `.msh` so CI does not need gmsh, and
   couples mesh regeneration to the Python runtime in a way
   that complicates ADR 0007 (pure-Python core; gmsh-occ
   imports compile a sizeable C++ extension).
4. **Drop the gate oxide as a separate volume and represent it
   via a Robin BC on the channel surface (lumped-cap shortcut).**
   Held in reserve as a stop-condition fallback if the OCC
   boolean intersection between Si and SiO2 misbehaves at the
   channel / oxide interface. Not the shipped default; a Robin
   shortcut would prevent M19 proper from exercising the
   multi-region Poisson coupling that mosfet_2d already covers
   in 2D.

## References

- ADR 0006 (V&V strategy): runner-driver / geometry PRs ship
  audit-only V&V via benchmark verifiers; precedent for the no-
  MMS-variant exemption here.
- ADR 0007 (pure-Python core): the new `semi/mesh/quality.py`
  module is pure-Python except for lazy dolfinx imports inside
  function bodies.
- ADR 0015 (Schottky Robin BC): audit-only V&V precedent.
- ADR 0017 (adaptive transient dt): runner-driver V&V precedent.
- ADR 0018 (SNES line-search type): solver-driver V&V precedent.
- `docs/IMPROVEMENT_GUIDE.md` M19 entry: full-M19 acceptance
  tests that this precursor defers.
- `docs/ROADMAP.md` M19 row: Pao-Sah 25 %, velsat 30 %, 5x GPU
  speedup at 500k DOFs (all deferred to M19 proper).
- `benchmarks/mosfet_2d/mosfet_2d.json`: doping-only S/D
  differentiation precedent.
- `benchmarks/resistor_3d/fixtures/box.geo`: 3D gmsh-ingest
  precedent; commits both `.geo` and `.msh`.
- `benchmarks/moscap_axisym_2d/moscap_axisym.geo`: hand-authored
  `.geo` style reference.
