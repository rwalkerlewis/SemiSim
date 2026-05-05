# Physics study guide for kronos-semi

A pedagogical, first-principles walkthrough of every line of physics in
the kronos-semi engine (v0.16.0+, schema v2.0.0, M16.1 shipped). After
working through this guide, you should be able to open any file in
`semi/`, `kronos_server/`, or `benchmarks/` and explain the physical
meaning of every term, parameter, boundary condition, scaling factor,
time-integrator coefficient, and design choice.

## Who this guide is for

You should arrive with:

- Undergraduate electromagnetism (Maxwell's equations, vector calculus).
- Undergraduate quantum mechanics (Schrödinger equation, free-particle states).
- Working knowledge of PDEs, ODEs, linear algebra.
- Familiarity with Python.

You do **not** need:

- Prior semiconductor device-physics background (TCAD, COMSOL, etc.).
- Prior finite-element method experience.
- Prior numerical-solver experience beyond textbook Newton.

## How this guide relates to the existing docs

The repo already has substantial documentation under `docs/`:

- [`docs/PHYSICS.md`](../PHYSICS.md) — reference: the canonical
  governing equations, scaling conventions, boundary conditions.
- [`docs/PHYSICS_INTRO.md`](../PHYSICS_INTRO.md) — narrative on-ramp
  for programmers without TCAD background. **If you want a fast 30-min
  read, start there.**
- [`docs/theory/`](../theory/) — focused theory notes.
- [`docs/ARCHITECTURE.md`](../ARCHITECTURE.md), [`docs/WALKTHROUGH.md`](../WALKTHROUGH.md), [`docs/ROADMAP.md`](../ROADMAP.md), [`docs/IMPROVEMENT_GUIDE.md`](../IMPROVEMENT_GUIDE.md), [`docs/adr/`](../adr/), [`docs/schema/reference.md`](../schema/reference.md), [`docs/gpu.md`](../gpu.md).

This study guide does **not** duplicate any of that. It is the
**learning path** that turns those reference documents into something
a beginner can absorb in a defined order. Cross-references into
`docs/` are for going *deeper*, not for filling gaps. Each chapter
explicitly cites which existing docs to consult next.

## Reading order

The chapters are ordered so each chapter only depends on physics
introduced in earlier chapters.

| # | Chapter | Why now |
|---|---|---|
| 1 | [Classical EM and Poisson](01_classical_em_and_poisson.md) | The PDE that everything is built on. |
| 2 | [Solid state and band theory](02_solid_state_and_band_theory.md) | Where carriers come from. |
| 3 | [Carrier statistics](03_carrier_statistics.md) | Boltzmann, mass action, $V_t$, $n_i$. |
| 4 | [Doping and charge neutrality](04_doping_and_charge_neutrality.md) | The bulk equilibrium. |
| 5 | [Drift-diffusion transport](05_transport_drift_diffusion.md) | Currents and continuity. |
| 6 | [Generation–recombination](06_generation_recombination.md) | SRH, Auger preview, radiative preview. |
| 7 | [pn junction physics](07_pn_junction_physics.md) | The first actual device. |
| 8 | [Metal-semiconductor and ohmic contacts](08_metal_semiconductor_and_ohmic_contacts.md) | What a contact is. |
| 9 | [MOS capacitor physics](09_mos_capacitor_physics.md) | Band bending under a gate. |
| 10 | [MOSFET physics](10_mosfet_physics.md) | The device of the digital era. |
| 11 | [Quasi-Fermi and Slotboom](11_quasi_fermi_and_slotboom.md) | The variable change that makes Galerkin work. |
| 12 | [Nondimensionalization](12_nondimensionalization.md) | Why the engine doesn't blow up at $10^{30}$ condition. |
| 13 | [FEM weak forms](13_fem_weak_forms.md) | Strong → weak via test function. |
| 14 | [Multi-region and interfaces](14_multi_region_and_interfaces.md) | Si/SiO₂ flux continuity. |
| 15 | [Axisymmetric formulations](15_axisymmetric_formulations.md) | r-weighted weak forms. |
| 16 | [Nonlinear solvers and continuation](16_nonlinear_solvers_and_continuation.md) | Newton, line search, bias ramp. |
| 17 | [Transient time integration](17_transient_time_integration.md) | BDF1 / BDF2 with Slotboom. |
| 18 | [Small-signal AC analysis](18_small_signal_ac_analysis.md) | Frequency-domain linearization. |
| 19 | [Linear solvers and GPU backends](19_linear_solvers_and_gpu_backends.md) | MUMPS, AMGX, hypre BoomerAMG. |
| 20 | [Material parameter database](20_material_parameter_database.md) | Si, Ge, GaAs; SiO₂, HfO₂, Si₃N₄. |
| 21 | [Verification and benchmarks](21_verification_and_benchmarks.md) | MMS, conservation, the eight benchmarks. |

Appendices:

- [Appendix A — Constants and units](appendix_A_constants_and_units.md). Quick reference.
- [Appendix B — Full derivations](appendix_B_full_derivations.md). Long algebraic detours referenced from chapters.
- [Appendix C — Glossary and symbols](appendix_C_glossary_and_symbols.md).
- [Appendix D — Milestone physics map](appendix_D_milestone_physics_map.md). M1 retrospective + M16.x preview.
- [References](references.md). Full bibliography.

The [00_inventory.md](00_inventory.md) file is the concept-to-code-to-docs
map; consult it if you want to see where any specific concept lives in
the codebase.

## Reading paths

**For a fast 30-min read:** [`docs/PHYSICS_INTRO.md`](../PHYSICS_INTRO.md)
in the parent docs/. This study guide is for the deep dive.

**For the device physics path** (if you don't know semiconductors):
read 1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10 in order. By the end
you'll understand every shipped benchmark physically.

**For the numerics path** (if you know semiconductors but not FEM):
read 11 → 12 → 13 → 14 → 15 → 16 → 17 → 18 → 19. You can skip 1-10
and refer back as needed for specific physics.

**For the V&V path** (if you want to verify the engine):
read 21 first, then back-fill from chapters cited in its Code Map.

**For a specific milestone:** consult [Appendix D](appendix_D_milestone_physics_map.md)
for the chapter pointers, then read those chapters.

## How the code-map tables work

Every chapter has a **Code map** section: a table with three columns
`Concept`, `Equation in this chapter`, `Code location` (file:line).
The code locations cite specific lines; if a line moves, the chapter
text describes the function so you can grep for it. For
forward-referenced (planned) milestones, the code location points to
[`docs/IMPROVEMENT_GUIDE.md`](../IMPROVEMENT_GUIDE.md) instead of a
file path — fabricated paths would be wrong.

## How the cross-references work

- Within the guide: `[Ch. 11 §3](11_quasi_fermi_and_slotboom.md#slotboom-transformation)`.
- Into `docs/`: `[scaling theory](../theory/scaling.md)`.
- Into source: `[`semi/physics/poisson.py:60-66`](../../semi/physics/poisson.py)`.
- Into ADRs: `[ADR 0011](../adr/0011-ac-small-signal.md)`.

The existing `docs/` are authoritative on notation and code anchors;
this guide reconciles to them. If a derivation in this guide implies
a different design choice than the one locked in by an ADR, the ADR
wins; the chapter explains why or notes the disagreement honestly.

## Conventions

- All math in LaTeX.
- All quantities in SI by default. Where the JSON schema and the code
  use cm⁻³ for densities or cm²/(V·s) for mobilities, the chapter
  flags it explicitly; conversions live in
  [`semi/constants.py`](../../semi/constants.py).
- Equations are numbered per chapter, e.g. (7.5) is the fifth equation
  of Chapter 7. Cross-chapter references repeat the chapter number.
- Code locations are file paths from the repo root, with
  `:line` or `:line-line` ranges. Markdown links resolve relative to
  this directory.
- Existing-docs cross-references use relative paths
  (`[../theory/...]`) so the guide is portable.

## What to do if you find an error

This guide should be self-contained for someone with the listed
prerequisites. If you hit a step you can't follow:

1. Check the chapter's *Common pitfalls* and *Further reading*.
2. Check the Code Map for the function in question and read the
   docstring; engine docstrings are typically expansive.
3. Check the relevant existing-docs cross-reference; the chapter
   identifies it.
4. Open an issue. The guide is meant to be living documentation; gaps
   in the *GAPS* section of [00_inventory.md](00_inventory.md) are the
   roadmap for revisions.

## Length and depth

Most chapters are 2000-3000 words and contain:

1. Learning objectives.
2. Physical motivation (why this concept exists in device simulation).
3. First-principles derivation (no "it can be shown that").
4. Key results (numbered equations, unit-checked).
5. Worked numerical example (using benchmark values).
6. Code map.
7. Existing-docs cross-reference.
8. Common pitfalls (≥3).
9. Exercises (3-6, with solutions).
10. Further reading (specific section/page references).

Depth is the priority over brevity. If you want shallow, read PHYSICS_INTRO.md.

## Status (as of v0.17.0)

- Physics shipped: equilibrium Poisson, coupled Slotboom drift-diffusion,
  SRH recombination, ohmic / gate contacts, multi-region (Si/SiO₂),
  axisymmetric (cylindrical) 2D, BDF1/BDF2 transient, AC small-signal,
  Caughey-Thomas mobility (M16.1).
- Physics planned: Lombardi mobility (M16.2), Auger (M16.3),
  Fermi-Dirac (M16.4), Schottky contacts (M16.5), BBT/TAT (M16.6),
  time-varying $V(t)$ (M16.7), heterojunctions (M17), 3D MOSFET (M19),
  MPI (M19.1), HTTP hardening (M20).

This guide covers all shipped physics in full and forward-references
the planned physics with milestone tags so the reader can follow
upcoming work without re-learning context.
