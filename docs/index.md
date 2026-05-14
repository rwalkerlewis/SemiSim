# Documentation

Top-level table of contents for the SemiSim documentation tree.
For the project overview, install instructions, and quick start, see
the repository [README](../README.md).

## Project planning and history

- [PLAN.md](../PLAN.md) — current state, single in-flight task, invariants. Single source of truth for what to work on.
- [CHANGELOG.md](../CHANGELOG.md) — Keep-a-Changelog formatted release notes.
- [CONTRIBUTING.md](../CONTRIBUTING.md) — how to set up, what to read first, conventions.
- [LICENSE](../LICENSE) — MIT.
- [ROADMAP.md](ROADMAP.md) — per-milestone delivery history (M1 through current).
- [IMPROVEMENT_GUIDE.md](IMPROVEMENT_GUIDE.md) — milestone roadmap with explicit acceptance tests.

## Architecture and physics reference

- [ARCHITECTURE.md](ARCHITECTURE.md) — five-layer component design and per-layer import rules.
- [WALKTHROUGH.md](WALKTHROUGH.md) — end-to-end trace of `semi.run.run(cfg)` with file:line anchors.
- [PHYSICS_INTRO.md](PHYSICS_INTRO.md) — physics tutorial for programmers without a TCAD background.
- [PHYSICS.md](PHYSICS.md) — governing equations, scaling, BCs (reference).
- [PHYSICS_AUDIT.md](PHYSICS_AUDIT.md) — phase 1 cross-runner consistency audit results.
- [adr/](adr/) — locked architecture decisions. Open a new ADR before changing any invariant.

## Theory notes (extracted from the README during the post-merge cleanup)

- [theory/scaling.md](theory/scaling.md) — why nondimensional scaling, how it's wired through the engine.
- [theory/slotboom.md](theory/slotboom.md) — why Slotboom variables for drift-diffusion.
- [theory/dolfinx_choice.md](theory/dolfinx_choice.md) — why dolfinx 0.10 and `NonlinearProblem`.
- [theory/axisymmetric.md](theory/axisymmetric.md) — derivation of the r-weighted weak forms, axis BC, domain truncation.
- [theory/moscap_cv.md](theory/moscap_cv.md) — LF and HF C–V derivation; Hu Ch. 5 cross-references; honest note on the depletion-clamp HF stand-in.

## Schema and JSON contract

- [schema/reference.md](schema/reference.md) — full JSON input schema reference (extends the inline README example).
- [`schemas/input.v1.json`](../schemas/input.v1.json) — the on-disk Draft-07 schema.
- [`schemas/manifest.v1.json`](../schemas/manifest.v1.json) — the run-artifact manifest schema.

## Benchmarks

- [benchmarks/pn_junction_1d.md](benchmarks/pn_junction_1d.md) — equilibrium 1D pn junction.
- [benchmarks/moscap_axisym_2d.md](benchmarks/moscap_axisym_2d.md) — axisymmetric 2D MOSCAP, LF/HF C–V.
- [`benchmarks/`](../benchmarks/) — full benchmark directory (one subdir per device).

## Per-device derivations

- [mos_derivation.md](mos_derivation.md) — multi-region MOS Poisson + C-V derivation.
- [resistor_derivation.md](resistor_derivation.md) — 3D doped resistor.
- [mms_dd_derivation.md](mms_dd_derivation.md) — MMS for coupled drift-diffusion.

## Tasks and prompts

- [tasks/](tasks/) — historical task prompts kept as artifacts of the two-agent workflow described in [`CLAUDE.md`](../CLAUDE.md).

## Notebooks

- [`notebooks/01_pn_junction_1d.ipynb`](../notebooks/01_pn_junction_1d.ipynb) — equilibrium Poisson, 1D pn junction vs. depletion approximation.
- [`notebooks/02_pn_junction_bias.ipynb`](../notebooks/02_pn_junction_bias.ipynb) — forward Shockley + reverse SNS bias sweep.
- [`notebooks/03_mos_cv.ipynb`](../notebooks/03_mos_cv.ipynb) — 2D MOS capacitor C–V sweep.
- [`notebooks/04_resistor_3d.ipynb`](../notebooks/04_resistor_3d.ipynb) — 3D doped bar, builtin and gmsh meshes.
- [`notebooks/05_moscap_axisym_cv.ipynb`](../notebooks/05_moscap_axisym_cv.ipynb) — axisymmetric 2D MOSCAP, LF/HF C–V (reproduces Hu Fig. 5-18).

## Issues and follow-ups

Open follow-ups are tracked on GitHub:
[rwalkerlewis/SemiSim/issues](https://github.com/rwalkerlewis/SemiSim/issues).
