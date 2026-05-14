# §8 — SemiSim vs. COMSOL Semiconductor Module

**Suggested slide title:** "How Does SemiSim Compare to COMSOL's Semiconductor Module?"
**Target time:** 8–12 minutes (stand-alone section; trim to 4–5 min if inserting into a shorter talk)

---

## Slide 8.0 — Why This Comparison Matters

**[optional intro]**

COMSOL Multiphysics with the Semiconductor Module is the most widely
cited GUI-driven device simulator in academic and industrial settings. It
is expensive, mature, and feature-rich. Understanding where SemiSim
agrees with it, where it intentionally differs, and where it still falls
short is the most honest framing you can give an audience that has used
COMSOL.

This section is organized into seven dimensions: physics models, solver
strategy, geometry and meshing, input/output contract, extensibility,
access and cost, and verification practice. A summary table closes the
section.

---

## Slide 8.1 — Physics Models

**[0:00 relative to this section]**

Both SemiSim and COMSOL Semiconductor share the same underlying
classical physics: Poisson's equation coupled with drift-diffusion
continuity equations for electrons and holes. They agree on the
governing PDE set and on the form of the Shockley-Read-Hall
recombination kernel.

Where they diverge is in the breadth of the **physics catalog**:

| Physics feature | COMSOL Semiconductor | SemiSim (v0.16) |
|---|---|---|
| Drift-diffusion (Poisson + 2 continuity) | ✓ | ✓ |
| Boltzmann carrier statistics | ✓ | ✓ |
| Fermi-Dirac carrier statistics | ✓ | planned M16.4 |
| SRH recombination (mid-gap trap) | ✓ | ✓ |
| Auger recombination | ✓ | planned M16.3 |
| Radiative recombination | ✓ | planned |
| Caughey-Thomas field-dependent mobility | ✓ | planned M16.1 |
| Lombardi surface scattering mobility | ✓ | planned M16.2 |
| Constant mobility (Sze defaults) | ✓ | ✓ |
| Ohmic contacts (Dirichlet) | ✓ | ✓ |
| Schottky contacts | ✓ | planned M16.5 |
| Gate-oxide / MOS interface | ✓ | ✓ |
| Transient (time-domain) solves | ✓ | ✓ BDF1/BDF2 |
| AC small-signal / impedance sweep | ✓ | ✓ (linearised) |
| Quantum correction (density-gradient) | ✓ | not planned |
| Lattice heat equation (electrothermal) | ✓ | not planned |
| Optical generation (ray tracing) | ✓ add-on | not planned |
| Band-to-band tunneling | ✓ | planned M16.6 |
| Trap-assisted tunneling (TAT) | ✓ | planned M16.6 |
| Heterostructures / position-dep. band gap | ✓ | planned M17 |

The honest summary: COMSOL has a wider physics catalog. SemiSim
covers the classical drift-diffusion core that describes most
silicon devices operating below 10¹⁹ cm⁻³ doping and at moderate
fields. The M16 roadmap closes the most important gaps; quantum effects
and electrothermal are out of scope for the current version.

**Key points**
- Both solve the same classical drift-diffusion PDE set.
- COMSOL's catalog includes advanced recombination, quantum corrections,
  and electrothermal; these are planned or out-of-scope for SemiSim.
- The classical silicon core (pn junction, MOSFET, MOSCAP) is fully
  covered and verified.

---

## Slide 8.2 — Solver Strategy and Numerical Methods

**[2:00]**

This is where the approaches diverge most significantly, and where
SemiSim makes a different set of tradeoffs.

**COMSOL Semiconductor Module solver strategy:**
- Proprietary PARDISO direct solver (default) with algebraic
  multigrid preconditioned iterative options.
- COMSOL uses a **generalized Scharfetter-Gummel (SG) discretization**
  for the carrier continuity equations in its semiconductor module:
  a finite-volume scheme on element edges that guarantees positive
  carrier densities and is exact for 1D exponential profiles.
- The MUMPS and PARDISO direct solvers are used for small/medium
  problems; GMRES + ILU for larger ones.
- Newton damping is handled internally by COMSOL's nonlinear solver;
  users can tune under-relaxation and continuation parameters through
  the GUI.
- Convergence monitoring is exposed through the GUI but not
  programmatically inspectable at each Newton step.

**SemiSim solver strategy:**
- PETSc SNES Newton with MUMPS LU (CPU) or AMGX/hypre BoomerAMG
  (GPU), all configurable through the JSON `solver` block.
- **Slotboom (quasi-Fermi) transformation** instead of
  Scharfetter-Gummel. The SG scheme is natural on FV meshes but
  requires special quadrature on element edges in FEM, does not
  generalize to 3D unstructured meshes in a mesh-independent way, and
  degrades to first-order convergence in 2D. The Slotboom
  transformation re-expresses the continuity equations as pure
  gradients of quasi-Fermi potentials, restoring standard Galerkin
  stability and second-order L² convergence (confirmed by MMS).
- **Adaptive continuation** with halving-and-growth logic exposed in
  JSON as `solver.continuation.*`. Every step is logged; convergence
  histories are artifacts.
- Full transparency: every Newton residual norm, every continuation
  step, every linear solve time is written to the artifact tree and
  visible in the logs.

| Solver dimension | COMSOL | SemiSim |
|---|---|---|
| Primary unknowns | ψ, n, p (or Boltzmann form) | ψ, Φ_n, Φ_p (Slotboom) |
| DD discretization | Scharfetter-Gummel (FV-on-FEM) | Galerkin FEM on pure-gradient form |
| Convergence order (L², DD) | ~1st (SG) | 2nd (MMS confirmed) |
| Newton solver | COMSOL built-in | PETSc SNES |
| Linear solver | PARDISO / MUMPS | MUMPS / AMGX / hypre |
| GPU linear-solver path | no (CPU only in SM module) | yes (PETSc-CUDA/HIP) |
| Continuation control | GUI, some JSON export | JSON `solver.continuation` |
| Convergence artifacts | GUI plot, no export | JSON + CSV artifacts |

**Key points**
- SG is optimal for 1D FV; Slotboom is optimal for FEM in 3D.
- MMS confirms 2nd-order L² convergence for SemiSim DD; this is
  higher than typical SG in 2D/3D.
- COMSOL's nonlinear solver is opaque; SemiSim exposes every
  step for inspection and scripting.
- COMSOL has no GPU linear-solver path in the Semiconductor Module;
  SemiSim does (M15).

---

## Slide 8.3 — Geometry and Meshing

**[4:00]**

**COMSOL geometry workflow:**
- Full parametric CAD kernel (built on the COMSOL geometry engine and
  optionally LiveLink for SolidWorks/CATIA).
- Meshing is driven by the GUI with automatic refinement near
  boundaries, junctions, and high-field regions.
- Physical regions (semiconductor, oxide, metal) are defined by
  geometric entities in the COMSOL model tree.
- Mesh export is possible (NASTRAN, COMSOL MPH format) but round-tripping
  to external tools is non-trivial.
- Typical geometry + mesh + physics setup: 30–60 minutes GUI work for
  a new device.

**SemiSim geometry workflow:**
- **Built-in generators** for axis-aligned 1D/2D/3D boxes (adequate
  for benchmarks and most first-principles studies).
- **gmsh .msh import**: for arbitrary geometry, create the mesh in
  gmsh (open source) or any tool that exports .msh; physical groups
  propagate verbatim as cell and facet tags.
- **XDMF import**: for meshes generated by FEniCSx or dolfinx
  directly.
- Geometry is specified entirely in the JSON file (box dimensions,
  region bounds, doping bounds) or referenced as a path to a .msh
  file. No GUI needed to set up a new device.
- Physical groups in gmsh are mapped to material regions and contact
  facets by integer tag in the JSON `regions` block.

| Geometry dimension | COMSOL | SemiSim |
|---|---|---|
| CAD kernel | Yes (built-in + LiveLink) | No (use gmsh externally) |
| GUI-driven setup | Yes | No (JSON only) |
| Auto mesh refinement at junctions | Yes | No (manual in gmsh) |
| Unstructured 2D/3D meshes | Yes | Yes (via gmsh or XDMF) |
| Physical group / material tag propagation | Yes | Yes (gmsh integer tags) |
| Mesh export for external tools | Limited | gmsh .msh, XDMF |
| Setup time for a new device | 30–60 min (GUI) | 10–20 min (JSON + gmsh) |

**Key points**
- COMSOL's CAD kernel is a significant advantage for complex geometry.
- For rectangular/cylindrical devices (the vast majority of benchmarks),
  SemiSim's built-in generators require zero external tools.
- For non-rectangular geometry, gmsh + .msh import is the recommended
  workflow; gmsh is free and its Python API is scriptable.

---

## Slide 8.4 — Input/Output Contract and Reproducibility

**[6:00]**

This is the dimension where SemiSim makes the strongest claim of
superiority over COMSOL.

**COMSOL I/O:**
- A simulation lives in a `.mph` file — a proprietary binary format.
  Opening it requires a COMSOL license and a compatible version.
- Batch runs can be scripted via MATLAB LiveLink or COMSOL's Java API,
  but both require the full COMSOL installation.
- Results are exportable to CSV, text, or MATLAB `.mat` via the GUI
  or the COMSOL Java API. This requires deliberate action; it is not
  automatic.
- Reproducibility depends on the COMSOL version and any GUI settings
  not captured in export scripts. A `.mph` file from COMSOL 6.0 may
  not open identically in COMSOL 6.2.

**SemiSim I/O:**
- A simulation is **fully determined** by a JSON file and a package
  version. No binary file, no GUI state.
- The JSON schema is versioned (`v2.0.0`), validated against JSON
  Schema Draft-07, and enforced with `additionalProperties: false`.
  A typo in a field name is a validation error before any FEM code
  runs.
- Results land in `runs/<run_id>/` as a schema-validated
  `manifest.json`, XDMF field files, IV CSVs, and convergence logs.
  Any tool that can read JSON and XDMF can consume them.
- The HTTP API (`POST /solve`, `GET /runs/{id}`) makes SemiSim
  callable from any language, any UI, any AI assistant — without a
  Python runtime on the client side.

| I/O dimension | COMSOL | SemiSim |
|---|---|---|
| Simulation file format | `.mph` (proprietary binary) | JSON (open, text) |
| Version-compatible replay | Fragile across COMSOL versions | JSON + package version = full spec |
| Schema validation at input | No | Yes (JSON Schema Draft-07) |
| Results format | GUI export, CSV/MAT | XDMF + CSV + JSON (automatic) |
| Scriptable without installation | No (requires COMSOL) | Yes (any HTTP client) |
| AI assistant / LLM compatible | No | Yes (JSON readable/writable) |
| Reproducibility guarantee | COMSOL version + .mph | JSON file + `pip install SemiSim==x.y.z` |

**Key points**
- The JSON-as-contract approach is the single biggest architectural
  differentiator.
- A `.mph` file requires COMSOL to open; a JSON file requires nothing.
- AI assistants (LLMs, Copilot agents) can inspect and generate
  SemiSim device specs without executing any code.
- Reproducibility in COMSOL is a documentation burden; in SemiSim
  it is automatic.

---

## Slide 8.5 — Extensibility and Open-Source Access

**[8:00]**

**COMSOL extensibility:**
- COMSOL is closed source. The semiconductor physics kernel is not
  inspectable.
- Users can add custom PDEs via the COMSOL PDE Interfaces, but this
  requires building on top of a proprietary abstraction layer.
- Physics extensions are not unit-testable outside the full COMSOL
  environment.
- A new physics module (e.g., a custom recombination term) requires
  COMSOL expertise and cannot be unit-tested with MMS outside COMSOL.

**SemiSim extensibility:**
- MIT license. Every line of the solver, weak form, and boundary
  condition is inspectable and modifiable.
- A new physics term is a new UFL weak-form expression in
  `semi/physics/`. It can be unit-tested with MMS using
  `pytest` without any GUI or license.
- The five-layer architecture means a new recombination model goes in
  Layer 4 (FEM), has no effect on Layers 1–3, and is gated by a new
  MMS test.
- Runnable on Google Colab with zero local install (Docker image or
  the Colab-FEniCSx setup script).

| Extensibility dimension | COMSOL | SemiSim |
|---|---|---|
| Source code inspectable | No | Yes (MIT) |
| Custom PDE terms | Yes (proprietary interface) | Yes (standard UFL Python) |
| MMS-testable physics extensions | Not outside COMSOL | Yes (pytest) |
| Zero-install execution | No | Yes (Colab / Docker) |
| License cost | $5k–$25k/seat/year (typical) | Free |
| CI integration | Not practical | Yes (GitHub Actions) |
| AI agent can write/verify extensions | No | Yes |

**Key points**
- For research groups, the cost differential is often $20k+/year.
- For teaching, Colab + zero-install is the enabling feature.
- For AI-assisted development, an LLM can read the weak form and
  propose extensions; this is not possible with COMSOL's closed kernel.

---

## Slide 8.6 — Verification and Validation Practice

**[10:00]**

This is a dimension where SemiSim makes a stronger engineering
claim than COMSOL's default user workflow.

**COMSOL V&V practice:**
- COMSOL ships a library of "Application Library" examples with
  published reference results.
- Convergence studies can be run manually through the GUI but are not
  automated.
- There is no built-in MMS (Method of Manufactured Solutions)
  framework — users must construct and inject forcing terms manually.
- Coverage metrics for the solver source are not published.
- A change to COMSOL's solver is not accompanied by a public regression
  test suite.

**SemiSim V&V practice:**
- Every shipped device capability has a registered analytical verifier
  in the test suite. Verification is not optional.
- The MMS suite covers Poisson 1D/2D, multi-region Poisson, and three
  variants of the coupled drift-diffusion system. Finest-pair L² rates
  land at the theoretical 2.000.
- Conservation checks run at every bias step: total current is constant
  in space; total charge integrates to zero at equilibrium.
- 95% code coverage gate enforced by CI on every push to main.
- The V&V suite is public, runnable, and reproducible from the
  repository without any license.

| V&V dimension | COMSOL | SemiSim |
|---|---|---|
| Bundled benchmark library | Yes (App Library) | Yes (benchmarks/) |
| Automated convergence (MMS) suite | No | Yes |
| Conservation checks per solve | No | Yes (automated in CI) |
| Coverage metric | Not published | 95% gate, CI-enforced |
| Reproducible from public repo | No (requires license) | Yes |
| Regression-tested on every commit | No | Yes (GitHub Actions) |

**Key points**
- COMSOL's App Library is validation; SemiSim's MMS suite is
  verification. Both are needed; only SemiSim runs both
  automatically.
- The 95% coverage gate means a change that breaks an untested code
  path must add a test before merging.
- Every result in this talk can be reproduced by cloning the repo and
  running `pytest`.

---

## Slide 8.7 — Summary: When to Use Which Tool

**[11:30]**

Both tools are legitimate; the right choice depends on the problem.

**Use COMSOL Semiconductor Module when:**
- You need quantum corrections (density-gradient) or
  electrothermal coupling — these are mature in COMSOL and not yet
  in SemiSim.
- You are working with complex 3D device geometry and need a full
  parametric CAD kernel with GUI meshing.
- Your team already has COMSOL licenses and institutional support.
- You are working on a proprietary product where the cost of COMSOL
  is acceptable and the source code opacity is not a concern.
- You need advanced optical or RF co-simulation through COMSOL's
  multiphysics coupling.

**Use SemiSim when:**
- You need full transparency into the solver: every weak form, every
  Newton step, every residual is inspectable and modifiable.
- You are a graduate student or research group without COMSOL funding.
- You need to script and automate device studies from a web API,
  CI pipeline, or AI assistant.
- You are teaching semiconductor device physics and want students to
  run simulations on Colab with zero install.
- You want your simulations to be version-controlled, reproducible
  from a plain JSON file, and regression-tested.
- You are implementing a new physics model and need MMS-level
  verification against a known exact solution.
- You want a GPU linear-solver path that operates transparently
  without changing input files.

| Decision factor | Choose COMSOL | Choose SemiSim |
|---|---|---|
| Quantum corrections / electrothermal | ✓ | — |
| Complex CAD geometry (GUI-driven) | ✓ | — |
| Classical silicon DD physics | ✓ | ✓ |
| Open-source, inspectable solver | — | ✓ |
| Zero-cost access / Colab | — | ✓ |
| Scriptable from HTTP / AI | — | ✓ |
| Automated MMS verification | — | ✓ |
| GPU linear-solve path | — | ✓ |
| Version-controlled JSON reproducibility | — | ✓ |

**Key points**
- COMSOL's advantages are its physics breadth, CAD kernel, and
  institutional maturity.
- SemiSim's advantages are transparency, cost, reproducibility,
  extensibility, and open-source AI-compatible workflow.
- These tools are complementary, not mutually exclusive. A research
  group might prototype in SemiSim (free, scriptable) and validate
  against COMSOL for the final physical result.

**Transition:** For a quick-reference glossary of all terms and
acronyms used in this talk, see the Glossary section.
