# §1 — Opening: What Problem Are We Solving?

**Suggested slide title:** "Building a Semiconductor Device Simulator from Scratch"
**Target time:** 2–3 minutes

---

## Slide 1.1 — Hook

**[0:00]** Open with a physical question, not a technology statement.

Every chip you have ever used depends on a tiny piece of silicon where
electrons and holes are pushed around by electric fields in precisely
controlled ways. Before you tape out a chip you need to know what that
silicon is going to do. You can't build a thousand chips and measure which
one works. You simulate.

The tools the industry uses for that simulation — COMSOL Semiconductor
Module, Sentaurus, Silvaco — are expensive, closed-source, and designed
for expert TCAD engineers. If you are a graduate student, an open-source
hardware researcher, or an AI assistant writing device models, you either
pay license fees or work around the tools.

**Key points**
- Device simulation is a necessary step between design and fabrication.
- Commercial TCAD tools are expensive and opaque.
- There is a gap for a transparent, extensible, JSON-driven alternative.

**Transition:** Let me tell you what SemiSim actually is, and then
explain why building it requires solving three genuinely hard problems.

---

## Slide 1.2 — What is SemiSim?

**[1:30]**

SemiSim is a finite-element semiconductor device simulator built on
FEniCSx — the leading open-source FEM framework — and driven entirely by
JSON input files. You describe your device: the geometry, the doping
profile, the contacts, the physics options. The engine produces a
self-consistent solution for the electric potential, electron density, and
hole density everywhere in the device, plus the currents flowing through
every contact. Everything is reproducible from the JSON file and the
package version. The output lands on disk as schema-validated artifacts
that any tool can read.

Today, version 0.16, the engine ships:

- Equilibrium Poisson in 1D, 2D, and 3D
- Coupled drift-diffusion bias sweeps
- A 2D MOSFET with analytical verification
- Transient (time-dependent) solves with BDF1 and BDF2
- Small-signal AC analysis for frequency-dependent capacitance
- An axisymmetric 2D MOS capacitor
- An optional GPU linear-solver path (PETSc CUDA/HIP via AMGX or
  hypre BoomerAMG)
- An HTTP API so a UI or another tool can drive solves remotely

All of this verified against analytical references and Method of
Manufactured Solutions, with a CI pipeline that enforces a 95% coverage
gate on every commit.

**Key points**
- JSON in → validated artifact tree out. One file fully determines a simulation.
- Ships 1D, 2D, 3D devices; equilibrium, bias sweep, transient, AC, GPU.
- Open source, MIT license, runnable on Colab with zero local install.

**Transition:** Before we get to the code, let me give you just enough
physics to understand why the solver looks the way it does.
