# Architecture

SemiSim is a FEniCSx-based finite-element semiconductor device simulator driven by a single JSON file (optionally referencing an external geometry/mesh artifact). It is organized as five layers. Higher layers depend on lower layers; lower layers must never import from higher layers. The pure-Python core (Layer 3) has an additional constraint: it must not import dolfinx, so that schema, material, and doping logic can be tested in a lightweight CI environment and used standalone.

For the project mission and scope relative to the COMSOL Semiconductor Module, see [docs/ROADMAP.md](ROADMAP.md) and [README.md](../README.md).

## Layer diagram

```
+-----------------------------------------------------------------------+
| Layer 5 : Delivery surface                                            |
|   benchmarks/  scripts/  notebooks/  tests/  Dockerfile               |
|   docker-compose.yml  .devcontainer/  .github/workflows/              |
+-----------------------------^-----------------------------------------+
                              |
+-----------------------------+-----------------------------------------+
| Layer 4 : FEM                                                         |
|   semi/mesh.py   semi/physics/*   semi/solver.py   semi/run.py        |
|   depends on: dolfinx 0.10, petsc4py, ufl, mpi4py + Layer 3           |
+-----------------------------^-----------------------------------------+
                              |
+-----------------------------+-----------------------------------------+
| Layer 3 : Pure-Python core                                            |
|   semi/constants.py   semi/materials.py   semi/scaling.py             |
|   semi/doping.py                                                      |
|   depends on: numpy only                                              |
+-----------------------------^-----------------------------------------+
                              |
+-----------------------------+-----------------------------------------+
| Layer 2 : Schema                                                      |
|   semi/schema.py                                                      |
|   depends on: jsonschema, Layer 3 (for defaults)                      |
+-----------------------------^-----------------------------------------+
                              |
+-----------------------------+-----------------------------------------+
| Layer 1 : JSON input                                                  |
|   benchmarks/*/*.json, user-provided files                            |
|   depends on: nothing                                                 |
+-----------------------------------------------------------------------+
```

## Layer-by-layer

### Layer 1: JSON input

- **Provides:** a declarative device specification (mesh, regions,
  doping, contacts, physics options, solver options, output).
- **Depends on:** nothing. JSON is a data format, not code.
- **Must never depend on:** Python, imports, or side effects. A JSON
  file must be consumable by non-Python tooling (UI, JS, other
  solvers' converters).

### Layer 2: Schema (`semi.schema`)

- **Provides:** `validate(cfg)` and `load(path)` that produce a
  normalized config dict with defaults filled in and unknown keys
  rejected.
- **Depends on:** jsonschema, Layer 3 constants for unit conversions
  used in defaults.
- **Must never depend on:** dolfinx, petsc4py, matplotlib, numpy for
  numerics (only numpy types are fine, but no heavy computation).

### Layer 3: Pure-Python core

Modules: `semi.constants`, `semi.materials`, `semi.scaling`,
`semi.doping`.

- **Provides:**
  - `constants`: SI physical constants and cm/m conversion helpers.
  - `materials`: material parameter database (Si, Ge, GaAs, SiO2,
    HfO2, Si3N4) as `Material` dataclasses.
  - `scaling`: `Scaling` class with `L0`, `C0`, `V0`, `lambda2`,
    `debye_length`, etc. Derives scales from the config.
  - `doping`: callable doping-profile evaluators
    (uniform, step, gaussian) that take an `(dim, N)` array and return
    net doping in $\mathrm{m}^{-3}$.
- **Depends on:** numpy only.
- **Must never depend on:** dolfinx, petsc4py, mpi4py, ufl, basix,
  matplotlib, jupyter. These imports are lazy in Layer 4 specifically
  so that Layer 3 can be imported and unit-tested in a plain Python
  environment.

Rationale: fast CI, no Docker required for schema/materials
development, and users with only the pure-Python install can still
drive the schema and material DB.

### Layer 4: FEM

Modules: `semi.mesh`, `semi.physics.*`, `semi.solver`, `semi.run`.

- **Provides:**
  - `mesh`: builtin interval/rectangle/box mesh construction plus
    region and facet tagging.
  - `physics.poisson`: equilibrium Poisson UFL residual builder (single-region and multi-region).
  - `physics.drift_diffusion`: Slotboom continuity forms (electron and hole continuity as separate form builders, plus a combined block residual).
  - `physics.recombination`: SRH recombination kernel as a UFL expression builder.
  - `solver`: `solve_nonlinear` wrapping `dolfinx.fem.petsc.NonlinearProblem`.
  - `runners/equilibrium.py`, `runners/bias_sweep.py`, `runners/_common.py`: per-solver-type orchestration.
  - `bcs.py`: contact BC resolution and Dirichlet BC builders (pure-Python tier; no dolfinx at module scope).
  - `postprocess.py`: current extraction and plotting helpers.
  - `run`: top-level `run(cfg) -> SimulationResult` dispatcher (74 lines).
- **Depends on:** dolfinx 0.10, petsc4py, mpi4py, ufl, basix, plus
  Layers 1-3.
- **Must never depend on:** matplotlib, jupyter, IPython. Plotting and
  notebook concerns live in Layer 5.

### Layer 5: Delivery surface

- **Provides:**
  - `benchmarks/`: self-contained JSON inputs with a README per benchmark.
  - `scripts/run_benchmark.py`: CLI that loads a benchmark, runs it, and
    verifies it.
  - `notebooks/`: thin Colab drivers that clone the repo and import the
    package.
  - `tests/`: pytest suite covering Layer 2-3 (offline) and Layer 4
    (in-Docker).
  - `Dockerfile`, `docker-compose.yml`, `.devcontainer/`: reproducible
    dev environment.
  - `.github/workflows/ci.yml`: lint and test the pure-Python layers on
    Python 3.10 / 3.11 / 3.12 (no dolfinx required).
- **Depends on:** everything below.
- **Must never depend on:** nothing depends on Layer 5; it is the leaf.

## Test matrix

- **Offline (GitHub Actions, no dolfinx):** Layers 2-3 with the pure-Python pytest suite, plus ruff lint. Python 3.10, 3.11, and 3.12.
- **In-Docker (GHA with dolfinx image):** full suite (206 tests at 95.58% coverage) plus all five benchmark verifiers and the V&V suite (`scripts/run_verification.py all`, 62/62 PASS) inside the `ghcr.io/fenics/dolfinx/dolfinx:stable` image. CI timeout is 15 minutes once the image is cached.

## Module boundaries at a glance

| Concern              | Where it lives         | Forbidden imports                       |
|----------------------|------------------------|------------------------------------------|
| Unit conversions     | `semi.constants`       | everything except stdlib                 |
| Material DB          | `semi.materials`       | dolfinx                                  |
| Scaling              | `semi.scaling`         | dolfinx                                  |
| Doping profiles      | `semi.doping`          | dolfinx                                  |
| JSON schema          | `semi.schema`          | dolfinx, matplotlib                      |
| Mesh assembly        | `semi.mesh`            | matplotlib                               |
| UFL residuals        | `semi.physics.*`       | matplotlib                               |
| SNES driver          | `semi.solver`          | matplotlib                               |
| Top-level runner     | `semi.run`             | matplotlib (plotting is Layer 5)         |
| Plotting             | `scripts/run_benchmark.py`, notebooks | ...                        |
