# Contributing

## Before you start

Read these documents before writing code. They are the source of truth
for what is done, what is next, and which design decisions are locked.

**Required reading (~45 minutes):**

- [PLAN.md](PLAN.md) — current state, single in-flight task, and invariants.
- [docs/IMPROVEMENT_GUIDE.md](docs/IMPROVEMENT_GUIDE.md) — M9+ roadmap with explicit acceptance tests per milestone. Do not start a milestone without reading its section.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — five-layer design and per-layer import rules.

**Depending on what you're doing:**

- **New physics or solver code:** [docs/PHYSICS.md](docs/PHYSICS.md) for equations and scaling conventions; then [docs/adr/](docs/adr/) for locked decisions. If your change conflicts with any invariant in `PLAN.md`, open an ADR first.
- **Code navigation and the data flow:** [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md) traces `semi.run.run(cfg)` end-to-end with file:line anchors.
- **No device-physics background:** [docs/PHYSICS_INTRO.md](docs/PHYSICS_INTRO.md) is a tutorial that assumes calculus and programming but not TCAD. Read this before `PHYSICS.md`.
- **Picking up an M9+ milestone:** the matching section of `IMPROVEMENT_GUIDE.md` has your acceptance tests. For M9 specifically, `docs/M9_STARTER_PROMPT.md` is a ready-to-paste prompt for coding agents.

## Development setup

```bash
git clone https://github.com/rwalkerlewis/SemiSim.git
cd SemiSim
# Option A: conda for dolfinx (full FEM capability)
conda create -n SemiSim -c conda-forge python=3.12 fenics-dolfinx mpich
conda activate SemiSim
pip install -e ".[dev]"

# Option B: pure-Python only (no FEM, just schema/materials/scaling/doping)
pip install -e ".[dev]"
```

## Before committing

```bash
pytest tests/ -v               # ~237 pure-Python tests; FEM-backed tests skip without dolfinx
ruff check semi/ tests/
python tests/check_analytical_math.py
```

If you added FEM code, also run:

```bash
docker compose run --rm test                     # full pytest including FEM tests
python scripts/run_verification.py               # V&V suite (MMS, conservation)
```

## Structure

- `semi/` — the package. Pure-Python modules (`constants`, `materials`, `scaling`, `doping`, `schema`, `cv`) never import dolfinx. FEM modules (`mesh`, `physics`, `solver`, `run`, `runners`, `bcs`, `postprocess`) import dolfinx lazily inside function bodies, not at module top level.
- `benchmarks/` — self-contained JSON inputs plus a per-benchmark README.
- `notebooks/` — thin Colab drivers that clone the repo and import the package.
- `tests/` — pytest suite. Tests at `tests/*.py` may be either pure-Python or FEM-gated; FEM-gated tests use `pytest.importorskip("dolfinx")` at module top so they skip cleanly without dolfinx. Tests under `tests/fem/*.py` always require dolfinx.
- `scripts/` — notebook generators, benchmark runner CLI, V&V driver. Not installed.
- `docs/` — physics reference, walkthrough, architecture, ADRs, theory notes, schema reference, per-benchmark landing pages, historical task prompts under `docs/tasks/`. See [docs/index.md](docs/index.md).
- `schemas/` — JSON schemas for input and manifest contracts (`input.v1.json`, `manifest.v1.json`).

## Conventions

- Internal units are SI. JSON input follows device-physics tradition where doping is in cm⁻³ and mobility in cm²/(V·s), converted at ingest.
- Physics-style variable names (`N_A`, `V`, `F`, `psi_L`, `eps_Si`) are fine; PEP 8 naming is explicitly relaxed via ruff config.
- Any new physics module requires (a) a unit test in `tests/` if it has pure-Python logic, (b) an MMS verifier following ADR 0006, (c) an analytical verification in its benchmark.
- Never skip `make_scaling_from_config`. The raw equations have a Jacobian condition number > 10³⁰.
- For M9+ work, don't cross milestone boundaries in one PR. The acceptance tests in `IMPROVEMENT_GUIDE.md` are non-negotiable gates.
- No em dashes in prose or code comments. Use commas, periods, parentheses, or colons.

## Adding a new benchmark

Each benchmark is a self-contained directory under `benchmarks/<name>/`:

```
benchmarks/<name>/
  <name>.json           # required: schema-validated input config
  README.md             # required: brief device description, expected results
  reference_*.csv       # optional: analytical reference for the verifier
  fixtures/*.msh        # optional: committed gmsh fixtures (see Mesh files below)
  *.geo                 # optional: gmsh source for the mesh (see below)
```

The JSON file must declare `schema_version` (currently `2.0.0`
strict; v1.x.y is accepted but emits a `DeprecationWarning`,
deprecated for one minor cycle); see
[docs/schema/reference.md](docs/schema/reference.md) for the full
contract. Wire the benchmark into the CLI driver
[`scripts/run_benchmark.py`](scripts/run_benchmark.py) by adding a
verifier and (optionally) plotters under the benchmark name, and
register it in [`PLAN.md`](PLAN.md) and the
[CHANGELOG](CHANGELOG.md).

A matching notebook under `notebooks/` and a landing page under
`docs/benchmarks/<name>.md` are required for any benchmark that is
called out in the README's notebook catalog.

## Coordinate systems and axisymmetric benchmarks

Schema 1.3.0 added a top-level `coordinate_system` field with values
`"cartesian"` (default) or `"axisymmetric"`. To write an axisymmetric
benchmark:

1. Set `"coordinate_system": "axisymmetric"` and `"dimension": 2` in the JSON.
2. The first mesh coordinate is interpreted as `r ≥ 0`; the second is `z`.
3. The symmetry axis `r = 0` is a natural (no-flux) boundary; do **not** put a Dirichlet contact there. Cross-field schema validation will reject it.
4. The outer radial wall `r = R` defaults to homogeneous Neumann; choose `R` large enough that the solution is insensitive to the cutoff (`R ≳ 5 W_dmax` for MOS-type devices).
5. Use the axisymmetric weak-form builders in
   [`semi/physics/axisymmetric.py`](semi/physics/axisymmetric.py)
   (the runners dispatch automatically based on `coordinate_system`).

Theory: [docs/theory/axisymmetric.md](docs/theory/axisymmetric.md).
Reference benchmark: [`benchmarks/moscap_axisym_2d/`](benchmarks/moscap_axisym_2d/).

## Mesh files: the gmsh `.geo` template convention

For non-trivial 2D / 3D devices, the meridian / structured mesh is
generated from a committed gmsh source file:

- Place the gmsh source as `<benchmark>.geo` next to the JSON config.
  Example:
  [`benchmarks/moscap_axisym_2d/moscap_axisym.geo`](benchmarks/moscap_axisym_2d/moscap_axisym.geo).
- The committed `.geo` is the source of truth; per-benchmark `.msh`
  files are regenerable and are not committed by default. Generate
  them locally with:

  ```bash
  cd benchmarks/<name>
  gmsh -2 <name>.geo -o <name>.msh   # 2D; use -3 for 3D
  ```

- Small fixtures used by loader tests (e.g.
  `benchmarks/resistor_3d/fixtures/box.msh`) are an exception and
  *are* committed. `.gitignore` whitelists them via
  `!benchmarks/**/fixtures/*.msh`.
- dolfinx 0.10's `dolfinx.io.gmsh.read_from_msh` reads the file
  directly and propagates physical groups as cell and facet tags.

## Pure-Python vs FEM test split

The test suite is split deliberately so that contributors without
dolfinx can still run a meaningful subset:

- **Pure-Python tests** at `tests/*.py` cover schema, materials,
  scaling, doping, MOSCAP analytics (`semi/cv.py`), continuation, and
  the offline math sanity checks.
- **FEM-gated tests** at `tests/fem/*.py` always require dolfinx; if
  dolfinx is not importable, the whole module skips.
- **Mixed-mode tests** that live at `tests/*.py` but exercise dolfinx
  internally use the canonical guard at the top of the file:

  ```python
  import pytest
  pytest.importorskip("dolfinx")
  ```

  Examples:
  [`tests/test_moscap_axisym_cv.py`](tests/test_moscap_axisym_cv.py),
  [`tests/test_moscap_axisym_cv_fem.py`](tests/test_moscap_axisym_cv_fem.py).

When you add a new test, decide deliberately which bucket it belongs
in. Pure-Python tests must not import `semi.runners`, `semi.physics`,
`semi.mesh`, `semi.solver`, `semi.run`, `semi.bcs`, or
`semi.postprocess` (those modules import dolfinx lazily but their
imports still pull dolfinx into the test process when the *modules*
are imported at test-collection time).
