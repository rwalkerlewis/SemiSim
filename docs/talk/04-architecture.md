# §4 — Architecture: Five Layers and a JSON Contract

**Suggested slide title:** "Architecture: Five Layers, One Contract"
**Target time:** 4–5 minutes

---

## Slide 4.1 — The Five-Layer Design

**[16:00]**

SemiSim is organized as five layers. Higher layers depend on lower
ones; lower layers never import from higher ones. The pure-Python core
(Layer 3) has an additional constraint: it must not import dolfinx. That
lets schema, material, and doping logic run in any Python environment
without a FEM toolchain.

```
Layer 5: Delivery surface
  benchmarks/  scripts/  notebooks/  tests/  Dockerfile

Layer 4: FEM
  semi/mesh.py  semi/physics/*  semi/solver.py  semi/run.py
  depends on: dolfinx 0.10, PETSc, UFL, mpi4py + Layer 3

Layer 3: Pure-Python core
  semi/constants.py  semi/materials.py  semi/scaling.py  semi/doping.py
  depends on: numpy only

Layer 2: Schema
  semi/schema.py
  depends on: jsonschema, Layer 3

Layer 1: JSON input
  benchmarks/*/*.json, user-provided files
  depends on: nothing
```

The import rule is enforced by the test suite: any test import of a Layer
3 module that accidentally pulls in dolfinx will fail in the GitHub
Actions "offline" CI job, which runs without a FEM container.

**Key points**
- Five layers; strict one-way dependency direction.
- Layer 3 (constants, materials, scaling, doping) is pure Python + numpy only.
- This means schema validation and material DB work standalone without FEM.
- The layering prevents circular dependencies and lets CI test each layer separately.

---

## Slide 4.2 — JSON as the Contract

**[18:00]**

The top-level design decision is that JSON is the only supported input
format (ADR 0001). A single JSON file fully determines a simulation. Why
JSON over Python DSL, TOML, or YAML?

- **Schema validation at the boundary.** `jsonschema` Draft-07 catches
  malformed inputs before any FEM code runs. The strict-mode schema
  (`additionalProperties: false`) means a typo in a field name is a
  validation error, not a silently dropped option.
- **Round-trip through any tool.** JSON goes cleanly through web APIs,
  JavaScript UIs, AI assistants, and test fixtures without a Python
  interpreter.
- **Reproducibility.** A JSON file plus a package version fully
  determines the simulation. No hidden mutable state.
- **AI assistant compatibility.** AI agents can inspect and generate
  device specs statically without executing code.

The active schema is `schemas/input.v2.json` (strict, v2.0.0). The legacy
`schemas/input.v1.json` is accepted with a `DeprecationWarning` for one
minor release cycle. Both use JSON Schema Draft-07 with `description`
annotations on every field so a UI can auto-generate forms.

**Key points**
- JSON contract: one file per simulation; schema-validated before any FEM code runs.
- Strict mode (`additionalProperties: false`) catches typos at input time.
- Schema versioned separately from package (v2.0.0 strict; v1.x deprecated).
- Every field annotated with `description` for UI form generation.

---

## Slide 4.3 — The HTTP API and Why the Server Never Imports dolfinx

**[20:00]**

Beyond the Python API, SemiSim ships an HTTP server
(`semisim_server/`) powered by FastAPI. You `POST /solve` with a JSON body
and get back a `run_id`. You poll `GET /runs/{id}` for status and fetch
the manifest, IV CSVs, and field files when complete. Progress streams
over a WebSocket.

One architectural curiosity: the server process **never imports dolfinx,
UFL, or PETSc at module scope**. All FEM work happens in worker
subprocesses spawned via `ProcessPoolExecutor`. This means:

- The server starts and serves `/health`, `/ready`, `/capabilities`, and
  `/schema` in a pure-Python environment without the FEM stack installed.
- UI developers can run a mock server on their laptop without building
  FEniCSx.
- The `GET /capabilities` endpoint reports available GPU backends so a UI
  can gate the GPU option dynamically on the host's hardware.

**Key points**
- HTTP API: POST /solve → run_id; GET /runs/{id}/manifest; WebSocket progress.
- Server process is pure Python; no dolfinx at module scope.
- FEM work in spawned subprocesses — same isolation pattern as a job queue.
- `/capabilities` reports GPU availability; UI can gate options dynamically.

**Transition:** Now let me walk you through what the engine can actually
do today.
