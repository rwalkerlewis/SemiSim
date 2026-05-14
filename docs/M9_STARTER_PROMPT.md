# Starter prompt: M9 — Result Artifact Writer

Paste this into Claude Code, Cursor, or whatever agent you use. It is
self-contained; the agent should not need to ask clarifying questions.

---

You are working on `SemiSim` (https://github.com/rwalkerlewis/SemiSim), a
FEniCSx-based JSON-driven finite-element semiconductor device simulator.

Before writing any code, read the following repository files in this order:

1. `PLAN.md`
2. `docs/ROADMAP.md`
3. `docs/ARCHITECTURE.md`
4. `docs/IMPROVEMENT_GUIDE.md` — **your task is M9 in that document**
5. `semi/run.py`
6. `semi/runners/bias_sweep.py`
7. `semi/runners/equilibrium.py`
8. `semi/runners/mos_cv.py`
9. `semi/schema.py`

Do not modify any physics kernels or existing tests. Your changes are purely
additive: a new module, a new CLI entry point, a new schema file, and a new
test file.

## Task

Implement M9 as specified in `docs/IMPROVEMENT_GUIDE.md` section 4:

1. Create `schemas/manifest.v1.json` — a JSON-schema Draft-07 file describing
   the manifest contract in section 3 of the improvement guide. Every property
   must have `description`. Required fields: `schema_version`, `engine`,
   `run_id`, `status`, `wall_time_s`, `input_sha256`, `solver`, `fields`,
   `mesh`. Optional: `sweeps`, `warnings`.

2. Create `semi/io/__init__.py` and `semi/io/artifact.py` exposing:

   ```python
   def write_artifact(
       result: "SimulationResult",
       out_dir: Path,
       run_id: str | None = None,
       input_json_path: Path | None = None,
   ) -> Path:
       """
       Write the full result tree under out_dir/<run_id>/ and return the path.
       If run_id is None, generate one as
           f"{iso_timestamp}_{cfg['name']}_{short_sha}"
       where short_sha is the first 7 chars of sha256(input_json_bytes).
       """

   def read_manifest(run_dir: Path) -> dict:
       """Load and JSON-schema-validate manifest.json from a run dir."""
   ```

   Output layout exactly as specified in section 3 of the improvement guide.

3. Field export rules:
   - Write `psi`, `n`, `p`, `N_net` always (all solver types produce these).
   - Write `E = -grad(psi)` as a DG0 vector field.
   - For `drift_diffusion` and `bias_sweep` runs, additionally write `phi_n`,
     `phi_p`, and per-step `J_n`, `J_p` (these can be time-series in a single
     BP file).
   - Use `dolfinx.io.VTXWriter` with ADIOS2 BP5 backend. If ADIOS2 is not
     available, fall back to `XDMFFile` and note the fallback in
     `manifest.warnings`.

4. IV export: if `result.iv` is non-empty, write `iv/<contact>.csv` with header
   `V,J_n,J_p,J_total`. For now, populate `J_total = result.iv[i]["J"]` and
   leave `J_n`, `J_p` as `NaN` — splitting those is M10's problem.

5. Convergence export: write `convergence/snes.csv` with columns
   `bias_step,V_applied,iterations,reason,converged` from the per-step info
   the runners already record. For equilibrium runs, this is a single row.

6. New CLI `semi-run`, registered via `pyproject.toml`:

   ```
   semi-run <path-to-input.json> [--out runs/] [--run-id <id>]
   ```

   Prints the run directory path on stdout on success. Nonzero exit on any
   failure.

7. Tests in `tests/test_artifact.py`:
   - For each of the five benchmarks in `benchmarks/*/`, run it through
     `semi.run.run`, call `write_artifact`, then `read_manifest`, and assert
     the manifest validates against `schemas/manifest.v1.json`.
   - Assert every field listed in `manifest.fields[].path` exists on disk.
   - Assert every sweep file listed in `manifest.sweeps[].path` is readable
     as CSV and has `n_steps + 1` rows (header + data).
   - Assert `manifest.input_sha256` matches the actual SHA256 of the input
     JSON file.
   - Use `pytest.mark.skipif` to skip the FEM-heavy tests when dolfinx is
     unavailable; the pure-Python tests (schema validation round-trip) must
     still run in the pure-Python CI job.

## Constraints

- No changes to `semi/physics/`, `semi/mesh.py`, `semi/bcs.py`,
  `semi/continuation.py`, or any existing test.
- No `import dolfinx` in `semi/io/artifact.py` at module top level. Do the
  imports inside functions, matching the existing convention in `semi/run.py`.
- The manifest JSON must be deterministic given identical input — sort keys,
  fix float precision at 10 significant figures, use ISO-8601 UTC for the
  `run_id` timestamp.
- The `semisim_server/` package from M10 will consume this artifact. Design
  `read_manifest` so it can be called from a subprocess that has no dolfinx
  or numpy installed (pure-stdlib JSON only). If you need to, split the
  reader out into `semi/io/reader.py` that has only stdlib dependencies.

## Acceptance criteria (gate for merge)

Run these commands; all must succeed:

```bash
# 1. New schema is valid JSON-schema Draft-07
python -c "import json, jsonschema; s = json.load(open('schemas/manifest.v1.json')); jsonschema.Draft7Validator.check_schema(s)"

# 2. Existing tests still pass
pytest -q

# 3. New artifact tests pass
pytest -q tests/test_artifact.py

# 4. CLI works end-to-end
semi-run benchmarks/pn_1d/pn_junction.json --out /tmp/runs
ls /tmp/runs/*/manifest.json

# 5. The manifest validates against its own schema
python -c "
import json, jsonschema
from pathlib import Path
run_dir = sorted(Path('/tmp/runs').iterdir())[-1]
m = json.load(open(run_dir / 'manifest.json'))
s = json.load(open('schemas/manifest.v1.json'))
jsonschema.Draft7Validator(s).validate(m)
print('manifest OK:', run_dir)
"

# 6. Lint/format clean
ruff check semi/ tests/ schemas/
```

## What to do when finished

Update `PLAN.md`:
- Move the "M9" line from "Next task" to "Completed work log".
- Set "Current state" to reflect artifact writer shipped.
- Set "Next task" to "M10: HTTP server (see docs/IMPROVEMENT_GUIDE.md §4)".

Update `CHANGELOG.md` with a `[0.9.0] - M9: Result artifact writer` entry
following the format of the existing entries.

Do not bump the version number yourself in `semi/__init__.py`. The maintainer
will cut the release.

## What not to do

- Do not start M10. The HTTP server is a separate PR.
- Do not add new physics, new solver backends, new benchmarks, or new
  materials.
- Do not refactor `semi/run.py` or the runners. If you find you need to, stop
  and ask first.
- Do not change the existing JSON input schema in this PR. Any schema changes
  belong in M11.

When you are done, post a single summary message listing: files added, files
modified (should be only `PLAN.md`, `CHANGELOG.md`, `pyproject.toml`),
acceptance-test output, and any open questions you had to resolve.
