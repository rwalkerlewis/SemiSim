"""
Tests for the M9 result artifact writer.

FEM-heavy tests (requiring dolfinx) are skipped when dolfinx is
unavailable, matching the pattern in tests/fem/conftest.py.
Pure-Python tests (schema validation round-trip) always run.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from semi.results import AcSweepResult, TransientResult

try:
    import dolfinx  # noqa: F401
    HAS_DOLFINX = True
except ImportError:
    HAS_DOLFINX = False

REPO_ROOT = Path(__file__).parent.parent
SCHEMAS_DIR = REPO_ROOT / "schemas"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"

BENCHMARK_JSONS = [
    BENCHMARKS_DIR / "pn_1d" / "pn_junction.json",
    BENCHMARKS_DIR / "pn_1d_bias" / "pn_junction_bias.json",
    BENCHMARKS_DIR / "pn_1d_bias_reverse" / "pn_junction_bias_reverse.json",
    BENCHMARKS_DIR / "mos_2d" / "mos_cap.json",
    BENCHMARKS_DIR / "resistor_3d" / "resistor.json",
]


# ---- Pure-Python tests (always run) ----

def test_manifest_schema_valid_draft7():
    import jsonschema
    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())
    jsonschema.Draft7Validator.check_schema(schema)


def test_manifest_schema_round_trip():
    """A minimal valid manifest validates without error."""
    import jsonschema
    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())
    manifest = {
        "schema_version": "1.0.0",
        "engine": {"name": "kronos-semi", "version": "0.9.0", "commit": "abc1234abc1234"},
        "run_id": "2026-04-23T12-00-00Z_test_abc1234",
        "status": "completed",
        "wall_time_s": 1.5,
        "input_sha256": "a" * 64,
        "solver": {
            "type": "equilibrium",
            "backend": "petsc-cpu-mumps",
            "n_dofs": 401,
            "n_steps": 1,
            "converged": True,
        },
        "fields": [
            {"name": "psi", "units": "V", "path": "fields/psi.bp", "rank": 0},
        ],
        "mesh": {
            "path": "mesh/mesh.xdmf",
            "dimension": 1,
            "n_cells": 400,
            "n_vertices": 401,
            "regions": [{"name": "silicon", "material": "Si", "tag": 1}],
        },
        "warnings": [],
    }
    jsonschema.Draft7Validator(schema).validate(manifest)


def test_manifest_schema_optional_sweeps():
    """Manifest with sweeps validates without error."""
    import jsonschema
    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())
    manifest = {
        "schema_version": "1.0.0",
        "engine": {"name": "kronos-semi", "version": "0.9.0", "commit": "abc1234"},
        "run_id": "2026-04-23T12-00-00Z_test_abc1234",
        "status": "completed",
        "wall_time_s": 5.0,
        "input_sha256": "b" * 64,
        "solver": {
            "type": "bias_sweep",
            "backend": "petsc-cpu-mumps",
            "n_dofs": 801,
            "n_steps": 13,
            "converged": True,
        },
        "fields": [
            {"name": "psi", "units": "V", "path": "fields/psi.bp", "rank": 0},
        ],
        "mesh": {
            "path": "mesh/mesh.xdmf",
            "dimension": 1,
            "n_cells": 800,
            "n_vertices": 801,
            "regions": [{"name": "silicon", "material": "Si", "tag": 1}],
        },
        "sweeps": [
            {
                "kind": "voltage",
                "contact": "anode",
                "path": "iv/anode.csv",
                "n_steps": 13,
                "bipolar": False,
            }
        ],
        "warnings": ["test warning"],
    }
    jsonschema.Draft7Validator(schema).validate(manifest)


def test_write_transient_artifact_writes_per_contact_csv_and_manifest(tmp_path):
    from semi.io.artifact import write_transient_artifact

    input_json = tmp_path / "input.json"
    input_cfg = {"name": "tiny_transient"}
    input_json.write_text(json.dumps(input_cfg))

    result = TransientResult(
        t=[0.0, 1.0e-12],
        iv=[
            {"t": 0.0, "contact": "anode", "V": 0.0, "J_n": 1.0, "J_p": -0.2, "J": 0.8},
            {
                "t": 1.0e-12,
                "contact": "anode",
                "V": 0.1,
                "J_n": 1.1,
                "J_p": -0.1,
                "J_total": 1.0,
            },
            {"t": 0.0, "contact": "cathode", "V": 0.0, "J_n": -1.0, "J_p": 0.2, "J": -0.8},
        ],
        meta={"order": 2, "dt": 1.0e-12, "n_steps_taken": 2},
    )

    run_dir = write_transient_artifact(
        result=result,
        out_dir=tmp_path / "runs",
        run_id="transient_test",
        input_json_path=input_json,
    )

    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["solver"]["type"] == "transient"
    assert manifest["solver"]["order"] == 2
    assert manifest["solver"]["n_steps"] == 2
    assert manifest["solver"]["dt"] == pytest.approx(1.0e-12)
    assert manifest["input_sha256"] == hashlib.sha256(input_json.read_bytes()).hexdigest()
    assert len(manifest["sweeps"]) == 2

    anode_csv = run_dir / "iv" / "anode.csv"
    cathode_csv = run_dir / "iv" / "cathode.csv"
    assert anode_csv.exists()
    assert cathode_csv.exists()

    with open(anode_csv, newline="") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["t", "V", "J_n", "J_p", "J_total"]
    # Row 2 comes from J_total fallback when J is not present.
    assert rows[2][-1] == "1.0"


def test_write_ac_sweep_artifact_writes_csv_and_manifest(tmp_path):
    from semi.io.artifact import write_ac_sweep_artifact

    input_json = tmp_path / "input_ac.json"
    input_cfg = {"name": "tiny_ac"}
    input_json.write_text(json.dumps(input_cfg))

    result = AcSweepResult(
        frequencies=[1.0, 10.0],
        Y=[1.0 + 2.0j, 2.0 + 4.0j],
        Z=[0.0j, 0.0j],
        C=[3.0, 6.0],
        G=[1.0, 2.0],
        dc_bias={"contact": "gate", "voltage": 0.25},
    )

    run_dir = write_ac_sweep_artifact(
        result=result,
        out_dir=tmp_path / "runs",
        run_id="ac_test",
        input_json_path=input_json,
    )

    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["solver"]["type"] == "ac_sweep"
    assert manifest["solver"]["n_freqs"] == 2
    assert manifest["solver"]["dc_contact"] == "gate"
    assert manifest["solver"]["dc_voltage"] == pytest.approx(0.25)
    assert manifest["sweeps"][0]["path"] == "iv/gate.csv"

    ac_csv = run_dir / "iv" / "gate.csv"
    with open(ac_csv, newline="") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["f_Hz", "C", "G", "Re_Y", "Im_Y"]
    assert rows[1] == ["1.0", "3.0", "1.0", "1.0", "2.0"]


def test_write_ac_sweep_artifact_empty_frequency_has_no_sweeps(tmp_path):
    from semi.io.artifact import write_ac_sweep_artifact

    result = AcSweepResult(frequencies=[], Y=[], Z=[], C=[], G=[], dc_bias={})
    run_dir = write_ac_sweep_artifact(
        result=result,
        out_dir=tmp_path / "runs",
        run_id="ac_empty_test",
        input_json_path=None,
    )

    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["solver"]["n_freqs"] == 0
    assert "sweeps" not in manifest


# ---- FEM-heavy tests (require dolfinx) ----

@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx not available")
@pytest.mark.parametrize(
    "bench_json",
    BENCHMARK_JSONS,
    ids=[p.parent.name for p in BENCHMARK_JSONS],
)
def test_write_read_artifact(bench_json, tmp_path):
    """Write artifact for a benchmark and assert manifest validity."""
    import jsonschema

    from semi.io.artifact import write_artifact
    from semi.io.reader import read_manifest
    from semi.run import run
    from semi.schema import load as schema_load

    schema = json.loads((SCHEMAS_DIR / "manifest.v1.json").read_text())

    cfg = schema_load(str(bench_json))
    result = run(cfg)
    run_dir = write_artifact(result, tmp_path / "runs", input_json_path=bench_json)

    # 1. Manifest validates against the schema
    manifest = read_manifest(run_dir)
    jsonschema.Draft7Validator(schema).validate(manifest)

    # 2. All field files listed in manifest exist on disk
    for field_entry in manifest["fields"]:
        field_path = run_dir / field_entry["path"]
        assert field_path.exists(), f"Field file not found: {field_entry['path']}"

    # 3. Sweep CSVs have n_steps + 1 rows (header + data)
    for sweep in manifest.get("sweeps", []):
        sweep_path = run_dir / sweep["path"]
        assert sweep_path.exists(), f"Sweep CSV not found: {sweep['path']}"
        with open(sweep_path, newline="") as f:
            rows = list(csv.reader(f))
        assert len(rows) == sweep["n_steps"] + 1, (
            f"Expected {sweep['n_steps'] + 1} rows in {sweep['path']}, got {len(rows)}"
        )

    # 4. input_sha256 matches the actual file
    expected = hashlib.sha256(bench_json.read_bytes()).hexdigest()
    assert manifest["input_sha256"] == expected, "input_sha256 mismatch"


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx not available")
def test_write_artifact_cfg_sha_when_no_input_path(tmp_path):
    """When input_json_path is omitted, input_sha256 is derived from result.cfg."""
    from semi.io.artifact import write_artifact
    from semi.run import run
    from semi.schema import load as schema_load

    bench_json = BENCHMARKS_DIR / "pn_1d" / "pn_junction.json"
    cfg = schema_load(str(bench_json))
    result = run(cfg)
    run_dir = write_artifact(result, tmp_path / "runs", input_json_path=None)

    manifest = json.loads((run_dir / "manifest.json").read_text())
    expected = hashlib.sha256(json.dumps(result.cfg, sort_keys=True).encode()).hexdigest()
    assert manifest["input_sha256"] == expected
