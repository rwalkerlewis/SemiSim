"""
Tests for the M11 schema-versioning contract.

These are pure-Python tests (no dolfinx needed) covering:
    - The standalone JSON schema file is a valid Draft-07 schema.
    - Every `type: object` node in the input schema has a description.
    - Every benchmark JSON declares schema_version.
    - validate() rejects mismatched major versions.
    - validate() accepts minor/patch skew silently.
    - validate() rejects missing schema_version.
    - The semi.schema loader caches the parsed schema.
    - read_manifest rejects mismatched manifest major versions.
"""
from __future__ import annotations

import glob
import json
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
SCHEMAS_DIR = REPO_ROOT / "schemas"
INPUT_SCHEMA_PATH = SCHEMAS_DIR / "input.v1.json"
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"


def _base_cfg() -> dict:
    """A minimal cfg that passes jsonschema validation; callers tweak schema_version."""
    return {
        "schema_version": "1.0.0",
        "name": "sv_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [100],
            "facets_by_plane": [
                {"name": "L", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "R", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {"si": {"material": "Si", "tag": 1, "role": "semiconductor"}},
        "doping": [
            {"region": "si", "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "L", "facet": "L", "type": "ohmic", "voltage": 0.0},
            {"name": "R", "facet": "R", "type": "ohmic", "voltage": 0.0},
        ],
    }


def test_schema_file_is_valid_draft7():
    import jsonschema

    schema = json.loads(INPUT_SCHEMA_PATH.read_text())
    jsonschema.Draft7Validator.check_schema(schema)


def test_every_object_node_has_description():
    schema = json.loads(INPUT_SCHEMA_PATH.read_text())
    offenders: list[str] = []

    def walk(node, path: str) -> None:
        if not isinstance(node, dict):
            return
        if node.get("type") == "object":
            if not node.get("description"):
                offenders.append(path)
        for key in ("properties", "patternProperties"):
            children = node.get(key)
            if isinstance(children, dict):
                for child_name, child in children.items():
                    walk(child, f"{path}.{key}.{child_name}")
        ap = node.get("additionalProperties")
        if isinstance(ap, dict):
            walk(ap, f"{path}.additionalProperties")
        for combinator in ("oneOf", "anyOf", "allOf"):
            branches = node.get(combinator)
            if isinstance(branches, list):
                for i, branch in enumerate(branches):
                    walk(branch, f"{path}.{combinator}[{i}]")
        items = node.get("items")
        if isinstance(items, dict):
            walk(items, f"{path}.items")
        elif isinstance(items, list):
            for i, sub in enumerate(items):
                walk(sub, f"{path}.items[{i}]")

    walk(schema, "(root)")
    assert offenders == [], (
        "Every type:object node must have a description. Missing: "
        + ", ".join(offenders)
    )


def test_schema_version_required_in_benchmarks():
    pattern = re.compile(r"^\d+\.\d+\.\d+$")
    paths = sorted(glob.glob(str(BENCHMARKS_DIR / "*" / "*.json")))
    assert paths, "no benchmark JSONs found under benchmarks/*/*.json"
    for p in paths:
        cfg = json.loads(Path(p).read_text())
        assert "schema_version" in cfg, f"{p} is missing schema_version"
        sv = cfg["schema_version"]
        assert isinstance(sv, str), f"{p}: schema_version must be a string, got {type(sv)}"
        assert pattern.match(sv), (
            f"{p}: schema_version {sv!r} must match MAJOR.MINOR.PATCH"
        )


def test_validate_rejects_major_mismatch():
    from semi import schema as schema_mod

    cfg = _base_cfg()
    # v1 and v2 are both supported as of M14.3; pick a major outside the
    # supported range to exercise the rejection path.
    cfg["schema_version"] = "3.0.0"
    with pytest.raises(schema_mod.SchemaError, match=r"major"):
        schema_mod.validate(cfg)


def test_validate_accepts_minor_skew():
    from semi import schema as schema_mod

    cfg = _base_cfg()
    cfg["schema_version"] = "1.99.3"
    result = schema_mod.validate(cfg)
    assert result["schema_version"] == "1.99.3"


def test_validate_rejects_missing_schema_version():
    from semi import schema as schema_mod

    cfg = _base_cfg()
    del cfg["schema_version"]
    with pytest.raises(schema_mod.SchemaError):
        schema_mod.validate(cfg)


def test_schema_loader_caches():
    import semi.schema as schema_mod

    a = schema_mod.SCHEMA
    b = schema_mod.SCHEMA
    assert a is b


def test_read_manifest_rejects_major_mismatch(tmp_path):
    from semi.io.reader import read_manifest

    manifest = {
        "schema_version": "2.0.0",
        "engine": {"name": "SemiSim", "version": "0.11.0", "commit": "abc1234"},
        "run_id": "2026-04-23T00-00-00Z_test_abc1234",
        "status": "completed",
        "wall_time_s": 0.0,
        "input_sha256": "a" * 64,
        "solver": {
            "type": "equilibrium",
            "backend": "petsc-cpu-mumps",
            "n_dofs": 1,
            "n_steps": 1,
            "converged": True,
        },
        "fields": [
            {"name": "psi", "units": "V", "path": "fields/psi.bp", "rank": 0},
        ],
        "mesh": {
            "path": "mesh/mesh.xdmf",
            "dimension": 1,
            "n_cells": 1,
            "n_vertices": 2,
            "regions": [{"name": "silicon", "material": "Si", "tag": 1}],
        },
        "warnings": [],
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match=r"major"):
        read_manifest(tmp_path)
