"""
Tests for the M18.1 schema addition `solver.diagnostics`.

Pure-Python; mirrors the M18 adaptive-dt schema tests:
  - the new field default-fills on bias_sweep,
  - non-bias_sweep solver types reject it,
  - existing benchmark and example JSONs validate unchanged, and
  - the current v2.10.0 minor version accepts legacy inputs unchanged
    aside from the version string.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from semi import schema

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"
EXAMPLES_DIR = REPO_ROOT / "examples"


@pytest.fixture
def bias_sweep_cfg():
    return {
        "schema_version": "2.10.0",
        "name": "diagnostics_schema_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [20],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "anode", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0},
            }
        ],
        "contacts": [
            {
                "name": "anode",
                "facet": "anode",
                "type": "ohmic",
                "voltage": 0.0,
                "voltage_sweep": {"start": 0.0, "stop": 0.1, "step": 0.1},
            },
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "solver": {"type": "bias_sweep"},
    }


def _all_existing_configs():
    for root in (BENCHMARKS_DIR, EXAMPLES_DIR):
        for json_path in sorted(root.rglob("*.json")):
            if json_path.name == "manifest.json":
                continue
            yield json_path


@pytest.mark.parametrize(
    "json_path",
    list(_all_existing_configs()),
    ids=lambda p: p.parent.name + "/" + p.name,
)
def test_existing_configs_validate_unchanged(json_path):
    with json_path.open() as f:
        cfg = json.load(f)
    if "schema_version" not in cfg:
        pytest.skip(f"{json_path} has no schema_version")
    schema.validate(cfg)


@pytest.mark.parametrize(
    "json_path",
    list(_all_existing_configs()),
    ids=lambda p: p.parent.name + "/" + p.name,
)
def test_existing_configs_validate_as_v210(json_path):
    with json_path.open() as f:
        cfg = json.load(f)
    if "schema_version" not in cfg:
        pytest.skip(f"{json_path} has no schema_version")
    cfg["schema_version"] = "2.10.0"
    schema.validate(cfg)


def test_diagnostics_default_fills_false(bias_sweep_cfg):
    result = schema.validate(copy.deepcopy(bias_sweep_cfg))
    assert result["solver"]["diagnostics"] is False


def test_diagnostics_true_validates(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["diagnostics"] = True
    result = schema.validate(bias_sweep_cfg)
    assert result["solver"]["diagnostics"] is True


def test_diagnostics_on_non_bias_sweep_rejected(bias_sweep_cfg):
    cfg = copy.deepcopy(bias_sweep_cfg)
    cfg["solver"] = {"type": "equilibrium", "diagnostics": True}
    cfg["contacts"][0].pop("voltage_sweep", None)
    with pytest.raises(schema.SchemaError, match="only consumed by the bias_sweep runner"):
        schema.validate(cfg)
