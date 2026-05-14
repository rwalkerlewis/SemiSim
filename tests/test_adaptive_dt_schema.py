"""
Tests for the M18 schema additions: solver.adaptive for the transient
runner. Pure-Python; no dolfinx required.

Mirrors the M16.5 / M16.6 / M16.7 schema-test pattern: every existing
benchmark JSON keeps validating, the new field is optional, and the
loader rejects mutually exclusive or solver-mismatched configs at
validate time so users see the failure before the FEM path.
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
def transient_cfg():
    return {
        "schema_version": "2.9.0",
        "name": "adaptive_dt_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [10],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
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
            {"name": "anode",   "facet": "anode",   "type": "ohmic",
             "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "solver": {
            "type": "transient",
            "t_end": 1.0e-9,
            "dt": 1.0e-10,
        },
    }


def _all_existing_benchmark_configs():
    """Yield every shipped benchmark / example JSON path."""
    for d in (BENCHMARKS_DIR, EXAMPLES_DIR):
        if not d.exists():
            continue
        for json_path in sorted(d.rglob("*.json")):
            if json_path.name == "manifest.json":
                continue
            yield json_path


@pytest.mark.parametrize(
    "json_path",
    list(_all_existing_benchmark_configs()),
    ids=lambda p: p.parent.name + "/" + p.name,
)
def test_existing_benchmark_validates_against_v29(json_path):
    """Every shipped benchmark / example JSON validates unchanged on v2.9.0."""
    with json_path.open() as f:
        cfg = json.load(f)
    if "schema_version" not in cfg:
        pytest.skip(f"{json_path} has no schema_version (probably not an input config)")
    schema.validate(cfg)


def test_adaptive_full_block_validates(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    result = schema.validate(transient_cfg)
    adp = result["solver"]["adaptive"]
    assert adp["enabled"] is True
    assert adp["dt_min"] == pytest.approx(1.0e-12)
    assert adp["dt_max"] == pytest.approx(5.0e-9)


def test_adaptive_absent_validates(transient_cfg):
    """Configs without solver.adaptive validate (default fixed-dt path)."""
    result = schema.validate(transient_cfg)
    assert result["solver"].get("adaptive") is None


def test_adaptive_disabled_validates_without_dt_bounds(transient_cfg):
    """When enabled is false, dt_min/dt_max are not required."""
    transient_cfg["solver"]["adaptive"] = {"enabled": False}
    result = schema.validate(transient_cfg)
    assert result["solver"]["adaptive"]["enabled"] is False


def test_adaptive_enabled_without_dt_min_rejected(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_max": 5.0e-9,
    }
    with pytest.raises(schema.SchemaError, match="dt_min is required"):
        schema.validate(transient_cfg)


def test_adaptive_enabled_without_dt_max_rejected(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
    }
    with pytest.raises(schema.SchemaError, match="dt_max is required"):
        schema.validate(transient_cfg)


def test_adaptive_dt_min_exceeds_dt_max_rejected(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-8,
        "dt_max": 1.0e-10,
    }
    with pytest.raises(schema.SchemaError, match="must not exceed"):
        schema.validate(transient_cfg)


def test_adaptive_initial_dt_below_dt_min_rejected(transient_cfg):
    transient_cfg["solver"]["dt"] = 1.0e-13
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
    }
    with pytest.raises(schema.SchemaError, match="initial dt"):
        schema.validate(transient_cfg)


def test_adaptive_initial_dt_above_dt_max_rejected(transient_cfg):
    transient_cfg["solver"]["dt"] = 1.0e-7
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
    }
    with pytest.raises(schema.SchemaError, match="initial dt"):
        schema.validate(transient_cfg)


def test_adaptive_grow_factor_one_rejected(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
        "grow_factor": 1.0,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_adaptive_easy_iter_threshold_zero_rejected(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
        "easy_iter_threshold": 0,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_adaptive_max_consecutive_failures_zero_rejected(transient_cfg):
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
        "max_consecutive_failures": 0,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_adaptive_unknown_field_rejected(transient_cfg):
    """additionalProperties: false on the adaptive block."""
    transient_cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
        "this_field_does_not_exist": 1.0,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(transient_cfg)


def test_adaptive_on_bias_sweep_rejected(transient_cfg):
    """Setting solver.adaptive on bias_sweep is rejected at validate time."""
    cfg = copy.deepcopy(transient_cfg)
    cfg["solver"]["type"] = "bias_sweep"
    cfg["solver"].pop("t_end", None)
    cfg["solver"].pop("dt", None)
    cfg["contacts"][0]["voltage_sweep"] = {"start": 0.0, "stop": 0.5, "step": 0.1}
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
    }
    with pytest.raises(schema.SchemaError, match="only consumed by the transient runner"):
        schema.validate(cfg)


def test_adaptive_on_equilibrium_rejected(transient_cfg):
    cfg = copy.deepcopy(transient_cfg)
    cfg["solver"]["type"] = "equilibrium"
    cfg["solver"].pop("t_end", None)
    cfg["solver"].pop("dt", None)
    cfg["solver"]["adaptive"] = {
        "enabled": False,
    }
    with pytest.raises(schema.SchemaError, match="only consumed by the transient runner"):
        schema.validate(cfg)
