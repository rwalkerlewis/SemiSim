"""Phase C tests: schema v2.11.0 mesh.quality_gate / quality_thresholds.

Covers:
- default-fill of mesh.quality_gate to false
- override-threshold parsing and per-group nesting
- strict-mode rejection of unknown keys under quality_thresholds
- every existing benchmark JSON validates as v2.11.0 unchanged
  (existing-benchmark byte-identity at the schema layer)
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from semi.schema import SCHEMA_SUPPORTED_MINOR, SchemaError, validate

BENCHMARKS_DIR = Path(__file__).parent.parent / "benchmarks"
EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


def _minimal_v2(extra_mesh: dict | None = None) -> dict:
    """Smallest config that validates against v2."""
    return {
        "name": "schema_quality_gate_smoke",
        "schema_version": "2.11.0",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1e-6]],
            "resolution": [10],
            **(extra_mesh or {}),
        },
        "regions": {"si": {"material": "Si"}},
        "doping": [
            {"region": "si",
             "profile": {"type": "uniform", "N_D": 1e16, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "left",  "facet": "left",  "type": "ohmic",
             "voltage": 0.0},
            {"name": "right", "facet": "right", "type": "ohmic",
             "voltage": 0.0},
        ],
        "physics": {},
        "solver": {"type": "equilibrium"},
    }


def test_schema_minor_advanced_to_11():
    assert SCHEMA_SUPPORTED_MINOR >= 11


def test_quality_gate_defaults_to_false():
    cfg = _minimal_v2()
    out = validate(cfg)
    assert out["mesh"]["quality_gate"] is False


def test_quality_gate_explicit_true_validates():
    cfg = _minimal_v2({"quality_gate": True})
    out = validate(cfg)
    assert out["mesh"]["quality_gate"] is True


def test_quality_thresholds_override_parses():
    cfg = _minimal_v2({
        "quality_gate": True,
        "quality_thresholds": {
            "skewness_max": 0.9,
            "aspect_ratio_max": 15.0,
            "cell_count_min": 100,
            "cell_count_max": 1000,
            "edge_min_global": 5e-10,
            "per_group": {
                "1": {"name": "si", "edge_max": 2.0e-7},
                "4": {"edge_max": 4.0e-9}
            }
        }
    })
    out = validate(cfg)
    qt = out["mesh"]["quality_thresholds"]
    assert qt["skewness_max"] == 0.9
    assert qt["per_group"]["1"]["edge_max"] == 2.0e-7
    assert qt["per_group"]["4"]["edge_max"] == 4.0e-9


def test_quality_thresholds_unknown_key_rejected():
    cfg = _minimal_v2({
        "quality_gate": True,
        "quality_thresholds": {"unknown_key": 42},
    })
    with pytest.raises(SchemaError):
        validate(cfg)


def test_quality_thresholds_unknown_per_group_subkey_rejected():
    cfg = _minimal_v2({
        "quality_gate": True,
        "quality_thresholds": {
            "per_group": {"1": {"bogus": 1}}
        },
    })
    with pytest.raises(SchemaError):
        validate(cfg)


def test_quality_gate_field_in_file_mesh_branch():
    """quality_gate is available on the file-source mesh branch too."""
    cfg = {
        "name": "file_mesh_quality_gate",
        "schema_version": "2.11.0",
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": "ignored.msh",
            "format": "gmsh",
            "quality_gate": True,
        },
        "regions": {"si": {"material": "Si"}},
        "doping": [
            {"region": "si",
             "profile": {"type": "uniform", "N_D": 1e16, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "a", "facet": "source", "type": "ohmic", "voltage": 0.0},
            {"name": "b", "facet": "drain",  "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {},
        "solver": {"type": "equilibrium"},
    }
    out = validate(cfg)
    assert out["mesh"]["quality_gate"] is True


def _all_benchmark_configs() -> list[Path]:
    paths: list[Path] = []
    for root in (BENCHMARKS_DIR, EXAMPLES_DIR):
        if not root.is_dir():
            continue
        for p in root.rglob("*.json"):
            if "/fixtures/" in str(p):
                continue
            paths.append(p)
    return sorted(paths)


@pytest.mark.parametrize("cfg_path", _all_benchmark_configs(),
                         ids=lambda p: str(p.relative_to(p.parent.parent)))
def test_existing_benchmark_validates_unchanged(cfg_path):
    """Every existing benchmark and example JSON validates as v2.11.0
    without source edits (additive bump byte-identity)."""
    with cfg_path.open() as f:
        cfg = json.load(f)
    # Skip configs that target schema major 1 (deprecated; one
    # minor cycle window per ADR M14.3).
    schema_version = cfg.get("schema_version", "2.0.0")
    if schema_version.startswith("1."):
        pytest.skip(f"v1 schema benchmark: {cfg_path}")
    # Deep-copy because validate fills defaults in place.
    original = copy.deepcopy(cfg)
    out = validate(cfg)
    # quality_gate either absent (then default-filled to false) or
    # explicitly set (e.g. mosfet_3d_eq sets true).
    assert out["mesh"]["quality_gate"] in (True, False)
    # Existing configs that did not set quality_gate must end up
    # at the default (False) so byte-identity holds on the no-gate
    # branch.
    if "quality_gate" not in original.get("mesh", {}):
        assert out["mesh"]["quality_gate"] is False
