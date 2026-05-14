"""
Tests for the M18.1 schema addition: solver.snes.line_search for the
bias_sweep coupled drift-diffusion block. Pure-Python; no dolfinx
required.

Mirrors the M18 schema-test pattern: every existing benchmark JSON
keeps validating, the new field is optional, and the loader rejects
out-of-enum values at validate time so users see the failure before
the FEM path.
"""
from __future__ import annotations

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
        "name": "snes_ls_test",
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
             "voltage": 0.0,
             "voltage_sweep": {"start": 0.0, "stop": 0.3, "step": 0.1}},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "solver": {
            "type": "bias_sweep",
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
def test_existing_benchmark_validates_against_v210(json_path):
    """Every shipped benchmark / example JSON validates unchanged on v2.10.0.

    Forward-compatibility gate: v2.0.0 through v2.9.0 inputs continue to
    parse without modification under the M18.1 schema bump.
    """
    with json_path.open() as f:
        cfg = json.load(f)
    if "schema_version" not in cfg:
        pytest.skip(f"{json_path} has no schema_version")
    schema.validate(cfg)


def test_line_search_absent_validates(bias_sweep_cfg):
    """Configs without solver.snes.line_search validate (default `bt`)."""
    result = schema.validate(bias_sweep_cfg)
    snes_block = result["solver"].get("snes") or {}
    # Default fill happens for snes.line_search only when the snes block
    # exists; otherwise the runner applies its own "bt" default. The
    # bias_sweep runner reads snes_opts.get("line_search", "bt").
    assert snes_block.get("line_search", "bt") == "bt"


def test_line_search_bt_validates(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["snes"] = {"line_search": "bt"}
    result = schema.validate(bias_sweep_cfg)
    assert result["solver"]["snes"]["line_search"] == "bt"


def test_line_search_nleqerr_validates(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["snes"] = {"line_search": "nleqerr"}
    result = schema.validate(bias_sweep_cfg)
    assert result["solver"]["snes"]["line_search"] == "nleqerr"


@pytest.mark.parametrize("ls", ["cp", "l2", "basic"])
def test_line_search_other_valid_types_validate(bias_sweep_cfg, ls):
    bias_sweep_cfg["solver"]["snes"] = {"line_search": ls}
    result = schema.validate(bias_sweep_cfg)
    assert result["solver"]["snes"]["line_search"] == ls


def test_line_search_unknown_rejected(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["snes"] = {"line_search": "not-a-real-type"}
    with pytest.raises(schema.SchemaError):
        schema.validate(bias_sweep_cfg)


def test_line_search_alongside_existing_snes_overrides(bias_sweep_cfg):
    """The line_search field coexists with the M12 rtol / atol overrides."""
    bias_sweep_cfg["solver"]["snes"] = {
        "rtol": 1.0e-12,
        "atol": 1.0e-9,
        "line_search": "nleqerr",
    }
    result = schema.validate(bias_sweep_cfg)
    snes_block = result["solver"]["snes"]
    assert snes_block["rtol"] == pytest.approx(1.0e-12)
    assert snes_block["atol"] == pytest.approx(1.0e-9)
    assert snes_block["line_search"] == "nleqerr"


def test_line_search_v29_input_still_validates_and_defaults_to_bt(
    bias_sweep_cfg,
):
    """A v2.9.0 input without line_search reads bt at runtime."""
    bias_sweep_cfg["schema_version"] = "2.9.0"
    result = schema.validate(bias_sweep_cfg)
    snes_block = result["solver"].get("snes") or {}
    assert snes_block.get("line_search", "bt") == "bt"
