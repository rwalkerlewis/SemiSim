"""
Tests for the M18.1 schema addition `solver.line_search`.

Pure-Python; verifies enum widening, default-fill, non-bias_sweep
rejection, and forward-compat validation on v2.10.0 inputs.
"""
from __future__ import annotations

import copy

import pytest

from semi import schema


@pytest.fixture
def bias_sweep_cfg():
    return {
        "schema_version": "2.10.0",
        "name": "line_search_schema_test",
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


def test_line_search_default_fills_bt(bias_sweep_cfg):
    result = schema.validate(copy.deepcopy(bias_sweep_cfg))
    assert result["solver"]["line_search"] == "bt"


@pytest.mark.parametrize("line_search", ["bt", "nleqerr", "cp", "basic", "l2"])
def test_line_search_enum_values_validate(bias_sweep_cfg, line_search):
    bias_sweep_cfg["solver"]["line_search"] = line_search
    result = schema.validate(bias_sweep_cfg)
    assert result["solver"]["line_search"] == line_search


def test_line_search_unknown_value_rejected(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["line_search"] = "not-a-real-linesearch"
    with pytest.raises(schema.SchemaError):
        schema.validate(bias_sweep_cfg)


def test_line_search_on_non_bias_sweep_rejected(bias_sweep_cfg):
    cfg = copy.deepcopy(bias_sweep_cfg)
    cfg["solver"] = {"type": "transient", "line_search": "nleqerr", "t_end": 1.0e-9, "dt": 1.0e-10}
    cfg["contacts"][0].pop("voltage_sweep", None)
    with pytest.raises(schema.SchemaError, match="only consumed by the bias_sweep runner"):
        schema.validate(cfg)
