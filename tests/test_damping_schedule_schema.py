"""
Tests for the M18.1 schema addition `solver.damping_schedule`.
"""
from __future__ import annotations

import copy

import pytest

from semi import schema


@pytest.fixture
def bias_sweep_cfg():
    return {
        "schema_version": "2.10.0",
        "name": "damping_schedule_schema_test",
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


def test_damping_schedule_disabled_validates_without_parameters(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["damping_schedule"] = {"enabled": False}
    result = schema.validate(copy.deepcopy(bias_sweep_cfg))
    assert result["solver"]["damping_schedule"]["enabled"] is False


def test_damping_schedule_enabled_validates(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["damping_schedule"] = {
        "enabled": True,
        "lambda_start": 0.25,
        "lambda_end": 1.0,
        "decay_iters": 6,
    }
    result = schema.validate(copy.deepcopy(bias_sweep_cfg))
    assert result["solver"]["damping_schedule"]["lambda_start"] == pytest.approx(0.25)


def test_damping_schedule_missing_decay_iters_default_fills(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["damping_schedule"] = {
        "enabled": True,
        "lambda_start": 0.25,
        "lambda_end": 1.0,
    }
    result = schema.validate(copy.deepcopy(bias_sweep_cfg))
    assert result["solver"]["damping_schedule"]["decay_iters"] == 1


def test_damping_schedule_rejects_lambda_out_of_range(bias_sweep_cfg):
    bias_sweep_cfg["solver"]["damping_schedule"] = {
        "enabled": True,
        "lambda_start": 0.0,
        "lambda_end": 1.0,
        "decay_iters": 6,
    }
    with pytest.raises(schema.SchemaError, match="less than or equal to the minimum"):
        schema.validate(bias_sweep_cfg)


def test_damping_schedule_on_non_bias_sweep_rejected(bias_sweep_cfg):
    cfg = copy.deepcopy(bias_sweep_cfg)
    cfg["solver"] = {
        "type": "equilibrium",
        "damping_schedule": {"enabled": False},
    }
    cfg["contacts"][0].pop("voltage_sweep", None)
    with pytest.raises(schema.SchemaError, match="only consumed by the bias_sweep runner"):
        schema.validate(cfg)
