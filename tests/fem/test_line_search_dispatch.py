"""
FEM integration tests for M18.1 line-search dispatch in bias_sweep.
"""
from __future__ import annotations

import copy

import pytest


def _tiny_bias_sweep_cfg() -> dict:
    return {
        "schema_version": "2.10.0",
        "name": "line_search_dispatch",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 2.0e-5]],
            "resolution": [80],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 2.0e-5]]},
            ],
            "facets_by_plane": [
                {"name": "anode", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 2.0e-5},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {
                    "type": "step",
                    "axis": 0,
                    "location": 1.0e-5,
                    "N_D_left": 0.0,
                    "N_A_left": 1.0e17,
                    "N_D_right": 1.0e17,
                    "N_A_right": 0.0,
                },
            }
        ],
        "contacts": [
            {
                "name": "anode",
                "facet": "anode",
                "type": "ohmic",
                "voltage": 0.0,
                "voltage_sweep": {"start": 0.0, "stop": 0.05, "step": 0.05},
            },
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
            "recombination": {"srh": True, "tau_n": 1.0e-8, "tau_p": 1.0e-8, "E_t": 0.0},
        },
        "solver": {
            "type": "bias_sweep",
            "continuation": {"min_step": 0.001, "max_halvings": 6, "max_step": 0.05},
        },
        "output": {"directory": "/tmp/line_search_dispatch", "fields": []},
    }


@pytest.mark.parametrize("line_search", ["bt", "nleqerr", "cp", "basic", "l2"])
def test_bias_sweep_line_search_reaches_petsc(line_search):
    from semi import schema
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = _tiny_bias_sweep_cfg()
    cfg["solver"]["line_search"] = line_search
    cfg = schema.validate(copy.deepcopy(cfg))
    result = run_bias_sweep(cfg)

    assert result.solver_info["converged"], result.solver_info
    ls_type = result.solver_info["problem"].solver.getLineSearch().getType()
    assert ls_type == line_search
