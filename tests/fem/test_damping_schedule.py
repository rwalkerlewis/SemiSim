"""
FEM integration test for the M18.1 damping-schedule path in bias_sweep.
"""
from __future__ import annotations

import copy


def _tiny_bias_sweep_cfg() -> dict:
    return {
        "schema_version": "2.10.0",
        "name": "damping_schedule",
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
            "line_search": "nleqerr",
            "diagnostics": True,
            "damping_schedule": {
                "enabled": True,
                "lambda_start": 0.25,
                "lambda_end": 1.0,
                "decay_iters": 4,
            },
            "continuation": {"min_step": 0.001, "max_halvings": 6, "max_step": 0.05},
        },
        "output": {"directory": "/tmp/damping_schedule", "fields": []},
    }


def test_bias_sweep_damping_schedule_executes():
    from semi import schema
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = schema.validate(copy.deepcopy(_tiny_bias_sweep_cfg()))
    result = run_bias_sweep(cfg)

    assert result.solver_info["converged"], result.solver_info
    assert result.snes_diagnostics["damping_schedule"]["enabled"] is True
    assert result.snes_diagnostics["damping_schedule"]["lambda_start"] == 0.25
