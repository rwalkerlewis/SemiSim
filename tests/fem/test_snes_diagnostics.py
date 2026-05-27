"""
FEM integration tests for M18.1 bias_sweep SNES diagnostics.
"""
from __future__ import annotations

import copy
import json


def _tiny_bias_sweep_cfg() -> dict:
    return {
        "schema_version": "2.10.0",
        "name": "snes_diagnostics",
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
        "output": {"directory": "/tmp/snes_diagnostics", "fields": []},
    }


def test_bias_sweep_diagnostics_write_artifact_and_manifest(tmp_path):
    from semi import schema
    from semi.io.artifact import write_artifact
    from semi.io.reader import read_manifest
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = _tiny_bias_sweep_cfg()
    cfg["solver"]["diagnostics"] = True
    cfg = schema.validate(copy.deepcopy(cfg))
    result = run_bias_sweep(cfg)
    run_dir = write_artifact(result, tmp_path / "runs")
    manifest = read_manifest(run_dir)

    diagnostics_path = run_dir / "snes_diagnostics.json"
    assert diagnostics_path.exists()
    assert manifest["logs"] == [{"name": "snes_diagnostics", "path": "snes_diagnostics.json"}]

    payload = json.loads(diagnostics_path.read_text())
    assert payload["line_search_type"] == "bt"
    assert payload["solves"]
    first_solve = payload["solves"][0]
    assert "snes_converged_reason" in first_solve
    assert "snes_linesearch_reason" in first_solve
    assert first_solve["newton_iterations"]
    assert "residual_norm" in first_solve["newton_iterations"][0]
    assert "line_search_reason" in first_solve["newton_iterations"][0]


def test_bias_sweep_diagnostics_off_writes_no_extra_artifact(tmp_path):
    from semi import schema
    from semi.io.artifact import write_artifact
    from semi.io.reader import read_manifest
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = schema.validate(copy.deepcopy(_tiny_bias_sweep_cfg()))
    result = run_bias_sweep(cfg)
    run_dir = write_artifact(result, tmp_path / "runs")
    manifest = read_manifest(run_dir)

    assert not (run_dir / "snes_diagnostics.json").exists()
    assert "logs" not in manifest
