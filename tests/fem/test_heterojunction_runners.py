"""
M17 runner-level smoke tests with heterojunction configs.

These tests exercise the `if cfg_uses_heterojunction(...): het_fields =
build_dg0_material_fields(...)` branch in each runner that integrates
the M17 heterojunction path. The simulation uses a coarse 1D Si pn
diode with a `material_overrides: {chi_eV, ...}` flag on the n side
so the runner takes the heterojunction-fields code path even though
both regions are nominally silicon. The goal is integration coverage
of the runner-level branch, not a physics gate (the physics gates
live in `tests/fem/test_heterojunction_assembly.py` and the
`hemt_2d` benchmark).

Runs only inside the Dockerized FEM environment; tests/fem/conftest.py
skips this directory when dolfinx is unavailable.
"""
from __future__ import annotations


def _heterojunction_diode_cfg() -> dict:
    """Coarse 1D pn diode with a `material_overrides` chi shift on the
    n side to opt the runner into the heterojunction code path."""
    return {
        "schema_version": "2.8.0",
        "name": "het_runner_smoke",
        "description": "M17 runner smoke (Si pn diode with chi override).",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-5]],
            "resolution": [40],
            "regions_by_box": [
                {"name": "p_side", "tag": 1, "bounds": [[0.0, 5.0e-6]]},
                {"name": "n_side", "tag": 2, "bounds": [[5.0e-6, 1.0e-5]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-5},
            ],
        },
        "regions": {
            "p_side": {"material": "Si", "tag": 1, "role": "semiconductor"},
            "n_side": {
                "material": "Si",
                "tag": 2,
                "role": "semiconductor",
                # A small chi shift opts the cfg into the heterojunction
                # code path (cfg_uses_heterojunction returns True). The
                # shift is tiny so the solver still converges quickly.
                "material_overrides": {"chi_eV": 4.10},
            },
        },
        "doping": [
            {"region": "p_side", "profile": {"type": "uniform", "N_D": 0.0, "N_A": 1.0e17}},
            {"region": "n_side", "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "physics": {
            "temperature": 300.0,
            "statistics": "boltzmann",
            "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
        },
        "solver": {"type": "equilibrium"},
        "output": {"directory": "/tmp/het_runner_smoke", "fields": []},
    }


def test_run_equilibrium_heterojunction_path_executes():
    """`run_equilibrium` on a heterojunction cfg builds DG0 material
    fields and threads them through the equilibrium-Poisson form
    builder. The smoke gate is convergence + a populated psi field;
    the physics gate lives in the assembly-level tests."""
    from semi import schema
    from semi.runners.equilibrium import run_equilibrium

    cfg = schema.validate(_heterojunction_diode_cfg())
    result = run_equilibrium(cfg)
    assert result.solver_info.get("converged", False), (
        f"heterojunction equilibrium did not converge: {result.solver_info}"
    )
    assert result.psi.x.array.size > 0


def test_run_bias_sweep_heterojunction_path_executes():
    """`run_bias_sweep` on a heterojunction cfg takes the het_fields
    branch in the main DD setup. Use a single bias step (V=0 -> 0.05 V)
    to keep wallclock low."""
    from semi import schema
    from semi.runners.bias_sweep import run_bias_sweep

    cfg = _heterojunction_diode_cfg()
    cfg["solver"] = {
        "type": "bias_sweep",
        "snes": {"rtol": 1.0e-10, "atol": 1.0e-8, "stol": 1.0e-14, "max_it": 80},
    }
    cfg["contacts"][0]["voltage_sweep"] = {"start": 0.0, "stop": 0.05, "step": 0.05}
    cfg = schema.validate(cfg)
    result = run_bias_sweep(cfg)
    assert len(result.iv) >= 1
