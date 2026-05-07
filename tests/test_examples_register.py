"""
Coverage gate for the examples catalogue verifier registrations.

Each example shipped under `examples/` registers a smoke verifier in
`scripts/run_benchmark.py` via the `@register("name")` decorator. This
test file asserts each example's verifier name is present in the
`_VERIFIERS` registry so the registrations land in the gated suite
coverage report; without this gate, a typo in the registration would
go unnoticed until CI runs the example end-to-end (which only happens
under the docker-fem-benchmarks job, not the docker-fem-tests
coverage job).

The registrations themselves are pure-Python lookups and do not pull
in dolfinx, so this test file is safe to run in any environment that
can import the script module.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_BENCHMARK_PATH = REPO_ROOT / "scripts" / "run_benchmark.py"


def _load_run_benchmark_module():
    spec = importlib.util.spec_from_file_location(
        "run_benchmark_under_test", RUN_BENCHMARK_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def verifiers():
    return _load_run_benchmark_module()._VERIFIERS


@pytest.mark.parametrize(
    "example_name",
    [
        "nmos_idvgs",
        "schottky_iv_temperature",
        "power_diode_reverse_recovery",
    ],
)
def test_example_verifier_registered(verifiers, example_name):
    assert example_name in verifiers, (
        f"example {example_name!r} has no verifier registered in "
        f"scripts/run_benchmark.py; add @register({example_name!r}) "
        f"to the corresponding verify_<name> function"
    )
