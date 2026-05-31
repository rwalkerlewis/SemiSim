"""End-to-end integration test for the M19 precursor benchmark.

Loads benchmarks/mosfet_3d_eq/mosfet_3d_eq.json, runs the
equilibrium solver, and asserts that every gate in
semi/verification/mosfet_3d_eq.py passes. This is also the CI
matrix gate (no allow-failure; ADR 0019).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from semi import run as run_mod
from semi import schema
from semi.verification.mosfet_3d_eq import verify_mosfet_3d_eq

BENCH = Path(__file__).parent.parent.parent / "benchmarks" / "mosfet_3d_eq"
CFG_PATH = BENCH / "mosfet_3d_eq.json"


@pytest.fixture(scope="module")
def equilibrium_result():
    if not (BENCH / "fixtures" / "mosfet_3d.msh").exists():
        pytest.skip("mosfet_3d.msh not present")
    cfg = schema.load(str(CFG_PATH))
    return run_mod.run(cfg)


def test_equilibrium_solve_converges(equilibrium_result):
    info = equilibrium_result.solver_info or {}
    assert info.get("converged", True), info


def test_all_verifier_gates_pass(equilibrium_result):
    checks = verify_mosfet_3d_eq(equilibrium_result)
    failures = [(label, detail) for label, ok, detail in checks if not ok]
    assert not failures, (
        "verify_mosfet_3d_eq failures:\n  - "
        + "\n  - ".join(f"{lab}: {det}" for lab, det in failures)
    )
    # ADR 0019 wants all five gates exercised.
    assert len(checks) == 5
