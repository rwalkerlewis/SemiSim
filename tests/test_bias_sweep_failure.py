"""
Unit tests for the BiasSweepFailure exception payload (M18.1).

Pure-Python; no FEM/PETSc dependency. The bias_sweep runner attaches a
JSON-serialisable SNES diagnostics payload to the exception when the
coupled DD continuation gives up, so the run artifact can record the
failing solve. These tests pin that contract without driving a solve.
"""
from __future__ import annotations

from semi.runners.bias_sweep import BiasSweepFailure


def test_failure_carries_diagnostics_payload():
    diag = {"line_search_type": "nleqerr", "solves": []}
    err = BiasSweepFailure("did not converge", snes_diagnostics=diag)
    assert isinstance(err, RuntimeError)
    assert str(err) == "did not converge"
    assert err.snes_diagnostics is diag


def test_failure_diagnostics_default_none():
    err = BiasSweepFailure("did not converge")
    assert err.snes_diagnostics is None
