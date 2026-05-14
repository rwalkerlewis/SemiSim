"""
Tests for semi.solver.apply_snes_line_search (M18.1 / ADR 0018).

Pure-Python; no dolfinx required, but petsc4py is. The helper is
exercised against a freshly-created SNES so we can verify the line
search type round-trips through ``SNESLineSearch.setType`` /
``getType``.

The dolfinx 0.10 ``NonlinearProblem`` configures its SNES with
``SNESSetFromOptions`` before the line-search context exists, so
``snes_linesearch_*`` options pushed through ``petsc_options`` are
silently dropped. The helper exists specifically to set the line-search
type directly on the SNES object after construction.
"""
from __future__ import annotations

import pytest


def _make_snes():
    """Create a bare PETSc SNES for line-search-type round-tripping."""
    pytest.importorskip("petsc4py")
    from petsc4py import PETSc

    snes = PETSc.SNES().create()
    return snes


def test_apply_snes_line_search_bt():
    from semi.solver import apply_snes_line_search

    snes = _make_snes()
    apply_snes_line_search(snes, "bt")
    assert snes.getLineSearch().getType() == "bt"


def test_apply_snes_line_search_nleqerr():
    from semi.solver import apply_snes_line_search

    snes = _make_snes()
    apply_snes_line_search(snes, "nleqerr")
    assert snes.getLineSearch().getType() == "nleqerr"


@pytest.mark.parametrize("ls", ["cp", "l2", "basic"])
def test_apply_snes_line_search_other_types(ls):
    from semi.solver import apply_snes_line_search

    snes = _make_snes()
    apply_snes_line_search(snes, ls)
    assert snes.getLineSearch().getType() == ls


def test_apply_snes_line_search_noop_on_none():
    """Passing None is a no-op so callers can pass the cfg field
    unconditionally without an upstream branch."""
    from semi.solver import apply_snes_line_search

    snes = _make_snes()
    before = snes.getLineSearch().getType()
    apply_snes_line_search(snes, None)
    assert snes.getLineSearch().getType() == before


def test_apply_snes_line_search_noop_on_empty_string():
    from semi.solver import apply_snes_line_search

    snes = _make_snes()
    before = snes.getLineSearch().getType()
    apply_snes_line_search(snes, "")
    assert snes.getLineSearch().getType() == before
