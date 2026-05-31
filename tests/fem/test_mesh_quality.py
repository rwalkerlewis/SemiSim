"""FEM-side integration tests for `semi.mesh_quality`.

Exercises `check_mesh_quality` on the committed mosfet_3d.msh,
plus the failure path triggered by deliberately-tight thresholds
(no separate broken .msh is committed; tightening the thresholds
exercises the same code path as a coarse-mesh failure).

The ingest-hook test confirms that `semi.mesh.build_mesh` runs
the gate when `mesh.quality_gate: true` and raises on violations.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from semi.mesh_quality import (
    MeshQualityError,
    check_mesh_quality,
)

BENCH = (Path(__file__).parent.parent.parent
         / "benchmarks" / "mosfet_3d_eq" / "fixtures")
MSH = BENCH / "mosfet_3d.msh"


@pytest.fixture(scope="module")
def mosfet_mesh():
    if not MSH.exists():
        pytest.skip(f"mesh not present: {MSH}")
    from dolfinx.io import gmsh as gmshio
    from mpi4py import MPI
    return gmshio.read_from_msh(str(MSH), MPI.COMM_WORLD, gdim=3)


def test_check_mesh_quality_passes_on_committed_mesh(mosfet_mesh):
    report = check_mesh_quality(mosfet_mesh.mesh, mosfet_mesh.cell_tags)
    assert report.passed, (
        f"committed mosfet_3d.msh fails default gate: "
        f"{report.violations}"
    )
    assert 200_000 < report.total_cells < 500_000
    assert report.skewness_max < 0.995
    assert report.aspect_ratio_max < 20.0
    names = {pg.name for pg in report.per_group}
    assert "si_body" in names and "sio2_gate" in names


def test_check_mesh_quality_per_group_stats(mosfet_mesh):
    report = check_mesh_quality(mosfet_mesh.mesh, mosfet_mesh.cell_tags)
    by_name = {pg.name: pg for pg in report.per_group}
    assert by_name["sio2_gate"].cell_count > 10_000
    assert by_name["si_body"].cell_count > 10_000
    # Min edge in the oxide is finer than min edge in the Si bulk
    # because the oxide carries the tightest size cap.
    assert by_name["sio2_gate"].edge_max < 7.0e-9
    assert by_name["si_body"].edge_max < 3.5e-7


def test_check_mesh_quality_fails_tight_oxide_threshold(mosfet_mesh):
    """Tighten the oxide max-edge below the as-shipped value;
    the same code path that catches a broken mesh fires here."""
    tight = {"per_group": {4: {"edge_max": 2.0e-9}}}
    report = check_mesh_quality(
        mosfet_mesh.mesh, mosfet_mesh.cell_tags, thresholds=tight,
    )
    assert not report.passed
    assert any("sio2_gate" in v for v in report.violations)


def test_check_mesh_quality_fails_tight_skewness(mosfet_mesh):
    """Tighten skewness below the worst tet's skewness."""
    tight = {"skewness_max": 0.0}
    report = check_mesh_quality(
        mosfet_mesh.mesh, mosfet_mesh.cell_tags, thresholds=tight,
    )
    assert not report.passed
    assert any("skewness" in v for v in report.violations)


def test_check_mesh_quality_fails_cell_count_band(mosfet_mesh):
    """Cell count outside the [min, max] band must fire."""
    report = check_mesh_quality(
        mosfet_mesh.mesh, mosfet_mesh.cell_tags,
        thresholds={"cell_count_min": 1_000_000, "cell_count_max": 2_000_000},
    )
    assert not report.passed
    assert any("below threshold" in v for v in report.violations)


def test_build_mesh_quality_gate_passes(tmp_path):
    """End-to-end: build_mesh with quality_gate=True returns
    successfully on the committed mesh."""
    from semi.mesh import build_mesh

    cfg = {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": str(MSH),
            "format": "gmsh",
            "quality_gate": True,
        },
        "_source_dir": str(MSH.parent),
    }
    msh, ct, ft = build_mesh(cfg)
    assert msh.topology.dim == 3
    assert ct is not None


def test_build_mesh_quality_gate_raises_on_tight_thresholds():
    """End-to-end: build_mesh with quality_gate=True and tight
    thresholds raises MeshQualityError."""
    from semi.mesh import build_mesh

    cfg = {
        "dimension": 3,
        "mesh": {
            "source": "file",
            "path": str(MSH),
            "format": "gmsh",
            "quality_gate": True,
            "quality_thresholds": {
                "per_group": {4: {"edge_max": 1.0e-12}},
            },
        },
        "_source_dir": str(MSH.parent),
    }
    with pytest.raises(MeshQualityError, match="sio2_gate"):
        build_mesh(cfg)
