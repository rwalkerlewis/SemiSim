"""dolfinx-side validation of the M19 precursor mosfet_3d.msh.

Confirms the .msh ingests cleanly via dolfinx.io.gmsh, has the
expected physical-group counts, and produces a mesh inside the
expected geometric bounding box (meters).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI

BENCH_DIR = Path(__file__).parent.parent.parent / "benchmarks" / "mosfet_3d_eq"
MSH = BENCH_DIR / "fixtures" / "mosfet_3d.msh"


@pytest.fixture(scope="module")
def mesh_data():
    if not MSH.exists():
        pytest.skip(f"mosfet_3d.msh not present at {MSH}")
    return gmshio.read_from_msh(str(MSH), MPI.COMM_WORLD, gdim=3)


def test_physical_group_names_and_tags(mesh_data):
    pg = mesh_data[5]
    assert pg["si_body"].tag == 1 and pg["si_body"].dim == 3
    assert pg["sio2_gate"].tag == 4 and pg["sio2_gate"].dim == 3
    assert pg["source"].tag == 10 and pg["source"].dim == 2
    assert pg["drain"].tag == 11 and pg["drain"].dim == 2
    assert pg["gate"].tag == 12 and pg["gate"].dim == 2
    assert pg["body"].tag == 13 and pg["body"].dim == 2


def test_mesh_bbox_in_meters(mesh_data):
    mesh = mesh_data.mesh
    pts = mesh.geometry.x
    assert pts[:, 0].min() == pytest.approx(0.0, abs=1e-15)
    assert pts[:, 0].max() == pytest.approx(1.5e-6, rel=1e-9)
    assert pts[:, 1].min() == pytest.approx(0.0, abs=1e-15)
    assert pts[:, 1].max() == pytest.approx(5.0e-7, rel=1e-9)
    assert pts[:, 2].min() == pytest.approx(0.0, abs=1e-15)
    assert pts[:, 2].max() == pytest.approx(2.005e-6, rel=1e-9)


def test_cell_counts_per_region(mesh_data):
    mesh = mesh_data.mesh
    ct = mesh_data.cell_tags
    total = mesh.topology.index_map(mesh.topology.dim).size_local
    # Cell budget set by ADR 0019: 200k to 500k total tets.
    assert 200_000 < total < 500_000, f"unexpected total cells: {total}"
    assert ct is not None
    si = (ct.values == 1).sum()
    ox = (ct.values == 4).sum()
    assert si + ox == total, "every cell must carry a physical-volume tag"
    # Both regions must be non-empty; the oxide carries the bulk
    # of the cells because of the 5 nm thickness and the 3 nm
    # max-edge cap (ADR 0019).
    assert si > 10_000
    assert ox > 10_000


def test_facet_counts_per_contact(mesh_data):
    ft = mesh_data.facet_tags
    assert ft is not None
    counts = {int(tag): int((ft.values == tag).sum())
              for tag in np.unique(ft.values)}
    for tag in (10, 11, 12, 13):
        assert counts.get(tag, 0) > 0, f"contact tag {tag} has no facets"
    # The gate contact (top of oxide) has the most facets because
    # the oxide region carries the finest mesh; body bottom has
    # the fewest because it sits at h_bulk far from any anchor.
    assert counts[12] > counts[10] > counts[13]
    assert counts[12] > counts[11] > counts[13]
