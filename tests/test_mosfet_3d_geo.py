"""Pure-Python smoke tests for the M19 precursor mosfet_3d.msh artefact.

These checks run in the pure-Python CI environment (no dolfinx) and
verify that the shipped mesh file exists, has plausible size, and
that the source .geo file is present alongside it for regeneration.
The richer parse-and-inspect tests live in tests/fem/ where gmsh
is available.
"""
from __future__ import annotations

from pathlib import Path

import pytest

BENCH_DIR = Path(__file__).parent.parent / "benchmarks" / "mosfet_3d_eq" / "fixtures"


def test_geo_file_exists():
    geo = BENCH_DIR / "mosfet_3d.geo"
    assert geo.exists(), f"missing {geo}"
    text = geo.read_text()
    assert 'SetFactory("OpenCASCADE")' in text
    assert 'Mesh.ScalingFactor = 1e-6' in text, (
        ".geo must scale to meters on output (ADR 0019)"
    )
    for tag in ('si_body', 'sio2_gate', 'source', 'drain', 'gate', 'body'):
        assert f'"{tag}"' in text, f"physical group {tag} missing in .geo"


def test_msh_file_exists():
    msh = BENCH_DIR / "mosfet_3d.msh"
    assert msh.exists(), (
        f"missing {msh}; regenerate via "
        f"`docker compose run --rm dev gmsh -3 -format msh4 "
        f"benchmarks/mosfet_3d_eq/fixtures/mosfet_3d.geo "
        f"-o benchmarks/mosfet_3d_eq/fixtures/mosfet_3d.msh`"
    )
    size = msh.stat().st_size
    # 200k-500k tets at msh4 binary produces ~5-25 MB. Smaller
    # than 1 MB indicates an unrefined / failed mesh; larger than
    # 50 MB indicates the bulk Si max-edge was over-tightened.
    assert 1_000_000 < size < 50_000_000, (
        f"mosfet_3d.msh size {size} bytes is out of plausible range"
    )


def test_msh_binary_header():
    msh = BENCH_DIR / "mosfet_3d.msh"
    if not msh.exists():
        pytest.skip("mesh file not generated yet")
    with msh.open("rb") as f:
        head = f.read(128)
    # gmsh msh files start with "$MeshFormat\n".
    assert head.startswith(b"$MeshFormat\n"), (
        "mosfet_3d.msh does not start with the gmsh MeshFormat header"
    )
