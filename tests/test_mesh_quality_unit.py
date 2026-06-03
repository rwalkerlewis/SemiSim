"""Pure-Python unit tests for `semi.mesh_quality`.

Exercises the dataclass shape, the threshold-merge logic, and the
geometric helpers (`_tet_edge_lengths`, `_tet_volumes`,
`_tet_skewness`, `_tet_aspect_ratio`) on hand-built tetrahedra.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.mesh_quality import (
    DEFAULT_THRESHOLDS,
    MeshQualityReport,
    PerGroupStats,
    _merge_thresholds,
    _tet_aspect_ratio,
    _tet_edge_lengths,
    _tet_skewness,
    _tet_volumes,
)


def _regular_tet(L: float = 1.0):
    """One regular tetrahedron with edge L."""
    return np.array([[
        [0.0, 0.0, 0.0],
        [L,   0.0, 0.0],
        [L * 0.5, L * np.sqrt(3) / 2.0, 0.0],
        [L * 0.5, L * np.sqrt(3) / 6.0, L * np.sqrt(2) / np.sqrt(3)],
    ]])


def _sliver_tet(t: float = 1e-3):
    """A degenerate sliver tetrahedron: one edge orders of magnitude
    shorter than the others, so aspect ratio L_max / L_min is large
    and skewness approaches 1."""
    return np.array([[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 1.0, t],   # last vertex sits just above the second one
    ]])


def test_regular_tet_metrics():
    coords = _regular_tet(1.0)
    edges = _tet_edge_lengths(coords)
    assert edges.shape == (1, 6)
    np.testing.assert_allclose(edges, 1.0, atol=1e-12)
    vol = _tet_volumes(coords)
    np.testing.assert_allclose(vol, np.sqrt(2.0) / 12.0, atol=1e-12)
    sk = _tet_skewness(coords)
    np.testing.assert_allclose(sk, 0.0, atol=1e-12)
    ar = _tet_aspect_ratio(coords)
    np.testing.assert_allclose(ar, 1.0, atol=1e-12)


def test_sliver_tet_metrics():
    coords = _sliver_tet(1e-3)
    sk = _tet_skewness(coords)
    ar = _tet_aspect_ratio(coords)
    assert sk[0] > 0.9, f"sliver skewness should approach 1, got {sk[0]}"
    assert ar[0] > 100.0, f"sliver aspect ratio should be large, got {ar[0]}"


def test_skewness_clipping_on_zero_volume():
    """Degenerate (coplanar) tet must give skewness exactly 1.0."""
    flat = np.array([[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ]])
    sk = _tet_skewness(flat)
    np.testing.assert_allclose(sk, 1.0, atol=1e-12)


def test_default_thresholds_shape():
    assert DEFAULT_THRESHOLDS["skewness_max"] > 0
    assert DEFAULT_THRESHOLDS["aspect_ratio_max"] > 1
    assert DEFAULT_THRESHOLDS["cell_count_min"] > 0
    assert (DEFAULT_THRESHOLDS["cell_count_max"]
            > DEFAULT_THRESHOLDS["cell_count_min"])
    assert 1 in DEFAULT_THRESHOLDS["per_group"]
    assert 4 in DEFAULT_THRESHOLDS["per_group"]


def test_merge_thresholds_overlays_per_group():
    merged = _merge_thresholds({
        "skewness_max": 0.5,
        "per_group": {4: {"edge_max": 2.0e-9}},
    })
    assert merged["skewness_max"] == 0.5
    # Defaults for other keys are preserved.
    assert (merged["aspect_ratio_max"]
            == DEFAULT_THRESHOLDS["aspect_ratio_max"])
    # The per-group override stacks on top of the default per_group;
    # tag 1 (si_body) is untouched, tag 4 (sio2_gate) gets the
    # override on edge_max.
    assert merged["per_group"][1]["edge_max"] == (
        DEFAULT_THRESHOLDS["per_group"][1]["edge_max"]
    )
    assert merged["per_group"][4]["edge_max"] == 2.0e-9
    assert merged["per_group"][4]["name"] == "sio2_gate"


def test_merge_thresholds_none_returns_copy():
    merged = _merge_thresholds(None)
    assert merged == DEFAULT_THRESHOLDS
    merged["skewness_max"] = 0.0
    assert DEFAULT_THRESHOLDS["skewness_max"] != 0.0


def test_per_group_stats_immutable():
    s = PerGroupStats(tag=1, name="x", cell_count=10,
                      edge_min=1e-9, edge_max=1e-8)
    with pytest.raises((AttributeError, Exception)):
        s.cell_count = 11


def test_mesh_quality_report_passed_property():
    r = MeshQualityReport(
        total_cells=10, skewness_max=0.5, aspect_ratio_max=5.0,
        per_group=(), thresholds={}, violations=(),
    )
    assert r.passed
    r2 = MeshQualityReport(
        total_cells=10, skewness_max=0.5, aspect_ratio_max=5.0,
        per_group=(), thresholds={}, violations=("oh no",),
    )
    assert not r2.passed
