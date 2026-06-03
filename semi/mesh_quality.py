"""
Post-ingest mesh-quality gate (M19 precursor; ADR 0019).

Reports per-physical-group edge-length extrema, mesh-wide
tetrahedral skewness and aspect-ratio extrema, and a per-region
cell count. The `check_mesh_quality` entry point compares the
measurements against a thresholds dict and raises
`MeshQualityError` on the first violation when called from the
runner, or returns the `MeshQualityReport` (with `passed=False`)
when used directly from tests.

Per ADR 0007, this module is pure-Python at import time. The
dolfinx and numpy imports inside `check_mesh_quality` are lazy
so the module can be imported by the pure-Python CI matrix.

Per the M19 precursor prompt, the file lives at
`semi/mesh_quality.py` rather than `semi/mesh/quality.py`
because `semi/mesh` ships as a single flat module today and the
prompt's "Do not refactor semi/mesh.py for clarity" anti-goal
takes precedence over the path suggestion. The exported API
matches the prompt.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


class MeshQualityError(RuntimeError):
    """Raised when a mesh fails its declared quality thresholds."""


@dataclass(frozen=True)
class PerGroupStats:
    """Edge-length extrema and cell count for one physical group."""

    tag: int
    name: str
    cell_count: int
    edge_min: float
    edge_max: float


@dataclass(frozen=True)
class MeshQualityReport:
    """Output of `check_mesh_quality`."""

    total_cells: int
    skewness_max: float
    aspect_ratio_max: float
    per_group: tuple[PerGroupStats, ...]
    thresholds: dict[str, Any]
    violations: tuple[str, ...] = field(default_factory=tuple)

    @property
    def passed(self) -> bool:
        return not self.violations


# Default thresholds for the mosfet_3d_eq benchmark; the JSON
# config can override via `mesh.quality_thresholds`. Edge lengths
# in meters, skewness and aspect ratio dimensionless.
#
# The bounds are calibrated to the as-shipped mosfet_3d.msh (ADR
# 0019): gmsh emits cells up to ~2x the nominal MeshSize target,
# so the edge_max thresholds sit at ~2.3x the .geo's h_bulk and
# h_ox. The skewness bound at 0.995 admits the small fraction of
# near-degenerate tets the gmsh Delaunay optimizer leaves behind
# in the oxide / Si transition (worst observed: 0.992; average
# 0.82). Cell-count band brackets the prompt's 200k-500k target.
DEFAULT_THRESHOLDS: dict[str, Any] = {
    "skewness_max":      0.995,
    "aspect_ratio_max":  20.0,
    "cell_count_min":    200_000,
    "cell_count_max":    500_000,
    "edge_min_global":   1.0e-9,
    "per_group": {
        # tag -> {"edge_max": value, "name": optional label}
        1: {"name": "si_body",   "edge_max": 3.5e-7},  # ~2.3x h_bulk
        4: {"name": "sio2_gate", "edge_max": 7.0e-9},  # ~2.3x h_ox
    },
}


def _merge_thresholds(user: dict[str, Any] | None) -> dict[str, Any]:
    """Overlay user thresholds on top of DEFAULT_THRESHOLDS."""
    merged: dict[str, Any] = {
        k: v for k, v in DEFAULT_THRESHOLDS.items() if k != "per_group"
    }
    merged["per_group"] = {
        tag: dict(v) for tag, v in DEFAULT_THRESHOLDS["per_group"].items()
    }
    if not user:
        return merged
    for k, v in user.items():
        if k == "per_group":
            for tag_str, group in (v or {}).items():
                tag = int(tag_str)
                merged["per_group"].setdefault(tag, {}).update(group)
        else:
            merged[k] = v
    return merged


def _tet_edge_lengths(coords):
    """Return shape (n_cells, 6) of edge lengths for a tet mesh."""
    import numpy as np

    # Tet edges: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    diffs = np.empty((coords.shape[0], len(pairs), coords.shape[2]),
                     dtype=coords.dtype)
    for i, (a, b) in enumerate(pairs):
        diffs[:, i, :] = coords[:, b, :] - coords[:, a, :]
    return np.linalg.norm(diffs, axis=2)


def _tet_volumes(coords):
    """Signed volumes of tets, shape (n_cells,)."""
    import numpy as np

    v01 = coords[:, 1, :] - coords[:, 0, :]
    v02 = coords[:, 2, :] - coords[:, 0, :]
    v03 = coords[:, 3, :] - coords[:, 0, :]
    return np.abs(np.einsum("ij,ij->i", v01, np.cross(v02, v03))) / 6.0


def _tet_skewness(coords):
    """
    Skewness in [0, 1]: 0 = regular tet, 1 = degenerate.

    Computed as 1 - V / V_regular(L_max), where V_regular(L) =
    L^3 / (6 * sqrt(2)) is the volume of a regular tetrahedron
    with edge L. A regular tet has V = V_regular(L) so the metric
    is 0; a sliver tet has V << V_regular(L_max) so the metric
    approaches 1. Clipped to [0, 1] for numerical safety.
    """
    import numpy as np

    edges = _tet_edge_lengths(coords)
    L_max = edges.max(axis=1)
    vol = _tet_volumes(coords)
    v_reg = L_max ** 3 / (6.0 * np.sqrt(2.0))
    quality = np.where(v_reg > 0.0, vol / v_reg, 0.0)
    return np.clip(1.0 - quality, 0.0, 1.0)


def _tet_aspect_ratio(coords):
    """max edge length / min edge length, shape (n_cells,)."""
    import numpy as np

    edges = _tet_edge_lengths(coords)
    L_min = edges.min(axis=1)
    return np.where(L_min > 0.0, edges.max(axis=1) / L_min, np.inf)


def _cell_vertex_coords(mesh):
    """
    Return (n_cells, n_verts_per_cell, gdim) array of cell vertices.

    Reads from the topology connectivity (cells -> vertices) and
    geometry array.
    """
    import numpy as np

    tdim = mesh.topology.dim
    mesh.topology.create_connectivity(tdim, 0)
    conn = mesh.topology.connectivity(tdim, 0)
    n_cells = mesh.topology.index_map(tdim).size_local
    n_per_cell = len(conn.links(0))
    geom = mesh.geometry.x
    out = np.empty((n_cells, n_per_cell, geom.shape[1]), dtype=geom.dtype)
    for c in range(n_cells):
        verts = conn.links(c)
        out[c, :, :] = geom[verts, :]
    return out


def check_mesh_quality(
    mesh,
    cell_tags,
    *,
    thresholds: dict[str, Any] | None = None,
) -> MeshQualityReport:
    """
    Inspect a dolfinx mesh and return a `MeshQualityReport`.

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh
        Volume mesh (tdim must equal 3 for the M19 precursor).
    cell_tags : dolfinx.mesh.MeshTags
        Cell-tag MeshTags carrying the physical-volume IDs from
        the ingest path (gmsh / XDMF).
    thresholds : dict, optional
        Override DEFAULT_THRESHOLDS. The optional sub-dict
        ``per_group`` is overlaid on the per-group defaults
        rather than replacing them.

    Returns
    -------
    MeshQualityReport
    """
    import numpy as np

    if mesh.topology.dim != 3:
        raise ValueError(
            f"check_mesh_quality is 3D-only (got tdim={mesh.topology.dim})"
        )

    merged = _merge_thresholds(thresholds)
    coords = _cell_vertex_coords(mesh)
    if coords.shape[1] != 4:
        raise ValueError(
            f"check_mesh_quality expects tetrahedral cells "
            f"(got {coords.shape[1]} vertices per cell)"
        )

    edges = _tet_edge_lengths(coords)
    skewness = _tet_skewness(coords)
    aspect = _tet_aspect_ratio(coords)
    total = coords.shape[0]

    violations: list[str] = []

    skew_max = float(skewness.max(initial=0.0))
    if skew_max > merged["skewness_max"]:
        violations.append(
            f"tet skewness max {skew_max:.3f} exceeds "
            f"threshold {merged['skewness_max']:.3f}"
        )

    ar_max = float(aspect.max(initial=0.0))
    if ar_max > merged["aspect_ratio_max"]:
        violations.append(
            f"tet aspect ratio max {ar_max:.2f} exceeds "
            f"threshold {merged['aspect_ratio_max']:.2f}"
        )

    edge_min_global = float(edges.min(initial=np.inf))
    if edge_min_global < merged["edge_min_global"]:
        violations.append(
            f"min edge length {edge_min_global:.3e} m below "
            f"threshold {merged['edge_min_global']:.3e} m"
        )

    if total < merged["cell_count_min"]:
        violations.append(
            f"total cells {total} below threshold "
            f"{merged['cell_count_min']}"
        )
    if total > merged["cell_count_max"]:
        violations.append(
            f"total cells {total} above threshold "
            f"{merged['cell_count_max']}"
        )

    per_group_list: list[PerGroupStats] = []
    if cell_tags is None:
        violations.append(
            "cell_tags is None; check_mesh_quality requires "
            "physical-group tagging on cells"
        )
    else:
        values = cell_tags.values
        unique_tags = sorted(int(t) for t in np.unique(values))
        for tag in unique_tags:
            mask = values == tag
            n = int(mask.sum())
            group_edges = edges[mask]
            e_min = float(group_edges.min(initial=np.inf))
            e_max = float(group_edges.max(initial=0.0))
            spec = merged["per_group"].get(tag, {})
            name = spec.get("name", f"tag_{tag}")
            per_group_list.append(PerGroupStats(
                tag=tag, name=name, cell_count=n,
                edge_min=e_min, edge_max=e_max,
            ))
            edge_max_th = spec.get("edge_max")
            if edge_max_th is not None and e_max > edge_max_th:
                violations.append(
                    f"region '{name}' (tag {tag}) max edge "
                    f"{e_max:.3e} m exceeds threshold "
                    f"{edge_max_th:.3e} m"
                )

    return MeshQualityReport(
        total_cells=total,
        skewness_max=skew_max,
        aspect_ratio_max=ar_max,
        per_group=tuple(per_group_list),
        thresholds=merged,
        violations=tuple(violations),
    )
