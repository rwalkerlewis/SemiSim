"""GET /runs, GET /runs/{id}, plus field / iv / log / input / manifest downloads."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse

from ..jobs import read_status
from ..models import RunListResponse, RunSummary

router = APIRouter()


def _created_at_from_run_id(run_id: str) -> str:
    """Parse the ISO-8601 timestamp prefix from `<ts>_<name>_<sha>` IDs.

    The run_id encodes the timestamp as `YYYY-MM-DDTHH-MM-SSZ` (hyphens in
    the time portion so the string is filename-safe). Convert back to the
    standard ISO-8601 `YYYY-MM-DDTHH:MM:SSZ`.
    """
    stamp = run_id.split("_", 1)[0]
    if "T" not in stamp:
        return run_id
    date_part, time_part = stamp.split("T", 1)
    return f"{date_part}T{time_part.replace('-', ':')}"


def _status_for(run_dir: Path) -> str:
    manifest_path = run_dir / "manifest.json"
    if manifest_path.exists():
        try:
            m = json.loads(manifest_path.read_text())
            status = m.get("status")
            if status in ("completed", "failed", "running"):
                return status
        except Exception:
            pass
    status_rec = read_status(run_dir)
    if status_rec and isinstance(status_rec.get("status"), str):
        return status_rec["status"]
    return "unknown"


def _run_summary_payload(run_id: str, run_dir: Path) -> dict[str, Any]:
    return {
        "run_id": run_id,
        "status": _status_for(run_dir),
        "created_at": _created_at_from_run_id(run_id),
    }


@router.get("/runs", response_model=RunListResponse)
def list_runs(request: Request) -> RunListResponse:
    storage = request.app.state.storage
    entries: list[RunSummary] = []
    for rid in storage.list_run_ids():
        rd = storage.run_dir(rid)
        entries.append(RunSummary(**_run_summary_payload(rid, rd)))
    return RunListResponse(runs=entries)


def _resolve_run_dir(request: Request, run_id: str) -> Path:
    storage = request.app.state.storage
    run_dir = storage.run_dir(run_id)
    if not run_dir.exists() or not run_dir.is_dir():
        raise HTTPException(status_code=404, detail=f"Run {run_id!r} not found.")
    return run_dir


@router.get("/runs/{run_id}")
def get_run(run_id: str, request: Request) -> JSONResponse:
    run_dir = _resolve_run_dir(request, run_id)
    manifest_path = run_dir / "manifest.json"
    status = _status_for(run_dir)
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
        except Exception as exc:  # pragma: no cover
            raise HTTPException(
                status_code=500, detail=f"Corrupt manifest: {exc}",
            ) from exc
        manifest["status"] = status
        return JSONResponse(content=manifest)
    payload = _run_summary_payload(run_id, run_dir)
    status_rec = read_status(run_dir) or {}
    for key in ("input_sha256", "wall_time_s", "error"):
        if key in status_rec:
            payload[key] = status_rec[key]
    return JSONResponse(content=payload)


@router.get("/runs/{run_id}/manifest")
def get_manifest(run_id: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "manifest.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Manifest not yet written.")
    return FileResponse(str(path), media_type="application/json")


@router.get("/runs/{run_id}/input")
def get_input(run_id: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "input.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Input JSON not found.")
    return FileResponse(str(path), media_type="application/json")


@router.get("/runs/{run_id}/fields/{name}")
def get_field(run_id: str, name: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    fields_dir = run_dir / "fields"
    bp_dir = fields_dir / f"{name}.bp"
    xdmf_path = fields_dir / f"{name}.xdmf"
    if bp_dir.exists() and bp_dir.is_dir():
        # VTX BP directories are multi-file; stream a tar.
        import tarfile
        import tempfile
        fd, tmp = tempfile.mkstemp(prefix=f"{name}_", suffix=".tar")
        import os
        os.close(fd)
        with tarfile.open(tmp, "w") as tar:
            tar.add(str(bp_dir), arcname=f"{name}.bp")
        return FileResponse(
            tmp,
            media_type="application/x-tar",
            filename=f"{name}.bp.tar",
        )
    if xdmf_path.exists():
        return FileResponse(str(xdmf_path), media_type="application/xml")
    raise HTTPException(status_code=404, detail=f"Field {name!r} not found.")


@router.get("/runs/{run_id}/iv/{contact}")
def get_iv(run_id: str, contact: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "iv" / f"{contact}.csv"
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"IV file for {contact!r} not found.")
    return FileResponse(str(path), media_type="text/csv")


@router.get("/runs/{run_id}/logs")
def get_logs(run_id: str, request: Request) -> FileResponse:
    run_dir = _resolve_run_dir(request, run_id)
    path = run_dir / "logs" / "engine.log"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Engine log not found.")
    return FileResponse(str(path), media_type="text/plain")


@router.get("/runs/{run_id}/fields/{name}/profile")
def get_field_profile(run_id: str, name: str, request: Request):
    """Return a 2-column CSV (x_m,value) for a 1-D scalar field stored in
    a VTX/ADIOS2 .bp directory.  Raises 404 if the field does not exist and
    422 if the mesh is not 1-D (profile extraction only makes sense for 1-D).
    """
    run_dir = _resolve_run_dir(request, run_id)
    bp_dir = run_dir / "fields" / f"{name}.bp"
    if not bp_dir.exists():
        raise HTTPException(status_code=404, detail=f"Field {name!r} not found.")

    try:
        import adios2  # type: ignore
        import numpy as np
        from fastapi.responses import Response

        with adios2.FileReader(str(bp_dir)) as f:
            vars_ = f.available_variables()
            if "geometry" not in vars_ or name not in vars_:
                raise HTTPException(
                    status_code=422,
                    detail=f"Field {name!r} or geometry not found in BP file.",
                )
            geom = f.read("geometry")   # (N, 3) -- VTX always writes xyz
            vals = f.read(name).flatten()

        if geom.ndim != 2 or geom.shape[1] != 3:
            raise HTTPException(
                status_code=422,
                detail="Profile export only supported for 1-D meshes.",
            )

        # Decide which axis has non-trivial spread.
        spans = geom.max(axis=0) - geom.min(axis=0)
        axis = int(np.argmax(spans))
        x = geom[:, axis]
        order = np.argsort(x)
        x = x[order]
        vals = vals[order]

        lines = [f"x_m,{name}"]
        lines += [f"{xi:.8e},{vi:.8e}" for xi, vi in zip(x, vals)]
        return Response("\n".join(lines), media_type="text/csv")

    except HTTPException:
        raise
    except ImportError:
        raise HTTPException(
            status_code=503,
            detail="adios2 is not available on this server.",
        )
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Could not read field: {exc}",
        ) from exc


@router.get("/runs/{run_id}/mesh")
def get_mesh(run_id: str, request: Request, max_nodes: int = 4000) -> JSONResponse:
    """Return mesh geometry and topology as JSON for visualization.

    Returns ``{dim, nodes: [[x,y,...]], cells: [[i,j,k],...], downsampled}``.
    For large meshes, nodes and cells are downsampled to ``max_nodes`` nodes.
    """
    run_dir = _resolve_run_dir(request, run_id)
    mesh_dir = run_dir / "mesh"
    xdmf_path = mesh_dir / "mesh.xdmf"
    if not xdmf_path.exists():
        raise HTTPException(status_code=404, detail="Mesh file not found for this run.")

    try:
        import mpi4py.MPI as MPI  # type: ignore
        import numpy as np
        from dolfinx.io import XDMFFile  # type: ignore

        with XDMFFile(MPI.COMM_WORLD, str(xdmf_path), "r") as fh:
            msh = fh.read_mesh()

        gdim = msh.geometry.dim
        tdim = msh.topology.dim
        coords = msh.geometry.x  # (N, 3)

        msh.topology.create_connectivity(tdim, 0)
        conn = msh.topology.connectivity(tdim, 0)
        num_cells = msh.topology.index_map(tdim).size_local

        cells_list = [conn.links(i).tolist() for i in range(num_cells)]

        # Downsample large meshes for payload size
        downsampled = False
        if coords.shape[0] > max_nodes:
            downsampled = True
            keep = np.random.default_rng(0).choice(
                coords.shape[0], size=max_nodes, replace=False
            )
            keep_set = set(keep.tolist())
            coords = coords[keep]
            # Remap cell vertices; drop cells with missing vertices
            remap = {old: new for new, old in enumerate(sorted(keep_set))}
            new_cells = []
            for c in cells_list:
                if all(v in remap for v in c):
                    new_cells.append([remap[v] for v in c])
            cells_list = new_cells

        # Only return coordinates in the active spatial dims (drop zero-span axes)
        spans = coords.max(axis=0) - coords.min(axis=0)
        active_axes = [i for i in range(3) if spans[i] > 1e-30]
        if not active_axes:
            active_axes = [0]
        nodes = coords[:, active_axes].tolist()

        return JSONResponse(content={
            "dim": tdim,
            "gdim": gdim,
            "nodes": nodes,
            "cells": cells_list,
            "downsampled": downsampled,
        })

    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(
            status_code=500,
            detail=f"Could not read mesh: {exc}",
        ) from exc
