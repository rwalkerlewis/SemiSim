"""
Result artifact writer for kronos-semi simulations.

Output layout under out_dir/<run_id>/:
    manifest.json     schema-versioned metadata
    input.json        copy of input JSON
    mesh/mesh.xdmf    mesh topology + geometry
    mesh/mesh.h5      mesh HDF5 companion
    fields/<name>.bp  scalar/vector fields (ADIOS2 BP5, or .xdmf fallback)
    iv/<contact>.csv  V,J_n,J_p,J_total per sweep step
    convergence/snes.csv  per-step Newton convergence data
    logs/             placeholder directory

No dolfinx import at module top level (Architecture Invariant 4).
"""
from __future__ import annotations

import csv
import datetime
import hashlib
import json
import subprocess
import time
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..run import SimulationResult

from .reader import read_manifest  # re-export

_SCHEMA_VERSION = "1.1.0"

__all__ = ["write_artifact", "read_manifest"]


def write_artifact(
    result: SimulationResult,
    out_dir: Path,
    run_id: str | None = None,
    input_json_path: Path | None = None,
) -> Path:
    """
    Write a versioned run artifact directory and return the run_dir Path.

    Parameters
    ----------
    result:
        Completed SimulationResult from semi.run.run().
    out_dir:
        Parent directory under which the run subdirectory is created.
    run_id:
        Override for the auto-generated run ID string. If None, one is
        derived from the current UTC timestamp, the config name, and a
        short SHA of the input bytes.
    input_json_path:
        Path to the original input JSON file. When supplied, the file bytes
        are used verbatim for the SHA-256 digest and for writing input.json.
        When None, result.cfg is serialised to JSON bytes (sort_keys=True).

    Returns
    -------
    Path
        Absolute path to the run directory (out_dir / run_id).
    """
    t_start = time.monotonic()

    # 1. Compute input bytes and sha256
    if input_json_path is not None:
        input_bytes = Path(input_json_path).read_bytes()
    else:
        input_bytes = json.dumps(result.cfg, sort_keys=True).encode()
    input_sha256 = hashlib.sha256(input_bytes).hexdigest()

    # 2. Build run_id if not supplied
    if run_id is None:
        short_sha = input_sha256[:7]
        ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H-%M-%SZ")
        name = result.cfg.get("name", "run")
        run_id = f"{ts}_{name}_{short_sha}"

    # 3. Create run_dir and subdirectories
    out_dir = Path(out_dir)
    run_dir = out_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    for subdir in ("fields", "mesh", "iv", "convergence", "logs"):
        (run_dir / subdir).mkdir(exist_ok=True)

    # 4. Write input.json
    (run_dir / "input.json").write_bytes(input_bytes)

    # 5. Write fields
    warnings: list[str] = []
    field_entries = _write_fields(result, run_dir, warnings)

    # 6. Write mesh XDMF
    _write_mesh(result, run_dir, warnings)

    # 7. Write IV CSV
    sweep_entries = _write_iv_csv(result, run_dir)

    # 8. Write convergence CSV
    _write_convergence_csv(result, run_dir)

    # 9. Build and write manifest
    wall_time_s = time.monotonic() - t_start

    from .. import __version__

    engine = {
        "name": "kronos-semi",
        "version": __version__,
        "commit": _git_commit(),
    }

    stype = result.cfg.get("solver", {}).get("type", "equilibrium")
    n_steps = max(1, len(result.iv))
    converged = bool((result.solver_info or {}).get("converged", True))
    n_dofs = int(result.V.dofmap.index_map.size_global)

    # M15 backend metadata (optional; populated when the runner ran a
    # solve through semi.solver.solve_nonlinear[_block]). Falls back to
    # the pre-M15 "petsc-cpu-mumps" string when absent so old reader
    # paths still find a sensible value in `solver.backend`.
    info = result.solver_info or {}
    backend_resolved = info.get("backend_resolved", "cpu-mumps")
    backend_requested = info.get("backend_requested", backend_resolved)
    device = info.get("device", "cpu" if backend_resolved == "cpu-mumps" else "auto")
    backend_str = f"petsc-{backend_resolved}"

    dim = result.mesh.topology.dim
    n_cells = int(result.mesh.topology.index_map(dim).size_global)
    n_vertices = int(result.mesh.topology.index_map(0).size_global)

    # Build regions list from cfg
    regions = []
    for region_name, region_data in result.cfg.get("regions", {}).items():
        entry: dict = {
            "name": region_name,
            "material": region_data.get("material", ""),
        }
        if "tag" in region_data:
            entry["tag"] = int(region_data["tag"])
        regions.append(entry)

    manifest: dict = {
        "schema_version": _SCHEMA_VERSION,
        "engine": engine,
        "run_id": run_id,
        "status": "completed",
        "wall_time_s": wall_time_s,
        "input_sha256": input_sha256,
        "solver": {
            "type": stype,
            "backend": backend_str,
            "backend_requested": backend_requested,
            "backend_resolved": backend_resolved,
            "device": device,
            "linear_solver": info.get("linear_solver", "auto"),
            "preconditioner": info.get("preconditioner", "auto"),
            "ksp_iters": int(info.get("ksp_iters", 0)),
            "linear_solve_wall_s": float(info.get("linear_solve_wall_s", 0.0)),
            "n_dofs": n_dofs,
            "n_steps": n_steps,
            "converged": converged,
        },
        "fields": field_entries,
        "mesh": {
            "path": "mesh/mesh.xdmf",
            "dimension": dim,
            "n_cells": n_cells,
            "n_vertices": n_vertices,
            "regions": regions,
        },
        "warnings": warnings,
    }

    if sweep_entries:
        manifest["sweeps"] = sweep_entries

    manifest_text = json.dumps(manifest, sort_keys=True, indent=2) + "\n"
    (run_dir / "manifest.json").write_text(manifest_text)

    return run_dir


def _write_fields(result, run_dir, warnings):
    """Write all exported fields; return list of field-entry dicts."""
    # dolfinx imports inside
    import ufl
    from dolfinx import fem

    msh = result.mesh
    sc = result.scaling
    fields_dir = run_dir / "fields"
    stype = result.cfg.get("solver", {}).get("type", "equilibrium")
    is_dd = stype in ("drift_diffusion", "bias_sweep")

    field_entries = []
    use_bp = [True]  # mutable list so inner closure can flip it

    def _write(fn, name):
        if use_bp[0]:
            bp_path = str(fields_dir / f"{name}.bp")
            try:
                _vtx_write(msh, fn, bp_path)
                return f"fields/{name}.bp"
            except Exception as exc:  # pragma: no cover
                warnings.append(f"ADIOS2/VTX unavailable ({exc}); fields written as XDMF")
                use_bp[0] = False
        xdmf_path = str(fields_dir / f"{name}.xdmf")  # pragma: no cover
        try:  # pragma: no cover
            _xdmf_write(msh, fn, xdmf_path)
            return f"fields/{name}.xdmf"
        except Exception as exc:
            warnings.append(f"Field {name} write failed (XDMF): {exc}")
            return None

    # psi (physical)
    psi_fn = fem.Function(result.V, name="psi")
    psi_fn.x.array[:] = result.psi.x.array * sc.V0
    psi_fn.x.scatter_forward()
    p = _write(psi_fn, "psi")
    if p:
        field_entries.append({"name": "psi", "units": "V", "path": p, "rank": 0})

    # n
    n_fn = fem.Function(result.V, name="n")
    n_fn.x.array[:] = result.n_phys
    n_fn.x.scatter_forward()
    p = _write(n_fn, "n")
    if p:
        field_entries.append({"name": "n", "units": "m^-3", "path": p, "rank": 0})

    # p
    p_fn = fem.Function(result.V, name="p")
    p_fn.x.array[:] = result.p_phys
    p_fn.x.scatter_forward()
    path = _write(p_fn, "p")
    if path:
        field_entries.append({"name": "p", "units": "m^-3", "path": path, "rank": 0})

    # N_net (physical)
    N_net_fn = fem.Function(result.V, name="N_net")
    N_net_fn.x.array[:] = result.N_hat.x.array * sc.C0
    N_net_fn.x.scatter_forward()
    path = _write(N_net_fn, "N_net")
    if path:
        field_entries.append({"name": "N_net", "units": "m^-3", "path": path, "rank": 0})

    # E = -grad(psi_phys) as DG0 vector
    try:
        dim = msh.topology.dim
        V_dg = fem.functionspace(msh, ("DG", 0, (dim,)))
        E_fn = fem.Function(V_dg, name="E")
        E_expr = fem.Expression(
            -sc.V0 * ufl.grad(result.psi),
            V_dg.element.interpolation_points(),
        )
        E_fn.interpolate(E_expr)  # pragma: no cover
        path = _write(E_fn, "E")  # pragma: no cover
        if path:  # pragma: no cover
            field_entries.append({"name": "E", "units": "V/m", "path": path, "rank": 1})
    except Exception as exc:
        warnings.append(f"E field export failed: {exc}")

    # phi_n, phi_p, J_n, J_p for DD runs
    if is_dd and result.phi_n is not None and result.phi_p is not None:
        phys_cfg = result.cfg.get("physics", {})
        mob = phys_cfg.get("mobility", {})
        mu_n_SI = float(mob.get("mu_n", 1400.0)) * 1.0e-4
        mu_p_SI = float(mob.get("mu_p", 450.0)) * 1.0e-4

        phi_n_fn = fem.Function(result.phi_n.function_space, name="phi_n")
        phi_n_fn.x.array[:] = result.phi_n.x.array * sc.V0
        phi_n_fn.x.scatter_forward()
        path = _write(phi_n_fn, "phi_n")
        if path:
            field_entries.append({"name": "phi_n", "units": "V", "path": path, "rank": 0})

        phi_p_fn = fem.Function(result.phi_p.function_space, name="phi_p")
        phi_p_fn.x.array[:] = result.phi_p.x.array * sc.V0
        phi_p_fn.x.scatter_forward()
        path = _write(phi_p_fn, "phi_p")
        if path:
            field_entries.append({"name": "phi_p", "units": "V", "path": path, "rank": 0})

        try:
            from petsc4py import PETSc

            from ..constants import Q
            from ..physics.slotboom import n_from_slotboom, p_from_slotboom
            dim = msh.topology.dim
            V_dg_vec = fem.functionspace(msh, ("DG", 0, (dim,)))
            ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))

            Jn_ufl = (
                Q * mu_n_SI
                * n_from_slotboom(result.psi, result.phi_n, ni_hat)
                * sc.C0
                * sc.V0 * ufl.grad(result.phi_n)
            )
            Jn_fn = fem.Function(V_dg_vec, name="J_n")
            Jn_fn.interpolate(fem.Expression(Jn_ufl, V_dg_vec.element.interpolation_points()))
            path = _write(Jn_fn, "J_n")  # pragma: no cover
            if path:  # pragma: no cover
                field_entries.append({"name": "J_n", "units": "A/m^2", "path": path, "rank": 1})

            Jp_ufl = (  # pragma: no cover
                Q * mu_p_SI
                * p_from_slotboom(result.psi, result.phi_p, ni_hat)
                * sc.C0
                * sc.V0 * ufl.grad(result.phi_p)
            )
            Jp_fn = fem.Function(V_dg_vec, name="J_p")  # pragma: no cover
            Jp_fn.interpolate(fem.Expression(Jp_ufl, V_dg_vec.element.interpolation_points()))  # pragma: no cover
            path = _write(Jp_fn, "J_p")  # pragma: no cover
            if path:  # pragma: no cover
                field_entries.append({"name": "J_p", "units": "A/m^2", "path": path, "rank": 1})
        except Exception as exc:
            warnings.append(f"J_n/J_p field export failed: {exc}")

    return field_entries


def _write_mesh(result, run_dir, warnings):
    """Export mesh topology and geometry to mesh/mesh.xdmf."""
    try:
        from dolfinx.io import XDMFFile
        mesh_path = str(run_dir / "mesh" / "mesh.xdmf")
        with XDMFFile(result.mesh.comm, mesh_path, "w") as xdmf:
            xdmf.write_mesh(result.mesh)
    except Exception as exc:  # pragma: no cover
        warnings.append(f"Mesh XDMF write failed: {exc}")


def _vtx_write(msh, fn, path_str: str) -> None:
    from dolfinx.io import VTXWriter
    with VTXWriter(msh.comm, path_str, [fn], engine="BP5") as vtx:
        vtx.write(0.0)


def _xdmf_write(msh, fn, path_str: str) -> None:  # pragma: no cover
    from dolfinx.io import XDMFFile
    with XDMFFile(msh.comm, path_str, "w") as xdmf:
        xdmf.write_mesh(msh)
        xdmf.write_function(fn)


def _write_iv_csv(result, run_dir) -> list[dict]:
    if not result.iv or result.bias_contact is None:
        return []
    contact = result.bias_contact
    iv_path = run_dir / "iv" / f"{contact}.csv"
    with open(iv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["V", "J_n", "J_p", "J_total"])
        for row in result.iv:
            writer.writerow([row["V"], float("nan"), float("nan"), row.get("J", 0.0)])

    stype = result.cfg.get("solver", {}).get("type", "equilibrium")
    kind = "gate" if stype == "mos_cv" else "voltage"
    vs = [row["V"] for row in result.iv]
    bipolar = any(v < -1.0e-12 for v in vs) and any(v > 1.0e-12 for v in vs)
    return [{
        "kind": kind,
        "contact": contact,
        "path": f"iv/{contact}.csv",
        "n_steps": len(result.iv),
        "bipolar": bipolar,
    }]


def _write_convergence_csv(result, run_dir) -> None:
    conv_path = run_dir / "convergence" / "snes.csv"
    info = result.solver_info or {}
    last_iters = int(info.get("iterations", -1))
    last_reason = info.get("reason", "unknown")
    last_conv = bool(info.get("converged", True))

    rows = []
    if result.iv:
        for i, iv_row in enumerate(result.iv):
            is_last = i == len(result.iv) - 1
            rows.append({
                "bias_step": i,
                "V_applied": iv_row["V"],
                "iterations": last_iters if is_last else -1,
                "reason": last_reason if is_last else "unknown",
                "converged": last_conv if is_last else True,
            })
    else:
        rows.append({
            "bias_step": 0,
            "V_applied": 0.0,
            "iterations": last_iters,
            "reason": last_reason,
            "converged": last_conv,
        })

    with open(conv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["bias_step", "V_applied", "iterations", "reason", "converged"],
        )
        writer.writeheader()
        writer.writerows(rows)


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
            cwd=Path(__file__).parent.parent.parent,
        ).strip()
    except Exception:  # pragma: no cover
        return "unknown"


def write_transient_artifact(
    result,
    out_dir: Path,
    run_id: str | None = None,
    input_json_path: Path | None = None,
) -> Path:
    """Write a minimal versioned artifact for a ``TransientResult``.

    Writes:
      input.json, manifest.json, iv/<contact>.csv (t,V,J_n,J_p,J_total)

    Returns the run directory Path.
    """
    from .. import __version__

    t_start = time.monotonic()

    if input_json_path is not None:
        input_bytes = Path(input_json_path).read_bytes()
        cfg = json.loads(input_bytes)
    else:
        cfg = {}
        input_bytes = b"{}"
    input_sha256 = hashlib.sha256(input_bytes).hexdigest()

    if run_id is None:
        short_sha = input_sha256[:7]
        ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H-%M-%SZ")
        name = cfg.get("name", "run")
        run_id = f"{ts}_{name}_{short_sha}"

    run_dir = Path(out_dir) / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    for subdir in ("iv", "logs"):
        (run_dir / subdir).mkdir(exist_ok=True)

    (run_dir / "input.json").write_bytes(input_bytes)

    # Collect unique contacts from iv rows.
    iv_rows = result.iv or []
    contacts = list(dict.fromkeys(r.get("contact", "unknown") for r in iv_rows))

    sweep_entries = []
    for contact in contacts:
        rows = [r for r in iv_rows if r.get("contact") == contact]
        iv_path = run_dir / "iv" / f"{contact}.csv"
        with open(iv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["t", "V", "J_n", "J_p", "J_total"])
            for row in rows:
                writer.writerow([
                    row.get("t", 0.0),
                    row.get("V", 0.0),
                    row.get("J_n", float("nan")),
                    row.get("J_p", float("nan")),
                    row.get("J", row.get("J_total", float("nan"))),
                ])
        sweep_entries.append({
            "kind": "transient",
            "contact": contact,
            "path": f"iv/{contact}.csv",
            "n_steps": len(rows),
            "bipolar": False,
        })

    meta = result.meta or {}
    manifest = {
        "schema_version": _SCHEMA_VERSION,
        "engine": {
            "name": "kronos-semi",
            "version": __version__,
            "commit": _git_commit(),
        },
        "run_id": run_id,
        "status": "completed",
        "wall_time_s": time.monotonic() - t_start,
        "input_sha256": input_sha256,
        "solver": {
            "type": "transient",
            "n_steps": int(meta.get("n_steps_taken", len(result.t))),
            "dt": float(meta.get("dt", 0.0)),
            "order": int(meta.get("order", 1)),
        },
        "fields": [],
        "warnings": [],
    }
    if sweep_entries:
        manifest["sweeps"] = sweep_entries

    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True, indent=2) + "\n"
    )
    return run_dir


def write_ac_sweep_artifact(
    result,
    out_dir: Path,
    run_id: str | None = None,
    input_json_path: Path | None = None,
) -> Path:
    """Write a minimal versioned artifact for an ``AcSweepResult``.

    Writes:
      input.json, manifest.json, iv/<contact>.csv (f_Hz,C,G,Re_Y,Im_Y)

    Returns the run directory Path.
    """
    from .. import __version__

    t_start = time.monotonic()

    if input_json_path is not None:
        input_bytes = Path(input_json_path).read_bytes()
        cfg = json.loads(input_bytes)
    else:
        cfg = {}
        input_bytes = b"{}"
    input_sha256 = hashlib.sha256(input_bytes).hexdigest()

    if run_id is None:
        short_sha = input_sha256[:7]
        ts = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H-%M-%SZ")
        name = cfg.get("name", "run")
        run_id = f"{ts}_{name}_{short_sha}"

    run_dir = Path(out_dir) / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    for subdir in ("iv", "logs"):
        (run_dir / subdir).mkdir(exist_ok=True)

    (run_dir / "input.json").write_bytes(input_bytes)

    # Derive contact name from dc_bias metadata.
    dc = result.dc_bias or {}
    contact = dc.get("contact", "ac_port")
    iv_path = run_dir / "iv" / f"{contact}.csv"
    sweep_entries = []
    if result.frequencies:
        with open(iv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["f_Hz", "C", "G", "Re_Y", "Im_Y"])
            for freq, Y_val, C_val, G_val in zip(
                result.frequencies, result.Y, result.C, result.G, strict=False
            ):
                writer.writerow([freq, C_val, G_val, Y_val.real, Y_val.imag])
        sweep_entries.append({
            "kind": "ac_sweep",
            "contact": contact,
            "path": f"iv/{contact}.csv",
            "n_steps": len(result.frequencies),
            "bipolar": False,
        })

    manifest = {
        "schema_version": _SCHEMA_VERSION,
        "engine": {
            "name": "kronos-semi",
            "version": __version__,
            "commit": _git_commit(),
        },
        "run_id": run_id,
        "status": "completed",
        "wall_time_s": time.monotonic() - t_start,
        "input_sha256": input_sha256,
        "solver": {
            "type": "ac_sweep",
            "n_freqs": len(result.frequencies),
            "dc_contact": contact,
            "dc_voltage": float(dc.get("voltage", 0.0)),
        },
        "fields": [],
        "warnings": [],
    }
    if sweep_entries:
        manifest["sweeps"] = sweep_entries

    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True, indent=2) + "\n"
    )
    return run_dir
