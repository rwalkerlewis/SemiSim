"""Job queue + worker function.

The server process never imports dolfinx. `_solve_job` runs in a fresh
subprocess spawned by `ProcessPoolExecutor` (spawn context, so PETSc state
is not shared across jobs, per IMPROVEMENT_GUIDE §2 P3). Progress events
are streamed to disk via `progress.append_event`; the WebSocket route tails
that file.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import multiprocessing
import subprocess
import time
import traceback
from concurrent.futures import Future, ProcessPoolExecutor
from pathlib import Path
from typing import Any

from .config import Settings
from .progress import append_event
from .storage import LocalStorage


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
            cwd=Path(__file__).parent.parent,
        ).strip()
    except Exception:  # pragma: no cover
        return "unknown"


def _get_version() -> str:
    try:
        from semi import __version__
        return str(__version__)
    except Exception:  # pragma: no cover
        return "unknown"


def _status_path(run_dir: Path) -> Path:
    return run_dir / "status.json"


def read_status(run_dir: Path) -> dict[str, Any] | None:
    path = _status_path(run_dir)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _utc_now_iso() -> str:
    now = datetime.datetime.now(datetime.timezone.utc)
    return now.strftime("%Y-%m-%dT%H:%M:%SZ")


def write_status(run_dir: Path, status: str, **extra: Any) -> None:
    path = _status_path(run_dir)
    payload = {"status": status, "updated_at": _utc_now_iso()}
    payload.update(extra)
    path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")


def _solve_job(run_id: str, run_dir_str: str, input_path_str: str) -> dict[str, Any]:
    """Subprocess worker entry point. Imports dolfinx only inside."""
    run_dir = Path(run_dir_str)
    input_path = Path(input_path_str)
    t_start = time.monotonic()

    def progress_cb(event: dict[str, Any]) -> None:
        append_event(run_dir, event)

    write_status(run_dir, "running", run_id=run_id)
    progress_cb({"type": "run_started", "run_id": run_id})

    try:
        from semi.io import write_artifact, write_ac_sweep_artifact, write_transient_artifact
        from semi.runners import run_ac_sweep, run_bias_sweep, run_equilibrium, run_mos_cap_ac, run_mos_cv, run_transient
        from semi.schema import load as schema_load

        cfg = schema_load(str(input_path))
        stype = cfg.get("solver", {}).get("type", "equilibrium")

        if stype == "equilibrium":
            result = run_equilibrium(cfg, progress_callback=progress_cb)
        elif stype in ("drift_diffusion", "bias_sweep"):
            result = run_bias_sweep(cfg, progress_callback=progress_cb)
        elif stype == "mos_cv":
            result = run_mos_cv(cfg, progress_callback=progress_cb)
        elif stype == "mos_cap_ac":
            result = run_mos_cap_ac(cfg, progress_callback=progress_cb)
        elif stype == "ac_sweep":
            result = run_ac_sweep(cfg, progress_callback=progress_cb)
            final_dir = write_ac_sweep_artifact(
                result,
                out_dir=run_dir.parent,
                run_id=run_id,
                input_json_path=input_path,
            )
            if str(final_dir) != str(run_dir):  # pragma: no cover
                raise RuntimeError(
                    f"write_ac_sweep_artifact wrote to {final_dir}, expected {run_dir}",
                )
            wall = time.monotonic() - t_start
            write_status(run_dir, "completed", run_id=run_id, wall_time_s=wall)
            progress_cb({"type": "run_done", "status": "completed", "wall_time_s": wall})
            return {"status": "completed", "wall_time_s": wall}
        elif stype == "transient":
            result = run_transient(cfg, progress_callback=progress_cb)
            final_dir = write_transient_artifact(
                result,
                out_dir=run_dir.parent,
                run_id=run_id,
                input_json_path=input_path,
            )
            if str(final_dir) != str(run_dir):  # pragma: no cover
                raise RuntimeError(
                    f"write_transient_artifact wrote to {final_dir}, expected {run_dir}",
                )
            wall = time.monotonic() - t_start
            write_status(run_dir, "completed", run_id=run_id, wall_time_s=wall)
            progress_cb({"type": "run_done", "status": "completed", "wall_time_s": wall})
            return {"status": "completed", "wall_time_s": wall}
        else:
            raise ValueError(f"Unknown solver.type {stype!r}")

        final_dir = write_artifact(
            result,
            out_dir=run_dir.parent,
            run_id=run_id,
            input_json_path=input_path,
        )
        if str(final_dir) != str(run_dir):  # pragma: no cover
            raise RuntimeError(
                f"write_artifact wrote to {final_dir}, expected {run_dir}",
            )

        wall = time.monotonic() - t_start
        write_status(run_dir, "completed", run_id=run_id, wall_time_s=wall)
        progress_cb({"type": "run_done", "status": "completed", "wall_time_s": wall})
        return {"status": "completed", "wall_time_s": wall}
    except Exception as exc:
        wall = time.monotonic() - t_start
        tb = traceback.format_exc()
        progress_cb({
            "type": "run_done",
            "status": "failed",
            "wall_time_s": wall,
            "error": str(exc),
        })
        _write_failure_manifest(run_dir, run_id, input_path, exc, tb, wall)
        write_status(run_dir, "failed", run_id=run_id, wall_time_s=wall, error=str(exc))
        return {"status": "failed", "error": str(exc)}


def _write_failure_manifest(
    run_dir: Path,
    run_id: str,
    input_path: Path,
    exc: Exception,
    tb: str,
    wall: float,
) -> None:
    """Emit a best-effort manifest on failure so `GET /runs/{id}` still works."""
    try:
        input_bytes = input_path.read_bytes() if input_path.exists() else b""
        sha = hashlib.sha256(input_bytes).hexdigest() if input_bytes else "0" * 64
        manifest = {
            "schema_version": "1.0.0",
            "engine": {"name": "kronos-semi", "version": _get_version(), "commit": _git_commit()},
            "run_id": run_id,
            "status": "failed",
            "wall_time_s": float(wall),
            "input_sha256": sha,
            "error": {"message": str(exc), "traceback": tb},
        }
        (run_dir / "manifest.json").write_text(
            json.dumps(manifest, sort_keys=True, indent=2) + "\n",
        )
    except Exception:  # pragma: no cover
        pass


class JobManager:
    """Process-pool-backed job dispatcher."""

    def __init__(self, settings: Settings, storage: LocalStorage):
        self.settings = settings
        self.storage = storage
        self.executor: ProcessPoolExecutor | None = None
        self._futures: dict[str, Future] = {}

    def start(self) -> None:
        ctx = multiprocessing.get_context("spawn")
        self.executor = ProcessPoolExecutor(
            max_workers=self.settings.workers,
            mp_context=ctx,
        )

    def shutdown(self) -> None:
        if self.executor is not None:
            self.executor.shutdown(wait=False, cancel_futures=True)
            self.executor = None

    def is_ready(self) -> bool:
        return self.executor is not None

    def submit(self, input_bytes: bytes, cfg: dict[str, Any]) -> str:
        if self.executor is None:
            raise RuntimeError("JobManager not started")
        sha = hashlib.sha256(input_bytes).hexdigest()
        short = sha[:7]
        name = cfg.get("name", "run")
        ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H-%M-%SZ")
        # Add a short counter suffix if this exact run_id is already taken,
        # to keep concurrent submissions from colliding within a 1-second slot.
        base = f"{ts}_{name}_{short}"
        run_id = base
        suffix = 1
        while self.storage.run_dir(run_id).exists():
            run_id = f"{base}_{suffix}"
            suffix += 1

        run_dir = self.storage.run_dir(run_id)
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "input.json").write_bytes(input_bytes)
        write_status(run_dir, "queued", run_id=run_id, input_sha256=sha)

        fut = self.executor.submit(_solve_job, run_id, str(run_dir), str(run_dir / "input.json"))
        self._futures[run_id] = fut
        return run_id

    def future(self, run_id: str) -> Future | None:
        return self._futures.get(run_id)
