"""WebSocket progress stream: /runs/{id}/stream."""
from __future__ import annotations

import asyncio
import json
from pathlib import Path

from fastapi import APIRouter, WebSocket
from starlette.websockets import WebSocketDisconnect

from ..jobs import read_status

router = APIRouter()


_POLL_INTERVAL_S = 0.05


def _read_new_lines(path: Path, pos: int) -> tuple[list[str], int]:
    if not path.exists():
        return [], pos
    with open(path) as f:
        f.seek(pos)
        lines = f.readlines()
        new_pos = f.tell()
    stripped = [ln.strip() for ln in lines if ln.strip()]
    return stripped, new_pos


def _is_terminal(status: str) -> bool:
    return status in ("completed", "failed")


@router.websocket("/runs/{run_id}/stream")
async def ws_stream(websocket: WebSocket, run_id: str) -> None:
    storage = websocket.app.state.storage
    run_dir = storage.run_dir(run_id)
    if not run_dir.exists():
        await websocket.close(code=4404)
        return

    await websocket.accept()

    progress_path = run_dir / "progress.ndjson"
    pos = 0
    saw_run_done = False

    try:
        while True:
            lines, pos = _read_new_lines(progress_path, pos)
            for line in lines:
                await websocket.send_text(line)
                try:
                    evt = json.loads(line)
                    if evt.get("type") == "run_done":
                        saw_run_done = True
                except json.JSONDecodeError:  # pragma: no cover
                    pass
            if saw_run_done:
                break

            # Also honour the status file: if the worker is finished but hasn't
            # emitted run_done (e.g. hard crash), close cleanly rather than hang.
            status_rec = read_status(run_dir)
            if status_rec and _is_terminal(status_rec.get("status", "")):
                # One more pass to catch any race-written lines.
                lines, pos = _read_new_lines(progress_path, pos)
                for line in lines:
                    await websocket.send_text(line)
                if not saw_run_done:
                    await websocket.send_text(
                        json.dumps({
                            "type": "run_done",
                            "status": status_rec["status"],
                            "wall_time_s": status_rec.get("wall_time_s", 0.0),
                        }),
                    )
                break

            await asyncio.sleep(_POLL_INTERVAL_S)
    except WebSocketDisconnect:  # pragma: no cover
        return

    try:
        await websocket.close(code=1000)
    except Exception:  # pragma: no cover
        pass
