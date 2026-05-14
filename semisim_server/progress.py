"""File-backed progress event stream.

The worker appends newline-delimited JSON to `<run_dir>/progress.ndjson`;
the WebSocket route tails that file. Using the filesystem as the transport
avoids shared-queue plumbing across a ProcessPoolExecutor boundary and keeps
the run directory self-describing.
"""
from __future__ import annotations

import datetime
import json
from pathlib import Path
from typing import Any


def iso_now() -> str:
    now = datetime.datetime.now(datetime.timezone.utc)
    return now.strftime("%Y-%m-%dT%H:%M:%S.") + f"{now.microsecond // 1000:03d}Z"


def append_event(run_dir: Path, event: dict[str, Any]) -> None:
    """Append a single event dict as one JSON line."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    payload = dict(event)
    payload.setdefault("timestamp", iso_now())
    line = json.dumps(payload, sort_keys=True)
    path = run_dir / "progress.ndjson"
    with open(path, "a") as f:
        f.write(line + "\n")
