"""Local-filesystem backend for run directories.

The on-disk run directory is the source of truth; the server holds no durable
state beyond what lives there. An S3 backend is a future milestone
(TODO(M10+)).
"""
from __future__ import annotations

from pathlib import Path


class LocalStorage:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def run_dir(self, run_id: str) -> Path:
        return self.root / run_id

    def list_run_ids(self) -> list[str]:
        if not self.root.exists():
            return []
        return sorted(p.name for p in self.root.iterdir() if p.is_dir())
