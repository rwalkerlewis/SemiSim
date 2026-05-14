"""Runtime settings sourced from environment variables."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


def _default_workers() -> int:
    cpus = os.cpu_count() or 4
    return max(1, cpus // 4)


def _parse_origins(raw: str) -> list[str]:
    return [o.strip() for o in raw.split(",") if o.strip()]


@dataclass
class Settings:
    host: str = field(
        default_factory=lambda: os.environ.get("SEMISIM_SERVER_HOST", "127.0.0.1"),
    )
    port: int = field(
        default_factory=lambda: int(os.environ.get("SEMISIM_SERVER_PORT", "8000")),
    )
    workers: int = field(
        default_factory=lambda: int(
            os.environ.get("SEMISIM_SERVER_WORKERS", str(_default_workers())),
        ),
    )
    runs_dir: Path = field(
        default_factory=lambda: Path(os.environ.get("SEMISIM_SERVER_RUNS_DIR", "runs")),
    )
    cors_origins: list[str] = field(
        default_factory=lambda: _parse_origins(
            os.environ.get(
                "SEMISIM_SERVER_CORS_ORIGINS",
                "http://localhost:3000,http://localhost:5173",
            ),
        ),
    )
