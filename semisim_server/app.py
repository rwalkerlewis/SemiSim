"""FastAPI application factory + `semisim-server` entry point.

The server process deliberately does NOT import dolfinx, UFL, or PETSc at
module scope. FEM-heavy imports happen only inside worker subprocesses
spawned by `jobs.JobManager`, so a developer without dolfinx can still run
`semisim-server` and get a useful 503 on `/ready` instead of an ImportError
at startup.
"""
from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import Settings
from .jobs import JobManager
from .storage import LocalStorage


def build_app(settings: Settings | None = None) -> FastAPI:
    if settings is None:
        settings = Settings()

    storage = LocalStorage(settings.runs_dir)
    job_mgr = JobManager(settings, storage)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        job_mgr.start()
        try:
            yield
        finally:
            job_mgr.shutdown()

    app = FastAPI(
        title="SemiSim HTTP server",
        version="0.10.0",
        description="HTTP API for the SemiSim FEM semiconductor device simulator.",
        lifespan=lifespan,
    )

    # TODO(M10+): add auth, rate limiting, API keys before any public deployment.

    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cors_origins,
        allow_credentials=False,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    app.state.settings = settings
    app.state.storage = storage
    app.state.job_mgr = job_mgr

    from .routes import health, runs, schema, solve, stream

    app.include_router(health.router)
    app.include_router(schema.router)
    app.include_router(solve.router)
    app.include_router(runs.router)
    app.include_router(stream.router)

    return app


def main() -> None:
    import uvicorn

    settings = Settings()
    uvicorn.run(build_app(settings), host=settings.host, port=settings.port)


if __name__ == "__main__":  # pragma: no cover
    main()
