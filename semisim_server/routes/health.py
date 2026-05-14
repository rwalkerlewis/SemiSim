"""Liveness and readiness probes."""
from __future__ import annotations

from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

from ..models import HealthResponse, ReadyResponse

router = APIRouter()


@router.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    return HealthResponse(status="ok")


@router.get(
    "/ready",
    response_model=ReadyResponse,
    responses={503: {"model": ReadyResponse}},
)
def ready(request: Request) -> JSONResponse:
    job_mgr = request.app.state.job_mgr
    is_ready = job_mgr.is_ready()
    return JSONResponse(
        content={"ready": bool(is_ready)},
        status_code=200 if is_ready else 503,
    )
