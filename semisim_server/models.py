"""Pydantic request/response models surfaced in the OpenAPI spec."""
from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class SolveResponse(BaseModel):
    run_id: str = Field(..., description="Identifier for the newly-enqueued run.")
    status: str = Field(..., description="Queued/running/completed/failed.")
    status_url: str = Field(..., description="Relative URL for polling run status.")


class RunSummary(BaseModel):
    run_id: str
    status: str
    created_at: str = Field(..., description="ISO-8601 UTC timestamp parsed from run_id.")


class RunListResponse(BaseModel):
    runs: list[RunSummary]


class HealthResponse(BaseModel):
    status: str


class ReadyResponse(BaseModel):
    ready: bool


class CapabilitiesEngine(BaseModel):
    name: str
    version: str


class CapabilitiesPhysics(BaseModel):
    boltzmann: bool
    fermi_dirac: bool
    srh_recombination: bool
    auger_recombination: bool
    field_dependent_mobility: bool
    transient: bool
    ac_small_signal: bool
    schottky_contacts: bool
    heterojunctions: bool
    tunneling: bool


class CapabilitiesBackends(BaseModel):
    cpu_mumps: bool
    gpu_amgx: bool
    gpu_hypre: bool


class CapabilitiesDeviceInfo(BaseModel):
    """Runtime probe of the linked PETSc build (M15 Phase C).

    Mirrors :func:`semi.compute.device_info`. Reported by
    ``GET /capabilities`` so a UI can decide whether to expose GPU
    backend choices.
    """
    engine_version: str
    petsc_version: str | None
    petsc_complex: bool
    petsc_int64: bool
    backends_available: list[str]
    device_count: int | None
    device_name: str | None


class CapabilitiesResponse(BaseModel):
    engine: CapabilitiesEngine
    schema_version: str
    solver_types: list[str]
    dimensions: list[int]
    physics: CapabilitiesPhysics
    mesh_sources: list[str]
    backends: CapabilitiesBackends
    device_info: CapabilitiesDeviceInfo


class ErrorResponse(BaseModel):
    detail: Any
