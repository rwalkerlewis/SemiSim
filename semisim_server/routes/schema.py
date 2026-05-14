"""Static metadata endpoints: /schema, /materials, /capabilities."""
from __future__ import annotations

from typing import Any

from fastapi import APIRouter

from ..models import (
    CapabilitiesBackends,
    CapabilitiesDeviceInfo,
    CapabilitiesEngine,
    CapabilitiesPhysics,
    CapabilitiesResponse,
)

router = APIRouter()


@router.get("/schema")
def get_schema() -> dict[str, Any]:
    """The current input JSON schema plus version metadata.

    UI form-builders consume this. Returns
    ``{"schema": <Draft-07 schema>, "version": "2.0.0",
    "supported_major": 2, "supported_majors": [1, 2]}``. The advertised
    ``schema`` is the strict v2 schema (current); v1 callers can still
    submit inputs and will get a DeprecationWarning back from the
    engine.
    """
    from semi.schema import (
        ENGINE_SUPPORTED_SCHEMA_MAJOR,
        ENGINE_SUPPORTED_SCHEMA_MAJORS,
        get_schema as _get_schema,
    )

    schema = _get_schema(ENGINE_SUPPORTED_SCHEMA_MAJOR)
    version = (
        schema.get("properties", {})
        .get("schema_version", {})
        .get("examples", ["2.0.0"])[0]
    )
    return {
        "schema": schema,
        "version": version,
        "supported_major": ENGINE_SUPPORTED_SCHEMA_MAJOR,
        "supported_majors": list(ENGINE_SUPPORTED_SCHEMA_MAJORS),
    }


@router.get("/materials")
def get_materials() -> dict[str, dict[str, Any]]:
    """Dump the material database: {name: parameters}."""
    from dataclasses import asdict

    from semi.materials import MATERIALS
    return {name: asdict(mat) for name, mat in MATERIALS.items()}


@router.get("/capabilities", response_model=CapabilitiesResponse)
def get_capabilities() -> CapabilitiesResponse:
    """Feature flags describing what this engine build supports.

    The static physics flags are hard-coded at the source level (see
    docs/IMPROVEMENT_GUIDE.md §4). The ``backends`` and ``device_info``
    fields are produced at request time by ``semi.compute.device_info``,
    so a UI can hide GPU backend choices on a CPU-only PETSc build.
    """
    try:
        from semi import __version__
        version = str(__version__)
    except Exception:  # pragma: no cover
        version = "unknown"

    from semi.compute import device_info as _device_info
    from semi.schema import ENGINE_SUPPORTED_SCHEMA_MAJOR, get_schema as _get_schema
    info = _device_info()
    available = set(info["backends_available"])
    schema = _get_schema(ENGINE_SUPPORTED_SCHEMA_MAJOR)
    schema_version = (
        schema.get("properties", {})
        .get("schema_version", {})
        .get("examples", ["2.0.0"])[-1]
    )

    return CapabilitiesResponse(
        engine=CapabilitiesEngine(name="SemiSim", version=version),
        schema_version=schema_version,
        solver_types=[
            "equilibrium",
            "drift_diffusion",
            "bias_sweep",
            "mos_cv",
            "mos_cap_ac",
            "transient",
            "ac_sweep",
        ],
        dimensions=[1, 2, 3],
        physics=CapabilitiesPhysics(
            boltzmann=True,
            fermi_dirac=False,
            srh_recombination=True,
            auger_recombination=False,
            field_dependent_mobility=False,
            transient=True,
            ac_small_signal=True,
            schottky_contacts=False,
            heterojunctions=False,
            tunneling=False,
        ),
        mesh_sources=["builtin", "gmsh"],
        backends=CapabilitiesBackends(
            cpu_mumps="cpu-mumps" in available,
            gpu_amgx="gpu-amgx" in available,
            gpu_hypre="gpu-hypre" in available,
        ),
        device_info=CapabilitiesDeviceInfo(**info),
    )
