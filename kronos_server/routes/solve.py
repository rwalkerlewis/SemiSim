"""POST /solve: validate JSON and enqueue a solve."""
from __future__ import annotations

import json

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import JSONResponse

from ..models import SolveResponse

router = APIRouter()


@router.post(
    "/solve",
    response_model=SolveResponse,
    status_code=202,
    responses={400: {"description": "Invalid JSON or schema validation failure."}},
)
async def submit_solve(request: Request) -> JSONResponse:
    raw = await request.body()
    if not raw:
        raise HTTPException(status_code=400, detail="Empty request body.")
    try:
        cfg = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise HTTPException(
            status_code=400,
            detail={"error": "invalid_json", "message": str(exc)},
        ) from exc

    if not isinstance(cfg, dict):
        raise HTTPException(
            status_code=400,
            detail={"error": "invalid_input", "message": "Top-level JSON must be an object."},
        )

    # Validate against the input schema before enqueuing.
    from semi.schema import SchemaError, get_schema, validate

    # Pick the schema matching the submitted major version so that v2 inputs
    # (which may use caughey_thomas mobility, ac_sweep, etc.) are validated
    # against input.v2.json rather than the legacy v1 schema.
    raw_version = cfg.get("schema_version", "1.0.0")
    try:
        schema_major = int(str(raw_version).split(".")[0])
    except (ValueError, AttributeError):
        schema_major = 1
    versioned_schema = get_schema(schema_major)

    try:
        import jsonschema
        validator = jsonschema.Draft7Validator(versioned_schema)
        errors = sorted(validator.iter_errors(cfg), key=lambda e: list(e.path))
        if errors:
            raise HTTPException(
                status_code=400,
                detail={
                    "error": "schema_validation_failed",
                    "errors": [
                        {
                            "path": ".".join(str(p) for p in e.absolute_path) or "<root>",
                            "message": e.message,
                        }
                        for e in errors
                    ],
                },
            )
    except ImportError:  # pragma: no cover
        try:
            validate(cfg)
        except SchemaError as exc:
            raise HTTPException(
                status_code=400,
                detail={"error": "schema_validation_failed", "message": str(exc)},
            ) from exc

    # Fill defaults (mirrors what `semi-run` does on disk) before handing to the worker.
    # The worker reloads the file via schema.load() which re-applies defaults; pass the
    # raw bytes unchanged so input.json on disk matches what the caller posted.
    cfg_with_defaults = validate(cfg)

    job_mgr = request.app.state.job_mgr
    if not job_mgr.is_ready():
        raise HTTPException(status_code=503, detail="Worker pool not ready.")

    run_id = job_mgr.submit(raw, cfg_with_defaults)
    response = SolveResponse(
        run_id=run_id,
        status="queued",
        status_url=f"/runs/{run_id}",
    )
    return JSONResponse(
        status_code=202,
        content=response.model_dump(),
    )
