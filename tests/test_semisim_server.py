"""
Tests for the M10 semisim_server HTTP API.

FEM-heavy tests (real /solve end-to-end) are skipped when dolfinx is
unavailable, matching the M9 test convention. Pure-Python tests
(schema, capabilities, CORS, 404) always run.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

try:
    import dolfinx  # noqa: F401
    HAS_DOLFINX = True
except ImportError:
    HAS_DOLFINX = False

try:
    import fastapi  # noqa: F401
    import httpx  # noqa: F401
    HAS_SERVER_DEPS = True
except ImportError:
    HAS_SERVER_DEPS = False

pytestmark = pytest.mark.skipif(
    not HAS_SERVER_DEPS,
    reason="semisim_server requires fastapi + httpx (install with pip install -e '.[server]')",
)

REPO_ROOT = Path(__file__).parent.parent
BENCHMARKS_DIR = REPO_ROOT / "benchmarks"
PN_1D_JSON = BENCHMARKS_DIR / "pn_1d" / "pn_junction.json"
PN_1D_BIAS_JSON = BENCHMARKS_DIR / "pn_1d_bias" / "pn_junction_bias.json"


@pytest.fixture
def server_app(tmp_path, monkeypatch):
    """Build a fresh app rooted at a temp runs dir and exercise its lifespan."""
    monkeypatch.setenv("SEMISIM_SERVER_RUNS_DIR", str(tmp_path))
    monkeypatch.setenv("SEMISIM_SERVER_WORKERS", "2")
    monkeypatch.setenv(
        "SEMISIM_SERVER_CORS_ORIGINS",
        "http://localhost:3000,http://example.com",
    )
    from semisim_server.app import build_app
    from semisim_server.config import Settings
    settings = Settings()
    app = build_app(settings)
    return app


@pytest.fixture
def client(server_app):
    from fastapi.testclient import TestClient
    with TestClient(server_app) as c:
        yield c


def _load_body(path: Path) -> bytes:
    return path.read_bytes()


def _poll_until_done(client, run_id: str, timeout_s: float = 60.0) -> dict:
    deadline = time.monotonic() + timeout_s
    last: dict = {}
    while time.monotonic() < deadline:
        r = client.get(f"/runs/{run_id}")
        assert r.status_code == 200, r.text
        last = r.json()
        if last.get("status") in ("completed", "failed"):
            return last
        time.sleep(0.25)
    raise AssertionError(f"run {run_id} did not finish in {timeout_s}s; last={last}")


# ---------- Pure-Python tests (no dolfinx required) ----------


def test_health(client):
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}


def test_ready(client):
    r = client.get("/ready")
    assert r.status_code == 200
    assert r.json() == {"ready": True}


def test_capabilities(client):
    r = client.get("/capabilities")
    assert r.status_code == 200
    body = r.json()
    assert body["engine"]["name"] == "SemiSim"
    # M13 (transient) and M14 (ac_sweep) are shipped as of v0.14.x.
    assert "equilibrium" in body["solver_types"]
    assert "transient" in body["solver_types"]
    assert "ac_sweep" in body["solver_types"]
    assert body["dimensions"] == [1, 2, 3]
    assert body["mesh_sources"] == ["builtin", "gmsh"]
    # M15 Phase B/C: schema_version and the device_info probe are now
    # part of the response.
    assert isinstance(body["schema_version"], str)
    info = body["device_info"]
    assert info["engine_version"] == body["engine"]["version"]
    assert "cpu-mumps" in info["backends_available"]
    phys = body["physics"]
    # Hard-coded booleans reflecting the current build:
    assert phys["boltzmann"] is True
    assert phys["srh_recombination"] is True
    assert phys["transient"] is True
    assert phys["ac_small_signal"] is True
    assert phys["fermi_dirac"] is False
    assert phys["auger_recombination"] is False
    assert phys["field_dependent_mobility"] is False
    assert phys["schottky_contacts"] is False
    assert phys["heterojunctions"] is False
    assert phys["tunneling"] is False
    backends = body["backends"]
    # cpu-mumps must be reported true; GPU flags are set from the
    # runtime probe and are typically False on a CPU-only CI.
    assert backends["cpu_mumps"] is True


def test_schema(client):
    from semi.schema import (
        ENGINE_SUPPORTED_SCHEMA_MAJOR,
        ENGINE_SUPPORTED_SCHEMA_MAJORS,
        get_schema,
    )
    r = client.get("/schema")
    assert r.status_code == 200
    body = r.json()
    # /schema serves the current major's schema (strict v2 as of M14.3).
    assert body["schema"] == get_schema(ENGINE_SUPPORTED_SCHEMA_MAJOR)
    assert body["supported_major"] == ENGINE_SUPPORTED_SCHEMA_MAJOR
    assert body["supported_majors"] == list(ENGINE_SUPPORTED_SCHEMA_MAJORS)
    assert isinstance(body["version"], str)
    major, _, _ = body["version"].partition(".")
    assert int(major) == ENGINE_SUPPORTED_SCHEMA_MAJOR


def test_openapi_and_docs(client):
    r = client.get("/openapi.json")
    assert r.status_code == 200
    spec = r.json()
    assert spec["info"]["title"] == "SemiSim HTTP server"
    assert "/solve" in spec["paths"]
    assert "/runs/{run_id}" in spec["paths"]
    assert "/capabilities" in spec["paths"]


def test_solve_validates_input(client):
    # Bad JSON
    r = client.post("/solve", content=b"{ this is not json", headers={"Content-Type": "application/json"})
    assert r.status_code == 400
    body = r.json()
    assert body["detail"]["error"] == "invalid_json"

    # Schema violation: missing required `name`
    bad_cfg = {"dimension": 1}
    r = client.post("/solve", json=bad_cfg)
    assert r.status_code == 400
    body = r.json()
    assert body["detail"]["error"] == "schema_validation_failed"
    assert isinstance(body["detail"]["errors"], list)
    assert len(body["detail"]["errors"]) >= 1


def test_run_not_found(client):
    r = client.get("/runs/nonexistent_run_id_abcdef")
    assert r.status_code == 404


def test_materials(client):
    r = client.get("/materials")
    assert r.status_code == 200
    body = r.json()
    assert "Si" in body
    assert body["Si"]["role"] == "semiconductor"


def test_cors_headers(client):
    r = client.options(
        "/capabilities",
        headers={
            "Origin": "http://example.com",
            "Access-Control-Request-Method": "GET",
        },
    )
    # FastAPI/Starlette's CORSMiddleware responds 200 on preflight for allowed origins.
    assert r.status_code == 200
    assert r.headers.get("access-control-allow-origin") == "http://example.com"


# ---------- FEM-backed end-to-end tests ----------


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx required")
def test_solve_enqueues_and_completes(client, tmp_path):
    body = _load_body(PN_1D_JSON)
    r = client.post("/solve", content=body, headers={"Content-Type": "application/json"})
    assert r.status_code == 202, r.text
    data = r.json()
    run_id = data["run_id"]
    assert data["status"] == "queued"
    assert data["status_url"] == f"/runs/{run_id}"

    final = _poll_until_done(client, run_id, timeout_s=60.0)
    assert final["status"] == "completed", final

    # Manifest on disk matches what the HTTP response returned.
    runs_root = Path(client.app.state.settings.runs_dir)
    from semi.io import read_manifest
    on_disk = read_manifest(runs_root / run_id)
    for key in ("run_id", "engine", "solver", "fields", "mesh", "input_sha256"):
        assert on_disk[key] == final[key], f"mismatch on {key}"
    assert final["status"] == "completed"


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx required")
def test_field_download(client):
    body = _load_body(PN_1D_JSON)
    r = client.post("/solve", content=body, headers={"Content-Type": "application/json"})
    run_id = r.json()["run_id"]
    _poll_until_done(client, run_id)

    r = client.get(f"/runs/{run_id}/fields/psi")
    assert r.status_code == 200
    # BP directories come back as a tar archive; just assert non-empty bytes.
    assert len(r.content) > 0
    # Also verify the file really exists on disk.
    runs_root = Path(client.app.state.settings.runs_dir)
    assert (runs_root / run_id / "fields" / "psi.bp").exists()


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx required")
def test_iv_download(client):
    body = _load_body(PN_1D_BIAS_JSON)
    r = client.post("/solve", content=body, headers={"Content-Type": "application/json"})
    run_id = r.json()["run_id"]
    _poll_until_done(client, run_id, timeout_s=120.0)

    r = client.get(f"/runs/{run_id}/iv/anode")
    assert r.status_code == 200
    runs_root = Path(client.app.state.settings.runs_dir)
    disk_csv = (runs_root / run_id / "iv" / "anode.csv").read_bytes()
    assert r.content == disk_csv


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx required")
def test_input_download(client):
    body = _load_body(PN_1D_JSON)
    r = client.post("/solve", content=body, headers={"Content-Type": "application/json"})
    run_id = r.json()["run_id"]
    _poll_until_done(client, run_id)

    r = client.get(f"/runs/{run_id}/input")
    assert r.status_code == 200
    assert r.content == body


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx required")
def test_websocket_progress(client):
    """Connect during a bias sweep and assert at least n_steps messages arrive."""
    from starlette.websockets import WebSocketDisconnect

    body = _load_body(PN_1D_BIAS_JSON)
    r = client.post("/solve", content=body, headers={"Content-Type": "application/json"})
    run_id = r.json()["run_id"]

    messages: list[dict] = []
    with client.websocket_connect(f"/runs/{run_id}/stream") as ws:
        while True:
            try:
                text = ws.receive_text()
            except WebSocketDisconnect:
                break
            try:
                evt = json.loads(text)
            except json.JSONDecodeError:
                continue
            messages.append(evt)
            if evt.get("type") == "run_done":
                break

    step_done = [m for m in messages if m.get("type") == "step_done"]
    # pn_1d_bias ramps 0 -> 0.6 V in 0.05 V steps: expect well over 5 accepted
    # steps including the seed solve at V=0. At least 5 is a comfortable floor.
    assert len(step_done) >= 5, f"too few step_done messages: {messages}"
    run_done = [m for m in messages if m.get("type") == "run_done"]
    assert len(run_done) == 1
    assert run_done[0]["status"] == "completed"


@pytest.mark.skipif(not HAS_DOLFINX, reason="dolfinx required")
def test_concurrent_solves(client):
    """Fire three solves; confirm they all complete with distinct run_ids."""
    body = _load_body(PN_1D_JSON)
    ids: list[str] = []
    for _ in range(3):
        r = client.post("/solve", content=body, headers={"Content-Type": "application/json"})
        assert r.status_code == 202, r.text
        ids.append(r.json()["run_id"])
        # Force a millisecond of separation so the `ts_name_sha` prefix differs;
        # the fallback suffix handles collisions within the same second.
        time.sleep(0.01)

    assert len(set(ids)) == 3, f"run_ids collided: {ids}"

    for rid in ids:
        final = _poll_until_done(client, rid, timeout_s=120.0)
        assert final["status"] == "completed", f"{rid} -> {final}"

    # Listing surface should reflect all three runs.
    r = client.get("/runs")
    assert r.status_code == 200
    listed = {e["run_id"] for e in r.json()["runs"]}
    assert set(ids).issubset(listed)
