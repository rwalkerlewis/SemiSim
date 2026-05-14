"""Tests for semi.compute (M15 Phase C).

Pure-Python; the only dependency is whatever PETSc the test env happens
to have. The resolve_backend logic is tested by passing an explicit
``available`` list, which avoids any dependency on the host PETSc build.
"""
from __future__ import annotations

import os
from unittest import mock

import pytest

from semi import compute


def test_known_backends_constant():
    assert compute._KNOWN_BACKENDS == ("cpu-mumps", "gpu-amgx", "gpu-hypre")


def test_auto_preference_order():
    assert compute._AUTO_PREFERENCE == ("gpu-amgx", "gpu-hypre", "cpu-mumps")


def test_available_backends_always_includes_cpu_mumps():
    backends = compute.available_backends()
    assert "cpu-mumps" in backends


def test_available_backends_returns_list_subset_of_known():
    backends = compute.available_backends()
    assert isinstance(backends, list)
    for b in backends:
        assert b in compute._KNOWN_BACKENDS


def test_device_info_keys():
    info = compute.device_info()
    expected_keys = {
        "engine_version",
        "petsc_version",
        "petsc_complex",
        "petsc_int64",
        "backends_available",
        "device_count",
        "device_name",
    }
    assert set(info.keys()) == expected_keys


def test_device_info_engine_version_matches_package():
    from semi import __version__
    info = compute.device_info()
    assert info["engine_version"] == __version__


def test_device_info_backends_consistent_with_available():
    info = compute.device_info()
    assert info["backends_available"] == compute.available_backends()


def test_resolve_cpu_mumps_always_succeeds():
    # cpu-mumps is the V&V reference and must resolve regardless of
    # what the available list looks like (it's defensive: even an empty
    # list resolves cpu-mumps because the runner re-validates anyway).
    assert compute.resolve_backend("cpu-mumps", available=["cpu-mumps"]) == "cpu-mumps"
    assert compute.resolve_backend("cpu-mumps", available=[]) == "cpu-mumps"


def test_resolve_gpu_amgx_when_available():
    assert (
        compute.resolve_backend("gpu-amgx", available=["cpu-mumps", "gpu-amgx"])
        == "gpu-amgx"
    )


def test_resolve_gpu_amgx_when_unavailable_raises():
    with pytest.raises(compute.ConfigError, match=r"gpu-amgx"):
        compute.resolve_backend("gpu-amgx", available=["cpu-mumps"])


def test_resolve_gpu_hypre_when_unavailable_raises():
    with pytest.raises(compute.ConfigError, match=r"gpu-hypre"):
        compute.resolve_backend("gpu-hypre", available=["cpu-mumps"])


def test_resolve_auto_prefers_amgx_over_hypre_over_cpu():
    assert (
        compute.resolve_backend(
            "auto", available=["cpu-mumps", "gpu-amgx", "gpu-hypre"]
        )
        == "gpu-amgx"
    )
    assert (
        compute.resolve_backend("auto", available=["cpu-mumps", "gpu-hypre"])
        == "gpu-hypre"
    )
    assert compute.resolve_backend("auto", available=["cpu-mumps"]) == "cpu-mumps"


def test_resolve_auto_with_env_override_cpu():
    with mock.patch.dict(os.environ, {"SEMISIM_BACKEND": "cpu-mumps"}):
        assert (
            compute.resolve_backend(
                "auto", available=["cpu-mumps", "gpu-amgx"]
            )
            == "cpu-mumps"
        )


def test_resolve_auto_with_env_override_unavailable_raises():
    with mock.patch.dict(os.environ, {"SEMISIM_BACKEND": "gpu-amgx"}):
        with pytest.raises(compute.ConfigError, match=r"SEMISIM_BACKEND"):
            compute.resolve_backend("auto", available=["cpu-mumps"])


def test_resolve_auto_with_env_override_unknown_backend():
    with mock.patch.dict(os.environ, {"SEMISIM_BACKEND": "tpu-special"}):
        with pytest.raises(compute.ConfigError, match=r"SEMISIM_BACKEND"):
            compute.resolve_backend("auto", available=["cpu-mumps"])


def test_resolve_unknown_backend_raises():
    with pytest.raises(compute.ConfigError, match=r"recognised"):
        compute.resolve_backend("tpu-special", available=["cpu-mumps"])


def test_resolve_default_available_uses_runtime_probe():
    # Smoke test: call without explicit `available` so it goes through
    # available_backends(); cpu-mumps must resolve.
    assert compute.resolve_backend("cpu-mumps") == "cpu-mumps"


def test_probe_petsc_handles_missing_module():
    """When petsc4py is not importable, _probe_petsc returns the
    conservative no-GPU answer."""
    with mock.patch.dict("sys.modules", {"petsc4py": None}):
        info = compute._probe_petsc()
    assert info["petsc_version"] is None
    assert info["mat_types"] == frozenset()
    assert info["vec_types"] == frozenset()
    assert info["pc_types"] == frozenset()
    assert info["hypre_gpu"] is False


def test_available_backends_with_full_gpu_probe_mock():
    """If the probe reports GPU mat/vec types and AMGX PC, gpu-amgx is
    advertised. If hypre PC is present and pc_hypre_use_gpu is set,
    gpu-hypre is also advertised."""
    fake = {
        "petsc_version": "3.21.0",
        "petsc_complex": False,
        "petsc_int64": False,
        "mat_types": frozenset({"aij", "aijcusparse"}),
        "vec_types": frozenset({"standard", "cuda"}),
        "pc_types": frozenset({"lu", "amgx", "hypre"}),
        "hypre_gpu": True,
    }
    with mock.patch.object(compute, "_probe_petsc", return_value=fake):
        backends = compute.available_backends()
    assert backends == ["cpu-mumps", "gpu-amgx", "gpu-hypre"]


def test_available_backends_hypre_without_gpu_flag_excluded():
    """hypre PC alone is not enough; HYPRE_USE_GPU must be detected."""
    fake = {
        "petsc_version": "3.21.0",
        "petsc_complex": False,
        "petsc_int64": False,
        "mat_types": frozenset({"aij", "aijcusparse"}),
        "vec_types": frozenset({"standard", "cuda"}),
        "pc_types": frozenset({"lu", "hypre"}),
        "hypre_gpu": False,
    }
    with mock.patch.object(compute, "_probe_petsc", return_value=fake):
        backends = compute.available_backends()
    assert "gpu-hypre" not in backends
    assert backends == ["cpu-mumps"]


def test_available_backends_no_gpu_mat_excludes_all_gpu():
    """Without aijcusparse/aijhipsparse, no GPU backend is offered even
    if the PC types are present."""
    fake = {
        "petsc_version": "3.21.0",
        "petsc_complex": False,
        "petsc_int64": False,
        "mat_types": frozenset({"aij"}),
        "vec_types": frozenset({"standard"}),
        "pc_types": frozenset({"lu", "amgx", "hypre"}),
        "hypre_gpu": True,
    }
    with mock.patch.object(compute, "_probe_petsc", return_value=fake):
        backends = compute.available_backends()
    assert backends == ["cpu-mumps"]
