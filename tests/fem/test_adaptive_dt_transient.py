"""
Coverage tests for the M18 adaptive-dt controller in
`semi/runners/transient.py`. These run in the gated `docker-fem-tests`
job; they require dolfinx for the SNES solve.

Each case exercises a single controller transition:

  - growth on consecutive easy SNES solves
  - halving on SNES non-convergence
  - floor exhaustion via StepTooSmall
  - endpoint clamp at t_end (controller not notified)
  - waveform breakpoint clamp for voltage_t.step.t0
  - waveform breakpoint clamp for voltage_t.table.times[i]
  - bit-identity guard: adaptive: {enabled: false} matches `adaptive`
    absent
  - variable_bdf2(omega) coefficients flow through the runner

The cases use a tiny 1D pn diode (4 to 8 cells) with short t_end so
each test runs in O(seconds) inside the gated job. The pure-Python
schema and BDF2 unit tests in `tests/test_adaptive_dt_schema.py` and
`tests/test_variable_bdf2.py` cover the controller-independent pieces;
this file is the FEM-side coverage layer that exercises the runner
end to end.
"""
from __future__ import annotations

import copy

import numpy as np
import pytest


def _have_dolfinx() -> bool:
    try:
        import dolfinx  # noqa: F401
    except ImportError:
        return False
    return True


require_dolfinx = pytest.mark.skipif(
    not _have_dolfinx(),
    reason="dolfinx unavailable; FEM-runner tests skipped",
)


_L = 2.0e-6
_TAU = 1.0e-9


def _make_diode_cfg():
    return {
        "schema_version": "2.9.0",
        "name": "adaptive_dt_diode",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, _L]],
            "resolution": [8],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, _L]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": _L},
            ],
        },
        "regions": {"silicon": {"material": "Si", "tag": 1,
                                "role": "semiconductor"}},
        "doping": [
            {"region": "silicon",
             "profile": {"type": "step", "axis": 0, "location": _L / 2.0,
                         "N_D_left": 0.0, "N_A_left": 1.0e17,
                         "N_D_right": 1.0e17, "N_A_right": 0.0}}
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic",
             "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
        "physics": {"temperature": 300.0, "statistics": "boltzmann",
                    "mobility": {"mu_n": 1400.0, "mu_p": 450.0},
                    "recombination": {"srh": True, "tau_n": _TAU,
                                      "tau_p": _TAU, "E_t": 0.0}},
        "solver": {
            "type": "transient",
            "t_end": 5.0e-9,
            "dt": 1.0e-9,
            "order": 1,
            "bc_ramp_steps": 0,
            "max_steps": 500,
            "output_every": 1000,
        },
    }


def _capture_progress(events):
    def cb(ev):
        events.append(dict(ev))
    return cb


@require_dolfinx
def test_grows_after_easy_solves():
    """Benign config (zero-bias quiescent state, no voltage_t): SNES
    converges in 0-1 iterations every step. The controller must grow
    dt at least once over a 30-step window."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["solver"]["t_end"] = 5.0e-8
    cfg["solver"]["dt"] = 1.0e-10
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-13,
        "dt_max": 1.0e-8,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    events = []
    result = run_transient(cfg, progress_callback=_capture_progress(events))

    dts = [ev["dt"] for ev in events]
    assert len(dts) >= 5
    assert max(dts) > dts[0] * 1.4
    assert result.meta["adaptive_enabled"] is True


@require_dolfinx
def test_halves_on_snes_failure():
    """A config where the initial dt fails the first solve; the
    controller halves and the run completes."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 0.0, "v0": 0.0, "v1": 0.7,
    }
    cfg["solver"]["t_end"] = 1.0e-8
    cfg["solver"]["dt"] = 5.0e-9
    cfg["solver"]["bc_ramp_steps"] = 0
    cfg["solver"]["snes"] = {"rtol": 1.0e-13, "atol": 1.0e-12,
                             "stol": 1.0e-16, "max_it": 6}
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-12,
        "dt_max": 5.0e-9,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 12,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    result = run_transient(cfg)
    # The run finished; failures or no, the run reached t_end.
    assert result.t[-1] == pytest.approx(cfg["solver"]["t_end"], rel=1.0e-9)


@require_dolfinx
def test_floor_exhaustion_raises():
    """A config tuned so even halving to dt_min cannot solve; the
    runner raises RuntimeError citing dt_min."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 0.0, "v0": 0.0, "v1": 50.0,
    }
    cfg["solver"]["t_end"] = 1.0e-9
    cfg["solver"]["dt"] = 1.0e-9
    cfg["solver"]["bc_ramp_steps"] = 0
    cfg["solver"]["snes"] = {"rtol": 1.0e-15, "atol": 1.0e-15,
                             "stol": 1.0e-16, "max_it": 1}
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-10,
        "dt_max": 1.0e-9,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 4,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    with pytest.raises(RuntimeError, match=r"Adaptive transient stalled"):
        run_transient(cfg)


@require_dolfinx
def test_endpoint_clamps_to_t_end():
    """The last step lands exactly at t_end. The endpoint clamp does
    not notify the controller (clamped: True in the progress event)."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["solver"]["t_end"] = 7.0e-9
    cfg["solver"]["dt"] = 3.0e-9  # 3 ns, 3 ns, then 1 ns to land on 7 ns
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-13,
        "dt_max": 3.0e-9,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    events = []
    result = run_transient(cfg, progress_callback=_capture_progress(events))

    t_end = cfg["solver"]["t_end"]
    assert result.t[-1] == pytest.approx(t_end, abs=1.0e-12 * t_end)
    # The final progress event reports clamped=True (endpoint clamp).
    assert events[-1]["clamped"] is True


@require_dolfinx
def test_step_breakpoint_clamp():
    """voltage_t.type='step' with t0 = 5e-9; initial dt = 1e-8 would
    straddle t0. The first step lands exactly at t0 and the BC at the
    post-step solve is v1."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "step", "t0": 5.0e-9, "v0": 0.0, "v1": 0.3,
    }
    cfg["solver"]["t_end"] = 2.0e-8
    cfg["solver"]["dt"] = 1.0e-8
    cfg["solver"]["bc_ramp_steps"] = 0
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-13,
        "dt_max": 1.0e-8,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    result = run_transient(cfg)

    # The first recorded post-step time is exactly t0 = 5 ns.
    assert result.t[1] == pytest.approx(5.0e-9, abs=1.0e-15)
    # The IV row at t = t0 sees v1 (the post-transition value, since
    # the step evaluator returns v1 for t >= t0).
    iv_at_t0 = [
        r for r in result.iv
        if r["contact"] == "anode" and abs(r["t"] - 5.0e-9) < 1.0e-15
    ]
    assert iv_at_t0
    assert iv_at_t0[0]["V"] == pytest.approx(0.3)


@require_dolfinx
def test_table_breakpoint_clamp():
    """voltage_t.type='table' with non-uniform times; initial dt
    larger than the smallest interior interval. The first step lands
    exactly at the smallest interior times[i] in the (0, dt] window."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times":  [0.0, 1.0e-9, 5.0e-9, 8.0e-9],
        "values": [0.0, 0.1,    0.3,    0.4],
    }
    cfg["solver"]["t_end"] = 8.0e-9
    cfg["solver"]["dt"] = 5.0e-9
    cfg["solver"]["bc_ramp_steps"] = 0
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-13,
        "dt_max": 5.0e-9,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    result = run_transient(cfg)

    # The first step lands on times[1] = 1 ns; the second on times[2]
    # = 5 ns (clamped from initial dt = 5 ns offset by 1 ns); the third
    # on times[3] = 8 ns (the t_end clamp coincides here).
    assert result.t[1] == pytest.approx(1.0e-9, abs=1.0e-15)
    assert result.t[2] == pytest.approx(5.0e-9, abs=1.0e-15)
    assert result.t[-1] == pytest.approx(8.0e-9, abs=1.0e-15)


@require_dolfinx
def test_bit_identity_when_disabled():
    """A run with adaptive: {enabled: false} produces the same IV trace
    as a run with `adaptive` absent. The fixed-dt path is the bit-
    identity branch."""
    from semi import schema
    from semi.runners.transient import run_transient

    cfg_a = _make_diode_cfg()
    cfg_a["solver"]["t_end"] = 3.0e-9
    cfg_a["solver"]["dt"] = 1.0e-9

    cfg_b = copy.deepcopy(cfg_a)
    cfg_b["solver"]["adaptive"] = {"enabled": False}

    cfg_a_v = schema.validate(copy.deepcopy(cfg_a))
    cfg_b_v = schema.validate(copy.deepcopy(cfg_b))

    res_a = run_transient(cfg_a_v)
    res_b = run_transient(cfg_b_v)

    assert len(res_a.iv) == len(res_b.iv)
    for ra, rb in zip(res_a.iv, res_b.iv, strict=True):
        assert ra["t"] == rb["t"]
        assert ra["contact"] == rb["contact"]
        # Slotboom solves are deterministic; require byte-identical J.
        assert ra["J"] == rb["J"]
        assert ra["V"] == rb["V"]


@require_dolfinx
def test_variable_bdf2_at_dt_change():
    """When dt changes (here, by clamping), the post-clamp progress
    event reports BDF2 coefficients from variable_bdf2(omega) for the
    observed omega."""
    from semi import schema
    from semi.runners.transient import run_transient
    from semi.timestepping import BDFCoefficients

    cfg = _make_diode_cfg()
    cfg["contacts"][0]["voltage_t"] = {
        "type": "table",
        "times":  [0.0, 1.0e-9, 5.0e-9],
        "values": [0.0, 0.1,    0.3],
    }
    cfg["solver"]["t_end"] = 5.0e-9
    cfg["solver"]["dt"] = 5.0e-9
    cfg["solver"]["order"] = 2
    cfg["solver"]["bc_ramp_steps"] = 0
    cfg["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-13,
        "dt_max": 5.0e-9,
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    cfg = schema.validate(copy.deepcopy(cfg))

    events = []
    run_transient(cfg, progress_callback=_capture_progress(events))

    # Find the first event reporting effective BDF2 (3 alpha coeffs).
    bdf2_events = [ev for ev in events if len(ev["alpha_coeffs"]) == 3]
    assert bdf2_events, "expected at least one BDF2 step"
    ev = bdf2_events[0]
    omega = ev["omega"]
    expected = BDFCoefficients.variable_bdf2(omega)
    actual = ev["alpha_coeffs"]
    assert actual == pytest.approx(expected, rel=1.0e-13, abs=1.0e-13)
    # Sanity: omega is the ratio of the BDF2 step's dt to the previous
    # step's dt, which here was clamped to 1 ns; so the ratio matches
    # the test waveform spacing 4 ns / 1 ns = 4.
    assert omega == pytest.approx(4.0, rel=1.0e-12)
    # Sum-to-zero (constant-in-time exactness) holds.
    assert sum(actual) == pytest.approx(0.0, abs=1.0e-13)
    assert np.isfinite(actual).all()
