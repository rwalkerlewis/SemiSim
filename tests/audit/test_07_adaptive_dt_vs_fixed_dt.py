"""
Audit case 7 (M18): adaptive-dt run matches the fixed-dt reference
within 1 %.

Loads `benchmarks/pn_1d_turnon` and runs the transient runner twice:

  (a) The shipped fixed-dt config (the v0.24.0 reference). This is
      the baseline I(t) trace.
  (b) The same config with `solver.adaptive.enabled = true,
      dt_min = 1e-13, dt_max = solver.dt`. Constraining dt_max to
      solver.dt isolates the controller's halving response: the
      adaptive run can never exceed the fixed-dt resolution, so any
      drift in I(t) comes from a dt halve in the shape-changing part
      of the transient (or from the variable-step BDF2 path picking
      up a non-unit omega).

The acceptance gate is `|I_adaptive(t) - I_fixed(t)| / |I_fixed(t)|
< 0.01` at every recorded sample. Subsamples the adaptive trace
onto the fixed-dt grid via linear interpolation before comparing,
since the adaptive run may have additional records when dt halves
during the transient. (M18)
"""
from __future__ import annotations

import copy

import numpy as np
import pytest

from ._helpers import (
    load_benchmark,
    require_dolfinx,
    write_csv,
    write_markdown,
)

CASE = "07_adaptive_dt_vs_fixed_dt"


@pytest.mark.audit
def test_adaptive_dt_vs_fixed_dt():
    require_dolfinx()

    from semi.runners.transient import run_transient

    cfg = load_benchmark("pn_1d_turnon/config.json")

    cfg_fixed = copy.deepcopy(cfg)
    res_fixed = run_transient(cfg_fixed)

    cfg_adapt = copy.deepcopy(cfg)
    cfg_adapt["solver"]["adaptive"] = {
        "enabled": True,
        "dt_min": 1.0e-13,
        "dt_max": float(cfg_adapt["solver"]["dt"]),
        "easy_iter_threshold": 4,
        "grow_factor": 1.5,
        "max_consecutive_failures": 6,
    }
    res_adapt = run_transient(cfg_adapt)

    fixed_anode = [r for r in res_fixed.iv if r["contact"] == "anode"]
    adapt_anode = [r for r in res_adapt.iv if r["contact"] == "anode"]
    fixed_anode.sort(key=lambda r: r["t"])
    adapt_anode.sort(key=lambda r: r["t"])

    t_fixed = np.array([r["t"] for r in fixed_anode])
    j_fixed = np.array([r["J"] for r in fixed_anode])
    t_adapt = np.array([r["t"] for r in adapt_anode])
    j_adapt = np.array([r["J"] for r in adapt_anode])

    # Subsample adaptive onto fixed grid via linear interpolation.
    j_adapt_on_fixed = np.interp(t_fixed, t_adapt, j_adapt)

    # Pointwise relative disagreement, denominator-floored against tiny
    # currents so a near-zero forward turn-on tail does not blow up.
    j_ref_scale = max(np.max(np.abs(j_fixed)), 1.0e-30)
    rel_err = np.abs(j_adapt_on_fixed - j_fixed) / np.maximum(
        np.abs(j_fixed), 1.0e-6 * j_ref_scale
    )
    worst = float(np.max(rel_err))
    worst_idx = int(np.argmax(rel_err))
    worst_t = float(t_fixed[worst_idx])
    worst_j_fixed = float(j_fixed[worst_idx])
    worst_j_adapt = float(j_adapt_on_fixed[worst_idx])

    n_steps_fixed = int(res_fixed.meta.get("n_steps_taken", 0))
    n_steps_adapt = int(res_adapt.meta.get("n_steps_taken", 0))
    n_failed_adapt = int(res_adapt.meta.get("n_failed_steps", 0))

    write_csv(
        CASE,
        [
            "t_worst", "J_fixed_worst", "J_adapt_worst",
            "rel_err_worst",
            "n_steps_fixed", "n_steps_adapt", "n_failed_adapt",
            "status",
        ],
        [[
            worst_t, worst_j_fixed, worst_j_adapt,
            worst, n_steps_fixed, n_steps_adapt, n_failed_adapt,
            "passed" if worst < 0.01 else "failed",
        ]],
    )
    write_markdown(
        CASE,
        "Case 07 - adaptive dt vs fixed dt (M18)",
        "On `benchmarks/pn_1d_turnon` the adaptive-dt transient run "
        "(M18; `solver.adaptive.enabled = true` with "
        "`dt_max = solver.dt` so the controller can only halve, not "
        "exceed the fixed-dt resolution) matches the shipped fixed-dt "
        f"reference within the 1 % audit gate. Worst-case relative "
        f"disagreement on I(t) is {worst:.2e} at t = {worst_t:.3e} s "
        f"(see `/tmp/audit/{CASE}.csv`). The adaptive run took "
        f"{n_steps_adapt} accepted steps "
        f"({n_failed_adapt} failed and halved); the fixed-dt run took "
        f"{n_steps_fixed} steps. The variable-step BDF2 path used "
        "when dt halves comes from "
        "`BDFCoefficients.variable_bdf2(omega)` and is bit-identical "
        "to uniform BDF2 at omega = 1.\n\n"
        f"CSV: `/tmp/audit/{CASE}.csv`",
    )

    assert worst < 0.01, (
        f"Worst-case |I_adapt - I_fixed| / |I_fixed| = {worst:.3e} "
        f"exceeds 1% gate at t = {worst_t:.3e} s "
        f"(I_fixed = {worst_j_fixed:.3e}, I_adapt = {worst_j_adapt:.3e})"
    )
