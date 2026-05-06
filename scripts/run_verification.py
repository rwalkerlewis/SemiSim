#!/usr/bin/env python3
"""
V&V (Verification & Validation) runner CLI.

Sister script to `scripts/run_benchmark.py`. Each subcommand runs one
verification activity, writes CSV+PNG artifacts under `results/`, and
exits 0 on success.

Subcommands (Phases 1 and 2 implement `mms_poisson`, `mesh_convergence`,
and `all`):

    python scripts/run_verification.py mms_poisson
    python scripts/run_verification.py mesh_convergence
    python scripts/run_verification.py conservation        # Phase 3 (placeholder)
    python scripts/run_verification.py mms_dd              # Phase 4 (placeholder)
    python scripts/run_verification.py all
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results"


def _gate_finest_pair_rate(rows, key: str, threshold: float, label: str) -> tuple[bool, str]:
    """Pass-fail check on the finest-pair rate column."""
    rate = rows[-1].get(key)
    if rate is None or rate != rate:  # NaN check without numpy
        return False, f"{label}: finest-pair {key} is NaN"
    ok = rate >= threshold
    marker = "PASS" if ok else "FAIL"
    return ok, f"[{marker}] {label}: finest-pair {key} = {rate:.3f} (threshold {threshold:.2f})"


def cmd_mms_poisson(args) -> int:
    """Run all MMS-Poisson studies and gate the finest-pair convergence rates."""
    from semi.verification.mms_poisson import report_table, run_cli_study

    out_dir = Path(args.out) if args.out else RESULTS_DIR / "mms_poisson"
    print(f"[run_verification] mms_poisson -> {out_dir.relative_to(REPO_ROOT)}")
    studies = run_cli_study(out_dir)

    all_ok = True
    for label in ("1d_linear", "1d_nonlinear", "2d_triangles"):
        rows = studies[label]
        print()
        print(report_table(rows, header=f"=== mms_poisson {label} ==="))
        ok_l2, msg_l2 = _gate_finest_pair_rate(rows, "rate_L2", 1.85, f"{label} L2")
        ok_h1, msg_h1 = _gate_finest_pair_rate(rows, "rate_H1", 0.85, f"{label} H1")
        print(msg_l2)
        print(msg_h1)
        all_ok = all_ok and ok_l2 and ok_h1

    smoke = studies["2d_quad_smoke"][0]
    print()
    print("=== mms_poisson 2d_quad_smoke (single mesh, N=64) ===")
    print(
        f"  triangle e_L2 = {smoke['e_L2_triangle']:.3e}, "
        f"quad e_L2 = {smoke['e_L2_quad']:.3e}, "
        f"ratio = {smoke['ratio_L2_quad_over_triangle']:.3f}"
    )
    smoke_ok = 0.1 < smoke["ratio_L2_quad_over_triangle"] < 10.0
    print(("[PASS]" if smoke_ok else "[FAIL]") + " quad/triangle ratio in [0.1, 10]")
    all_ok = all_ok and smoke_ok

    # M6: multi-region coefficient-jump guard on the Si/SiO2 assembly.
    # The M6 acceptance criterion from docs/mos_derivation.md is
    # rate_L2 >= 1.99 on the finest pair (theoretical 2.00). The CLI
    # enforces that stricter floor here; pytest uses a looser 1.85
    # floor to absorb solver-tolerance noise.
    mr_rows = studies["2d_multiregion"]
    print()
    print(report_table(mr_rows, header="=== mms_poisson 2d_multiregion (Si/SiO2) ==="))
    ok_mr_l2, msg_mr_l2 = _gate_finest_pair_rate(
        mr_rows, "rate_L2", 1.99, "2d_multiregion L2",
    )
    ok_mr_h1, msg_mr_h1 = _gate_finest_pair_rate(
        mr_rows, "rate_H1", 0.95, "2d_multiregion H1",
    )
    print(msg_mr_l2)
    print(msg_mr_h1)
    all_ok = all_ok and ok_mr_l2 and ok_mr_h1

    return 0 if all_ok else 5


def cmd_mesh_convergence(args) -> int:
    """Run the pn_1d mesh-convergence sweep and gate monotone error reduction."""
    from semi.verification.mesh_convergence import (
        CLI_NS,
        report_table,
        run_cli_study,
    )

    out_dir = Path(args.out) if args.out else RESULTS_DIR / "mesh_convergence"
    print(f"[run_verification] mesh_convergence -> {out_dir.relative_to(REPO_ROOT)}")
    print(f"[run_verification] sweeping N = {CLI_NS}")
    rows = run_cli_study(out_dir)

    print()
    print(report_table(rows, header="=== mesh_convergence pn_1d ==="))
    print()
    print("Honest flag: the depletion approximation is a model, not the")
    print("truth. The FEM solver converges to the full Poisson-Boltzmann")
    print("solution, so relative errors against the depletion references")
    print("plateau on fine meshes at the physics-model gap, not at zero.")
    print("See semi/verification/mesh_convergence.py docstring.")

    # Acceptance gates.
    # V_bi is set exactly by the Ohmic BCs and is mesh-independent to
    # machine epsilon (not a meaningful convergence indicator), so it is
    # reported only, not gated. E_peak and W against the depletion
    # approximation are gated for strict monotone reduction over the
    # first four refinements (before the physics-model plateau). The
    # Cauchy self-convergence ratio is gated at >= 1.8x per doubling on
    # the same range, which isolates FEM discretization error from the
    # plateau and is the mathematically meaningful convergence rate.
    ok = True

    def _monotone(key, label):
        errs = [r[key] for r in rows[:4]]
        strictly_decreasing = all(
            errs[i + 1] < errs[i] for i in range(len(errs) - 1)
        )
        print(
            f"[{'PASS' if strictly_decreasing else 'FAIL'}] {label} monotone "
            f"on N={CLI_NS[:4]}: {[f'{e:.3e}' for e in errs]}"
        )
        return strictly_decreasing

    def _cauchy_rate(key, label, floor=1.8):
        ratios = [r[key + "_ratio"] for r in rows[1:4]]
        finite = [r for r in ratios if r == r and r > 0.0]
        worst = min(finite) if finite else float("nan")
        passing = finite and worst >= floor
        print(
            f"[{'PASS' if passing else 'FAIL'}] {label} Cauchy rate >= "
            f"{floor:.1f}x per doubling on N={CLI_NS[1:4]}: "
            f"ratios={[f'{r:.2f}' for r in ratios]}"
        )
        return bool(passing)

    print(
        f"[info] V_bi err vs depletion at finest N={CLI_NS[-1]}: "
        f"{rows[-1]['err_Vbi_rel']:.2e} (set exactly by BCs, informational only)"
    )
    ok = _monotone("err_Epeak_rel", "E_peak err vs depletion") and ok
    ok = _monotone("err_W_rel", "W err vs depletion") and ok
    ok = _cauchy_rate("err_Epeak_cauchy", "E_peak") and ok
    ok = _cauchy_rate("err_W_cauchy", "W") and ok

    finest = rows[-1]
    print(
        f"[info] finest level N={finest['N']}: "
        f"solve_time = {finest['solve_time_s']:.3f} s, "
        f"newton_iters = {finest['newton_iters']}"
    )
    return 0 if ok else 6


def cmd_conservation(args) -> int:
    """
    Phase 3: run conservation checks at specific bias points.

    Two benchmarks drive the report:

      - `pn_1d` equilibrium: integrate rho(x) and assert
        |Q_net| < 1e-10 * q * max|N_net| * L_device.
      - `pn_1d_bias` forward (targets 0.3, 0.45, 0.6 V) and
        `pn_1d_bias_reverse` (targets -0.5, -1.0, -2.0 V): at each
        target sample J_total at 10 interior facets and assert
        max |J - mean| / |mean| < 5% (forward) or 15% (reverse).

    The forward and reverse studies each run a single bias sweep (not
    one per target) with a `post_step_hook` that records the interior
    continuity samples on every converged bias, then the report reads
    back the rows at the target biases. That keeps runtime close to
    the production benchmark path.
    """
    import time

    from semi import run as semi_run
    from semi import schema
    from semi.verification.conservation import (
        charge_conservation_from_result,
        run_bias_sweep_with_continuity,
    )

    out_dir = Path(args.out) if args.out else RESULTS_DIR / "conservation"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[run_verification] conservation -> {out_dir.relative_to(REPO_ROOT)}")

    all_ok = True

    # ---- Charge conservation on pn_1d equilibrium ----
    pn_1d_cfg = REPO_ROOT / "benchmarks" / "pn_1d" / "pn_junction.json"
    cfg = schema.load(pn_1d_cfg)
    t0 = time.perf_counter()
    result = semi_run.run(cfg)
    dt_eq = time.perf_counter() - t0
    Q = charge_conservation_from_result(result)
    print()
    print("=== charge conservation (pn_1d equilibrium) ===")
    print(f"  Q_net     = {Q.Q_net:.3e}  C/m^2")
    print(f"  Q_ref     = {Q.Q_ref:.3e}  C/m^2   (= q * max|N_net| * L)")
    print(f"  threshold = {Q.Q_ref * 1.0e-10:.3e}  C/m^2   (1e-10 * Q_ref)")
    print(f"  rel       = {Q.rel_error:.3e}")
    print(f"  solve_time= {dt_eq:.3f} s")
    q_ok = Q.rel_error < 1.0e-10
    print(f"  [{'PASS' if q_ok else 'FAIL'}] |Q_net| < 1e-10 * Q_ref")
    all_ok = all_ok and q_ok

    # ---- Current continuity: forward sweep ----
    print()
    print("=== current continuity: pn_1d_bias (forward, tol 5%) ===")
    fwd_cfg_path = REPO_ROOT / "benchmarks" / "pn_1d_bias" / "pn_junction_bias.json"
    fwd_cfg = schema.load(fwd_cfg_path)
    t0 = time.perf_counter()
    fwd_result = run_bias_sweep_with_continuity(fwd_cfg, n_samples=10)
    dt_fwd = time.perf_counter() - t0
    fwd_targets = [0.30, 0.45, 0.60]
    fwd_ok = _report_continuity_table(fwd_result.iv, fwd_targets, tol=0.05)
    print(f"  sweep_time = {dt_fwd:.3f} s")
    all_ok = all_ok and fwd_ok

    # ---- Current continuity: reverse sweep ----
    print()
    print("=== current continuity: pn_1d_bias_reverse (tol 15%) ===")
    rev_cfg_path = (
        REPO_ROOT / "benchmarks" / "pn_1d_bias_reverse" / "pn_junction_bias_reverse.json"
    )
    rev_cfg = schema.load(rev_cfg_path)
    t0 = time.perf_counter()
    rev_result = run_bias_sweep_with_continuity(rev_cfg, n_samples=10)
    dt_rev = time.perf_counter() - t0
    rev_targets = [-0.50, -1.00, -2.00]
    rev_ok = _report_continuity_table(rev_result.iv, rev_targets, tol=0.15)
    print(f"  sweep_time = {dt_rev:.3f} s")
    all_ok = all_ok and rev_ok

    # ---- CSV artifact ----
    import csv
    csv_path = out_dir / "conservation.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["benchmark", "V", "mean_J", "max_abs_dev", "max_rel_dev",
                    "n_pts", "pass", "tol"])
        for V in fwd_targets:
            row = _closest_iv_row(fwd_result.iv, V)
            w.writerow([
                "pn_1d_bias", row.get("V", ""), row.get("continuity_mean_J", ""),
                row.get("continuity_max_dev", ""), row.get("continuity_max_rel", ""),
                row.get("continuity_n_pts", ""),
                bool(row.get("continuity_max_rel", float("inf")) < 0.05), 0.05,
            ])
        for V in rev_targets:
            row = _closest_iv_row(rev_result.iv, V)
            w.writerow([
                "pn_1d_bias_reverse", row.get("V", ""), row.get("continuity_mean_J", ""),
                row.get("continuity_max_dev", ""), row.get("continuity_max_rel", ""),
                row.get("continuity_n_pts", ""),
                bool(row.get("continuity_max_rel", float("inf")) < 0.15), 0.15,
            ])
        w.writerow(["pn_1d_equilibrium_charge",
                    "", Q.Q_net, "", Q.rel_error, "",
                    Q.rel_error < 1.0e-10, 1.0e-10])
    print()
    print(f"[run_verification] wrote {csv_path.relative_to(REPO_ROOT)}")
    print(f"[run_verification] total runtime: "
          f"{dt_eq + dt_fwd + dt_rev:.1f} s "
          f"(eq {dt_eq:.1f}s, fwd {dt_fwd:.1f}s, rev {dt_rev:.1f}s)")
    return 0 if all_ok else 7


def _closest_iv_row(iv_rows, V_target):
    """Return the iv_row whose V is closest to V_target."""
    if not iv_rows:
        return {}
    return min(iv_rows, key=lambda r: abs(float(r.get("V", 0.0)) - V_target))


def _report_continuity_table(iv_rows, targets, *, tol):
    """Print a pass/fail table at `targets` and return overall pass/fail."""
    print(f"  {'V':>8}  {'mean_J [A/m^2]':>16}  {'max_dev [A/m^2]':>16}  "
          f"{'max_rel':>10}  {'pts':>4}  result")
    all_pass = True
    for V in targets:
        row = _closest_iv_row(iv_rows, V)
        if not row or "continuity_max_rel" not in row:
            print(f"  {V:>8.3f}  {'(no row)':>16}")
            all_pass = False
            continue
        rel = float(row["continuity_max_rel"])
        passed = rel < tol
        all_pass = all_pass and passed
        marker = "PASS" if passed else "FAIL"
        print(
            f"  {float(row['V']):>8.3f}  "
            f"{float(row['continuity_mean_J']):>16.3e}  "
            f"{float(row['continuity_max_dev']):>16.3e}  "
            f"{rel*100:>9.2f}%  "
            f"{int(row['continuity_n_pts']):>4d}  [{marker}]"
        )
    return all_pass


def cmd_mms_dd(args) -> int:
    """
    Phase 4: MMS for the coupled drift-diffusion block residual.

    Runs studies for every variant in {A, B, C, D, E, F, G, H} across
    {1D linear, 1D nonlinear, 2D} and gates the finest-pair L^2 and
    H^1 rates of every meaningful block per variant (psi only for
    Variant A; psi + phi_n + phi_p otherwise). Variants D, E, F, G,
    and H carry the M16.x milestone-acceptance floor L^2 >= 1.99 and
    H^1 >= 0.99 (D: M16.1, E: M16.2, F: M16.3, G: M16.4, H: M16.6).
    Variants A, B, and C stay at the existing 1.75 / 0.80 floor (the
    pytest gate uses the same floor; see docs/mms_dd_derivation.md).

    Variant H asymmetry: the Kane / Hurkx kernels depend only on the
    field magnitude L_0 |grad(psi_hat)|, so the BBT/TAT contribution
    to the (phi_n) / (phi_p) row residuals is decoupled from those
    unknowns. With the small continuity-row weight (L_0^2 mu_hat
    n_hat ~ 1e-18 in scaled units) the phi-block discretization
    error sits near the SNES absolute-residual floor and the textbook
    P1 rates are not recoverable on those blocks. The CLI gate
    therefore restricts Variant H rate gating to the psi block in 1D
    and drops rate gating entirely in 2D (where the integration
    weight L_0^2 ~ 4e-12 lands every block below machine precision
    before Newton can resolve the discretization envelope). This
    matches the pytest gate in tests/fem/test_mms_tunneling.py; the
    Kane and Hurkx physics themselves are independently verified by
    the closed-form NumPy unit tests in tests/test_recombination.py
    and by the zener_1d Kane analytical benchmark.
    """
    from semi.verification.mms_dd import report_table, run_cli_study

    out_dir = Path(args.out) if args.out else RESULTS_DIR / "mms_dd"
    print(f"[run_verification] mms_dd -> {out_dir.relative_to(REPO_ROOT)}")
    studies = run_cli_study(out_dir)

    all_ok = True
    for label, rows in studies.items():
        if "_A_" in label or label.endswith("_A"):
            variant = "A"
        elif "_B_" in label or label.endswith("_B"):
            variant = "B"
        elif "_D_" in label or label.endswith("_D"):
            variant = "D"
        elif "_E_" in label or label.endswith("_E"):
            variant = "E"
        elif "_F_" in label or label.endswith("_F"):
            variant = "F"
        elif "_G_" in label or label.endswith("_G"):
            variant = "G"
        elif "_H_" in label or label.endswith("_H"):
            variant = "H"
        else:
            variant = "C"
        if variant == "A":
            blocks = ("psi",)
        elif variant == "H":
            # 1D: gate only the psi block. 2D: smoke (no rate gate).
            # See docstring above and tests/fem/test_mms_tunneling.py.
            blocks = () if label.startswith("2d_") else ("psi",)
        else:
            blocks = ("psi", "phi_n", "phi_p")
        # Variants D, E, F, G, H are M16.x acceptance gates: tighter floor.
        if variant in ("D", "E", "F", "G", "H"):
            l2_floor = 1.99
            h1_floor = 0.99
        else:
            l2_floor = 1.75
            h1_floor = 0.80
        print()
        print(report_table(rows, header=f"=== mms_dd {label} ==="))
        for b in blocks:
            ok_l2, msg_l2 = _gate_finest_pair_rate(
                rows, f"rate_L2_{b}", l2_floor, f"{label} L2[{b}]",
            )
            ok_h1, msg_h1 = _gate_finest_pair_rate(
                rows, f"rate_H1_{b}", h1_floor, f"{label} H1[{b}]",
            )
            print(msg_l2)
            print(msg_h1)
            all_ok = all_ok and ok_l2 and ok_h1

    return 0 if all_ok else 8


def cmd_all(args) -> int:
    """Run every implemented V&V activity in sequence."""
    rc = cmd_mms_poisson(args)
    if rc != 0:
        return rc
    rc = cmd_mesh_convergence(args)
    if rc != 0:
        return rc
    rc = cmd_conservation(args)
    if rc != 0:
        return rc
    rc = cmd_mms_dd(args)
    if rc != 0:
        return rc
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run a kronos-semi V&V study and write CSV+PNG artifacts."
    )
    parser.add_argument("--out", help="output root directory (default: results/<study>)")
    sub = parser.add_subparsers(dest="study", required=True)

    sub.add_parser("mms_poisson", help="MMS for equilibrium Poisson (Phase 1)")
    sub.add_parser("mesh_convergence", help="pn_1d mesh convergence (Phase 2)")
    sub.add_parser("conservation", help="current and charge conservation (Phase 3)")
    sub.add_parser("mms_dd", help="MMS for coupled drift-diffusion (Phase 4)")
    sub.add_parser("all", help="run every implemented study")

    args = parser.parse_args(argv)
    if args.study == "mms_poisson":
        return cmd_mms_poisson(args)
    if args.study == "mesh_convergence":
        return cmd_mesh_convergence(args)
    if args.study == "conservation":
        return cmd_conservation(args)
    if args.study == "mms_dd":
        return cmd_mms_dd(args)
    if args.study == "all":
        return cmd_all(args)
    parser.error(f"Unknown study {args.study!r}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
