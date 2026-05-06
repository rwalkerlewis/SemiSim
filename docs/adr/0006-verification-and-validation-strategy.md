# 0006. Verification & Validation strategy

- Status: Accepted
- Date: 2026-04-21

> Note (2026-05-06, M16.5): the per-physics-module MMS rule below
> applies to *domain* physics (residual additions over `dx`).
> Boundary-physics milestones (Schottky thermionic-emission Robin
> BC, etc.) use analytical-benchmark plus byte-identity gates
> instead of MMS rate gates, because the existing MMS-DD harness
> uses Dirichlet BCs everywhere by construction and substituting a
> Robin BC into the manufactured-source residual either cancels the
> Robin constraint into the source or produces a BC-dominated error
> that does not converge at the optimal rate. See ADR 0015.

## Context

Through M3: Adaptive continuation the project had three benchmark "verifiers" under
`scripts/run_benchmark.py`: `pn_1d` (equilibrium depletion-approximation
check), `pn_1d_bias` (forward Shockley / SNS check), and
`pn_1d_bias_reverse` (SRH generation check). Each runs at a single
mesh resolution against an analytical reference.

In the Roache / Oberkampf-Roy sense these are **validation** exercises
(are the equations we wrote down a reasonable model of a diode?) rather
than **verification** exercises (did we solve those equations
correctly?). They do not prove convergence of the FEM discretization,
do not exercise the solver away from the junction, and do not catch
sign or coefficient errors that happen to cancel at a single
resolution. Several failure modes are invisible to them by
construction:

- A Poisson coefficient off by a factor (the M1 `lambda2` vs.
  `L_D^2` bug went undetected for the first pass because V_bi is set
  by the BCs alone).
- A residual block whose sign is flipped: the coupled Newton may still
  converge to a fixed point that happens to pass the single-point
  check.
- Anisotropic mesh error (e.g. a 2D form that is correct only on
  quads but wrong on triangles).
- Discretization that loses order silently (Galerkin degradation
  under high Peclet, or a P1 form that is accidentally P0-like).

M4: V&V suite was originally scoped as a refactor pass. We re-scoped it to
build a proper V&V suite because the cost of shipping a silently
non-convergent simulator to the evaluation reviewer was higher than
the cost of a refactor delay.

The candidate verification activities considered:

- **Method of Manufactured Solutions (MMS).** A smooth exact
  solution is injected as a forcing term; the code must recover
  theoretical FE convergence rates on a mesh sweep. Catches most
  sign/coefficient bugs and all rate-losing bugs.
- **Mesh convergence studies.** Re-run a physical benchmark at
  increasing resolutions; each derived quantity must converge
  monotonically.
- **Discrete conservation checks.** At each bias, the total
  current J_n + J_p is a constant in space (continuity); the
  space-integral of charge density is zero at equilibrium (charge
  neutrality). Both are exact consequences of the equations and
  so are easy to gate tightly.
- **Code-to-code comparison** against devsim / Sentaurus. Rejected
  for M4: introduces an external-tool dependency and a
  licensing gate, and adds little over MMS for a short-horizon
  project.
- **Experimental validation.** Out of scope; no lab data.

## Decision

Ship a four-phase V&V suite as the M4 deliverable. The suite runs
automatically inside the `docker-fem` CI job via
`scripts/run_verification.py all` and uploads every CSV and PNG
artifact under `results/`.

1. **MMS for equilibrium Poisson** (Phase 1). Module
   `semi/verification/mms_poisson.py`. 1D and 2D smooth
   manufactured solutions (sin and sin-sin with distinct
   wavenumbers). UFL weak-form forcing (do not form
   `ufl.div(ufl.grad(...))` directly; the constant prefactor on
   a coarse mesh can collapse to numerical zero). Acceptance:
   finest-pair L^2 rate >= 1.85, H^1 rate >= 0.85, monotone error
   reduction across the sweep. A `2d_quad_smoke` sanity run
   checks that the quad-element path is within one order of
   magnitude of the triangle path.

2. **Mesh convergence on `pn_1d`** (Phase 2). Module
   `semi/verification/mesh_convergence.py`. Sweep N in
   [50, 100, 200, 400, 800, 1600]; record V_bi, peak |E|, W,
   Newton iterations, solve time. Report error both against the
   depletion approximation (physics-model check) and as Cauchy
   self-convergence ratios between consecutive meshes (pure
   discretization check). Gate monotone error reduction on the
   depletion errors for the first four levels (before the
   physics-model plateau), and gate Cauchy ratios at >= 1.8x per
   doubling over the same range. V_bi is set by the Ohmic BCs
   alone, so we report it but do not gate it.

3. **Conservation checks** (Phase 3). Module
   `semi/verification/conservation.py`. Two gates:
   - Charge conservation on `pn_1d` equilibrium:
     |integrate(q*(p - n + N_D - N_A))| < 1e-10 * q * max|N_net| * L.
   - Current continuity at forward bias targets (0.30, 0.45,
     0.60 V) with `pn_1d_bias` and reverse bias targets (-0.50,
     -1.00, -2.00 V) with `pn_1d_bias_reverse`: sample J_total
     on ten interior facets per target, require
     max |J - mean| / |mean| < 5% forward, < 15% reverse.

4. **MMS for coupled drift-diffusion** (Phase 4). Module
   `semi/verification/mms_dd.py`. Derivation
   `docs/mms_dd_derivation.md`. Three variants exercise
   progressively more of the residual:
   - **Variant A** (psi-only) holds `phi_n_e = phi_p_e = 0` so
     the continuity rows are rate-tested at machine noise and
     only the psi block is gated. Purpose: prove the Poisson
     row of the block residual assembles correctly when driven
     by the coupled nonlinear solver, which the standalone
     Poisson MMS cannot see.
   - **Variant B** exercises the full three-block coupling with
     lifetimes set to 1e+20 so R_e is numerically negligible.
   - **Variant C** runs realistic Si lifetimes (1e-7 s) so the
     SRH coupling between the two continuity blocks is exercised.

   Acceptance: finest-pair L^2 >= 1.75, H^1 >= 0.80 on every
   meaningful block per variant (psi-only on A; psi, phi_n,
   phi_p on B and C). The 0.25 / 0.20 rate headroom below
   theoretical 2 / 1 absorbs block-coupling bleed and the
   exponential SRH numerator.

**CI integration.** The `docker-fem` job adds a step that runs
`python scripts/run_verification.py all` inside the dolfinx image.
The job timeout tightens from 30 to 15 minutes. The existing
`benchmark-plots` artifact is renamed `fem-results` and its path
stays `results/`, so MMS CSVs, mesh-convergence CSVs, conservation
CSVs, and all PNGs are uploaded together.

**Documentation.** A new "Verification & Validation" section in
`docs/PHYSICS.md` records the current finest-pair rates and the
honest flags (depletion-model plateau, Variant A scope, block
residual scale disparity). `docs/mms_dd_derivation.md` stays the
standalone mathematical gate artifact.

## Consequences

Easier:

- Any future change that silently breaks convergence order (sign
  flip, coefficient bug, lost nonlinear term) fails CI before
  reviewer sees it.
- Three distinct verification activities (MMS, mesh convergence,
  conservation) cross-check each other: a bug that hides behind
  one usually surfaces in at least one of the other two.
- The derivation-before-code workflow that produced
  `mms_dd_derivation.md` is now a template. Future V&V extensions
  (2D MOS MMS, 3D resistor convergence) follow the same pattern:
  derivation artifact first, then implementation, then gate.
- The `pn_1d*` benchmarks keep their role as validation exercises
  against analytical physics; V&V and validation are now
  visibly separate activities.

Harder:

- CI wall-clock pressure. The 15-minute budget is tight once the
  dolfinx image is cached; adding a 2D MOS MMS or a 3D resistor
  study will need either per-mesh caching or a split job. This
  is acceptable for the submission timeline but will need
  revisiting in M6-M7.
- Variant A is only gated on psi because the continuity-row
  errors collapse to machine-noise when the quasi-Fermis are
  identically zero. This is documented but could confuse a
  future reader who expects "all three blocks gated on every
  variant." The PHYSICS V&V section and the ADR both call it
  out; the csv rows still report the numeric values.
- The rate thresholds (L^2 >= 1.85 Poisson, >= 1.75 DD; H^1
  >= 0.85 Poisson, >= 0.80 DD) are tighter than standard FE
  theory strictly requires on a single mesh pair. This is
  intentional: it keeps the gate sensitive to half-order losses
  on smooth problems, where there is no mechanism for a real
  rate reduction. Future stabilization work (SUPG / SG) that
  introduces a legitimate half-order loss will need an ADR
  update and a corresponding threshold relaxation.
- The existing `benchmark-plots` CI artifact is renamed to
  `fem-results`. Any downstream tooling or reviewer that fetched
  the old name needs to update; at the time of writing there is
  no such tooling.

## Rejected alternatives

- **Code-to-code comparison against devsim.** Adds an install
  path and a license check without addressing the convergence-order
  question that MMS settles definitively. Deferred to post-submission.
- **Single MMS case covering both Poisson and DD.** Rejected
  because a regression in the DD row would then swamp the Poisson
  signal and vice versa; splitting by physics block keeps the
  signal localized when a gate fails.
- **Only gate L^2 rate, not H^1.** Rejected because a discretization
  that loses one order of H^1 regularity can still pass L^2 by
  luck; H^1 is a strictly stronger gate on smooth problems.
- **Run V&V on every pytest invocation.** Rejected because the
  full sweep is ~1-2 minutes of wall clock; pytest stays fast by
  running a reduced three-level sweep in-process, and the full
  CLI sweep runs only in CI and on demand.

## Related

- `docs/mms_dd_derivation.md` (Phase 4 mathematical gate).
- `docs/PHYSICS.md` Section 5 (current V&V results).
- `semi/verification/{mms_poisson, mesh_convergence, conservation, mms_dd}.py`.
- `scripts/run_verification.py` (CLI entry point).
- `.github/workflows/ci.yml` (CI integration).
- Invariants 3, 4, 5, 6 in `PLAN.md` remain load-bearing across
  the new verification activities.
