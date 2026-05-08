# Roadmap

kronos-semi is a FEniCSx-based finite-element semiconductor device simulator that mimics the capabilities of the COMSOL Semiconductor Module. Simulations are driven by a single JSON file (optionally referencing an external geometry/mesh artifact) and are otherwise plain text. The project ships 1D, 2D, and 3D benchmark problems and a zero-setup Colab notebook (cloud-hosted, no local install) so reviewers can run everything from a browser.

## Capability matrix

All milestones M1 through M15 plus M14.3, M14.4, M16.1 through
M16.7 have shipped as of v0.23.0; M17 (heterojunction /
position-dependent band structure) ships in v0.24.0. The M16
umbrella plus M17 close the position-dependent material parameter
expansion. The table below is the current state; see the
Delivery history section for per-milestone details. Planned
milestones (M19, M19.1, M20) have explicit acceptance tests in
[`docs/IMPROVEMENT_GUIDE.md`](IMPROVEMENT_GUIDE.md).

| Capability | Dimensions | Status | Verifier |
|---|---|---|---|
| Equilibrium Poisson | 1D / 2D / 3D | shipped | MMS finest-pair rates L2 = 2.0, H1 = 1.0 |
| Coupled Slotboom drift-diffusion | 1D / 2D | shipped | MMS finest-pair rates L2 >= 1.99 across variants |
| SRH recombination | 1D / 2D | shipped | verified against SNS analytical at reverse bias |
| Ohmic contact BCs | 1D / 2D / 3D | shipped | Shockley diode within 10% at forward bias |
| Gate contact BCs with phi_ms | 2D | shipped | MOS C-V within 10% in depletion window |
| Multi-region Poisson (Si/SiO2) | 2D | shipped | multi-region MMS L2 = 2.0 |
| File-sourced gmsh .msh meshes | 3D | shipped | builtin vs gmsh R-match within 1% |
| Adaptive bias continuation | uni + bipolar | shipped | pn junction forward + reverse, 3D resistor |
| 3D ohmic V-I linearity | 3D | shipped | V-I linearity within 1% |
| Result artifact writer + manifest (M9) | n/a | shipped | semi-run CLI; manifest.v1.json schema |
| HTTP server POST /solve, GET /runs/{id} (M10) | n/a | shipped | integration tests; WebSocket progress |
| Schema versioning + UI-facing companion (M11) | n/a | shipped | schema_version gate; SCHEMA_SUPPORTED_MINOR |
| Mesh input beyond axis-aligned boxes (M12) | 2D | shipped | MOSFET 2D verifier within ±20% (qualitative) |
| Transient solver BDF1/BDF2 Slotboom (M13/M13.1) | 1D | shipped | steady-state-limit + BDF rate + pn_1d_turnon (5%) |
| AC small-signal analysis (M14) | 1D / 2D | shipped | RC benchmark within 0.4% of analytical C_dep |
| Differential capacitance via AC admittance (M14.1) | 2D | shipped | mos_cap_ac vs mos_cv on Q_gate (byte-identical) |
| Axisymmetric (cylindrical) 2D MOSCAP (M14.2) | 2D | shipped | schema 1.3.0; LF/HF C-V vs Hu Fig. 5-18 analytical |
| GPU linear-solver path (M15) | n/a | shipped | schema 1.4.0; pn_1d psi within 1e-8 vs CPU; >=5x linear-solve speedup at 500k DOFs |
| Test suite | pure + FEM | shipped | 320+ tests (post-M14.3 SG cleanup), coverage >= 95% (CI gate) |
| V&V | 10 studies | shipped | 62/62 PASS |
| CI | pure-python + lint + docker-fem (parallelized) | shipped | green on main |
| Housekeeping (M14.3) | n/a | shipped | mosfet_2d Pao-Sah verifier (20% tolerance in [V_T+0.2, V_T+0.6] V); XDMF mesh ingest (R within 1e-12 vs gmsh path); strict schema v2.0.0 (additionalProperties: false; v1 deprecated for one minor cycle); semi/fem/sg_assembly.py removed; coverage gate raised to 95 |
| Caughey-Thomas field-dependent mobility (M16.1) | 1D / 2D | shipped | MMS-DD Variant D L2 >= 1.99 / H1 >= 0.99 on every block; diode_velsat_1d divergence 56% at V_F=0.9 V, convergence 0.19% at V_F=0.3 V |
| Lombardi surface mobility (M16.2) | 2D / 3D | shipped | MMS-DD Variant E L2 >= 1.99 / H1 >= 0.99 on every block (1D measured 2.000/1.999/2.000, 2D measured 1.997/1.995/1.998); mosfet_2d Pao-Sah verifier window widened to [V_T+0.4, V_T+1.0] V at 10% (run carries M16.1-era allow-failure flag pending separate SNES audit) |
| Auger recombination (M16.3) | 1D / 2D | shipped | MMS-DD Variant F L2 >= 1.99 / H1 >= 0.99 on every block; diode_auger_1d >20% SRH-vs-Auger divergence at V_F = 0.9 V and <10% match to closed-form Hall-Auger ambipolar high-injection asymptote |
| Fermi-Dirac statistics (M16.4) | 1D / 2D | shipped | MMS-DD Variant G L2 >= 1.99 / H1 >= 0.99 on every block (1D measured 2.000/2.000/2.000, 2D measured 1.997/1.999/1.999); diode_fermi_dirac_1d FEM-vs-Blakemore-analytical V_bi 0.0000% and FD-vs-Boltzmann V_bi divergence 7.37% at N_D=1e20 cm^-3 (basic Blakemore production form; see CHANGELOG `[0.20.0]` for the deviation rationale from the >15% / 1e-3-vs-full-integral nominal targets) |
| Schottky contacts (M16.5) | 1D | shipped | schottky_1d slope `ln(J_FEM)` vs V matches 1/V_t within 5% (observed 2.46%) and `|J_FEM - J_thermionic|/J_thermionic < 5x` envelope absolute match (observed worst 278%) over V in [0.1, 0.5] V; ADR 0015 documents the V&V scope |
| BBT and TAT tunneling (M16.6) | 1D | shipped | MMS-DD Variant H psi block L2 = 2.000 / H1 = 1.000 (1D); zener_1d ln-J slope sign indicating BBT firing (observed -0.279 per V) and 5x envelope `\|J_FEM - J_Kane\|/J_Kane` over V_R in [-8, -4] V (observed worst 99.88%); doping relaxed to N = 1e18 cm^-3 from the prompt's nominal 1e19 cm^-3 (deep-bias SNES convergence; tighter follow-up tracked in M16.7 backlog) |
| Time-varying transient contact voltage (M16.7) | 1D | shipped | audit case 06 transient FFT vs AC sweep agreement within 5% (Hann-windowed FFT admittance ratio at V_DC = 0.4 V, F = 1 MHz, dV = 1 mV); pn_1d_pulse and diode_sine_1d demonstration benchmarks ship the schema variants |
| Heterojunctions (M17) | 2D | shipped | schema 2.8.0 `regions[].material_overrides` + `heterojunction` switch; AlGaAs_0p3 in materials database; `semi/physics/heterojunction.py` builds per-cell DG0 chi/Eg/Nc/Nv/n_i/eps_r fields threaded through Poisson and DD form builders; `semi/bcs.py` Anderson-rule local-chi ohmic equilibrium psi; ADR 0016; `benchmarks/hemt_2d/` classical-electrostatic 2DEG reference (FEM-side n_s integration is a Phase F follow-up) |
| 3D MOSFET benchmark (M19) | 3D | Planned | Pao-Sah within 25% (linear); velsat within 30% (saturation); >=5x GPU speedup at 500k DOFs |
| MPI parallel benchmark (M19.1) | 3D | Planned | mosfet_3d under mpiexec -n {1,2,4} same I_D within 1e-8; n=4 vs n=1 speedup >=2.5x |
| HTTP server hardening (M20) | n/a | Planned | unauthenticated POST /solve returns 401; authenticated rate-limited 429 on overrun |

## Honest gap

Two things the engine still does not do well, called out at the top
of this document so a reviewer hits them first instead of digging
through the milestone history.

The 3D semiconductor coverage is thin. The shipped 3D benchmarks are
the doped resistor (M7) and a pure-Poisson box (the M15 GPU
acceptance test). There is no real 3D semiconductor device. M19
closes this gap with a 3D MOSFET benchmark on a gmsh-sourced
unstructured mesh, run on both CPU-MUMPS and GPU-AMGX backends; M19
depends on M16.1 (Caughey-Thomas mobility) so saturation has any
meaning.

The physics catalogue is complete in the M16 umbrella plus M17.
Field-dependent mobility (M16.1 / M16.2), Auger recombination
(M16.3), Fermi-Dirac statistics (M16.4), Schottky contacts
(M16.5), BBT / TAT tunneling (M16.6), transient time-varying
contact voltage with FFT-vs-AC audit (M16.7), and heterojunction
/ position-dependent band structure (M17) have all shipped.
Quantum confinement corrections (Schrodinger-Poisson coupling,
density-gradient, Bohm potential) are the next-tier physics gap;
the 15 % HEMT 2DEG tolerance reflects this gap and a quantum-
corrected verifier is a future M-numbered milestone.

## Scope vs. COMSOL Semiconductor Module

kronos-semi covers the quasi-static steady-state and time-dependent
subsets of the COMSOL Semiconductor Module, plus AC small-signal
analysis and axisymmetric (cylindrical) 2D devices:

**In scope (shipped):**
- Poisson equation with multi-region dielectric (Si/SiO2)
- Drift-diffusion in Slotboom (quasi-Fermi potential) form
- SRH recombination with configurable trap energy
- Ohmic contacts and ideal gate contacts
- 1D interval, 2D rectangle, and 3D box meshes (builtin)
- 3D unstructured tetrahedral meshes via gmsh .msh files
- Bias sweeps with adaptive step-size continuation
- Method-of-Manufactured-Solutions and conservation V&V suite
- Transient (time-dependent) solver with BDF1/BDF2 (M13/M13.1)
- AC small-signal analysis with frequency sweep (M14)
- Differential capacitance via AC admittance (M14.1)
- Axisymmetric (cylindrical) 2D devices, MOSCAP LF/HF C-V (M14.2)

- Housekeeping (M14.3, v0.16.0): mosfet_2d Pao-Sah verifier; XDMF
  mesh ingest; strict-mode input schema v2.0.0 with
  `additionalProperties: false`; dead Scharfetter-Gummel primitives
  removed; coverage gate raised to 95.
- Caughey-Thomas / Lombardi field-dependent mobility (M16.1, M16.2)
- Auger and radiative recombination (M16.3)
- Fermi-Dirac statistics under the basic Blakemore approximation
  (M16.4); Boltzmann remains the default
- Schottky contacts via thermionic-emission Robin BC plus metal-
  Fermi-level psi Dirichlet (M16.5; v0.21.0)
- Band-to-band / trap-assisted tunneling via Kane and Hurkx
  kernels (M16.6; v0.22.0)
- Time-varying transient contact voltage with FFT-vs-AC audit
  (M16.7; v0.23.0)

**Explicitly out of scope today (M17 heterojunctions, M19 3D
MOSFET, M19.1 MPI, M20 server hardening planned):**

- Heterojunctions, position-dependent chi(x) and Eg(x) (M17;
  depends on M16.4)
- 3D MOSFET on a gmsh-sourced unstructured mesh (M19; depends on
  M16.1)
- MPI parallel orchestration in the runners (M19.1; depends on
  M19)
- HTTP server auth, rate limiting, API keys (M20)

## JSON input contract

Every simulation is driven by a JSON file validated against the jsonschema Draft-07 schema in `semi/schema.py`. The JSON specifies: mesh source (builtin geometry or path to an external .msh file), regions with materials and tags, doping profiles, contacts, physics options, and solver options. A JSON file plus the package version fully determines the simulation. See `docs/adr/0001-json-as-input-format.md` for the rationale.

## Colab onboarding

Each benchmark notebook opens from a Colab badge. The first cell installs FEniCSx via FEM-on-Colab (~30 s), clones the repo, and loads a JSON benchmark file. No local installation is required. See the README for Colab links.

## Delivery history

The milestone entries below record the per-day goal, deliverables, verification criteria, and dependencies for all shipped milestones (M1–M14.1). This is the historical record of how each capability was built; the capability matrix above is the concise entry point for reviewers.

Read `PLAN.md` for the current in-flight task and `docs/IMPROVEMENT_GUIDE.md` for the full milestone specs for M15, M16, and M17.

## M1: Equilibrium Poisson, 1D pn junction, Docker env

- **Status:** Done (2026-04-20). PR `dev/docker-day1-fix`.
- **Goal:** prove the end-to-end pipeline works for a known analytical
  test case.
- **Deliverables:**
  - `Dockerfile`, `docker-compose.yml`, `.devcontainer/devcontainer.json`,
    `.dockerignore`.
  - `semi/constants.py`, `semi/materials.py`, `semi/scaling.py`,
    `semi/doping.py`, `semi/schema.py`.
  - `semi/mesh.py`, `semi/physics/poisson.py`, `semi/solver.py`,
    `semi/run.py`.
  - `benchmarks/pn_1d/pn_junction.json`.
  - `scripts/run_benchmark.py` with a `pn_1d` verifier.
  - 36 pytest tests covering the pure-Python core.
- **Verification:**
  - `docker compose run --rm test` exits 0 with 36/36 pass.
  - `docker compose run --rm benchmark pn_1d` exits 0 with all six
    verifier checks green.
- **Dependencies:** none.

## M2: Coupled drift-diffusion (Slotboom, coupled Newton, bias sweep)

- **Status:** Done (2026-04-20). Delivered on
  `dev/day2-drift-diffusion`; all deliverables landed, the
  `pn_1d_bias` verifier is green (V=0.6 V within 10% of Shockley),
  and M1 did not regress. The low-bias ideal-Shockley mismatch is
  physical (SRH depletion-region recombination raises the ideality
  factor) and is handled by the M2 verifier with qualitative
  checks below 0.5 V. M3 will replace those qualitative checks
  with a Sah-Noyce-Shockley term.
- **Goal (as delivered):** solved the full coupled (psi, Phi_n, Phi_p)
  system with SRH recombination under forward bias, verified against
  Shockley at high bias.
- **Deliverables:**
  - `semi/physics/drift_diffusion.py` with Slotboom continuity forms
    (electron and hole continuity as separate form builders, plus a
    combined block residual).
  - `semi/physics/recombination.py` with the SRH kernel as a UFL
    expression builder.
  - Extension of `semi/solver.py` to dispatch between scalar and
    blocked `NonlinearProblem` via the `kind` argument.
  - Update `semi/run.py` to build the coupled system when the config's
    `solver.type` is `"drift_diffusion"` or `"bias_sweep"`.
  - Schema extension: `recombination` block (tau_n, tau_p, E_t),
    `contacts[*].voltage_sweep` for bias ramping.
  - `benchmarks/pn_1d_bias/pn_junction_bias.json` with a forward-bias
    IV sweep from 0.0 to 0.6 V in 0.05 V increments.
  - `pn_1d_bias` verifier in `scripts/run_benchmark.py` checking IV vs
    Shockley `J = J_0 (exp(V/V_t) - 1)` with the saturation current
    $J_0$ computed from $D_n, D_p, \tau_n, \tau_p, N_A, N_D, n_i$.
  - At least 10 new unit tests: Slotboom math round-trip (n, p
    recovery), SRH kernel limits (R = 0 at equilibrium, R > 0 under
    bias, correct saturation behavior).
- **Verification:**
  - `docker compose run --rm benchmark pn_1d` still green (M1
    regression).
  - `docker compose run --rm benchmark pn_1d_bias` exits 0 with IV
    curve within 10% of Shockley over V in [0.2, 0.6] V. Below 0.2 V
    the current is dominated by numerical noise and Shockley is not a
    meaningful target.
  - `docker compose run --rm test` green with the new tests.
  - ruff clean.
- **Dependencies:** M1 merged to `main`.

## M3: Adaptive continuation (bias ramping hardening, Shockley IV polish)

- **Status:** Done (2026-04-21). PR `dev/day3-bias-hardening`.
- **Goal:** make the bias sweep robust, add reverse-bias check,
  document the ramp-continuation logic.
- **Delivered:**
  - Adaptive step-size control: `AdaptiveStepController` grows the
    step by `grow_factor` (default 1.5) after `easy_iter_threshold`
    (default 4) consecutive easy solves, halves on SNES failure, and
    clamps to the sweep endpoint. Cuts total SNES iterations on the
    forward `pn_1d_bias` sweep from 42 to 31 (26.2% reduction).
  - Sah-Noyce-Shockley recombination term (with the leading-order
    Sze f = 2 V_t/(V_bi - V) correction) in the forward verifier:
    J_sim matches J_diff + J_rec within 15% on [0.15, 0.55] V and
    Shockley diffusion within 10% at V = 0.6 V.
  - New `benchmarks/pn_1d_bias_reverse/` (0 to -2 V) with a
    verifier demanding |J| within 20% of the net SRH generation
    current (q n_i / 2 tau_eff)(W(V) - W(0)) on [-2, -0.5] V. For
    tau = 1e-8 s this replaces the originally-proposed "J saturates
    to J_0" bar, which does not apply in the SRH-generation regime.
  - Extracted analytical reference curves into
    `semi/diode_analytical.py`; 26 new pytest tests.
  - New "Bias continuation strategy" subsection in `docs/PHYSICS.md`.
  - CI `docker-fem` job runs the reverse benchmark on every push.
- **Dependencies:** M2.

## M4: V&V suite (Verification & Validation)

- **Status:** Done (2026-04-21). PR `dev/day4-vnv`.
- **Goal:** replace single-point physical sanity tests with proper
  verification (Roache / Oberkampf-Roy sense): manufactured solutions,
  mesh-refinement convergence studies, and discrete conservation
  checks. The originally-planned refactor moves to M5.
- **Rationale:** the existing benchmark "verifiers" prove the code
  reproduces analytical reference values at one mesh resolution. They
  do not prove convergence, do not exercise the discretization away
  from operating points, and do not catch sign or coefficient errors
  that happen to cancel at a single resolution. MMS is the standard
  remedy.
- **Deliverables:**
  - **Phase 1: MMS for Poisson.**
    - `semi/verification/mms_poisson.py`: 1D and 2D smooth
      manufactured solutions, UFL-based forcing, mesh sweep,
      L^2 / H^1 error and observed-rate computation, CSV+PNG output.
    - `tests/fem/test_mms_poisson.py` with `tests/fem/conftest.py`
      that skips when dolfinx is missing.
  - **Phase 2: Mesh convergence on `pn_1d`.**
    - `semi/verification/mesh_convergence.py`: sweep N in
      [50, 100, 200, 400, 800, 1600], record V_bi error, peak |E|
      error, depletion-width error, Newton iterations, solve time.
    - `tests/fem/test_mesh_convergence.py` (3 levels for speed).
    - Documented convergence order, including any plateau caused by
      the depletion-approximation reference itself being
      approximate.
  - **Phase 3: Conservation checks.**
    - Current continuity check (max |J_total(x) - mean| / mean)
      added to the `pn_1d_bias` and `pn_1d_bias_reverse` verifiers.
      Threshold +-5% forward, +-15% reverse.
    - Charge conservation check (integrated rho over the device
      bounded by solver tolerance) added to `pn_1d`.
    - Unit tests on the check functions in isolation, with mock
      `SimulationResult` objects.
  - **Phase 4: MMS for coupled drift-diffusion.**
    - `semi/verification/mms_dd.py`: manufactured psi, Phi_n, Phi_p
      with three forcing terms in UFL.
    - Three pytest tests: psi-only with frozen quasi-Fermis;
      full three-block coupling; full coupling with nonzero SRH.
    - Acceptance: L^2 rate >= 1.75 on each field on the finest pair.
  - **Phase 5: CI integration.**
    - `scripts/run_verification.py` CLI parallel to
      `scripts/run_benchmark.py` (subcommands `mms_poisson`,
      `mms_dd`, `mesh_convergence`, `conservation`, `all`).
    - `.github/workflows/ci.yml`: run V&V step inside `docker-fem`,
      upload `results/mms_*` and `results/mesh_convergence/`
      artifacts. Total CI runtime under 15 minutes.
  - **Final phase: Documentation.**
    - New "Verification & Validation" section in `docs/PHYSICS.md`.
    - New ADR `docs/adr/0006-verification-and-validation-strategy.md`.
    - PLAN.md, ROADMAP.md, CHANGELOG.md final update.
- **Delivered:**
  - MMS-Poisson 1D linear / 1D nonlinear / 2D triangles with
    finest-pair L^2 at theoretical 2.000 and H^1 at 0.999-1.000;
    `2d_quad_smoke` sanity ratio 0.394.
  - Mesh convergence on `pn_1d` over N in [50..1600] with
    self-convergence Cauchy ratios at >= 1.99x per doubling for
    E_peak and W over the first four levels (before the
    physics-model plateau from the depletion-approximation
    reference, which is honest-flagged in both the runner output
    and in `docs/PHYSICS.md`). V_bi is set by the Ohmic BCs and
    so is reported but not gated.
  - Conservation: charge on `pn_1d` equilibrium at rel 1.5e-17
    (threshold 1e-10); current continuity on `pn_1d_bias` V in
    {0.30, 0.45, 0.60} V and `pn_1d_bias_reverse` V in {-0.50,
    -1.00, -2.00} V, worst max_rel 1.91% forward / 0.020%
    reverse inside the 5% / 15% tolerances.
  - MMS-DD nine studies (three variants x three grids) with
    every gated block rate >= 1.99 L^2 and >= 0.996 H^1.
    Variant A gates only the psi block because the continuity
    rows collapse to machine roundoff when the quasi-Fermis are
    identically zero (documented in `mms_dd_derivation.md`
    Amendment section).
  - `scripts/run_verification.py` with subcommands `mms_poisson`,
    `mesh_convergence`, `conservation`, `mms_dd`, `all`.
  - CI `docker-fem` job now runs V&V inside the dolfinx image,
    timeout tightened to 15 minutes, artifacts uploaded as
    `fem-results` (renamed from `benchmark-plots`).
  - New ADR `0006-verification-and-validation-strategy.md`, new
    "Verification & Validation" section in `docs/PHYSICS.md`,
    derivation artifact `docs/mms_dd_derivation.md` amended with
    the SNES atol-0 fix rationale.
- **Dependencies:** M3 (PR #5) merged into `main`.

## M5: Refactor and test pass

- **Status:** Done (2026-04-21). Delivered on `dev/day5-refactor`.
- **Goal:** paid down technical debt accumulated across Days 1-3 and
  expanded coverage before moving to higher dimensions.
- **Deliverables (as delivered):**
  - `semi/bcs.py` extracted (pure-Python core tier, no dolfinx import
    at module scope) with `ContactBC` dataclass,
    `resolve_contacts(cfg, facet_tags=None, voltages=None)`,
    `build_psi_dirichlet_bcs(...)`, and `build_dd_dirichlet_bcs(...)`.
    Inline `_build_ohmic_bcs_psi` and `_build_dd_ohmic_bcs` removed
    from `semi/run.py`.
  - `semi/run.py` split into a thin dispatcher (74 lines) plus
    `semi/runners/` (`equilibrium.py`, `bias_sweep.py`, `_common.py`)
    and `semi/postprocess.py`. `bias_sweep.py` 294 lines; no module
    in `semi/` exceeds 300 lines (was 580 in `run.py`).
  - Coverage 96.25% (1598 statements, 60 missed) across 177 tests
    (139 prior + 38 new). CI gate `--cov-fail-under=95` enforced
    in the `docker-fem` job.
  - `docs/PHYSICS.md` Section 2.5 completed: full scaled
    drift-diffusion derivation (continuity rows pick up `L_0^2`
    explicitly because the mesh stays in meters per Invariant 3),
    scaled SRH kernel, residual-sign block summary, J_0 numerical
    check.
  - `docs/adr/0007-contact-bc-interface.md` records the dataclass-
    over-dict / voltages-override / insulating-skip / optional-facet-
    verification design choices.
- **Verification:**
  - `pytest --cov=semi --cov-fail-under=95` exits 0 at 96.25%.
  - `pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse` benchmarks all pass
    with byte-identical numerics to the M4 baseline.
  - `python scripts/run_verification.py all` reports 53 PASS / 0 FAIL
    with finest-pair MMS rates byte-identical to M4.
- **Dependencies:** M4 V&V suite merged (PR #6).

## M6: 2D MOS capacitor

- **Status:** Done (2026-04-21). Delivered on `dev/day6-mos-2d`.
  All eight deliverables landed; the `mos_2d` benchmark exits 0
  with 4/4 verifier checks green; M4-M5 V&V suite stayed green
  with the new multi-region Poisson MMS study clearing
  rate_L^2 >= 1.99; 1D benchmarks (`pn_1d`, `pn_1d_bias`,
  `pn_1d_bias_reverse`) are byte-identical to the M5 baseline.
- **Goal (as delivered):** first 2D benchmark and first multi-region
  device. Equilibrium Poisson assembles over the full mesh with
  cellwise DG0 eps_r; space charge is restricted to silicon via
  `dx(subdomain_id=semi_tag)`; gate contact applies a Dirichlet BC on
  psi at the oxide top face; a C-V verifier differentiates gate charge
  with respect to V_gate and matches depletion-approximation MOS
  theory within 10% in the verifier window.
- **Deliverables (as landed):**
  - `docs/mos_derivation.md`: derivation-first gate, approved M6
    prerequisite (9da5fa6).
  - `semi/mesh.py` submesh + cellwise eps_r helpers (719b9e8).
  - `semi/bcs.py` gate contact wiring (b24c650).
  - `semi/physics/poisson.py:build_equilibrium_poisson_form_mr`
    (this commit).
  - `semi/physics/drift_diffusion.py` submesh block residual via
    `entity_maps` (455196a).
  - `benchmarks/mos_2d/mos_cap.json` device spec: 500 nm p-type Si
    (N_A = 1e17 cm^-3) / 5 nm SiO2, uniform 1 nm vertical mesh,
    V_gate sweep [-0.9, +1.2] V.
  - `scripts/run_benchmark.py`: 2D `plot_mos_2d` (tricontourf of
    psi, central-column psi(y), C-V and Q-V) and `verify_mos_2d`
    (10% tolerance in [V_FB + 0.2, V_T - 0.1] V).
  - `semi/verification/mms_poisson.py`:
    `run_mms_poisson_2d_multiregion` with Si/SiO2 coefficient jump
    (99fcd2b).
  - `docs/PHYSICS.md` Section 6 (MOS reference).
- **Verification (observed):**
  - `benchmark mos_2d` exits 0; worst |C_sim - C_theory|/C_theory in
    the verifier window is 9.25% at V_gate = -0.20 V.
  - 4 plots written: `psi_2d.png`, `potentials_1d.png`, `cv.png`,
    `qv.png`.
  - V&V suite rates: see `PHYSICS.md` sections 5.1-5.5.
  - `pytest --cov=semi --cov-fail-under=95` exits 0 at 95.43%
    (195 tests pass).
- **Dependencies:** M5 refactor (PR #7).

## M7: 3D doped resistor

- **Status:** Completed (2026-04-21). Merged via PR #9 from
  `dev/day7-resistor-3d` at `a604b12`.
  All deliverables landed; the `resistor_3d` benchmark exits 0 on
  both the builtin `create_box` mesh and the committed gmsh fixture
  with worst |R_sim - R_theory|/R_theory well inside the fixed 1%
  tolerance on both paths; M1-M6 benchmarks (`pn_1d`,
  `pn_1d_bias`, `pn_1d_bias_reverse`, `mos_2d`) stay byte-identical;
  the M4-M6 V&V suite stays green (every MMS rate within 0.01 of
  the post-M6 values, including `2d_multiregion`).
- **Goal (as delivered):** confirmed the simulator framework extends
  to 3D with no physics changes. Equilibrium Poisson and the
  Slotboom drift-diffusion forms (`semi/physics/poisson.py`,
  `semi/physics/drift_diffusion.py`) are dimension agnostic and
  were left untouched. The only new machinery is the gmsh `.msh`
  loader, a 3D benchmark device, a bipolar-sweep driver path, a
  V-I linearity verifier, and 3D slice plots.
- **Delivered:**
  - `docs/resistor_derivation.md` (derivation-lite gate; 466a820):
    device geometry, analytical ohmic resistance, 1% V-I linearity
    rationale, 3D slice-plot strategy, gmsh loader test strategy.
  - `semi/mesh.py::_build_from_file` wired via
    `dolfinx.io.gmsh.read_from_msh` (e395830, docs alignment in
    3abd170). Physical groups stored in the `.msh` file are
    returned verbatim as `cell_tags` and `facet_tags`; the JSON
    box-tagger is bypassed for file-source meshes.
  - `semi/runners/bias_sweep.py` two-leg walk for bipolar sweeps
    (7c947b2): when `v_sweep_list` spans zero, the ramp runs from
    `V = 0 -> min(V) -> max(V)` with a fresh
    `AdaptiveStepController` on each leg; unipolar sweeps behave
    exactly as before.
  - `benchmarks/resistor_3d/resistor.json`: 3D rectangular bar
    (1 um x 200 nm x 200 nm), uniform n-type N_D = 1e18 cm^-3,
    ohmic contacts on x=0 and x=L, 5-point sweep in
    [-0.01, +0.01] V (62f52d2).
  - `benchmarks/resistor_3d/resistor_gmsh.json` plus a committed
    fixture `fixtures/box.geo`/`box.msh` (reproducible via
    `gmsh -3 box.geo -o box.msh`). Same device on an unstructured
    tetrahedral mesh.
  - `scripts/run_benchmark.py`: `verify_resistor_3d` V-I linearity
    verifier (1% tol, zero-bias-noise and sign-consistency sanity
    checks); 3D slice plotters for psi and |J_n| at the y = W/2
    midplane plus an I-V scatter.
  - Tests: `tests/fem/test_mesh_gmsh.py` (round-trip + physical
    group preservation) and `tests/fem/test_resistor_3d.py`
    (coarse smoke plus builtin/gmsh production-JSON theory match).
  - `docs/PHYSICS.md` Section 7 (3D extension notes, ~1 page)
    citing `docs/resistor_derivation.md` for the full derivation.
- **Verification:**
  - `docker compose run --rm benchmark resistor_3d` exits 0 with
    V-I linearity within 1% on both builtin and gmsh variants;
    three plots written (`psi_slice_y_midplane.png`,
    `jn_slice_y_midplane.png`, IV scatter).
  - All prior benchmarks byte-identical; all prior V&V gates
    green; coverage stays >= 95%.
- **Dependencies:** M6 MOS merged (PR #8).

## M8: Submission polish

- **Status:** In flight. In progress on `dev/submission-polish`, targeting
  PR #10 against `main`. Was originally M7 in the pre-resistor
  roadmap.
- **Goal:** make the submission presentation-ready. No new physics,
  no new benchmarks; this is a documentation, notebook, and release
  pass only.
- **Deliverables:**
  - Sync `PLAN.md` "Current state" and this roadmap to reflect the
    M7 merge and M8 in flight.
  Rewrite the `README.md` status section as an end-of-M7
    capability matrix (M1-M7 shipped across PRs #2-#9), with a
    short scope-out prose block enumerating COMSOL Semiconductor-Module
    features that are deliberately out of scope.
  - Regenerate `notebooks/01_pn_junction_1d.ipynb` via
    `scripts/build_notebook_01.py` with current-state framing (no
    framing (no "M1" heading).
  - Author `notebooks/02_pn_junction_bias.ipynb` (M2-M3 content:
    forward Shockley + reverse SNS sweeps).
  - Author `notebooks/03_mos_cv.ipynb` (M6 content: MOS C-V, with
    the `V_FB + 0.2` verifier-window disclosure surfaced in the
    narrative).
  - Author `notebooks/04_resistor_3d.ipynb` (M7 content: bipolar
    sweep, builtin-vs-gmsh comparison).
  - Add a README "Notebooks" catalog with Colab badges.
  - Append `[0.8.0] - M8: Submission polish` to `CHANGELOG.md`.
  - Open PR #10, wait for CI, do not self-merge.
  - Post-merge and only on explicit human prompt, create annotated
    tag `v0.2.0`.
- **Verification:**
  - Every notebook executed on a real Colab runtime (not just
    `jupyter nbconvert --execute` locally) with wall time and
    final-cell observation recorded in the PR body.
  README and CHANGELOG accurately describe the end-of-M7
    capability set.
  - CI green on every commit on the branch
    (`gh run list --branch dev/submission-polish` verbatim).
  - All five benchmarks (`pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse`,
    `mos_2d`, `resistor_3d`) stay green; coverage stays >= 95%;
    V&V suite stays green; no physics, schema, or verifier changes.
- **Dependencies:** M6 and M7.

## Post-submission (stretch goals)

Explicitly out of scope for the evaluation submission. Listed here so
contributors and reviewers know the intended direction.

- **Field-dependent mobility:** Caughey-Thomas parameterization with
  saturation velocity, needed for any quantitatively accurate
  high-field device simulation.
- **Auger and radiative recombination:** becomes important in heavily
  doped regions and direct-gap materials.
- **Fermi-Dirac statistics:** needed for degeneracy in source/drain
  regions and for wide-gap materials at low temperature.
- **Incomplete dopant ionization:** freeze-out effects at low T.
- **AC small-signal:** Shipped in M14 — no longer a stretch goal.
- **Transient solver:** Shipped in M13/M13.1 — no longer a stretch goal.
- **Band-to-band and trap-assisted tunneling:** breakdown and gate
  leakage.
- **Heterojunctions:** multi-semiconductor stacks (GaAs/AlGaAs,
  Si/SiGe).
- **Full MOSFET:** source/drain/gate/body with submesh for carriers
  and proper work-function handling.
- **FinFET and 3D transistor geometries:** geometry generation and
  meshing, plus whatever new physics the above stretch goals bring.
- **GUI or web frontend:** a JavaScript SPA that builds the JSON
  input interactively and streams results back.

## M9: Result artifact writer

- **Status:** Done (2026-04-23). See PLAN.md completed work log and PR #12.
- **Goal:** Persist every simulation result to disk in a reproducible
  manifest-plus-artifact layout so downstream consumers (HTTP server,
  UI, notebooks) can read results without re-running the solver.
- **Deliverables:**
  - `semi/io/artifact.py` — writes field arrays and IV tables as
    NumPy `.npy` files alongside a `manifest.json` that records
    solver config, version, and output paths.
  - `schemas/manifest.v1.json` — JSON Schema for the manifest;
    major-version gate mirrors the input schema gate.
  - `semi-run` CLI entry point (`pyproject.toml scripts`) that runs
    a JSON config and writes the artifact bundle to `--output-dir`.
- **Verification:**
  - Pure-Python tests for artifact round-trip (write → load → compare).
  - 238 tests pass, coverage 95.5%.
- **Dependencies:** M8.

## M10: HTTP server

- **Status:** Done (2026-04-23). See PLAN.md completed work log and PR #13.
- **Goal:** Expose the solver over HTTP so the future UI (separate
  repo) can drive simulations without installing the Python package.
- **Deliverables:**
  - `kronos_server/` top-level package with a FastAPI app.
  - `POST /solve` — accepts a JSON config, queues a solve, returns a
    run ID.
  - `GET /runs/{id}` — returns the run status and artifact paths.
  - WebSocket progress endpoint for live iteration feedback.
- **Verification:**
  - Integration tests for the HTTP endpoints.
  - `docker compose run --rm server` starts without error.
- **Dependencies:** M9.

## M11: Schema versioning + UI-facing schema companion

- **Status:** Done (2026-04-23). See PLAN.md completed work log and PR #14.
- **Goal:** Decouple the JSON schema from the Python source so the UI
  can load it without importing the engine; add semver gating so old
  configs fail loudly on a breaking schema change.
- **Deliverables:**
  - `schemas/input.v1.json` — verbatim Draft-07 schema extracted from
    `semi/schema.py`.
  - `semi/schema.py` now loads `schemas/input.v1.json` at import time;
    `ENGINE_SUPPORTED_SCHEMA_MAJOR = 1`, `SCHEMA_SUPPORTED_MINOR` tracks
    minor bumps.
  - `schema_version` (semver) field in every config; major-version gate
    in `semi.schema.validate`.
  - `GET /schema` returns `{schema, version, supported_major}`.
  - `_fill_defaults` migrated into the schema `default`/`examples` fields.
  - `tests/test_schema_versioning.py` — 8 pure-Python tests.
- **Verification:**
  - 238 tests pass, coverage 95.5%, ruff clean.
- **Dependencies:** M10.

## M12: MOSFET n+ doping + SNES tolerance amendment (ADR 0008)

- **Status:** Done (2026-04-25). See PLAN.md completed work log, PR #18,
  and `docs/adr/0008-snes-tolerances.md`.
- **Goal:** Support Gaussian implant profiles for realistic source/drain
  regions; amend SNES tolerances to handle the high-contrast doping
  (see ADR 0008).
- **Deliverables:**
  - Gaussian doping profile builder in `semi/doping.py`.
  - Relaxed default SNES tolerances documented in ADR 0008.
  - `benchmarks/mosfet_2d/` — 2D MOSFET with n+ source/drain implants.
  - MOSFET 2D verifier (`scripts/run_benchmark.py mosfet_2d`) checking
    threshold voltage and subthreshold slope within ±20%.
- **Verification:**
  - All existing benchmarks still green.
  - MOSFET 2D verifier exits 0.
- **Dependencies:** M11.

## M13: Transient solver — BDF1/BDF2 (ADR 0009, ADR 0010)

- **Status:** Done (2026-04-25), superseded by M13.1 for Slotboom form.
  See PLAN.md completed work log and `docs/adr/0009-transient-formulation.md`
  (superseded), `docs/adr/0010-bdf-time-integration.md`.
- **Goal:** Time-dependent drift-diffusion with BDF1/BDF2 time integration
  for diode turn-on and switching analysis.
- **Deliverables:**
  - `semi/runners/transient.py` — BDF1/BDF2 time-stepper with
    configurable dt, t_end, and step output.
  - `benchmarks/pn_1d_turnon/` — 1D diode turn-on from 0 to 0.6 V
    over 500 ns; reference vs analytical exponential approach to SS.
  - `semi/fem/mass.py` — lumped mass matrix for the capacitive term.
  - BDF rate tests in `tests/fem/test_bdf_rates.py`.
- **Verification:**
  - Steady-state-limit test: transient at t → ∞ matches bias_sweep
    within tolerance.
  - BDF1/BDF2 rate tests: temporal convergence 1st / 2nd order.
- **Dependencies:** M12.

## M13.1: Transient close-out — Slotboom primary unknowns (ADR 0012, 0013, 0014)

- **Status:** Done (2026-04-27). PRs #44, #45, #46, #47, #54; tag v0.14.1.
  See `docs/adr/0012-sg-primitives.md`, `docs/adr/0013-bc-ramp-continuation.md`,
  `docs/adr/0014-slotboom-transient.md`.
- **Goal:** Close the xfailed `test_transient_steady_state_limit` and BDF
  rate tests by switching the transient runner to Slotboom primary unknowns
  (psi, phi_n, phi_p) and fixing the MUMPS rank-deficiency failure.
- **Deliverables:**
  - `semi/fem/scharfetter_gummel.py`, `semi/fem/sg_assembly.py` — SG flux
    primitives (ADR 0012).
  - BC-ramp continuation in `run_transient` (ADR 0013): smooth initial
    condition ramp from equilibrium to the target voltage.
  - Slotboom transient formulation in `run_transient` (ADR 0014);
    supersedes ADR 0009.
  - `solver.jacobian_shift` config field + direct petsc4py MUMPS workspace
    bump (`mat_mumps_icntl_14=200`) for the rank-deficient minority-side
    row in deep p-bulk.
- **Verification:**
  - `test_transient_steady_state_limit`: bias_sweep vs transient at deep SS
    within 1e-4 relative error.
  - BDF1/BDF2 rate tests: active (no xfail), pass.
  - `pn_1d_turnon` benchmark within 5%.
- **Dependencies:** M13.

## M14: AC small-signal analysis (ADR 0011)

- **Status:** Done (2026-04-26). See PLAN.md completed work log and
  `docs/adr/0011-ac-small-signal.md`.
- **Goal:** Linearised frequency-domain response (J + jωM) δu = −dF/dV δV
  solved via real 2×2 block factorisation; enable C-V analysis without
  running a time-domain simulation.
- **Deliverables:**
  - `semi/runners/ac_sweep.py` — AC sweep runner.
  - `schemas/input.v1.json` updated to schema version 1.2.0
    (`solver.type = "ac_sweep"`, `solver.dc_bias`, `solver.ac`).
  - `benchmarks/rc_ac_sweep/` — 1D pn diode at V_DC = -1 V; depletion
    capacitance vs analytical within 0.4%.
  - `tests/mms/test_ac_consistency.py` — three AC formulation invariants.
- **Verification:**
  - RC benchmark depletion capacitance within 0.4% of analytical.
  - AC MMS consistency tests pass.
  - 249 tests (pre-audit), coverage 95.x%, ruff clean.
- **Dependencies:** M13.1.

## M14.1: Differential capacitance via AC admittance (PR #38)

- **Status:** Done (2026-04-26). PR #38.
- **Goal:** Add a `mos_cap_ac` runner that computes MOS differential
  capacitance from Im(Y) / (2πf) so the mos_2d C-V benchmark can be
  verified by the AC path as well as the numerical dQ/dV path.
- **Deliverables:**
  - `semi/runners/mos_cap_ac.py` — `run_mos_cap_ac` runner.
  - Integration of `mos_cap_ac` into the `mos_2d` benchmark verifier.
  - `tests/fem/test_mos_cap_ac.py` — byte-identity claim: mos_cv and
    mos_cap_ac agree on Q_gate to machine precision.
- **Verification:**
  - Phase 1 audit case 03 confirms byte-identity (rel_err = 0.000 at all
    42 gate voltages).
- **Dependencies:** M14.

## M14.2: Axisymmetric (cylindrical) 2D MOSCAP (schema 1.3.0)

- **Status:** Done (2026-04-30). PR #64 (`a4649be`). PR #65
  (`docs/post-merge-cleanup`) follows up with the runner dispatch and
  Colab self-contained notebook.
- **Goal:** Add an axisymmetric (cylindrical) coordinate-system path so
  rotationally symmetric devices (MOSCAPs, vertical resistors, junction
  cylinders) can be solved on a 2D meridian mesh with the radial weight
  carried by the weak forms.
- **Deliverables:**
  - **Schema 1.3.0:** new top-level `coordinate_system` field accepts
    `"cartesian"` (default, unchanged behaviour) or `"axisymmetric"`.
    Cross-field validation enforces dimension == 2, non-negative
    radial extent, and rejects Dirichlet contacts on the symmetry
    axis r = 0. `SCHEMA_SUPPORTED_MINOR` bumped 2 → 3.
  - `semi/physics/axisymmetric.py`: r-weighted equilibrium Poisson
    (`build_equilibrium_poisson_form_axisym_mr`) and Slotboom
    drift-diffusion forms on the meridian half-plane.
  - `semi/cv.py` (pure-Python, no dolfinx): MOSCAP analytical
    parameters (V_fb, V_t, |phi_B|, W_dmax, C_ox, C_min) plus LF
    (quasi-static, brentq) and HF (depletion-clamp) C-V helpers; FEM
    postprocessors `compute_lf_cv_fem` and
    `compute_hf_cv_depletion_clamp`. scipy promoted to base
    dependencies because `lf_cv_quasistatic` uses `scipy.optimize.brentq`.
  - `semi/runners/mos_cap_ac.py` axisymmetric branch (PR #65):
    replaces `W_lat = extents[0][1] - extents[0][0]` with
    `L_gate = ∫_gate r ds_meridian` (assembled via
    `fem.assemble_scalar(fem.form(r * ds_gate))` + MPI allreduce);
    r-weights `charge_form` and `sensitivity_form` by
    `r = ufl.SpatialCoordinate(msh)[0]`.
  - `benchmarks/moscap_axisym_2d/`: gmsh `.geo` template for the
    meridian mesh, `moscap_axisym.json` config (Hu Fig. 5-18
    parameters: p-type body N_A = 5e16 cm^-3, T_ox = 10 nm,
    phi_ms = -0.95 V, gate radius 50 μm), `reference_cv.csv`
    analytical baseline (V_fb = -0.950 V, V_t = +0.181 V,
    |phi_B| = 0.399 V, W_dmax = 144 nm, C_min/C_ox = 0.173).
  - `notebooks/05_moscap_axisym_cv.ipynb`: end-to-end Colab notebook
    reproducing Hu Fig. 5-18. Self-contained bootstrap (FEM-on-Colab
    installer, gmsh CLI + Python bindings, repository clone with
    branch pin, scipy install fallback). 22 cells: analytical LF/HF
    plot, mesh visualisation (region-shaded tripcolor, three panels),
    equilibrium psi(r,z) contourf, full 81-point V_g sweep, FEM-vs-
    analytical C-V plot, CSV export.
  - Tests: `tests/check_axisym_moscap_math.py` (15/15 analytical
    anchors green), `tests/test_coordinate_system.py` (schema gate),
    `tests/test_moscap_axisym_cv.py` (dolfinx FEM smoke),
    `tests/test_moscap_axisym_cv_fem.py` (FEM-vs-reference
    regression; skip-clean when `fem_cv.csv` is absent).
  - `docs/theory/axisymmetric.md` and `docs/theory/moscap_cv.md`:
    derivations and conventions.
- **Verification:**
  - Pure-Python: 237 passed / 22 skipped (M14.2 baseline); 15/15
    MOSCAP analytical anchors green.
  - FEM (Colab): notebook 05 runs end-to-end; FEM C(V) tracks
    analytical LF/HF curves through accumulation, depletion, and
    inversion regimes.
- **Dependencies:** M14.1.

## Deferred

These items remain in-tree as named follow-ups but are not on the
critical path. They are documented here, with the milestone that
supersedes them, so that a reviewer scanning for "what is left of
M14.2" finds the answer in one place.

- **M14.2.x: Cartesian-2D MOSCAP variant and rigorous gate-driven
  HF C-V method.** Originally tracked in
  [`PLAN.md`](../PLAN.md) §Backlog and on the README capability
  matrix. Superseded for practical purposes by M16.4 (Fermi-Dirac;
  the rigorous HF C-V formulation depends on FD-corrected statistics
  at the high doping that makes HF dispersion physically meaningful)
  and M19 (3D MOSFET; the same multi-region MOSCAP infrastructure
  exercised on a more important device than a 2D MOSCAP variant).
  If a contributor still wants to ship the Cartesian-2D MOSCAP
  variant in isolation, the work is bounded (~2 days) and can be
  picked up at any time; the rigorous HF C-V should not be
  attempted before M16.4.
- **Physics validation suite, Phase 2 (Sze, Nicollian-Brews).**
  Phase 1 is closed. Case 06 (transient FFT vs AC sweep) is the
  M16.7 deliverable; the rest (cases 01-05) are
  internal-consistency checks that pass and were never load-bearing
  on external references.

## Forward-looking

See [`PLAN.md`](../PLAN.md) for the active "Next task" and milestone
progression. See
[`docs/IMPROVEMENT_GUIDE.md`](IMPROVEMENT_GUIDE.md) for the full
milestone specs:

- M16.1 through M16.7: physics-completeness sub-milestones, each its
  own PR. M16.1 (Caughey-Thomas) is the next shipping PR; pick up via
  [`docs/M16_1_STARTER_PROMPT.md`](M16_1_STARTER_PROMPT.md).
- M17: heterojunctions (depends on M16.4).
- M19: 3D MOSFET capstone (depends on M16.1; the M14.3 XDMF / strict
  schema work is now landed).
- M19.1: MPI parallel orchestration (depends on M19).
- M20: HTTP server hardening (independent).

This file (`ROADMAP.md`) is the "what does this project do"
document and is updated when milestones merge. `PLAN.md` tracks
active work. `IMPROVEMENT_GUIDE.md` defines the acceptance criteria
for future milestones.
