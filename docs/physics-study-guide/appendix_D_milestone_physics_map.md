# Appendix D — Milestone physics map

This appendix maps every shipped milestone (M1–M16.1) and every open
milestone (M16.2 onward) to the physics it required or will require,
with cross-references into the chapters of this study guide.

Source: [`docs/ROADMAP.md`](../ROADMAP.md), [`docs/IMPROVEMENT_GUIDE.md`](../IMPROVEMENT_GUIDE.md).

## Shipped milestones

### M1 — Equilibrium Poisson, 1D pn junction (2026-04-20)
**Physics introduced:**
- Poisson's equation in matter (Ch. 1).
- Boltzmann carrier statistics (Ch. 3).
- Charge neutrality, equilibrium $\psi$ (Ch. 4).
- Built-in voltage (Ch. 7).
- Depletion approximation, peak field (Ch. 7).
- Nondimensionalization, $\lambda^2$ (Ch. 12).
- Ohmic Dirichlet BC (Ch. 1, Ch. 8).

**Verifier:** `pn_1d` benchmark vs depletion-approximation analytics (Ch. 21).

### M2 — Coupled drift-diffusion, Slotboom, SRH (2026-04-20)
**Physics introduced:**
- Drift-diffusion currents and continuity (Ch. 5).
- Einstein relation (Ch. 5).
- Slotboom transformation (Ch. 11).
- SRH recombination kernel (Ch. 6).
- Forward-bias diode (Ch. 7).
- Bias continuation (Ch. 16; refined in M3).

**Verifier:** `pn_1d_bias` vs Shockley diode equation.

### M3 — Adaptive bias-step continuation (2026-04-21)
**Physics:** No new physics; numerical refinement of M2.
**Numerics introduced:**
- `AdaptiveStepController` (Ch. 16).
- Sah-Noyce-Shockley recombination current correction (Ch. 7, Ch. 21).
- Reverse-bias generation current (Ch. 6, Ch. 7).

**Verifiers:** Forward (`pn_1d_bias`) tightened with SNS; new `pn_1d_bias_reverse`.

### M4 — V&V suite (2026-04-21)
**No new physics**; verification infrastructure.
**Introduced:**
- MMS for Poisson (Ch. 21).
- MMS for DD with three variants (Ch. 21).
- Mesh-convergence study (Ch. 21).
- Conservation gates (Ch. 21).
- ADR 0006.

### M5 — Refactor and test pass (2026-04-21)
**No new physics**; code reorganization and `bcs.py` extraction
(Ch. 8 code anchors).

### M6 — 2D MOS capacitor, multi-region (2026-04-21)
**Physics introduced:**
- Multi-region dielectric (Si/SiO₂) coupling, flux continuity (Ch. 14).
- MOS capacitor band bending: accumulation, depletion, inversion (Ch. 9).
- Flat-band, threshold, $W_\mathrm{dmax}$, $C_{ox}$, $C_\mathrm{min}$ (Ch. 9).
- Gate contact with $\phi_{ms}$ (Ch. 9, Ch. 8).
- Submesh for carriers in multi-region (Ch. 14).

**Verifier:** `mos_2d` C-V vs depletion-approximation theory; multi-region MMS (Ch. 21).

### M7 — 3D doped resistor (2026-04-21)
**Physics:** No new physics (dimension extension).
**Introduced:**
- gmsh `.msh` ingest with physical groups (Ch. 13).
- Bipolar bias sweep (Ch. 16).
- 3D V-I linearity verifier.

**Verifier:** `resistor_3d` vs $R = L/(qN_D\mu_n A)$ within 1%.

### M8 — Submission polish (2026-04-23)
**No new physics**; documentation and notebook deliverables.

### M9 — Result artifact writer (2026-04-23)
**No new physics**; persistence layer.
**Introduced:**
- `manifest.json` schema (Ch. 19 for KSP/wall-time fields).
- `semi-run` CLI.

### M10 — HTTP server (2026-04-23)
**No new physics**; service layer.
**Introduced:**
- FastAPI app, WebSocket progress (Ch. 19 for `GET /capabilities`).

### M11 — Schema versioning (2026-04-23)
**No new physics**; schema infrastructure.

### M12 — Mesh input + 2D MOSFET, ADR 0008 (2026-04-25)
**Physics introduced:**
- Gaussian source/drain implants (Ch. 4, Ch. 10).
- 2D MOSFET inversion-layer formation (Ch. 10).
- SNES tolerance amendments (Ch. 16).

**Verifier:** Initial 2D MOSFET I-V (refined in M14.3 with Pao-Sah).

### M13 — Transient solver, BDF1/BDF2 (2026-04-25, superseded by M13.1)
**Physics introduced:**
- BDF1, BDF2 time integration (Ch. 17).
- Lumped mass matrix (Ch. 17).
- ADR 0009 (transient form, superseded), ADR 0010 (BDF choice).

### M13.1 — Slotboom transient close-out (2026-04-27)
**Physics introduced:**
- Slotboom transient with chain-rule mass (Ch. 17).
- BC-ramp continuation IC (Ch. 16, Ch. 17).
- Jacobian shift for minority-side stability (Ch. 17).
- ADR 0014 (Slotboom transient, supersedes ADR 0009).
- ADR 0013 (BC-ramp).

**Verifier:** `pn_1d_turnon` benchmark within 5% of $\tau_p$.

### M14 — AC small-signal analysis (2026-04-26)
**Physics introduced:**
- Small-signal AC: $(J + j\omega M)\delta u$ system (Ch. 18).
- Real 2×2 block reformulation (Ch. 18).
- Displacement vs conduction current (Ch. 18).
- Sign convention errata (Ch. 18, ADR 0011).
- ADR 0011.

**Verifier:** `rc_ac_sweep` within 0.4% of analytical depletion C.

### M14.1 — MOSCAP differential capacitance via AC (2026-04-26)
**Physics introduced:**
- Poisson sensitivity / analytic $dQ/dV$ (Ch. 18).
- `mos_cap_ac` runner.

**Verifier:** Byte-identity with `mos_cv` on $Q_\mathrm{gate}$ at all 42 $V_g$ points.

### M14.2 — Axisymmetric 2D MOSCAP (2026-04-30)
**Physics introduced:**
- r-weighted Poisson (Ch. 15).
- r-weighted Slotboom DD (Ch. 15).
- Axis-symmetry natural BC (Ch. 1, Ch. 15).
- Hu Fig. 5-18 reproduction (Ch. 9, Ch. 15).
- LF / HF C-V helpers (Ch. 9).

**Verifier:** `moscap_axisym_2d` vs Hu Fig. 5-18 analytical reference.

### M14.3 — Housekeeping (v0.16.0, 2026-05-22)
**Physics introduced:**
- Pao-Sah linear-regime MOSFET verifier (Ch. 10).
- XDMF mesh ingest (Ch. 13).

**Verifier:** `mosfet_2d` Pao-Sah within 20% in $[V_T+0.2, V_T+0.6]$ V.

### M15 — GPU linear-solver path (v0.15.0)
**Physics:** No new physics.
**Numerics introduced:**
- `solver.backend` field (Ch. 19).
- AMGX, hypre BoomerAMG (Ch. 19).
- PETSc-CUDA / PETSc-HIP (Ch. 19).
- `GET /capabilities` endpoint (Ch. 19).

**Verifier:** `poisson_3d_gpu` with $5\times$ linear-solve speedup acceptance.

### M16.1 — Caughey-Thomas field-dependent mobility (v0.17.0, 2026-05-01)
**Physics introduced:**
- Caughey-Thomas closed-form velocity saturation (Ch. 5).
- Velocity-saturation in MOSFET (Ch. 10).
- MMS-DD Variant D acceptance gate (Ch. 21).

**Verifier:** `diode_velsat_1d` divergence/convergence anchors; MMS-DD Variant D
$L^2 \ge 1.99$ on every block.

## Open milestones

### M16.2 — Lombardi surface mobility (planned)
**Physics:** Surface scattering for inversion-layer mobility in MOSFETs (Ch. 5).

**Acceptance:** `mosfet_2d` Pao-Sah verifier tightened from 20% to 10% in
the inversion-strong-field window.

**Depends on:** M16.1.

### M16.3 — Auger recombination (planned)
**Physics:** Three-particle Auger recombination $R = (C_n n + C_p p)(np - n_i^2)$ (Ch. 6).

**Acceptance:** `diode_auger_1d` vs analytical high-injection long-diode reference
within 5%; existing benchmarks bit-identical with auger=false.

**Depends on:** M14.3.

### M16.4 — Fermi-Dirac statistics (planned)
**Physics:**
- Full Fermi-Dirac integral $F_{1/2}$ for degenerate doping (Ch. 3).
- Blakemore approximation for production (Ch. 3).
- Replaces `exp(±ψ)` calls in Slotboom builders.

**Acceptance:** `diode_fermi_dirac_1d` (1D pn, $N_D = 10^{20}\,\mathrm{cm^{-3}}$):
Boltzmann vs FD diverge by >15% on $V_{bi}$; FD matches scipy reference within $10^{-3}$.

**Depends on:** M14.3, M16.1. Hard prerequisite for M16.6 and M17.

### M16.5 — Schottky contacts (planned)
**Physics:**
- Thermionic-emission flux at metal-semiconductor contact (Ch. 8).
- Robin-style boundary form on continuity rows.
- Schema `type: "schottky"` with `barrier_height_eV` parameter.

**Acceptance:** `schottky_1d` vs thermionic-emission analytical I-V within 10%
on $V_F \in [0.1, 0.5]\,\mathrm{V}$.

**Depends on:** M14.3.

### M16.6 — BBT and TAT tunneling (planned)
**Physics:**
- Kane band-to-band tunneling rate (Ch. 6).
- Hurkx trap-assisted tunneling rate (Ch. 6).
- Generation kernel additions to `recombination.py`.

**Acceptance:** `zener_1d` reverse-bias breakdown vs Kane analytical reference
within 20% on $V_R \in [4, 8]\,\mathrm{V}$.

**Depends on:** M14.3, M16.4 (FD-corrected DOS).

### M16.7 — Time-varying transient contact voltage (planned)
**Physics:** No new physics; runner extension only.
**Numerics introduced:**
- Per-timestep contact voltage callable / table (Ch. 17).
- Closes audit case 06 (transient FFT ↔ AC sweep agreement).

**Acceptance:** Audit case 06 transient FFT vs `run_ac_sweep`'s $Y(\omega)$
within 5%.

**Depends on:** M14.3.

### M17 — Heterojunctions (planned)
**Physics:**
- Position-dependent electron affinity $\chi(\mathbf{x})$ (Ch. 2).
- Position-dependent bandgap $E_g(\mathbf{x})$ (Ch. 2).
- Cellwise DG0 fields for $\chi$, $E_g$ (Ch. 14).
- Heterojunction band-edge discontinuities.

**Acceptance:** `hemt_2d` AlGaAs/GaAs benchmark, 2DEG sheet density within
15% of self-consistent Poisson-Schrödinger reference.

**Depends on:** M16.4 (Fermi-Dirac for the heterojunction barrier).

### M19 — 3D MOSFET capstone (planned, slotted between M16 and M17)
**Physics:**
- 3D MOSFET inversion-layer formation (Ch. 10 extended).
- Saturation regime at $V_{DS} = 1\,\mathrm{V}$.

**Acceptance:** `mosfet_3d` Pao-Sah linear within 25%, velocity-saturation
within 30%; CPU/GPU speedup $\ge 5\times$ at 500 kDOFs.

**Depends on:** M16.1.

### M19.1 — MPI parallel benchmark (planned)
**No new physics**; orchestration audit.
**Acceptance:** `mosfet_3d` under `mpiexec -n {1,2,4}` same I_D within $10^{-8}$;
$n=4$ vs $n=1$ speedup $\ge 2.5\times$.

**Depends on:** M19.

### M20 — HTTP server hardening (planned)
**No physics**; auth/rate-limit/API-key infrastructure (Ch. 19 §HTTP).

## Summary table

| Milestone | New physics? | Primary chapter(s) | Status |
|---|---|---|---|
| M1 | Poisson + Boltzmann | 1, 3, 4, 7, 12 | Done |
| M2 | DD + Slotboom + SRH | 5, 6, 11 | Done |
| M3 | (numerical) | 16 | Done |
| M4 | (V&V) | 21 | Done |
| M5 | – | – | Done |
| M6 | Multi-region MOS | 9, 14 | Done |
| M7 | (3D extension) | 13, 14, 21 | Done |
| M8 | – | – | Done |
| M9 | – | – | Done |
| M10 | – | – | Done |
| M11 | – | – | Done |
| M12 | Gaussian implants, MOSFET | 4, 10 | Done |
| M13/M13.1 | Transient (Slotboom) | 17 | Done |
| M14 | AC small-signal | 18 | Done |
| M14.1 | MOSCAP $dQ/dV$ via Poisson sensitivity | 18 | Done |
| M14.2 | Axisymmetric 2D | 15 | Done |
| M14.3 | Pao-Sah verifier | 10, 21 | Done |
| M15 | GPU backends (numerics) | 19 | Done |
| M16.1 | Caughey-Thomas $\mu(F)$ | 5 | Done |
| M16.2 | Lombardi surface $\mu$ | 5 (forward) | Planned |
| M16.3 | Auger | 6 (forward) | Planned |
| M16.4 | Fermi-Dirac | 3 (forward) | Planned |
| M16.5 | Schottky contacts | 8 (forward) | Planned |
| M16.6 | BBT, TAT | 6 (forward) | Planned |
| M16.7 | Time-varying $V(t)$ | 17 (forward) | Planned |
| M17 | Heterojunctions | 2 (forward) | Planned |
| M19 | 3D MOSFET capstone | 10 (forward) | Planned |
| M19.1 | MPI orchestration | – | Planned |
| M20 | HTTP hardening | – | Planned |

## See also

- [`docs/ROADMAP.md`](../ROADMAP.md) — full delivery history with per-milestone deliverable / verification / dependencies.
- [`docs/IMPROVEMENT_GUIDE.md`](../IMPROVEMENT_GUIDE.md) — scope and acceptance tests for open milestones.
- [`docs/PLAN.md`](../../PLAN.md) — current task and immediate next steps.
- Chapter 21 §"Tolerance philosophy" — why each verifier accepts the
  tolerance it does, and how M16.x will tighten.
