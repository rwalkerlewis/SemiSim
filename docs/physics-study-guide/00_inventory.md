# 00 — Physics-concept inventory

This inventory maps every physical concept, equation, parameter, boundary
condition, time integrator, linearization, or numerical method that
appears anywhere in the kronos-semi repository (v0.16.0 / schema v2.0.0
/ M16.1 shipped) to:

1. **Code location** — file path with line range or function name
2. **Schema field(s)** — JSON keys that expose the concept
3. **Existing-doc coverage** — where it is already documented
4. **Gap** — what a beginner needs that the existing docs don't supply,
   and which chapter of this study guide fills the gap

The inventory is organized by chapter bucket. For the chapter content
itself, see `01_classical_em_and_poisson.md` through
`appendix_D_milestone_physics_map.md`.

---

## A. Classical electromagnetism and Poisson's equation (Ch. 1)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Maxwell's equations → electrostatic limit | `semi/physics/poisson.py:25-79` (`build_equilibrium_poisson_form`) | n/a | `docs/PHYSICS.md` §1.1 (final form only) | First-principles derivation of why we drop $\partial \mathbf{B}/\partial t$ and $\partial \mathbf{D}/\partial t$ |
| Permittivity and displacement field | `semi/materials.py:38-40` (`Material.epsilon`) | `regions[*].material` | `docs/PHYSICS.md` §1.1 | Microscopic origin of ε; why insulators carry only $\varepsilon_r$ |
| Poisson PDE in matter $-\nabla\!\cdot\!(\varepsilon\nabla\psi)=\rho$ | `semi/physics/poisson.py:73-79` | n/a | `docs/PHYSICS.md` §1.1, `docs/PHYSICS_INTRO.md` §2.1 | Why we solve for $\psi$, not $\mathbf{E}$ |
| Dirichlet BC (ohmic, gate) | `semi/bcs.py:145-196` | `contacts[*].type == "ohmic"`, `"gate"` | `docs/PHYSICS.md` §3.1, §3.2 | First-principles motivation for "infinite recombination velocity" interpretation of ohmic |
| Homogeneous Neumann BC (insulating) | `semi/bcs.py:34, 100` (skipped); natural in weak form | `contacts[*].type == "insulating"` | `docs/PHYSICS.md` §3.3 | Why this is the natural BC of the variational principle |
| Interface flux continuity (Si/SiO₂) | `semi/physics/poisson.py:82-118` (`build_equilibrium_poisson_form_mr`) | `regions` with multiple `material` | `docs/PHYSICS.md` §6.2 | Derivation that the Galerkin form encodes $[\![\,\varepsilon\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!]=0$ for free |
| Symmetry-axis natural BC ($r=0$) | `semi/physics/axisymmetric.py:24-27`, schema cross-validator | `coordinate_system: "axisymmetric"` | `docs/theory/axisymmetric.md` | Forward reference from Ch. 1 to Ch. 15 |

## B. Solid state and band theory (Ch. 2)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Bandgap $E_g$, electron affinity $\chi$ | `semi/materials.py:31-32` | `regions[*].material` (looked up in DB) | `docs/PHYSICS.md` §4.2 (table only) | Where these come from physically |
| Effective DOS $N_c$, $N_v$ | `semi/materials.py:31-32, 56-57` | n/a (material DB) | `docs/PHYSICS.md` §4.2, Sze ref | 3D parabolic-band integral derivation |
| Direct vs indirect gap | (mentioned in `docs/talk/02-physics.md`) | n/a | `docs/talk/02-physics.md` | Not yet exposed in code; relevant for radiative recombination M16.3 |
| Work-function offset $\phi_{ms}$ | `semi/bcs.py:186-188` | `contacts[*].workfunction` | `docs/PHYSICS.md` §3.2 | Derivation of $V_{fb}$ from band alignment |

## C. Carrier statistics (Ch. 3)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Boltzmann approximation $n=n_i\exp((\psi-\Phi_n)/V_t)$ | `semi/physics/slotboom.py:27-42` | n/a | `docs/PHYSICS.md` §1.2 | Reduction Fermi–Dirac → Boltzmann; validity boundary |
| Mass action $np = n_i^2$ | `tests/check_analytical_math.py:78-80` | n/a | `docs/PHYSICS_INTRO.md` §3.2 | Proof from Fermi-level constancy at equilibrium |
| Intrinsic concentration $n_i$ | `semi/materials.py:33, 58, 70, 82` | n/a | `docs/PHYSICS.md` §4.2 | Why $n_i = \sqrt{N_c N_v}\exp(-E_g/2kT)$ |
| Thermal voltage $V_t = kT/q$ | `semi/constants.py:24-26`, `tests/check_analytical_math.py:16` | `physics.temperature` | `docs/PHYSICS.md` §4.1 | Dimensional argument; numerical value at 300 K |
| Fermi–Dirac (planned, M16.4) | not yet implemented | `physics.statistics: "fermi_dirac"` (planned) | `docs/IMPROVEMENT_GUIDE.md` §M16.4 | Forward reference; degenerate-doping motivation |

## D. Doping and charge neutrality (Ch. 4)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Donors/acceptors, complete ionization | `semi/doping.py:70-107` | `doping[*].profile` | `docs/PHYSICS.md` §1.1 | Why we assume complete ionization at 300 K |
| Net doping $N_D - N_A$ | `semi/doping.py:23-56` | implicit (signed sum) | `docs/PHYSICS_INTRO.md` §2.1 | Sign convention rationale |
| Uniform / step / Gaussian profiles | `semi/doping.py:70-107` | `doping[*].profile.type` | `docs/schema/reference.md` §doping | Physical meaning (implant, diffusion, abrupt junction) |
| Charge neutrality in bulk | `tests/check_analytical_math.py:82-93` | n/a | `docs/PHYSICS.md` §3.1 | Derivation of $\psi_\mathrm{eq}=V_t\,\mathrm{asinh}(N_\mathrm{net}/2n_i)$ |

## E. Drift-diffusion transport (Ch. 5)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Boltzmann transport equation | conceptual; not in code | n/a | (mentioned) | Reduction BTE → moment equations → DD |
| Einstein relation $D = \mu V_t$ | `semi/diode_analytical.py:32-33` | n/a | `docs/PHYSICS.md` §1.3 | Derivation from detailed balance |
| Continuity equations $\nabla\!\cdot\!\mathbf{J}_n = qR$ | `semi/physics/drift_diffusion.py:158-169` | n/a | `docs/PHYSICS.md` §1.3 | Sign convention proof |
| Drift-diffusion currents | `semi/postprocess.py:97-99` | n/a | `docs/PHYSICS_INTRO.md` §2.3 | Péclet number / non-monotonicity discussion |
| Mobility (constant) | `semi/physics/mobility.py:87-93` | `physics.mobility.mu_n`, `mu_p` | `docs/PHYSICS.md` §4.2 | Microscopic origin (scattering rates) |
| Mobility (Caughey–Thomas, M16.1) | `semi/physics/mobility.py:57-84, 119-217` | `physics.mobility.model: "caughey_thomas"` | `docs/M16_1_STARTER_PROMPT.md` | Velocity-saturation derivation |

## F. Generation–recombination (Ch. 6)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| SRH kernel $R = (np-n_i^2)/(\tau_p(n+n_1)+\tau_n(p+p_1))$ | `semi/physics/recombination.py:35-78` | `physics.recombination.{tau_n,tau_p,E_t}` | `docs/PHYSICS.md` §1.4 | Trap-state-balance derivation |
| Trap energy $E_t$, $n_1$, $p_1$ | `semi/physics/recombination.py:58-59` | `physics.recombination.E_t` | `docs/PHYSICS.md` §1.4 | Why $n_1 p_1 = n_i^2$ |
| Equilibrium $R = 0$ | `tests/test_recombination.py` | n/a | `docs/PHYSICS.md` §1.4 | Algebraic proof |
| Auger (planned, M16.3) | not implemented | `physics.recombination.auger` | `docs/IMPROVEMENT_GUIDE.md` §M16.3 | Three-particle motivation |
| Radiative (planned) | not implemented | (planned) | `docs/IMPROVEMENT_GUIDE.md` (post-submission) | Forward reference |

## G. pn junction physics (Ch. 7)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Built-in voltage (ln form) | `semi/diode_analytical.py:50, 78, 112` | n/a | `docs/PHYSICS_INTRO.md` §5 | Derivation from $\Phi_n^L = \Phi_n^R$ |
| Built-in voltage (asinh form) | `semi/bcs.py:183`, `semi/physics/slotboom.py:73-83` | n/a | `docs/PHYSICS.md` §3.1 | Equivalence proof for $N\gg n_i$ |
| Depletion width $W$, $x_n$, $x_p$ | `semi/diode_analytical.py:40-53` | n/a | `tests/check_analytical_math.py:50-55` | Charge-balance + double integration of Poisson |
| Peak field $|E_\mathrm{max}|$ | `tests/check_analytical_math.py:57-61` | n/a | values in tests | Derivation; numerical check 113 kV/cm |
| Shockley diffusion current | `semi/diode_analytical.py:19-37` | n/a | `docs/PHYSICS.md` §3.1 | Long-diode minority-carrier derivation |
| SRH generation in reverse bias | `semi/diode_analytical.py:96-121` | n/a | `docs/PHYSICS.md` §5.3 | Why $J_\mathrm{gen}\propto W(V)-W(0)$ |

## H. Metal–semiconductor and ohmic contacts (Ch. 8)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Ohmic contact idealization | `semi/bcs.py:181-189, 251-266` | `contacts[*].type == "ohmic"` | `docs/PHYSICS.md` §3.1 | "Infinite recombination velocity" interpretation |
| Gate contact idealization | `semi/bcs.py:186-189, 236-249` | `contacts[*].type == "gate"` | `docs/PHYSICS.md` §3.2 | Forward-link Schottky |
| Schottky thermionic emission (planned, M16.5) | `semi/bcs.py:190-193` raises | `contacts[*].type == "schottky"` (planned) | `docs/IMPROVEMENT_GUIDE.md` §M16.5 | Forward reference; full derivation in chapter |

## I. MOS capacitor physics (Ch. 9)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Flat-band voltage $V_{fb}$ | `semi/cv.py:127-128` (`MoscapAnalytic.V_fb`) | – | `docs/theory/moscap_cv.md` | Derivation from $\phi_{ms}$ and oxide charge |
| Threshold voltage $V_t^*$ | `semi/cv.py:135-139` | – | `docs/theory/moscap_cv.md` | Strong-inversion condition |
| Oxide capacitance $C_{ox}$ | `semi/cv.py:126` | – | `docs/theory/moscap_cv.md` | Parallel-plate derivation |
| Maximum depletion $W_\mathrm{dmax}$ | `semi/cv.py:131` | – | `docs/theory/moscap_cv.md` | Pinned at $\psi_s = 2\phi_B$ |
| LF / HF C–V | `semi/cv.py:159-360` | – | `docs/theory/moscap_cv.md` | Pedagogical contrast |
| BC-convention shift $V_{fb} = \phi_{ms} - \phi_F$ | `docs/PHYSICS.md` §6.3 | – | `docs/PHYSICS.md` §6.3 | Why kronos differs from textbook |

## J. MOSFET physics (Ch. 10)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Inversion-layer formation | implicit (Poisson + Slotboom DD) | – | `docs/PHYSICS.md` §6.6 | Sketch from band bending |
| Pao–Sah linear-regime expression | `tests/test_mosfet_2d_verifier.py` | – | `docs/PHYSICS.md` §6.6 | Derivation |
| Gaussian source/drain implants | `semi/doping.py:92-107` | `doping[*].profile.type == "gaussian"` | `docs/schema/reference.md` | Physical meaning |
| 2D MOSFET benchmark | `benchmarks/mosfet_2d/` | – | `benchmarks/mosfet_2d/README.md` | Verifier window rationale |

## K. Quasi-Fermi potentials and Slotboom transformation (Ch. 11)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Quasi-Fermi $\Phi_n, \Phi_p$ | `semi/physics/slotboom.py:55-70` | n/a | `docs/PHYSICS.md` §1.2-1.3 | Physical meaning ("Fermi level shift under bias") |
| Slotboom transformation | `semi/physics/slotboom.py:27-52` | – (chosen for engine) | `docs/theory/slotboom.md`, ADR 0004 | Why $\mathbf{J}_n$ becomes a pure gradient |
| Galerkin without SUPG | `docs/theory/slotboom.md` | – | `docs/theory/slotboom.md` | Why the discrete form is well-posed |
| Equilibrium $\Phi_n=\Phi_p=0$ | `semi/runners/bias_sweep.py:72-73` | n/a | `docs/PHYSICS.md` §1.2 | Why this is exact |

## L. Nondimensionalization (Ch. 12)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| $10^{30}$ condition number | `semi/scaling.py:1-10` | – | `docs/theory/scaling.md`, `docs/PHYSICS.md` §2 | Numerical example |
| Scales $L_0, V_t, C_0, \mu_0, t_0, J_0$ | `semi/scaling.py:32-71` | implicit (derived from cfg) | `docs/PHYSICS.md` §2.1 | Choice motivation |
| $\lambda^2$ definition | `semi/scaling.py:62-70` | – | `docs/PHYSICS.md` §2.3 | Singular perturbation interpretation |
| Debye length | `semi/scaling.py:72-80` | – | `docs/PHYSICS.md` §2.3 | Physical meaning |
| $L_D^2 = \lambda^2 L_0^2$ trap | `semi/physics/poisson.py:60-66` | – | `docs/PHYSICS.md` §2.3 (cautionary note) | Spell out the bug + lesson |
| Bias continuation | `semi/continuation.py`, `semi/runners/bias_sweep.py` | `solver.continuation` | `docs/PHYSICS.md` §2.4 | Homotopy framing |

## M. FEM weak forms (Ch. 13)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Strong→weak via test function | every `physics/*.py` form builder | – | none in this repo | First-principles intro for FEM novices |
| $H^1$ / $H^1_0$ Lagrange spaces | `semi/physics/drift_diffusion.py:55-63` | – | dolfinx tutorial | Conforming-space sketch |
| `dx(tag)` / `ds(tag)` measures | `semi/physics/poisson.py:107-110`, `semi/postprocess.py:84-86` | `mesh.regions_by_box`, `mesh.facets_by_plane` | `docs/schema/reference.md` | Subdomain integration explained |
| P1 vertex DOFs | `semi/physics/drift_diffusion.py:57-59` | – | – | Lagrange interpolation primer |

## N. Multi-region and interfaces (Ch. 14)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Cellwise $\varepsilon_r$ via DG0 | `semi/mesh.py::build_eps_r_function` | `regions` | `docs/PHYSICS.md` §6.2 | What "DG0" means physically |
| Submesh for carriers | `semi/physics/drift_diffusion.py:173-216` | `regions` with role | `docs/PHYSICS_INTRO.md` §4.2 | Why $\Phi_{n,p}$ are ill-defined in oxide |
| `entity_maps` | `semi/physics/drift_diffusion.py:219-327` | – | code comment | Conceptual explanation |
| Flux continuity | `semi/physics/poisson.py:82-118` | – | `docs/PHYSICS.md` §6.2 | Derivation that Galerkin gives it for free |

## O. Axisymmetric formulations (Ch. 15)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Cylindrical coordinates | `semi/physics/axisymmetric.py:1-46` | `coordinate_system: "axisymmetric"` | `docs/theory/axisymmetric.md` | Geometric picture |
| Meridian half-plane | `semi/physics/axisymmetric.py:5-13` | – | `docs/theory/axisymmetric.md` | What it means |
| r-weighted Poisson | `semi/physics/axisymmetric.py:55-103` | – | `docs/theory/axisymmetric.md` | Derivation from $dV_{3D} = 2\pi r\,dr\,dz$ |
| r-weighted DD | `semi/physics/axisymmetric.py:142-219` | – | `docs/theory/axisymmetric.md` | Same |
| Forbidden Dirichlet at $r=0$ | `semi/schema.py::_validate_coordinate_system` | – | `docs/theory/axisymmetric.md` | Geometric reason |
| MOSCAP benchmark (Hu Fig. 5-18) | `benchmarks/moscap_axisym_2d/` | `coordinate_system: "axisymmetric"` | `docs/benchmarks/moscap_axisym_2d.md` | Reproduction strategy |

## P. Nonlinear solvers and continuation (Ch. 16)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Newton on residual | `semi/solver.py:186-249` | – | `docs/theory/dolfinx_choice.md` | Textbook Newton primer |
| Backtracking line search | `semi/solver.py:21` (`snes_linesearch_type: bt`) | `solver.snes` | `docs/adr/0008-snes-tolerances.md` | Why pure Newton diverges |
| `NonlinearProblem` (PETSc SNES wrap) | `semi/solver.py:217-249` | – | `docs/theory/dolfinx_choice.md`, ADR 0003 | Why dolfinx 0.10 over earlier |
| `AdaptiveStepController` | `semi/continuation.py:32-133` | `solver.continuation.{max_step,min_step,max_halvings,easy_iter_threshold,grow_factor}` | `docs/PHYSICS.md` §2.4 | Homotopy interpretation |
| Bipolar sweep | `semi/runners/bias_sweep.py:347-361` | `voltage_sweep.{start,stop,step}` spanning 0 | `docs/PHYSICS.md` §7.3 | Why two legs |
| Jacobian shift | `semi/solver.py:121-149` | `solver.jacobian_shift` | code comment | Why minority-side rows are rank-deficient |

## Q. Transient time integration (Ch. 17)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| BDF1 (backward Euler) | `semi/timestepping.py:54-65` | `solver.order: 1` | ADR 0010 | Stability argument |
| BDF2 | `semi/timestepping.py:54-65` | `solver.order: 2` | ADR 0010 | Truncation-error sketch |
| Slotboom transient | `semi/runners/transient.py:528-668` | `solver.type: "transient"` | ADR 0014 | Why ADR 0014 supersedes ADR 0009 |
| Implicit-only choice | `semi/runners/transient.py:1-62` | – | ADR 0010 | Stiffness argument |
| BC-ramp continuation IC | `semi/runners/transient.py:671-789` | `solver.bc_ramp_steps` | ADR 0013 | Why two-stage IC works |
| `pn_1d_turnon` benchmark | `benchmarks/pn_1d_turnon/` | – | benchmark README | Step-bias-at-$t=0$ scenario |
| Time-varying $V(t)$ (planned, M16.7) | not implemented | `contacts[*].voltage_t` (planned) | `docs/IMPROVEMENT_GUIDE.md` §M16.7 | Forward reference |

## R. Small-signal AC analysis (Ch. 18)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Linearization $(J + j\omega M)\delta u = -dF/dV\,\delta V$ | `semi/runners/ac_sweep.py:151-499` | `solver.type: "ac_sweep"` | ADR 0011 | Phasor derivation |
| Real 2×2 block | `semi/runners/ac_sweep.py:518-591` | – | ADR 0011 | Why we don't use complex PETSc |
| Mass matrix in Slotboom | `semi/runners/ac_sweep.py:326-387` | – | ADR 0011 Errata #2 | Chain-rule mass |
| Sign convention (current INTO device) | ADR 0011 Errata | – | ADR 0011 | Why the audit flipped sign |
| `mos_cap_ac` analytic dQ/dV | `semi/runners/mos_cap_ac.py:59-329` | `solver.type: "mos_cap_ac"` | `docs/theory/moscap_cv.md` | Sensitivity derivation |
| `rc_ac_sweep` benchmark | `benchmarks/rc_ac_sweep/` | – | benchmark README | 0.4% acceptance |

## S. Linear solvers and GPU backends (Ch. 19)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| MUMPS direct LU | `semi/solver.py:20-30` | `solver.backend: "cpu-mumps"` | ADR 0008 | What sparse direct factorization does |
| AMGX preconditioner | `semi/compute.py` (M15) | `solver.backend: "gpu-amgx"` | `docs/gpu.md` | Algebraic multigrid sketch |
| hypre BoomerAMG | `semi/compute.py` | `solver.backend: "gpu-hypre"` | `docs/gpu.md` | Parallel AMG sketch |
| `auto` fallback | `semi/compute.py` | `solver.backend: "auto"` | `docs/gpu.md` | Capability negotiation |
| `GET /capabilities` | `kronos_server/routes/health.py` | – | `docs/gpu.md` | UI integration |
| Bit-identity guarantee | `docs/gpu.md` | – | `docs/gpu.md` | Why CPU is reference |

## T. Material parameter database (Ch. 20)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| Si parameters | `semi/materials.py:50-61` | `regions[*].material: "Si"` | `docs/PHYSICS.md` §4.2 | Provenance (Sze Table 7, Altermatt $n_i$) |
| Ge, GaAs | `semi/materials.py:62-85` | `material: "Ge"`/`"GaAs"` | – | Same |
| SiO₂, HfO₂, Si₃N₄ | `semi/materials.py:86-88` | `material: "SiO2"`/`"HfO2"`/`"Si3N4"` | `docs/PHYSICS.md` §4.3 | Why insulators have only $\varepsilon_r$ |
| `MATERIALS` dict | `semi/materials.py:49-89` | `GET /materials` endpoint | `docs/ARCHITECTURE.md` | UI integration |

## U. Verification and benchmarks (Ch. 21)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| MMS philosophy | `semi/verification/mms_poisson.py`, `mms_dd.py` | – | ADR 0006 | What "manufactured" means |
| Conservation checks | `semi/verification/conservation.py` | – | `docs/PHYSICS.md` §5.3 | Why these are independent witnesses |
| Mesh convergence | `semi/verification/mesh_convergence.py` | – | `docs/PHYSICS.md` §5.2 | Cauchy-ratio interpretation |
| `tests/check_analytical_math.py` | `tests/check_analytical_math.py` | – | – | Reproduce every assertion |
| pn_1d benchmark | `benchmarks/pn_1d/` | – | benchmark README | Equilibrium target |
| pn_1d_bias / reverse | `benchmarks/pn_1d_bias{,_reverse}/` | – | benchmark READMEs | Shockley + SNS verification |
| mos_2d, mos_cv, mos_cap_ac | `benchmarks/mos_2d/` | – | benchmark README | C–V validation |
| mosfet_2d | `benchmarks/mosfet_2d/` | – | benchmark README | Pao–Sah window 20% |
| resistor_3d | `benchmarks/resistor_3d/` | – | benchmark README | 1% V–I linearity |
| moscap_axisym_2d | `benchmarks/moscap_axisym_2d/` | `coordinate_system: "axisymmetric"` | benchmark README | Hu Fig. 5-18 |
| rc_ac_sweep | `benchmarks/rc_ac_sweep/` | `solver.type: "ac_sweep"` | benchmark README | 0.4% AC depletion |
| pn_1d_turnon | `benchmarks/pn_1d_turnon/` | `solver.type: "transient"` | benchmark README | 5% transient target |

## V. Constants, units, and unit conversions (Appendix A)

| Concept | Code | Schema | Existing docs | Gap |
|---|---|---|---|---|
| 2019-redefinition $q$, $k_B$, $\varepsilon_0$, $\hbar$, $m_0$ | `semi/constants.py:14-21` | – | `docs/PHYSICS.md` §4.1 | Provenance |
| cm⁻³ ↔ m⁻³ | `semi/constants.py:29-46` | – | `docs/PHYSICS.md` §1 (preamble) | Mixed-unit pitfalls |

## W. Milestone physics map (Appendix D)

Pulled from `docs/ROADMAP.md` and `docs/IMPROVEMENT_GUIDE.md`:

| Milestone | New physics | Status | Chapter |
|---|---|---|---|
| M1 | Equilibrium Poisson, Boltzmann | Done | Ch. 1, 3, 7 |
| M2 | Coupled Slotboom DD, SRH, ohmic BCs | Done | Ch. 5, 6, 8, 11 |
| M3 | Adaptive bias-ramp continuation | Done | Ch. 16 |
| M4 | V&V suite (MMS, conservation) | Done | Ch. 21 |
| M5 | Refactor (no new physics) | Done | – |
| M6 | 2D MOS capacitor, multi-region | Done | Ch. 9, 14 |
| M7 | 3D resistor, gmsh ingest | Done | Ch. 13, 14, 21 |
| M8 | Submission polish | Done | – |
| M9 | Result artifact writer | Done | – |
| M10 | HTTP server | Done | – |
| M11 | Schema versioning | Done | – |
| M12 | Gaussian implants, MOSFET 2D | Done | Ch. 4, 10 |
| M13 / M13.1 | Transient BDF1/BDF2 (Slotboom) | Done | Ch. 17 |
| M14 | AC small-signal | Done | Ch. 18 |
| M14.1 | mos_cap_ac analytic dQ/dV | Done | Ch. 18 |
| M14.2 | Axisymmetric 2D MOSCAP | Done | Ch. 15 |
| M14.3 | Strict schema, XDMF, Pao-Sah verifier | Done | Ch. 10, 13, 21 |
| M15 | GPU linear-solver path | Done | Ch. 19 |
| M16.1 | Caughey–Thomas mobility | Done | Ch. 5 |
| M16.2 | Lombardi surface mobility | Planned | Ch. 5 (forward) |
| M16.3 | Auger recombination | Planned | Ch. 6 (forward) |
| M16.4 | Fermi–Dirac statistics | Planned | Ch. 3 (forward) |
| M16.5 | Schottky contacts | Planned | Ch. 8 (forward) |
| M16.6 | BBT and TAT tunneling | Planned | Ch. 6 (forward) |
| M16.7 | Time-varying $V(t)$ | Planned | Ch. 17 (forward) |
| M17 | Heterojunctions | Planned | Ch. 2 (forward) |
| M19 | 3D MOSFET capstone | Planned | Ch. 10 (forward) |
| M19.1 | MPI orchestration | Planned | – |
| M20 | HTTP server hardening | Planned | – |

---

## GAPS

The following items are referenced in code or docs but not yet given a
chapter section. They are flagged here so future revisions of this guide
can fill them:

- **Lombardi surface mobility (M16.2).** Full surface-roughness /
  Coulomb / phonon decomposition. Currently only a forward reference
  in Ch. 5; the full derivation will land when M16.2 ships.
- **Hurkx trap-assisted tunneling (M16.6).** Forward reference only.
- **Kane band-to-band tunneling (M16.6).** Forward reference only.
- **Heterojunction quantum confinement.** Out of scope for the DD
  engine; M17 captures only the classical band-bending piece.
- **MPI collective-communication audit (M19.1).** Numerical concern,
  not a physics concern; intentionally not chapter-mapped.

These are documented as future work in `docs/IMPROVEMENT_GUIDE.md`.
