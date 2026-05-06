# PHYSICS_INTRO: A Programmer's Walkthrough of What kronos-semi Actually Computes

**Audience.** You can program, you know calculus, you have maybe taken one
undergraduate physics or EE course a long time ago. You have never opened
a TCAD manual. You want to contribute to kronos-semi or write code that
consumes its output, and you need to understand what the solver is doing
without taking a graduate semester on semiconductor device physics first.

**Scope.** This document is the on-ramp. It tells you *what* the equations
mean and *why* they exist in the shape they do. For exact equations,
scaling conventions, and boundary-condition variants, go to
[PHYSICS.md](PHYSICS.md) — that's the reference. This is the tutorial.

**Time to read:** about 30 minutes if you take the math seriously, 10 if
you skim.

---

## 1. What is a semiconductor device simulator, really

Three things happen inside a piece of silicon when you apply a voltage:

1. **An electric field exists** because of charges and applied voltages.
2. **Charges move** because of that field (drift) and because charges
   diffuse from high concentration to low (diffusion).
3. **Electron–hole pairs appear and disappear** because of thermal
   generation and various recombination processes.

A device simulator solves a small set of coupled partial differential
equations (PDEs) describing these three phenomena simultaneously, on a
mesh that represents the geometry of your device. The output is a
self-consistent picture: at every point in the device, what is the
electrostatic potential ψ, what is the electron density n, what is the
hole density p, and what currents are flowing through the contacts.

kronos-semi is that kind of simulator, built on the FEniCSx finite-element
framework and driven by JSON input files.

---

## 2. The three equations, in English first

### 2.1 Poisson's equation: "charges make fields"

If you've seen electrostatics, this is the one you already know:

```
∇ · (ε ∇ψ) = -ρ
```

Translated: the divergence of the electric displacement field equals the
charge density. In a semiconductor, the charge density is built up from
four things:

- Positive ionized donors (N_D⁺, atoms that gave up an electron)
- Negative ionized acceptors (N_A⁻, atoms that grabbed an electron)
- Mobile electrons (n, negative)
- Mobile holes (p, positive — a hole is the absence of an electron in a
  crystal lattice, treated as a positive charge carrier)

So the right side becomes q(p − n + N_D⁺ − N_A⁻), where q is the
elementary charge. In equilibrium, "ionized" just means "all of them,"
which is what kronos-semi assumes.

**Code pointer:** `semi/physics/poisson.py::build_equilibrium_poisson_form`
for the equilibrium case, and the first term `F_psi` in
`semi/physics/drift_diffusion.py::build_dd_block_residual` for the
biased case.

### 2.2 The continuity equations: "charge is conserved"

Electrons and holes cannot appear or disappear without a paired partner
(you can't destroy charge). So the rate at which electrons flow into a
region must equal the rate at which electrons accumulate there, minus
the rate at which they recombine with holes:

```
∇ · J_n = +q·R     (electron continuity)
∇ · J_p = -q·R     (hole continuity)
```

Here J_n is the electron current density, J_p is the hole current
density, and R is the net recombination rate (positive means electrons
and holes are disappearing in pairs). The signs are opposite because
electrons are negative and holes are positive — a recombination event
removes one of each, but the divergence of their currents has opposite
signs.

In steady state there is no time derivative on either side, so a solve
is just "find n(x), p(x), ψ(x) such that charges balance and currents
are conservative everywhere."

**Code pointer:** `F_phi_n` and `F_phi_p` in
`semi/physics/drift_diffusion.py::build_dd_block_residual`.

### 2.3 Drift-diffusion currents: "how carriers move"

Each current is the sum of a drift term (carriers pushed by the field)
and a diffusion term (carriers spreading from high concentration to low):

```
J_n = q μ_n n E + q D_n ∇n     = -q μ_n n ∇ψ + q D_n ∇n
J_p = q μ_p p E - q D_p ∇p     = -q μ_p p ∇ψ - q D_p ∇p
```

where μ_n, μ_p are mobilities (how fast carriers drift per unit field),
D_n, D_p are diffusivities, and E = -∇ψ is the electric field. The
Einstein relation D = μ·V_t, where V_t = kT/q ≈ 25.85 mV at room
temperature, links them. Substituting everything into the continuity
equations gives you the drift-diffusion system.

That is literally the whole physics of a classical semiconductor
simulator. Everything else — band bending, depletion regions, MOS
inversion, diode IV curves, MOSFET saturation — *falls out* of those
three equations plus the right boundary conditions.

---

## 3. Why we don't solve those equations directly

If you implemented the equations above naively in FEniCSx, it would
diverge. Three reasons:

### 3.1 The numbers span 40 orders of magnitude

In SI units, a typical 1 µm device at 10¹⁷ cm⁻³ doping has:

- Permittivity ε ≈ 10⁻¹¹ F/m
- Elementary charge q ≈ 10⁻¹⁹ C
- Density N ≈ 10²³ m⁻³
- Length L ≈ 10⁻⁶ m

Plugging those into the Jacobian of Poisson's equation gives entries
that span from 10⁻¹¹ to 10²³. Newton's method on that produces a matrix
with a condition number around 10³⁰. Floating-point arithmetic laughs
at you and returns nonsense. The fix is nondimensionalization: scale
every variable so the scaled versions are all O(1). See PHYSICS.md §2
for the full recipe. The short version:

- Potentials are measured in units of V_t ≈ 25.85 mV.
- Densities are measured in units of C₀ (peak doping, typically ~10¹⁷ cm⁻³).
- Lengths stay in meters. This is on purpose — it lets JSON users
  specify geometry in SI, which is less error-prone than a custom
  length scale.

After scaling, the only small number left is λ² = ε V_t / (q C₀ L₀²),
which equals the squared ratio of the Debye length to the device
length. For a 1 µm device at 10¹⁷ cm⁻³, λ² ≈ 10⁻⁴. That's a singular
perturbation — small λ² means thin boundary layers at junctions and
interfaces, which is exactly the depletion-region physics. The solver
handles it via the scaling, not by ignoring it.

**Code pointer:** `semi/scaling.py`.

### 3.2 Carrier densities vary exponentially

Electron density n varies as exp(ψ/V_t), which at V_t = 25.85 mV means
a 1 V swing changes n by a factor of e^40 ≈ 10¹⁷. The Jacobian of
Newton's method gets entries that vary by the same ratio, which
destroys conditioning even worse than the raw-units problem.

The standard fix is the **Slotboom transformation**. Instead of solving
for n and p directly, we solve for the "quasi-Fermi potentials" Φ_n and
Φ_p, defined by:

```
n = n_i · exp((ψ - Φ_n) / V_t)
p = n_i · exp((Φ_p - ψ) / V_t)
```

where n_i is the intrinsic carrier density (about 10¹⁰ cm⁻³ in silicon
at room temperature). Φ_n is roughly "how far above the Fermi level are
the electrons," in volts. At thermal equilibrium Φ_n = Φ_p = 0, and n·p
equals n_i² identically — this is the mass-action law.

When you rewrite the continuity equations in terms of Φ_n and Φ_p, a
miracle happens: the currents become pure gradients of the Φ's:

```
J_n = -q μ_n n ∇Φ_n
J_p = -q μ_p p ∇Φ_p
```

which means the weak form is well-posed under Galerkin FEM without any
upwind stabilization. This is why every modern device simulator uses
some form of Slotboom or Scharfetter-Gummel discretization. See ADR
0004 (`docs/adr/0004-slotboom-variables-for-dd.md`) for why we picked
Slotboom over SG.

**Code pointer:** `semi/physics/slotboom.py` and
`semi/physics/drift_diffusion.py`.

### 3.3 Forward bias breaks Newton from a cold start

At forward bias of 0.6 V on a silicon diode, carrier injection raises n
and p by a factor of about exp(0.6/0.0258) ≈ 10¹⁰ over equilibrium. If
you try to Newton-solve for the biased state starting from equilibrium,
you overshoot by many orders of magnitude on the first step and never
recover.

The fix is **continuation**: don't jump to V = 0.6 V directly. Start at
V = 0, walk up in small steps (say 0.02 V), solve at each one using the
previous solution as the initial guess. If a step fails, halve it and
retry. If three steps in a row converge in under four Newton iterations,
grow the step size. This is standard arc-length-style homotopy.

**Code pointer:** `semi/continuation.py` for the controller,
`semi/runners/bias_sweep.py` for the walk.

---

## 4. Boundary conditions, a.k.a. "where contacts happen"

A device is a chunk of semiconductor with metal contacts attached at
specific faces. Each contact imposes a boundary condition on the PDEs,
and different physical types of contact give different BCs.

### 4.1 Ohmic contacts

An ohmic contact is one with no rectifying behavior — it passes current
both ways linearly. Mathematically, we treat the metal as a perfect
reservoir that holds the semiconductor at a fixed potential and
maintains charge neutrality (N_net + p - n = 0) right at the contact.
That gives us Dirichlet conditions on ψ, Φ_n, and Φ_p:

- ψ = V_applied + ψ_built-in, where ψ_built-in is the potential that
  satisfies charge neutrality at the doping profile present at the
  contact (n-type contact → positive built-in, p-type contact →
  negative built-in).
- Φ_n = Φ_p = V_applied.

**Code pointer:** `semi/bcs.py::build_psi_dirichlet_bcs` and
`build_dd_dirichlet_bcs`.

### 4.2 Gate contacts

A gate is a metal electrode separated from the semiconductor by an
insulator (for silicon MOS, the insulator is SiO₂). The gate holds its
own potential but there's no carrier exchange. In kronos-semi:

- The oxide is a meshed region with its own permittivity, no doping, no
  mobile carriers.
- ψ is continuous across the Si/SiO₂ interface with a jump in ε·∇ψ (the
  displacement field is continuous, not the field itself).
- Φ_n and Φ_p live only on the semiconductor — we use a submesh for
  them, since quasi-Fermi potentials are ill-defined in an ideal
  insulator.
- The gate contact applies ψ_gate = V_applied + φ_ms (work-function
  difference between gate metal and semiconductor) on the top surface
  of the oxide. No BC is applied to Φ_n, Φ_p at the gate because they
  don't exist there.

**Code pointer:** `semi/mesh.py::build_submesh_by_role`,
`semi/physics/drift_diffusion.py::build_dd_block_residual_mr`.

### 4.3 Insulating / reflective boundaries

Everywhere else on the device boundary — the "sides" of the mesh that
aren't contacts — we apply homogeneous Neumann (zero-gradient) on all
three fields. This means "no current, no field leaking out." It's
sometimes called a reflecting boundary and is the default when you
don't say otherwise.

Natural to the weak form: if you look at the UFL in
`build_dd_block_residual`, you'll see only volume integrals over `dx`
and no boundary integrals. Homogeneous Neumann is what you get when you
integrate by parts and then drop the boundary term. If you ever see
someone adding a `ds` integral, that's a non-trivial flux BC.

---

## 5. A toy example: the 1D pn junction

Let's walk through the simplest realistic case.

**Device.** A 2 µm-long piece of silicon, with the left micron doped
p-type (N_A = 10¹⁷ cm⁻³) and the right micron doped n-type (N_D = 10¹⁷
cm⁻³). Contacts at both ends.

**What should happen physically.** At the junction (x = 1 µm), there's
a sudden step change in doping. Electrons on the n-side diffuse left
into the p-side, where they recombine with holes. This leaves behind
positive ionized donors on the n-side and negative ionized acceptors
on the p-side. The resulting charge distribution sets up an electric
field that opposes further diffusion. Equilibrium is reached when the
drift current (field pushing carriers) exactly cancels the diffusion
current (concentration gradient pushing carriers the other way).

The region where the charge imbalance exists is called the **depletion
region**. The built-in potential across it is V_bi = V_t · ln(N_A N_D /
n_i²) ≈ 0.83 V for symmetric 10¹⁷/10¹⁷ doping. The depletion width is
about 150 nm, and the peak field is about 113 kV/cm. These are
textbook formulas, and they fall out of a zero-bias Poisson solve with
no continuity equations needed (at equilibrium, Φ_n = Φ_p = 0, so the
continuity equations are trivially satisfied).

**What kronos-semi does.**

1. Parse `benchmarks/pn_1d/pn_junction.json`.
2. Build a 1D interval mesh with 400 cells between 0 and 2 µm.
3. Tag the left facet as "anode," the right as "cathode."
4. Compute the doping function N_net(x) = N_D(x) - N_A(x) (negative on
   the left, positive on the right).
5. Set up the scaled Poisson equation with Dirichlet BCs at both contacts.
6. Start with the initial guess ψ(x) = V_t · asinh(N_net(x) / 2n_i),
   which is the local-charge-neutrality potential (this is the
   equilibrium ψ in the *bulk*, far from the junction; using it as an
   initial guess gives Newton a head start).
7. Solve with PETSc SNES using MUMPS LU. Converges in 4-6 Newton
   iterations.
8. Read off ψ, n, p at every mesh node. Compare V_bi, depletion width,
   and peak field against the analytical formulas. Done.

**Walk the code.** Start at `benchmarks/pn_1d/pn_junction.json`, then
`semi/schema.py::load`, then `semi/run.py::run`, which dispatches to
`semi/runners/equilibrium.py::run_equilibrium`. Everything else is
details. The whole solve, including mesh build, takes about a second
on a laptop.

---

## 6. What else can the solver do today

Beyond equilibrium, kronos-semi handles:

- **Bias sweeps.** Apply a non-zero voltage, solve the full coupled
  (psi, phi_n, phi_p) system. `benchmarks/pn_1d_bias/` produces a
  forward-bias I-V curve that matches the Shockley diode equation
  within 10 %.
- **Reverse bias.** Same solver, negative voltage. Recovers the
  reverse saturation current from SRH generation. See
  `benchmarks/pn_1d_bias_reverse/`.
- **MOS capacitors.** Multi-region (silicon + oxide), gate contact,
  C-V sweep. Both the Cartesian path (`benchmarks/mos_2d/`) and the
  axisymmetric path (`benchmarks/moscap_axisym_2d/`, M14.2). The
  submesh for phi_n, phi_p is the trickiest part of the codebase.
- **MOSFETs.** 2D MOSFET with Pao-Sah analytical anchor; the
  Lombardi surface-mobility benchmark widens the verifier window
  into the inversion regime. See `benchmarks/mosfet_2d/`.
- **3D resistors.** Uniform doped bar, unstructured gmsh mesh,
  bipolar sweep. `benchmarks/resistor_3d/`.
- **Transient (BDF1/BDF2).** Time-dependent solves with adaptive
  step control. `semi/runners/transient.py` plus the 1D pn turn-on
  benchmark `benchmarks/pn_1d_turnon/`.
- **AC small-signal.** Linearized admittance around a DC operating
  point. `semi/runners/ac_sweep.py` plus the AC sensitivity path
  in `semi/runners/mos_cap_ac.py`. RC test case in
  `benchmarks/rc_ac_sweep/`.
- **Field-dependent mobility.** Caughey-Thomas velocity-saturation
  (M16.1, `benchmarks/diode_velsat_1d/`) and Lombardi composite
  surface mobility (M16.2; the mosfet_2d Pao-Sah window widens
  accordingly).
- **Auger recombination.** Closed-form additive kernel
  R = (C_n n + C_p p)(np - n_i^2) (M16.3,
  `benchmarks/diode_auger_1d/`).
- **Fermi-Dirac statistics.** Generalized-Slotboom substitution
  under the basic Blakemore approximation; the production path
  under `physics.statistics: "fermi_dirac"` (M16.4,
  `benchmarks/diode_fermi_dirac_1d/`). Boltzmann remains the
  default.
- **Schottky contacts.** Thermionic-emission Robin BC on the
  electron continuity row plus a metal-Fermi-level psi Dirichlet
  (M16.5, `benchmarks/schottky_1d/`). New in v0.21.0.

Verification & validation: every domain-physics module has an MMS
variant in `semi/verification/mms_dd.py` (Variants A through G).
Boundary-physics milestones use analytical-benchmark plus byte-
identity gates instead; see ADR 0015.

---

## 7. What it cannot do (yet)

If your mental model of "semiconductor simulator" is COMSOL
Semiconductor Module, the remaining gaps are:

- **No tunneling.** No band-to-band (Kane / Zener / Esaki), no
  trap-assisted (Hurkx). Required for any non-toy reverse-bias
  diode and for floating-gate flash modeling. M16.6 closes this.
- **No impact ionization.** Avalanche multiplication is not
  modeled; reverse-breakdown benchmarks rely on tunneling instead
  and are bounded above the avalanche regime. Not scheduled in
  M16; tracked in IMPROVEMENT_GUIDE.md as a post-M19 deliverable.
- **No heterojunctions.** Position-dependent chi and Eg are not
  yet supported. M17 closes this; depends on M16.4 (Fermi-Dirac)
  because heterojunctions break the non-degenerate approximation
  at the barrier.
- **No 3D MOSFET / FinFET / planar transistor.** The 3D coverage
  today is the doped resistor and a pure-Poisson box. M19 closes
  this with a 3D MOSFET on a gmsh-sourced unstructured mesh.
- **No transient FFT vs AC sweep validation.** The transient and
  AC paths are independently verified, but the cross-check that
  an FFT of the transient response matches the AC sweep at the
  same operating point is deferred. M16.7 closes this.

These are all on the roadmap. See
[IMPROVEMENT_GUIDE.md](IMPROVEMENT_GUIDE.md) § 4 for the per-
milestone numerical acceptance thresholds.

---

## 8. What to read next

Depending on what you want to do:

- **Contribute a new physics model.** Read PHYSICS.md in full, then the
  existing derivation documents (`mms_dd_derivation.md`,
  `mos_derivation.md`, `resistor_derivation.md`) to see how new models
  are introduced with derivations checked into the repo. Then read
  `docs/adr/0004-slotboom-variables-for-dd.md` and
  `docs/adr/0006-verification-and-validation-strategy.md` because
  anything new needs an MMS verifier.
- **Write a UI or API consumer.** Read [WALKTHROUGH.md](WALKTHROUGH.md)
  for how a JSON turns into a result, then IMPROVEMENT_GUIDE.md §3 and
  §4 (M9, M10) for the artifact and HTTP contracts.
- **Port to GPU or optimize the linear solver.** Read
  [IMPROVEMENT_GUIDE.md](IMPROVEMENT_GUIDE.md) §5. You will not be
  happy if you skip it.
- **Add a new benchmark.** Read `benchmarks/pn_1d/` as a template, then
  `scripts/run_benchmark.py` to see how verifiers are registered, then
  ADR 0006 for the standard of proof.

## 9. Recommended external references

If you want to go deeper on the device physics, two books are enough
to carry you through anything in kronos-semi:

- **Sze, Ng — Physics of Semiconductor Devices (3rd ed., 2007).** The
  standard reference. Read chapter 2 for carrier statistics and
  chapter 3 for pn junctions. Material parameters in kronos-semi's
  `materials.py` cite this book.
- **Selberherr — Analysis and Simulation of Semiconductor Devices
  (1984).** Old but still the clearest derivation of the
  drift-diffusion system and the numerical methods for it, including
  Slotboom and Scharfetter-Gummel.

For the finite-element side:

- **Logg, Mardal, Wells — Automated Solution of Differential Equations
  by the Finite Element Method (2012).** Free online. The FEniCS book,
  pre-FEniCSx but most concepts carry over. Chapter 1 is enough for
  what this repo does.
- **FEniCSx tutorial:** https://jsdokken.com/dolfinx-tutorial/ — the
  most up-to-date practical guide for dolfinx 0.10+.
