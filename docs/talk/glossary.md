# Glossary of Terms and Acronyms

This glossary covers every technical term and acronym used across the
SemiSim talk slides (§1–§8). Entries are grouped by category for
easy review at the podium; an alphabetical quick-reference index
follows at the end.

---

## 1. Physics Terms

**Acceptor (N_A)**
A p-type dopant atom (e.g., boron in silicon) that introduces a hole
into the valence band when ionized. Acceptor concentration N_A sets
the hole-dominated carrier density in p-type material.

**Band gap (E_g)**
The energy separation between the valence band maximum and the
conduction band minimum. For silicon at 300 K, E_g ≈ 1.12 eV. Charge
carriers must receive at least this energy to be excited across the
gap.

**Boltzmann statistics**
An approximation to Fermi-Dirac statistics valid when the Fermi level
is at least 3 k_B T below the conduction band and above the valence
band (non-degenerate semiconductor). Carrier densities become
n = N_C exp((E_F − E_C)/k_B T). Valid for doping below ~10¹⁹ cm⁻³.

**Built-in potential (V_bi)**
The electrostatic potential difference that develops across a pn
junction at thermal equilibrium due to the charge redistribution.
V_bi ≈ V_t ln(N_A N_D / n_i²). For symmetric 10¹⁷ cm⁻³ silicon
doping, V_bi ≈ 0.83 V.

**Carrier density (n, p)**
The number of electrons (n) or holes (p) per unit volume in the
semiconductor, measured in cm⁻³. At equilibrium, n·p = n_i².

**Continuity equation**
A PDE expressing conservation of carriers: ∇·J_n = qR (electrons),
∇·J_p = −qR (holes). States that the divergence of current equals net
generation minus recombination.

**Debye length (L_D)**
The characteristic length over which electrostatic screening occurs.
L_D = sqrt(ε V_t / q N), where N is the doping density. Sets the
width of transition regions at heterojunctions and contacts.

**Depletion approximation**
An analytical model for the space-charge region in a pn junction that
assumes the carrier density is zero in the depletion region. Gives
simple formulas for depletion width and built-in field. Accurate to a
few percent for silicon at moderate doping.

**Depletion region**
The region near a pn junction or MOS interface that is depleted of
mobile carriers due to the internal electric field. Contains ionized
dopant atoms whose charge gives rise to the space-charge field.

**Donor (N_D)**
An n-type dopant atom (e.g., phosphorus, arsenic in silicon) that
contributes a free electron to the conduction band when ionized.

**Drift current**
The component of carrier current due to carrier motion in an electric
field: J_drift = qμnE (electrons). Proportional to mobility and field.

**Drift-diffusion (DD) model**
The standard macroscopic model for semiconductor carrier transport:
current = drift + diffusion, coupled to Poisson's equation. Valid for
device feature sizes above ~50 nm where quantum and ballistic effects
are negligible.

**Diffusion current**
The component of carrier current due to concentration gradients:
J_diff = qD ∇n (electrons). Proportional to diffusivity and gradient.

**Diffusivity (D)**
Carrier diffusion coefficient. Linked to mobility by the Einstein
relation: D = μ V_t.

**Einstein relation**
D_n = μ_n k_B T / q = μ_n V_t. Links the diffusion coefficient to
mobility via the thermal voltage. Follows from requiring that drift
and diffusion balance at equilibrium.

**Electric field (E)**
E = −∇ψ. Points from high to low potential; drives carrier drift.

**Electrostatic potential (ψ)**
Scalar potential defined so that E = −∇ψ. One of the three primary
unknowns in the drift-diffusion system.

**Fermi-Dirac statistics**
The exact quantum-mechanical distribution function for electrons.
Required when the semiconductor is degenerate (heavily doped,
> ~10¹⁹ cm⁻³ in silicon). Planned in SemiSim M16.4.

**Fermi level (E_F)**
The chemical potential of electrons. At thermal equilibrium, a single
E_F characterizes the entire device. Under bias, splits into quasi-
Fermi levels Φ_n, Φ_p for electrons and holes.

**Generation-recombination (G-R)**
Processes that create (generation) or annihilate (recombination)
electron-hole pairs. Principal mechanism in SemiSim: SRH.

**Hole**
A missing electron in the valence band, behaving as a positive charge
carrier. Holes have their own effective mass and mobility.

**Intrinsic carrier density (n_i)**
The carrier density in undoped silicon at thermal equilibrium:
n_i ≈ 1.5 × 10¹⁰ cm⁻³ at 300 K. Sets the scale for minority carrier
injection.

**Mass-action law**
At thermal equilibrium, n·p = n_i² everywhere in the semiconductor,
independent of doping. A fundamental constraint verified in the
`pn_1d` benchmark.

**Minority carrier**
Electrons in p-type material, or holes in n-type material. Their
concentration is ~n_i²/N. Minority carrier injection under forward
bias drives diode current.

**Mobility (μ)**
The proportionality constant between carrier drift velocity and
electric field: v_d = μE. Units: cm² V⁻¹ s⁻¹. In v0.16, constant
(Sze defaults). Field-dependent Caughey-Thomas mobility planned M16.1.

**n-type semiconductor**
Silicon doped with donor atoms. Majority carriers are electrons.

**Ohmic contact**
A metal-semiconductor contact that maintains charge neutrality and
thermal equilibrium at the contact surface. Implemented in SemiSim
as Dirichlet BCs on ψ, Φ_n, and Φ_p equal to the applied voltage plus
the built-in offset.

**p-type semiconductor**
Silicon doped with acceptor atoms. Majority carriers are holes.

**Permittivity (ε)**
ε = ε_0 ε_r. Describes how strongly a material polarizes in an
electric field. ε_r(Si) ≈ 11.7; ε_r(SiO₂) ≈ 3.9.

**pn junction**
The interface between p-type and n-type semiconductor. The simplest
and most fundamental device; the basis of diodes, BJTs, and MOSFETs.

**Poisson's equation**
−∇·(ε∇ψ) = ρ = q(p − n + N_D − N_A). One of the three governing PDEs.
Relates electrostatic potential to charge density.

**Quasi-Fermi potential (Φ_n, Φ_p)**
Generalizations of the equilibrium Fermi level to non-equilibrium
conditions. Electrons follow Φ_n; holes follow Φ_p. The Slotboom
transformation uses these as primary unknowns.

**Recombination (R)**
Net rate at which electron-hole pairs annihilate. In SemiSim:
SRH recombination through mid-gap traps. R > 0 near forward-biased
junctions; R < 0 (generation) in reverse-biased depletion regions.

**SRH recombination (Shockley-Read-Hall)**
Indirect recombination mediated by deep trap levels (impurities or
defects). Rate R = (np − n_i²) / [τ_p(n + n_1) + τ_n(p + p_1)].
Dominant recombination mechanism in indirect-gap semiconductors such
as silicon.

**Space-charge region**
See *depletion region*. Region containing net fixed ionized dopant
charge and depleted of mobile carriers.

**Thermal voltage (V_t)**
V_t = k_B T / q ≈ 25.85 mV at 300 K. The natural voltage scale for
semiconductor physics. Used as the potential normalization unit in
SemiSim's nondimensionalization.

---

## 2. Numerical Methods Terms

**Adaptive continuation**
A solver strategy for nonlinear PDEs: incrementally walk from a
known solution (e.g., V = 0) to the target (e.g., V = 0.6 V) in
small steps, using each converged solution as the initial guess for
the next. Step size is halved on failure and grown on easy convergence.
Implemented in `semi/continuation.py`.

**BDF (Backward Differentiation Formula)**
A family of implicit multistep methods for time integration of ODEs
and DAEs. BDF1 = backward Euler; BDF2 uses the two most recent time
levels. Unconditionally A-stable up to order 2. Used in SemiSim
for transient drift-diffusion (M13).

**BDF1 (Backward Euler)**
First-order BDF: (u^{n+1} − u^n) / Δt = f(u^{n+1}). Unconditionally
stable, first-order accurate in time.

**BDF2**
Second-order BDF: (3u^{n+1} − 4u^n + u^{n-1}) / (2Δt) = f(u^{n+1}).
Unconditionally stable, second-order accurate in time.

**Condition number**
For a matrix A: κ(A) = ‖A‖ · ‖A⁻¹‖. Measures sensitivity of the
linear system to perturbations. κ ~ 10³⁰ for raw DD Jacobian in SI
units; reduced to ~1–100 after nondimensionalization.

**Continuation**
See *adaptive continuation*.

**Convergence order (rate)**
The exponent p such that ‖error‖ ~ h^p as mesh spacing h → 0. Second-
order (p = 2) in L² norm is the theoretical rate for P1 Lagrange FEM
on smooth problems. Verified by the MMS suite in SemiSim.

**Direct solver**
A linear-algebra algorithm that computes an exact factorization of the
system matrix (e.g., LU decomposition). MUMPS is the direct solver in
SemiSim. Exact but memory-intensive for large 3D problems.

**DOF (degree of freedom)**
A scalar unknown in the discretized PDE system. For P1 elements, one
DOF per mesh vertex per field. Total DOFs ≈ N_vertices × N_fields.

**FEM (Finite Element Method)**
A numerical PDE discretization based on a weak (variational) form.
The domain is divided into elements; the solution is approximated by
piecewise polynomial basis functions. SemiSim uses FEniCSx/dolfinx
as its FEM backend.

**Galerkin FEM**
The standard FEM formulation where test functions are from the same
space as trial functions. Galerkin FEM on the Slotboom form of the
drift-diffusion equations is stable and gives second-order convergence
without stabilization.

**H¹ norm (H^1, H-one)**
Sobolev norm including both the L² norm of the function and its
gradient: ‖u‖²_{H¹} = ‖u‖²_{L²} + ‖∇u‖²_{L²}. Theoretical FEM
convergence for the gradient (potential, field) is first-order in H¹.

**Iterative solver**
A linear algebra algorithm that produces successive approximations to
the exact solution (e.g., GMRES, CG). Requires a preconditioner for
efficiency on ill-conditioned systems.

**Jacobian**
The matrix of partial derivatives ∂F_i/∂u_j of the nonlinear residual
F with respect to the unknowns u. Newton's method requires a Jacobian
solve at each iteration.

**L² norm (L-two)**
The root-mean-square error over the domain: ‖e‖_{L²} = sqrt(∫e²dΩ).
Theoretical FEM convergence for the solution itself is second-order in
L² for smooth problems on quasi-uniform meshes.

**MMS (Method of Manufactured Solutions)**
A V&V technique: choose an arbitrary "manufactured" exact solution
u_exact, compute the body force f that makes it satisfy the PDE
exactly, solve the modified problem numerically, and measure the error
‖u_h − u_exact‖ as the mesh is refined. The convergence rate (p ≈ 2
for P1 FEM) confirms correct implementation. MMS catches sign errors,
coefficient errors, and coupling bugs that single-point benchmarks miss.

**Multigrid**
An iterative solver that uses a hierarchy of coarser meshes to
efficiently damp errors at all spatial scales. Algebraic multigrid
(AMG) builds the hierarchy algebraically from the matrix. Used in
hypre BoomerAMG (GPU path).

**MUMPS**
MUltifrontal Massively Parallel sparse direct Solver. The default
linear solver backend in SemiSim for CPU solves. Exact LU
factorization; robust on ill-conditioned systems.

**Newton's method (Newton-Raphson)**
An iterative algorithm for solving nonlinear systems F(u) = 0 by
repeated linearization: u^{k+1} = u^k − J^{−1} F(u^k). Quadratic
convergence near the root; requires a good initial guess.

**Nondimensionalization**
The process of scaling physical variables so that they are O(1) in
the computation. In SemiSim: potentials scaled by V_t, densities
by C_0 (peak doping), spatial coordinates kept in meters. Reduces
Jacobian condition number from ~10³⁰ to manageable levels.

**Péclet number (local)**
For convection-diffusion problems: Pe = |velocity| · h / diffusivity.
Pe >> 1 means convection dominates and standard Galerkin FEM produces
spurious oscillations. The Slotboom transformation eliminates
high-Péclet stiffness in the DD continuity equations.

**PETSc**
Portable, Extensible Toolkit for Scientific Computation. A library of
data structures and solvers for large-scale linear and nonlinear
algebra. SemiSim uses PETSc via petsc4py for Newton (SNES),
linear solvers (KSP), and GPU backends.

**P1 Lagrange elements**
Piecewise linear (first-degree polynomial) finite elements. One DOF
per vertex; continuous across element boundaries. Standard choice for
Poisson and drift-diffusion in FEniCSx.

**Preconditioner**
A matrix M ≈ A^{−1} used to improve the convergence of iterative
linear solvers. AMGX and hypre BoomerAMG are AMG preconditioners used
in the SemiSim GPU path.

**Residual**
The vector F(u) measuring how far the current solution u is from
satisfying the PDE discretization F(u) = 0. Newton converges when
‖F‖ < tolerance.

**Scharfetter-Gummel (SG) discretization**
A finite-volume edge-based discretization for the drift-diffusion
current that is exact for exponential profiles in 1D. Guarantees
positive carrier densities. Used in COMSOL and many TCAD codes.
SemiSim uses the Slotboom transformation instead, which achieves
equivalent stability on FEM meshes with higher-order convergence.

**Singular perturbation**
A PDE whose small parameter multiplies the highest derivative, causing
the solution to have thin boundary layers. The Debye-length parameter
λ² in the scaled Poisson equation is a singular perturbation;
physically this is why space-charge layers are narrow.

**Slotboom transformation**
Change of primary unknowns from (ψ, n, p) to (ψ, Φ_n, Φ_p). Rewrites
drift-diffusion currents as pure gradients, making the weak form
coercive and Galerkin-stable without upwinding. Locked in ADR 0004.

**SNES (Scalable Nonlinear Equations Solver)**
The PETSc component for Newton-type nonlinear solvers. SemiSim
uses SNES for the outer Newton loop and KSP for the inner linear
solves.

**SUPG (Streamline Upwind Petrov-Galerkin)**
A stabilization technique for convection-dominated PDEs that adds
artificial diffusion along streamlines. Not used in SemiSim —
the Slotboom transformation eliminates the need for stabilization.

**UFL (Unified Form Language)**
A domain-specific language for expressing weak forms of PDEs in
FEniCSx. UFL expressions are compiled to efficient C++ kernels by
FFCx (FEniCSx Form Compiler). Weak forms in `semi/physics/` are
written in UFL.

**Weak form**
The variational (integral) form of a PDE obtained by multiplying by
a test function and integrating by parts. The FEM is based on the
weak form; Neumann boundary conditions appear naturally as boundary
integrals.

---

## 3. Device and Component Terms

**BJT (Bipolar Junction Transistor)**
A three-terminal device using two pn junctions. Not modeled in
v0.16 SemiSim.

**Contact**
A boundary of the device where an electrical terminal is attached.
In SemiSim: ohmic contacts (Dirichlet BC on all fields) or
gate contacts (Dirichlet on ψ only).

**C-V (Capacitance-Voltage)**
A characterization technique that measures the small-signal
capacitance of a device as a function of DC bias. Key measurement for
MOS capacitors. SemiSim's AC solver reproduces LF and HF C-V
curves (M14.2).

**HEMT (High Electron Mobility Transistor)**
A FET using a heterojunction to confine carriers. Requires
position-dependent band structure (planned M17).

**IV curve (Current-Voltage)**
A plot of terminal current vs. applied voltage. The primary output of
a bias sweep solve. Compared against the Shockley equation for pn
diode verification.

**MOS (Metal-Oxide-Semiconductor)**
A structure consisting of a metal gate, an insulating oxide, and a
semiconductor body. The basis of MOSFETs and MOS capacitors.

**MOSCAP (MOS Capacitor)**
A two-terminal MOS structure used to characterize the oxide and
semiconductor interface. Produces C-V curves with distinct
accumulation, depletion, and inversion regimes.

**MOSFET (Metal-Oxide-Semiconductor Field-Effect Transistor)**
A four-terminal transistor controlled by the gate voltage. The
dominant device in modern CMOS technology. Modeled in SemiSim as
a 2D multi-region structure.

**n+ region**
Heavily doped (>10¹⁸ cm⁻³) n-type semiconductor. Forms the source and
drain extensions of a MOSFET.

**Oxide (SiO₂)**
Silicon dioxide, the gate insulator in silicon MOSFETs. ε_r ≈ 3.9,
E_g ≈ 9 eV. Modeled as a passive dielectric region in SemiSim.

**pn diode**
A two-terminal device with a pn junction. The simplest device in
SemiSim; used for bias-sweep and IV verification.

**TCAD (Technology Computer-Aided Design)**
The field of simulation tools used to model semiconductor fabrication
processes and device behavior. Sentaurus, Silvaco, COMSOL Semiconductor
Module are commercial TCAD tools; SemiSim is an open-source TCAD
solver.

**Threshold voltage (V_T)**
The gate voltage at which the surface of the MOSFET channel transitions
from depletion to strong inversion. Used as a reference point in the
Pao-Sah linear-regime verification for the 2D MOSFET.

---

## 4. Software and Architecture Terms

**ADR (Architecture Decision Record)**
A short document (200–300 words) capturing a consequential design
decision, its context, and its consequences. Immutable once accepted;
superseded by a new ADR, not edited. ADRs live in `docs/adr/`.

**AMGX**
NVIDIA's Algebraic Multigrid solver library for GPU-accelerated sparse
linear systems. One of the GPU linear-solver backends in SemiSim.
Used via PETSc-CUDA.

**API (Application Programming Interface)**
A defined interface for programmatic interaction. SemiSim exposes
a Python API (direct function calls) and an HTTP API (REST endpoints).

**BDF (see Numerical Methods section)**

**CI (Continuous Integration)**
Automated build, test, and verification runs triggered on every git
push. SemiSim's CI runs on GitHub Actions and enforces the 95%
coverage gate and the full V&V suite.

**CUDA**
NVIDIA's parallel computing platform and programming model for GPU
computation. PETSc-CUDA enables GPU linear solves in SemiSim.

**DOF (see Numerical Methods section)**

**dolfinx**
The computational backend of FEniCSx. Provides the mesh, function
space, assembly, and linear algebra abstractions used in SemiSim.
Layer 4 of the five-layer architecture.

**FEniCSx**
An open-source finite-element computing platform. Consists of UFL
(form language), FFCx (form compiler), basix (element library), and
dolfinx (computational backend). SemiSim is built on FEniCSx 0.10.

**GPU (Graphics Processing Unit)**
A massively parallel processor. Used in SemiSim for accelerated
sparse linear solves (M15). Accessed via PETSc-CUDA or PETSc-HIP.

**HIP**
AMD's GPU computing framework (analogous to CUDA). PETSc-HIP enables
GPU linear solves on AMD hardware.

**HTTP (Hypertext Transfer Protocol)**
The protocol used by the SemiSim HTTP API (`semisim_server/`).
REST endpoints: `POST /solve`, `GET /runs/{id}`, `GET /capabilities`.

**hypre**
A library of parallel high-performance preconditioners, including
BoomerAMG (algebraic multigrid). Used as the GPU linear-solver backend
in SemiSim when PETSc-HIP/CUDA is available.

**JSON (JavaScript Object Notation)**
The only supported input format for SemiSim. Text-based, human-
readable, schema-validatable. Locked as the input contract in ADR 0001.

**JSON Schema**
A vocabulary for describing and validating JSON documents. SemiSim
uses JSON Schema Draft-07. The strict schema (`additionalProperties:
false`) rejects unknown fields at the boundary.

**KSP (Krylov Subspace Solver)**
The PETSc component for linear iterative solvers (GMRES, CG, BiCGSTAB,
etc.) and direct solvers (MUMPS). Used for inner linear solves in
Newton.

**MIT license**
A permissive open-source license. Users may use, modify, and
redistribute the code, including in commercial products, without
royalties.

**MPI (Message Passing Interface)**
A standard for distributed-memory parallel computing. FEniCSx is
MPI-parallel. MPI-parallel runs planned in SemiSim M19.1.

**MUMPS (see Numerical Methods section)**

**PETSc (see Numerical Methods section)**

**REST (Representational State Transfer)**
An architectural style for web APIs using standard HTTP verbs (GET,
POST). The SemiSim HTTP API follows REST conventions.

**run_id**
A UUID assigned to each simulation run by the SemiSim server.
Results are stored under `runs/<run_id>/`. Enables parallel runs and
unambiguous artifact tracking.

**Schema versioning**
The practice of giving the JSON input schema an explicit version
string. SemiSim v2.0.0 schema is strict; the v1.x legacy schema
is accepted with a deprecation warning.

**UFL (see Numerical Methods section)**

**V&V (Verification and Validation)**
- **Verification**: Are we solving the equations correctly? (MMS,
  convergence rates, conservation checks)
- **Validation**: Are we solving the right equations? (comparison to
  analytical formulas, experimental data)
Both are required and automated in SemiSim's test suite.

**WebSocket**
A full-duplex communication protocol over a single TCP connection.
Used by the SemiSim HTTP server for streaming solve progress
updates to a client UI.

**XDMF (eXtensible Data Model and Format)**
A lightweight XML + HDF5 format for scientific mesh and field data.
SemiSim writes all field output (potential, carrier densities) in
XDMF. Readable by ParaView, VisIt, and any HDF5 tool.

---

## 5. Alphabetical Quick-Reference Index

| Term / Acronym | Category | Section |
|---|---|---|
| AC (small-signal) | Physics | §5.2 |
| ADR | Architecture | §6.2 |
| AMGX | Software | §5.3, §8.2 |
| API | Software | §4.3 |
| BDF / BDF1 / BDF2 | Numerical | §5.2, §3.3 |
| BJT | Device | §8.1 |
| Boltzmann statistics | Physics | §2.1, §3.2 |
| Built-in potential | Physics | §2.3 |
| C-V curve | Device | §5.2, §8.1 |
| CI | Software | §6.1, §8.6 |
| Condition number | Numerical | §3.1 |
| Continuation | Numerical | §3.3 |
| CUDA | Software | §5.3, §8.2 |
| DD (drift-diffusion) | Physics | §2.1–2.3 |
| Debye length | Physics | §3.1, §7.4 |
| Depletion approximation | Physics | §2.3, §5.1 |
| Depletion region | Physics | §2.3 |
| DOF | Numerical | §5.3 |
| dolfinx | Software | §4.1 |
| Einstein relation | Physics | §2.1 |
| FEM | Numerical | §3.2 |
| FEniCSx | Software | §1.2, §4 |
| Fermi-Dirac | Physics | §7.1, §8.1 |
| Galerkin | Numerical | §3.2 |
| GPU | Software | §5.3, §8.2 |
| HEMT | Device | §7.2, §8.1 |
| HIP | Software | §8.2 |
| HTTP | Software | §4.3 |
| hypre BoomerAMG | Software | §5.3, §8.2 |
| IV curve | Device | §5.1, §2.3 |
| Jacobian | Numerical | §3.1 |
| JSON | Software | §1.2, §4.2 |
| JSON Schema | Software | §4.2 |
| KSP | Software | §8.2 |
| L² norm | Numerical | §3.2, §5.1 |
| Mass-action law | Physics | §2.3 |
| MMS | Numerical | §6.3, §8.6 |
| MOSCAP | Device | §5.2 |
| MOSFET | Device | §5.1 |
| MOS | Device | §2.2 |
| MPI | Software | §7.2 |
| MUMPS | Numerical | §4.3, §8.2 |
| Newton's method | Numerical | §3.1, §3.3 |
| n_i (intrinsic density) | Physics | §2.3 |
| Nondimensionalization | Numerical | §3.1 |
| Ohmic contact | Physics | §2.2 |
| Péclet number | Numerical | §3.2 |
| PETSc | Software | §4.1, §8.2 |
| P1 Lagrange | Numerical | §3.2 |
| pn junction | Device | §2.3 |
| Poisson equation | Physics | §2.1 |
| Quasi-Fermi potential | Physics | §2.2, §3.2 |
| Recombination | Physics | §2.1, §2.3 |
| REST | Software | §4.3 |
| run_id | Software | §4.3 |
| Scharfetter-Gummel | Numerical | §7.4, §8.2 |
| Singular perturbation | Numerical | §3.1 |
| Slotboom | Numerical | §3.2, §8.2 |
| SNES | Software | §8.2 |
| SRH | Physics | §5.1, §8.1 |
| SUPG | Numerical | §3.2 |
| TCAD | Device | §1.1 |
| Thermal voltage V_t | Physics | §3.1 |
| Threshold voltage V_T | Device | §5.1 |
| UFL | Software | §4.1 |
| V&V | Numerical | §6.3, §8.6 |
| V_bi | Physics | §2.3 |
| WebSocket | Software | §4.3 |
| XDMF | Software | §5.3 |
