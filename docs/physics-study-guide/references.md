# References

A bibliography of the textbooks, original papers, and software docs
that this study guide draws from. Each chapter's *Further reading*
section points into this list with chapter/section/equation precision
where possible.

## Textbooks

### Semiconductor device physics

- **Sze, S. M., and Ng, K. K. (2007).** *Physics of Semiconductor Devices*,
  3rd ed. Wiley-Interscience.
  The standard reference. Chapter 1 (carriers, statistics), Chapter 2
  (pn junctions), Chapter 3 (Schottky contacts), Chapter 4 (MOSCAP),
  Chapter 6 (MOSFET), Appendix C (material parameters). The numerical
  values in [`semi/materials.py`](../../semi/materials.py) come from
  this book's Table 7.

- **Sze, S. M. (1981).** *Physics of Semiconductor Devices*, 2nd ed.
  Wiley. The earlier edition is still cited for some derivations.

- **Pierret, R. F. (1996).** *Semiconductor Device Fundamentals*.
  Addison-Wesley. More approachable than Sze; chapter 3 for
  doping/charge-neutrality, chapter 5 for pn junctions. Modular Vols.
  I-VI cover individual topics in more depth.

- **Pierret, R. F. (2002).** *Advanced Semiconductor Fundamentals*,
  2nd ed. Prentice Hall. Quantum-mechanical foundations of band theory.

- **Hu, C. (2010).** *Modern Semiconductor Devices for Integrated
  Circuits*. Pearson. Chapter 5 is the C-V chapter that the kronos-semi
  M14.2 axisymmetric MOSCAP benchmark reproduces (Hu Fig. 5-18).
  Chapter 4 covers the pn-junction picture used by the M3 reverse-bias
  verifier.

- **Tsividis, Y. (1999).** *Operation and Modeling of the MOS
  Transistor*, 2nd ed. McGraw-Hill. The MOSFET-modeling reference;
  detailed compact-model derivations.

- **Nicollian, E. H., and Brews, J. R. (1982).** *MOS (Metal Oxide
  Semiconductor) Physics and Technology*. Wiley. The C-V bible.

- **Ashcroft, N. W., and Mermin, N. D. (1976).** *Solid State Physics*.
  Saunders. Chapter 8 (Bloch theorem), Chapter 12 (semiconductor band
  structure).

- **Kittel, C. (2004).** *Introduction to Solid State Physics*, 8th ed.
  Wiley. Chapters 7-8 for band theory at undergraduate level.

### Computational and numerical methods

- **Selberherr, S. (1984).** *Analysis and Simulation of Semiconductor
  Devices*. Springer. The classic device-simulation reference.
  Chapter 3 derives DD from the BTE moments; Chapter 5 §5.4 covers
  Slotboom variables; Chapter 7 covers nondimensionalization
  (matches kronos-semi); Chapter 8 covers nonlinear solvers and
  continuation.

- **Vasileska, D., Goodnick, S. M., and Klimeck, G. (2010).**
  *Computational Electronics*. CRC Press. Chapter 1 (BTE → DD),
  Chapter 3 (transient methods), Chapter 4 (Slotboom vs Scharfetter-
  Gummel), Chapter 5 (Caughey-Thomas, Lombardi mobility).

- **Markowich, P. A. (1986).** *The Stationary Semiconductor Device
  Equations*. Springer. Mathematical analysis: singular perturbation,
  depletion approximation, well-posedness.

- **Brezzi, F., and Boffi, D. (2003).** *Mixed and Hybrid Finite
  Element Methods*. Springer. For mixed-formulation FEM; relevant
  background even though kronos-semi uses pure $H^1$-conforming
  Lagrange.

- **Brenner, S. C., and Scott, L. R. (2008).** *The Mathematical
  Theory of Finite Element Methods*, 3rd ed. Springer. Rigorous
  analysis of FE convergence rates: the source for the
  $L^2 = O(h^2)$, $H^1 = O(h)$ rates that the MMS gates target.

- **Logg, A., Mardal, K.-A., and Wells, G. N. (2012).** *Automated
  Solution of Differential Equations by the Finite Element Method*.
  Springer. The FEniCS book. Free online. Chapter 1 is enough for
  what kronos-semi does.

- **Hairer, E., Nørsett, S. P., and Wanner, G. (1996).** *Solving
  Ordinary Differential Equations II: Stiff and Differential-
  Algebraic Problems*, 2nd ed. Springer. The reference for BDF
  stability theory used by Ch. 17.

- **Saad, Y. (2003).** *Iterative Methods for Sparse Linear Systems*,
  2nd ed. SIAM. Free online. Krylov methods, AMG; relevant to
  Ch. 19 GPU backends.

- **Trottenberg, U., Oosterlee, C. W., and Schuller, A. (2001).**
  *Multigrid*. Academic Press. The classic AMG reference.

- **Allgower, E. L., and Georg, K. (2003).** *Numerical Continuation
  Methods*. SIAM. Continuation strategies; relevant to the engine's
  bias-continuation logic.

- **Nocedal, J., and Wright, S. J. (2006).** *Numerical Optimization*,
  2nd ed. Springer. Chapter 11 covers Newton with line-search.

- **Dennis, J. E., and Schnabel, R. B. (1996).** *Numerical Methods
  for Unconstrained Optimization and Nonlinear Equations*. SIAM
  Classics in Applied Mathematics.

### Verification and validation

- **Roache, P. J. (1998).** *Verification and Validation in
  Computational Science and Engineering*. Hermosa.

- **Oberkampf, W. L., and Roy, C. J. (2010).** *Verification and
  Validation in Scientific Computing*. Cambridge University Press.

- **Salari, K., and Knupp, P. (2000).** *Code Verification by the
  Method of Manufactured Solutions*. Sandia report SAND2000-1444.

### Electromagnetism (prerequisite background)

- **Griffiths, D. J. (2013).** *Introduction to Electrodynamics*,
  4th ed. Pearson. Standard undergraduate text; Ch. 1 prerequisite.

- **Jackson, J. D. (1998).** *Classical Electrodynamics*, 3rd ed.
  Wiley. Graduate-level reference.

### Process technology

- **Plummer, J. D., Deal, M. D., and Griffin, P. B. (2000).** *Silicon
  VLSI Technology*. Prentice Hall. For Gaussian-implant rationale and
  the diffusion-driven profiles the engine's `gaussian` doping mimics.

## Original papers

### Slotboom (drift-diffusion variable change)

- **Slotboom, J. W. (1973).** "Computer-aided two-dimensional analysis
  of bipolar transistors." *IEEE Trans. Electron Devices* 20, 669-679.
  The original Slotboom paper.

### Scharfetter-Gummel (alternative discretization)

- **Scharfetter, D. L., and Gummel, H. K. (1969).** "Large-signal
  analysis of a silicon Read diode oscillator." *IEEE Trans. Electron
  Devices* 16, 64-77. The original Scharfetter-Gummel discretization;
  kronos-semi's M13.1 ships Scharfetter-Gummel primitives in
  [`semi/fem/scharfetter_gummel.py`](../../semi/fem/scharfetter_gummel.py)
  for future use though Slotboom is the active path.

- **Brezzi, F., Marini, L. D., and Pietra, P. (1989).** "Two-dimensional
  exponential fitting and applications to drift-diffusion models."
  *SIAM J. Numer. Anal.* 26, 1342-1355. Mathematical justification
  for SG and Slotboom-equivalent schemes.

### Shockley-Read-Hall

- **Shockley, W., and Read, W. T., Jr. (1952).** "Statistics of the
  recombinations of holes and electrons." *Phys. Rev.* 87, 835-842.

- **Hall, R. N. (1952).** "Electron-hole recombination in germanium."
  *Phys. Rev.* 87, 387-392.

### Pao-Sah MOSFET

- **Pao, H. C., and Sah, C. T. (1966).** "Effects of diffusion current
  on characteristics of metal-oxide(insulator)-semiconductor
  transistors." *Solid State Electron.* 9, 927-937.

### Diode I-V

- **Shockley, W. (1949).** "The theory of p-n junctions in
  semiconductors and p-n junction transistors." *Bell System Tech.
  J.* 28, 435-489.

- **Sah, C. T., Noyce, R. N., and Shockley, W. (1957).** "Carrier
  generation and recombination in p-n junctions and p-n junction
  characteristics." *Proc. IRE* 45, 1228-1243. The SNS recombination
  current the engine's reverse-bias verifier uses.

### Schottky contacts

- **Crowell, C. R., and Sze, S. M. (1966).** "Current transport in
  metal-semiconductor barriers." *Solid State Electron.* 9, 1035-1048.

- **Bardeen, J. (1947).** "Surface states and rectification at a
  metal-semiconductor contact." *Phys. Rev.* 71, 717-727.

### Mobility models

- **Caughey, D. M., and Thomas, R. E. (1967).** "Carrier mobilities in
  silicon empirically related to doping and field." *Proc. IEEE* 55,
  2192-2193. The kronos-semi M16.1 closed-form.

- **Lombardi, C., Manzini, S., Saporito, A., and Vanzi, M. (1988).**
  "A physically based mobility model for numerical simulation of
  nonplanar devices." *IEEE Trans. CAD* 7, 1164-1171. Planned in
  M16.2.

- **Klaassen, D. B. M. (1992).** "A unified mobility model for device
  simulation - I." *Solid-State Electron.* 35, 953-959. Doping-
  dependent mobility model; not yet in the engine.

### Tunneling

- **Kane, E. O. (1961).** "Theory of tunneling." *J. Appl. Phys.* 32,
  83-91. The original BBT paper for M16.6.

- **Hurkx, G. A. M., Klaassen, D. B. M., and Knuvers, M. P. G. (1992).**
  "A new recombination model for device simulation including
  tunneling." *IEEE Trans. Electron Devices* 39, 331-338.
  The trap-assisted tunneling reference for M16.6.

### Intrinsic carrier density

- **Altermatt, P. P., et al. (2003).** "A simulation model for the
  density of states and for incomplete ionization in crystalline
  silicon." *J. Appl. Phys.* 93, 1598-1604. The source for the modern
  $n_i^\mathrm{Si} = 1.0\times 10^{10}\,\mathrm{cm^{-3}}$ value.

- **Misiakos, K., and Tsamakis, D. (1993).** "Accurate measurements of
  the silicon intrinsic carrier density from 78 to 340 K." *J. Appl.
  Phys.* 74, 3293-3297. Experimental basis.

### Bandgap narrowing

- **Slotboom, J. W., and de Graaff, H. C. (1976).** "Measurements of
  bandgap narrowing in Si bipolar transistors." *Solid-State
  Electron.* 19, 857-862. Heavy-doping bandgap narrowing; not yet in
  kronos-semi.

### AC small-signal analysis

- **Schenk, A. (1998).** *Advanced Physical Models for Silicon Device
  Simulation*. Springer. Chapter 5 §"Linearised continuity equations
  for AC analysis."

- **Snowden, C. (1989).** *Semiconductor Device Modelling*. Peter
  Peregrinus. §"Small-signal extraction from drift-diffusion."

### Velocity saturation in short channels

- **Sodini, C. G., Ko, P. K., and Moll, J. L. (1984).** "The effect of
  high fields on MOS device and circuit performance." *IEEE Trans.
  Electron Devices* 31, 1386-1393. Velocity-saturation in short-channel
  MOSFETs; M16.1 motivation.

## Software documentation

- **FEniCSx project** (Wells et al.). https://fenicsproject.org and
  GitHub at https://github.com/FEniCS/dolfinx. The dolfinx 0.10
  documentation is the runtime reference.

- **Dokken, J. S.** *The dolfinx tutorial.* https://jsdokken.com/dolfinx-tutorial/.
  Practical guide for dolfinx 0.10+.

- **PETSc** (Balay et al.). https://petsc.org. Including the SNES,
  KSP, PC, and Mat manual sections.

- **MUMPS** (Amestoy et al.). https://mumps-solver.org. User guide.

- **AMGX** (NVIDIA). https://github.com/NVIDIA/AMGX.

- **hypre** (Lawrence Livermore). https://hypre.readthedocs.io.

- **gmsh** (Geuzaine and Remacle). https://gmsh.info.

## Engine-internal references (this repo)

For convenience, the most-cited internal documents:

- [`docs/PHYSICS.md`](../PHYSICS.md) — definitive physics reference.
- [`docs/PHYSICS_INTRO.md`](../PHYSICS_INTRO.md) — narrative on-ramp.
- [`docs/theory/`](../theory/) — focused theory notes (slotboom, scaling, axisymmetric, moscap_cv, dolfinx_choice).
- [`docs/adr/`](../adr/) — locked architectural / numerical decisions.
- [`docs/ROADMAP.md`](../ROADMAP.md) — milestone delivery history.
- [`docs/IMPROVEMENT_GUIDE.md`](../IMPROVEMENT_GUIDE.md) — open-milestone scopes and acceptance tests.
- [`docs/ARCHITECTURE.md`](../ARCHITECTURE.md) — five-layer design.
- [`docs/WALKTHROUGH.md`](../WALKTHROUGH.md) — end-to-end `semi.run.run(cfg)` trace.
- [`docs/schema/reference.md`](../schema/reference.md) — JSON contract.
- [`docs/gpu.md`](../gpu.md) — GPU backend reference.
- [`docs/mms_dd_derivation.md`](../mms_dd_derivation.md), [`docs/mos_derivation.md`](../mos_derivation.md), [`docs/resistor_derivation.md`](../resistor_derivation.md) — derivation-first artifacts.
