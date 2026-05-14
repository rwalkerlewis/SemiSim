# §2 — The Physics: What SemiSim Actually Computes

**Suggested slide title:** "Three Equations That Describe Every Classical Semiconductor Device"
**Target time:** 5–6 minutes

---

## Slide 2.1 — The Three Equations in English

**[3:00]** This is the physics section. You don't need a semiconductor
background to follow this. You need calculus and the concept of a PDE.

Three things happen inside silicon when you apply a voltage:

1. **Charges make fields.** An electric potential ψ exists because of
   fixed ions and mobile carriers. This is described by Poisson's
   equation: the divergence of the displacement field equals the charge
   density.

2. **Charge is conserved.** Electrons and holes can't appear or disappear
   without a paired partner. The rate at which electrons flow into a
   volume equals the rate at which they accumulate there minus the rate at
   which they recombine with holes. This is the continuity equation —
   there is one for electrons and one for holes.

3. **Carriers move by drift and diffusion.** Electrons are pushed by the
   electric field (drift) and spread from regions of high concentration to
   low (diffusion). The current density is the sum of both.

That is literally the whole physics. Everything you have ever seen in a
semiconductor datasheet — diode IV curves, transistor transfer
characteristics, MOSFET saturation — falls out of those three equations
plus the right boundary conditions.

**Key points**
- Poisson: ∇·(ε∇ψ) = −ρ; charge density drives the potential.
- Continuity: ∇·J_n = +qR; ∇·J_p = −qR; carriers conserved.
- Drift-diffusion: J = drift + diffusion; Einstein relation links them.
- Three coupled nonlinear PDEs; three unknowns (ψ, n, p) at every mesh point.

**Transition:** Now here is where it gets interesting. Those equations are
well-defined physically. But if you try to solve them naively in a
finite-element code, they fail immediately. Let me show you why.

---

## Slide 2.2 — Boundary Conditions: Where Contacts Happen

**[5:30]**

A semiconductor device is a chunk of silicon with metal contacts attached
at specific faces. Each contact imposes a boundary condition on the PDEs.

**Ohmic contacts** are perfect current sources. The metal holds the
semiconductor at a fixed potential and maintains charge neutrality right
at the contact surface. In SemiSim, this means Dirichlet conditions on
all three fields: ψ = V_applied + ψ_builtin, and both quasi-Fermi
potentials equal V_applied.

**Gate contacts** in a MOS structure sit behind an insulator — silicon
dioxide. The gate holds a potential but there is no carrier exchange
across the oxide. In the code this means: Poisson is solved in the whole
domain (silicon + oxide), but the continuity equations live on a submesh
that is the semiconductor only. The oxide just has permittivity.

**Insulating boundaries** — the "sides" of the device that aren't contacts
— apply zero-flux Neumann conditions. No current leaks out the walls. This
is free in the weak form: integrating by parts and dropping the boundary
term gives you Neumann automatically.

**Key points**
- Ohmic: three Dirichlet BCs; Φ_n = Φ_p = V_applied at the contact.
- Gate: ψ Dirichlet on the oxide top face; no BC on quasi-Fermis (they don't live in oxide).
- Neumann on everything else: implicit in the weak form.

---

## Slide 2.3 — A Concrete Example: 1D pn Junction

**[7:00]** Let's walk through the simplest realistic case.

Take a 2 µm strip of silicon. The left micron is p-type (doped with
acceptors at 10¹⁷ cm⁻³). The right micron is n-type (doped with donors at
10¹⁷ cm⁻³). Contacts at both ends.

At the junction, electrons from the n-side diffuse left and recombine with
holes. This leaves behind positive ionized donors on the n-side and
negative ionized acceptors on the p-side. That charge distribution creates
a field that opposes further diffusion. Equilibrium is the point where
drift exactly cancels diffusion.

The resulting "depletion region" has a built-in potential of about 0.83 V,
a width of about 150 nm, and a peak electric field of about 113 kV/cm.
These are textbook formulas. The Poisson solver recovers them to better
than 1% in a few Newton iterations.

**Key points**
- pn junction is the simplest real device; all fundamentals are visible.
- Built-in voltage V_bi ≈ V_t · ln(N_A · N_D / n_i²) ≈ 0.83 V at symmetric 10¹⁷ doping.
- Depletion width ~150 nm; peak field ~113 kV/cm.
- Solver recovers these from first principles, not from fitted formulas.

**Transition:** Good. Now let me tell you what goes wrong if you don't
handle the numerics carefully.
