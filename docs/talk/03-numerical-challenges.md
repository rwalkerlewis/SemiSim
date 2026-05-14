# §3 — Numerical Challenges: Why Naïve FEM Fails

**Suggested slide title:** "Three Reasons a Naïve Implementation Diverges Immediately"
**Target time:** 5–6 minutes

---

## Slide 3.1 — Problem 1: The Numbers Span 40 Orders of Magnitude

**[8:30]**

Take a typical 1 µm device doped at 10¹⁷ cm⁻³. The quantities that go
into the Jacobian of Poisson's equation are:

| Quantity | Typical SI value |
|---|---|
| Permittivity ε | ~10⁻¹¹ F/m |
| Elementary charge q | ~10⁻¹⁹ C |
| Carrier density N | ~10²³ m⁻³ |
| Device length L | ~10⁻⁶ m |

Plugging those into the stiffness matrix gives entries that span from
10⁻¹¹ to 10²³. The condition number of Newton's Jacobian is around 10³⁰.
Floating-point arithmetic breaks down. Newton diverges or stagnates
immediately.

The fix is nondimensionalization. Scale every variable so the scaled
versions are all O(1):

- Potentials in units of the thermal voltage V_t = kT/q ≈ 25.85 mV
- Densities in units of C₀, the peak doping (typically ~10¹⁷ cm⁻³)
- Spatial coordinates stay in meters — this is a deliberate asymmetric
  choice so users can specify geometry in SI in the JSON file

After scaling, the only "small" coefficient left is λ² ≈ 10⁻⁴ — the
squared ratio of the Debye length to the device length. That small
parameter is not a problem; it is the physics. It is exactly why
space-charge regions (depletion zones) are narrow compared to the device.

**Key points**
- Raw Jacobian condition number ~10³⁰ → Newton fails immediately.
- Fix: scale ψ by V_t, densities by C₀, keep x in meters.
- After scaling, λ² ≈ 10⁻⁴ is a singular perturbation that describes depletion physics.
- Code: `semi/scaling.py`; locked in ADR 0002.

---

## Slide 3.2 — Problem 2: Carrier Densities Are Exponential in the Potential

**[11:00]**

Even after scaling, there is a second problem. Electron density n varies
as exp(ψ/V_t). At 1 V forward bias, that is a factor of e⁴⁰ ≈ 10¹⁷ over
equilibrium. The Newton Jacobian has entries that vary by the same ratio.
Standard Galerkin FEM with raw (ψ, n, p) as unknowns is oscillatory in
the drift-dominated regime — the local Péclet number is >> 1 nearly
everywhere in a real device.

The standard fix in device simulation is the **Slotboom transformation**.
Instead of solving for n and p directly, solve for the quasi-Fermi
potentials Φ_n and Φ_p, defined by:

```
n = n_i · exp((ψ − Φ_n) / V_t)
p = n_i · exp((Φ_p − ψ) / V_t)
```

When you rewrite the continuity equations in terms of Φ_n and Φ_p, the
currents become pure gradients:

```
J_n = −q μ_n n ∇Φ_n
J_p = −q μ_p p ∇Φ_p
```

A pure gradient of a scalar: this is what standard Galerkin FEM handles
well. No upwind stabilization, no SUPG parameters, no mesh-dependent
tuning. The weak form is coercive. MMS confirms second-order L² convergence.

The boundary conditions on Φ_n and Φ_p at ohmic contacts are also natural:
they simply equal the applied voltage. At equilibrium, Φ_n = Φ_p = 0
everywhere, and n·p = n_i² (the mass-action law) is satisfied identically.

**Key points**
- Raw drift-diffusion: Galerkin unstable at high Péclet (everywhere in real devices).
- Slotboom transform: (ψ, n, p) → (ψ, Φ_n, Φ_p); currents become pure gradients.
- No stabilization parameters to tune; FEM coercive; second-order convergence confirmed by MMS.
- Code: `semi/physics/drift_diffusion.py`; locked in ADR 0004.

---

## Slide 3.3 — Problem 3: Forward Bias Breaks Newton From a Cold Start

**[13:30]**

There is a third problem specific to biased solves. At 0.6 V forward bias
on a silicon pn diode, minority carrier injection raises n and p by a
factor of exp(0.6 / 0.0258) ≈ 10¹⁰ over equilibrium. If you try to
Newton-solve for V = 0.6 V starting from the V = 0 equilibrium solution,
the first Newton step overshoots by many orders of magnitude. The solver
never recovers.

The fix is **continuation**: walk up from zero voltage in small steps.
Solve at V = 0, then at V = ΔV, using the previous solution as the
initial guess. If a step fails to converge, halve it and retry. If three
consecutive steps converge in very few Newton iterations, grow the step.
This is adaptive arc-length homotopy.

In SemiSim this logic lives in `semi/continuation.py` as
`AdaptiveStepController`. The bias sweep runner in
`semi/runners/bias_sweep.py` drives the controller, records an IV pair
at each successful step, and stops at the target voltage.

The interplay of continuation and Slotboom is the reason the solver can
reach 0.6 V forward bias on a silicon diode and recover IV curves that
match the Shockley equation within 10%. Either ingredient alone is
insufficient.

**Key points**
- At 0.6 V forward, carrier injection is 10¹⁰× equilibrium — Newton can't jump there.
- Fix: adaptive continuation; halve step on failure, grow step on easy convergence.
- Code: `semi/continuation.py`, `semi/runners/bias_sweep.py`.
- These three fixes — scaling, Slotboom, continuation — are the minimum viable solver.

**Transition:** With those three problems solved, the code has a stable
foundation. Now let me show you how it is organized.
