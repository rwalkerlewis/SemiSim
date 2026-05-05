# 11 — Quasi-Fermi potentials and the Slotboom transformation

## Learning objectives

- Define the electron and hole quasi-Fermi potentials $\Phi_n, \Phi_p$
  and connect them to the equilibrium Fermi level under bias.
- Derive the Slotboom expressions $n = n_i\exp((\psi-\Phi_n)/V_t)$,
  $p = n_i\exp((\Phi_p-\psi)/V_t)$.
- Show that the electron current density becomes a *pure gradient* of
  $\Phi_n$: $\mathbf{J}_n = -q\mu_n n\,\nabla\Phi_n$, and explain why
  this is what makes Galerkin discretization of the continuity equations
  well-posed.
- State the equilibrium identity $\Phi_n = \Phi_p = 0$ globally (under
  the engine's convention) and recognize how the engine pins them at
  ohmic contacts under bias.
- Locate the Slotboom helpers in `semi/physics/slotboom.py` and the
  block residual in `semi/physics/drift_diffusion.py`.

## Physical motivation

Without the Slotboom transformation, the drift-diffusion residual is
unstable under Galerkin discretization at the Péclet numbers typical
of real devices (Ch. 5 §Péclet). With the Slotboom transformation, the
residual becomes coercive and standard P1 Lagrange elements work
without any stabilization (no SUPG, no Scharfetter–Gummel edge fluxes).
The cost is a change of primary unknown — instead of solving for
$(\psi, n, p)$, the engine solves for $(\psi, \Phi_n, \Phi_p)$ —
and a chain-rule mass matrix that becomes relevant for the transient
runner (Ch. 17).

This is the most important variable change in the entire engine. It is
locked in by ADR 0004; the M13 transient runner originally tried to use
$(\psi, n, p)$ primary form (ADR 0009) and was eventually superseded
by ADR 0014 to put the transient runner on Slotboom too.

## Derivation from first principles

### Quasi-Fermi levels under bias

At thermal equilibrium, $E_F$ is a single global constant. The carrier
densities follow Boltzmann (3.4) with this single $E_F$.

Under bias, the electron and hole populations are no longer in
equilibrium with each other. Define **quasi-Fermi levels** $E_{F,n}(\mathbf{x})$
and $E_{F,p}(\mathbf{x})$ such that the carrier densities still take the
Boltzmann form, but with each carrier-type's own quasi-Fermi level:

$$
n(\mathbf{x}) = N_c\,\exp\!\left(\frac{E_{F,n}(\mathbf{x}) - E_c(\mathbf{x})}{kT}\right),
\quad
p(\mathbf{x}) = N_v\,\exp\!\left(\frac{E_v(\mathbf{x}) - E_{F,p}(\mathbf{x})}{kT}\right).
\tag{11.1}
$$

In equilibrium, both $E_{F,n}$ and $E_{F,p}$ collapse to the global
$E_F$. Under bias, they generally split: $E_{F,n} \neq E_{F,p}$ in any
region where the carriers are out of equilibrium with each other (e.g.
the depletion region of a forward-biased pn junction).

### Quasi-Fermi *potentials*

The engine works with **potentials** rather than energies, defined by

$$
\Phi_n(\mathbf{x}) \equiv -\frac{E_{F,n}(\mathbf{x})}{q},
\qquad
\Phi_p(\mathbf{x}) \equiv -\frac{E_{F,p}(\mathbf{x})}{q}.
\tag{11.2}
$$

(Negative sign so that high $\Phi_n$ corresponds to low electron density
and vice versa; matches the standard sign convention for the Fermi
*potential*.)

Substitute (3.6) — $\psi = -E_i/q$ — and rearrange (11.1) using
$E_c - E_i = E_g/2$ at the band-edge convention:

$$
\boxed{
n = n_i\,\exp\!\left(\frac{\psi - \Phi_n}{V_t}\right),
\qquad
p = n_i\,\exp\!\left(\frac{\Phi_p - \psi}{V_t}\right).
}
\tag{11.3}
$$

This is the engine's working form, in [`semi/physics/slotboom.py:27-42`](../../semi/physics/slotboom.py)
and [`semi/physics/drift_diffusion.py:144-145`](../../semi/physics/drift_diffusion.py)
(in scaled units where the $V_t$ is absorbed).

### Equilibrium identity

At thermal equilibrium $E_{F,n} = E_{F,p} = E_F$, so $\Phi_n = \Phi_p$
globally. The engine pins this single value by convention to zero —
the choice of energy zero at the intrinsic level is the convention. So:

$$
\Phi_n = \Phi_p = 0\quad\text{at thermal equilibrium}.
\tag{11.4}
$$

Substituting (11.4) into (11.3): $n = n_i\exp(\psi/V_t)$, $p = n_i\exp(-\psi/V_t)$,
which is the equilibrium expression used in [`semi/physics/poisson.py:73`](../../semi/physics/poisson.py).
Mass action $np = n_i^2$ falls out trivially from (11.3) when
$\Phi_n = \Phi_p$.

### Under bias: ohmic contact pinning

At an ohmic contact under applied bias $V_\mathrm{applied}$, the carriers
are at equilibrium with the contact metal, but the contact metal's Fermi
level is at $-qV_\mathrm{applied}$ (negative because $\Phi$ is the negative
Fermi *potential*). So both $\Phi_n$ and $\Phi_p$ at the contact equal
the applied bias:

$$
\Phi_n|_\mathrm{contact} = \Phi_p|_\mathrm{contact} = V_\mathrm{applied}.
\tag{11.5}
$$

This is the **Shockley boundary condition** (cf. (8.5)). Inside the
device, $\Phi_n$ and $\Phi_p$ vary smoothly between the contacts; inside
a depletion region they can split by up to $V_\mathrm{applied}$.

### The Slotboom transformation: currents become pure gradients

Substitute (11.3) into the drift-diffusion current (5.5) for electrons:

$$
\mathbf{J}_n = q\mu_n n\mathbf{E} + qD_n\nabla n
   = -q\mu_n n\nabla\psi + q\mu_n V_t\,\nabla n.
$$

Compute $\nabla n$ from (11.3):

$$
\nabla n = n_i\,e^{(\psi - \Phi_n)/V_t}\cdot \frac{1}{V_t}(\nabla\psi - \nabla\Phi_n)
        = \frac{n}{V_t}(\nabla\psi - \nabla\Phi_n).
$$

Substitute and simplify:

$$
\mathbf{J}_n = -q\mu_n n\,\nabla\psi + q\mu_n n\,(\nabla\psi - \nabla\Phi_n)
   = -q\mu_n n\,\nabla\Phi_n.
\tag{11.6}
$$

The drift and diffusion terms have *cancelled* against each other; the
remaining structure is

$$
\boxed{
\mathbf{J}_n = -q\,\mu_n n\,\nabla\Phi_n,
\qquad
\mathbf{J}_p = -q\,\mu_p p\,\nabla\Phi_p,
}
\tag{11.7}
$$

a **pure gradient times a coefficient**. This is the magic of Slotboom.

### Why this fixes the discretization

The continuity equation $\nabla\!\cdot\!\mathbf{J}_n = qR$ with (11.7)
becomes

$$
-\nabla\!\cdot\!(q\mu_n n\,\nabla\Phi_n) = qR,
\tag{11.8}
$$

a second-order *elliptic* PDE in $\Phi_n$ with a positive coefficient
$q\mu_n n > 0$. Galerkin Lagrange discretization of (11.8) is *coercive*
in $H^1$ — there is no Péclet number to worry about, no upwind needed,
no instability under refinement. Sign-symmetry and coupling structure
of the resulting block residual are documented in
[`docs/theory/slotboom.md`](../theory/slotboom.md).

The price paid: the carrier density $n$ inside the coefficient is itself
nonlinear in the unknowns ($n = n_i\exp((\psi - \Phi_n)/V_t)$), so the
discrete system is still nonlinear. Newton handles this. The
discretization is now well-posed at every iterate.

## Key results

- Slotboom expressions: (11.3).
- Equilibrium $\Phi_n = \Phi_p = 0$: (11.4).
- Shockley boundary at ohmic contact: (11.5).
- Pure-gradient current form: (11.7).
- Coercive elliptic continuity: (11.8).

## Worked numerical example

**M2 forward bias.** $V_\mathrm{applied} = 0.6\,\mathrm{V}$ on the cathode.
At the cathode (n-side ohmic, n-type):
- $\psi_R = \psi_R^\mathrm{eq} + V_\mathrm{applied} = 0.4167 + 0.6 = 1.0167\,\mathrm{V}$
- $\Phi_n = \Phi_p = 0.6\,\mathrm{V}$
- $n_R = n_i\exp((\psi_R - \Phi_n)/V_t) = 10^{16}\exp((1.0167 - 0.6)/0.02585)
       = 10^{16}\exp(16.118) = 10^{16}\cdot 10^7 = 10^{23}\,\mathrm{m^{-3}}$
       $= 10^{17}\,\mathrm{cm^{-3}}$ (= $N_D$, ✓).
- $p_R = n_i\exp((\Phi_p - \psi_R)/V_t) = 10^{16}\exp(-16.118) = 10^{16}\cdot 10^{-7}
       = 10^9\,\mathrm{m^{-3}} = 10^3\,\mathrm{cm^{-3}}$ (minority).

Mass action in the forward-biased ohmic contact: $n_R p_R = 10^{17}\cdot 10^3 = 10^{20}$.
But $n_i^2 = 10^{20}$. ✓ — at the contact, mass action *holds* because
$\Phi_n = \Phi_p$, even though the contact is biased.

**Inside the depletion region under forward bias.** $\Phi_n$ and $\Phi_p$
split. At the metallurgical junction with $\Phi_n - \Phi_p \approx V$,
$np = n_i^2\exp((\Phi_p - \Phi_n)/V_t) = n_i^2\exp(-V/V_t)\cdot\exp(2V/V_t) = n_i^2\exp(V/V_t)$.

Wait — let me redo the algebra. $np = n_i\exp((\psi-\Phi_n)/V_t)\cdot n_i\exp((\Phi_p - \psi)/V_t)
= n_i^2\exp((\Phi_p - \Phi_n)/V_t)$. Under bias, $\Phi_n - \Phi_p = -V$
in the depletion region (so $\Phi_p - \Phi_n = +V$); $np = n_i^2\exp(V/V_t)$.
At $V = 0.6\,\mathrm{V}$: $np/n_i^2 = e^{23.21} = 10^{10.08}$, an
enormous excess. This is the **mass-action breaking** that drives the
diffusion current.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| Slotboom UFL helpers | (11.3) | `semi/physics/slotboom.py:27-42` (`n_from_slotboom`, `p_from_slotboom`) |
| Slotboom numpy helpers | (11.3) | `semi/physics/slotboom.py:45-52` |
| Recover $\Phi_n$ from $(n, \psi)$ | (11.3) inverted | `semi/physics/slotboom.py:55-70` |
| Equilibrium $\Phi_n = \Phi_p = 0$ | (11.4) | `semi/runners/bias_sweep.py:72-73`, `transient.py:243-244` |
| Shockley boundary at contact | (11.5) | `semi/bcs.py:255-266` |
| Block residual in Slotboom form | (11.8) | `semi/physics/drift_diffusion.py:158-169` |
| Current in pure-gradient form | (11.7) | `semi/postprocess.py:97-99` |
| Sign-convention for residual | – | `docs/PHYSICS.md` §1.3, ADR 0004 |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §1.3](../PHYSICS.md) — Slotboom continuity equations as shipped.
- [`docs/PHYSICS.md` §2.5](../PHYSICS.md) — scaled DD residual.
- [`docs/theory/slotboom.md`](../theory/slotboom.md) — why Slotboom over SG, including the SG flux primitives in M13.1.
- [`docs/adr/0004-slotboom-variables-for-dd.md`](../adr/0004-slotboom-variables-for-dd.md) — the locked decision.
- [`docs/adr/0014-slotboom-transient.md`](../adr/0014-slotboom-transient.md) — transient runner alignment with the steady-state choice.

## Common pitfalls

1. **Sign convention for $\Phi_n$.** The engine uses
   $E_F = -q\Phi$, so high $n$ means low $\Phi_n$ (more negative), and
   low $n$ (e.g. minority side) means high $\Phi_n$. Some textbooks use
   the opposite sign. Match the engine's convention or your residual
   will be off.
2. **$\Phi_n \neq V_\mathrm{applied}$ inside the device.** $\Phi_n$
   varies smoothly between contacts; only at the *contact facets* is
   it pinned to the applied bias. In the bulk of a forward-biased
   diode, $\Phi_n$ is roughly constant across the n-bulk and through
   the depletion region (because $\nabla\Phi_n$ is small except in the
   depletion edge), then transitions to the p-side bias.
3. **Underflow in Slotboom densities.** $n = n_i\exp((\psi-\Phi_n)/V_t)$
   underflows to zero when the argument is below $-700/\ln(10) \approx -304$
   (about $-7.86\,\mathrm{V}$). On a 2 µm device at $V_F = 0.6\,\mathrm{V}$,
   the minority-side $\Phi_n - \psi$ can hit this floor. The engine's
   bias-continuation strategy (Ch. 16) keeps Newton iterates within the
   safe range; if you skip continuation, expect SNES failures.
4. **Equilibrium initial guess for $\Phi_n, \Phi_p$.** Setting them to
   zero is correct at $V = 0$ but seeds Newton at the wrong basin under
   bias; the bias-continuation walk is what relocates the iterates
   smoothly into the bias-converged state.
5. **Caughey–Thomas's $F$ uses $\nabla\Phi_n$, not $\nabla\psi$.**
   In the Slotboom formulation the *physically correct* field magnitude
   for velocity-saturation purposes is $|\nabla\Phi_n|$ — the gradient
   of the quasi-Fermi potential, not of the electrostatic potential.
   See [`semi/physics/mobility.py:96-116`](../../semi/physics/mobility.py)
   and the M16.1 starter prompt.

## Exercises

**Exercise 11.1.** Show that the hole current also reduces to a pure
gradient: $\mathbf{J}_p = -q\mu_p p\nabla\Phi_p$. Reproduce the algebra.

**Exercise 11.2.** At thermal equilibrium ($\Phi_n = \Phi_p = 0$),
verify that $np = n_i^2$ identically by direct substitution into (11.3).

**Exercise 11.3.** Under reverse bias $V = -1\,\mathrm{V}$ on the M3
device, what are $\Phi_n$ and $\Phi_p$ at the swept ohmic contact?
What does this imply about $n$ and $p$ at the contact?

**Exercise 11.4.** Read [`semi/physics/drift_diffusion.py:144-145`](../../semi/physics/drift_diffusion.py).
Why is `n_hat = n_from_slotboom(psi, phi_n, ni_hat)` evaluated in the
form expression rather than precomputed? What does this buy?

**Exercise 11.5.** A naive $(\psi, n, p)$ Galerkin discretization
gives $\nabla\!\cdot\!(q\mu_n n\mathbf{E})$ as a primary term, which has a
mixed-derivative structure. Show that this is *not* coercive in $H^1$,
i.e. the discrete operator can have negative or zero eigenvalues at
high Péclet number.

### Solutions

**11.1.** Same algebra: $\nabla p = (p/V_t)(\nabla\Phi_p - \nabla\psi)$.
Substitute into $\mathbf{J}_p = q\mu_p p(-\nabla\psi) - qD_p\nabla p
= -q\mu_p p\nabla\psi - q\mu_p V_t\cdot(p/V_t)(\nabla\Phi_p - \nabla\psi)
= -q\mu_p p\nabla\psi - q\mu_p p\nabla\Phi_p + q\mu_p p\nabla\psi
= -q\mu_p p\nabla\Phi_p$. ✓

**11.2.** $np = n_i\exp((\psi-\Phi_n)/V_t)\cdot n_i\exp((\Phi_p-\psi)/V_t)
= n_i^2\exp((\Phi_p-\Phi_n)/V_t)$. With $\Phi_n = \Phi_p = 0$, the
exponent vanishes, so $np = n_i^2$.

**11.3.** $\Phi_n = \Phi_p = -1\,\mathrm{V}$ at the swept (cathode)
contact. (Or +1 if it's the anode.) At a cathode (n-type) with
$\psi_R = 0.4167 - 1 = -0.5833\,\mathrm{V}$ (the equilibrium value
shifted by the applied bias): $n_R = n_i\exp((\psi_R - \Phi_n)/V_t)
= 10^{16}\exp((-0.5833 + 1)/0.02585) = 10^{16}\exp(16.12) = 10^{23}\,\mathrm{m^{-3}}$
— still equal to $N_D$, as expected for an ideal ohmic contact.

**11.4.** UFL evaluates the form expression at quadrature points, so
$n_\mathrm{hat}$ is computed *as a function of the unknown*, not a
precomputed coefficient. This means UFL can compute $\partial n_\mathrm{hat}/\partial\psi$
and $\partial n_\mathrm{hat}/\partial\Phi_n$ automatically when forming
the Jacobian. Precomputing $n$ would freeze it and break the
chain-rule pickup that UFL needs for SNES.

**11.5.** Take a 1D test problem $-(D u_x)_x + (b u)_x = f$ with
constant $D, b$ and $\mathrm{Pe} = bL/D$. The bilinear form is
$a(u, v) = \int (Du_x v_x - b u v_x)\,dx$. The first term is coercive
($\geq c\|u_x\|^2$); the second is skew-symmetric (its symmetric part
is zero) and adds nothing to coercivity. Galerkin on this gives a
discrete operator with eigenvalues $D - b\,h$, which can be negative
when $b\,h > D$, i.e. $\mathrm{Pe}\cdot h/L > 1$ — the cell Péclet
condition. SUPG and SG cure this; Slotboom changes the variable to a
problem without an advection term at all.

## Further reading

- **Slotboom, J. W. (1973).** "Computer-aided two-dimensional analysis
  of bipolar transistors." *IEEE Trans. Electron Devices* 20, 669.
  The original.
- **Selberherr, *Analysis and Simulation of Semiconductor Devices*
  (1984).** Chapter 5 §5.4. Slotboom variables in the FEM context;
  shows the coercivity proof.
- **Vasileska, Goodnick, and Klimeck (2010), Chapter 4.** Slotboom vs
  Scharfetter–Gummel comparison.
- **Brezzi, Marini, and Pietra (1989).** "Two-dimensional exponential
  fitting and applications to drift-diffusion models." *SIAM J. Numer.
  Anal.* 26, 1342. The mathematical justification for why Slotboom
  works.
- **ADR 0004** in this repo for the engine's locked decision.
