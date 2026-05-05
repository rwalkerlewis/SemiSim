# 5 — Transport: drift-diffusion

## Learning objectives

- Sketch the Boltzmann transport equation and identify the moments that
  produce the drift-diffusion (DD) currents.
- State the Einstein relation $D = \mu V_t$ and derive it from detailed
  balance at thermal equilibrium.
- Write the steady-state continuity equations
  $\nabla\!\cdot\!\mathbf{J}_n = qR$, $\nabla\!\cdot\!\mathbf{J}_p = -qR$
  and explain the sign convention.
- Recognize why naive Galerkin discretization of $\mathbf{J} = \mu n\mathbf{E}
  + D\nabla n$ fails when drift dominates diffusion (the Péclet
  problem) and why Slotboom (Ch. 11) cures it.
- Locate the constant- and Caughey–Thomas-mobility branches in
  `semi/physics/mobility.py` and explain when each is appropriate.

## Physical motivation

Carriers move because of two distinct physical drivers. **Drift** is
the response to an applied electric field: an electron in a field
$\mathbf{E}$ accelerates between scattering events and reaches an
average velocity $\mathbf{v}_n = -\mu_n\mathbf{E}$ (the minus sign
because electrons are negatively charged). **Diffusion** is the
response to a concentration gradient: carriers spread from regions of
high concentration to low without needing a field, by random thermal
motion. Both processes happen simultaneously, and at a pn junction in
equilibrium they cancel each other exactly (Ch. 7).

Capturing both in one equation requires writing the *current density*
$\mathbf{J}$ as a sum of a drift term proportional to $\mathbf{E}$ and
a diffusion term proportional to $\nabla n$. The continuity equations
then conserve carriers locally, with sources and sinks coming from
generation–recombination (Ch. 6). Together — Poisson, two
continuities, and constitutive expressions for $\mathbf{J}$ —
the system has three coupled PDEs in three unknowns $(\psi, n, p)$,
which the engine reformulates in Slotboom variables before solving
(Ch. 11).

## Derivation from first principles

### Boltzmann transport equation

The fundamental kinetic description of an electron gas is the
Boltzmann transport equation (BTE) for the distribution function
$f(\mathbf{r}, \mathbf{k}, t)$:

$$
\frac{\partial f}{\partial t}
+ \mathbf{v}_\mathbf{k}\!\cdot\!\nabla_\mathbf{r} f
+ \frac{q}{\hbar}\,\mathbf{E}\!\cdot\!\nabla_\mathbf{k} f
= \left.\frac{\partial f}{\partial t}\right|_\mathrm{coll}.
\tag{5.1}
$$

The left side is convection in real space and acceleration in
$\mathbf{k}$-space; the right is the collision operator capturing
phonon, impurity, and carrier–carrier scattering. The drift-diffusion
limit comes from taking moments of (5.1) under the relaxation-time
approximation $(\partial f/\partial t)_\mathrm{coll} \approx -(f - f_0)/\tau_m$
and assuming the distribution stays close to local equilibrium $f_0$.

### Reduction to drift-diffusion

The zeroth $\mathbf{k}$-moment of (5.1) gives the continuity equation:

$$
\frac{\partial n}{\partial t} + \nabla\!\cdot\!(n\mathbf{v}_n) = G - R,
\tag{5.2}
$$

with $G$, $R$ generation and recombination from the collision operator.
The first moment, after the relaxation-time and gradient-expansion
approximations, gives the **drift-diffusion current**

$$
\mathbf{J}_n = q\,\mu_n n\,\mathbf{E} + q\,D_n\,\nabla n,
\qquad
\mathbf{J}_p = q\,\mu_p p\,\mathbf{E} - q\,D_p\,\nabla p,
\tag{5.3}
$$

where $\mu_n = q\tau_m/m_n^*$ is the electron mobility (sign absorbed
into the velocity convention $\mathbf{v}_n = -\mu_n\mathbf{E}$) and
$D_n$ is the diffusivity. The signs differ between (5.3a) and (5.3b)
because electrons drift opposite to $\mathbf{E}$ and holes drift along
it.

The full derivation from (5.1) to (5.3) is in Selberherr §3 or in
Vasileska, Goodnick, and Klimeck, *Computational Electronics*, Chapter 1.
We accept (5.3) as the working starting point.

### Einstein relation

At thermal equilibrium with no applied bias, the net electron current
is zero everywhere: drift and diffusion cancel locally. Setting
$\mathbf{J}_n = 0$ in (5.3a) and using the equilibrium Boltzmann form
$n(\psi) = n_i\,e^{\psi/V_t}$ from (3.7) with $\Phi_n = 0$:

$$
0 = q\,\mu_n n\,(-\nabla\psi) + q D_n\,\nabla n
   = q n\,(-\mu_n\nabla\psi + D_n\,\nabla\psi/V_t),
$$

so $D_n = \mu_n V_t$. This is the **Einstein relation**:

$$
\boxed{
D_n = \mu_n V_t,
\qquad
D_p = \mu_p V_t,
\qquad
V_t = kT/q.
}
\tag{5.4}
$$

Substituting (5.4) into (5.3) and using $\mathbf{E} = -\nabla\psi$:

$$
\mathbf{J}_n = -q\mu_n n\,\nabla\psi + q\mu_n V_t\,\nabla n
   = q\mu_n V_t\,\bigl(\nabla n - (n/V_t)\nabla\psi\bigr).
\tag{5.5}
$$

This is the form that goes into the engine's terminal-current evaluator
([`semi/postprocess.py:97-99`](../../semi/postprocess.py)). For the
*solver* itself, the Slotboom transformation (Ch. 11) rewrites (5.5)
as a pure gradient of $\Phi_n$, which is what makes Galerkin
discretization stable.

### Continuity equations and sign convention

In steady state, (5.2) drops the time derivative, and combining with the
hole equation:

$$
\boxed{
\nabla\!\cdot\!\mathbf{J}_n = +q\,R(n,p),
\qquad
\nabla\!\cdot\!\mathbf{J}_p = -q\,R(n,p),
}
\tag{5.6}
$$

with $R$ the *net* recombination (positive when electrons and holes
disappear in pairs). The signs differ between (5.6a) and (5.6b) because
recombination removes one electron and one hole simultaneously, but
their respective charges have opposite sign — divergence of $\mathbf{J}_n$
matches the rate of electron *removal* with a $+q$ sign, while
divergence of $\mathbf{J}_p$ matches the rate of hole *removal* with a
$-q$ sign.

Adding (5.6a) and (5.6b) gives $\nabla\!\cdot\!(\mathbf{J}_n + \mathbf{J}_p) = 0$,
which is the **total-current conservation** law that the V&V suite
checks per-target-bias to within 5% forward / 15% reverse
([`docs/PHYSICS.md` §5.3](../PHYSICS.md)).

### Why naive Galerkin fails: Péclet number

If you discretize (5.5) with standard Galerkin Lagrange elements on
$(n, p, \psi)$ as primary unknowns, the convective term
$-q\mu_n n\,\nabla\psi$ dominates the diffusive $q\mu_n V_t\,\nabla n$ wherever
$|\nabla\psi|$ is large compared with $V_t/h$, where $h$ is the mesh
spacing. The dimensionless ratio is the cell **Péclet number**

$$
\mathrm{Pe} = \frac{|\nabla\psi|\,h}{V_t} = \frac{|E|\,h}{V_t}.
\tag{5.7}
$$

A 1 V drop over a 1 µm cell at $V_t = 25\,\mathrm{mV}$ gives $\mathrm{Pe} \approx 40$.
At Pe $\gtrsim 1$ the discrete operator stops being coercive and
solutions develop nonphysical oscillations ("wiggles"). Stabilization
schemes like SUPG and Scharfetter–Gummel exist to cure this; the
engine instead changes primary unknowns (Slotboom, Ch. 11) so the
underlying PDE has no Péclet problem at all. See ADR 0004 for the
rationale.

### Mobility

Constant mobility is the simplest model: $\mu_n$ and $\mu_p$ are
material parameters, independent of field. The shipped engine defaults
to $\mu_n = 1400\,\mathrm{cm^2 V^{-1} s^{-1}}$ and
$\mu_p = 450\,\mathrm{cm^2 V^{-1} s^{-1}}$ for silicon at 300 K
([`semi/materials.py:59-60`](../../semi/materials.py)). These are the
*low-field bulk* values; mobility falls at high field (carriers
saturate at $v_\mathrm{sat} \approx 10^7\,\mathrm{cm/s}$ in silicon)
and at high doping (impurity scattering).

**Caughey–Thomas (M16.1, shipped).** Velocity saturation at high field:

$$
\mu(F) = \frac{\mu_0}{\bigl(1 + (\mu_0 F / v_\mathrm{sat})^\beta\bigr)^{1/\beta}},
\qquad F = |\nabla\Phi_n|\,(\text{or }|\nabla\Phi_p|).
\tag{5.8}
$$

Limits: $F\to 0$ gives $\mu \to \mu_0$ (low-field); $F\to\infty$ gives
$\mu F \to v_\mathrm{sat}$ (saturation). Default $\beta_n = 2$,
$\beta_p = 1$. Implementation:
[`semi/physics/mobility.py:57-84, 119-217`](../../semi/physics/mobility.py).
The Caughey–Thomas branch keys on `physics.mobility.model: "caughey_thomas"`
and is verified by MMS-DD Variant D at L² ≥ 1.99 on every block.

**Lombardi (M16.2, planned).** Surface-scattering composite needed for
inversion-layer mobility in MOSFETs:

$$
\frac{1}{\mu} = \frac{1}{\mu_\mathrm{Coulomb}}
              + \frac{1}{\mu_\mathrm{phonon}}
              + \frac{1}{\mu_\mathrm{surface}}.
\tag{5.9}
$$

Forward reference to [`docs/IMPROVEMENT_GUIDE.md` §M16.2](../IMPROVEMENT_GUIDE.md).
Will tighten the `mosfet_2d` Pao–Sah verifier from 20% to 10% once
shipped.

## Key results

- Drift-diffusion currents: (5.3).
- Einstein relation: (5.4).
- Steady-state continuity with sign convention: (5.6).
- Péclet failure mode: (5.7).
- Caughey–Thomas mobility: (5.8).

## Worked numerical example

For the M2 forward-bias pn junction at $V = 0.6\,\mathrm{V}$, the
saturation current is given by the Shockley long-diode formula:

$$
J_s = q n_i^2\left(\frac{D_n}{L_n N_A} + \frac{D_p}{L_p N_D}\right),
\qquad L_n = \sqrt{D_n\tau_n},\quad L_p = \sqrt{D_p\tau_p}.
$$

Plug in (using SI throughout):
- $n_i = 10^{16}\,\mathrm{m^{-3}}$, so $n_i^2 = 10^{32}\,\mathrm{m^{-6}}$.
- $\mu_n = 1400\,\mathrm{cm^2/Vs} = 0.14\,\mathrm{m^2/Vs}$,
  $\mu_p = 450\,\mathrm{cm^2/Vs} = 0.045\,\mathrm{m^2/Vs}$.
- Einstein: $D_n = 0.14 \cdot 0.02585 = 3.62\times 10^{-3}\,\mathrm{m^2/s}$,
  $D_p = 0.045 \cdot 0.02585 = 1.16\times 10^{-3}\,\mathrm{m^2/s}$.
- Lifetimes $\tau_n = \tau_p = 10^{-7}\,\mathrm{s}$ (default).
  $L_n = \sqrt{3.62\times 10^{-3} \cdot 10^{-7}} = 1.90\times 10^{-5}\,\mathrm{m} = 19\,\mu\mathrm{m}$.
  $L_p = \sqrt{1.16\times 10^{-3} \cdot 10^{-7}} = 1.08\times 10^{-5}\,\mathrm{m} = 10.8\,\mu\mathrm{m}$.
- $N_A = N_D = 10^{23}\,\mathrm{m^{-3}}$.
- $J_s = 1.602\times 10^{-19} \cdot 10^{32}\cdot (3.62\times 10^{-3}/(1.90\times 10^{-5}\cdot 10^{23})
   + 1.16\times 10^{-3}/(1.08\times 10^{-5}\cdot 10^{23}))$.
- First term: $3.62\times 10^{-3}/1.90\times 10^{18} = 1.91\times 10^{-21}$.
- Second: $1.16\times 10^{-3}/1.08\times 10^{18} = 1.07\times 10^{-21}$.
- Sum: $2.98\times 10^{-21}$.
- $J_s = 1.602\times 10^{-19} \cdot 10^{32} \cdot 2.98\times 10^{-21}
   = 4.78\times 10^{-8}\,\mathrm{A/m^2}$.

At $V = 0.6\,\mathrm{V}$, $\exp(V/V_t) = e^{23.21} = 1.20\times 10^{10}$.
$J_\mathrm{Shockley} = J_s\cdot(e^{V/V_t} - 1) \approx J_s\cdot 1.20\times 10^{10}
= 5.74\times 10^2\,\mathrm{A/m^2}$.

The actual benchmark reports $\sim 1.6\times 10^3\,\mathrm{A/m^2}$ at
$V = 0.6\,\mathrm{V}$ ([`docs/PHYSICS.md` §2.5 close](../PHYSICS.md)),
about 2.8× higher. The factor-of-3 discrepancy is real and is the
standard "ideality factor n=1.0 vs effective n>1" mismatch: the M2
device is short enough ($L=2\,\mu\mathrm{m}$ vs $L_n=19\,\mu\mathrm{m}$,
$L_p=10.8\,\mu\mathrm{m}$) that the **short-base** approximation
applies and the long-diode $L_{n,p}$ in the saturation current should be
replaced by the metallurgical region width. Within an order of magnitude,
agreement is the right call. The benchmark verifier accepts 10%
agreement on $V \geq 0.5\,\mathrm{V}$ specifically because of effects
like this.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| DD current $\mathbf{J}_n$ | (5.3) | `semi/postprocess.py:97-99` (`evaluate_current_at_contact`) |
| Einstein $D = \mu V_t$ | (5.4) | implicit; `D = sc.V0 * sc.mu0` in `Scaling.D0` (`semi/scaling.py:48-50`) |
| Continuity in Slotboom | (5.6) | `semi/physics/drift_diffusion.py:158-169` |
| Constant mobility branch | $\mu = \mu_0$ | `semi/physics/mobility.py:87-93` (`constant_mu`) |
| Caughey–Thomas mobility | (5.8) | `semi/physics/mobility.py:57-84` (`caughey_thomas_mu`) |
| Mobility dispatcher | – | `semi/physics/mobility.py:119-217` (`build_mobility_expressions`) |
| Saturation velocity (scaled) | – | `semi/physics/mobility.py:96-116` (`caughey_thomas_vsat_for_form`) |
| Shockley saturation current | – | `semi/diode_analytical.py:19-37` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §1.3](../PHYSICS.md) — Slotboom-form continuity equations.
- [`docs/PHYSICS.md` §2.5](../PHYSICS.md) — scaled DD residual.
- [`docs/PHYSICS_INTRO.md` §2.3, §3.1](../PHYSICS_INTRO.md) — narrative version.
- [`docs/theory/slotboom.md`](../theory/slotboom.md) — why Slotboom over SG.
- [`docs/adr/0004-slotboom-variables-for-dd.md`](../adr/0004-slotboom-variables-for-dd.md).
- [`docs/M16_1_STARTER_PROMPT.md`](../M16_1_STARTER_PROMPT.md) — Caughey–Thomas details.

## Common pitfalls

1. **Drift sign for electrons.** In the convention where
   $\mathbf{J}_n = q\mu_n n\mathbf{E} + qD_n\nabla n$ (positive q,
   positive μ), the *current* flows along $\mathbf{E}$, but electron
   *velocity* is $\mathbf{v}_n = -\mu_n\mathbf{E}$. The negative-charge
   sign has been absorbed into the current expression.
2. **Mobility units.** JSON inputs are in $\mathrm{cm^2 V^{-1} s^{-1}}$;
   conversion via `cm2_to_m2` ([`semi/constants.py:39-41`](../../semi/constants.py))
   produces $\mathrm{m^2 V^{-1} s^{-1}}$. A factor $10^4$ separates the
   two; using the wrong units gives currents $10^4$ too small or large.
3. **Péclet wiggles look like physics.** A drift-dominated DD solve on
   a coarse mesh with primary $(n, \psi)$ produces a smoothly oscillating
   $n(x)$ with hundreds of decades amplitude. It looks like a converged
   solution by SNES residual norms, because the discrete equations are
   what they are. The fix is not to "tune the solver"; it is to change
   primary variables (Slotboom) or stabilize (SG / SUPG).
4. **Constant mobility is wrong at high field.** A 0.5 V drop over a
   100 nm channel gives $5\times 10^6\,\mathrm{V/m}$, well into
   saturation: real silicon has $\mu(F) \approx \mu_0/3$ at this field.
   The constant-mobility benchmarks accept 20% agreement specifically
   because of this. M16.1 Caughey–Thomas tightens the gap at high field.
5. **Mobility models couple to Slotboom.** In the Caughey–Thomas branch
   of `build_mobility_expressions`, $F$ is taken as
   $|\nabla\Phi_n|$ or $|\nabla\Phi_p|$ — the gradient of the
   *quasi-Fermi potential*, not of $\psi$. This is consistent with the
   Slotboom flux form (Ch. 11) and matches the dimensional-analysis note
   in [`semi/physics/mobility.py:96-116`](../../semi/physics/mobility.py).

## Exercises

**Exercise 5.1.** Derive the Einstein relation for holes by setting
$\mathbf{J}_p = 0$ at equilibrium with $p = n_i\exp(-\psi/V_t)$.

**Exercise 5.2.** For the M1 device (1 µm half, 10¹⁷ doping), estimate
the cell Péclet number on a 100-cell mesh at $V_{bi} = 0.83\,\mathrm{V}$.
Conclude that naive Galerkin would produce wiggles.

**Exercise 5.3.** Compute the diffusion length $L_n = \sqrt{D_n\tau_n}$
in silicon for $\tau_n = 10^{-7}\,\mathrm{s}$. Repeat for $\tau_n =
10^{-9}\,\mathrm{s}$ (heavily damaged silicon). Comment on which case
is "short-base" for a 2 µm device.

**Exercise 5.4.** At what field does Caughey–Thomas (5.8) predict
$\mu = \mu_0/2$ for electrons in silicon, with $\beta_n = 2$ and
$v_\mathrm{sat} = 10^7\,\mathrm{cm/s}$?

**Exercise 5.5.** Show that adding the two continuities (5.6) gives
$\nabla\!\cdot\!(\mathbf{J}_n + \mathbf{J}_p) = 0$ — total current is
divergence-free, i.e. has the same value on every cross-section of a
1D device.

### Solutions

**5.1.** $0 = q\mu_p p\,\mathbf{E} - qD_p\,\nabla p
= q\mu_p p(-\nabla\psi) - qD_p\,(-(p/V_t)\nabla\psi)
= q p \nabla\psi(-\mu_p + D_p/V_t)$. So $D_p = \mu_p V_t$, same as electrons.

**5.2.** Mesh spacing $h = 2\,\mu\mathrm{m}/100 = 20\,\mathrm{nm}$.
Field at the junction is $\sim V_{bi}/W \approx 0.83/146\,\mathrm{nm}
\approx 5.7\times 10^6\,\mathrm{V/m}$. So $\mathrm{Pe} = Eh/V_t
= 5.7\times 10^6\cdot 20\times 10^{-9}/0.02585 = 4.4$. Well above 1;
naive Galerkin would wiggle. Slotboom sidesteps it.

**5.3.** $D_n = 3.62\times 10^{-3}\,\mathrm{m^2/s}$. $L_n =
\sqrt{D_n\cdot 10^{-7}} = 1.90\times 10^{-5}\,\mathrm{m} = 19\,\mu\mathrm{m}$
for $\tau = 10^{-7}$. $L_n = \sqrt{D_n\cdot 10^{-9}} = 1.90\times 10^{-6}\,\mathrm{m}
= 1.9\,\mu\mathrm{m}$ for $\tau = 10^{-9}$. The 2 µm device is
"short-base" relative to the first ($L\approx W$); "long-base" relative
to the second ($L \gg W$ doesn't hold here either; both are
intermediate). Short-base reduces $J_s$ by a factor $W/L_n$; this is the
correction missing from the naive long-diode formula.

**5.4.** Setting $\mu(F)/\mu_0 = 1/2$: $1/(1 + (\mu_0 F/v_s)^2)^{1/2} = 1/2$,
so $1 + (\mu_0 F/v_s)^2 = 4$, $\mu_0 F/v_s = \sqrt 3$, $F = \sqrt 3\,v_s/\mu_0$.
With $v_s = 10^5\,\mathrm{m/s}$ and $\mu_0 = 0.14\,\mathrm{m^2/Vs}$:
$F = 1.732\cdot 10^5/0.14 = 1.24\times 10^6\,\mathrm{V/m} = 12.4\,\mathrm{kV/cm}$.

**5.5.** Add (5.6a) and (5.6b): $\nabla\!\cdot\!(\mathbf{J}_n + \mathbf{J}_p)
= +qR + (-qR) = 0$. The recombination is equal-and-opposite on the two
rows, so it cancels. In 1D, this means $J_n(x) + J_p(x)$ is constant in
$x$ across the device — the **total-current conservation** that the
V&V conservation gate verifies on ten interior facets per bias point.

## Further reading

- **Selberherr, *Analysis and Simulation of Semiconductor Devices*
  (1984).** Chapter 3 derives the DD model from the BTE moments. Best
  reference for *why* DD is the limit it is.
- **Vasileska, Goodnick, and Klimeck, *Computational Electronics*
  (2010).** Chapter 1 covers BTE → DD with worked examples; chapter 5
  covers Caughey–Thomas and Lombardi mobility models.
- **Sze and Ng (2007), §1.6** for mobility data on Si, Ge, GaAs.
- **Caughey, D. M., and Thomas, R. E. (1967).** "Carrier mobilities in
  silicon empirically related to doping and field." *Proc. IEEE* 55,
  2192. The original paper for (5.8).
- **Lombardi, C., et al. (1988).** "A physically based mobility model for
  numerical simulation of nonplanar devices." *IEEE Trans. CAD* 7, 1164.
  Source for (5.9), planned in M16.2.
- **Scharfetter, D. L., and Gummel, H. K. (1969).** "Large-signal
  analysis of a silicon Read diode oscillator." *IEEE Trans. Electron
  Devices* 16, 64. The discretization that solved the Péclet problem
  before Slotboom.
