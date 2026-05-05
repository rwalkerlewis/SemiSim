# 1 — Classical electromagnetism and Poisson's equation

## Learning objectives

By the end of this chapter you will be able to:

- Derive Poisson's equation $-\nabla\!\cdot\!(\varepsilon\nabla\psi) = \rho$ in matter from Maxwell's equations and the electrostatic limit.
- Explain why kronos-semi solves for the electrostatic potential $\psi$ rather than the electric field $\mathbf{E}$.
- Identify the four standard boundary conditions used in the engine — Dirichlet (ohmic, gate), homogeneous Neumann (insulating), interface flux continuity (Si/SiO₂), and the natural axis-of-symmetry BC at $r=0$ — and recognize them in the JSON schema and in `semi/bcs.py`.
- Recognize the dielectric materials shipped in `semi/materials.py` and connect their $\varepsilon_r$ to a microscopic picture of polarization.
- Locate every term of the dimensional Poisson equation in
  `semi/physics/poisson.py`.

## Physical motivation

A semiconductor device is, electrostatically, a chunk of polarizable
matter with mobile charges (electrons and holes) and fixed charges
(ionized donors and acceptors). When you bias contacts or apply a gate
voltage, charge redistributes until the resulting electric field is
self-consistent with the charge density. Poisson's equation is the
*equilibrium* statement of that self-consistency: at every interior
point, the divergence of the displacement field equals the local charge
density. Solve it on the device geometry with the right boundary
conditions and you have $\psi(\mathbf{x})$; differentiate to get
$\mathbf{E}(\mathbf{x}) = -\nabla\psi$, and from $\psi$ via Boltzmann
statistics (Ch. 3) you recover $n$ and $p$.

Why isn't the magnetic field part of the picture? Inside a typical
device the operating-frequency wavelengths ($\sim$10 cm at 1 GHz) are
many orders of magnitude longer than the device dimensions ($\sim$1 µm
to 1 mm), so $\partial\mathbf{B}/\partial t$ contributes negligibly to
$\nabla\!\times\!\mathbf{E}$. We will make this rigorous below; the
upshot is the **quasi-electrostatic** approximation, in which $\mathbf{E}$
is curl-free and admits a scalar potential.

## Derivation from first principles

Start from Maxwell's equations in matter, in SI units:

$$
\nabla\!\cdot\!\mathbf{D} = \rho_\mathrm{free},
\qquad \nabla\!\cdot\!\mathbf{B} = 0,
\tag{1.1}
$$

$$
\nabla\!\times\!\mathbf{E} = -\,\frac{\partial \mathbf{B}}{\partial t},
\qquad \nabla\!\times\!\mathbf{H} = \mathbf{J}_\mathrm{free}
                                   + \frac{\partial\mathbf{D}}{\partial t}.
\tag{1.2}
$$

The constitutive relations for a linear isotropic dielectric are
$\mathbf{D} = \varepsilon_0\varepsilon_r\mathbf{E}$ and
$\mathbf{B} = \mu_0\mu_r\mathbf{H}$. The free-charge density
$\rho_\mathrm{free}$ in a doped semiconductor is

$$
\rho_\mathrm{free} = q\,(p - n + N_D^+ - N_A^-),
\tag{1.3}
$$

with $q$ the elementary charge, $p$ and $n$ the (positive) hole and
electron densities, and $N_D^+$, $N_A^-$ the ionized donor and acceptor
densities (Ch. 4). All four are number densities in $\mathrm{m}^{-3}$.

### Quasi-electrostatic limit

Let $L$ be a characteristic device length (say 1 µm) and $T$ a
characteristic timescale (say $1/(2\pi f)$ for AC at frequency $f$).
The order-of-magnitude balance in Faraday's law is

$$
|\nabla\!\times\!\mathbf{E}| \;\sim\; \frac{|\mathbf{E}|}{L},
\qquad \left|\frac{\partial\mathbf{B}}{\partial t}\right|
   \;\sim\; \frac{|\mathbf{B}|}{T}.
$$

Inside a passive device $|\mathbf{B}|/|\mathbf{E}| \sim 1/c$ where
$c$ is the speed of light, so the magnetic-induction term is smaller
than the geometric term by a factor

$$
\frac{|\partial\mathbf{B}/\partial t|}{|\mathbf{E}|/L}
   \;\sim\; \frac{L}{cT}
   \;=\; \frac{L}{c}\cdot 2\pi f.
$$

For $L = 10\,\mu\mathrm{m}$ and $f = 1\,\mathrm{GHz}$ this ratio is
about $2\times 10^{-4}$. We drop $\partial\mathbf{B}/\partial t$ from
$\nabla\!\times\!\mathbf{E}$, which gives

$$
\nabla\!\times\!\mathbf{E} = 0 \;\;\Longrightarrow\;\;
\mathbf{E} = -\nabla\psi
\tag{1.4}
$$

for some scalar potential $\psi(\mathbf{x},t)$. The displacement-current
term in Ampère's law cannot be dropped equally cheaply — that is what
generates the small-signal AC mass matrix in Ch. 18 — but in the
electrostatic / DC steady-state regime that drives the equilibrium and
bias-sweep runners, $\partial\mathbf{D}/\partial t = 0$ as well.

### Substituting the constitutive relation

Combining $\mathbf{D} = \varepsilon_0\varepsilon_r\mathbf{E}$ with (1.4)
and the divergence equation $\nabla\!\cdot\!\mathbf{D} = \rho_\mathrm{free}$:

$$
-\nabla\!\cdot\!\bigl(\varepsilon_0\varepsilon_r(\mathbf{x})\,\nabla\psi\bigr)
= q\,(p - n + N_D^+ - N_A^-).
\tag{1.5}
$$

Equation (1.5) is the **dimensional Poisson equation** as it appears in
[`docs/PHYSICS.md` §1.1](../PHYSICS.md). Note that $\varepsilon_r$
depends on position because the device is multi-region (silicon, oxide,
etc.) — that piecewise dependence is exactly what
`build_equilibrium_poisson_form_mr` ([`semi/physics/poisson.py:82-118`](../../semi/physics/poisson.py))
handles.

### Why $\psi$ and not $\mathbf{E}$?

Two reasons. First, $\psi$ is a scalar field; $\mathbf{E}$ is a vector
field. The number of unknowns in the discrete linear system is smaller
by a factor of the spatial dimension. Second, the curl-free condition
$\nabla\!\times\!\mathbf{E}=0$ is automatic if $\mathbf{E}$ is the
gradient of a single-valued scalar; if you discretize $\mathbf{E}$ as
the primary unknown you have to enforce the curl-free constraint
separately (e.g. via Nédélec elements), which is heavier machinery than
required for elliptic problems with continuous $\psi$. Hodge theory
gives a clean justification when interfaces or topology demand
$\mathbf{H}(\mathrm{curl})$-conforming spaces; for the simply connected
domains of a typical device, $H^1$-conforming Lagrange elements on
$\psi$ are exactly the right tool (Ch. 13).

### Boundary conditions

The Poisson PDE on a bounded domain $\Omega$ has a one-parameter family
of solutions until you fix the value of $\psi$ somewhere on $\partial\Omega$.
The four BCs the engine uses are:

1. **Ohmic Dirichlet.** At an ohmic contact, $\psi$ is fixed to a
   value that satisfies local charge neutrality plus the applied bias
   (see Ch. 8 for the derivation):
   $$
   \psi_\mathrm{ohmic} \;=\;
       V_t\,\mathrm{asinh}\!\left(\frac{N_D - N_A}{2\,n_i}\right)
       + V_\mathrm{applied}.
   \tag{1.6}
   $$
   Code: [`semi/bcs.py:181-189`](../../semi/bcs.py).

2. **Gate Dirichlet.** At a gate-over-oxide contact, $\psi$ is fixed to
   the applied gate voltage minus the metal–semiconductor work-function
   difference $\phi_{ms}$:
   $$
   \psi_\mathrm{gate} \;=\; V_\mathrm{gate} - \phi_{ms}.
   \tag{1.7}
   $$
   Code: [`semi/bcs.py:186-189`](../../semi/bcs.py). The work function
   originates in the band alignment at the metal–semiconductor interface
   (Ch. 2 for $\chi$, Ch. 9 for $\phi_{ms}$).

3. **Insulating (homogeneous Neumann).** No current and no field through
   the boundary: $\nabla\psi\cdot\hat{\mathbf{n}} = 0$. This is the
   *natural* boundary condition of the Galerkin weak form — it falls out
   for free if you don't impose anything else; see Ch. 13. Code: BCs of
   type `"insulating"` are skipped during BC resolution
   ([`semi/bcs.py:34, 100`](../../semi/bcs.py)).

4. **Interface flux continuity.** Across a Si/SiO₂ interface,
   $\psi$ is continuous and the normal component of the displacement
   field $\mathbf{D}$ is continuous:
   $$
   [\![\,\psi\,]\!] = 0,
   \qquad [\![\,\varepsilon_0\varepsilon_r\,\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!] = 0.
   \tag{1.8}
   $$
   The jump bracket $[\![\,f\,]\!]$ denotes $f|_+ - f|_-$ across the
   interface. Equation (1.8) is the local form of $\nabla\!\cdot\!\mathbf{D}
   = \rho$ when there is no surface charge sheet at the interface; see
   Ch. 14 for the proof that the Galerkin form encodes (1.8) automatically.

5. **Axis-of-symmetry (forward reference).** In an axisymmetric problem
   on the meridian half-plane $(r,z)$, the symmetry axis $r=0$ is a
   *natural* boundary because the volume measure $dV = 2\pi r\,dr\,dz$
   weights test-function contributions by $r$, which vanishes there.
   Dirichlet BCs at $r=0$ are forbidden by the schema cross-validator
   ([`semi/schema.py::_validate_coordinate_system`](../../semi/schema.py));
   see Ch. 15 for the geometric picture.

## Key results

$$
\boxed{
-\nabla\!\cdot\!\bigl(\varepsilon_0\varepsilon_r(\mathbf{x})\,\nabla\psi\bigr)
   = q\,(p - n + N_D^+ - N_A^-)
}
\tag{1.9}
$$

Units: both sides have $[\mathrm{C/m^3}]$. The left side is
$[\mathrm{F/m}] \cdot [\mathrm{V/m^2}] = [\mathrm{C/m^3}]$; the right
side is $[\mathrm{C}] \cdot [\mathrm{m^{-3}}] = [\mathrm{C/m^3}]$. ✓

$$
\boxed{
\mathbf{E} = -\nabla\psi
}
\tag{1.10}
$$

Units: $[\mathrm{V/m}]$. ✓

## Worked numerical example

Take the M1 benchmark, [`benchmarks/pn_1d/pn_junction.json`](../../benchmarks/pn_1d/pn_junction.json):
silicon ($\varepsilon_r = 11.7$), $L = 2\,\mu\mathrm{m}$, doping
$10^{17}\,\mathrm{cm}^{-3}$ on each side of a step junction at the
midpoint. We will not solve the PDE here — that is Ch. 7 — but we will
unit-check the right-hand side of (1.9) at a quasi-neutral point on the
n-side.

Far inside the n-bulk, charge neutrality requires
$p - n + N_D - N_A = 0$. Convert the doping from cm⁻³ to m⁻³:
$N_D - N_A = 10^{17}\,\mathrm{cm}^{-3} = 10^{23}\,\mathrm{m}^{-3}$.

Charge neutrality is enforced by the equilibrium-Poisson initial
guess in [`semi/runners/equilibrium.py:51`](../../semi/runners/equilibrium.py):
the asinh formula (1.6) sets $\psi$ to the value at which the Boltzmann
expressions for $n$ and $p$ exactly cancel the doping. The unit check
uses $q = 1.6022\times 10^{-19}\,\mathrm{C}$:

$$
q \cdot (N_D - N_A) = 1.6022\times 10^{-19}\,\mathrm{C} \cdot 10^{23}\,\mathrm{m}^{-3}
   = 1.6022\times 10^{4}\,\mathrm{C/m^3}.
$$

That is the ionized-donor space charge that the Poisson equation has to
balance with mobile holes. The same number appears as the right-hand
side at a depletion-region point where $n,p \ll N_D - N_A$.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| Dimensional Poisson, single region | (1.9) | `semi/physics/poisson.py:25-79` (`build_equilibrium_poisson_form`) |
| Dimensional Poisson, multi-region | (1.9) with cellwise $\varepsilon_r$ | `semi/physics/poisson.py:82-118` (`build_equilibrium_poisson_form_mr`) |
| $\varepsilon$ in matter | $\mathbf{D} = \varepsilon_0\varepsilon_r\mathbf{E}$ | `semi/materials.py:38-40` (`Material.epsilon`) |
| Free-charge density | (1.3) | `semi/physics/drift_diffusion.py:156` (`rho_hat = p - n + N`) |
| Ohmic Dirichlet | (1.6) | `semi/bcs.py:181-189` (`build_psi_dirichlet_bcs`) |
| Gate Dirichlet | (1.7) | `semi/bcs.py:186-189` |
| Insulating (natural) | $\nabla\psi\cdot\hat{\mathbf{n}}=0$ | `semi/bcs.py:34` (kind skipped) |
| Interface flux continuity | (1.8) | encoded by piecewise $\varepsilon_r$ in the bilinear form (Ch. 14) |
| Axis-of-symmetry natural | r-weighted measure | `semi/physics/axisymmetric.py:24-27` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §1.1](../PHYSICS.md) — final form of (1.9) and the engine's notation.
- [`docs/PHYSICS.md` §3.1–3.3](../PHYSICS.md) — BC catalogue.
- [`docs/PHYSICS.md` §6.2](../PHYSICS.md) — interface flux continuity in the multi-region MOS context.
- [`docs/PHYSICS_INTRO.md` §2.1](../PHYSICS_INTRO.md) — Poisson in plain English.
- [`docs/theory/axisymmetric.md`](../theory/axisymmetric.md) — natural BC at $r=0$.

## Common pitfalls

1. **Densities in cm⁻³ vs. m⁻³.** The JSON schema takes doping and
   intrinsic densities in cm⁻³ because that is the device-physics
   convention; everything inside the engine is SI ($\mathrm{m^{-3}}$).
   The conversion happens once at ingest in
   [`semi/constants.py:29-31`](../../semi/constants.py)
   (`cm3_to_m3`). If you insert a number directly into one of the UFL
   forms without the helper you will be off by $10^6$.
2. **Sign of $\rho$ in the displayed equation.** Some textbooks write
   $\nabla\!\cdot\!(\varepsilon\nabla\psi) = -\rho$ (with the minus on the
   right), others write $-\nabla\!\cdot\!(\varepsilon\nabla\psi) = \rho$
   (with the minus on the left). Both are correct; they differ by a
   notational choice for which side absorbs the negation. kronos-semi
   uses the second form, matching `build_equilibrium_poisson_form`'s
   residual `+stiffness − rho_hat`.
3. **$\mathbf{E}$ vs $\mathbf{D}$ continuity.** At a dielectric interface,
   $\mathbf{D}\cdot\hat{\mathbf{n}}$ is continuous; $\mathbf{E}\cdot\hat{\mathbf{n}}$
   is *not* (it jumps by the ratio of permittivities). $\psi$ is
   continuous; $\nabla\psi$ is not. A common beginner mistake is to
   demand $\mathbf{E}$-continuity at Si/SiO₂ — the Galerkin form
   correctly does not.
4. **Why insulators carry only $\varepsilon_r$.** In the SiO₂ region of
   a MOS capacitor there are no mobile carriers and no doping, so
   $\rho = 0$ inside the oxide. The Poisson equation reduces to
   Laplace's equation $-\nabla\!\cdot\!(\varepsilon_r\nabla\psi)=0$,
   and the only material datum the oxide needs to provide is its
   relative permittivity. This is why `Material` instances with
   `role="insulator"` set every other field to zero
   ([`semi/materials.py:86-88`](../../semi/materials.py)).

## Exercises

**Exercise 1.1.** Show that for a uniform doping $N_D - N_A = N$ with
zero applied bias, the Poisson equation has a constant solution
$\psi(\mathbf{x}) = \psi_\mathrm{eq}$ with $n(\psi_\mathrm{eq}) -
p(\psi_\mathrm{eq}) = N$. (Hint: bulk charge neutrality.)

**Exercise 1.2.** Verify that for a thin capacitor of thickness $d$,
plate area $A$, dielectric $\varepsilon_r$, and applied voltage $V$,
the integral $-\int_\Omega \nabla\psi\cdot\nabla\psi\,dV / V^2$ equals
$\varepsilon_r\varepsilon_0 A / d$. This is the parallel-plate
capacitance, which we will revisit when computing $C_{ox}$ in Ch. 9.

**Exercise 1.3.** Take the silicon entry from
[`semi/materials.py`](../../semi/materials.py): $\varepsilon_r = 11.7$.
Compute the absolute permittivity $\varepsilon = \varepsilon_0\varepsilon_r$
in F/m and check it matches `Material(...).epsilon`. Repeat for SiO₂
($\varepsilon_r = 3.9$) and HfO₂ ($\varepsilon_r = 25.0$).

**Exercise 1.4.** Show that the homogeneous Neumann condition
$\nabla\psi\cdot\hat{\mathbf{n}}=0$ on a portion of $\partial\Omega$
falls out of the Galerkin form
$\int\varepsilon\nabla\psi\cdot\nabla v\,dV = \int\rho v\,dV$
for any test function $v$ that is allowed to be nonzero on that
portion. (Hint: integrate by parts and ask which boundary integral
must vanish for the equation to hold.)

**Exercise 1.5.** At a sharp Si/SiO₂ interface with no surface charge,
prove that $\nabla\psi$ has a jump $[\![\,\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!]
= -[\![\,\varepsilon_r\,]\!]/\varepsilon_r^+\cdot \nabla\psi^+\!\cdot\!\hat{\mathbf{n}}$
even though $\psi$ is continuous. Estimate the magnitude of this jump
for the M6 MOS capacitor at $V_g = 1\,\mathrm{V}$.

### Solutions

**1.1.** With $V_\mathrm{applied} = 0$ and $N$ uniform, the asinh
formula (1.6) gives a single value $\psi_\mathrm{eq} = V_t\,\mathrm{asinh}(N/2n_i)$
that is independent of $\mathbf{x}$. Substitution into the Boltzmann
expressions (Ch. 3) gives $n - p = 2n_i\sinh(\psi_\mathrm{eq}/V_t) = N$
exactly. So the constant $\psi_\mathrm{eq}$ satisfies (1.9) trivially
(both sides are zero) and the homogeneous Neumann BC on the side walls.

**1.2.** $\psi$ varies linearly through the dielectric:
$\psi(z) = V z/d$ for $z \in [0, d]$. Then $\nabla\psi = V/d\,\hat{\mathbf{z}}$,
$|\nabla\psi|^2 = V^2/d^2$. The volume integral evaluates to
$\varepsilon_r\varepsilon_0 \cdot V^2/d^2 \cdot Ad = \varepsilon_r\varepsilon_0 A V^2 / d$;
dividing by $V^2$ leaves $\varepsilon_r\varepsilon_0 A/d = C$. ✓

**1.3.** Si: $11.7 \times 8.854\times 10^{-12} = 1.036\times 10^{-10}\,\mathrm{F/m}$.
SiO₂: $3.9 \times 8.854\times 10^{-12} = 3.453\times 10^{-11}\,\mathrm{F/m}$.
HfO₂: $25.0 \times 8.854\times 10^{-12} = 2.214\times 10^{-10}\,\mathrm{F/m}$.

**1.4.** Integrating $-\nabla\!\cdot\!(\varepsilon\nabla\psi)v$ by parts:
$\int_\Omega \varepsilon\nabla\psi\!\cdot\!\nabla v\,dV
- \int_{\partial\Omega} v\,\varepsilon\,\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,dS
= \int_\Omega \rho v\,dV$. The Galerkin form drops the surface term,
so unless we replace it with an explicit `ds` integral, we are
implicitly demanding $\varepsilon\nabla\psi\!\cdot\!\hat{\mathbf{n}}=0$
on the part of $\partial\Omega$ where $v$ is free. ✓

**1.5.** Continuity of $\psi$ across the interface is imposed
by the conforming function space; continuity of
$\varepsilon\nabla\psi\cdot\hat{\mathbf{n}}$ is encoded by the
bilinear form $\int\varepsilon_r\nabla\psi\cdot\nabla v\,dx$ with the
piecewise constant $\varepsilon_r(\mathbf{x})$. Solving for the field
jump: $\varepsilon_r^+(\nabla\psi^+\!\cdot\!\hat{\mathbf{n}}) =
\varepsilon_r^-(\nabla\psi^-\!\cdot\!\hat{\mathbf{n}})$, so
$\nabla\psi^-\!\cdot\!\hat{\mathbf{n}} - \nabla\psi^+\!\cdot\!\hat{\mathbf{n}}
= (1 - \varepsilon_r^+/\varepsilon_r^-)\,\nabla\psi^+\!\cdot\!\hat{\mathbf{n}}$.
At Si/SiO₂ with $\varepsilon_r^+ = 3.9$ (oxide), $\varepsilon_r^- = 11.7$
(silicon), the field is $11.7/3.9 = 3$ times larger in the oxide than in
the silicon at the interface; this is the quantitative origin of the
"oxide field" in MOS physics (Ch. 9).

## Further reading

- **Griffiths, *Introduction to Electrodynamics*, 4th ed. (2013).**
  Chapters 2 and 4. Maxwell's equations in matter, dielectrics, BCs at
  dielectric interfaces. Standard undergraduate text; this is exactly
  the level you should be at coming into Ch. 1.
- **Jackson, *Classical Electrodynamics*, 3rd ed. (1998).** §1.5 for
  Poisson, §1.9 for boundary conditions at a dielectric interface, §6.4
  for the quasi-electrostatic limit.
- **Sze and Ng, *Physics of Semiconductor Devices*, 3rd ed. (2007).**
  Appendix C tabulates dielectric constants for Si, Ge, GaAs, SiO₂, and
  Si₃N₄ (the same numbers appear in [`semi/materials.py`](../../semi/materials.py)).
- **Hu, *Modern Semiconductor Devices for Integrated Circuits* (2010).**
  Chapter 2 for Maxwell-to-Poisson in a device-physics context.
