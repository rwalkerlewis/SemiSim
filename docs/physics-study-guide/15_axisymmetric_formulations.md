# 15 — Axisymmetric formulations

## Learning objectives

- Recognize when a 3D device has rotational symmetry and can be reduced
  to a 2D meridian computation.
- Derive the r-weighted Galerkin weak forms for Poisson and the
  Slotboom drift-diffusion system on the meridian half-plane $(r, z)$.
- Explain why the symmetry axis $r = 0$ is a natural (Neumann) boundary
  in the r-weighted form, and why Dirichlet contacts at $r = 0$ are
  forbidden by the schema cross-validator.
- Recognize the M14.2 axisymmetric MOSCAP benchmark and its Hu Fig. 5-18
  reproduction.

## Physical motivation

Many real devices have cylindrical symmetry: vertical p-n junction
photodiodes, MOSCAPs with circular gates, gate-all-around nanowire
transistors with circular cross-section, vertical resistors. A full 3D
mesh of such a device would be wasteful — the physical fields depend
only on $(r, z)$, not on the azimuthal angle $\theta$. Solving on the
2D **meridian half-plane** $(r, z)$ with $r \ge 0$ recovers the full 3D
solution by revolution. The price is that the volume measure becomes
$dV = 2\pi r\,dr\,dz$, and the weak forms gain an explicit factor of
$r$ as an integration weight.

The M14.2 milestone added this path. Schema 1.3.0 introduced the
top-level field `coordinate_system: "axisymmetric"` (default `"cartesian"`),
and a new module [`semi/physics/axisymmetric.py`](../../semi/physics/axisymmetric.py)
ships r-weighted variants of the equilibrium Poisson and Slotboom DD
form builders. The M14.2 benchmark
([`benchmarks/moscap_axisym_2d/`](../../benchmarks/moscap_axisym_2d/))
reproduces Hu Fig. 5-18 — a textbook MOSCAP C–V curve with circular
gate of 50 µm radius — using a 2D meridian mesh.

## Derivation from first principles

### From 3D Cartesian to 2D meridian

For a problem with $\theta$-independent fields, define cylindrical
coordinates $(r, \theta, z)$ with $r = \sqrt{x^2 + y^2}$,
$\tan\theta = y/x$, $z = z$. The volume element is

$$
dV = r\,dr\,d\theta\,dz.
\tag{15.1}
$$

Integrating an $\theta$-independent integrand over the full 3D domain
collapses the $\theta$ integral to $2\pi$:

$$
\int_\Omega f(r, z)\,dV = 2\pi\int_{\Omega_\mathrm{merid}} f(r,z)\,r\,dr\,dz.
\tag{15.2}
$$

The factor $2\pi$ cancels on both sides of any equation where every
term is a 3D volume integral, so we drop it. The remaining $r$ is the
**radial weight** in the meridian-mesh weak form.

### Gradient and divergence in cylindrical coords

For a scalar $\psi(r, z)$ (no $\theta$ dependence):

$$
\nabla\psi = \frac{\partial\psi}{\partial r}\hat{\mathbf{e}}_r
           + \frac{\partial\psi}{\partial z}\hat{\mathbf{e}}_z,
\qquad
\nabla\!\cdot\!\mathbf{F} = \frac{1}{r}\frac{\partial(rF_r)}{\partial r}
                          + \frac{\partial F_z}{\partial z}.
\tag{15.3}
$$

The Laplacian is $\nabla^2\psi = (1/r)\partial_r(r\partial_r\psi) + \partial_z^2\psi$.

### r-weighted Poisson weak form

Strong form: $-\nabla\!\cdot\!(\varepsilon\nabla\psi) = \rho$ on the
meridian half-plane (with the cylindrical $\nabla\!\cdot$ from (15.3)).
Multiply by a test function $v(r, z)$ and integrate against $r\,dr\,dz$:

$$
-\int_\Omega \nabla\!\cdot\!(\varepsilon\nabla\psi)\,v\,r\,dr\,dz
= \int_\Omega \rho v\,r\,dr\,dz.
$$

Use the cylindrical divergence identity to integrate by parts. The
result is the **r-weighted Galerkin weak form**:

$$
\boxed{
\int_\Omega \varepsilon\,\nabla\psi\!\cdot\!\nabla v\,r\,dr\,dz
= \int_\Omega \rho\,v\,r\,dr\,dz.
}
\tag{15.4}
$$

(With $\nabla$ now the meridian-plane gradient $(\partial_r, \partial_z)$.)
Compared to the Cartesian form (13.1), the *only* difference is the
factor $r$ in the integrands. UFL handles this via
`r = ufl.SpatialCoordinate(msh)[0]` and multiplying the integrands
by `r` ([`semi/physics/axisymmetric.py:49-52`](../../semi/physics/axisymmetric.py)
defines the helper `_radial_coord(msh)`).

### Boundary conditions on the meridian

The meridian half-plane has four kinds of boundary:

1. **Symmetry axis $r = 0$.** Natural BC: the r-weighted measure has
   $r \to 0$, so the test-function weight vanishes there automatically.
   No Dirichlet row is needed *or allowed*. The schema cross-validator
   [`semi/schema.py::_validate_coordinate_system`](../../semi/schema.py)
   rejects any `contacts[*].facet` that maps to a facet on the $r = 0$
   plane.
2. **Outer radial wall $r = R$.** Treated as homogeneous Neumann
   (no-flux). The user must choose $R$ large enough that the solution
   near the device is insensitive to the cutoff. For a MOSCAP,
   $R \gtrsim 5\,W_\mathrm{dmax}$ is the rule of thumb.
3. **Top and bottom planar boundaries** ($z = 0$ and $z = H$): standard
   Dirichlet (ohmic, gate) or Neumann (insulating).
4. **Interior interfaces** (e.g. Si/SiO₂): the Galerkin form encodes
   flux continuity automatically (Ch. 14).

### r-weighted DD

The same recipe applies to the continuity equations:

$$
\int_\Omega L_0^2\,\hat\mu_n\,\hat n\,\nabla\hat\Phi_n\!\cdot\!\nabla v_n\,r\,dr\,dz
- \int_\Omega \hat R\,v_n\,r\,dr\,dz = 0,
\tag{15.5}
$$

and similarly for holes with sign flip. See [`semi/physics/axisymmetric.py:142-219`](../../semi/physics/axisymmetric.py)
(`build_dd_block_residual_axisym`).

### r-weighted boundary measures

For the AC sweep / `mos_cap_ac` runner, the gate area integral becomes

$$
A_\mathrm{gate}^{3D} = 2\pi\int_\mathrm{gate facet} r\,ds_\mathrm{merid}.
\tag{15.6}
$$

The $\int r\,ds$ in the meridian computes the per-unit-radian gate
"length" that, multiplied by $2\pi$, gives the actual gate area on the
revolved 3D device. [`semi/runners/mos_cap_ac.py:188-198`](../../semi/runners/mos_cap_ac.py)
performs this assembly to convert the integrated charge $Q_\mathrm{semi}^{3D}
= 2\pi q\int\rho r\,dA$ into a per-unit-area $Q_\mathrm{gate}/A_\mathrm{gate}$.

### When axisymmetric is wrong

If the device has a non-radial feature — a finger gate, a square
contact, an off-axis implant — the axisymmetric reduction loses
information. The user must drop back to a full 2D Cartesian or 3D
solve. The schema validator does not catch this; it only enforces the
geometric invariants (`dimension == 2`, $r \geq 0$, no Dirichlet on $r=0$).
Choosing axisymmetric is a modeling decision, not just a runtime flag.

## Key results

- 3D-to-meridian volume reduction: (15.2).
- r-weighted Poisson weak form: (15.4).
- r-weighted DD block: (15.5).
- Natural BC at $r = 0$: from $r \to 0$ in the integrand.
- Outer radial wall: homogeneous Neumann; choose $R \gtrsim 5\,W_\mathrm{dmax}$.

## Worked numerical example

The M14.2 MOSCAP benchmark targets Hu Fig. 5-18 parameters:
$N_a = 5\times 10^{16}\,\mathrm{cm^{-3}}$, $T_{ox} = 10\,\mathrm{nm}$,
$\phi_{ms} = -0.95\,\mathrm{V}$. From Ch. 9, $W_\mathrm{dmax} = 144\,\mathrm{nm}$.

Choose the meridian-mesh radial extent $R = 50\,\mu\mathrm{m}$ (much
larger than $5W_\mathrm{dmax} = 720\,\mathrm{nm}$). The gate is a circle
of radius $R_\mathrm{gate} = 50\,\mu\mathrm{m}$ at the top — *exactly*
the meridian's radial extent in this benchmark, because the entire top
face is the gate. (More general benchmarks could have a smaller gate
inside a larger meridian.)

Gate area in 3D: $A_\mathrm{gate} = \pi R_\mathrm{gate}^2 = \pi(5\times 10^{-5})^2
= 7.85\times 10^{-9}\,\mathrm{m^2} = 7.85\times 10^{-5}\,\mathrm{cm^2}$.

Integrate $r$ over the gate facet in the meridian:
$L_\mathrm{gate}^\mathrm{merid} = \int_0^{R_\mathrm{gate}} r\,dr = R_\mathrm{gate}^2/2$.
Multiplied by $2\pi$: $2\pi\cdot R_\mathrm{gate}^2/2 = \pi R_\mathrm{gate}^2$, ✓.
The integral check at [`semi/runners/mos_cap_ac.py:191-192`](../../semi/runners/mos_cap_ac.py)
exercises exactly this piece.

C-V curve at $V_g = -0.7\,\mathrm{V}$ (depletion regime): the engine
solves the r-weighted Poisson + sensitivity, integrates the body charge,
and reports $C/C_{ox}$ on the curve. Comparison with Hu's Fig. 5-18
analytical line (using the LF and HF helpers in
[`semi/cv.py`](../../semi/cv.py)) shows agreement to a few percent
through accumulation, depletion, and inversion regions; see the
benchmark output and notebook 05 for the visual reproduction.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| `_radial_coord(msh)` | $r$ via `SpatialCoordinate` | `semi/physics/axisymmetric.py:49-52` |
| Axisymmetric Poisson form | (15.4) | `semi/physics/axisymmetric.py:55-103` |
| Multi-region axisym Poisson | (15.4) + cellwise $\varepsilon_r$ | `semi/physics/axisymmetric.py:106-139` |
| Axisymmetric DD block | (15.5) | `semi/physics/axisymmetric.py:142-219` |
| Schema cross-validator | $r=0$ axis-Dirichlet rejection | `semi/schema.py::_validate_coordinate_system` |
| `mos_cap_ac` axisym branch | (15.6) | `semi/runners/mos_cap_ac.py:88-198` |
| Hu Fig. 5-18 benchmark | – | `benchmarks/moscap_axisym_2d/` |
| MOSCAP analytical helpers | – | `semi/cv.py:62-360` |
| `coordinate_system` schema field | – | `schemas/input.v2.json`, `docs/schema/reference.md` |

## Existing-docs cross-reference

- [`docs/theory/axisymmetric.md`](../theory/axisymmetric.md) — definitive theory note.
- [`docs/theory/moscap_cv.md`](../theory/moscap_cv.md) — LF/HF C–V, with the depletion clamp rationale.
- [`docs/benchmarks/moscap_axisym_2d.md`](../benchmarks/moscap_axisym_2d.md) — benchmark landing page.
- [`docs/schema/reference.md` §coordinate_system](../schema/reference.md) — JSON contract.
- `notebooks/05_moscap_axisym_cv.ipynb` — end-to-end Colab demonstration.

## Common pitfalls

1. **Dirichlet on the symmetry axis.** Forbidden by the schema
   cross-validator. The natural condition is no-flux; trying to impose
   a Dirichlet there would over-constrain the system and make the
   discrete operator non-symmetric in unexpected ways.
2. **Choosing $R$ too small.** If the outer radial wall is too close
   to the device active region, the no-flux BC there shadows the
   real-device solution. The rule of thumb is $R \gtrsim 5\,W_\mathrm{dmax}$
   for MOSCAPs; for vertical pn junctions, $R \gtrsim 10\,L_n$ ensures
   diffusion currents have decayed.
3. **Forgetting the $r$ factor.** A common bug when adding a new term
   to the axisymmetric residual: write $\hat\rho\,v\,dx$ instead of
   $\hat\rho\,v\,r\,dx$. The form will assemble fine and Newton will
   converge to a *wrong* solution (the Cartesian solution scaled by
   $1/\bar r$ or some such). The MMS verification catches this; the
   M14.2 axisymmetric MMS (planned) is the way to gate against it.
4. **Boundary measures in axisymmetric.** Just like the volume measure,
   the boundary measure on a meridian-plane face must include $r$ for
   3D-equivalent integrals. [`semi/runners/mos_cap_ac.py:191-192`](../../semi/runners/mos_cap_ac.py)
   is the canonical example: `fem.assemble_scalar(fem.form(r * ds_gate))`.
5. **The $2\pi$ does *not* always cancel.** It cancels in any equation
   that's homogeneous in 3D-volume integrals. It does *not* cancel
   when comparing per-unit-area quantities (e.g. $C/A_\mathrm{gate}$)
   — there you have to multiply through by $2\pi$ explicitly to get
   the right normalization, *or* express both numerator and denominator
   in meridian-mesh-equivalent units. The mos_cap_ac runner picks the
   second route.

## Exercises

**Exercise 15.1.** A vertical p-n junction has a circular contact of
radius $R = 1\,\mu\mathrm{m}$ on top and a planar contact at the bottom.
What does the meridian-mesh look like? Where is the symmetry axis?

**Exercise 15.2.** Show that the r-weighted weak form (15.4) reduces to
the Cartesian form (13.1) far from the symmetry axis (i.e. $r$ is
roughly constant across an element). What is the validity range of
this reduction?

**Exercise 15.3.** A 50 µm-radius MOSCAP has $W_\mathrm{dmax} = 144\,\mathrm{nm}$.
Should the meridian-mesh radial extent $R$ be smaller, equal, or larger
than 50 µm? Why? Compute the minimum $R$ that respects the rule of
thumb.

**Exercise 15.4.** Read the gate-facet integration at
[`semi/runners/mos_cap_ac.py:188-192`](../../semi/runners/mos_cap_ac.py).
What does $L_\mathrm{gate}^\mathrm{merid} = \int r\,ds_\mathrm{merid}$
represent physically? Why isn't it just the radial extent of the gate
facet?

**Exercise 15.5.** The schema rejects Dirichlet on $r = 0$. Why doesn't
it also reject Neumann (e.g. an applied field) there?

### Solutions

**15.1.** The meridian half-plane is a rectangle $(r, z)$ with $0 \le r \le R_\mathrm{outer}$
and $0 \le z \le L$. The circular top contact is the segment
$\{0 \le r \le 1\,\mu\mathrm{m}, z = L\}$; the planar bottom contact is
the entire $z = 0$ edge. The symmetry axis is the left edge $r = 0$.

**15.2.** Far from $r = 0$, $r$ varies slowly over a P1 element.
Approximate $r$ by its element midpoint $\bar r$; the integrand factors
as $\bar r\,\varepsilon\nabla\psi\!\cdot\!\nabla v$ for the stiffness term,
and similarly for the source. The $\bar r$ cancels on both sides of
any *local* (per-element) test of the equation, leaving the Cartesian
form. Validity: $r$ must vary by $\ll 1$ across an element, i.e.
$h_r \ll \bar r$. Near the axis ($\bar r \to 0$), this fails and the
$r$ weight matters.

**15.3.** $R \ge 5W_\mathrm{dmax} = 720\,\mathrm{nm}$ is sufficient for
the *radial* boundary not to influence the field profile near the
gate. The gate radius is 50 µm — the outer radial wall must extend at
least to 50 µm just to contain the gate, plus some margin. The M14.2
benchmark uses $R = 50\,\mu\mathrm{m}$ exactly (the gate is the entire
top face), accepting the 1D-like geometry that this implies far from
$r = 0$. A more general benchmark with a smaller gate would extend $R$
substantially beyond the gate to allow the field to fall off radially.

**15.4.** $L_\mathrm{gate}^\mathrm{merid}$ is *half* the gate area
divided by $\pi$. That is: when you revolve the meridian around the
$z$-axis, a meridian-line segment at radius $r$ traces out an annulus
of circumference $2\pi r$ on the 3D gate; the meridian-line element
$ds$ traces out an area element $2\pi r\,ds$. So $\int r\,ds_\mathrm{merid}
= A_\mathrm{gate}^{3D}/(2\pi)$. The radial extent of the gate facet is
just $R_\mathrm{gate} - r_\mathrm{inner}$; the area picks up the $r$
factor.

**15.5.** A nonzero Neumann condition (i.e. nonzero radial flux into
the device at $r = 0$) is unphysical for a rotationally-symmetric
field — a field that points outward from the axis would have a
discontinuity at the axis itself. So the only physically meaningful
condition at $r = 0$ is *zero* flux, which is exactly what the
homogeneous Neumann case gives. The schema doesn't need to forbid
nonzero Neumann because the kronos-semi schema doesn't expose nonzero
Neumann at all (the only options are Dirichlet via `ohmic`/`gate` or
homogeneous Neumann via `insulating` / unspecified).

## Further reading

- **Sze and Ng (2007), §6.5** for cylindrical-symmetry device examples
  (gate-all-around nanowires, vertical photodetectors).
- **Hu (2010), Chapter 5 §5-9** for the MOSCAP C–V reproduction target
  the M14.2 benchmark uses.
- **`docs/theory/axisymmetric.md`** — the engine's locked theory note.
- **Cliffe, K. A. (1991).** "Numerical methods for axisymmetric
  problems." Workshop notes; classic FEM reference for cylindrical
  problems.
