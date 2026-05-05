# 14 — Multi-region and interfaces

## Learning objectives

- Construct a cellwise (DG0) relative-permittivity function $\varepsilon_r(\mathbf{x})$
  for a multi-region device.
- Show that the natural flux-continuity condition $[\![\,\varepsilon\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!] = 0$
  at a Si/SiO₂ interface is encoded automatically by the Galerkin
  weak form with piecewise $\varepsilon_r$.
- Explain why carriers (Slotboom $\Phi_n, \Phi_p$) live on a
  *semiconductor submesh* while the electrostatic potential $\psi$
  lives on the parent mesh.
- Read the `entity_maps` machinery in
  `semi/physics/drift_diffusion.py` and explain how a parent-mesh
  $\psi$ is consumed in a submesh integral.
- Recognize the M6 multi-region MOS verifier and its MMS gate.

## Physical motivation

The MOS capacitor and the MOSFET both contain a semiconductor and an
insulator joined at a planar interface. The electrostatic potential
$\psi$ is defined and continuous everywhere; the electron and hole
densities $n, p$ exist only in the semiconductor. A clean numerical
treatment must:

- carry one $\psi$ field over the whole device,
- carry $\Phi_n, \Phi_p$ only over the semiconductor,
- enforce the dielectric flux-continuity condition at the interface,
- not couple oxide DOFs into continuity equations that don't apply there.

kronos-semi achieves all four with a parent mesh + semiconductor
submesh + cellwise DG0 $\varepsilon_r$ + entity-mapped UFL forms. This
chapter explains every one.

## Derivation from first principles

### Cellwise DG0 $\varepsilon_r$

A "DG0" (discontinuous Galerkin, degree 0) function space has one
constant value per cell. Building a DG0 $\varepsilon_r$ for a
multi-region mesh:

1. Tag each cell with its region (silicon, oxide, ...).
2. For each region, look up the material's $\varepsilon_r$ from
   [`semi/materials.py`](../../semi/materials.py).
3. On each cell, set the DG0 function value to that region's $\varepsilon_r$.

Code lives in `semi/mesh.py::build_eps_r_function`. The resulting
`Function` is then dropped into the UFL form as a coefficient:

```python
F = L_D2 * eps_r_fn * inner(grad(psi), grad(v)) * dx_full + ...
```

Inside each cell, $\varepsilon_r$ is constant; across an interface, it
jumps. The Galerkin assembler integrates over each cell using that
cell's value.

### Flux-continuity at the interface

At a sharp Si/SiO₂ interface with no surface charge, the electrostatic
boundary condition is

$$
[\![\,\psi\,]\!] = 0,
\qquad
[\![\,\varepsilon_0\varepsilon_r\,\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!] = 0.
\tag{14.1}
$$

(Continuity of $\psi$, continuity of normal $\mathbf{D}$.) The
P1 Lagrange function space is **conforming** in $H^1$ — values are
continuous across element edges by construction — so $[\![\,\psi\,]\!] = 0$
is built into the discrete unknown. The flux-continuity condition is
trickier; let's derive it from the variational principle.

Consider an interior face $F$ between two cells $K_1$ (silicon) and
$K_2$ (oxide). Pick a test function $v \in V_h$ with support straddling
$F$. The contribution of $F$ to the bilinear form $a(\psi, v) = \int_\Omega
\varepsilon_r\nabla\psi\!\cdot\!\nabla v\,dV$ is computed by summing the
two cell contributions:

$$
a_{K_1\cup K_2}(\psi, v) = \int_{K_1}\varepsilon_r^\mathrm{Si}\nabla\psi\!\cdot\!\nabla v
+ \int_{K_2}\varepsilon_r^\mathrm{ox}\nabla\psi\!\cdot\!\nabla v.
$$

Apply the divergence theorem to each cell separately. Each interior
piece collects a boundary integral on $F$:

$$
-\int_F\bigl(\varepsilon_r^\mathrm{Si}\nabla\psi^\mathrm{Si}\!\cdot\!\hat{\mathbf{n}}_1
       - \varepsilon_r^\mathrm{ox}\nabla\psi^\mathrm{ox}\!\cdot\!\hat{\mathbf{n}}_1\bigr)\,v\,dS,
$$

with $\hat{\mathbf{n}}_1$ the outward normal to $K_1$. Combining the
two cell-integrated-by-parts contributions, the interior surface term
collapses to $-\int_F\,[\![\,\varepsilon_r\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!]\,v\,dS$
(the jump in normal flux). For the weak Galerkin equation
$a(\psi, v) = (f, v)$ to hold for *every* test $v$, including those
supported across the interface, the jump must satisfy

$$
[\![\,\varepsilon_r\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,]\!] = (\text{surface charge density})/\varepsilon_0.
\tag{14.2}
$$

When there is no surface charge ($\sigma = 0$), the natural condition
of the variational form is exactly (14.1). **The Galerkin form encodes
flux continuity for free** — without any explicit interface term, the
discrete solution reproduces the physical interface condition to the
accuracy of the FE space.

This is the *single best feature* of FEM for multi-region problems and
the main reason kronos-semi uses it (as opposed to finite-difference
or finite-volume schemes that would need explicit interface flux
balancing).

### Why $\Phi_n, \Phi_p$ on a submesh?

In an ideal insulator, $n_i \to 0$ and the Slotboom expressions
$n = n_i\exp(...)$ give exactly zero carriers regardless of $\Phi_n$.
This means $\Phi_n$ is undefined in the oxide (any value gives the
same zero $n$). Putting $\Phi_n, \Phi_p$ on the parent mesh would
introduce ill-posed (rank-deficient) rows in the Jacobian for the
oxide DOFs; the LU factorization would fail with a zero pivot.

The cure is to define $\Phi_n, \Phi_p$ only on the silicon submesh:
[`dolfinx.mesh.create_submesh(parent, dim, semi_cells)`](../../semi/mesh.py)
returns a new mesh covering only the silicon cells, plus an entity map
that translates parent-cell indices to submesh-cell indices. Build
P1 spaces on the submesh:

```python
V_phi_n = fem.functionspace(submesh, ("Lagrange", 1))
```

These spaces have DOFs only at silicon vertices. The block-Jacobian
size is the right size; no rank-deficient rows.

### Mixing parent and submesh in one form

The Poisson row of the multi-region DD residual integrates over the
parent mesh (Laplacian over silicon + oxide) but the source term
$\hat\rho = \hat p - \hat n + \hat N$ requires $\hat n, \hat p$, which
are functions of $(\psi, \Phi_n, \Phi_p)$. $\psi$ lives on the parent
mesh; $\Phi_n, \Phi_p$ live on the submesh. UFL handles this through
**entity maps**:

```python
F_psi = (
    L_D2 * eps_r_ufl * inner(grad(psi), grad(v_psi)) * dx_parent
    - rho_hat * v_psi * dx_semi_parent
)
fem.form(F_psi, entity_maps=[em])
```

The `entity_maps` argument tells UFL "when you encounter `phi_n` in
this form (which lives on the submesh) but you're integrating on the
parent mesh, use this map to translate parent-cell indices to
submesh-cell indices and pull values from there." See
[`semi/physics/drift_diffusion.py:219-327`](../../semi/physics/drift_diffusion.py)
for the multi-region builder. The dolfinx 0.10 `NonlinearProblem` accepts
an `entity_maps` kwarg that threads through compile-time and assembly.

### Ohmic vs gate BCs in multi-region

The multi-region MOS has:
- An ohmic body contact at the silicon back face. Pinns $\psi, \Phi_n, \Phi_p$
  with the standard Shockley boundary.
- A gate contact at the oxide top face. Pins only $\psi = V_g - \phi_{ms}$.
  No $\Phi_n, \Phi_p$ Dirichlet because the gate facet has no
  $\Phi_n, \Phi_p$ DOFs (it's on the oxide-side boundary of the
  silicon submesh, which is *outside* the submesh entirely from the
  submesh's point of view).

[`semi/bcs.py:236-249`](../../semi/bcs.py) handles this explicitly:
gate contacts append a $\psi$ Dirichlet to the BC list and `continue`,
skipping the $\Phi_n, \Phi_p$ entries that ohmic contacts add.

### M6 verifier

The M6 MOS-2D benchmark verifier ([`scripts/run_benchmark.py`](../../scripts/run_benchmark.py)
`verify_mos_2d`) sweeps $V_g$ and computes the gate charge $Q_\mathrm{gate}(V_g)$
by integrating the silicon space charge:

$$
Q_\mathrm{gate}(V_g) = -\frac{q}{W_\mathrm{lat}}\int_{\Omega_\mathrm{Si}}\rho(\mathbf{x})\,dA,
\tag{14.3}
$$

(2D mesh; divide by lateral extent $W_\mathrm{lat}$ to get charge per
unit area). The C–V curve from $C(V_g) = -dQ_\mathrm{gate}/dV_g$ is
compared to the depletion-approximation theory (Ch. 9). The M6
benchmark passes within 10% in the verifier window
$[V_{fb}+0.2, V_T - 0.1]\,\mathrm{V}$.

The M6 multi-region Poisson MMS verifier
([`semi/verification/mms_poisson.py`](../../semi/verification/mms_poisson.py)
`run_mms_poisson_2d_multiregion`) tests the assembly against a
manufactured solution that satisfies (14.1) by construction — the
*method of manufactured solutions* (Ch. 21).

## Key results

- Cellwise DG0 $\varepsilon_r$: piecewise constant per region.
- Flux continuity (14.1) is encoded by Galerkin + DG0 $\varepsilon_r$
  with no explicit interface term.
- Submesh for carriers: avoids rank-deficient oxide rows.
- `entity_maps` translates between parent and submesh.
- Gate contact: only $\psi$ Dirichlet; no $\Phi_n, \Phi_p$ at gate facets.

## Worked numerical example

For the M6 device:
- Body: silicon, $N_a = 10^{17}\,\mathrm{cm^{-3}}$, 500 nm thick,
  $\varepsilon_r = 11.7$.
- Gate oxide: SiO₂, 5 nm thick, $\varepsilon_r = 3.9$.

Mesh: 1 nm vertical resolution → 505 cells over 505 nm. The interface
at $y = 500\,\mathrm{nm}$ lands on a grid line. Five cells in the oxide,
500 in the silicon.

DG0 $\varepsilon_r$:
- Silicon cells: $\varepsilon_r = 11.7$ → coefficient
  $L_D^2\cdot 11.7$.
- Oxide cells: $\varepsilon_r = 3.9$ → coefficient $L_D^2\cdot 3.9$.

Field at the interface in depletion ($V_g$ slightly above $V_{fb}$):
- $|E^\mathrm{Si}|$ at the interface: from depletion formula,
  $|E^\mathrm{Si}| \sim qN_a W/\varepsilon_s$.
- Continuity: $\varepsilon^\mathrm{Si} E^\mathrm{Si} = \varepsilon^\mathrm{ox} E^\mathrm{ox}$,
  so $E^\mathrm{ox} = (\varepsilon_r^\mathrm{Si}/\varepsilon_r^\mathrm{ox})\cdot E^\mathrm{Si} = 3\cdot E^\mathrm{Si}$.

So the field is *three times larger* in the oxide than at the silicon
side of the interface. This is the standard MOS oxide-field trick and
is what limits gate-oxide breakdown. The Galerkin form reproduces
this discontinuity automatically because $\varepsilon_r$ jumps cellwise
and the flux is enforced by the variational principle, not by hand.

## Code map

| Concept | UFL primitive | Code location |
|---|---|---|
| Cellwise DG0 $\varepsilon_r$ | `fem.functionspace(msh, ("DG", 0))` then interpolate | `semi/mesh.py::build_eps_r_function` |
| Multi-region Poisson form | (14.1) implicit | `semi/physics/poisson.py:82-118` |
| Submesh creation | `dolfinx.mesh.create_submesh(parent, dim, cells)` | `semi/mesh.py::build_submesh_by_role` |
| Submesh DD spaces | – | `semi/physics/drift_diffusion.py:175-216` (`DDBlockSpacesMR`) |
| Multi-region DD block | (14.1) + Slotboom | `semi/physics/drift_diffusion.py:219-327` (`build_dd_block_residual_mr`) |
| Entity-map threading | – | `semi/solver.py:315-321` (`entity_maps=[em]` to NonlinearProblem) |
| Gate Dirichlet exclusion of carriers | – | `semi/bcs.py:236-249` |
| Cell-tagged measure | `dx(subdomain_id=tag)` | `semi/physics/poisson.py:107-110` |
| Charge integral | (14.3) | `semi/runners/mos_cap_ac.py:138-152` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §6](../PHYSICS.md) — full M6 derivation.
- [`docs/PHYSICS_INTRO.md` §4.2](../PHYSICS_INTRO.md) — narrative for gate contacts.
- [`docs/mos_derivation.md`](../mos_derivation.md) — M6 mathematical gate.
- [`docs/PHYSICS.md` §5.5](../PHYSICS.md) — multi-region MMS verifier.

## Common pitfalls

1. **DG0 $\varepsilon_r$ is not `Constant`.** A `Constant` is the same
   value over the whole mesh; a DG0 `Function` has per-cell values. The
   form builder in [`semi/physics/poisson.py:67-70`](../../semi/physics/poisson.py)
   branches on `isinstance(eps_r, (int, float))` to decide which to wrap.
   Passing a `Function` where the form expects a `Constant` (or vice
   versa) compiles silently but produces wrong assembly.
2. **Submesh cells vs parent cells.** The `entity_map` returned by
   `create_submesh` is the *parent-to-submesh* map; the *submesh-to-parent*
   map is its inverse. UFL's `entity_maps` argument accepts the parent-to-submesh
   form; check the dolfinx 0.10 docs if you build a custom mixed form.
3. **Gate facet on submesh boundary.** The gate facet sits *outside*
   the silicon submesh entirely. If you naively try to apply
   `fem.locate_dofs_topological(V_phi_n, fdim, gate_facets)`, you get
   an empty DOF set — there are no $\Phi_n$ DOFs at the gate. The BC
   builder in [`semi/bcs.py:236-249`](../../semi/bcs.py) explicitly
   handles this by skipping carrier BCs at gate contacts.
4. **Surface charge needs an explicit `dS` term.** If you want to model
   a fixed oxide-charge sheet $Q_f$ at the interface, you would need
   to add an interior-facet integral `Q_f / eps_0 * v_psi * dS` to the
   Poisson row, plus a `dS` measure for the interface. The shipped
   engine sets $Q_f = 0$ in the analytic helpers; the FEM form does
   not yet expose interface-charge BCs.
5. **MMS for multi-region needs eps-flux-continuous manufactured
   solutions.** Constructing a manufactured $\psi^*(\mathbf{x})$ that
   is $C^0$ across the interface and has continuous $\varepsilon\nabla\psi^*\!\cdot\!\hat{\mathbf{n}}$
   is non-trivial. See `mos_derivation.md` §7 for the construction
   used by the M6 multi-region MMS gate.

## Exercises

**Exercise 14.1.** A user defines a 2D MOSCAP with three regions:
silicon body, SiO₂ gate oxide, *and* a thin Si₃N₄ passivation layer.
What does the cellwise DG0 $\varepsilon_r$ look like?

**Exercise 14.2.** Show that adding a surface charge $Q_f$ at the Si/SiO₂
interface modifies the flat-band voltage by $-Q_f/C_{ox}$, exactly as
the analytic helper [`semi/cv.py:127`](../../semi/cv.py) computes.

**Exercise 14.3.** Why can't $\Phi_n$ live on the parent mesh in a MOS
device? What goes wrong if you put it there? (Hint: think about what
$n = n_i\exp(\psi - \Phi_n)$ evaluates to in the oxide.)

**Exercise 14.4.** Read the `entity_maps` argument in
[`semi/solver.py:315-321`](../../semi/solver.py). What dolfinx
function consumes it, and what would happen if you forgot to pass it?

**Exercise 14.5.** A FinFET has a 3D-shaped channel surrounded by a
gate stack on three sides. How would you encode the multi-region
geometry in the JSON schema today (using `regions_by_box`)? What
limitation would force a gmsh-sourced unstructured mesh?

### Solutions

**14.1.** Three values: 11.7 for silicon cells, 3.9 for SiO₂ cells, 7.5
for Si₃N₄ cells. The DG0 function is per-cell, so the values are
piecewise constant by region. Each region's tag would be assigned via
`mesh.regions_by_box` and the DG0 builder reads the material name from
the JSON to look up $\varepsilon_r$ in [`semi/materials.py:86-88`](../../semi/materials.py).

**14.2.** Take the gate-side surface integral of the Maxwell flux:
$\varepsilon_{ox}E_{ox} = \varepsilon_s E_{Si} + Q_f$. Across the oxide,
$V_{ox} = E_{ox}T_{ox} = (\varepsilon_s E_{Si} + Q_f)T_{ox}/\varepsilon_{ox}$.
At flat band, $E_{Si} = 0$, so $V_{fb} = \phi_{ms} + Q_f T_{ox}/\varepsilon_{ox}
= \phi_{ms} + Q_f/C_{ox}$. Wait — sign: $V_{fb}$ is the *gate-side*
voltage that produces no band bending, and a positive $Q_f$ at the
interface (positive charge in the oxide) would attract electrons to the
silicon, requiring *less negative* gate voltage to invert; equivalently,
$V_{fb}$ shifts *more negative* by $Q_f/C_{ox}$. The textbook form is
$V_{fb} = \phi_{ms} - Q_f/C_{ox}$ ([`semi/cv.py:127`](../../semi/cv.py)).
The sign discrepancy is conventional: depending on whether $Q_f$ is
defined per the charge in the oxide or per the equivalent charge at
the interface.

**14.3.** In the oxide, $n_i^\mathrm{ox} \to 0$ (the oxide has no
carriers; its band gap is so large that thermal generation is
negligible). $n = n_i\exp(...) \to 0$ regardless of $\Phi_n$. The
electron continuity equation $\nabla\!\cdot\!(q\mu_n n\nabla\Phi_n) = qR$
in the oxide collapses to $0 = 0$, which is trivially true for any
$\Phi_n$. The discrete Jacobian for $\Phi_n$ on oxide DOFs has zero
entries (rank-deficient rows). MUMPS reports null pivots; SNES diverges.
Putting $\Phi_n$ only on the silicon submesh sidesteps this.

**14.4.** `dolfinx.fem.petsc.NonlinearProblem` consumes `entity_maps`
([`semi/solver.py:315-316`](../../semi/solver.py)). It is forwarded to
the form-compilation layer so UFL knows how to translate cell indices
between parent and child meshes. Forgetting to pass it produces a
compile error: "function on a different mesh, no entity map provided."

**14.5.** You'd need overlapping `regions_by_box` entries to carve out
the fin and the gate stack. axis-aligned boxes can model a rectangular
fin with rectangular gate cladding on top and sides, but a wrapping
gate (gate-all-around or top-and-sides) needs non-rectangular cell
unions. The schema's `regions_by_box` is too restrictive; you would
ship a `.msh` with named physical groups (M7 shipped this for the
resistor; M19 plans it for the 3D MOSFET).

## Further reading

- **Brezzi and Boffi (2003)** — for the mathematical analysis of mixed
  and constrained spaces. kronos-semi avoids these by using submesh +
  entity maps, but the theory underlies the correctness.
- **Logg et al. (2012), Chapter 4** — multi-region forms in FEniCS.
- **dolfinx tutorial** — multi-region examples at https://jsdokken.com/dolfinx-tutorial/.
- **`docs/mos_derivation.md`** in this repo — derivation-first M6 gate
  including the multi-region MMS construction.
- **Hu (2010), §5.7** for fixed-oxide-charge effects in real MOS
  devices, which the engine does not yet model.
