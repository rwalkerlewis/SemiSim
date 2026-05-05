# 13 — FEM weak forms

## Learning objectives

- Convert a strong-form PDE into its Galerkin weak form by multiplying
  with a test function and integrating by parts.
- Recognize the conforming Lagrange function spaces $H^1$ and $H^1_0$
  and what "P1" means in code.
- Map the JSON schema's `mesh.regions_by_box` and `mesh.facets_by_plane`
  to UFL `dx(tag)` and `ds(tag)` integration measures.
- Identify the natural-vs-essential boundary conditions in a UFL form.
- Read the kronos-semi residual builders and explain every term.

## Physical motivation

A device PDE in strong form ("at every interior point, this differential
equation holds; on the boundary, this condition is imposed") is the
natural physics statement. To solve it numerically on a complex domain
with arbitrary BCs, we convert to a **weak form**: integrate the strong
PDE against a test function over the domain, integrate by parts to
shift derivatives onto the test function, and discretize the resulting
integral equation by replacing the unknown and the test functions with
linear combinations of finite-element basis functions on a mesh.

The resulting nonlinear algebraic system has familiar structure (sparse
Jacobian, BCs lift to identity rows, integration is element-by-element)
and decades of mature solver infrastructure. dolfinx + UFL + PETSc
exposes exactly this pipeline. This chapter explains what every line
of `semi/physics/poisson.py` and `drift_diffusion.py` is doing in FEM
terms.

## Derivation from first principles

### Strong form to weak form: Poisson example

Strong form (1.5): $-\nabla\!\cdot\!(\varepsilon\nabla\psi) = \rho$ on
$\Omega$, with $\psi = g$ on $\Gamma_D$ and $\nabla\psi\!\cdot\!\hat{\mathbf{n}} = 0$
on $\Gamma_N$ (homogeneous Neumann).

Multiply by a test function $v$ that vanishes on $\Gamma_D$ and
integrate over $\Omega$:

$$
-\int_\Omega \nabla\!\cdot\!(\varepsilon\nabla\psi)\,v\,dV = \int_\Omega \rho\,v\,dV.
$$

Integrate the left side by parts:

$$
\int_\Omega \varepsilon\nabla\psi\!\cdot\!\nabla v\,dV
- \int_{\partial\Omega} \varepsilon\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,v\,dS
= \int_\Omega \rho\,v\,dV.
$$

The boundary integral splits into $\Gamma_D \cup \Gamma_N$. On
$\Gamma_D$, $v = 0$ (test function vanishes); the integral drops. On
$\Gamma_N$, $\nabla\psi\!\cdot\!\hat{\mathbf{n}} = 0$ (homogeneous Neumann);
the integral is zero. Both pieces vanish; the surface term *disappears*.
The Galerkin **weak form** is

$$
\boxed{
\int_\Omega \varepsilon\nabla\psi\!\cdot\!\nabla v\,dV
= \int_\Omega \rho\,v\,dV
\quad\forall v \in H^1_0(\Omega).
}
\tag{13.1}
$$

### Function spaces

The weak form lives in a **function space** of allowed $\psi, v$. The
choice that makes (13.1) well-posed is the Sobolev space $H^1(\Omega)$
of functions with square-integrable gradients. The test space adds the
constraint $v|_{\Gamma_D} = 0$, giving $H^1_0(\Gamma_D)$.

The discrete approximation replaces $H^1$ with a finite-dimensional
subspace, the **conforming Lagrange space $V_h \subset H^1$**:

- Triangulate $\Omega$ into elements (triangles, quads, tets).
- On each element, $\psi$ is a polynomial of degree $k$.
- Across elements, $\psi$ is continuous (and so is its definition; gradients
  jump but each pair of adjacent values agrees).

For $k = 1$ (P1 Lagrange), each element has one DOF per vertex and the
basis functions are the "tent functions" that equal 1 at one vertex
and 0 at all other vertices. kronos-semi uses P1 throughout
([`semi/physics/drift_diffusion.py:55-63`](../../semi/physics/drift_diffusion.py)
sets `("Lagrange", 1)`).

### Discretization

Expand $\psi(\mathbf{x}) = \sum_i u_i\phi_i(\mathbf{x})$ in the P1
basis $\{\phi_i\}$, and require (13.1) to hold for $v = \phi_j$ at each
interior DOF $j$:

$$
\sum_i u_i \underbrace{\int_\Omega \varepsilon\,\nabla\phi_i\!\cdot\!\nabla\phi_j\,dV}_{K_{ji}}
   = \underbrace{\int_\Omega \rho\,\phi_j\,dV}_{f_j}.
$$

The matrix $K$ is the **stiffness matrix**; the vector $f$ is the
**load vector**. For a nonlinear PDE (Poisson + Boltzmann), $\rho$ depends
on $\psi$, so the system is $K(\mathbf{u})\mathbf{u} = \mathbf{f}(\mathbf{u})$
or — better — written as a residual $F(\mathbf{u}) = \mathbf{0}$ that
Newton's method handles (Ch. 16).

### Dirichlet ("essential") BCs

$\psi = g$ on $\Gamma_D$ is not built into the test space; it is
imposed *after* assembly by:

1. Identify DOFs $j$ on $\Gamma_D$ via mesh-topology lookup
   ([`semi/bcs.py:180`](../../semi/bcs.py): `fem.locate_dofs_topological`).
2. Set $u_j = g_j$ in the solution vector.
3. Replace row $j$ of the matrix with the identity row $u_j = g_j$
   (via `dolfinx.fem.dirichletbc`).

This is the standard "lift-and-replace" treatment for essential BCs in
FEM.

### Neumann ("natural") BCs

$\nabla\psi\!\cdot\!\hat{\mathbf{n}} = h$ on $\Gamma_N$ is *not* an extra
constraint to impose; instead, it modifies the boundary integral term:

$$
\int_\Omega \varepsilon\nabla\psi\!\cdot\!\nabla v\,dV
- \int_{\Gamma_N} \varepsilon h\,v\,dS = \int_\Omega\rho v\,dV.
$$

For $h = 0$ (homogeneous Neumann, the engine's insulating BC), the
extra term is zero and you don't add it. Look at the UFL form: if you
see only `dx` integrals and no `ds`, you have implicitly assumed
homogeneous Neumann everywhere on $\Gamma_N$ ([`docs/PHYSICS.md` §3.3](../PHYSICS.md)).

### Region tags and the `dx(tag)` measure

For a multi-region device (silicon + oxide), the integration measure
splits per region:

```python
dx_full = ufl.Measure("dx", domain=msh)
dx_semi = ufl.Measure("dx", domain=msh,
                       subdomain_data=cell_tags,
                       subdomain_id=int(semi_tag))
```

Cells in the mesh carry an integer tag (`cell_tags` from `dolfinx.mesh.MeshTags`).
The schema's `mesh.regions_by_box` maps axis-aligned sub-boxes to tags;
[`semi/mesh.py`](../../semi/mesh.py) builds the `MeshTags` object from
that JSON. Then `dx_semi` integrates only over the cells with the
silicon tag; oxide cells contribute zero. This is exactly how the
multi-region MOS Poisson form
([`semi/physics/poisson.py:107-117`](../../semi/physics/poisson.py))
restricts space charge to silicon while letting Laplacian stiffness
integrate over the full mesh.

### Facet tags and the `ds(tag)` measure

Boundary facets carry tags from the schema's `mesh.facets_by_plane`
list. Each entry specifies an axis-aligned plane and a tag:

```json
{"name": "anode", "tag": 1, "axis": 0, "value": 0.0}
```

[`semi/mesh.py`](../../semi/mesh.py) collects facets with their
midpoints on the named plane and builds a `MeshTags` over the facets.
The runner uses these tags to (a) locate DOFs for Dirichlet BCs
([`semi/bcs.py:180`](../../semi/bcs.py)), (b) integrate terminal
currents over the facet via a `ds` measure
([`semi/postprocess.py:84-86`](../../semi/postprocess.py)).

### Submesh and `entity_maps`

For the multi-region MOS, $\psi$ lives on the parent mesh (silicon +
oxide); $\Phi_n, \Phi_p$ live only on the silicon submesh
([`semi/physics/drift_diffusion.py:173-216`](../../semi/physics/drift_diffusion.py)).
The submesh is created via `dolfinx.mesh.create_submesh(parent, dim, cells)`,
returning the new mesh plus an **entity map** that translates parent-mesh
cell indices to submesh cell indices. UFL forms that mix the two
(Poisson row reading $n, p$ from the submesh) compile by passing
`entity_maps=[em]` to `fem.form` and `NonlinearProblem`. Ch. 14 covers
this in depth; here the takeaway is "submesh + entity map = mixed-domain
assembly".

### Reading the engine's residual

Take the equilibrium-Poisson form
([`semi/physics/poisson.py:73-79`](../../semi/physics/poisson.py)):

```python
F = (
    L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v)) * ufl.dx
    - rho_hat * v * ufl.dx
)
```

Term by term:
- `L_D2 * eps_r_ufl * ufl.inner(grad(psi), grad(v)) * dx`:
  the stiffness term $\int L_D^2\varepsilon_r\nabla\hat\psi\!\cdot\!\nabla v\,dV$
  from (12.1).
- `- rho_hat * v * ufl.dx`: the source term, sign-flipped because the
  engine writes the form as $F = 0$ residual rather than $a = L$
  bilinear-linear pair.
- The `dx` is the default volume measure over the full mesh; no region
  restriction (single-region path).

Compare to the multi-region path (`build_equilibrium_poisson_form_mr`,
[`semi/physics/poisson.py:107-117`](../../semi/physics/poisson.py)):
- Stiffness: $\int_{\Omega} L_D^2\varepsilon_r(\mathbf{x})\nabla\hat\psi\!\cdot\!\nabla v\,dV$
  with $\varepsilon_r$ a cellwise DG0 Function (Ch. 14).
- Source: $\int_{\Omega_\mathrm{Si}} \hat\rho\,v\,dV$, restricted to
  silicon via `dx_semi`.
- Oxide cells contribute only the Laplacian (no charge), encoding
  Laplace's equation in the insulator.

## Key results

- Galerkin weak form: (13.1).
- Conforming Lagrange spaces: $V_h \subset H^1$, P1 = piecewise linear
  with one DOF per vertex.
- Stiffness + load decomposition: $K\mathbf{u} = \mathbf{f}$.
- Region/facet tags map JSON to UFL measures.
- Natural BCs are dropped surface integrals; Dirichlet BCs are imposed
  by row replacement.

## Worked numerical example

For the M1 device (1D, $L = 2\,\mu\mathrm{m}$, 400 cells), build the
P1 stiffness matrix $K_{ij} = \int_0^L (L_D^2\varepsilon_r)\,\phi_i'\,\phi_j'\,dx$
on a uniform mesh with $h = L/400 = 5\,\mathrm{nm}$.

The P1 hat function at vertex $i$ has $\phi_i' = +1/h$ on the left
element, $-1/h$ on the right element, and $0$ elsewhere. So
$\int \phi_i'\,\phi_i'\,dx = 2/h$ and $\int \phi_i'\,\phi_{i+1}'\,dx = -1/h$.

Per-element coefficient: $L_D^2\varepsilon_r = 1.67\times 10^{-16}$
(M1 worked example, Ch. 12). Multiplied:
- $K_{ii}^\mathrm{interior} = 2\cdot 1.67\times 10^{-16}/h = 3.34\times 10^{-16}/5\times 10^{-9}
  = 6.68\times 10^{-8}\,\mathrm{m^{-1}}$.
- $K_{i,i+1} = -1.67\times 10^{-16}/h = -3.34\times 10^{-8}$.

(These are the dimensional entries of $K$; the load vector has matching
dimensions.) The matrix is tridiagonal, symmetric, positive-definite,
and very well-conditioned with the M1 scaling — the M1 bug (Ch. 12)
showed how badly things go when this conditioning is lost.

## Code map

| Concept | UFL primitive | Code location |
|---|---|---|
| `Lagrange P1` space | `fem.functionspace(msh, ("Lagrange", 1))` | `semi/physics/drift_diffusion.py:57-59` |
| Test function | `ufl.TestFunction(V)` | `semi/physics/poisson.py:57` |
| Volume measure | `ufl.dx` | `semi/physics/poisson.py:75-77` |
| Region-restricted measure | `ufl.Measure("dx", subdomain_data=cell_tags, subdomain_id=tag)` | `semi/physics/poisson.py:107-110` |
| Facet measure for terminal current | `ufl.Measure("ds", subdomain_data=facet_mt, subdomain_id=tag)` | `semi/postprocess.py:84-86` |
| Stiffness assembly | `inner(grad(u), grad(v)) * dx` | `semi/physics/poisson.py:75` |
| Residual form $F$ | UFL composition | `semi/physics/poisson.py:73-78` |
| Dirichlet BC (essential) | `fem.dirichletbc(value, dofs, V)` | `semi/bcs.py:185` |
| Mesh + cell tags + facet tags builder | – | `semi/mesh.py` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §3.3](../PHYSICS.md) — "natural BC requires no explicit code" (the `ds`-less form).
- [`docs/PHYSICS.md` §6.2](../PHYSICS.md) — multi-region MOS Poisson with the `dx_full`/`dx_semi` split.
- [`docs/ARCHITECTURE.md`](../ARCHITECTURE.md) — Layer 4 (FEM) is where dolfinx imports live.
- [`docs/theory/dolfinx_choice.md`](../theory/dolfinx_choice.md) — why dolfinx 0.10's NonlinearProblem.
- **dolfinx tutorial** at https://jsdokken.com/dolfinx-tutorial/ — beginner-friendly intro to UFL.

## Common pitfalls

1. **Forgetting that $v = 0$ on $\Gamma_D$ in the test space.** When you
   apply a Dirichlet BC via `fem.dirichletbc`, the dolfinx assembler
   replaces the matrix row corresponding to that DOF with the identity
   row — equivalent to saying "the test function at this DOF is zero
   for the residual purposes; the unknown is fixed at $g$." If you
   manually try to add a constraint to the test function, you double-impose.
2. **Facet tag vs facet name.** The schema has `mesh.facets_by_plane[*].name`
   (string) and `mesh.facets_by_plane[*].tag` (int). [`semi/bcs.py:91-94`](../../semi/bcs.py)
   builds a `tag_by_name` map so contacts can refer to either by name
   or by integer. If you mix conventions, the resolver raises.
3. **`dx(tag)` vs `dx`.** Just `dx` integrates over the entire mesh.
   `dx(subdomain_id=tag)` integrates only over the tagged cells. Mixing
   them in the same form is fine — the multi-region MOS form uses both
   on the Poisson row.
4. **Mesh in physical meters vs scaled.** As Ch. 12 emphasized, the
   mesh stays in physical meters; UFL `grad(u)` returns the physical
   gradient (1/m). Forms must include explicit $L_0^2$ or $L_D^2$
   factors. This is captured in `semi/physics/`'s explicit constants.
5. **`fem.assemble_scalar` returns local, not global.** When integrating
   scalars (e.g. terminal currents), you have to `comm.allreduce` to
   get the MPI-aggregated value. [`semi/postprocess.py:103-106`](../../semi/postprocess.py)
   does this; [`semi/runners/mos_cap_ac.py:191-192`](../../semi/runners/mos_cap_ac.py)
   does too.

## Exercises

**Exercise 13.1.** Multiply the strong-form Poisson by a test function
$v$ that does *not* vanish on $\Gamma_D$. What surface term remains in
the weak form, and how does the dolfinx assembler handle it?

**Exercise 13.2.** Sketch the support of the P1 hat function $\phi_i$
in 1D and 2D. How many neighbouring DOFs is $\phi_i$ "non-orthogonal
to" in the stiffness matrix?

**Exercise 13.3.** A user wants to impose a non-homogeneous Neumann BC
$\nabla\psi\!\cdot\!\hat{\mathbf{n}} = h$ on a facet. Write the
modification to (13.1) and identify what UFL `ds` integral would need
to be added to the form.

**Exercise 13.4.** Read the multi-region Poisson form at
[`semi/physics/poisson.py:113-117`](../../semi/physics/poisson.py).
Why does the stiffness term integrate over `dx_full` while the source
term integrates over `dx_semi`? What does this enforce physically?

**Exercise 13.5.** Why does kronos-semi use *three separate* P1 scalar
spaces for $(\psi, \Phi_n, \Phi_p)$ rather than a single
`MixedElement`? See [`semi/physics/drift_diffusion.py:50-63`](../../semi/physics/drift_diffusion.py)
for the answer.

### Solutions

**13.1.** With $v$ free on $\Gamma_D$, the surface term
$\int_{\Gamma_D}\varepsilon\nabla\psi\!\cdot\!\hat{\mathbf{n}}\,v\,dS$
remains. The dolfinx assembler imposes the Dirichlet condition by row
replacement *after* assembly: the DOFs at $\Gamma_D$ are set to $g$
and the corresponding rows of $K$ become identity. The test space is
implicitly restricted to $\{v: v|_{\Gamma_D} = 0\}$ for the residual
purposes; the surface term is effectively zeroed even though it wasn't
formally dropped in the bilinear form.

**13.2.** 1D: $\phi_i$ is a tent over the two elements adjacent to vertex $i$.
Non-orthogonal to $\phi_{i-1}, \phi_i, \phi_{i+1}$ — three neighbours
(self-coupling counts).
2D triangles: support is a "star" of triangles around vertex $i$;
typically 6 neighbouring vertices on a uniform triangular mesh, plus
self.

**13.3.** Modified weak form:
$\int_\Omega \varepsilon\nabla\psi\!\cdot\!\nabla v\,dV
= \int_\Omega \rho\,v\,dV + \int_{\Gamma_N}\varepsilon h\,v\,dS$.
UFL: add `+ eps * h * v * ds(neumann_tag)` to the form. Currently
neither the schema nor `semi/bcs.py` exposes non-homogeneous Neumann;
adding it would be a small extension.

**13.4.** The Poisson stiffness term (Laplacian) acts on the *full*
mesh — both silicon and oxide regions carry $\psi$ and have a Laplacian.
The space-charge term acts only in *silicon* — the oxide has no carriers
and no doping, so $\hat\rho = 0$ there. Splitting the integration this
way encodes "Laplace's equation in the oxide, full Poisson in the
semiconductor" without any extra subdomain machinery on the LHS. See
[`docs/PHYSICS.md` §6.2](../PHYSICS.md).

**13.5.** Three independent scalar spaces simplify block assembly: the
Jacobian is naturally a 3×3 block, and the M1 equilibrium Poisson can
be re-used as the initial guess by extracting the $\psi$ subspace
without reformatting. A single MixedElement on the same mesh would
work too but couples the three unknowns at every DOF, complicating the
block-LU solve and the multi-region case where $\Phi_n, \Phi_p$ live
on a submesh and $\psi$ on the parent. The submesh case is impossible
to express as a MixedElement on a single mesh.

## Further reading

- **Logg, Mardal, and Wells, *Automated Solution of Differential
  Equations by the Finite Element Method* (2012).** The FEniCS book.
  Free online. Chapter 1 covers Galerkin from scratch and is exactly
  the right level for this guide. Pre-FEniCSx but most concepts carry
  over.
- **dolfinx tutorial: https://jsdokken.com/dolfinx-tutorial/**
  The current practical guide for dolfinx 0.10+; covers UFL
  expressions, forms, and the `MeshTags`/`Measure` API.
- **Brezzi and Boffi, *Mixed and Hybrid Finite Element Methods*
  (1991/2003).** For when MixedElement and the inf-sup condition
  matter. kronos-semi avoids these constructs by using $H^1$-conforming
  Lagrange on each scalar; this reference is for further study.
- **Brenner and Scott, *The Mathematical Theory of Finite Element
  Methods* (3rd ed., 2008).** Rigorous functional analysis of FEM.
  Heavier reading.
