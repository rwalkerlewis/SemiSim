# 12 — Nondimensionalization

## Learning objectives

- Explain the $10^{30}$-condition-number problem of the raw drift-diffusion
  Jacobian.
- Choose scales $L_0, V_t, C_0, \mu_0, t_0, J_0$ and apply them to derive
  the dimensionless Poisson equation with the small parameter
  $\lambda^2 = \varepsilon V_t/(qC_0L_0^2)$.
- Recognize $\lambda$ as the Debye-length-to-device ratio and interpret
  $\lambda^2 \ll 1$ as a singular perturbation.
- Reproduce the scaled drift-diffusion residual as it appears in
  `semi/physics/drift_diffusion.py:158-169`.
- Identify the M1 `lambda2` vs `L_D^2` trap and explain why the mesh
  staying in physical meters forces an explicit $L_0^2$ factor in
  every $\nabla\cdot\nabla$ term.

## Physical motivation

If you feed the raw drift-diffusion equations to Newton's method in SI
units, the Jacobian's largest entry has magnitude $\sim 10^{23}$ (carrier
density terms) and its smallest is $\sim 10^{-11}$ (permittivity). The
condition number is order $10^{30}$ and double-precision floating-point
arithmetic returns nonsense. Every device simulator ever written has
solved this with nondimensionalization: scale every variable so the
scaled versions are O(1). The kronos-semi engine's choice of scales is
the textbook one (Selberherr Ch. 7), with one twist: the spatial
coordinate stays in physical meters while everything else is scaled.
That choice has consequences which the M1 coefficient bug — preserved
as a cautionary note — illustrates vividly.

## Derivation from first principles

### The conditioning problem

Plug typical numbers into the dimensional Poisson equation
$-\nabla\!\cdot\!(\varepsilon\nabla\psi) = q(p - n + N)$:
- $\varepsilon \sim 10^{-11}\,\mathrm{F/m}$
- $q \sim 10^{-19}\,\mathrm{C}$
- $N \sim 10^{23}\,\mathrm{m^{-3}}$
- $L \sim 10^{-6}\,\mathrm{m}$

The discrete Laplacian has entries $\sim \varepsilon/h^2 \sim 10^{-11}/10^{-12}
= 10$. The discrete carrier-density rows (continuity equations) have entries
$\sim qN/V_t \sim 10^{-19}\cdot 10^{23}/0.025 = 4\times 10^5$. Block
ratio of the unknowns $\psi$ vs $n$ is $V_t/N \sim 10^{-25}$. Newton
sees a Jacobian whose diagonal entries span 30 orders of magnitude, and
both pivoting heuristics and ill-posedness destroy convergence.

### Scaling choice

The engine's scales (`Scaling` dataclass at [`semi/scaling.py:32-86`](../../semi/scaling.py)):

| Quantity | Scale | Default |
|---|---|---|
| Length | (kept in m) | – |
| Potential | $V_0 = V_t = kT/q$ | 25.85 mV at 300 K |
| Density | $C_0 = \max\|N\|$ | floor at $10^{16}\,\mathrm{cm^{-3}}$ |
| Mobility | $\mu_0 = \mu_n^\mathrm{ref}$ | reference electron mobility |
| Diffusivity | $D_0 = V_0\mu_0$ | Einstein at the reference |
| Time | $t_0 = L_0^2/D_0$ | sets the rate scale |
| Current density | $J_0 = qD_0C_0/L_0$ | natural for terminal currents |

Scaled variables are denoted with hats: $\hat\psi = \psi/V_t$,
$\hat n = n/C_0$, $\hat\mu = \mu/\mu_0$, etc. The mesh coordinate
remains $\mathbf{x}$ in meters; this is **Invariant 3** in `PLAN.md`.

### Scaled Poisson equation

Substitute $\psi = V_t\hat\psi$, $n = C_0\hat n$, etc. into (1.9):

$$
-\nabla\!\cdot\!(\varepsilon\,V_t\,\nabla\hat\psi) = q\,C_0\,(\hat p - \hat n + \hat N).
$$

Divide both sides by $qC_0$:

$$
\boxed{
-\nabla\!\cdot\!(L_D^2\,\varepsilon_r\,\nabla\hat\psi) = \hat p - \hat n + \hat N,
\qquad
L_D^2 \equiv \frac{\varepsilon_0 V_t}{q C_0}.
}
\tag{12.1}
$$

$L_D$ is the **extrinsic Debye length** at $C_0$: the length over which
the screening charge cloud around a unit perturbation decays. For Si at
$C_0 = 10^{17}\,\mathrm{cm^{-3}}$:

$$
L_D = \sqrt{\frac{\varepsilon_0 V_t}{q C_0}}
   = \sqrt{\frac{8.854\times 10^{-12}\cdot 0.02585}{1.602\times 10^{-19}\cdot 10^{23}}}
   = \sqrt{1.43\times 10^{-17}} = 3.78\,\mathrm{nm}.
$$

(Note: this is the bare $\sqrt{\varepsilon_0 V_t/(qC_0)}$; the often-quoted
"Debye length in silicon" of ~13 nm folds in $\varepsilon_r = 11.7$:
$\sqrt{\varepsilon_r}\cdot 3.78 = 12.9\,\mathrm{nm}$, matching
[`tests/check_analytical_math.py:30-33`](../../tests/check_analytical_math.py).)

### The dimensionless number $\lambda^2$

Define $\lambda^2 = L_D^2/L_0^2$. With $L_0 = 2\,\mu\mathrm{m}$
(M1 device length):

$$
\lambda^2 = \frac{L_D^2}{L_0^2}
   = \frac{(3.78\times 10^{-9})^2}{(2\times 10^{-6})^2}
   = \frac{1.43\times 10^{-17}}{4\times 10^{-12}}
   = 3.57\times 10^{-6}.
$$

Multiplied by $\varepsilon_r = 11.7$, this is $4.18\times 10^{-5}$;
matches [`tests/check_analytical_math.py:27-29`](../../tests/check_analytical_math.py).
The `Scaling.lambda2` property returns the bare $\lambda^2$ without
$\varepsilon_r$ ([`semi/scaling.py:62-70`](../../semi/scaling.py)).

### The mesh-stays-in-meters subtlety

If you naively rewrote (12.1) with the gradient operator scaled
($\nabla \to \nabla/L_0$), you would get $-\lambda^2\varepsilon_r\Delta\hat\psi
= \hat\rho$ with $\Delta = $ scaled Laplacian. But the kronos-semi
mesh stays in physical meters, so $\nabla$ in the UFL form is the
physical gradient (1/m). The coefficient is therefore $L_D^2 = \lambda^2 L_0^2$:

$$
\text{stiffness coefficient in code: }
   L_D^2 \cdot \varepsilon_r = \lambda^2 \cdot L_0^2 \cdot \varepsilon_r.
$$

Concretely: [`semi/physics/poisson.py:60-66`](../../semi/physics/poisson.py)
sets `L_D2 = sc.lambda2 * sc.L0**2` and uses `L_D2 * eps_r_ufl` in the
UFL form. Using `lambda2` directly would suppress the diffusion term
by $L_0^2 \sim 10^{-12}$ on a micron device — that is the M1
coefficient bug, preserved as a cautionary note in
[`docs/PHYSICS.md` §2.3](../PHYSICS.md). Newton happened to converge
on the *Laplace* solution (no nonlinear charge), $V_{bi}$ matched
because it is set by the BC, but peak field was off by 20×.

### Scaled drift-diffusion

Following the same recipe for the continuity rows, all the algebra in
[`docs/PHYSICS.md` §2.5](../PHYSICS.md):

$$
-\nabla\!\cdot\!(L_D^2\varepsilon_r\nabla\hat\psi) = \hat p - \hat n + \hat N
\tag{12.2a}
$$

$$
-\nabla\!\cdot\!(L_0^2\,\hat\mu_n\,\hat n\,\nabla\hat\Phi_n) = +\hat R
\tag{12.2b}
$$

$$
-\nabla\!\cdot\!(L_0^2\,\hat\mu_p\,\hat p\,\nabla\hat\Phi_p) = -\hat R
\tag{12.2c}
$$

with $\hat n, \hat p$ from Slotboom (11.3), and the scaled SRH kernel

$$
\hat R(\hat n, \hat p) = \frac{\hat n\hat p - \hat n_i^2}{\hat\tau_p(\hat n+\hat n_1) + \hat\tau_n(\hat p+\hat p_1)}.
\tag{12.2d}
$$

The $L_0^2$ in (12.2b)–(12.2c) is the same "mesh stays in meters" pickup;
each $\nabla\cdot\nabla$ contributes one $1/L_0^2$, which is multiplied
back to give the explicit $L_0^2$. Both $L_D^2$ and $L_0^2$ are computed
in the runner from `sc.lambda2 * sc.L0**2` and `sc.L0**2` respectively
([`semi/physics/drift_diffusion.py:131-132`](../../semi/physics/drift_diffusion.py)).

### Singular-perturbation interpretation

$\lambda^2 \ll 1$ means the diffusive term in (12.1) is small compared
with the source. This is a **singular perturbation**: in regions where
the source dominates ($\hat\rho = O(1)$), the diffusive term is
negligible and the solution is *locally constant*. In thin transition
regions of width $\sim L_D$, the diffusive term restores continuity
with the boundary conditions and the solution rapidly varies. Those
thin transition regions are the **depletion regions** at junctions and
contacts — exactly the regime that Ch. 7 derives in the depletion
approximation.

The singular-perturbation framing makes precise why depletion regions
are narrow: they have width $\sim L_D \sim 13\,\mathrm{nm}$ at $10^{17}$
doping, while quasi-neutral bulk regions can be hundreds of microns.

## Key results

- Scaled Poisson coefficient: $L_D^2\,\varepsilon_r$, with
  $L_D^2 = \lambda^2 L_0^2 = \varepsilon_0 V_t/(qC_0)$.
- Scaled DD residual: (12.2).
- Singular-perturbation bound: width of transition regions $\sim L_D$.
- Mesh-in-meters convention: every $\nabla\!\cdot\!\nabla$ picks up an
  explicit $L_0^2$.

## Worked numerical example

For the M1 `pn_1d` benchmark:
- $L_0 = 2\,\mu\mathrm{m}$, $V_t = 25.85\,\mathrm{mV}$,
  $C_0 = 10^{17}\,\mathrm{cm^{-3}} = 10^{23}\,\mathrm{m^{-3}}$.
- `lambda2` $= \varepsilon_0 V_t/(qC_0L_0^2)
  = 8.854\times 10^{-12}\cdot 0.02585/(1.602\times 10^{-19}\cdot 10^{23}\cdot 4\times 10^{-12})
  = 2.288\times 10^{-13}/6.41\times 10^{-8} = 3.57\times 10^{-6}$. ✓
- $L_D^2 = \lambda^2 L_0^2 = 3.57\times 10^{-6}\cdot 4\times 10^{-12} = 1.43\times 10^{-17}\,\mathrm{m^2}$.
- Stiffness coefficient: $L_D^2 \cdot \varepsilon_r = 1.67\times 10^{-16}$.

Compare to the pre-fix bug: using `lambda2` directly:
$\lambda^2 \cdot \varepsilon_r = 4.18\times 10^{-5}$, a factor of
$\sim 10^{12}$ too large for the *coefficient* (and equivalently,
diffusion is suppressed by $10^{12}$ relative to the correct value when
$\nabla$ is interpreted in physical units). That mismatch is the
M1 bug.

**$L_0^2$ in the continuity rows.** Same $L_0 = 2\,\mu\mathrm{m}$ gives
$L_0^2 = 4\times 10^{-12}\,\mathrm{m^2}$. The continuity stiffness in
(12.2b) is $L_0^2\hat\mu\hat n$; with $\hat\mu \approx 1$ and $\hat n
\approx 1$ at full doping, the coefficient is again $\sim 10^{-12}$,
balancing against the dimensionless $\hat R$ (which is normalized to
$C_0/t_0$).

**$t_0$ for the transient runner.** $D_0 = V_t\mu_0 = 0.02585\cdot 0.14 = 3.62\times 10^{-3}\,\mathrm{m^2/s}$.
$t_0 = L_0^2/D_0 = 4\times 10^{-12}/3.62\times 10^{-3} = 1.10\times 10^{-9}\,\mathrm{s} = 1.1\,\mathrm{ns}$.
A typical $\tau_n = 100\,\mathrm{ns}$ scaled is $\hat\tau_n = 100/1.1 = 90.9$.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| `Scaling` dataclass | – | `semi/scaling.py:32-86` |
| `Scaling.V0` | – | `semi/scaling.py:42-45` |
| `Scaling.D0` | – | `semi/scaling.py:48-50` |
| `Scaling.t0` | – | `semi/scaling.py:53-55` |
| `Scaling.J0` | – | `semi/scaling.py:58-60` |
| `Scaling.lambda2` | $\lambda^2$ | `semi/scaling.py:62-70` |
| `Scaling.debye_length` | $L_D$ | `semi/scaling.py:72-80` |
| `make_scaling_from_config` | – | `semi/scaling.py:89-103` |
| Stiffness coefficient $L_D^2\cdot\varepsilon_r$ | (12.1), (12.2a) | `semi/physics/poisson.py:60-66`, `drift_diffusion.py:131` |
| Continuity coefficient $L_0^2\cdot\hat\mu\cdot\hat n$ | (12.2b) | `semi/physics/drift_diffusion.py:132, 162-163` |
| M1 trap cautionary note | – | `docs/PHYSICS.md` §2.3 |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §2](../PHYSICS.md) — full scaled DD derivation including the trap note.
- [`docs/theory/scaling.md`](../theory/scaling.md) — concise "why".
- [`docs/adr/0002-nondimensionalization-mandatory.md`](../adr/0002-nondimensionalization-mandatory.md) — the locked decision.

## Common pitfalls

1. **Using `lambda2` where `L_D^2 = lambda2 * L0^2` is needed.** This is
   the M1 bug. Modern code should follow [`semi/physics/poisson.py:60-66`](../../semi/physics/poisson.py)
   verbatim: `L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0**2))`.
2. **Mixing scaled and dimensional quantities.** Inside `semi/physics/`
   everything is in scaled units (hats); inside `semi/postprocess/` and
   the runners' return values, everything is in SI. Crossing the boundary
   without rescaling produces silent bugs.
3. **`C_0` floor at $10^{16}\,\mathrm{cm^{-3}}$.** [`semi/scaling.py:113-115`](../../semi/scaling.py)
   floors `C_max` at $10^{16}$. This means a uniform-doping device with
   $N = 10^{15}$ scales to $\hat N = 0.1$, not 1. The floor exists so
   the engine does not divide by zero on undoped substrates.
4. **`L_0` from mesh extents.** [`semi/scaling.py:106-111`](../../semi/scaling.py)
   takes $L_0$ from the *largest* mesh extent for builtin meshes; for
   file-based meshes it falls back to $1\,\mu\mathrm{m}$. If your gmsh
   mesh has a 100 µm-long resistor, the fallback under-estimates $L_0$
   by 100×; expect conditioning issues. The right fix is to add an
   explicit `L0_override` schema field, which has not yet been
   prioritized.
5. **Time scaling has its own subtlety.** $t_0 = L_0^2/D_0$ is the
   diffusion time scale. Real device timescales (the RC charging time
   of an MOS capacitor, the SRH lifetime) are different orders of
   magnitude. The scaled lifetime $\hat\tau = \tau/t_0$ can be very
   large (hundreds), which is fine, but the scaled time step $\hat{dt}
   = dt/t_0$ has to fall inside a reasonable range or BDF1/BDF2 stalls.

## Exercises

**Exercise 12.1.** Compute $\lambda^2$ for the M6 MOSCAP body
($N_A = 10^{17}\,\mathrm{cm^{-3}}$) on a 500 nm-thick silicon body.

**Exercise 12.2.** Show that the dimensional Debye length
$\sqrt{\varepsilon V_t/(qC_0)}$ comes out 13 nm in silicon at
$10^{17}\,\mathrm{cm^{-3}}$, matching the [`tests/check_analytical_math.py:30-33`](../../tests/check_analytical_math.py)
assertion.

**Exercise 12.3.** Derive (12.2b) from the dimensional electron
continuity equation $\nabla\!\cdot\!(q\mu_n n\,\nabla\Phi_n) = qR$ by
substituting $\Phi_n = V_t\hat\Phi_n$, $n = C_0\hat n$, $\mu_n = \mu_0\hat\mu_n$,
and dividing by $qC_0/t_0$.

**Exercise 12.4.** What does $\lambda^2 \to 0$ mean physically? Sketch
the limit of (12.1) and connect it to the depletion-approximation
charge-sheet picture from Ch. 7.

**Exercise 12.5.** A user submits a JSON with a 10 µm-long device but
a doping of $10^{14}\,\mathrm{cm^{-3}}$. The engine floors $C_0$ at
$10^{16}\,\mathrm{cm^{-3}}$. What value of $\lambda^2$ does the scaling
produce, and why is it different from what you would compute from the
actual doping?

### Solutions

**12.1.** $L_0 = 500\,\mathrm{nm}$. $\lambda^2 = 8.854\times 10^{-12}\cdot 0.02585
/(1.602\times 10^{-19}\cdot 10^{23}\cdot (5\times 10^{-7})^2)
= 2.288\times 10^{-13}/(1.602\times 10^4\cdot 2.5\times 10^{-13})
= 2.288\times 10^{-13}/4.005\times 10^{-9}
= 5.71\times 10^{-5}$.
Smaller device → larger $\lambda^2$ at fixed doping. Multiplied by
$\varepsilon_r = 11.7$: $6.69\times 10^{-4}$.

**12.2.** $L_D^\mathrm{absolute} = \sqrt{\varepsilon_r\varepsilon_0V_t/(qC_0)}
= \sqrt{11.7\cdot 8.854\times 10^{-12}\cdot 0.02585/(1.602\times 10^{-19}\cdot 10^{23})}
= \sqrt{1.679\times 10^{-16}} = 1.296\times 10^{-8}\,\mathrm{m} = 12.96\,\mathrm{nm}$. ✓
The test asserts $10 < L_D/\mathrm{nm} < 16$.

**12.3.** Substitute: $\nabla\!\cdot\!(q\mu_0\hat\mu_n C_0\hat n V_t\nabla\hat\Phi_n)
= qR$. Pull constants out: $q\mu_0 V_t C_0\nabla\!\cdot\!(\hat\mu_n\hat n\nabla\hat\Phi_n) = qR$.
Divide by $qC_0/t_0 = qC_0D_0/L_0^2$: $\mu_0 V_t/(D_0/L_0^2)\cdot\nabla\!\cdot\!(\hat\mu_n\hat n\nabla\hat\Phi_n) = \hat R$.
With $D_0 = V_t\mu_0$: $L_0^2\nabla\!\cdot\!(\hat\mu_n\hat n\nabla\hat\Phi_n) = \hat R$. ✓

**12.4.** $\lambda^2 \to 0$ means the stiffness term in (12.1)
vanishes; the equation reduces to $\hat\rho = \hat p - \hat n + \hat N = 0$
— pointwise charge neutrality. This is the depletion-approximation
limit *outside* the transition layer: bulk regions are quasi-neutral.
Inside the transition layer (width $\sim L_D \to 0$), the diffusive
term restores continuity by enforcing rapid $\hat\psi$ variation.
The depletion-approximation charge-sheet picture is exact in the
$\lambda^2 \to 0$ limit.

**12.5.** $C_0 = 10^{16}\,\mathrm{cm^{-3}}$ (floored), $L_0 = 10\,\mu\mathrm{m}$.
$\lambda^2 = 8.854\times 10^{-12}\cdot 0.02585/(1.602\times 10^{-19}\cdot 10^{22}\cdot 10^{-10})
= 2.288\times 10^{-13}/(1.602\times 10^3\cdot 10^{-10}\cdot 10^{-12})$
Wait, recompute: $C_0L_0^2 = 10^{22}\cdot 10^{-10} = 10^{12}$;
$qC_0L_0^2 = 1.602\times 10^{-7}$.
$\lambda^2 = 2.288\times 10^{-13}/1.602\times 10^{-7} = 1.43\times 10^{-6}$.
With *actual* doping $10^{14}\,\mathrm{cm^{-3}} = 10^{20}\,\mathrm{m^{-3}}$:
$\lambda^2 = 2.288\times 10^{-13}/(1.602\times 10^{-19}\cdot 10^{20}\cdot 10^{-10})
= 2.288\times 10^{-13}/1.602\times 10^{-9} = 1.43\times 10^{-4}$.
Two orders of magnitude smaller because the floor over-rates $C_0$.
The scaled equation is solved with the floored $C_0$; the result is
correctly recovered in physical units after un-scaling, but the
condition number is sub-optimal during the solve.

## Further reading

- **Selberherr, *Analysis and Simulation of Semiconductor Devices*
  (1984), Chapter 7.** The standard reference for nondimensionalization
  in DD. The kronos-semi scaling matches Selberherr's choice.
- **Brezzi, Marini, Pietra (1989).** "Two-dimensional exponential
  fitting and applications to drift-diffusion models." Discusses why
  $\lambda^2 \ll 1$ is a singular perturbation and how that drives
  discretization choices.
- **Markowich, *The Stationary Semiconductor Device Equations* (1986).**
  Mathematical analysis of the singular-perturbation limit and the
  depletion-approximation derivation.
- **ADR 0002** in this repo for the engine's stance.
- [`docs/theory/scaling.md`](../theory/scaling.md) for a one-page
  summary aimed at code readers.
