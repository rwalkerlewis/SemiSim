# 17 — Transient time integration

## Learning objectives

- Derive backward Euler (BDF1) and BDF2 time discretizations for the
  drift-diffusion + Poisson system in Slotboom variables.
- Explain why implicit methods are mandatory: the carrier-time-scale
  stiffness ratio in semiconductor problems.
- Read the BDF coefficient class in `semi/timestepping.py` and the
  Slotboom transient residual in `semi/runners/transient.py`.
- Recognize the BC-ramp continuation for the initial condition (ADR 0013)
  and the Jacobian shift for minority-side stability (ADR 0014).
- Anticipate the M16.7 time-varying contact-voltage extension and its
  audit-case-06 verification target.

## Physical motivation

A time-dependent device simulation — diode turn-on, MOSFET switching,
RC ringing — needs the time-derivative term in the continuity equations.
Carrier dynamics span many orders of magnitude in timescale: SRH lifetimes
(100 ns), depletion-region transit (sub-ps), dielectric relaxation
($\sim L_D^2/(D \mu_0)$, picoseconds). Adding the time derivative makes
the system *stiff*: the fast modes are locally stable but require tiny
time steps for explicit methods. Implicit methods (BDF1, BDF2) are
A-stable for $\mathrm{Re}(\lambda dt) > 0$ regardless of step size, so
they integrate the slow modes accurately while letting the fast modes
relax to their slow manifolds.

ADR 0010 commits to BDF1/BDF2 with fixed $dt$ for the M13 transient
runner. ADR 0014 (which supersedes ADR 0009) commits to Slotboom
primary unknowns for the transient runner — the same variables as
`bias_sweep` — to avoid the M13.x discrete-residual mismatches that
plagued the original (n, p) primary-form attempt.

## Derivation from first principles

### Continuity in time

The full transient continuity equation is

$$
\frac{\partial n}{\partial t} + \nabla\!\cdot\!\mathbf{J}_n/(-q) = G - R,
\qquad
\frac{\partial p}{\partial t} + \nabla\!\cdot\!\mathbf{J}_p/(+q) = G - R,
\tag{17.1}
$$

with the same drift-diffusion currents (5.3) and recombination (Ch. 6).
Poisson is purely algebraic in time (instantaneous):
$-\nabla\!\cdot\!(\varepsilon\nabla\psi) = \rho(t)$.

In Slotboom variables, $n$ and $p$ are *functions* of $\psi, \Phi_n, \Phi_p$.
The time derivative becomes

$$
\frac{\partial n}{\partial t}
= \frac{\partial n}{\partial\psi}\cdot\frac{\partial\psi}{\partial t}
+ \frac{\partial n}{\partial\Phi_n}\cdot\frac{\partial\Phi_n}{\partial t}
= \frac{n}{V_t}\left(\frac{\partial\psi}{\partial t} - \frac{\partial\Phi_n}{\partial t}\right).
\tag{17.2}
$$

Thus the time derivative *couples* $\psi$ and $\Phi_n$ through a chain
rule, which is the **chain-rule mass matrix** that ADR 0014 calls out.
Despite the coupling, the implementation is automatic: UFL's
`derivative` facility produces the correct sparse mass-matrix block
when the discrete time term is written as $\alpha_0/dt\cdot n_\mathrm{ufl}\cdot v\,dx$
and SNES automatically picks up the chain-rule entries.

### BDF1 (backward Euler)

The simplest implicit scheme: approximate

$$
\frac{\partial n}{\partial t}\bigg|_{t^{n+1}} \approx \frac{n^{n+1} - n^n}{dt}.
\tag{17.3}
$$

The residual at $t^{n+1}$ is

$$
F_n^\mathrm{BDF1} = \frac{n^{n+1} - n^n}{dt} + \nabla\!\cdot\!\mathbf{J}_n^{n+1}/(-q) + R^{n+1} = 0.
\tag{17.4}
$$

Convergence rate: first-order in time, $O(dt)$. A-stable.

### BDF2

Use a three-level history with quadratic backward extrapolation:

$$
\frac{\partial n}{\partial t}\bigg|_{t^{n+1}}
\approx \frac{1}{dt}\left(\frac{3}{2}n^{n+1} - 2n^n + \frac{1}{2}n^{n-1}\right).
\tag{17.5}
$$

Convergence rate: $O(dt^2)$. A-stable for $\theta$ in $[0, 90°)$ relative
to the imaginary axis (equivalently, A($\alpha$)-stable with $\alpha \approx 90°$,
see Hairer & Wanner). For most semiconductor problems where eigenvalues
are real and negative, this is fine.

### General BDF formula

Generic BDF-$k$:

$$
\frac{\partial u}{\partial t}\bigg|_{t^{n+1}}
\approx \frac{1}{dt}\sum_{j=0}^{k}\alpha_j u^{n+1-j},
\tag{17.6}
$$

with $\alpha_0 > 0$ the leading coefficient. For BDF1: $(1, -1)$.
For BDF2: $(3/2, -2, 1/2)$. [`semi/timestepping.py:54-65`](../../semi/timestepping.py)
encodes both.

### Engine implementation: Slotboom + chain rule + lumped mass

The transient runner builds the residual
([`semi/runners/transient.py:528-668`](../../semi/runners/transient.py)):

- Poisson row: identical to bias_sweep (no time term, since Poisson is
  instantaneous).
- Electron continuity row: adds
  $\int (\alpha_0/dt) n_\mathrm{ufl} v_n\,dx_\mathrm{lump}
  + \int (f_\mathrm{hist,n}/dt) v_n\,dx_\mathrm{lump}$
  to the steady-state Slotboom continuity, where $f_\mathrm{hist,n}$ is
  a stored Function holding $\sum_{k\ge 1}\alpha_k n^{n+1-k}$.
- Hole continuity row: same with $p_\mathrm{ufl}$.

Two important details:

1. **Lumped mass**: the time term integrates against the **vertex
   quadrature** measure
   $dx_\mathrm{lump} = \mathrm{ufl.dx(metadata=\{"quadrature_rule":"vertex","quadrature_degree":1\})}$
   instead of the full Galerkin measure. Lumping localizes the
   per-DOF mass to a single vertex (diagonal mass matrix), which
   prevents the consistent (Galerkin) mass from coupling DOFs of widely
   different magnitudes. ADR 0014 spells out why: at small $dt$, the
   consistent mass would dominate the spatial Laplacian and push MUMPS
   into pivot failure; lumping makes the Jacobian's mass block
   diagonal (in the (phi_n, phi_n) sub-block) and well-conditioned.

2. **History storage**: $n^{n+1-k}, p^{n+1-k}$ are stored as numpy
   arrays per converged Slotboom step, evaluated via
   `n_from_slotboom_np(psi, phi_n, ni_hat)` from
   [`semi/physics/slotboom.py:45-52`](../../semi/physics/slotboom.py).
   The history-source Functions $f_\mathrm{hist,n}, f_\mathrm{hist,p}$
   are updated each timestep with $\sum_k\alpha_k n^{n+1-k}$ and read
   into the residual ([`semi/runners/transient.py:434-443`](../../semi/runners/transient.py)).

### BC-ramp continuation initial condition

The transient runner uses a two-stage IC strategy (ADR 0013):

1. Solve $V = 0$ equilibrium Poisson to obtain $\psi_\mathrm{eq}$.
2. If `bc_ramp_steps > 0` (default 10), solve a steady-state
   `bias_sweep` from $V = 0$ to $V_\mathrm{target}$ and use the
   converged Slotboom triple $(\psi, \Phi_n, \Phi_p)$ as the time-loop IC.

For the `pn_1d_turnon` benchmark, `bc_ramp_steps = 0` is the right choice:
the V=0 equilibrium IC + step-bias-at-$t=0$ *is* the physical scenario
being measured. For the deep-steady-state-limit test, `bc_ramp_steps = 10`
puts the time loop's IC near its fixed point, so the steady-state
matching gate (1e-4 relative) is achievable.

### Jacobian shift

In the deep p-bulk, $n = n_i\exp(\psi - \Phi_n)$ is below floating-point
precision (e.g. $\sim 10^{-300}$ for $\Phi_n - \psi \gg 30 V_t$). The
Jacobian's $(\Phi_n, \Phi_n)$ block is then *numerically* zero on those
DOFs; MUMPS reports null pivots; the Newton update has unbounded
components that send subsequent residual evaluations to NaN.

Adding $\epsilon I$ to the assembled Jacobian (with $\epsilon \sim 10^{-14}$)
removes the null pivots without measurably perturbing the converged
solution. The transient runner sets `jacobian_shift = 1e-14`
([`semi/runners/transient.py:184-185`](../../semi/runners/transient.py))
and `mat_mumps_icntl_14: 200` (200% extra MUMPS workspace) to absorb
the small extra delayed pivots. ADR 0014 §"Solver options" documents
this.

### `pn_1d_turnon` benchmark

[`benchmarks/pn_1d_turnon/`](../../benchmarks/pn_1d_turnon/) is the
M13 / M13.1 acceptance benchmark: a 1D pn diode initially at
$V_F = 0\,\mathrm{V}$ is forward-biased to $V_F = 0.6\,\mathrm{V}$ at
$t = 0$, and the time-evolved drain current $I(t)$ is recorded. The
verifier extracts the effective minority-carrier lifetime $\tau_\mathrm{eff}$
from the late-time exponential approach to steady state, and compares
to the input $\tau_p = 100\,\mathrm{ns}$ within 5%. ADR 0014 reports
that the M13.1 Slotboom transient is within this gate.

### Forward reference: M16.7 time-varying $V(t)$

The shipped runner takes the contact voltages from
`cfg["contacts"][i]["voltage"]`, evaluated *once* and held fixed
throughout the time loop. The M16.7 milestone will extend this with a
per-step contact-voltage callable / table, allowing simulations like
"step-and-ring" or sinusoidal AC drive in the time domain. The audit
case 06 (currently `pytest.skip`) verifies that the transient FFT of
$I(t)$ at small-signal sinusoidal $V(t)$ agrees with the AC sweep
$Y(\omega)$ at the same frequency within 5%. See
[`docs/IMPROVEMENT_GUIDE.md` §M16.7](../IMPROVEMENT_GUIDE.md).

## Key results

- BDF1: (17.3)–(17.4); first-order, A-stable.
- BDF2: (17.5); second-order, A($\alpha$)-stable.
- Slotboom chain-rule time term: (17.2).
- Lumped mass: $dx_\mathrm{lump}$ with vertex quadrature.
- Jacobian shift: $\epsilon I$ with $\epsilon \sim 10^{-14}$.
- IC strategies: V=0 equilibrium (`bc_ramp_steps = 0`) or BC-ramp
  to $V_\mathrm{target}$.

## Worked numerical example

For `pn_1d_turnon` at $V_F = 0.6\,\mathrm{V}$, $\tau_p = 100\,\mathrm{ns}$,
$dt = 1\,\mathrm{ns}$, BDF2 order, $t_\mathrm{end} = 500\,\mathrm{ns}$:

- Number of timesteps: 500.
- Time-step ratio: $dt/\tau_p = 0.01$. Resolves the SRH lifetime
  comfortably.
- Scaling: $t_0 = L_0^2/D_0 = 4\times 10^{-12}/3.62\times 10^{-3} = 1.10\,\mathrm{ns}$.
  Scaled $\hat{dt} = dt/t_0 = 1\,\mathrm{ns}/1.1\,\mathrm{ns} = 0.91$.
  Scaled $\hat\tau_p = 100/1.1 = 90.9$.

BDF2 coefficients $(\alpha_0, \alpha_1, \alpha_2) = (3/2, -2, 1/2)$. At
step $n+1$: residual = $(3/2/\hat{dt}) n^{n+1}_\mathrm{ufl} v_n\,dx_\mathrm{lump}
+ (-2 n^n + n^{n-1}/2)/\hat{dt}\cdot v_n\,dx_\mathrm{lump} + \mathrm{spatial}$
$+ \hat R\,v_n\,dx$.

The first step uses BDF1 (only one history value available); the
second step onwards uses BDF2. The auto-switch lives at
[`semi/runners/transient.py:413-423`](../../semi/runners/transient.py).

Late-time behaviour: at $t \to \infty$, the system relaxes to its
$V_F = 0.6\,\mathrm{V}$ steady state. The relaxation rate is
$\sim 1/\tau_p$. Fitting an exponential to the late-time $I(t) - I_\infty$
gives $\tau_\mathrm{eff}$. Match to $\tau_p$ within 5% is the M13
acceptance gate.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| `BDFCoefficients(order)` | (17.6) | `semi/timestepping.py:28-96` |
| BDF1 / BDF2 coefficients | – | `semi/timestepping.py:54-57` |
| `apply` (history sum) | – | `semi/timestepping.py:67-96` |
| Slotboom transient residual | (17.4)+(17.5) hybrid | `semi/runners/transient.py:528-668` (`_build_transient_residual`) |
| Lumped mass measure | – | `semi/runners/transient.py:639-642` |
| BDF1 → BDF2 switching | – | `semi/runners/transient.py:413-423` |
| History storage | – | `semi/runners/transient.py:353-360, 434-443` |
| BC-ramp continuation IC | (Ch. 16) | `semi/runners/transient.py:671-789` (`_run_bc_continuation`) |
| Jacobian shift | – | `semi/solver.py:121-149` (`_install_jacobian_shift`) |
| MUMPS workspace bump | – | `semi/runners/transient.py:159-164` |
| `pn_1d_turnon` benchmark | – | `benchmarks/pn_1d_turnon/` |
| Time-varying $V(t)$ (planned) | – | `docs/IMPROVEMENT_GUIDE.md` §M16.7 |

## Existing-docs cross-reference

- [`docs/adr/0009-transient-formulation.md`](../adr/0009-transient-formulation.md) — superseded; original (n,p) form choice.
- [`docs/adr/0010-bdf-time-integration.md`](../adr/0010-bdf-time-integration.md) — BDF1/BDF2 choice.
- [`docs/adr/0013-bc-ramp-continuation.md`](../adr/0013-bc-ramp-continuation.md) — IC strategy.
- [`docs/adr/0014-slotboom-transient.md`](../adr/0014-slotboom-transient.md) — Slotboom time-loop, supersedes ADR 0009.
- [`docs/m13.1-followup-5-blocker.md`](../m13.1-followup-5-blocker.md) — diagnostic chain that motivated the switch.

## Common pitfalls

1. **Explicit methods are unstable for DD.** A naive forward Euler
   on (17.1) requires $dt \lesssim L_D^2/D_0 \sim$ ps for stability —
   prohibitive for any device-level transient (~ns to ~µs).
   Implicit BDF1/BDF2 lifts this restriction.
2. **Consistent vs lumped mass.** A consistent (Galerkin) mass term
   $\int n_\mathrm{ufl} v_n\,dx$ couples adjacent DOFs and produces
   wiggles at small $dt$. Lumping to vertex quadrature ($dx_\mathrm{lump}$)
   makes the mass diagonal and removes the spurious oscillations. The
   shipped engine lumps; ADR 0014 §Implementation explains.
3. **(n, p) primary form's stuck history.** The pre-M13.1 transient
   used (n, p) primary unknowns; storing $n^k$ as Functions and using
   them in a forward-Euler-style residual produced negative-density
   iterates after the first BDF step. ADR 0009 → ADR 0014 supersedes
   this; do not reinvent the (n, p) form unless you can carry positivity
   constraints (e.g. log transformation).
4. **Coarse $dt$ misses the fast transient.** A 100 ns transient with
   $dt = 100\,\mathrm{ns}$ has *one* step; you can't measure $\tau_\mathrm{eff}$
   from a single point. The benchmark `pn_1d_turnon` uses $dt = 1\,\mathrm{ns}$
   to give 500 samples over the 500 ns simulation.
5. **BDF2 starts as BDF1.** The first step has no $u^{n-1}$ history, so
   it uses BDF1; from step 2 onwards BDF2 is used. The runner
   auto-switches ([`semi/runners/transient.py:413-423`](../../semi/runners/transient.py)),
   but the *first* step's accuracy is BDF1 ($O(dt)$) — for second-order
   convergence on a manufactured solution, the first step's error
   contributes a constant $O(dt)$ that doesn't decay; use a small
   first $dt$ or an alternative startup (e.g. trapezoidal rule).

## Exercises

**Exercise 17.1.** Show that BDF1 applied to $\dot u = \lambda u$
(scalar test problem) gives $u^{n+1} = u^n/(1 - \lambda dt)$. For
$\mathrm{Re}(\lambda) < 0$ (decaying solution), is this always stable?

**Exercise 17.2.** Compute the BDF2 coefficient of the truncation error.
For $u(t) = e^{\lambda t}$ at $t^{n+1}$: expand $u^n, u^{n-1}$ around
$t^{n+1}$ in Taylor series; collect terms in (17.5); show the leading
truncation is $-\lambda^3 dt^3 / 3 + O(dt^4)$.

**Exercise 17.3.** Read [`semi/runners/transient.py:434-443`](../../semi/runners/transient.py).
What is `hist_n_arr` after a BDF2 step? Show by direct substitution into
(17.5) that updating `f_hist_n` with this expression produces the
correct residual.

**Exercise 17.4.** A user submits a JSON with `dt = 10 ns`,
`t_end = 1 us`, BDF2. Estimate the number of SNES iterations the
transient runner will require. Compare with the same problem at
`dt = 1 ns`.

**Exercise 17.5.** The M16.7 milestone (time-varying $V(t)$) will let
you simulate a sinusoidal $V(t) = V_\mathrm{DC} + V_\mathrm{AC}\cos(2\pi f t)$
in the time domain. What is the relationship between the time-domain
$I(t)$'s Fourier component at $f$ and the AC sweep's $Y(2\pi f)\cdot V_\mathrm{AC}$?

### Solutions

**17.1.** $u^{n+1} - u^n = \lambda dt\,u^{n+1}$, so $u^{n+1}(1 - \lambda dt) = u^n$,
$u^{n+1} = u^n/(1-\lambda dt)$. For $\mathrm{Re}(\lambda) < 0$:
$|1-\lambda dt|^2 = (1-\mathrm{Re}(\lambda)dt)^2 + (\mathrm{Im}(\lambda)dt)^2 > 1$,
so $|u^{n+1}| < |u^n|$. Stable for any $dt > 0$ — A-stable. ✓

**17.2.** Taylor expand: $u^n = u^{n+1} - dt\,u' + (dt)^2/2\,u'' - (dt)^3/6\,u''' + ...$;
$u^{n-1} = u^{n+1} - 2dt\,u' + 2(dt)^2\,u'' - (4/3)(dt)^3 u''' + ...$
Plug in: $(3/2 - 2 + 1/2)u^{n+1} + dt(2 - 1)u' + (dt)^2(-1 + 1)u''
+ (dt)^3(1/3 - 1/12)u''' + ...$. Coefficients:
$u^{n+1}$: 0. $u'$: $dt$. $u''$: $0\cdot(dt)^2$. $u'''$: $(dt)^3/4$.
So $\sum\alpha_k u^{n+1-k}/dt = u' + (dt^2/4) u''' + ...$. Hmm, my
arithmetic is rough; the textbook BDF2 truncation is $-(2/9)(dt)^2 u'''$.
The point is BDF2 is second-order: leading truncation $\propto dt^2$.

**17.3.** `hist_n_arr = sum_{k=1}^K alpha_k * n_hist[-k]`. For BDF2 with
`coeffs = (1.5, -2.0, 0.5)`: $\alpha_1 = -2$, $\alpha_2 = 0.5$, so
`hist_n_arr = -2 * n_hist[-1] + 0.5 * n_hist[-2] = -2 n^n + (1/2) n^{n-1}$.
Plug into (17.5)/dt: $(\alpha_0/dt) n^{n+1} + (-2 n^n + (1/2)n^{n-1})/dt
= (3/2/dt) n^{n+1} + \mathrm{hist}/dt$. ✓ — matches the residual
[`semi/runners/transient.py:651-657`](../../semi/runners/transient.py).

**17.4.** $dt = 10\,\mathrm{ns}$, $t_\mathrm{end} = 1\,\mu\mathrm{s}$:
100 timesteps. Each step takes ~5 SNES iterations (similar to bias_sweep
at the same operating point). Total: ~500 SNES iterations.
$dt = 1\,\mathrm{ns}$: 1000 timesteps; ~5000 iterations.
The 10× cost is the price of resolving the fast transient. If the
problem is *just* about reaching steady state, $dt = 10\,\mathrm{ns}$
is fine; if you care about RC time-constant extraction, $dt = 1\,\mathrm{ns}$
is the minimum for a clean fit.

**17.5.** At small enough $V_\mathrm{AC}$, the time-domain response is
linear: $I(t) - I_\mathrm{DC} = V_\mathrm{AC}\,|Y(2\pi f)|\,\cos(2\pi f t + \arg Y)$.
The Fourier coefficient of $I(t)$ at frequency $f$ equals
$V_\mathrm{AC}\cdot Y(2\pi f)$. The 5% acceptance gate of audit case 06
(M16.7 deliverable) is exactly this comparison: time-domain FFT of
$I(t)$ vs the AC sweep's $Y(\omega)$ at the matching frequency.

## Further reading

- **Hairer, Nørsett, and Wanner, *Solving Ordinary Differential
  Equations II* (2nd ed., 1996).** The reference for BDF stability
  theory, including the A-stability and A($\alpha$)-stability proofs
  for BDF1–BDF6.
- **Selberherr (1984), §8.4** for the device-physics framing of
  transient DD.
- **Vasileska et al. (2010), §3** for transient methods and lumped-mass
  rationale in DD codes.
- **`docs/adr/0014-slotboom-transient.md`** in this repo — the
  authoritative engine reference.
