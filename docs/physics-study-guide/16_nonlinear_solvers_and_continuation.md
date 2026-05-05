# 16 — Nonlinear solvers and continuation

## Learning objectives

- Apply Newton's method to a vector residual $F(\mathbf{u}) = \mathbf{0}$.
- Explain why Newton diverges from a cold start at high bias and how
  backtracking line search rescues some — but not all — failure modes.
- Recognize the PETSc SNES interface that dolfinx 0.10's
  `NonlinearProblem` wraps, and understand the role of MUMPS LU as the
  per-step linear solver.
- Identify the bias-continuation strategy as a homotopy from the
  equilibrium solution to the high-bias solution.
- Read the `AdaptiveStepController` state machine in
  `semi/continuation.py` and trace its grow/halve/clamp logic on a
  failed bias step.

## Physical motivation

The drift-diffusion residual $F(\mathbf{u})$ is sharply nonlinear
because $\mathbf{u}$ enters through exponentials in the Slotboom
expressions for $n$ and $p$. A single Newton step linearizes around the
current iterate; if the iterate is far from the root, the linearization
is a bad model and the resulting step overshoots into a nonphysical
state, after which subsequent residual evaluations may overflow,
underflow, or NaN. Bias continuation walks the solution from a
known-good equilibrium starting point through a sequence of small bias
increments, each well-conditioned, to the target bias. The shipped
engine combines Newton (with line search) inside SNES with an outer
adaptive-bias loop that grows and halves the step based on Newton
convergence behaviour. ADR 0008 documents the SNES tolerance defaults;
ADR 0011 documents the use of SNES for the AC linearisation; ADR 0014
documents the Slotboom-transient time-loop choice and its Jacobian
shift regularisation.

## Derivation from first principles

### Newton's method

For $F: \mathbb{R}^N \to \mathbb{R}^N$ with $F(\mathbf{u}^*) = 0$,
Newton iterates

$$
\mathbf{u}^{(k+1)} = \mathbf{u}^{(k)} - J(\mathbf{u}^{(k)})^{-1}\,F(\mathbf{u}^{(k)}),
\tag{16.1}
$$

with $J = \partial F/\partial\mathbf{u}$ the Jacobian. The local
convergence is quadratic ($\|\mathbf{u}^{(k+1)} - \mathbf{u}^*\| \lesssim C\|\mathbf{u}^{(k)} - \mathbf{u}^*\|^2$)
when $\mathbf{u}^{(0)}$ is close enough to $\mathbf{u}^*$. The radius of
quadratic convergence depends on the Lipschitz constant of $J$ — for
DD with exponential nonlinearities, this radius is small.

### Backtracking line search

When the full Newton step $\delta\mathbf{u} = -J^{-1}F$ produces a
worse residual, *backtracking* line search reduces the step by a factor
$\alpha < 1$ and re-evaluates:

$$
\mathbf{u}^{(k+1)} = \mathbf{u}^{(k)} + \alpha\,\delta\mathbf{u},
\quad \alpha\in\{1, 1/2, 1/4, ...\}.
$$

The step is accepted as soon as $\|F(\mathbf{u}^{(k+1)})\|$ is reduced
by a sufficient fraction (the **Armijo condition** with a default
factor $\sim 10^{-4}$). PETSc's `snes_linesearch_type bt` implements
this; kronos-semi uses it as the default
([`semi/solver.py:21`](../../semi/solver.py)).

Line search rescues moderate overshoots but not all failure modes. If
the Newton step pushes $\mathbf{u}$ into a region where Slotboom
densities overflow (Ch. 11 §pitfalls), even a very small $\alpha$ may
not bring the residual back. The cure is to ensure the *initial*
iterate is close enough — that is what bias continuation does.

### PETSc SNES wraps Newton + line search

`PETSc.SNES` (Scalable Nonlinear Equations Solver) provides Newton with
several line-search options (`bt`, `cp`, `l2`, `nleqerr`) and a
trust-region variant (`tr`). dolfinx 0.10's `NonlinearProblem` wraps a
SNES instance and exposes its options:

```python
NonlinearProblem(F, u, bcs=bcs,
                 petsc_options_prefix=prefix,
                 petsc_options={"snes_type": "newtonls",
                                "snes_linesearch_type": "bt", ...})
```

Default options ([`semi/solver.py:20-30`](../../semi/solver.py)):
- `snes_type: newtonls` — Newton with line search.
- `snes_linesearch_type: bt` — backtracking.
- `snes_rtol: 1e-10`, `snes_atol: 1e-12`, `snes_max_it: 50` — convergence
  criteria.
- `ksp_type: preonly`, `pc_type: lu`, `pc_factor_mat_solver_type: mumps`
  — per-Newton-step linear solve via MUMPS direct LU.

ADR 0008 explains why these defaults; the `bias_sweep` runner adopts
slightly different SNES tolerances (`atol = 1e-7`) for the DD block solve
because the absolute-tolerance scale differs across blocks
([`semi/runners/bias_sweep.py:106-111`](../../semi/runners/bias_sweep.py)).

### Why dolfinx 0.10 over earlier

dolfinx's pre-0.10 `NewtonSolver` did its own quasi-Newton iteration
without exposing PETSc SNES; that meant no built-in line search and no
direct way to query convergence diagnostics. dolfinx 0.10's
`NonlinearProblem` wraps SNES directly. ADR 0003 records this
decision; [`docs/theory/dolfinx_choice.md`](../theory/dolfinx_choice.md)
summarizes.

### Bias continuation as homotopy

Define the parametric residual $F(\mathbf{u}, V)$ at a sequence of
biases $0 = V_0 < V_1 < ... < V_K$. Solve at $V_0$ (typically equilibrium,
which has a known-good initial guess from charge neutrality). Use the
converged $\mathbf{u}^{(0)}$ as the *initial guess* for the solve at
$V_1$. Repeat. This is **homotopy continuation** in $V$.

For small enough $|V_{k+1} - V_k|$, the converged $\mathbf{u}_k$ is
inside the radius of quadratic convergence at $V_{k+1}$. The bias step
size is therefore the key knob: too small wastes work; too large
divurges Newton.

### Adaptive step controller

[`semi/continuation.py:32-133`](../../semi/continuation.py) implements
the engine's adaptive controller. State:
- `step`: signed step size (sign = sweep direction).
- `easy_count`: number of consecutive easy SNES solves
  (iterations < `easy_iter_threshold`, default 4).

On success with $n_\mathrm{iter} < $ threshold:
- Increment `easy_count`.
- If `easy_count >= threshold`: multiply `|step|` by `grow_factor`
  (default 1.5), clamp to `max_step_abs`, reset counter.

On success with $n_\mathrm{iter} \ge$ threshold:
- Reset `easy_count`.

On failure:
- Halve `|step|`. If the halved value is below `min_step_abs`, raise
  `StepTooSmall` (the controller's signal that the regime is outside
  the model's reach).

`clamp_to_endpoint(remaining)` rounds the next step so the new bias
exactly lands on the endpoint, avoiding overshoot.

This is a textbook **arc-length-style adaptive continuation** without
the curvilinear arc-length (which would cross fold points; the engine
doesn't need that yet).

### Bipolar sweep

When the requested sweep has both negative and positive entries, the
sweep is split into two single-direction legs:
$V = 0 \to \min(V) \to \max(V)$, each with a fresh
`AdaptiveStepController` ([`semi/runners/bias_sweep.py:347-361`](../../semi/runners/bias_sweep.py)).
This avoids re-traversing the equilibrium point twice and keeps the
controller state coherent on each leg.

### Jacobian shift for transient stability

The transient runner uses a small Jacobian shift
$\epsilon I$ added to the assembled block Jacobian
([`semi/solver.py:121-149`](../../semi/solver.py), `_install_jacobian_shift`).
This is a regularization for *rank-deficient minority-side rows* in the
Slotboom continuity equation: in the deep p-bulk, the electron density
$n = n_i\exp(\psi - \Phi_n)$ is below floating-point precision, so the
$(\Phi_n)$ row of the Jacobian is numerically zero at those vertices.
MUMPS reports null pivots and the Newton update is unbounded. Adding
$\epsilon I$ with $\epsilon \sim 10^{-14}$ removes the null pivots
without measurably perturbing the converged solution. ADR 0014
documents this; the steady-state `bias_sweep` runner does *not* need
the shift because its bias-continuation walk keeps every iterate near
a converged neighbour.

## Key results

- Newton iteration: (16.1).
- Backtracking line search: $\alpha\in\{1, 1/2, 1/4, ...\}$ chosen to
  satisfy Armijo.
- SNES default options: see [`semi/solver.py:20-30`](../../semi/solver.py).
- Adaptive controller: grow on easy, halve on fail, clamp to endpoint.
- Bipolar sweep: two single-direction legs with fresh controllers.
- Jacobian shift: $\epsilon\sim 10^{-14}$ for transient.

## Worked numerical example

**M2 forward bias from $V = 0$ to $V = 0.6\,\mathrm{V}$** with default
controller settings ($V_{step} = 0.05\,\mathrm{V}$, threshold = 4,
grow = 1.5, max = 0.05 V, min = $10^{-4}\,\mathrm{V}$):

- Step 1: solve at $V = 0.05$. Newton converges in 3 iterations
  (easy). `easy_count = 1`.
- Steps 2, 3: similarly, 3-4 iter each. `easy_count = 4` triggers
  growth to $0.075\,\mathrm{V}$ — but `max_step_abs = 0.05`, so it
  clamps to 0.05.
- Result: a uniform 0.05 V step throughout. Total: 12 steps
  $\times$ 4 iter $\approx$ 48 SNES iterations.

If `max_step_abs` were $0.1\,\mathrm{V}$, growth would land at 0.075 V
after step 4; if iterations stay easy, step 5 grows to 0.075 V,
step 6 grows to 0.1 V (clamped). Sweep finishes faster (~9 steps).

**With harder physics (high doping, deep reverse bias):**

- Step 1: solve at $V = -0.05$. SNES diverges after 50 iterations
  (line-search fails to find an improvement).
- Halve `step` to $0.025\,\mathrm{V}$. Restore the previous solution.
  Try $V = -0.025\,\mathrm{V}$.
- This time SNES converges in 8 iter (not "easy" but converged).
  `easy_count = 0`.
- Continue at $0.025\,\mathrm{V}$ steps until `easy_count = 4` triggers
  growth.
- If a step at $0.05\,\mathrm{V}$ later fails, halve to 0.025; if that
  fails, halve to 0.0125; ...; if 0.05 V fails 6 times in a row, the
  controller raises `RuntimeError` (`max_halvings = 6`).

## Code map

| Concept | Equation | Code location |
|---|---|---|
| Newton iteration | (16.1) | inside `dolfinx.fem.petsc.NonlinearProblem` (PETSc SNES) |
| `solve_nonlinear` wrapper | – | `semi/solver.py:186-249` |
| `solve_nonlinear_block` | – | `semi/solver.py:252-343` |
| Default PETSc options | – | `semi/solver.py:20-30` (`DEFAULT_PETSC_OPTIONS`) |
| Backend resolution | – | `semi/solver.py:152-183` (`_resolve_backend_options`) |
| Jacobian shift install | – | `semi/solver.py:121-149` (`_install_jacobian_shift`) |
| MUMPS factor options | – | `semi/solver.py:74-118` (`_apply_factor_options`) |
| `AdaptiveStepController` | – | `semi/continuation.py:32-133` |
| `on_success`, `on_failure`, `clamp_to_endpoint` | – | `semi/continuation.py:82-124` |
| Bias-sweep ramp loop | – | `semi/runners/bias_sweep.py:248-326` |
| Bipolar legs | – | `semi/runners/bias_sweep.py:347-361` |
| BC-ramp continuation IC | – | `semi/runners/transient.py:671-789` (transient-only) |

## Existing-docs cross-reference

- [`docs/theory/dolfinx_choice.md`](../theory/dolfinx_choice.md) — why dolfinx 0.10's `NonlinearProblem`.
- [`docs/PHYSICS.md` §2.4](../PHYSICS.md) — bias-continuation strategy.
- [`docs/adr/0003-dolfinx-0-10-api.md`](../adr/0003-dolfinx-0-10-api.md).
- [`docs/adr/0008-snes-tolerances.md`](../adr/0008-snes-tolerances.md).
- [`docs/adr/0013-bc-ramp-continuation.md`](../adr/0013-bc-ramp-continuation.md) (transient IC ramp; not separately covered here).
- [`docs/adr/0014-slotboom-transient.md`](../adr/0014-slotboom-transient.md) — Jacobian shift, transient SNES options.

## Common pitfalls

1. **Cold-start bias.** Starting Newton at $V = 0.6\,\mathrm{V}$ from
   the equilibrium initial guess fails because the Slotboom carrier
   ratio $\hat n/\hat n^\mathrm{eq}$ is $e^{0.6/V_t} \sim 10^{10}$.
   The bias-continuation walk is *not* optional — skipping it produces
   `SNES_DIVERGED_LINE_SEARCH` with cryptic diagnostics.
2. **Tolerances too tight.** Setting `snes_atol: 1e-14` for the DD
   block (instead of $10^{-7}$) makes SNES report "not converged" even
   when the residual is at the discretization's noise floor. The
   `bias_sweep` default of $10^{-7}$ is calibrated for this; ADR 0008
   explains.
3. **MUMPS workspace.** The transient runner bumps
   `mat_mumps_icntl_14: 200` (200% extra workspace) to handle the
   small extra delayed pivots from the Jacobian shift. Without this
   bump, MUMPS errors out with `INFO(1)=-9` ([`semi/runners/transient.py:159-164`](../../semi/runners/transient.py)).
4. **`prefix` collision.** PETSc options are stored in a global
   options database keyed by prefix. Reusing a prefix across solves
   leaks state. Each `solve_nonlinear[_block]` call must use a unique
   prefix. The runners build the prefix from the config name and a
   per-step tag (`fmt_tag(V)` for the bias step) to ensure uniqueness.
5. **Jacobian shift is *not* a residual fix.** Adding $\epsilon I$
   regularizes ill-conditioning, not a wrong residual. If the shift
   makes a benchmark "converge" but the answer is wrong, the residual
   is the bug, not the conditioning.

## Exercises

**Exercise 16.1.** Why does Newton converge in $\sim 4$ iterations on
the M1 equilibrium Poisson but in 8–12 iterations on the M2 forward-bias
DD block? (Hint: the residual nonlinearity grows.)

**Exercise 16.2.** A user submits a JSON with `voltage_sweep: {start: 0,
stop: 1.0, step: 1.0}` (a single 1 V step). What happens with default
continuation settings? Why is this different from `{start: 0, stop: 1.0, step: 0.05}`?

**Exercise 16.3.** Read the bipolar-leg logic in
[`semi/runners/bias_sweep.py:347-361`](../../semi/runners/bias_sweep.py).
Why does each leg get a fresh `AdaptiveStepController`?

**Exercise 16.4.** Compute the Newton iteration count for a 12-step
forward sweep with `easy_iter_threshold = 4`, `grow_factor = 1.5`,
each step converging in 3 iterations. How many total SNES iterations?

**Exercise 16.5.** A Schottky-contact runner (M16.5, planned) will
introduce a Robin-style boundary form on the carriers (Ch. 8). Predict
how this affects the radius of Newton convergence — better, worse, or
unchanged compared with the ohmic case?

### Solutions

**16.1.** Equilibrium Poisson is a single-block scalar PDE with a soft
exponential nonlinearity. The DD block is three coupled equations with
exponential carrier expressions on every term. Each Newton step has to
balance the cross-block coupling, which doubles or triples the
iteration count. Adding bias makes the residual sharper: at $V = 0.6$,
the minority carriers are 10 orders of magnitude away from their
equilibrium values, and the iterate must traverse that range smoothly.

**16.2.** A single 1 V step has the equilibrium iterate as initial
guess, which is far outside the radius of quadratic convergence at
$V = 1\,\mathrm{V}$. SNES diverges; the controller halves the step to
0.5 V. If that also diverges, halve to 0.25 V, and so on. After 6
halvings (default `max_halvings`) the controller raises `RuntimeError`.
With `step: 0.05`, the controller starts with a 0.05 V step that
typically converges in 3-4 iterations from the start, and grows
modestly. Total work: ~12 steps × 3-4 iter = 36-48 iter; vs the
single-step path which fails.

**16.3.** State (`easy_count`, current `step`) accumulated on the
forward leg may not be appropriate for the reverse leg: a forward
sweep that just grew its step size to 0.1 V might immediately fail on
the first reverse-direction step because the physics differs. Fresh
state lets each leg adapt independently. The *initial* `step` for each
leg starts at the user-specified `step` size and signs by direction.

**16.4.** Step 1: 3 iter, `easy = 1`. Step 2: 3 iter, `easy = 2`.
Step 3: 3 iter, `easy = 3`. Step 4: 3 iter, `easy = 4` triggers growth
to 0.075 V (or to `max_step_abs` if smaller). Reset `easy = 0`.
Step 5: 3 iter (still easy), `easy = 1`. ... With max_step clamped at
0.05, the growth is a no-op. Total: 12 × 3 = 36 SNES iter. With
max_step = 0.1 V allowing growth, fewer total steps but possibly more
iter per step. The exact tradeoff depends on the residual nonlinearity.

**16.5.** A Robin BC on the carriers turns a Dirichlet (rank-1 BC)
into a flux relation that couples the carrier value to the residual at
the contact. The Jacobian gains an extra block at the contact rows but
remains well-conditioned (the Robin term is positive-definite). Newton
convergence radius is roughly *unchanged* for moderate barrier heights;
for very large barriers ($\phi_B \gg V$), the saturation current is
exponentially small and the contact rows are nearly decoupled — at
which point the regime is essentially reverse-biased Schottky and the
solver may need its own continuation in $\phi_B$ from a small reference
to the target. M16.5's acceptance benchmark will surface any such
issues.

## Further reading

- **Dennis and Schnabel, *Numerical Methods for Unconstrained
  Optimization and Nonlinear Equations* (1996).** Newton + line-search
  fundamentals.
- **Nocedal and Wright, *Numerical Optimization*, 2nd ed. (2006).**
  Chapter 11 covers Newton's method for nonlinear systems with
  algorithmic detail, including Armijo and Wolfe conditions.
- **Allgower and Georg, *Numerical Continuation Methods* (2003).**
  The reference for continuation methods, including arc-length and
  pseudo-arc-length variants.
- **PETSc SNES manual:** https://petsc.org/release/manualpages/SNES/
  for the runtime options and internals.
- **Selberherr (1984), Chapter 8** for continuation strategies in DD
  device simulators.
- **`docs/theory/dolfinx_choice.md`** in this repo.
