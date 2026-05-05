# 18 — Small-signal AC analysis

## Learning objectives

- Derive the small-signal AC system $(J + j\omega M)\delta\mathbf{u} =
  -\partial F/\partial V\,\delta V$ by linearizing the steady-state
  residual around a converged DC operating point.
- Explain the **real 2×2 block reformulation** that lets a real-PETSc
  build solve the complex linear system at doubled size.
- Recognize the displacement-current contribution $j\omega Q_\psi$ and
  the conduction-current contribution $\delta J_\mathrm{cond}$ in the
  terminal admittance.
- Distinguish the general AC-sweep runner (`run_ac_sweep`) from the
  analytic Poisson-sensitivity runner for MOSCAPs (`run_mos_cap_ac`),
  and identify when each applies.
- Locate the sign-convention errata in ADR 0011 and explain why the
  fix preserved bit-identity on the `rc_ac_sweep` benchmark.

## Physical motivation

Capacitance-voltage analysis, admittance spectroscopy, S-parameter
extraction at GHz, and noise modeling all reduce to "linearize the
device around its DC operating point and solve a complex frequency-domain
linear system." Compared with running a transient simulation and
computing an FFT, the AC small-signal solve is orders of magnitude
cheaper per frequency point: one linear solve per $\omega$ vs hundreds
or thousands of nonlinear time steps.

The shipped engine has two AC paths:

1. **`run_ac_sweep`** ([`semi/runners/ac_sweep.py`](../../semi/runners/ac_sweep.py)):
   the general two-terminal admittance $Y(\omega)$ at a chosen DC bias,
   for ohmic-only devices. Verified by the `rc_ac_sweep` benchmark
   to within 0.4% of analytical depletion-C over 1 Hz to 1 MHz.
2. **`run_mos_cap_ac`** ([`semi/runners/mos_cap_ac.py`](../../semi/runners/mos_cap_ac.py)):
   the MOSCAP-specific differential capacitance $dQ_\mathrm{gate}/dV_g$
   via Poisson sensitivity. Replaces the noisier `numpy.gradient(Q, V)`
   in `mos_cv` and is byte-identical on $Q_\mathrm{gate}$.

ADR 0011 is the load-bearing reference; this chapter walks the math.

## Derivation from first principles

### Linearization of the steady-state residual

Let $u_0$ be the converged steady-state Slotboom solution at DC bias
$V_\mathrm{DC}$: $F(u_0; V_\mathrm{DC}) = 0$. Under a small bias
perturbation $V(t) = V_\mathrm{DC} + \delta V\,e^{j\omega t}$, the
unknown moves to $u(t) = u_0 + \delta u\,e^{j\omega t}$. The transient
residual is

$$
F_\mathrm{trans}(u; V) = F(u; V) + M\,\frac{du}{dt},
\tag{18.1}
$$

with $M$ the mass matrix from Ch. 17. Linearize around $(u_0, V_\mathrm{DC})$:

$$
F_\mathrm{trans}(u_0 + \delta u; V_\mathrm{DC} + \delta V)
\approx \frac{\partial F}{\partial u}\delta u + \frac{\partial F}{\partial V}\delta V
+ M\,j\omega\,\delta u = 0.
$$

Defining $J = \partial F/\partial u$ (the steady-state Jacobian),
rearranging gives the **AC small-signal system**:

$$
\boxed{
(J + j\omega M)\,\delta u = -\frac{\partial F}{\partial V}\,\delta V.
}
\tag{18.2}
$$

This is the AC linearisation. Solving at each $\omega$ in a sweep
gives the small-signal response of the device.

### The right-hand side via finite-difference DC sensitivity

Computing $\partial F/\partial V$ analytically is messy because the bias
enters through the contact Dirichlet BC. The engine instead uses a
finite-difference trick. Solve the steady state at $V_\mathrm{DC}$
*and* at $V_\mathrm{DC} + \epsilon_V$:

$$
F(u_0; V_\mathrm{DC}) = 0,
\qquad F(u_0'; V_\mathrm{DC} + \epsilon_V) = 0.
$$

The DC sensitivity is

$$
\delta u_\mathrm{DC} \equiv \frac{u_0' - u_0}{\epsilon_V}
\approx \frac{\partial u}{\partial V}\bigg|_{V_\mathrm{DC}}.
\tag{18.3}
$$

Differentiate $F(u(V); V) = 0$: $J\delta u_\mathrm{DC} + \partial F/\partial V = 0$,
so $-\partial F/\partial V = J\delta u_\mathrm{DC}$. The AC RHS is
$b = J\delta u_\mathrm{DC}$ — one mat-vec, no analytical
differentiation. This is the construction at
[`semi/runners/ac_sweep.py:209-217`](../../semi/runners/ac_sweep.py)
with $\epsilon_V = 10^{-3}\,\mathrm{V}$.

### Real 2×2 block reformulation

PETSc supports both real and complex scalar builds; the dolfinx-real
build (the default in CI) uses real `PetscScalar`. Solving (18.2)
directly with complex matrices requires `PETSC_USE_COMPLEX=1`, which
the engine does not assume.

Standard workaround: write $\delta u = x + j y$ with $x, y$ real. The
complex equation $(J + j\omega M)(x + jy) = b_R + jb_I$ becomes the
**real 2×2 block system**:

$$
\begin{bmatrix} J & -\omega M \\ \omega M & J \end{bmatrix}
\begin{bmatrix} x \\ y \end{bmatrix}
= \begin{bmatrix} b_R \\ b_I \end{bmatrix}.
\tag{18.4}
$$

The block size is $2\cdot 3N$ where $N$ is the per-block DOF count. A
single direct LU (MUMPS) factorization handles this at each frequency.
For sub-1 kDOF 1D problems (the M14 benchmark sizes), this is sub-second
per frequency. For ~100 kDOF 2D / 3D problems, the doubled size
becomes painful and the M15 GPU path becomes attractive (Ch. 19).

Implementation: [`semi/runners/ac_sweep.py:518-591`](../../semi/runners/ac_sweep.py)
(`_solve_2x2_real_block`) builds the block matrix via scipy.sparse and
hands it to PETSc.

### Mass matrix in Slotboom form

The mass matrix on the continuity rows is, from (17.2):

$$
M_n = \int (n_\mathrm{ufl})\,v_n\,dx_\mathrm{lump}\;\text{evaluated as Jacobian wrt }(\psi, \Phi_n, \Phi_p).
\tag{18.5}
$$

UFL's `derivative` produces the chain-rule entries automatically:

$$
\frac{\partial(n_\mathrm{ufl} v_n)}{\partial\psi}\Big|_{\mathrm{lump}} = +n\,v_n\,\delta_{ij},
\qquad
\frac{\partial(n_\mathrm{ufl} v_n)}{\partial\Phi_n}\Big|_{\mathrm{lump}} = -n\,v_n\,\delta_{ij},
\qquad
\frac{\partial(n_\mathrm{ufl} v_n)}{\partial\Phi_p}\Big|_{\mathrm{lump}} = 0,
$$

with $\delta_{ij}$ the vertex-quadrature lumping kronecker. Same for the
hole row. The Poisson row has no time term, so the corresponding mass
block is zero. The full mass matrix is therefore a sparse 3×3 block with
entries on the (psi,phi_n), (phi_n,psi), (phi_n,phi_n), (psi,phi_p),
(phi_p,psi), (phi_p,phi_p) sub-blocks.

[`semi/runners/ac_sweep.py:326-388`](../../semi/runners/ac_sweep.py)
assembles $M$ via UFL `derivative` and standard `assemble_matrix`. BC
rows of $M$ are zeroed (`diag = 0.0`) so the BC perturbation is
imposed identically at every frequency.

### Terminal admittance

The total terminal current at the swept contact is

$$
I_\mathrm{total} = I_\mathrm{cond} + I_\mathrm{disp}
= \int_\Gamma \mathbf{J}_\mathrm{cond}\!\cdot\!\hat{\mathbf{n}}\,dS
+ \int_\Gamma \frac{\partial\mathbf{D}}{\partial t}\!\cdot\!\hat{\mathbf{n}}\,dS,
$$

with $\mathbf{D} = \varepsilon\mathbf{E} = -\varepsilon\nabla\psi$ the
displacement field. In the small-signal limit:

$$
\delta I_\mathrm{cond} = \int_\Gamma\frac{\partial\mathbf{J}_\mathrm{cond}}{\partial u}\,\delta u\,\hat{\mathbf{n}}\,dS,
\qquad
\delta I_\mathrm{disp} = -j\omega\int_\Gamma\varepsilon\nabla(\delta\psi)\!\cdot\!\hat{\mathbf{n}}\,dS.
\tag{18.6}
$$

The conduction part uses UFL `derivative` of the steady-state
terminal-current expression
$\mathbf{J}_n = -q\mu_n n\nabla\Phi_n$ (Ch. 11) at the converged $u_0$
([`semi/runners/ac_sweep.py:594-670`](../../semi/runners/ac_sweep.py)).
The displacement part is $-j\omega Q_\psi$ where $Q_\psi$ is the
gate-side charge integral
([`semi/runners/ac_sweep.py:708-757`](../../semi/runners/ac_sweep.py)).

The terminal admittance is

$$
Y(\omega) = \frac{\delta I_\mathrm{total}}{\delta V}.
\tag{18.7}
$$

The capacitance is $C(\omega) = +\mathrm{Im}(Y)/(2\pi f)$; the
conductance is $G(\omega) = \mathrm{Re}(Y)$.

### Sign convention (current INTO device)

ADR 0011 reports the original M14 implementation flipped sign:
$Y$ was reported with current OUT of device, while the prose claimed
current INTO. The fix made $Y$ match `bias_sweep`'s convention
(current INTO; $\mathrm{Re}(Y(\omega \to 0)) = dI/dV$). The fix
preserved bit-identity of $C$ on the `rc_ac_sweep` benchmark because
both $\mathrm{Im}(Y)$ and the formula's sign flipped together.
ADR 0011 §Errata documents the diagnosis.

### Slotboom vs (n,p)-primary AC

The original M14 implementation used (n, p)-primary form for AC. Audit
case 05 (forward bias, $V_\mathrm{DC} = +0.4\,\mathrm{V}$) showed a
~12% disagreement vs `bias_sweep`'s centred-difference $dI/dV$. Errata
#2 in ADR 0011 traced this to a *discrete-residual mismatch*: the
Slotboom-converged operating point had a non-zero (n,p)-form discrete
residual, so $J_{(n,p)}\delta u_\mathrm{DC}$ was systematically biased
away from the linearisation.

The fix was to rewrite `ac_sweep` in **Slotboom variables** throughout:
the Jacobian is now exactly the discrete `bias_sweep` Jacobian; the
mass matrix is the chain-rule Slotboom mass; the terminal-current form
is the `ufl.derivative` of the same Slotboom-form expression
`bias_sweep` uses. By Newton's lemma at the bias-sweep operating point,
$J_\mathrm{Slotboom}\delta u_\mathrm{DC} = -\partial F_\mathrm{Slotboom}/\partial V$
is exact with no conversion step. The cross-runner consistency is
restored.

### MOSCAP differential capacitance via Poisson sensitivity

For a MOSCAP — multi-region, gate-Dirichlet, equilibrium-Poisson-only —
the AC sweep machinery doesn't directly apply (it targets ohmic
two-terminal pn devices). The dedicated `mos_cap_ac` runner
([`semi/runners/mos_cap_ac.py`](../../semi/runners/mos_cap_ac.py))
implements the $\omega \to 0$ limit:

1. At each $V_g$, solve the equilibrium Poisson and build the Jacobian
   $K = \partial F/\partial\psi$ at the converged $\psi_0$.
2. The gate's Dirichlet condition shifts by $\delta V_g$ uniformly:
   $\partial F/\partial V_g$ is concentrated at the gate row.
3. Solve $K\,\delta\psi = -\partial F/\partial V_g$ with $\delta\psi_\mathrm{gate}
   = 1/V_t$ (per-unit-V_g BC perturbation) and homogeneous BCs elsewhere.
4. Compute the differential gate charge:

$$
\frac{dQ_\mathrm{gate}}{dV_g}
= -\frac{q}{W_\mathrm{lat}}\int_{\Omega_\mathrm{Si}}\frac{\partial\rho}{\partial\psi}\,\delta\psi\,dA
= +\frac{q}{W_\mathrm{lat}}\int_{\Omega_\mathrm{Si}} n_i\bigl(e^{-\psi_0/V_t} + e^{\psi_0/V_t}\bigr)\,\delta\psi\,dA.
\tag{18.8}
$$

(Differentiating $\rho = n_i(e^{-\psi/V_t} - e^{\psi/V_t}) + N$ wrt $\psi$
gives $\partial\rho/\partial\psi = -(n_i/V_t)(e^{-\psi/V_t} + e^{\psi/V_t})$.
Including the chain rule for $\delta\psi$ and per-area normalization
gives (18.8).)

The numerator is computed by `sensitivity_form` at
[`semi/runners/mos_cap_ac.py:158-165`](../../semi/runners/mos_cap_ac.py).
This avoids the noise of `numpy.gradient(Q, V)`. The audit case 03
verifies byte-identity on $Q_\mathrm{gate}$ (the integrand is the same;
only the differentiation pathway differs).

## Key results

- AC linearisation: (18.2).
- DC sensitivity via FD: (18.3).
- Real 2×2 block: (18.4).
- Slotboom mass: (18.5).
- Terminal admittance: (18.7).
- MOSCAP differential C via Poisson sensitivity: (18.8).
- Sign convention: $Y$ is current INTO device; $C = +\mathrm{Im}(Y)/(2\pi f)$.

## Worked numerical example

**`rc_ac_sweep` benchmark.** 1D pn diode at $V_\mathrm{DC} = -1\,\mathrm{V}$,
41 logspace frequencies from 1 Hz to 1 GHz. Analytical depletion
capacitance per area: $C_\mathrm{dep} = \varepsilon_s/W(V) = 1.036\times 10^{-10}/W$.
At $V = -1\,\mathrm{V}$, $W = 218\,\mathrm{nm}$ (Ch. 7 worked example),
so $C_\mathrm{dep} = 4.75\times 10^{-4}\,\mathrm{F/m^2}$.

The benchmark verifier asserts the simulated $C$ over 1 Hz to 1 MHz
matches this value within 5%. The actual match is **0.41%** worst-case
([`docs/adr/0011-ac-small-signal.md`](../adr/0011-ac-small-signal.md)
Validation section). The plateau flatness $\max(C)/\min(C) = 1.000$
across 27 in-band samples — the depletion capacitance is genuinely
frequency-independent in this regime.

At $f \to 1\,\mathrm{GHz}$, displacement current and minority-carrier
response are no longer negligible; the curve drops as the carriers
stop following.

**MOSCAP differential capacitance.** At $V_g = -0.7\,\mathrm{V}$ on
the M14.2 axisymmetric MOSCAP (Hu Fig. 5-18 reproduction): the engine
solves the sensitivity (18.8) and reports $C_\mathrm{ac}$ in $\mathrm{F/m^2}$.
Comparison with the LF analytical curve from `lf_cv_quasistatic` on the
same parameters shows agreement to a few percent through the
depletion regime. The C–V curve as a whole reproduces Hu Fig. 5-18 to
visual accuracy on the notebook 05 plot.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| `run_ac_sweep` | (18.2)+(18.4)+(18.7) | `semi/runners/ac_sweep.py:151-499` |
| DC sensitivity FD | (18.3) | `semi/runners/ac_sweep.py:209-217` |
| Steady-state Jacobian assembly | $J$ | `semi/runners/ac_sweep.py:285-324` |
| Mass matrix assembly | (18.5) | `semi/runners/ac_sweep.py:326-388` |
| Real 2×2 block solve | (18.4) | `semi/runners/ac_sweep.py:518-591` (`_solve_2x2_real_block`) |
| Linearised terminal current | (18.6) | `semi/runners/ac_sweep.py:594-705` |
| Displacement charge | $Q_\psi$ | `semi/runners/ac_sweep.py:708-757` |
| `AcSweepResult` | – | `semi/results.py:56-114` |
| `run_mos_cap_ac` | (18.8) | `semi/runners/mos_cap_ac.py:59-329` |
| Sensitivity form (MOSCAP) | (18.8) integrand | `semi/runners/mos_cap_ac.py:158-165` |
| `LinearProblem` for sensitivity | – | `semi/runners/mos_cap_ac.py:260-268` |
| `rc_ac_sweep` benchmark | – | `benchmarks/rc_ac_sweep/` |

## Existing-docs cross-reference

- [`docs/adr/0011-ac-small-signal.md`](../adr/0011-ac-small-signal.md) — full ADR with Decision, Validation, Errata sections.
- [`docs/theory/moscap_cv.md`](../theory/moscap_cv.md) — LF/HF discussion and the depletion-clamp HF C–V.
- [`docs/IMPROVEMENT_GUIDE.md` §M14, §M14.1](../IMPROVEMENT_GUIDE.md) — original deliverable.
- [`docs/PHYSICS_AUDIT.md`](../PHYSICS_AUDIT.md) — audit cases 02, 03, 05 referenced in the Errata.

## Common pitfalls

1. **Sign conventions are easy to get wrong.** ADR 0011 has *two*
   errata: the original sign convention (current INTO vs OUT) and the
   primary-unknown choice (Slotboom vs (n,p)-primary). Both are subtle
   mistakes that pass coarse plausibility checks (the magnitudes look
   right; only the cross-runner audit catches them). When adding a
   new AC code path, replicate the audit-case structure: compare to
   `bias_sweep` $dI/dV$ at the same operating point.
2. **Real 2×2 block doubles the linear-system size.** A 100k DOF
   problem becomes a 200k×200k system per frequency. Acceptable for
   sub-10k benchmarks; pricey for 3D problems. M15 GPU AMG handles the
   resulting scale; complex-PETSc would be cleaner if available
   ([`docs/gpu.md`](../gpu.md) lists this as deferred).
3. **Mass matrix at BC rows.** Setting `diag = 0.0` on the BC rows of
   $M$ (vs `diag = 1.0` on $J$) ensures the BC perturbation $\delta u_\mathrm{BC}
   = \delta V/V_t$ is imposed at every frequency identically. If you
   set `diag = 1.0` on $M$ too, you get spurious $j\omega$ contributions
   at BC DOFs.
4. **Finite-difference $\epsilon_V$ tradeoff.** Too small ($10^{-7}$):
   FD noise dominates the DC sensitivity; the AC RHS is corrupted.
   Too large ($10^{-1}$): linearization error dominates. The default
   $10^{-3}\,\mathrm{V}$ at [`semi/runners/ac_sweep.py:209`](../../semi/runners/ac_sweep.py)
   is a careful choice well below $V_t \sim 25\,\mathrm{mV}$ but above
   the SNES tolerance floor.
5. **MOSCAP sensitivity $K$ is not the same as bias_sweep's $J$.**
   The MOSCAP `mos_cap_ac` runner solves the equilibrium-Poisson
   sensitivity, not the full DD Jacobian. The Jacobian is constant in
   $V_g$ (only the BC depends), so this is a legitimate simplification —
   but if you wanted the full Slotboom AC for a MOSCAP, you would need
   to extend `run_ac_sweep` to handle the gate-Dirichlet path, which
   ADR 0011 explicitly defers.

## Exercises

**Exercise 18.1.** Show that for a pure capacitor ($J = 0$, $M = C\cdot I$),
$Y(\omega) = j\omega C$. What does $C(\omega)$ from $\mathrm{Im}(Y)/(2\pi f)$
return?

**Exercise 18.2.** Reproduce the analytical depletion capacitance per
area at $V = -1\,\mathrm{V}$ for the `rc_ac_sweep` device (M3 parameters):
$N_A = N_D = 10^{17}\,\mathrm{cm^{-3}}$, $\varepsilon_r = 11.7$.

**Exercise 18.3.** Why does the M14 implementation prefer the
finite-difference DC sensitivity over computing $\partial F/\partial V$
analytically? List two practical advantages.

**Exercise 18.4.** Read the mass-matrix assembly at
[`semi/runners/ac_sweep.py:326-388`](../../semi/runners/ac_sweep.py).
Why does the Poisson row's mass form start with `zero_psi_c * v_psi * ufl.dx`?

**Exercise 18.5.** A user wants to extract $Y(\omega)$ for a MOSFET at
GHz frequencies. Explain why `run_ac_sweep` is the right tool and
`run_mos_cap_ac` is not.

### Solutions

**18.1.** Pure capacitor: $J = 0$, $M = C I$. (18.2): $j\omega C\delta u
= -\partial F/\partial V \delta V$. With $-\partial F/\partial V = 1$
(unit BC perturbation at the contact), $\delta u = 1/(j\omega C)\delta V$.
$\delta I = j\omega C\delta u\cdot ?$... Actually, for a parallel-plate
capacitor in this idealization, $\delta I = j\omega Q' = j\omega C \delta V$,
so $Y = \delta I/\delta V = j\omega C$. $C(\omega) = \mathrm{Im}(Y)/(2\pi f)
= \omega C/(2\pi f) = C$. ✓ (independent of $\omega$, as expected for
an ideal capacitor).

**18.2.** $V_{bi} = 0.834\,\mathrm{V}$, $W(-1) = W(0)\sqrt{(0.834+1)/0.834}
= 146.9\,\mathrm{nm}\cdot\sqrt{2.20} = 218\,\mathrm{nm}$.
$C_\mathrm{dep} = \varepsilon_s/W = 11.7\cdot 8.854\times 10^{-12}/2.18\times 10^{-7}
= 4.75\times 10^{-4}\,\mathrm{F/m^2}$. Matches the ADR 0011 acceptance
quote.

**18.3.** (a) Analytical $\partial F/\partial V$ requires inspecting
which Dirichlet BC depends on $V$ and computing its derivative — error-prone
across multi-contact, multi-runner-shared-BC code. The FD trick reuses
the *same* BC-construction code path that `bias_sweep` uses. (b) The FD
trick goes through `bias_sweep` at $V_\mathrm{DC}$ and $V_\mathrm{DC} + \epsilon$,
which exercises the full discretized residual; it implicitly captures
any subtle BC handling that an analytical $\partial F/\partial V$
would have to re-derive.

**18.4.** UFL refuses to compile a bare `0 * v_psi * dx` form (the
type system needs a function-valued integrand). Multiplying by a zero
`Constant` produces a well-formed form with a zero contribution to the
matrix. The Poisson row's mass block is thereby zero, as it should be.

**18.5.** `run_mos_cap_ac` only handles the equilibrium-Poisson
sensitivity, not the time-derivative term — it returns the $\omega \to 0$
differential capacitance only. A MOSFET at GHz needs the full
$(J + j\omega M)$ system with the carrier mass matrix; that's
`run_ac_sweep`. (The `run_ac_sweep` runner targets ohmic two-terminal
devices; a MOSFET has four terminals plus a gate, so the runner would
also need extension to handle multi-terminal Y-parameter extraction —
deferred work.)

## Further reading

- **Selberherr (1984), Chapter 7 §"Small-signal AC analysis."**
- **Snowden, *Semiconductor Device Modelling* (1989).** §"Small-signal
  extraction from drift-diffusion."
- **Schenk, *Advanced Physical Models for Silicon Device Simulation*
  (1998).** Chapter 5 §"Linearised continuity equations for AC analysis."
- **Tsividis (1999), Chapter 4.** Linearization of the MOSFET
  small-signal admittance.
- **`docs/adr/0011-ac-small-signal.md`** in this repo — load-bearing.
