# 21 — Verification and benchmarks

## Learning objectives

- State the difference between **verification** ("did we solve the
  equations correctly?") and **validation** ("are the equations the
  right model?") in the Roache / Oberkampf-Roy sense.
- Sketch the Method of Manufactured Solutions (MMS) and explain why
  observed L² rates of 2.0 confirm the discrete operator's
  consistency.
- Identify the four conservation gates the V&V suite uses: charge
  conservation in equilibrium; total-current continuity in
  forward/reverse bias.
- Enumerate the eight shipped benchmarks and the analytical reference
  each verifier compares against.
- Reproduce every numerical assertion in
  `tests/check_analytical_math.py` from formulas earlier in the guide.

## Physical motivation

A device simulator that converges to a wrong answer is worse than no
simulator at all. The kronos-semi V&V suite was the M4 deliverable
(ADR 0006); it is the project's load-bearing answer to "how do you know
the engine is right?"

The suite lives under [`semi/verification/`](../../semi/verification/)
and runs in CI on every push via `scripts/run_verification.py all`
inside the `docker-fem` job. Acceptance is binary: every gate must be
green or the build fails. As of v0.16.0 the suite reports 62/62 PASS;
finest-pair MMS rates clear $L^2 \ge 1.99$ on every gated study.

This chapter is your roadmap to the V&V suite. Unlike the physics
chapters, the goal here is to *trust the engine* — to convince yourself
that the numbers it produces are correct.

## Verification vs validation

- **Verification** is a property of the *code*: did we solve the
  equations we wrote down correctly? MMS, mesh convergence, and
  conservation gates are verification activities. Failure of a
  verification gate is a code bug.
- **Validation** is a property of the *model*: are the equations a
  reasonable description of the physical device? Comparison against
  Shockley diode, Pao-Sah MOSFET, depletion-approximation MOSCAP are
  validation activities. Failure of a validation gate often means the
  model is missing physics (e.g. the M16.x catalogue), not that the
  code is wrong.

ADR 0006 commits to keeping these activities visibly separate. The
existing `pn_1d`, `pn_1d_bias`, `pn_1d_bias_reverse` benchmarks are
*validation* exercises (compare the engine to depletion-approximation
analytics); the new V&V suite is *verification* (MMS, mesh convergence,
conservation).

## The four V&V activities (ADR 0006)

### Phase 1 — MMS for equilibrium Poisson

[`semi/verification/mms_poisson.py`](../../semi/verification/mms_poisson.py).
Manufactured smooth solutions ($\sin(kx)$ in 1D, $\sin(k_xx)\sin(k_yy)$
in 2D) injected with UFL-weak forcing; mesh refined; observed L² and
H¹ rates checked.

Acceptance: finest-pair $L^2 \ge 1.85$, $H^1 \ge 0.85$, monotone error
reduction. Latest finest-pair rates:

| Study | $L^2$ rate | $H^1$ rate |
|---|---|---|
| 1D linear | 2.000 | 1.000 |
| 1D nonlinear | 2.000 | 1.000 |
| 2D triangles | 1.998 | 0.999 |

Theoretical P1 Lagrange rates are $L^2 = 2$, $H^1 = 1$ (Brenner-Scott).
Observed = theoretical to roundoff.

### Phase 2 — Mesh convergence on `pn_1d`

[`semi/verification/mesh_convergence.py`](../../semi/verification/mesh_convergence.py).
A mesh sweep at $N \in \{50, 100, 200, 400, 800, 1600\}$ on the M1
equilibrium pn junction. Records $V_{bi}$, peak $|E|$, depletion
width $W$. Two error metrics:

1. **Error vs depletion approximation** — physics-model gap; plateaus
   on fine meshes at the analytical-formula error.
2. **Cauchy self-convergence ratios** between consecutive meshes —
   pure FE discretization error; should improve as $h \to 0$.

Acceptance: Cauchy ratios $\ge 1.8\times$ per doubling in the first
four refinements (before the physics-model plateau). $V_{bi}$ is set
by the BCs and is mesh-independent; reported but not gated.

### Phase 3 — Conservation checks

[`semi/verification/conservation.py`](../../semi/verification/conservation.py).

- **Charge conservation on `pn_1d` equilibrium.** $|\int q(p-n+N_D-N_A)|/(qN_\mathrm{ref}L) < 10^{-10}$.
  Latest: $1.5\times 10^{-17}$. ✓
- **Current continuity on `pn_1d_bias` and `pn_1d_bias_reverse`.**
  At each target bias, sample $J_\mathrm{total} = J_n + J_p$ on ten
  interior facets; require $\max|J - \bar J|/|\bar J| < 5\%$ forward,
  $< 15\%$ reverse. Worst case observed: 1.91% forward, 0.020% reverse.

### Phase 4 — MMS for coupled drift-diffusion

[`semi/verification/mms_dd.py`](../../semi/verification/mms_dd.py).
Three (now four with M16.1) variants exercise progressively more of
the residual:

- **Variant A (psi-only).** $\Phi_n^* = \Phi_p^* = 0$. Continuity rows
  collapse to noise; only the psi block is rate-gated.
- **Variant B (full coupling, no R).** All three fields nontrivial,
  lifetimes set to $10^{20}\,\mathrm{s}$ to make $R$ negligible.
- **Variant C (full coupling with SRH).** Realistic Si lifetimes
  $10^{-7}\,\mathrm{s}$.
- **Variant D (Caughey-Thomas mobility, M16.1).** Full coupling +
  field-dependent $\mu$. Acceptance gate $L^2\ge 1.99$, $H^1\ge 0.99$.

Latest finest-pair rates: every gated study clears $L^2\ge 1.99$ at the
M16.1 acceptance gate; A/B/C variants clear $L^2\ge 1.85$ floor with
$\ge 0.19$ headroom; 1D and 2D paths consistent.

## Method of Manufactured Solutions

Pick an exact solution $u^*(\mathbf{x})$ that you *know* (e.g.
$u^*(x) = \sin(\pi x/L)$). Plug it into the strong-form residual:

$$
F_\mathrm{strong}[u^*] = -\nabla\!\cdot\!(\varepsilon\nabla u^*) + (\text{source}) - \rho^*
\equiv f^*(\mathbf{x}),
$$

with $f^*$ in general nonzero. Now redefine the right-hand side of
your PDE to be $\rho^*(u^*) + f^*$. By construction, $u^*$ is the exact
solution of the modified problem.

Discretize the modified problem on a sequence of meshes; compute the
error $\|u_h - u^*\|$ at each level. The error should scale as
$h^{p+1}$ in $L^2$ and $h^p$ in $H^1$ where $p$ is the polynomial
degree (P1 → $p=1$). Observed rate $\to 2$ in $L^2$ and $\to 1$ in $H^1$
with mesh refinement.

If the observed rate is *less than* the theoretical: the discrete
operator is missing a derivative or has a sign error. MMS catches
silent bugs that single-point validation misses.

The forcing term must be added *in weak form*: do not compute
$f^* = -\nabla\!\cdot\!(\varepsilon\nabla u^*)$ symbolically and then
integrate it against $v$ — at coarse mesh resolution, the symbolic
discretization can collapse to numerical zero. Instead, build $u^*$ as
a UFL expression and compose the strong form *as UFL*, which the
assembler then integrates by parts in weak form
([`semi/verification/mms_poisson.py`](../../semi/verification/mms_poisson.py)
follows this rule).

## Reproducing `tests/check_analytical_math.py`

This script asserts the basic physics math without dolfinx. Reproducing
each assertion from formulas in this guide:

**Thermal voltage** (Ch. 3).
$V_t = kT/q = 1.381\times 10^{-23}\cdot 300/1.602\times 10^{-19}
= 0.025852\,\mathrm{V}$. Asserted $|V_t - 0.02585| < 0.001$. ✓

**$\lambda^2$ bare** (Ch. 12).
$\lambda^2 = \varepsilon_0 V_t/(qC_0L_0^2)
= 8.854\times 10^{-12}\cdot 0.02585/(1.602\times 10^{-19}\cdot 10^{23}\cdot 4\times 10^{-12})
= 3.57\times 10^{-6}$. Asserted matches.

**Debye length in silicon** (Ch. 12).
$L_D = \sqrt{\varepsilon_r\varepsilon_0V_t/(qC_0)}
= \sqrt{11.7\cdot 8.854\times 10^{-12}\cdot 0.02585/(1.602\times 10^{-19}\cdot 10^{23})}
= 1.296\times 10^{-8}\,\mathrm{m} = 12.96\,\mathrm{nm}$. Asserted
$10 < L_D < 16\,\mathrm{nm}$. ✓

**Cross-check $(L_D/L_0)^2 = \lambda^2\varepsilon_r$** (Ch. 12 Exercise).
$(L_D/L_0)^2 = (1.296\times 10^{-8}/2\times 10^{-6})^2 = 4.20\times 10^{-5}$.
$\lambda^2\varepsilon_r = 3.57\times 10^{-6}\cdot 11.7 = 4.18\times 10^{-5}$.
Match within 0.5%. ✓ (Asserted $|...|/(...) < 0.01$.)

**$V_{bi}$ via asinh and via log** (Ch. 7).
asinh: $V_t[\mathrm{asinh}(N_A/(2n_i)) + \mathrm{asinh}(N_D/(2n_i))]
= V_t[\mathrm{asinh}(5\times 10^6) + \mathrm{asinh}(5\times 10^6)]
= V_t\cdot 2\cdot 16.118 = 0.8334\,\mathrm{V}$.
log: $V_t\ln(N_AN_D/n_i^2) = V_t\ln(10^{14}) = 0.8334\,\mathrm{V}$.
Relative diff $< 0.01$. ✓

**Depletion width $W$, $x_p$, $x_n$** (Ch. 7).
$N_\mathrm{eff} = N/2 = 5\times 10^{22}\,\mathrm{m^{-3}}$.
$W = \sqrt{2\varepsilon V_{bi}/(qN_\mathrm{eff})}
= \sqrt{2\cdot 1.036\times 10^{-10}\cdot 0.834/(1.602\times 10^{-19}\cdot 5\times 10^{22})}
= 1.469\times 10^{-7}\,\mathrm{m} = 146.9\,\mathrm{nm}$.
$x_p = x_n = W/2 = 73.4\,\mathrm{nm}$. Asserted at 1% relative.

**Peak $|E|$** (Ch. 7).
$|E_\mathrm{max}| = qN_Ax_p/\varepsilon
= 1.602\times 10^{-19}\cdot 10^{23}\cdot 7.34\times 10^{-8}/1.036\times 10^{-10}
= 1.135\times 10^{7}\,\mathrm{V/m} = 113.5\,\mathrm{kV/cm}$. Asserted in $[90, 130]$ kV/cm. ✓

**Charge balance** (Ch. 4).
$x_pN_A = 7.34\times 10^{-8}\cdot 10^{23} = 7.34\times 10^{15}$;
$x_nN_D$ same. Match to roundoff. ✓

**Bulk carriers** (Ch. 7).
$n^R = n_i\exp(\psi_R^\mathrm{eq}/V_t) = 10^{16}\cdot e^{16.12} = 10^{23}\,\mathrm{m^{-3}}$
($= N_D$, ✓). $p^R = n_i\exp(-\psi_R^\mathrm{eq}/V_t) = 10^{9}\,\mathrm{m^{-3}}
= 10^3\,\mathrm{cm^{-3}}$.

**Mass action $np = n_i^2$** (Ch. 3).
$10^{23}\cdot 10^9 = 10^{32}\,\mathrm{m^{-6}} = n_i^2$. ✓

**Charge neutrality in p-bulk** (Ch. 4).
$p^L - n^L - N_A = ?$ Substituting the asinh-derived $\psi_L$, the
identity $p^L - n^L = N_A$ holds exactly (Ch. 4 §charge neutrality).
Asserted residual / $N_A$ < $10^{-10}$. ✓

All ten checks reproduce.

## Per-device benchmarks

### `pn_1d` (M1)

Equilibrium Poisson, symmetric $10^{17}/10^{17}$ junction.
Verifier: $V_{bi}, |E_\mathrm{max}|, W$ vs depletion-approximation
formulas. Acceptance: 5-10% per quantity.

### `pn_1d_bias` (M2)

Forward-bias sweep, 0 to 0.6 V. Verifier: I-V vs Shockley long-diode
(plus SNS recombination correction at low bias). Acceptance: 10% on
$V \ge 0.5\,\mathrm{V}$.

### `pn_1d_bias_reverse` (M3)

Reverse-bias sweep, 0 to -2 V. Verifier: I-V vs SRH net generation
$J_\mathrm{gen}(V) = qn_i(W(V) - W(0))/(2\tau_\mathrm{eff})$. Acceptance: 20%
on $V \in [-2, -0.5]\,\mathrm{V}$.

### `mos_2d` (M6)

2D MOSCAP, multi-region (Si + SiO₂). Verifier: C-V via
$dQ_\mathrm{gate}/dV_g$ from `numpy.gradient(Q, V)` vs depletion-
approximation theory. Acceptance: 10% in $[V_{fb}+0.2, V_T-0.1]\,\mathrm{V}$.

### `resistor_3d` (M7)

3D doped resistor, 1 µm × 200 nm × 200 nm. Verifier: V-I linearity vs
$R_\mathrm{theory} = L/(qN_D\mu_n A)$. Acceptance: 1%. Run on both
builtin Cartesian mesh and gmsh-sourced unstructured mesh.

### `mosfet_2d` (M12 / M14.3)

2D n-channel MOSFET with Gaussian source/drain implants. Verifier:
$I_D/W$ vs Pao-Sah linear regime in $V_{GS}\in[V_T+0.2, V_T+0.6]\,\mathrm{V}$
at $V_{DS} = 0.05\,\mathrm{V}$. Acceptance: 20%.

### `rc_ac_sweep` (M14)

1D pn diode, AC sweep at $V_\mathrm{DC} = -1\,\mathrm{V}$, 41 frequencies
1 Hz to 1 GHz. Verifier: $C(\omega)$ plateau vs analytical
depletion C. Acceptance: 5% over [1 Hz, 1 MHz]; observed 0.41% worst-case.

### `pn_1d_turnon` (M13/M13.1)

Transient 1D pn-diode turn-on, $V_F: 0\to 0.6\,\mathrm{V}$ at $t=0$.
Verifier: $\tau_\mathrm{eff}$ extracted from late-time $I(t)$ vs input
$\tau_p$. Acceptance: 5%.

### `moscap_axisym_2d` (M14.2)

Axisymmetric 2D MOSCAP, 50 µm gate radius, Hu Fig. 5-18 parameters.
Verifier: LF and HF C-V curves vs analytical helpers in `semi/cv.py`.
Acceptance: visual + small percent on individual points; the notebook
05 plot is the deliverable.

### `diode_velsat_1d` (M16.1)

1D pn diode with Caughey-Thomas mobility. Verifier: divergence
between constant-mobility and Caughey-Thomas at $V_F = 0.9\,\mathrm{V}$
(56% expected) and convergence at $V_F = 0.3\,\mathrm{V}$ (0.19%
expected). Anchors that the M16.1 mobility model is wired correctly.

### `poisson_3d_gpu` (M15)

3D Poisson box on $80^3$ mesh, GPU-vs-CPU bit-identity check.
Acceptance: $\|\psi_\mathrm{gpu}-\psi_\mathrm{cpu}\|_\infty < 10^{-8}$;
$T_\mathrm{cpu}/T_\mathrm{gpu} \ge 5$ on the linear solve.

## Tolerance philosophy

Each benchmark's tolerance reflects the *physics-model gap* between the
engine's PDE and the reference's closed-form expression:

- 1% (resistor): closed-form is exact in the limit; tolerance just
  absorbs FE discretization.
- 5% (turn-on): exponential transient fit is sensitive to early-time
  noise.
- 10% (Shockley, MOSCAP): depletion approximation has ~5% bias; FEM
  matches to that limit.
- 20% (MOSFET Pao-Sah): long-channel + constant-mobility approximations
  are 10-15% off in the verifier window before discretization. M16.2
  Lombardi will tighten to 10%.

Accepting a tolerance is a statement about the *model gap*, not the
*code accuracy*. Pushing tolerances tighter without changing the model
would cause spurious CI failures.

## Code map

| Concept | Code location |
|---|---|
| MMS Poisson | `semi/verification/mms_poisson.py` |
| MMS DD | `semi/verification/mms_dd.py` |
| Mesh convergence | `semi/verification/mesh_convergence.py` |
| Conservation | `semi/verification/conservation.py` |
| Convergence helpers | `semi/verification/_convergence.py`, `_norms.py` |
| `run_verification.py` CLI | `scripts/run_verification.py` |
| `tests/check_analytical_math.py` | `tests/check_analytical_math.py` |
| `tests/test_axisym_moscap_math.py` | `tests/test_axisym_moscap_math.py` |
| `tests/test_diode_analytical.py` | `tests/test_diode_analytical.py` |
| Per-benchmark verifiers | `scripts/run_benchmark.py` |
| pn_1d_turnon verifier | `benchmarks/pn_1d_turnon/verify.py` |
| Mobility MMS Variant D | `semi/verification/mms_dd.py` (M16.1 acceptance) |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §5](../PHYSICS.md) — V&V suite results table.
- [`docs/adr/0006-verification-and-validation-strategy.md`](../adr/0006-verification-and-validation-strategy.md) — locked decisions.
- [`docs/mms_dd_derivation.md`](../mms_dd_derivation.md) — derivation of the MMS-DD forcing.
- [`docs/mos_derivation.md`](../mos_derivation.md) — multi-region MMS.
- [`docs/resistor_derivation.md`](../resistor_derivation.md) — 3D resistor verifier.

## Common pitfalls

1. **MMS forcing in strong form.** Computing $f^* = -\nabla\!\cdot\!(\varepsilon\nabla u^*)$
   symbolically and then integrating against $v$ in the discrete form
   can produce spurious zeros at coarse mesh. Always inject the
   manufactured solution into the weak form via UFL composition; let
   the assembler integrate by parts.
2. **Plateau in the depletion-approximation comparison.** The
   pn_1d mesh-convergence test at $N \ge 800$ shows the error
   *plateauing* at ~3-5% — that is the depletion-approximation
   reference's own error, not the engine's. The Cauchy ratio gate is
   the proper convergence check; the depletion-error comparison is
   informational only.
3. **Variant A's continuity rows are noise.** When $\Phi_n^* = \Phi_p^* = 0$,
   the continuity terms collapse to machine-roundoff residuals. The
   reported "rates" on those rows are meaningless. ADR 0006 documents
   that Variant A gates only the psi block.
4. **SNES atol matters at fine mesh.** At $N = 1600$, the absolute
   residual scale shrinks; the default `snes_atol = 1e-7` may declare
   convergence prematurely. The MMS-DD runner sets `atol = 0.0` with
   a tight `stol = 1e-12` ([`docs/mms_dd_derivation.md`](../mms_dd_derivation.md)
   Amendment).
5. **Validation $\neq$ verification.** A benchmark passing within
   tolerance does not prove the code is correct; it proves the code
   reproduces a closed-form approximation. MMS proves the discrete
   operator is consistent with the PDE; conservation proves global
   integral identities. Three witnesses agree → the code is correct.

## Exercises

**Exercise 21.1.** Pick a manufactured solution $\psi^*(x) = \sin(2\pi x/L)$
on a 1D Poisson problem with $\varepsilon = 1$. Compute the strong-form
residual $f^*$. Plug the modified problem into the engine and predict
the observed L² rate.

**Exercise 21.2.** Derive the SNS recombination current correction in
[`semi/diode_analytical.py:56-93`](../../semi/diode_analytical.py).
Why the factor $2V_t/(V_{bi} - V)$?

**Exercise 21.3.** A user's CI run fails with "MMS-DD Variant C 1D
$\phi_n$ rate 1.84 < 1.99". List three possible causes and one
diagnostic step for each.

**Exercise 21.4.** Reproduce the assertion at [`tests/check_analytical_math.py:46`](../../tests/check_analytical_math.py)
that the asinh- and log-form $V_{bi}$ agree to within 0.01%. Show the
algebra.

**Exercise 21.5.** A new milestone introduces a Schottky-contact
benchmark (M16.5). Predict what verification activities it will need:
MMS? mesh convergence? conservation? validation?

### Solutions

**21.1.** $\psi^* = \sin(2\pi x/L)$. $\nabla\psi^* = (2\pi/L)\cos(...)$.
$\nabla\cdot(\varepsilon\nabla\psi^*) = -(2\pi/L)^2\sin(...)$.
$f^* = -(-1\cdot(2\pi/L)^2\sin(...)) = (2\pi/L)^2\sin(2\pi x/L) = (2\pi/L)^2\psi^*$.
Inject into the runner. Observed L² rate should be 2.0 (P1 Lagrange,
smooth solution).

**21.2.** The Sah-Noyce-Shockley factor accounts for the carrier-density
gradient inside the depletion region. The naive
$J_\mathrm{rec} = qn_iW/(2\tau_\mathrm{eff})\cdot(\exp(V/(2V_t))-1)$
overestimates because it assumes uniform $n,p$ across the depletion
edges. The leading-order Sze correction divides by $1$ in deep depletion
and by $2V_t/(V_{bi}-V)$ in the transition; this brings the
recombination current into 10-15% agreement with the engine's full
Slotboom solution at the M3 device. See ADR 0006 for the engineering
context.

**21.3.** (a) **SNES tolerance too loose.** Diagnostic: lower
`snes_atol` to $0$ with tight `stol`; rerun. (b) **Mesh too coarse for
the manufactured wavenumber.** Diagnostic: refine by another factor of
2 and check rate convergence. (c) **Recently introduced sign error in
the SRH residual.** Diagnostic: run Variant B (no R) and confirm rates
clear; if Variant B is fine, the bug is in R; check
`semi/physics/recombination.py` for sign changes since last green CI.

**21.4.** $V_{bi}^\mathrm{asinh} = V_t[\mathrm{asinh}(N_D/(2n_i)) + \mathrm{asinh}(N_A/(2n_i))]$.
For $N \gg n_i$: $\mathrm{asinh}(N/(2n_i)) = \ln(N/(2n_i) + \sqrt{(N/(2n_i))^2 + 1})
\approx \ln(N/n_i)$. Sum: $V_t[\ln(N_D/n_i) + \ln(N_A/n_i)] = V_t\ln(N_AN_D/n_i^2)
= V_{bi}^\mathrm{log}$. The remainder is
$O((n_i/N)^2)$, of order $10^{-12}$ for $N/n_i = 10^7$ — well below 0.01%.

**21.5.** New benchmark `schottky_1d` will need:
- **Validation**: thermionic-emission analytical I-V from (8.3); the
  benchmark verifier per [`docs/IMPROVEMENT_GUIDE.md` §M16.5](../IMPROVEMENT_GUIDE.md).
- **MMS**: a manufactured solution that exercises the new Robin-style
  boundary form (8.7). MMS rate $\ge 1.99$ on the $\Phi_n$ block
  needed to prove the new BC machinery doesn't lose order.
- **Mesh convergence**: optional; the validation acceptance is the
  primary gate, and Schottky physics has its own depletion-region
  scale that mesh convergence would expose.
- **Conservation**: same total-current continuity check as
  pn_1d_bias; the Schottky contact's BC must produce conservative
  currents. This is automatic if the FEM form is consistent.

## Further reading

- **Roache, *Verification and Validation in Computational Science and
  Engineering* (1998).** The reference for V&V terminology.
- **Oberkampf and Roy, *Verification and Validation in Scientific
  Computing* (2010).** Modern textbook treatment.
- **Salari and Knupp, *Code Verification by the Method of Manufactured
  Solutions* (2000).** SAND2000-1444. The MMS practitioner's bible.
- **`docs/adr/0006-verification-and-validation-strategy.md`** in this
  repo — load-bearing.
- **`docs/mms_dd_derivation.md`** — MMS-DD derivation with all sign
  and scaling details.
- **Brenner and Scott, *The Mathematical Theory of Finite Element
  Methods* (3rd ed., 2008), Chapter 5.** Theoretical $L^2$ and $H^1$
  convergence rates for P1 Lagrange.
