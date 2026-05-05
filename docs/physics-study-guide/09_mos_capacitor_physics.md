# 9 — MOS capacitor physics

## Learning objectives

- Sketch the band diagram of a MOS capacitor in accumulation, flat-band,
  depletion, and inversion.
- Derive the flat-band voltage $V_{fb}$, the threshold voltage $V_T$,
  the maximum depletion width $W_\mathrm{dmax}$, the oxide capacitance
  $C_{ox}$, and the minimum HF capacitance $C_\mathrm{min}$.
- Distinguish low-frequency (LF) from high-frequency (HF) C–V curves
  and explain why kronos-semi computes HF via a depletion-approximation
  clamp rather than a true AC small-signal solve.
- Reproduce the analytical anchors at the Hu Fig. 5-18 reference parameters
  ($N_a = 5\times 10^{16}\,\mathrm{cm^{-3}}$, $T_{ox} = 10\,\mathrm{nm}$,
  $\phi_{ms} = -0.95\,\mathrm{V}$): $V_{fb} = -0.950\,\mathrm{V}$,
  $V_T = +0.181\,\mathrm{V}$, $|\phi_B| = 0.399\,\mathrm{V}$,
  $W_\mathrm{dmax} = 144\,\mathrm{nm}$, $C_\mathrm{min}/C_{ox} = 0.173$.
- Recognize the kronos-specific BC-convention shift $V_{fb} = \phi_{ms} - \phi_F$
  (vs the textbook $V_{fb} = \phi_{ms}$) and explain why it differs.

## Physical motivation

The MOS capacitor — metal / oxide / semiconductor — is the gate of every
MOSFET and the unit cell of every DRAM bit. Applying a voltage to the
metal modulates the carrier density at the oxide–semiconductor interface
through a purely electrostatic mechanism (no current flows in steady
state, since the oxide is insulating). The C(V) characteristic of a MOS
capacitor encodes the doping, the oxide thickness, the work-function
difference, and the interface trap density — making it the workhorse
diagnostic for MOSFET process technology.

kronos-semi's M6 benchmark and M14.2 axisymmetric refresh both target
the MOSCAP for two reasons: (1) it exercises the multi-region (Si/SiO₂)
infrastructure (Ch. 14); (2) the textbook closed-form C–V curve gives a
hard verification target. The shipped engine reproduces Hu's Fig. 5-18
HF C–V curve to within a few percent.

## Derivation from first principles

### Device and band diagram

A MOS capacitor is gate / oxide / body. We treat the n-MOS / p-body case
(positive gate voltage induces electron inversion); pMOS is symmetric
under $N \to N$, $V \to -V$. The body has acceptor doping $N_a$ and an
ohmic back contact. The gate is a metal with work function $\Phi_m$.

Define $\phi_{ms} = \Phi_m - \Phi_s$ where $\Phi_s$ is the semiconductor
work function (Ch. 2). For an n+ poly-silicon gate on p-Si at $N_a =
5\times 10^{16}\,\mathrm{cm^{-3}}$, $\phi_{ms} \approx -0.95\,\mathrm{V}$
(this is the value in [`semi/cv.py:69`](../../semi/cv.py)).

Three voltage regimes (taking V_g as the gate-to-body voltage):

```mermaid
flowchart LR
  A[Vg < Vfb : accumulation] --> B[Vg = Vfb : flat band]
  B --> C[Vfb < Vg < VT : depletion]
  C --> D[Vg > VT : inversion]
```

- **Accumulation.** Holes accumulate at the surface; the surface is
  more p-type than the bulk; the band bending is small and negative.
- **Flat band.** The applied voltage exactly cancels the
  work-function-induced band bending. Inside the silicon, $\psi$ is
  uniform.
- **Depletion.** The surface is depleted of holes; ionized acceptors are
  exposed; band bending grows as $V_g$ rises above $V_{fb}$.
- **Inversion.** The surface electron density exceeds the bulk hole
  density; an "inversion layer" forms at the oxide interface.

### Flat-band voltage

At flat band, the silicon's $\psi$ is uniform and equal to its bulk
equilibrium value $\psi_b$. The gate's potential is $\psi_b + \phi_{ms}$
(work-function offset from the silicon's). With the engine's BC convention
(Ch. 1), $\psi_b = -|\phi_B|$ for a p-body (negative, i.e. below the
intrinsic level). So

$$
V_{fb} = \phi_{ms} + (-|\phi_B|) - (-|\phi_B|) = \phi_{ms}\quad\text{(textbook)}.
\tag{9.1a}
$$

But the engine's body-ohmic BC pins $\psi_\mathrm{body} = \psi_b = -|\phi_B|$
on the back contact. The flat-band condition is then $\psi_\mathrm{gate}
= \psi_\mathrm{body}$, giving (with (8.6) for the gate):

$$
V_{fb} = \phi_{ms} - \phi_F,
\tag{9.1b}
$$

where $\phi_F = -|\phi_B|$ on a p-body, so $V_{fb} = \phi_{ms} + |\phi_B|$
in absolute terms. *In the engine's BC convention*, the form is (9.1b);
this differs from the textbook (9.1a) by $\phi_F$. See
[`docs/PHYSICS.md` §6.3](../PHYSICS.md) for the full discussion. The
analytical helper [`semi/cv.py:127`](../../semi/cv.py) uses the textbook
form internally (`V_fb = phi_ms - Q_f/C_ox`), and the engine uses the
shifted form for the BC; that is the convention disconnect that the
M6 derivation calls out.

### Bulk Fermi potential and surface potential

For a p-body, $|\phi_B| = V_t\,\ln(N_a/n_i)$. Strong inversion is
defined by the surface potential reaching $\psi_s = 2|\phi_B|$ (the
electron density at the surface equals the bulk hole density).

### Threshold voltage

To reach $\psi_s = 2|\phi_B|$, the gate must supply enough voltage to
support the depletion charge and bend the bands by $2|\phi_B|$:

$$
V_T = V_{fb} + 2|\phi_B| + \frac{\sqrt{2\varepsilon_s q N_a (2|\phi_B|)}}{C_{ox}},
\tag{9.2}
$$

with $C_{ox} = \varepsilon_{ox}/T_{ox}$ the per-area oxide capacitance.
The square-root term is the depletion charge per unit area at strong
inversion divided by $C_{ox}$ — the voltage drop across the oxide due
to the maximum depletion charge.

### Maximum depletion width

At strong inversion the depletion region freezes at:

$$
W_\mathrm{dmax} = \sqrt{\frac{2\varepsilon_s\cdot 2|\phi_B|}{q N_a}}.
\tag{9.3}
$$

(Same form as the pn-junction depletion-width formula, with
$V_{bi} \to 2|\phi_B|$ and $N_\mathrm{eff} \to N_a$.)

### Oxide capacitance and series stack

The oxide is a thin parallel-plate capacitor:

$$
C_{ox} = \varepsilon_{ox}/T_{ox}.
\tag{9.4}
$$

The depletion region acts like a second capacitor in series, with
voltage-dependent thickness $W_\mathrm{dep}(V_g)$:

$$
\frac{1}{C(V_g)} = \frac{1}{C_{ox}} + \frac{1}{C_\mathrm{dep}(V_g)},
\qquad
C_\mathrm{dep} = \frac{\varepsilon_s}{W_\mathrm{dep}}.
\tag{9.5}
$$

In accumulation, $W_\mathrm{dep} \to 0$ and $C \to C_{ox}$.
In strong depletion, $W_\mathrm{dep}$ grows; $C$ falls.
At inversion, the inversion layer at the surface either
- behaves like a metal plate (LF: minority carriers respond at the AC
  signal): $W_\mathrm{dep}$ effectively becomes zero, $C \to C_{ox}$;
- can't follow the AC signal (HF: minority generation is too slow):
  $W_\mathrm{dep}$ stays at $W_\mathrm{dmax}$, $C$ stays at
$$
C_\mathrm{min} = \frac{1}{1/C_{ox} + W_\mathrm{dmax}/\varepsilon_s}.
\tag{9.6}
$$

This is the LF–HF distinction that the C–V community lives by.

### LF and HF curves: kronos-semi's implementation

**Low-frequency (quasi-static).** The textbook formula in
[`semi/cv.py:210-285`](../../semi/cv.py) (`lf_cv_quasistatic`) inverts
the body-charge function $F(\psi_s)$ for $\psi_s(V_g)$ via Brent
root-finding, then differentiates $Q_s(V_g)$ via centered FD on a fine
grid. The FEM postprocessor `compute_lf_cv_fem` does the same FD on a
$Q_s$-vs-$V_g$ table from a swept FEM solve.

**High-frequency.** The strict definition is "AC capacitance with
frozen minority carriers." The proper way to compute it is an AC
small-signal sweep at high $\omega$ (Ch. 18). kronos-semi instead uses
a **depletion-approximation clamp**: for $\psi_s > 2|\phi_B|$, freeze
$W_\mathrm{dep} = W_\mathrm{dmax}$ and use (9.5). This matches the
textbook HF curve to a few percent on the standard reference parameters
and avoids implementing AC for the multi-region MOSCAP today. See
[`docs/theory/moscap_cv.md`](../theory/moscap_cv.md) for the discussion;
the rigorous AC HF C–V is deferred (post-M14.2; see ROADMAP §Deferred).

### Differential capacitance via AC admittance (M14.1)

The runner [`semi/runners/mos_cap_ac.py`](../../semi/runners/mos_cap_ac.py)
solves the *Poisson sensitivity* $\partial\psi/\partial V_g$ at each
operating point. The sensitivity satisfies a linearization of the
nonlinear Poisson equation:

$$
K(\psi_0)\,\delta\psi = -\frac{\partial F}{\partial V_g}\,\delta V_g,
\tag{9.7}
$$

with $K = \partial F/\partial\psi$ the discrete Jacobian at the
converged state and $\partial F/\partial V_g$ concentrated at the gate
Dirichlet row. Because the gate's BC is $\psi_g = V_g - \phi_{ms}$,
$\partial\psi_g/\partial V_g = 1$, and the sensitivity solve is exactly
the per-unit-V_g response of the field. The differential capacitance
is then

$$
\frac{dQ_\mathrm{gate}}{dV_g}
   = -\frac{q}{W_\mathrm{lat}}\int_\mathrm{Si}
       n_i\bigl(e^{-\psi_0/V_t} + e^{\psi_0/V_t}\bigr)\,\delta\psi\,dV.
\tag{9.8}
$$

(With an r-weighted integrand in the axisymmetric case; see Ch. 15.)
This is the M14.1 audit anchor that confirms byte-identity with
`mos_cv` on Q_gate at all 42 gate voltages (i.e. (9.8) is exact, modulo
discretization).

## Key results

- $|\phi_B| = V_t\,\ln(N_a/n_i)$.
- $V_{fb} = \phi_{ms} - Q_f/C_{ox}$ (textbook); engine BC shift in §6.3.
- $V_T$: (9.2).
- $W_\mathrm{dmax}$: (9.3).
- $C_{ox}$: (9.4); $C_\mathrm{min}$: (9.6).
- LF/HF distinction: (9.5) with frozen vs free $W_\mathrm{dep}$.
- Differential capacitance from sensitivity: (9.8).

## Worked numerical example

Hu Fig. 5-18 parameters: $N_a = 5\times 10^{16}\,\mathrm{cm^{-3}}$,
$T_{ox} = 10\,\mathrm{nm}$, $\phi_{ms} = -0.95\,\mathrm{V}$,
$Q_f = 0$, $T = 300\,\mathrm{K}$, $\varepsilon_r^\mathrm{Si} = 11.7$,
$\varepsilon_r^\mathrm{SiO_2} = 3.9$.

**$|\phi_B|$.** $V_t\ln(5\times 10^{16}/10^{10}) = 0.02585\cdot\ln(5\times 10^6)
= 0.02585\cdot 15.42 = 0.399\,\mathrm{V}$. ✓

**$C_{ox}$.** $\varepsilon_{ox}/T_{ox} = 3.9\cdot 8.854\times 10^{-12}/10^{-8}
= 3.45\times 10^{-3}\,\mathrm{F/m^2} = 345\,\mathrm{nF/cm^2}$. ✓

**$V_{fb}$.** $\phi_{ms} - Q_f/C_{ox} = -0.95 - 0 = -0.950\,\mathrm{V}$. ✓

**$W_\mathrm{dmax}$.** $\sqrt{2\cdot 11.7\cdot 8.854\times 10^{-12}\cdot
2\cdot 0.399 / (1.602\times 10^{-19}\cdot 5\times 10^{22})}
= \sqrt{1.654\times 10^{-10}/8.012\times 10^3}
= \sqrt{2.06\times 10^{-14}} = 1.44\times 10^{-7}\,\mathrm{m} = 144\,\mathrm{nm}$. ✓

**Depletion-charge term in $V_T$.**
$\sqrt{2\cdot 11.7\cdot 8.854\times 10^{-12}\cdot 1.602\times 10^{-19}\cdot 5\times 10^{22}\cdot 2\cdot 0.399}/C_{ox}$.
The numerator under the sqrt: $2\cdot 1.036\times 10^{-10}\cdot 1.602\times 10^{-19}\cdot 5\times 10^{22}\cdot 0.798
= 1.32\times 10^{-6}\,\mathrm{C^2/m^4}$. sqrt gives $1.15\times 10^{-3}\,\mathrm{C/m^2}$.
Divide by $C_{ox} = 3.45\times 10^{-3}\,\mathrm{F/m^2}$: $0.333\,\mathrm{V}$.

**$V_T$.** $V_{fb} + 2|\phi_B| + \mathrm{depterm} = -0.950 + 0.798 + 0.333 = +0.181\,\mathrm{V}$. ✓

**$C_\mathrm{min}/C_{ox}$.**
$C_\mathrm{dep,max} = \varepsilon_s/W_\mathrm{dmax} = 11.7\cdot 8.854\times 10^{-12}/1.44\times 10^{-7}
= 7.19\times 10^{-4}\,\mathrm{F/m^2}$.
$C_\mathrm{min} = 1/(1/C_{ox} + 1/C_\mathrm{dep,max}) = 1/(290 + 1391) = 1/1681 = 5.95\times 10^{-4}\,\mathrm{F/m^2}$.
Ratio: $5.95\times 10^{-4}/3.45\times 10^{-3} = 0.173$. ✓

All five anchors match the [`semi/cv.py`](../../semi/cv.py) helper to
three decimal places. These are the assertions in
[`tests/test_axisym_moscap_math.py`](../../tests/test_axisym_moscap_math.py).

## Code map

| Concept | Equation | Code location |
|---|---|---|
| `analytical_moscap_params` | all | `semi/cv.py:62-156` |
| $|\phi_B|$ | – | `semi/cv.py:123-124` |
| $V_{fb}$ (textbook) | (9.1a) | `semi/cv.py:127` |
| $V_T$ | (9.2) | `semi/cv.py:135-139` |
| $W_\mathrm{dmax}$ | (9.3) | `semi/cv.py:131` |
| $C_{ox}$ | (9.4) | `semi/cv.py:126` |
| $C_\mathrm{min}$ | (9.6) | `semi/cv.py:141-142` |
| HF C–V depletion clamp | (9.5)+(9.6) | `semi/cv.py:159-207` (`hf_cv_depletion_approximation`) |
| LF C–V quasistatic | – | `semi/cv.py:210-285` (`lf_cv_quasistatic`) |
| FEM LF post-processor | – | `semi/cv.py:288-310` (`compute_lf_cv_fem`) |
| FEM HF clamp | – | `semi/cv.py:313-360` (`compute_hf_cv_depletion_clamp`) |
| AC differential C runner | (9.7), (9.8) | `semi/runners/mos_cap_ac.py:59-329` |
| Engine BC convention shift | (9.1b) | `docs/PHYSICS.md` §6.3 |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §6](../PHYSICS.md) — full M6 reference for the 2D MOSCAP including the BC-convention shift discussion.
- [`docs/theory/moscap_cv.md`](../theory/moscap_cv.md) — LF/HF and the depletion-approximation clamp rationale.
- [`docs/benchmarks/moscap_axisym_2d.md`](../benchmarks/moscap_axisym_2d.md) — landing page for the M14.2 axisymmetric reproduction of Hu Fig. 5-18.

## Common pitfalls

1. **Engine BC vs textbook $V_{fb}$.** The engine uses
   $V_{fb} = \phi_{ms} - \phi_F$ in its body-ohmic BC convention; the
   textbook is $V_{fb} = \phi_{ms}$ when $\psi$ is referenced from the
   bulk Fermi level rather than from the intrinsic level. The
   $V_{fb}$ reported by `analytical_moscap_params` is the textbook
   form; if you compare numerically against the FEM solve, account for
   the $\phi_F$ offset.
2. **HF clamp is not AC.** The clamp gives the correct asymptote
   (frozen $W$ at $\psi_s = 2|\phi_B|$) but does *not* capture
   transition-region dispersion. Real C–V at MHz frequencies has a
   smooth knee around $V_T$; the clamp has a kink.
3. **Inversion charge follows the surface potential, not the body.**
   In strong inversion the inversion-layer electron density grows
   exponentially with $\psi_s$, so a small-signal $\delta\psi_s$ is
   strongly amplified; this is what makes the differential
   capacitance peak in inversion.
4. **Oxide thickness in nm, but stored in m.** [`semi/cv.py`](../../semi/cv.py)
   takes `T_ox_m` in meters; the JSON benchmarks use SI all the way.
   `T_ox = 10 nm = 1e-8 m`.
5. **Polarity for pMOS.** Equation (9.2) is for n-MOS / p-body. A
   pMOS device with an n-body has $\phi_F > 0$ and the depletion
   charge term flips sign; `analytical_moscap_params` handles this
   via the `body_dopant` parameter ([`semi/cv.py:135-139`](../../semi/cv.py)).

## Exercises

**Exercise 9.1.** A MOS capacitor with $T_{ox} = 5\,\mathrm{nm}$ instead
of 10 nm, all other parameters as Hu Fig. 5-18. Compute $C_{ox}$, $V_T$,
$C_\mathrm{min}$. What about $V_{fb}$?

**Exercise 9.2.** Show that as $N_a$ increases, $W_\mathrm{dmax}$
decreases (in fact $W_\mathrm{dmax} \propto 1/\sqrt{N_a}\cdot \sqrt{\ln(N_a/n_i)}$).
Compute $W_\mathrm{dmax}$ at $N_a = 10^{18}\,\mathrm{cm^{-3}}$.

**Exercise 9.3.** Replace SiO₂ with HfO₂ ($\varepsilon_r = 25$) at the
same physical thickness $T_\mathrm{ox} = 10\,\mathrm{nm}$. Compute
$C_{ox}$ and the *equivalent oxide thickness* (EOT, the SiO₂ thickness
that would give the same $C_{ox}$).

**Exercise 9.4.** At what $V_g$ does $\psi_s$ first reach $2|\phi_B|$ on
the Hu device, given the depletion-approximation surface-potential
equation $V_g - V_{fb} = \psi_s + \sqrt{2\varepsilon_s q N_a \psi_s}/C_{ox}$?

**Exercise 9.5.** Walk through the [`mos_cap_ac.py`](../../semi/runners/mos_cap_ac.py)
runner for the M14.1 byte-identity test. What does it assert? How can a
differential-capacitance solve match a finite-difference solve to
machine precision?

### Solutions

**9.1.** $C_{ox}(5\,\mathrm{nm}) = 3.9\cdot 8.854\times 10^{-12}/5\times 10^{-9}
= 6.91\times 10^{-3}\,\mathrm{F/m^2}$ (twice the original).
Depletion term in $V_T$: same depletion charge but divided by twice the
$C_{ox}$, so half: $0.166\,\mathrm{V}$.
$V_T = -0.950 + 0.798 + 0.166 = +0.014\,\mathrm{V}$.
$C_\mathrm{min} = 1/(1/6.91\times 10^{-3} + W_\mathrm{dmax}/\varepsilon_s)
= 1/(145 + 1391) = 6.51\times 10^{-4}\,\mathrm{F/m^2}$.
$V_{fb}$ unchanged: $-0.950\,\mathrm{V}$ (independent of $T_{ox}$ when
$Q_f = 0$).

**9.2.** $W_\mathrm{dmax} = \sqrt{2\varepsilon_s\cdot 2|\phi_B|/(qN_a)}$.
$|\phi_B| \propto \ln(N_a)$ grows with doping, so the numerator grows
slowly while the denominator grows linearly: net $W \propto \sqrt{\ln/N}$
decreases. At $N_a = 10^{18}\,\mathrm{cm^{-3}}$:
$|\phi_B| = V_t\ln(10^8) = 0.476\,\mathrm{V}$;
$W_\mathrm{dmax} = \sqrt{2\cdot 1.036\times 10^{-10}\cdot 0.952/(1.602\times 10^{-19}\cdot 10^{24})}
= \sqrt{1.97\times 10^{-10}/1.6\times 10^5} = \sqrt{1.23\times 10^{-15}} = 35\,\mathrm{nm}$.
About 1/4 of the $5\times 10^{16}$ value.

**9.3.** $C_{ox}^\mathrm{HfO_2} = 25\cdot 8.854\times 10^{-12}/10^{-8}
= 2.21\times 10^{-2}\,\mathrm{F/m^2}$ (about 6.4× larger).
EOT = $\varepsilon_\mathrm{SiO_2}/C_{ox}^\mathrm{HfO_2} = 3.9\cdot 8.854\times 10^{-12}/2.21\times 10^{-2}
= 1.56\times 10^{-9}\,\mathrm{m} = 1.56\,\mathrm{nm}$. So a 10 nm HfO₂
acts like a 1.56 nm SiO₂. This is why high-$\varepsilon_r$ ("high-k")
dielectrics replaced SiO₂ in advanced CMOS — they give thin-EOT gate
insulators while keeping a physically thicker layer that does not tunnel.

**9.4.** Solve $V_g - V_{fb} = 2|\phi_B| + \sqrt{2\varepsilon_s qN_a\cdot 2|\phi_B|}/C_{ox}$.
At $\psi_s = 2|\phi_B| = 0.798\,\mathrm{V}$:
$V_g = V_{fb} + 0.798 + 0.333 = -0.950 + 1.131 = +0.181\,\mathrm{V}$.
This is exactly $V_T$, by definition.

**9.5.** [`mos_cap_ac.py`](../../semi/runners/mos_cap_ac.py) solves the
sensitivity (9.7) at each gate bias. The byte-identity test in audit
case 03 verifies that the integrated charge $Q_\mathrm{gate}(V_g)$ from
this runner exactly matches `mos_cv`'s $Q_\mathrm{gate}$ at the same
$V_g$ — which makes sense because both are computed from the *same*
converged $\psi_0$. Where they differ is in the *derivative*: `mos_cv`
uses `numpy.gradient(Q, V)` (finite-difference noise), `mos_cap_ac`
uses (9.8) (analytic). The byte-identity is on $Q$, not on $C$.

## Further reading

- **Hu, *Modern Semiconductor Devices for Integrated Circuits* (2010),
  Chapter 5.** The reference for everything in this chapter,
  including Fig. 5-18.
- **Sze and Ng (2007), Chapter 4** for a more rigorous derivation
  including the inversion charge integrals.
- **Nicollian and Brews, *MOS (Metal Oxide Semiconductor) Physics and
  Technology* (1982).** The C–V bible. Chapter 4 covers the
  body-charge function and quasi-static C–V; Chapter 8 covers HF and
  the trap-state effects.
- **Lopez-Villanueva, J. A., et al. (2003).** "On the validity of the
  parabolic effective-mass approximation in the inversion layer of
  silicon nMOSFETs." Comparison of classical depletion-approximation
  HF C–V to a full Schrödinger–Poisson treatment; relevant if you
  need to interpret why the engine's HF clamp is a few percent off
  from a "true" measurement.
