# 7 — pn junction physics

## Learning objectives

- Derive the built-in voltage of an abrupt pn junction in two equivalent
  forms ($V_t\ln(N_A N_D/n_i^2)$ and the asinh-difference form) and
  show they agree.
- Apply the depletion approximation to derive the depletion width $W$,
  the per-side widths $x_n, x_p$, and the peak field $|E_\mathrm{max}|$
  in closed form.
- Reproduce the numerical assertions in `tests/check_analytical_math.py`
  for the $10^{17}/10^{17}$ silicon junction: $V_{bi} \approx 0.83\,\mathrm{V}$,
  $W \approx 146\,\mathrm{nm}$, $|E_\mathrm{max}| \approx 113\,\mathrm{kV/cm}$.
- Sketch the band diagram and quasi-Fermi-level structure under
  forward and reverse bias and connect to the Shockley diode equation.
- Trace the M1, M2, and M3 benchmarks through the engine and identify
  which physical regime each verifier targets.

## Physical motivation

The pn junction is the building block of every semiconductor device: a
diode is a pn junction; a bipolar transistor is two of them; a MOSFET
has two of them in series with a gate; a solar cell is a pn junction in
the dark + a generation source. The physics of carrier exchange across
the junction — drift balancing diffusion at equilibrium, exponential
response to bias, generation-driven leakage at reverse bias — is the
mental model that everything else builds on.

The shipped engine validates against this model in three benchmarks
(M1: equilibrium, M2: forward bias, M3: reverse bias) which between
them exercise every term in the drift-diffusion residual. This chapter
does the analytical groundwork, with worked numbers from `pn_1d`'s
$10^{17}/10^{17}$ symmetric junction.

## Derivation from first principles

### The device

An abrupt pn junction is two semi-infinite slabs of doped silicon
joined at a metallurgical interface $x = 0$:

- Left ($x < 0$): p-type, $N_A = $ const, $N_D = 0$.
- Right ($x > 0$): n-type, $N_D = $ const, $N_A = 0$.

The shipped M1 benchmark ([`benchmarks/pn_1d/pn_junction.json`](../../benchmarks/pn_1d/pn_junction.json))
puts the junction at $x = 1\,\mu\mathrm{m}$ with 1 µm of p-bulk and
1 µm of n-bulk on either side, $N_A = N_D = 10^{17}\,\mathrm{cm^{-3}}$,
ohmic contacts at both ends.

### Built-in voltage

At thermal equilibrium $\Phi_n = \Phi_p = 0$ globally (Ch. 11). From
the asinh form (Ch. 4), the equilibrium potential in the bulks is

$$
\psi_L^\mathrm{eq} = V_t\,\mathrm{asinh}(-N_A/(2n_i)),
\qquad
\psi_R^\mathrm{eq} = V_t\,\mathrm{asinh}(+N_D/(2n_i)).
$$

The built-in voltage is the difference:

$$
V_{bi} = \psi_R^\mathrm{eq} - \psi_L^\mathrm{eq}
       = V_t\bigl[\mathrm{asinh}(N_D/(2n_i)) + \mathrm{asinh}(N_A/(2n_i))\bigr].
\tag{7.1}
$$

For $N \gg n_i$, $\mathrm{asinh}(x) \approx \ln(2x)$, giving

$$
V_{bi} \approx V_t\bigl[\ln(N_D/n_i) + \ln(N_A/n_i)\bigr]
        = V_t\,\ln\!\left(\frac{N_A N_D}{n_i^2}\right).
\tag{7.2}
$$

Equation (7.2) is the textbook expression. The two forms agree to
better than $V_t \ln 2 \cdot \text{(small correction)}$ for $N \gg n_i$;
see Ch. 4 Exercise 3.2 for the algebra.

### Depletion approximation

Mobile carriers in a thin region around $x=0$ are swept out by the
built-in field, leaving the ionized dopants exposed. The **depletion
approximation** assumes:

- Outside the depletion region, charge neutrality holds and $\psi$ is
  flat at the bulk equilibrium value.
- Inside the depletion region $-x_p < x < x_n$, mobile carriers are
  zero and the charge density is just the ionized dopants.

Under these assumptions, the Poisson equation in the depletion region
becomes:

$$
\frac{d^2\psi}{dx^2} = -\frac{q}{\varepsilon}\,N_\mathrm{net}(x)
   = \begin{cases}
       +qN_A/\varepsilon & -x_p < x < 0 \\
       -qN_D/\varepsilon & 0 < x < x_n
     \end{cases}.
\tag{7.3}
$$

(Sign convention: I have used $-d^2\psi/dx^2 = \rho/\varepsilon$ from
Ch. 1, with $\rho_\mathrm{p-side} = -qN_A$ and $\rho_\mathrm{n-side}
= +qN_D$.)

### Charge balance

The total exposed charge on each side must be equal in magnitude:

$$
N_A x_p = N_D x_n.
\tag{7.4}
$$

(The charge per unit area is $qN x$; if they were unequal there would
be net charge sloshing around, and the depletion-edge condition would
move until they balance.) For a symmetric junction $N_A = N_D$,
$x_p = x_n = W/2$.

### Solving for the field and potential

Integrate (7.3) once with the boundary condition $E(\pm x_n, \pm x_p) = 0$
(the field vanishes at the depletion edges, where mobile carriers
shield the field):

$$
E(x) = -\frac{d\psi}{dx}
     = \begin{cases}
         -(qN_A/\varepsilon)(x + x_p) & -x_p < x < 0 \\
         -(qN_D/\varepsilon)(x_n - x) & 0 < x < x_n
       \end{cases}.
\tag{7.5}
$$

The peak field at $x = 0$:

$$
\boxed{
|E_\mathrm{max}| = \frac{q N_A x_p}{\varepsilon} = \frac{q N_D x_n}{\varepsilon}.
}
\tag{7.6}
$$

(Both expressions are equal by (7.4), and they have to be: $E$ is
continuous across the junction.)

Integrate (7.5) once more, demanding $\psi(-x_p) = \psi_L^\mathrm{eq}$
and $\psi(+x_n) = \psi_R^\mathrm{eq}$:

$$
V_{bi} - V = \frac{q}{2\varepsilon}\bigl(N_A x_p^2 + N_D x_n^2\bigr),
\tag{7.7}
$$

where the left-hand side is $V_{bi}$ at thermal equilibrium and
$V_{bi} - V$ under applied forward bias $V > 0$ (the bias subtracts from
the built-in barrier).

### Solving for the widths

Combining charge balance (7.4) and the integrated potential (7.7):

$$
\boxed{
W(V) = x_n + x_p
    = \sqrt{\frac{2\varepsilon (V_{bi} - V)(N_A + N_D)}{q\,N_A N_D}},
}
\tag{7.8}
$$

$$
x_n = W\,\frac{N_A}{N_A + N_D},
\qquad
x_p = W\,\frac{N_D}{N_A + N_D}.
\tag{7.9}
$$

Combining (7.6) and (7.9):

$$
|E_\mathrm{max}| = \sqrt{\frac{2 q (V_{bi} - V) N_\mathrm{eff}}{\varepsilon}},
\qquad
N_\mathrm{eff} \equiv \frac{N_A N_D}{N_A + N_D}.
\tag{7.10}
$$

For the symmetric case $N_A = N_D$, $N_\mathrm{eff} = N/2$.

### Forward and reverse bias

Under bias $V$ (positive forward, by the convention applied to the
n-side contact relative to the p-side):

- $V > 0$ (forward): $V_{bi} - V$ shrinks; $W$ shrinks; the barrier is
  lowered and minority carriers flood across the junction. The
  Shockley diode equation $J = J_s(\exp(V/V_t) - 1)$ describes the
  resulting current; $J_s$ comes from solving the minority-carrier
  diffusion equation in each bulk (Ch. 5 worked example).
- $V < 0$ (reverse): $V_{bi} - V$ grows; $W$ grows; the barrier is raised
  and forward injection is exponentially suppressed. The remaining
  current is the **SRH generation current** integrated over the
  expanded depletion region, which the M3 reverse-bias verifier checks
  against $J_\mathrm{gen} = (qn_i/2\tau_\mathrm{eff})(W(V) - W(0))$.
  See [`semi/diode_analytical.py:96-121`](../../semi/diode_analytical.py).

Quasi-Fermi-level picture under bias (Ch. 11): $\Phi_n - \Phi_p = -V$
in the depletion region (the ohmic contacts pin $\Phi_n = \Phi_p = 0$
on the cathode side and $= V$ on the anode side, but inside the
depletion region they remain split by $\sim V$ because the SRH rate is
slow enough to keep them out of equilibrium).

## Key results

- Built-in voltage: (7.1) and (7.2).
- Depletion width: (7.8).
- Per-side widths: (7.9).
- Peak field: (7.6) and (7.10).
- Charge balance: (7.4).
- Shockley diode (forward): $J = J_s(e^{V/V_t} - 1)$ from Ch. 5.
- SRH generation (reverse): $J_\mathrm{gen} = (qn_i/2\tau)(W(V)-W(0))$ from Ch. 6.

## Worked numerical example

Reproduce every check in
[`tests/check_analytical_math.py`](../../tests/check_analytical_math.py)
for the $10^{17}/10^{17}$ symmetric Si junction at 300 K.

**Inputs.** $\varepsilon_r = 11.7$, $n_i = 10^{16}\,\mathrm{m^{-3}}$,
$N_A = N_D = 10^{23}\,\mathrm{m^{-3}}$, $V_t = 25.85\,\mathrm{mV}$.
$\varepsilon = \varepsilon_r\varepsilon_0 = 1.036\times 10^{-10}\,\mathrm{F/m}$.

**$V_{bi}$.** Log form:
$V_{bi} = V_t\ln(10^{23}\cdot 10^{23}/(10^{16})^2) = V_t\ln(10^{14})
= 0.02585\cdot 32.236 = 0.8334\,\mathrm{V}$. ✓

**$W$ at $V = 0$.** $N_\mathrm{eff} = N/2 = 5\times 10^{22}\,\mathrm{m^{-3}}$.
$W = \sqrt{2\cdot 1.036\times 10^{-10}\cdot 0.8334 / (1.602\times 10^{-19}\cdot 5\times 10^{22})}
   = \sqrt{1.727\times 10^{-10}/8.012\times 10^{3}}
   = \sqrt{2.156\times 10^{-14}}
   = 1.469\times 10^{-7}\,\mathrm{m}
   = 146.9\,\mathrm{nm}$. ✓
The test asserts $W \approx 146\,\mathrm{nm}$.

**$x_n, x_p$.** Symmetric: $x_n = x_p = W/2 = 73.4\,\mathrm{nm}$. ✓

**$|E_\mathrm{max}|$.**
$|E_\mathrm{max}| = qN_A x_p/\varepsilon
   = 1.602\times 10^{-19}\cdot 10^{23}\cdot 7.34\times 10^{-8}/1.036\times 10^{-10}
   = 1.176\times 10^{-3}/1.036\times 10^{-10}
   = 1.135\times 10^{7}\,\mathrm{V/m}
   = 113.5\,\mathrm{kV/cm}$. ✓
The test gates $90 < |E|/(\mathrm{kV/cm}) < 130$.

**Charge balance check.** $x_p N_A = 7.34\times 10^{-8}\cdot 10^{23} = 7.34\times 10^{15}$.
$x_n N_D = $ same number; ratio $1\pm 10^{-6}$. ✓

**Bulk carrier check.** From Ch. 4, $n_R = 10^{17}\,\mathrm{cm^{-3}}$,
$p_R = 10^3\,\mathrm{cm^{-3}}$. Both correct to 1%. ✓

All five tests pass.

### Forward-bias estimate

At $V = 0.6\,\mathrm{V}$ on the M2 device (same parameters, plus $\tau_n = \tau_p
= 10^{-7}\,\mathrm{s}$): $\exp(V/V_t) = e^{0.6/0.02585} = e^{23.21}
= 1.20\times 10^{10}$. Saturation current $J_s = 4.78\times 10^{-8}\,\mathrm{A/m^2}$
(Ch. 5 worked example, long-diode formula). Predicted Shockley current
$\approx 5.74\times 10^2\,\mathrm{A/m^2}$. Engine reports
$\sim 1.6\times 10^3\,\mathrm{A/m^2}$. The factor-of-3 discrepancy is
consistent with short-base correction; see Ch. 5 worked example.

### Reverse-bias estimate

At $V = -1\,\mathrm{V}$: $W(-1) = \sqrt{2\varepsilon(V_{bi}+1) N_A+N_D)/(qN_AN_D)}
= \sqrt{(0.834+1)/0.834}\cdot W(0) = \sqrt{2.20}\cdot 146.9 = 218\,\mathrm{nm}$.
$\Delta W = 218 - 147 = 71\,\mathrm{nm}$.
$|J_\mathrm{gen}| = qn_i\Delta W /(2\tau)
= 1.602\times 10^{-19}\cdot 10^{16}\cdot 7.1\times 10^{-8}/(2\cdot 10^{-7})
= 5.7\times 10^{-4}\,\mathrm{A/m^2}$.
The M3 verifier accepts 20% agreement; the engine match is well within.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| $V_{bi}$ via asinh BCs | (7.1) | `semi/bcs.py:181-189`, `semi/physics/slotboom.py:73-83` |
| $V_{bi}$ via log form | (7.2) | `semi/diode_analytical.py:50, 78, 112` |
| Depletion width $W(V)$ | (7.8) | `semi/diode_analytical.py:40-53` |
| Peak field | (7.6), (7.10) | `tests/check_analytical_math.py:57-61` |
| Charge balance | (7.4) | `tests/check_analytical_math.py:64-65` |
| Shockley $J_s$ | – | `semi/diode_analytical.py:19-37` (`shockley_long_diode_saturation`) |
| Sah-Noyce-Shockley + recombination | – | `semi/diode_analytical.py:56-93` (`sns_total_reference`) |
| SRH reverse-bias generation | – | `semi/diode_analytical.py:96-121` (`srh_generation_reference`) |
| M1 verifier (equilibrium) | – | `scripts/run_benchmark.py` `verify_pn_1d` |
| M2 verifier (forward) | – | `scripts/run_benchmark.py` `verify_pn_1d_bias` |
| M3 verifier (reverse) | – | `scripts/run_benchmark.py` `verify_pn_1d_bias_reverse` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §3.1, §5.3](../PHYSICS.md) — ohmic-contact BC, current-continuity gates.
- [`docs/PHYSICS_INTRO.md` §5](../PHYSICS_INTRO.md) — pn junction in plain English.
- [`benchmarks/pn_1d/README.md`](../../benchmarks/pn_1d/README.md) and the `pn_1d_bias{,_reverse}` READMEs.

## Common pitfalls

1. **The depletion approximation is not exact.** It treats the
   depletion edge as a sharp boundary; in reality the transition is
   smeared by $\sim L_D \sim 13\,\mathrm{nm}$ at $10^{17}\,\mathrm{cm^{-3}}$.
   The full-Poisson FEM solution shows a smooth transition; depletion
   approximation results overestimate $W$ by a few percent. See the
   M4 mesh-convergence study ([`docs/PHYSICS.md` §5.2](../PHYSICS.md))
   for the *physics-model gap* honest flag.
2. **Sign of $V$.** kronos-semi sets the *cathode* (n-side) ohmic
   contact to a swept voltage; the anode is at zero. So a "forward"
   bias of +0.6 V on the cathode is actually reverse — be careful with
   per-benchmark conventions. The `pn_1d_bias_reverse` config inverts
   this; check the JSON for which contact is swept.
3. **Asymmetric junction $\neq$ symmetric.** The depletion width
   formula (7.8) holds for any $N_A, N_D$, but the peak field migrates
   off-center: in a $10^{17}/10^{15}$ junction, the lightly-doped side
   carries most of the depletion ($x \propto 1/N$) and the field is
   asymmetric across $x=0$.
4. **$V_{bi}$ depends on $n_i$.** Using the older textbook
   $n_i = 1.45\times 10^{10}\,\mathrm{cm^{-3}}$ instead of Altermatt's
   $1.0\times 10^{10}$ shifts $V_{bi}$ by $V_t\ln((1.45/1.0)^2) = 19\,\mathrm{mV}$.
   For most device-level numbers this is below tolerance; for
   sub-threshold MOSFET subthreshold-slope predictions it matters.
5. **The Shockley equation assumes long-base, low-injection, ideal
   contacts.** Real diodes have ideality factor $n > 1$ at low bias due
   to SRH-recombination current (the Sah-Noyce-Shockley correction in
   `sns_total_reference`), and current saturation at high bias due to
   series resistance. The M2 verifier accepts 10% agreement above
   $V \geq 0.5\,\mathrm{V}$ specifically to absorb these.

## Exercises

**Exercise 7.1.** Compute $W$ at $V_F = 0.5\,\mathrm{V}$ on the M2
device. Compare with $W(V=0)$.

**Exercise 7.2.** For an asymmetric junction $N_A = 10^{15}$,
$N_D = 10^{17}$ in silicon: compute $V_{bi}$, $W$, $x_n$, $x_p$, and the
peak field at $V = 0$.

**Exercise 7.3.** Show that the Shockley diode equation reduces to
$J \approx J_s\exp(V/V_t)$ at $V \gg V_t$ and to $J \approx -J_s$ at
$V \ll -V_t$. Combine the two limits to explain "saturation current."

**Exercise 7.4.** Derive (7.7) from (7.5) by integrating $\psi(x)$ once
more. (Watch the sign of the integration.)

**Exercise 7.5.** What does a 10 V reverse bias do to the M1 device?
Compute $W$, $|E_\mathrm{max}|$. Is the device safe (silicon breakdown
field is $\sim 3\times 10^5\,\mathrm{V/cm}$)?

### Solutions

**7.1.** $W(0.5) = W(0)\sqrt{(V_{bi}-0.5)/V_{bi}}
= 146.9\,\mathrm{nm}\cdot\sqrt{0.334/0.834} = 146.9\cdot 0.633 = 92.9\,\mathrm{nm}$.
The depletion region shrinks by 37%.

**7.2.** $V_{bi} = V_t\ln(10^{15}\cdot 10^{17}/(10^{10})^2) = V_t\ln(10^{12})
= 0.02585\cdot 27.63 = 0.714\,\mathrm{V}$.
$N_\mathrm{eff} = 10^{15}\cdot 10^{17}/(10^{15}+10^{17}) \approx 10^{15}\,\mathrm{cm^{-3}}
= 10^{21}\,\mathrm{m^{-3}}$.
$W = \sqrt{2\cdot 1.036\times 10^{-10}\cdot 0.714 / (1.602\times 10^{-19}\cdot 10^{21})}
= \sqrt{1.479\times 10^{-10}/1.602\times 10^2}
= \sqrt{9.23\times 10^{-13}} = 9.61\times 10^{-7}\,\mathrm{m} = 961\,\mathrm{nm}$.
$x_p = W\cdot N_D/(N_A+N_D) = 961\cdot 10^{17}/(10^{17}+10^{15}) = 951\,\mathrm{nm}$,
$x_n = W - x_p = 10\,\mathrm{nm}$. The depletion is 99% on the lightly-doped p-side.
$|E_\mathrm{max}| = qN_D x_n/\varepsilon = 1.602\times 10^{-19}\cdot 10^{23}\cdot 10^{-8}/1.036\times 10^{-10}
= 1.55\times 10^6\,\mathrm{V/m} = 15.5\,\mathrm{kV/cm}$.

**7.3.** At $V \gg V_t$: $\exp(V/V_t)$ dominates the $-1$, so $J \to J_s e^{V/V_t}$
(exponential growth). At $V \ll -V_t$: $\exp(V/V_t) \to 0$, so $J \to -J_s$
(constant). $J_s$ is the *reverse saturation current* — the small,
near-temperature-independent current that flows when the junction is
back-biased and the only available current is minority-carrier
diffusion from the bulks.

**7.4.** Integrate (7.5) for $-x_p < x < 0$ once:
$\psi(x) - \psi_L^\mathrm{eq} = -\int_{-x_p}^x E(x')\,dx'
= (qN_A/\varepsilon)\int_{-x_p}^x (x' + x_p)\,dx'
= (qN_A/(2\varepsilon))(x + x_p)^2$.
At $x = 0$: $\psi(0) - \psi_L = qN_A x_p^2/(2\varepsilon)$.
Similarly on the n-side: $\psi_R - \psi(0) = qN_D x_n^2/(2\varepsilon)$.
Sum: $V_{bi} = (q/2\varepsilon)(N_A x_p^2 + N_D x_n^2)$, replacing
$V_{bi}$ with $V_{bi}-V$ under applied bias. ✓

**7.5.** $W(-10) = W(0)\sqrt{(V_{bi}+10)/V_{bi}} = 146.9\sqrt{10.83/0.83}
= 146.9\cdot 3.61 = 530\,\mathrm{nm}$.
$|E_\mathrm{max}| = $ scales as $\sqrt{V_{bi}-V}$, so
$E_\mathrm{max}(-10) = E_\mathrm{max}(0)\cdot 3.61 = 410\,\mathrm{kV/cm}$.
This exceeds the silicon breakdown field of $\sim 300\,\mathrm{kV/cm}$;
the M1 device would avalanche. Real diodes designed to handle this
voltage have lower doping (smaller $|E_\mathrm{max}|$ at the same $W$)
and a guard ring to keep the field uniform. M16.6 BBT models the
pre-breakdown leakage; avalanche multiplication is post-M16.

## Further reading

- **Sze and Ng (2007), Chapter 2** for the full pn-junction derivation
  in textbook form, including the abrupt-junction-with-Gaussian
  doping correction.
- **Pierret (1996), §5.1–5.3** for the depletion approximation and
  short-base diode discussion.
- **Hu (2010), Chapter 4.** Compact, modern treatment with worked
  examples that match common SPICE compact-model parameters.
- **Shockley, W. (1949).** "The theory of p-n junctions in semiconductors
  and p-n junction transistors." *Bell System Tech. J.* 28, 435.
  The original.
- **Sah, C. T., Noyce, R. N., and Shockley, W. (1957).** "Carrier
  generation and recombination in p-n junctions and p-n junction
  characteristics." *Proc. IRE* 45, 1228. The SNS correction used by
  the M3 reverse-bias verifier.
