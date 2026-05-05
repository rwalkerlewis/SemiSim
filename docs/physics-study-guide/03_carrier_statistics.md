# 3 — Carrier statistics

## Learning objectives

- Reduce the Fermi–Dirac distribution to the Boltzmann approximation and
  state the validity boundary explicitly.
- Derive $n = N_c\exp(-(E_c - E_F)/kT)$ and its dual.
- Connect the engine's working expressions $n = n_i\exp((\psi - \Phi_n)/V_t)$
  and $p = n_i\exp((\Phi_p - \psi)/V_t)$ to the band-edge form via the
  intrinsic level.
- Use mass action $np = n_i^2$ as a sanity check, and reproduce the
  numerical assertions about $V_t$, $n_i$, and bulk $n$, $p$ in
  `tests/check_analytical_math.py`.
- Anticipate where Boltzmann breaks down (degenerate doping above
  $\sim 10^{19}\,\mathrm{cm^{-3}}$) and locate the M16.4 Fermi–Dirac
  extension in the milestone backlog.

## Physical motivation

Whether a current flows when you apply a voltage depends on how many
mobile carriers there are and how they are distributed in energy.
Carriers are fermions, so the *most fundamental* statement of "how many
carriers at energy $E$" comes from the Fermi–Dirac distribution. But
solving the drift-diffusion equations with that nonlinearity baked in
is expensive and unnecessary at the doping levels that dominate today's
silicon technology. The Boltzmann approximation captures the same
physics with a closed-form exponential, and its validity boundary is a
well-defined doping threshold. This chapter quantifies all of that.

## Derivation from first principles

### Fermi–Dirac distribution

For a fermion with energy $E$ in equilibrium with a reservoir at
temperature $T$ and Fermi level $E_F$, the occupancy probability is

$$
f_\mathrm{FD}(E; E_F, T)
   = \frac{1}{1 + \exp\!\big((E - E_F)/kT\big)}.
\tag{3.1}
$$

The number density of electrons in the conduction band is the integral
over the band's density of states $g_c(E)$ weighted by occupancy:

$$
n = \int_{E_c}^{\infty} g_c(E)\,f_\mathrm{FD}(E; E_F, T)\,dE.
\tag{3.2}
$$

Using the parabolic-band density of states (2.4), substituting
$\eta = (E - E_c)/kT$ and $\eta_F = (E_F - E_c)/kT$, gives

$$
n = N_c\,F_{1/2}(\eta_F),
\qquad
F_{1/2}(\eta_F) = \frac{2}{\sqrt\pi}\int_0^\infty \frac{\sqrt\eta}{1 + e^{\eta - \eta_F}}\,d\eta,
\tag{3.3}
$$

where $N_c$ is the effective DOS from (2.6) and $F_{1/2}$ is the
Fermi–Dirac integral of order $1/2$. Implementing $F_{1/2}$ correctly is
the main numerical task of M16.4.

### Boltzmann limit

For $\eta_F < -3$ (equivalently $E_c - E_F > 3kT$), the integrand in
(3.3) is exponentially small and the $1/(1 + e^{\eta - \eta_F})$ in the
denominator can be replaced by $e^{\eta_F - \eta}$. The integral
collapses to $F_{1/2}(\eta_F) \approx e^{\eta_F}$, giving

$$
\boxed{
n = N_c \exp\!\left(\frac{E_F - E_c}{kT}\right),
\qquad
p = N_v \exp\!\left(\frac{E_v - E_F}{kT}\right).
}
\tag{3.4}
$$

Equation (3.4) is **Boltzmann statistics**. The error in (3.4) relative
to (3.3) is below 2% at $\eta_F = -3$, below 10% at $\eta_F = -1$, and
explodes for $\eta_F \geq 0$ (degenerate doping). The numerical
boundary kronos-semi flags is around $n \sim 10^{19}\,\mathrm{cm^{-3}}$
in silicon, which corresponds to $\eta_F \approx 0$.

### Intrinsic level and the engine's working form

Define the **intrinsic Fermi level** $E_i$ as the value of $E_F$ at
which $n = p$ in an undoped sample. Setting $n = p$ in (3.4) and
solving:

$$
E_i = \frac{E_c + E_v}{2} + \frac{kT}{2}\ln\!\left(\frac{N_v}{N_c}\right),
\qquad
n_i = \sqrt{N_c N_v}\,\exp\!\left(-\frac{E_g}{2kT}\right).
\tag{3.5}
$$

Define the electrostatic potential measured from $E_i$:

$$
\psi(\mathbf{x}) \equiv -\frac{E_i(\mathbf{x})}{q} + \mathrm{const}.
\tag{3.6}
$$

(Add an arbitrary constant; the engine pins it via the asinh formula at
the contacts.) Then $E_c = -q\psi - q\chi - E_g/2 + \mathrm{const}$ in
the bulk, and substituting (3.6) into (3.4):

$$
n = N_c\,\exp\!\left(\frac{E_F - E_c}{kT}\right)
   = n_i\,\exp\!\left(\frac{q\psi - (E_i - E_F)}{kT}\right).
$$

Define the **electron quasi-Fermi potential** $\Phi_n$ by $E_F \equiv -q\Phi_n$
under bias (see Ch. 11 for the bias case; at equilibrium $\Phi_n$ is a
single global constant, and the engine's convention pins it to zero).
Then

$$
\boxed{
n = n_i\,\exp\!\left(\frac{\psi - \Phi_n}{V_t}\right),
\qquad
p = n_i\,\exp\!\left(\frac{\Phi_p - \psi}{V_t}\right),
}
\tag{3.7}
$$

with $V_t = kT/q$ the **thermal voltage**. Equation (3.7) is the
working form used throughout `semi/physics/`. At equilibrium
$\Phi_n = \Phi_p = 0$ identically, so $n p = n_i^2$ falls out — that is
the **mass-action law** in disguise:

$$
\boxed{
n p \;=\; n_i^2 \;=\; N_c N_v \exp(-E_g/kT)
\quad\text{at thermal equilibrium}.
}
\tag{3.8}
$$

Mass action is exact under Boltzmann statistics; it fails by the same
factor that the Boltzmann-vs-FD ratio fails (a few percent at
$10^{19}\,\mathrm{cm^{-3}}$, an order of magnitude at
$10^{20}\,\mathrm{cm^{-3}}$).

## Key results

- Boltzmann: (3.4).
- Engine's working form: (3.7).
- Mass action: (3.8).
- Thermal voltage: $V_t = kT/q$. At $T = 300\,\mathrm{K}$,
  $V_t = (1.381\times 10^{-23})(300) / (1.602\times 10^{-19})
   = 25.85\,\mathrm{mV}$.

## Worked numerical example

The first half of [`tests/check_analytical_math.py`](../../tests/check_analytical_math.py)
exercises exactly the formulas in this chapter. Reproducing it line by
line:

**Thermal voltage at 300 K.**
$V_t = kT/q = 1.381\times 10^{-23} \cdot 300 / 1.602\times 10^{-19}
       = 0.025852\,\mathrm{V} = 25.85\,\mathrm{mV}$. ✓
(The test asserts $|V_t - 0.02585| < 0.001$.)

**Bulk electron density on the n-side of the M1 pn junction.**
With $N_D = 10^{17}\,\mathrm{cm^{-3}} = 10^{23}\,\mathrm{m^{-3}}$ and
$n_i = 10^{10}\,\mathrm{cm^{-3}} = 10^{16}\,\mathrm{m^{-3}}$:
$\psi_R^\mathrm{eq} = V_t\,\mathrm{asinh}(N_D/2n_i)$.
The argument is $10^{23}/(2 \cdot 10^{16}) = 5\times 10^6$;
$\mathrm{asinh}(5\times 10^6) = \ln(5\times 10^6 + \sqrt{25\times 10^{12}+1})
\approx \ln(10^7) = 7\ln 10 = 16.118$.
So $\psi_R^\mathrm{eq} = 0.02585 \cdot 16.118 = 0.4167\,\mathrm{V}$.
Then $n^R = n_i \exp(\psi_R/V_t) = 10^{10} \cdot e^{16.118}
= 10^{10} \cdot 10^7 \approx 10^{17}\,\mathrm{cm^{-3}} \approx N_D$, ✓ (the test
asserts $n_R \approx N_D$ to 1%).

**Mass action.**
$n^R \cdot p^R = (10^{17})(10^3) = 10^{20}\,\mathrm{cm^{-6}}$ (the
minority hole density is $p^R = n_i^2/N_D = 10^3\,\mathrm{cm^{-3}}$).
$n_i^2 = 10^{20}\,\mathrm{cm^{-6}}$. ✓

**Charge neutrality.**
On the p-bulk: $p - n - N_A = 10^{17} - 10^3 - 10^{17} \approx 0$ to
14 digits. The test asserts the residual is below $10^{-10} \cdot N_A$. ✓
The exact identity is the asinh(1) trick:
$\exp(-\psi_L/V_t) - \exp(\psi_L/V_t) = -2\sinh(\psi_L/V_t)
= 2\,\mathrm{asinh}(N_A/(2n_i))$ when $\psi_L = V_t\,\mathrm{asinh}(-N_A/(2n_i))$,
giving $p^L - n^L = N_A$ exactly.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| $V_t = kT/q$ | – | `semi/constants.py:24-26` (`thermal_voltage`); `semi/scaling.py:42-45` (`Scaling.V0`) |
| Boltzmann $n,p$ | (3.7) | `semi/physics/slotboom.py:27-42` (UFL); `:45-52` (numpy) |
| $n_i$ from material DB | (3.5) (loaded as a constant) | `semi/materials.py:33, 58, 70, 82` |
| Mass action sanity check | (3.8) | `tests/check_analytical_math.py:78-80` |
| Equilibrium $\psi^\mathrm{eq}$ | $V_t\,\mathrm{asinh}(N/2n_i)$ | `semi/physics/slotboom.py:73-83` (`equilibrium_psi_hat`); `semi/bcs.py:183, 253` |
| Validity-boundary forward link | $\sim 10^{19}\,\mathrm{cm^{-3}}$ | `docs/IMPROVEMENT_GUIDE.md` §M16.4 (planned) |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §1.2](../PHYSICS.md) — Boltzmann form used by the engine.
- [`docs/PHYSICS.md` §4.1, §4.2](../PHYSICS.md) — $V_t$ and $n_i$ at 300 K.
- [`docs/PHYSICS_INTRO.md` §3.2](../PHYSICS_INTRO.md) — narrative description of Slotboom and mass action.
- [`docs/IMPROVEMENT_GUIDE.md` §M16.4](../IMPROVEMENT_GUIDE.md) — Fermi–Dirac extension scope.

## Common pitfalls

1. **Boltzmann breaks down quietly.** Equation (3.4) gives perfectly
   well-defined numbers at $N_D = 10^{20}\,\mathrm{cm^{-3}}$; the
   numbers are simply *wrong* by a factor of order 2 because $\eta_F$ is
   positive. Source/drain extensions in real MOSFETs sit at $10^{19}$ to
   $10^{21}\,\mathrm{cm^{-3}}$. M16.4 will add the Fermi–Dirac branch
   gated by `physics.statistics: "fermi_dirac"`.
2. **$n_i$ is not a fundamental constant.** Different references quote
   different values: 1.45×10¹⁰ in older textbooks, 1.0×10¹⁰ in modern
   TCAD (Altermatt 2003, what kronos-semi uses), 1.08×10¹⁰ from the
   Sze tabulation. The discrepancy is from temperature, bandgap
   narrowing, and effective-mass values; for built-in voltage
   comparisons, $\Delta V_{bi} = 2 V_t \ln(n_{i,1}/n_{i,2})$, about 19 mV
   between the 1.45 and 1.0 values. See [`docs/PHYSICS.md` §4.2](../PHYSICS.md)
   for the engine's stance.
3. **The engine works in scaled units, mostly.** Inside `semi/physics/`
   everything is dimensionless: $\hat\psi = \psi/V_t$,
   $\hat n = n/C_0$, $\hat n_i = n_i/C_0$. The Boltzmann form in
   `semi/physics/slotboom.py` reads $\hat n = \hat n_i\,\exp(\hat\psi - \hat\Phi_n)$
   without an explicit $V_t$ in the exponent — the $V_t$ has been
   absorbed by the scaling (Ch. 12).
4. **$\Phi_n \neq E_F/q$.** Because $\Phi_n$ is defined with a sign so
   that $n$ rises when $\Phi_n$ falls below $\psi$, and because we choose
   to measure $\psi$ from the intrinsic level $E_i$ rather than from
   vacuum, the relation $E_F = -q\Phi_n$ holds only with the engine's
   sign and reference conventions. A textbook that defines the
   "imref" with the opposite sign will have flipped the relations.
5. **At equilibrium $\Phi_n$ and $\Phi_p$ are not separately determined.**
   They are only determined modulo a single global additive constant
   (the choice of where to put the Fermi level). The engine pins them
   to zero by convention (`semi/runners/bias_sweep.py:72-73`); under
   bias the ohmic contacts pin them to $V_\mathrm{applied}/V_t$ in scaled
   units (Ch. 8).

## Exercises

**Exercise 3.1.** Show that for a non-degenerate Boltzmann limit, the
electron density doubles when $E_F$ rises by $kT\ln 2 = 17.9\,\mathrm{mV}$
at 300 K.

**Exercise 3.2.** Take silicon with $N_D = 10^{17}\,\mathrm{cm^{-3}}$.
Compute $\psi_\mathrm{eq}$ on the n-side using both the asinh form
$V_t\,\mathrm{asinh}(N/2n_i)$ and the simpler logarithmic form
$V_t\,\ln(N/n_i)$. Show they agree to 4 decimal places.

**Exercise 3.3.** Estimate the Boltzmann-vs-Fermi-Dirac error in
$N_c$ for silicon at $\eta_F = -1$ (about 26 meV below the conduction
band, corresponding roughly to $n = N_c/e \approx 10^{19}\,\mathrm{cm^{-3}}$).
Use $F_{1/2}(-1) \approx 0.347$ vs $\exp(-1) = 0.368$.

**Exercise 3.4.** Verify mass action numerically at the n-bulk
($n \approx N_D = 10^{17}$, $p = n_i^2/N_D = 10^3$) and at the p-bulk
($p \approx N_A = 10^{17}$, $n = n_i^2/N_A = 10^3$) of the M1 benchmark.
Reproduce the assertion in
[`tests/check_analytical_math.py:78-80`](../../tests/check_analytical_math.py).

**Exercise 3.5.** Why does (3.7) predict $n \to \infty$ as $\Phi_n \to -\infty$
or $\psi \to \infty$? In what physical regime would this be a problem,
and how does the engine sidestep it numerically? (Hint: look up the
Slotboom underflow / overflow note in
[`docs/adr/0004-slotboom-variables-for-dd.md`](../adr/0004-slotboom-variables-for-dd.md).)

### Solutions

**3.1.** $n \propto \exp(E_F/kT)$, so $n_2/n_1 = \exp((E_{F,2} - E_{F,1})/kT)$.
Setting the ratio to 2: $E_{F,2} - E_{F,1} = kT\ln 2 = 0.02585 \cdot 0.693
= 0.01791\,\mathrm{V}$. ✓

**3.2.** asinh form: $\mathrm{asinh}(10^{23}/(2 \cdot 10^{16}))
= \mathrm{asinh}(5\times 10^6) = 16.1180$, $V_t \cdot 16.1180 = 0.4167\,\mathrm{V}$.
log form: $\ln(10^{23}/10^{16}) = \ln(10^7) = 7\ln 10 = 16.1181$,
$V_t \cdot 16.1181 = 0.4167\,\mathrm{V}$. They agree to 4 decimal places
because $\mathrm{asinh}(x) \approx \ln(2x)$ for $x \gg 1$, and the
difference $\ln 2 = 0.6931$ vs $7\ln 10 - 6\ln 10 - \ln(2) = ... $; the
correct identity is $\mathrm{asinh}(x) = \ln(x + \sqrt{x^2+1}) \approx \ln(2x)$,
giving $V_t\ln(N/n_i)$ when you do the algebra carefully. The two
forms differ by $V_t\ln 2 = 17.9\,\mathrm{mV}$ in absolute terms, but
when computing the *built-in voltage* $V_{bi} = \psi^R - \psi^L$ the
log-2 cancels (the asinhs sum to give a single ln(2)) and the difference
collapses below 0.01% (see [`tests/check_analytical_math.py:39-47`](../../tests/check_analytical_math.py)).

**3.3.** Ratio $F_{1/2}(-1)/\exp(-1) = 0.347/0.368 = 0.943$, so Boltzmann
overestimates $n$ by 5.7% at this Fermi level. At $\eta_F = 0$ ($E_F =
E_c$, severe degeneracy), $F_{1/2}(0) \approx 0.765$ vs $\exp(0) = 1$,
a 23% overestimate.

**3.4.** n-bulk: $n p = 10^{17} \cdot 10^3 = 10^{20}\,\mathrm{cm^{-6}}
= n_i^2 = (10^{10})^2$. ✓ p-bulk: same number by symmetry. ✓

**3.5.** $\Phi_n$ very negative under heavy forward bias on the n-side
overflows $\exp(\psi - \Phi_n)$. The Slotboom path keeps the *quasi-Fermi*
potentials bounded (they sit between contact biases); it is $\psi$
that is allowed to vary across $V_{bi}$. Newton iterates pull $\psi$
into its valid range continuously; the failure mode is when an over-large
Newton step pushes an iterate outside this range and the next residual
evaluation overflows. ADR 0004 describes the bias-continuation cure: ramp
slowly enough that no Newton step exceeds the radius of quadratic
convergence.

## Further reading

- **Sze and Ng, *Physics of Semiconductor Devices*, 3rd ed. (2007).**
  §1.4: Fermi–Dirac, density of states, intrinsic concentration. The
  effective DOS values in [`semi/materials.py`](../../semi/materials.py)
  come from Table 7 of this chapter.
- **Pierret, *Semiconductor Device Fundamentals* (1996).** §2.5–2.6:
  Boltzmann limit and validity boundary, with worked examples.
- **Altermatt, P. P., et al. (2003).** "A simulation model for the
  density of states and for incomplete ionization in crystalline silicon."
  *J. Appl. Phys.* 93, 1598. The source of the modern $n_i = 1.0\times 10^{10}$
  value used by kronos-semi (vs the older textbook 1.45×10¹⁰).
- **Hu, *Modern Semiconductor Devices for Integrated Circuits* (2010).**
  §1.7: Fermi level, $n_i$ at 300 K, intuitive numerical examples.
