# 4 — Doping and charge neutrality

## Learning objectives

- State the donor / acceptor / complete-ionization model and recognize
  its assumptions in the engine's net-doping callable.
- Derive the bulk equilibrium potential $\psi_\mathrm{eq} =
  V_t\,\mathrm{asinh}(N_\mathrm{net}/(2n_i))$ from charge neutrality.
- Identify the three doping-profile types in the schema (uniform, step,
  Gaussian) and connect each to a physical fabrication step.
- Predict the sign of $\psi_\mathrm{eq}$ for n-type and p-type bulk and
  unit-check the engine's contact BCs against that prediction.
- Locate every doping evaluator in `semi/doping.py` and trace a
  benchmark JSON's `doping[]` array to the function it builds.

## Physical motivation

A pure (intrinsic) semiconductor at room temperature has very few
mobile carriers — about $10^{10}\,\mathrm{cm^{-3}}$ in silicon, ten
orders of magnitude below the atom density. To make a useful device
you replace a controlled fraction of host-lattice atoms with **donors**
(group-V atoms in silicon: P, As, Sb, contributing one electron each)
or **acceptors** (group-III atoms: B, Al, contributing one hole each).
Doping at $10^{17}\,\mathrm{cm^{-3}}$ is one impurity per $5\times 10^5$
silicon atoms — still very dilute on the chemistry scale, but enough
to change the carrier density by seven orders of magnitude.

The shape of the doping profile $N_D(\mathbf{x}), N_A(\mathbf{x})$
encodes the device's electrical structure. A pn junction is a step in
$N_D - N_A$ at a plane (Ch. 7). A MOSFET source/drain is a Gaussian
implant offset from the surface (Ch. 10). A homogeneous resistor is a
uniform $N_D$ everywhere (M7 benchmark). The engine accepts these as
JSON specifications and builds a callable `N_net(x)` that the rest of
the pipeline treats as a known coefficient field.

## Derivation from first principles

### Donor and acceptor energy levels

A substitutional phosphorus atom in silicon has one extra valence
electron compared with silicon's four. That electron is loosely bound
to the P⁺ core by hydrogenic coupling, with a binding energy of about
45 meV in silicon. At 300 K, $kT = 25.85\,\mathrm{meV}$, so
$\exp(-45/25.85) = 0.176$, meaning roughly 95% of donor electrons are
ionized into the conduction band at room temperature for a low doping
density. As doping increases, donor states broaden into a band that
overlaps the conduction band, and ionization approaches 100% — this is
the **complete ionization** assumption.

Symmetrically, a substitutional boron atom captures an extra electron
from the valence band, leaving behind a free hole and a B⁻ acceptor
core. The acceptor binding energy in silicon is about 45 meV; complete
ionization holds for the same reason.

kronos-semi assumes complete ionization at 300 K throughout. Freeze-out
effects below ~100 K are out of scope; cryogenic device simulation
would need to relax this.

### Charge neutrality in bulk

Far from any junction or contact, the electrostatic potential is flat
and the local volume is electrically neutral:

$$
p(\psi) - n(\psi) + N_D - N_A = 0,
\tag{4.1}
$$

with $n, p$ from Boltzmann (3.7) at $\Phi_n = \Phi_p = 0$:

$$
n_i\,e^{-\psi/V_t} - n_i\,e^{\psi/V_t} + (N_D - N_A) = 0.
$$

Using $\sinh x = (e^x - e^{-x})/2$, this rearranges to

$$
2\,n_i\,\sinh(\psi/V_t) = N_D - N_A,
$$

so

$$
\boxed{
\psi_\mathrm{eq} = V_t\,\mathrm{asinh}\!\left(\frac{N_D - N_A}{2\,n_i}\right).
}
\tag{4.2}
$$

This is the **equilibrium potential at local charge neutrality**.
Equation (4.2) is what `semi/physics/slotboom.py:73-83` (`equilibrium_psi_hat`)
returns and what the contact BC builder pulls in (1.6).

For $N_\mathrm{net} \gg n_i$ (heavily doped),
$\mathrm{asinh}(x) \approx \ln(2x)$, so

$$
\psi_\mathrm{eq} \approx V_t\,\ln\!\left(\frac{N_\mathrm{net}}{n_i}\right)
\quad (N_\mathrm{net} \gg n_i),
\tag{4.3}
$$

with sign matching $N_\mathrm{net}$. This is the textbook "log form"
that produces $V_{bi} = V_t\,\ln(N_A N_D / n_i^2)$ for a pn junction
when you take $\psi_\mathrm{eq}^R - \psi_\mathrm{eq}^L$.

### The three profile types in the schema

`semi/doping.py` ships three callable builders:

**Uniform.** $N_\mathrm{net}(\mathbf{x}) = N_D - N_A$ everywhere
in the region. Physical interpretation: uniform background doping over
a substrate. Used by the M7 resistor benchmark.

**Step.** $N_\mathrm{net}(\mathbf{x}) = N_\mathrm{left}$ for
$x_\mathrm{axis} < x_\mathrm{loc}$ and $N_\mathrm{right}$ otherwise.
Physical interpretation: an abrupt junction (idealized; real diffusion
profiles have finite width). Used by the M1 pn benchmark.

**Gaussian.** $N_\mathrm{net}(\mathbf{x}) = N_\mathrm{bg} \pm
N_\mathrm{peak}\,\exp(-\tfrac12 r^2)$ with $r^2 = \sum_d
((x_d - c_d)/\sigma_d)^2$. Physical interpretation: an ion implant
followed by drive-in diffusion. The plus sign is for donors, minus for
acceptors. Used by the M12 MOSFET benchmark for n+ source/drain
implants.

The three are summed when multiple `doping[]` entries cover the same
region — that lets you superpose a Gaussian implant on a uniform
background ([`semi/doping.py:23-56`](../../semi/doping.py)).

### Sign convention

The engine works with the **signed net doping** $N_\mathrm{net} = N_D - N_A$:
positive on n-side, negative on p-side. The Boltzmann right-hand side
of Poisson is $p - n + N_\mathrm{net}$, with the same sign convention.
This single signed number replaces the two separate $N_D, N_A$ inputs
inside `build_dd_block_residual`; the JSON exposes both because user
intent (donor implant vs acceptor implant) is clearer that way.

## Key results

- Net doping evaluator: `build_profile(doping_list) -> Callable`
  ([`semi/doping.py:23-56`](../../semi/doping.py)).
- Bulk equilibrium potential: (4.2).
- Heavy-doping log form: (4.3).
- Built-in voltage of an asymmetric pn junction (preview of Ch. 7):
  $V_{bi} = V_t\,\ln(N_A N_D / n_i^2)$.

## Worked numerical example

The M1 benchmark uses a step profile:
- left region ($x < 1\,\mu\mathrm{m}$): $N_A = 10^{17}\,\mathrm{cm^{-3}}$, p-type
- right region ($x > 1\,\mu\mathrm{m}$): $N_D = 10^{17}\,\mathrm{cm^{-3}}$, n-type

In SI: $N_A = 10^{23}\,\mathrm{m^{-3}}$, $N_D = 10^{23}\,\mathrm{m^{-3}}$,
$n_i = 10^{16}\,\mathrm{m^{-3}}$, $V_t = 0.02585\,\mathrm{V}$.

Left bulk: $N_\mathrm{net} = -10^{23}$ (p-type).
$\psi_L = V_t\,\mathrm{asinh}(-10^{23}/(2\cdot 10^{16})) = -V_t\,\mathrm{asinh}(5\times 10^6) = -V_t \cdot 16.118 = -0.4167\,\mathrm{V}$.

Right bulk: $N_\mathrm{net} = +10^{23}$.
$\psi_R = +V_t\,\mathrm{asinh}(5\times 10^6) = +0.4167\,\mathrm{V}$.

Built-in voltage: $V_{bi} = \psi_R - \psi_L = 0.8334\,\mathrm{V}$.

Cross-check via the log form:
$V_t\,\ln(N_A N_D/n_i^2) = 0.02585 \cdot \ln(10^{46}/10^{32}) = 0.02585 \cdot \ln(10^{14}) = 0.02585 \cdot 32.236 = 0.8333\,\mathrm{V}$. ✓

The 0.0001 V mismatch is the asinh/log algebraic remainder, well below
SNES tolerance and irrelevant for any device-level number. This same
match is asserted in `tests/check_analytical_math.py:46-47`.

Mass action at the n-bulk:
$n^R \approx N_D = 10^{17}\,\mathrm{cm^{-3}}$,
$p^R = n_i \exp(-\psi_R/V_t) = 10^{10}\cdot e^{-16.118} = 10^{10} \cdot 10^{-7} = 10^3\,\mathrm{cm^{-3}}$.
$n^R \cdot p^R = 10^{20}\,\mathrm{cm^{-6}} = n_i^2$. ✓
This is the assertion at [`tests/check_analytical_math.py:78-80`](../../tests/check_analytical_math.py).

## Code map

| Concept | Equation | Code location |
|---|---|---|
| `build_profile(doping_list)` | sums profiles, returns callable | `semi/doping.py:23-56` |
| Uniform profile | $N = N_D - N_A$ | `semi/doping.py:70-76` |
| Step profile | left/right step | `semi/doping.py:79-89` |
| Gaussian implant | $N_\mathrm{bg} \pm N_\mathrm{pk}\,e^{-r^2/2}$ | `semi/doping.py:92-107` |
| Equilibrium $\psi$ | (4.2) | `semi/physics/slotboom.py:73-83` |
| Charge neutrality test | (4.1) | `tests/check_analytical_math.py:82-93` |
| Net doping interpolation onto a Function | $\hat N = N/C_0$ | `semi/runners/equilibrium.py:36-40`, `bias_sweep.py:67-68` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §1.1](../PHYSICS.md) — Poisson with $N_D - N_A$.
- [`docs/PHYSICS.md` §3.1](../PHYSICS.md) — ohmic-contact equilibrium $\psi$ with the asinh form.
- [`docs/schema/reference.md` §doping](../schema/reference.md) — JSON contract for the three profile types.

## Common pitfalls

1. **Density units.** Doping in JSON is in cm⁻³ by device-physics
   convention; the engine converts to m⁻³ at ingest via `cm3_to_m3`
   ([`semi/constants.py:29-31`](../../semi/constants.py)). All three
   profile builders call this helper. If you bypass them and inject a
   raw 10¹⁷ into a UFL form, you have just doped the device $10^6$ times
   too lightly.
2. **Sign of "donor" in the Gaussian profile.** The `dopant: "donor"` /
   `"acceptor"` field on a Gaussian entry chooses the sign:
   donor → +peak (adds to n-type), acceptor → −peak (subtracts from
   net doping, making it more p-type) ([`semi/doping.py:96`](../../semi/doping.py)).
   Forgetting this when building a multi-region MOSFET implant on a
   p-type body inverts the device.
3. **Complete ionization is an assumption.** At doping above
   $\sim 10^{19}\,\mathrm{cm^{-3}}$, dopant–dopant interaction broadens
   the impurity band into the conduction band and the assumption is
   accurate; below $\sim 10^{15}$ at room temperature the assumption is
   fine because $kT$ is well above the binding energy. The risky regime
   is heavy doping at low $T$ (cryogenic), where freeze-out reduces the
   ionized fraction. kronos-semi targets 300 K and does not model this.
4. **`region` field is informational.** As of the current ingest, the
   doping callable is evaluated at every mesh point regardless of the
   `region` key in the JSON entry; the field is read but does not gate
   the evaluation ([`semi/doping.py:42`](../../semi/doping.py)). If
   region-restricted evaluation matters (heterojunctions, M17), it will
   need to be added.
5. **Multiple entries sum.** If you have a uniform background plus a
   Gaussian implant for a MOSFET source, the JSON has two entries; the
   `build_profile` call sums them ([`semi/doping.py:50-54`](../../semi/doping.py)).
   This is *additive*, not *override*; if you accidentally list the
   same uniform doping twice you have doubled the doping.

## Exercises

**Exercise 4.1.** Write down the explicit form of $\psi_\mathrm{eq}$ for
a p-type body with $N_A = 5\times 10^{16}\,\mathrm{cm^{-3}}$ (the M14.2
MOSCAP body). Compare with the $|\phi_B| = 0.399\,\mathrm{V}$ asserted
by `analytical_moscap_params` in [`semi/cv.py:123`](../../semi/cv.py).

**Exercise 4.2.** Show that for the symmetric pn junction (M1 benchmark),
$x_n = x_p = W/2$ where $W$ is the depletion width. Use charge balance
$N_D x_n = N_A x_p$.

**Exercise 4.3.** Construct the Gaussian-implant parameters for an n+
source contact in a 2D MOSFET: peak doping $10^{20}\,\mathrm{cm^{-3}}$,
center at $(x, y) = (0.25\,\mu\mathrm{m}, 0)$, $\sigma_x = 0.1\,\mu\mathrm{m}$,
$\sigma_y = 0.05\,\mu\mathrm{m}$, on a p-type background of
$10^{17}\,\mathrm{cm^{-3}}$. What is the net doping at the implant
center? What about 200 nm laterally and 100 nm deep?

**Exercise 4.4.** Verify the assertion at
[`tests/check_analytical_math.py:64-65`](../../tests/check_analytical_math.py):
for the symmetric junction, $x_p N_A = x_n N_D$ at machine precision.

**Exercise 4.5.** A p-type body is gradually doped from
$10^{17}\,\mathrm{cm^{-3}}$ at the surface to $10^{15}$ at 1 µm depth
(linearly in $\log_{10} N_A$). Sketch $\psi_\mathrm{eq}(z)$ on this
profile. Why does kronos-semi not yet support arbitrary table-based
doping profiles?

### Solutions

**4.1.** $N_\mathrm{net} = -5\times 10^{16}\,\mathrm{cm^{-3}}
= -5\times 10^{22}\,\mathrm{m^{-3}}$.
$\psi_\mathrm{eq} = V_t\,\mathrm{asinh}(-5\times 10^{22}/(2\cdot 10^{16}))
= -V_t\,\mathrm{asinh}(2.5\times 10^6)
= -V_t \cdot \ln(5\times 10^6)
= -0.02585 \cdot 15.425 = -0.3989\,\mathrm{V}$. The MOSCAP code reports
$|\phi_B| = 0.399\,\mathrm{V}$ — same number to three decimal places.
The sign convention shifts: the engine reports $|\phi_B|$ as a positive
magnitude; the equilibrium $\psi$ on a p-body is $-|\phi_B|$.

**4.2.** Charge balance in the depletion approximation: the total
positive charge on the n-side (ionized donors $N_D$ over width $x_n$)
equals the total negative charge on the p-side (ionized acceptors $N_A$
over $x_p$): $N_D x_n = N_A x_p$. With $N_D = N_A$,
$x_n = x_p = W/2$. ✓

**4.3.** At the center, both Gaussian arguments are zero:
$N_\mathrm{net} = -10^{17} + 10^{20} = +9.99\times 10^{19}\,\mathrm{cm^{-3}}$
(strongly n-type).
At $\Delta x = 200\,\mathrm{nm}$ (i.e. $r_x^2 = (0.2/0.1)^2 = 4$, $r_y = 0$):
$\exp(-r^2/2) = \exp(-2) = 0.135$, so the Gaussian contributes
$10^{20}\cdot 0.135 = 1.35\times 10^{19}$;
$N_\mathrm{net} = -10^{17} + 1.35\times 10^{19} = +1.34\times 10^{19}$
(still n-type, but less heavily).
At $\Delta y = 100\,\mathrm{nm}$ deep ($r_y^2 = (0.1/0.05)^2 = 4$):
$\exp(-r^2/2) = e^{-2} = 0.135$, same magnitude as the lateral case.

**4.4.** Symmetric junction: $N_A = N_D$, so $x_p = x_n$ from charge
balance. The test computes $x_p \cdot N_A$ and $x_n \cdot N_D$ from the
M1 device parameters and asserts the relative difference is below
$10^{-6}$. With $N_A = N_D$, this collapses to $x_p = x_n$, which is
true to roundoff from the symmetric initial guess.

**4.5.** $\psi_\mathrm{eq}(z) = V_t\,\mathrm{asinh}(-N_A(z)/(2n_i))$
varies from $-0.4167\,\mathrm{V}$ at $z = 0$ to $-0.2980\,\mathrm{V}$
at $z = 1\,\mu\mathrm{m}$ (where $N_A = 10^{15}$:
$V_t\,\mathrm{asinh}(-5\times 10^4) \approx -V_t \cdot 11.513 = -0.298\,\mathrm{V}$).
The graded profile is not yet supported because the schema does not
expose a `table` profile type ([`semi/doping.py:67`](../../semi/doping.py)
explicitly raises). Adding it is a small change but has not been
prioritized; the Gaussian profile covers most realistic implant cases.

## Further reading

- **Sze and Ng (2007), §1.5** for donor/acceptor energy levels and the
  ionization equation in detail.
- **Pierret (1996), §3.2–3.3** for the bulk equilibrium picture and
  charge neutrality with worked examples.
- **Selberherr (1984), §1.4** for the device-physics framing of doping
  profiles, including diffusion-driven Gaussian shapes — the same model
  kronos-semi uses for source/drain implants.
- **Plummer, Deal, and Griffin, *Silicon VLSI Technology* (2000),
  §7–8** for the fabrication context: how Gaussian (or
  Pearson-shape-corrected Gaussian) profiles arise from real ion
  implantation followed by drive-in.
