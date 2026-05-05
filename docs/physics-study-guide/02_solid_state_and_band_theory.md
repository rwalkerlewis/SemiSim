# 2 — Solid state and band theory

## Learning objectives

- Sketch the Bloch picture: how the periodicity of a crystal lattice
  forces electron states into bands and gaps.
- Explain the meaning of $E_c$, $E_v$, $E_g$, and the electron affinity $\chi$.
- State the parabolic-band / effective-mass approximation and write down
  the resulting expressions for the effective densities of states $N_c$
  and $N_v$.
- Distinguish direct- from indirect-gap semiconductors and connect the
  distinction to which recombination mechanisms (SRH vs Auger vs
  radiative) dominate.
- Recognize $\chi$ and $\phi_{ms}$ in `semi/materials.py` and `semi/bcs.py`,
  and predict where they will reappear in MOS physics (Ch. 9) and
  heterojunctions (forward reference to M17).

## Physical motivation

Why do some materials conduct (metals), some insulate (SiO₂, diamond),
and a third class do "almost neither" but with electrically tuneable
behaviour (Si, Ge, GaAs)? The answer requires quantum mechanics applied
to a crystal lattice. A free electron has a continuous spectrum of
allowed energies $E = \hbar^2 k^2 / 2m_0$. An electron in a crystal sees
a periodic potential $V(\mathbf{r}+\mathbf{a}_i) = V(\mathbf{r})$ for the
lattice vectors $\mathbf{a}_i$, and the resulting eigenstates have
allowed energies organized into **bands** separated by **gaps** —
forbidden energy ranges where no propagating eigenstate exists. Whether
the highest occupied band (the **valence band**) is full and the next
(the **conduction band**) is empty determines whether the material is an
insulator, a metal, or a semiconductor.

In silicon at room temperature, the valence band is essentially full
and the conduction band is essentially empty, separated by an energy
gap $E_g = 1.12\,\mathrm{eV}$. Promoting an electron across the gap
(thermally, optically, or by injection from a contact) leaves a
positively charged hole behind in the valence band. Both the electron
and the hole then act as mobile carriers, with effective masses
$m_n^*$ and $m_p^*$ that capture the curvature of the band dispersion
near its extremum. The drift-diffusion machinery from Ch. 5 onward
treats $n$, $p$, $E_c$, $E_v$, $\chi$, and $E_g$ as bulk material
properties whose origin is band structure but whose use is purely
classical.

## Derivation from first principles

### Bloch's theorem (sketch)

Write the single-electron Schrödinger equation in a periodic potential
$V(\mathbf{r})$:

$$
\left[-\,\frac{\hbar^2}{2m_0}\,\nabla^2 + V(\mathbf{r})\right]\psi(\mathbf{r})
   = E\,\psi(\mathbf{r}),
\qquad V(\mathbf{r}+\mathbf{R}) = V(\mathbf{r}).
\tag{2.1}
$$

Bloch's theorem says the eigenfunctions take the form

$$
\psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}}\,u_{n\mathbf{k}}(\mathbf{r}),
\qquad u_{n\mathbf{k}}(\mathbf{r}+\mathbf{R}) = u_{n\mathbf{k}}(\mathbf{r}),
\tag{2.2}
$$

where $n$ is a band index and $\mathbf{k}$ ranges over the first
Brillouin zone (the unit cell of the reciprocal lattice). Each band has
its own dispersion $E_n(\mathbf{k})$, periodic in $\mathbf{k}$. Across a
band, $\mathbf{k}$ varies continuously and so does $E$; *between*
bands, gaps open up — there exist energy ranges $E$ for which no
$\mathbf{k}$ in the Brillouin zone satisfies $E_n(\mathbf{k}) = E$. We
will not derive (2.2) here; see Ashcroft & Mermin §8 or Kittel §7 for the
full derivation.

### Effective mass and parabolic bands

Near the extremum of a band — the maximum of the valence band or the
minimum of the conduction band — Taylor-expand $E(\mathbf{k})$ to second
order. For an isotropic effective mass $m^*$:

$$
E(\mathbf{k}) \;\approx\; E_0 + \frac{\hbar^2 |\mathbf{k} - \mathbf{k}_0|^2}{2m^*}.
\tag{2.3}
$$

This is the **effective-mass approximation**: an electron near the
conduction-band minimum behaves like a free particle with mass $m_n^*$,
and a hole near the valence-band maximum behaves like a free particle
with mass $m_p^*$. The drift-diffusion equations are written for these
quasi-classical particles.

### Density of states from a parabolic band

Counting the number of $\mathbf{k}$-states per unit volume between $E$
and $E + dE$ gives, for a 3D parabolic band,

$$
g(E) = \frac{1}{2\pi^2}\left(\frac{2m^*}{\hbar^2}\right)^{3/2} \sqrt{E - E_0},
\qquad E \geq E_0.
\tag{2.4}
$$

Multiplying by 2 for spin and integrating against the Fermi–Dirac
distribution gives the carrier densities (Ch. 3). When the Fermi level
sits well below the conduction-band minimum (or well above the
valence-band maximum) — the **non-degenerate** regime — the
Fermi–Dirac integral collapses to a Boltzmann factor and the result is

$$
n = N_c\,\exp\!\left(-\frac{E_c - E_F}{kT}\right),
\qquad
p = N_v\,\exp\!\left(-\frac{E_F - E_v}{kT}\right),
\tag{2.5}
$$

with the **effective densities of states**

$$
\boxed{
N_c = 2\left(\frac{m_n^* kT}{2\pi\hbar^2}\right)^{3/2},
\qquad
N_v = 2\left(\frac{m_p^* kT}{2\pi\hbar^2}\right)^{3/2}.
}
\tag{2.6}
$$

The full algebra of (2.6) is in Appendix B §B.1; the result is what
matters for kronos-semi. At $T = 300\,\mathrm{K}$ in silicon,
$N_c = 2.86\times 10^{19}\,\mathrm{cm}^{-3}$ and
$N_v = 3.10\times 10^{19}\,\mathrm{cm}^{-3}$; these values come from
Sze 3rd ed. Table 7 and appear in [`semi/materials.py:56-57`](../../semi/materials.py).

### Direct vs indirect gap

In a **direct-gap** material (GaAs, InGaAs), the valence-band maximum
and the conduction-band minimum sit at the same $\mathbf{k}$ in the
Brillouin zone. An electron–hole recombination event can emit a photon
of energy $\sim E_g$ without violating momentum conservation — the
photon carries effectively zero momentum. GaAs LEDs and laser diodes
exploit this. Radiative recombination is therefore a first-order
process in direct-gap materials.

In an **indirect-gap** material (Si, Ge), the valence-band maximum and
conduction-band minimum sit at *different* $\mathbf{k}$. A direct
photon-emission process would violate momentum conservation; a phonon
must be involved as an intermediate momentum reservoir, which makes
radiative recombination a second-order process and very slow. In
silicon, recombination is dominated by Shockley–Read–Hall (Ch. 6) and,
at high carrier densities, Auger (M16.3). The shipped engine treats only
SRH; Auger and radiative are M16.3 (planned).

### Electron affinity $\chi$ and work function

The electron affinity $\chi$ is the energy required to remove an
electron from the conduction-band minimum to the vacuum level. For
silicon, $\chi = 4.05\,\mathrm{eV}$ (`semi/materials.py:55`). The Fermi
level $E_F$ sits inside the gap, so the work function — the energy
required to remove an electron from the Fermi level to vacuum — is

$$
\Phi_s = \chi + (E_c - E_F)
       = \chi + V_t\,\ln(N_c / n),
\tag{2.7}
$$

which depends on the doping. For a metal, the work function $\Phi_m$
is a single number characterizing the metal; the **metal–semiconductor
work-function difference** $\phi_{ms} = \Phi_m - \Phi_s$ shows up in
the gate boundary condition (1.7) and is the principal driver of MOS
flat-band voltage shifts (Ch. 9).

## Key results

$$
n = N_c\,\exp\!\left(-\frac{E_c - E_F}{kT}\right),
\qquad
p = N_v\,\exp\!\left(-\frac{E_F - E_v}{kT}\right)
\tag{2.8}
$$

$$
n_i^2 = N_c N_v\,\exp\!\left(-\frac{E_g}{kT}\right)
\tag{2.9}
$$

$$
\phi_{ms} = \Phi_m - \Phi_s = \Phi_m - \big[\chi + V_t\,\ln(N_c/n)\big]
\tag{2.10}
$$

Equation (2.9) follows from multiplying $n$ and $p$ from (2.8), using
$E_c - E_v = E_g$, and dropping the dependence on $E_F$ that exactly
cancels — this is the **mass-action law** that Ch. 3 puts to numerical
use.

## Worked numerical example

For silicon at $T = 300\,\mathrm{K}$, plug Sze's values into (2.6):

- $m_n^* = 1.08\,m_0$ (density-of-states effective mass for the six
  conduction-band valleys), $kT = 0.02585\,\mathrm{eV} =
  4.14\times 10^{-21}\,\mathrm{J}$.
- $\hbar = 1.0546\times 10^{-34}\,\mathrm{J\,s}$,
  $m_0 = 9.109\times 10^{-31}\,\mathrm{kg}$.
- $\frac{m_n^* kT}{2\pi\hbar^2} = \frac{1.08 \cdot 9.109\times 10^{-31} \cdot 4.14\times 10^{-21}}{2\pi\cdot 1.112\times 10^{-68}}
   = \frac{4.07\times 10^{-51}}{6.99\times 10^{-68}} = 5.83\times 10^{16}\,\mathrm{m^{-2}}$.
- $N_c = 2 \cdot (5.83\times 10^{16})^{3/2} = 2 \cdot 1.41\times 10^{25}
   = 2.82\times 10^{25}\,\mathrm{m^{-3}}$.
- Convert to cm⁻³: divide by $10^6$, get $2.82\times 10^{19}\,\mathrm{cm}^{-3}$.

This matches Sze's tabulated $N_c = 2.86\times 10^{19}\,\mathrm{cm}^{-3}$
to better than 2%. The remaining difference comes from the slightly
different $m_n^*$ values quoted in different references; the engine uses
Sze's value directly.

For (2.9) with $E_g = 1.12\,\mathrm{eV}$ at 300 K:
$\exp(-E_g/kT) = \exp(-1.12/0.02585) = \exp(-43.32) = 1.55\times 10^{-19}$.
Then $n_i^2 = 2.86\times 10^{19} \cdot 3.10\times 10^{19} \cdot 1.55\times 10^{-19}
= 1.37\times 10^{20}\,\mathrm{cm^{-6}}$, so $n_i \approx 1.17\times 10^{10}\,\mathrm{cm}^{-3}$.

The engine uses $n_i = 1.0\times 10^{10}\,\mathrm{cm}^{-3}$ (Altermatt's
2003 reassessment; see [`semi/materials.py:58`](../../semi/materials.py)
and the note in [`docs/PHYSICS.md` §4.2](../PHYSICS.md)). The
discrepancy of order $V_t\ln(1.17/1.0)\approx 4\,\mathrm{mV}$ on a
built-in voltage is well within practical tolerance and matches the
modern TCAD convention.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| Bandgap $E_g$ | – | `semi/materials.py:31` (`Material.Eg`) |
| Electron affinity $\chi$ | (2.7) | `semi/materials.py:32` (`Material.chi`) |
| Effective DOS $N_c$, $N_v$ | (2.6) | `semi/materials.py:31-32` (`Material.Nc`, `Nv`); values at lines 56-57 (Si), 68-69 (Ge), 80-81 (GaAs) |
| Boltzmann carrier expressions | (2.8) | `semi/physics/slotboom.py:27-42` (uses $\psi - \Phi_n$ and $\Phi_p - \psi$ in scaled units; see Ch. 11) |
| Mass action $n_i^2$ | (2.9) | `tests/check_analytical_math.py:78-80` |
| Work-function offset $\phi_{ms}$ | (2.10) | `semi/bcs.py:186-188` (`workfunction` field) |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §4.2](../PHYSICS.md) — material parameter table for Si.
- [`docs/PHYSICS.md` §4.3](../PHYSICS.md) — Ge, GaAs, and the insulators.
- [`docs/PHYSICS_INTRO.md` §3](../PHYSICS_INTRO.md) — non-rigorous overview of bands.
- [`docs/talk/02-physics.md`](../talk/02-physics.md) — talk-style overview (no derivations).

## Common pitfalls

1. **Effective mass is not unique.** There are at least three commonly
   quoted "electron effective masses" in silicon: the longitudinal
   $m_l$, the transverse $m_t$, and the density-of-states-effective
   $m_n^*$. The latter is what enters (2.6). The number Sze tabulates,
   and the number behind kronos-semi's $N_c$, is $m_n^* = 1.08\,m_0$,
   not the conductivity effective mass that enters mobility formulas.
2. **Boltzmann limit only.** Equation (2.8) assumes $E_F$ sits well
   below $E_c$ (and well above $E_v$). At doping above $\sim 10^{19}\,\mathrm{cm}^{-3}$,
   $E_F$ approaches $E_c$ and (2.8) overestimates $n$ by roughly a
   factor of 2; the correct expression uses the Fermi–Dirac integral
   $F_{1/2}$. This is the boundary that motivates **M16.4** (planned).
3. **$E_F$ is not constant under bias.** At thermal equilibrium one
   Fermi level is well-defined throughout the device. Under bias the
   electron and hole populations are out of equilibrium with each other
   and you need *two* quasi-Fermi levels $E_{F,n}$ and $E_{F,p}$, which
   become the engine's primary unknowns $\Phi_n$ and $\Phi_p$ once you
   factor out the sign and the unit conversion (Ch. 11).
4. **$E_g$ is temperature-dependent.** $E_g(T)$ in silicon decreases by
   about 0.4 meV/K above 100 K; at 300 K, $E_g = 1.12\,\mathrm{eV}$, but
   at 77 K it is about 1.16 eV. kronos-semi uses the 300 K value
   throughout; cryogenic device simulation would need a $T$-dependent
   $E_g$, which is M16+ work.
5. **Direct vs indirect distinction is invisible to the DD engine.** The
   drift-diffusion equations don't know whether a transition is direct
   or indirect; they only know the SRH and (planned) Auger / radiative
   rate constants. The distinction shows up in the *values* of those
   constants — silicon's radiative coefficient is $\sim 10^{-14}\,\mathrm{cm^3/s}$,
   GaAs's is $\sim 10^{-10}\,\mathrm{cm^3/s}$ — not in the equations.

## Exercises

**Exercise 2.1.** Repeat the $N_c$ calculation for germanium using the
Sze table value $m_n^* = 0.56\,m_0$ at 300 K. Compare with `semi/materials.py:68`.

**Exercise 2.2.** Show that $n_i$ in (2.9) is independent of which Fermi
level $E_F$ you use, despite (2.8) appearing to depend on it.

**Exercise 2.3.** A metal has work function $\Phi_m = 4.7\,\mathrm{eV}$
(close to typical poly-silicon). The semiconductor underneath is p-type
silicon with $N_a = 5\times 10^{16}\,\mathrm{cm^{-3}}$. Compute
$\phi_{ms}$. Hint: the bulk Fermi potential
$|\phi_B| = V_t\ln(N_a/n_i) = 0.399\,\mathrm{V}$ (this is the same number
that drives the MOSCAP analytic helpers in `semi/cv.py`); express $\Phi_s$
in terms of $\chi$, $E_g$, and $|\phi_B|$.

**Exercise 2.4.** Why does the engine treat insulators (SiO₂, HfO₂,
Si₃N₄) as having $E_g = 0$ and $\chi = 0$ in `semi/materials.py`?
Reconcile this with the fact that SiO₂ has a real bandgap of about
9 eV.

**Exercise 2.5.** GaAs has $E_g = 1.424\,\mathrm{eV}$ and a small $N_c$
($4.7\times 10^{17}\,\mathrm{cm}^{-3}$, see `semi/materials.py:80`).
Predict whether $n_i$ in GaAs is larger or smaller than in silicon at
300 K, and by approximately what factor.

### Solutions

**2.1.** $\frac{m_n^* kT}{2\pi\hbar^2} = \frac{0.56 \cdot 9.109\times 10^{-31} \cdot 4.14\times 10^{-21}}{2\pi\cdot 1.112\times 10^{-68}}
= 3.02\times 10^{16}\,\mathrm{m^{-2}}$. Then $N_c = 2(3.02\times 10^{16})^{1.5} = 2 \cdot 5.25\times 10^{24} = 1.05\times 10^{25}\,\mathrm{m^{-3}}
= 1.05\times 10^{19}\,\mathrm{cm^{-3}}$. Sze gives $1.04\times 10^{19}$; matches.

**2.2.** Multiply $n$ and $p$ in (2.8): the $\exp(\pm E_F/kT)$ factors
cancel, leaving $np = N_c N_v \exp(-(E_c-E_v)/kT) = N_c N_v \exp(-E_g/kT)$,
which is the right side of (2.9). The $E_F$-cancellation is the
mass-action law's microscopic origin: $n_i$ is an intrinsic property of
the band structure, independent of where the Fermi level sits.

**2.3.** $\Phi_s = \chi + (E_c - E_F)$. For p-type at the equilibrium
$E_F = E_v + |\phi_B| = E_c - E_g + |\phi_B|$ (positive $|\phi_B|$
distance below mid-gap on the p-side; conventions vary, but the
combination $E_c - E_F = E_g/2 + |\phi_B|$ for a non-degenerate p-body).
With $E_g = 1.12$, $|\phi_B| = 0.399$, $\chi = 4.05$:
$\Phi_s = 4.05 + 1.12/2 + 0.399 = 5.01\,\mathrm{eV}$.
Then $\phi_{ms} = 4.7 - 5.01 = -0.31\,\mathrm{V}$.
The MOSCAP benchmark uses $\phi_{ms} = -0.95\,\mathrm{V}$; the
discrepancy is because that benchmark uses the n+ poly-silicon gate
work function ($\Phi_m \approx \chi_\mathrm{Si} = 4.05$ eV, since the
Fermi level of n+ poly is near $E_c$), not 4.7 V. Recompute:
$\phi_{ms} = 4.05 - 5.01 = -0.96\,\mathrm{V}$, matching `semi/cv.py:69`.

**2.4.** The engine uses the SiO₂ region only as a Laplacian dielectric:
no carriers, no doping, no bandgap-driven recombination. The Poisson
equation in the oxide reduces to $-\nabla\!\cdot\!(\varepsilon_r\nabla\psi) = 0$
which depends only on $\varepsilon_r$. The 9 eV bandgap of SiO₂ matters
for tunneling (M16.6, planned) and breakdown, neither of which the
shipped engine models, so $E_g$ and $\chi$ for SiO₂ are not loaded into
any UFL form and are stored as zero by convention.

**2.5.** Compare $n_i$ via (2.9). For silicon: $n_i^2 = N_c N_v \exp(-E_g/kT)
= 2.86\times 10^{19} \cdot 3.10\times 10^{19} \cdot \exp(-1.12/0.02585)
= 8.87\times 10^{38} \cdot 1.55\times 10^{-19} = 1.37\times 10^{20}$, so
$n_i \approx 1.2\times 10^{10}\,\mathrm{cm^{-3}}$. For GaAs:
$n_i^2 = 4.7\times 10^{17} \cdot 7.0\times 10^{18} \cdot \exp(-1.424/0.02585)
= 3.29\times 10^{36} \cdot \exp(-55.08) = 3.29\times 10^{36} \cdot 1.16\times 10^{-24}
= 3.82\times 10^{12}$, so $n_i \approx 1.95\times 10^{6}\,\mathrm{cm^{-3}}$.
GaAs $n_i$ is about $10^4$ times smaller; this is a generic feature of
wide-gap materials and explains why GaAs devices are insensitive to
the SRH-generation reverse-leakage that plagues silicon (Ch. 6).

## Further reading

- **Ashcroft and Mermin, *Solid State Physics* (1976).** Chapter 8 (Bloch
  theorem), Chapter 12 (semiconductor band structure). Standard
  graduate text; rigorous and exhaustive.
- **Kittel, *Introduction to Solid State Physics*, 8th ed. (2004).**
  Chapters 7–8 are at the right level for this guide's audience.
- **Sze and Ng, *Physics of Semiconductor Devices*, 3rd ed. (2007).**
  Chapter 1 §§1.2–1.4: bands, density of states, intrinsic carriers, with
  the parameter tables that kronos-semi's material database draws from.
- **Pierret, *Advanced Semiconductor Fundamentals*, 2nd ed. (2002).**
  Chapter 4 for a more device-physics-oriented treatment.
