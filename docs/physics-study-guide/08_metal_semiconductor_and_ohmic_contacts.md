# 8 — Metal–semiconductor and ohmic contacts

## Learning objectives

- Sketch the band diagram of an ideal metal–semiconductor (Schottky)
  junction and identify the barrier height $\phi_B$.
- Derive the thermionic-emission current density
  $J = A^* T^2\exp(-\phi_B/V_t)\bigl(\exp(V/V_t)-1\bigr)$.
- Distinguish a Schottky-rectifying contact from an ohmic contact and
  state the engine's idealization of the latter.
- Recognize the ohmic Dirichlet construction in `semi/bcs.py` as the
  "infinite recombination velocity at the contact" limit.
- Locate the gate Dirichlet construction (with $\phi_{ms}$) and forward-
  reference the Schottky branch (M16.5, planned, currently raises
  `NotImplementedError` at `semi/bcs.py:190-193`).

## Physical motivation

Devices need contacts. Two distinct kinds matter for kronos-semi:

- **Ohmic contacts** pass current both ways linearly, behaving like a
  Thevenin-equivalent voltage source on the semiconductor's terminals.
  This is the engine's default contact and the only carrier-pinning
  contact today.
- **Gate contacts** sit on top of an oxide and electrically *don't*
  exchange carriers with the semiconductor — they impose only a
  Dirichlet condition on $\psi$, with the value offset by the
  metal–semiconductor work-function difference.

A third kind, **Schottky contacts**, would replace the ohmic Dirichlet
with a thermionic-emission boundary condition that introduces nonlinear
contact resistance. M16.5 is scoped to add this; the schema reserves
`type: "schottky"` and the BC builder raises `NotImplementedError`
([`semi/bcs.py:190-193`](../../semi/bcs.py)). This chapter covers all
three so the reader can tell what is shipped from what is forward-referenced.

## Derivation from first principles

### The Schottky barrier

A metal with work function $\Phi_m$ contacts an n-type semiconductor
with electron affinity $\chi$ and equilibrium Fermi level $E_F$. Before
contact, the metal Fermi level sits $\Phi_m$ below the local vacuum
level; the semiconductor Fermi level sits $\Phi_s = \chi + (E_c - E_F)$
below its local vacuum level. After contact (and a fast equilibration
of charge), the two Fermi levels must align — that is the definition
of equilibrium.

The result is an offset between the metal Fermi level and the
semiconductor conduction band right at the interface, called the
**Schottky barrier**:

$$
\phi_B = \Phi_m - \chi.
\tag{8.1}
$$

(Idealized; in real interfaces the barrier is also affected by surface
states and image-force lowering, but (8.1) is the first-order picture.)
For Pt on n-Si, $\Phi_m \approx 5.65\,\mathrm{eV}$ and
$\chi = 4.05\,\mathrm{eV}$, giving $\phi_B \approx 1.6\,\mathrm{eV}$.
For Al on n-Si, $\Phi_m \approx 4.10\,\mathrm{eV}$, giving
$\phi_B \approx 0.05\,\mathrm{eV}$ — barely a barrier, and Al-on-Si is
the canonical *ohmic* contact for moderately doped silicon.

The barrier creates a potential-energy hill that electrons must climb
to flow from the semiconductor into the metal. The width of the
depletion region under the contact is given by the same formula as a
pn junction (Ch. 7) with $V_{bi}$ replaced by $\phi_B - (E_c - E_F)/q
= \phi_B - V_t\ln(N_c/N_D)$.

### Thermionic emission

In the thermionic-emission model, electrons that have thermal kinetic
energy above $\phi_B$ are emitted from the semiconductor into the
metal. The rate is governed by the Richardson–Dushman law:

$$
J_\mathrm{s\to m} = A^* T^2\,\exp(-\phi_B/V_t),
\qquad
A^* = \frac{4\pi q m_n^* k^2}{h^3} \approx 110\,\mathrm{A/cm^2/K^2}\,\mathrm{(for\ Si)}.
\tag{8.2}
$$

Under applied bias $V$ with the metal positive, the barrier seen from
the metal side is unchanged ($\phi_B$ is a property of the interface),
but the barrier seen from the semiconductor side is lowered to
$\phi_B - V$. The reverse current $J_\mathrm{m\to s}$ stays at the value
in (8.2); the forward $J_\mathrm{s\to m}$ rises by $\exp(V/V_t)$. Net:

$$
\boxed{
J = A^* T^2\,\exp(-\phi_B/V_t)\bigl(\exp(V/V_t) - 1\bigr).
}
\tag{8.3}
$$

This is the **thermionic-emission diode equation** for a Schottky
contact. The shape mirrors the Shockley diode equation (Ch. 7), but
the saturation current $A^* T^2\exp(-\phi_B/V_t)$ depends on barrier
height rather than on minority-carrier diffusion in a bulk. Schottky
diodes typically have orders-of-magnitude larger reverse leakage and
forward currents than pn diodes of comparable size — the saturation
current's exponential is in $\phi_B$ (typically 0.5–1 eV) rather than
$E_g$ (1.12 eV in Si), so $\exp(-\phi_B/V_t)$ is much larger.

### Image-force lowering

Quantum-mechanically, an electron near the metal interface induces an
image charge in the metal which lowers the effective barrier:

$$
\Delta\phi_B = \sqrt{\frac{q|E|}{4\pi\varepsilon_s}},
\tag{8.4}
$$

where $|E|$ is the field at the interface. This adds a few tens of mV
of barrier reduction at typical fields. It is a small correction and
is usually folded into an *effective* $\phi_B$ extracted from
measurement.

### Ohmic contact: the engine's idealization

A real ohmic contact is achieved by a combination of: (a) heavy doping
near the contact so the depletion region is thin enough for tunneling,
turning the Schottky barrier into a tunneling resistor; (b) a metal
chosen so that $\phi_B$ is small ($<3 V_t$, allowing thermionic
emission to dominate). The result is a contact whose I–V is linear and
reciprocal.

Modeling (a) and (b) self-consistently is expensive. The standard
**ideal ohmic-contact** approximation pins both the potential *and* the
carriers at the contact:

$$
\psi|_\mathrm{contact} = V_t\,\mathrm{asinh}(N_\mathrm{net}/(2n_i)) + V_\mathrm{applied},
\qquad
\Phi_n|_\mathrm{contact} = \Phi_p|_\mathrm{contact} = V_\mathrm{applied}.
\tag{8.5}
$$

The first equation is local charge neutrality plus the applied bias,
i.e. the equilibrium $\psi$ shifted by $V_\mathrm{applied}$. The second
sets the quasi-Fermi potentials (Ch. 11) to the applied bias —
equivalently, it pins the carrier densities at their bulk equilibrium
values $n = N_D$, $p = n_i^2/N_D$ (n-side) etc. Substituted into
Boltzmann (3.7), the contact carriers come out the same as bulk
equilibrium values regardless of bias.

What this *means* physically is "the contact is a perfect reservoir
that absorbs any minority carrier that arrives and emits whatever the
quasi-equilibrium density demands." In recombination-velocity
language, the surface recombination velocity $S \to \infty$.

### Gate contact: the engine's idealization

A gate sits on an oxide; no carrier exchange. Only $\psi$ has a
boundary condition, and that condition is the applied bias offset by
the metal–semiconductor work function difference:

$$
\psi_\mathrm{gate} = V_g - \phi_{ms}.
\tag{8.6}
$$

Carriers are not defined in the oxide (Slotboom variables are
ill-defined where $n_i \to 0$, see Ch. 14), so $\Phi_n, \Phi_p$ have no
BC at the gate. The submesh handling in `build_dd_dirichlet_bcs`
([`semi/bcs.py:236-249`](../../semi/bcs.py)) explicitly skips
quasi-Fermi BCs for gate contacts.

### The Schottky branch (planned)

M16.5's deliverable adds a Robin-style contact form to the continuity
rows:

$$
\mathbf{J}_n\!\cdot\!\hat{\mathbf{n}}\,\big|_\mathrm{contact}
   = q v_R\,(n - n_0\exp(V/V_t)),
\tag{8.7}
$$

with $v_R = A^*T^2/(qN_c)$ the Richardson recombination velocity. The
schema will accept `type: "schottky"`; the contact pins neither $\psi$
nor carriers but adds a flux through the boundary integral. Forward
reference: [`docs/IMPROVEMENT_GUIDE.md` §M16.5](../IMPROVEMENT_GUIDE.md);
acceptance test is the `schottky_1d` benchmark, planned.

## Key results

- Schottky barrier: (8.1).
- Thermionic-emission diode: (8.3).
- Image-force lowering: (8.4).
- Ohmic contact pinning: (8.5).
- Gate contact $\psi$ BC: (8.6).
- Schottky thermionic flux BC (planned): (8.7).

## Worked numerical example

**Pt on n-Si Schottky barrier.** $\Phi_m \approx 5.65\,\mathrm{eV}$,
$\chi = 4.05\,\mathrm{eV}$. $\phi_B = 5.65 - 4.05 = 1.60\,\mathrm{eV}$.

Thermionic saturation current at 300 K, $A^* = 110\,\mathrm{A/cm^2/K^2}$:
$J_s = 110\cdot 300^2\cdot \exp(-1.60/0.02585)
= 9.9\times 10^6\cdot \exp(-61.89)
= 9.9\times 10^6\cdot 1.07\times 10^{-27}
= 1.06\times 10^{-20}\,\mathrm{A/cm^2}$.

This is a vanishingly small reverse leakage — Pt-Si is essentially a
perfect rectifier in this idealization. Real measurements give
$\sim 10^{-7}\,\mathrm{A/cm^2}$ at room temperature; the discrepancy is
image-force lowering plus interface-state effects, and it shows up as
an *effective* barrier of about 0.85 eV.

**Aluminium on n-Si.** $\Phi_m = 4.10\,\mathrm{eV}$,
$\phi_B = 0.05\,\mathrm{eV}$. The barrier is comparable to $V_t$
($25.85\,\mathrm{mV}$) — almost no barrier. The contact is essentially
ohmic. Saturation current: $J_s = 9.9\times 10^6\cdot \exp(-0.05/0.02585)
= 9.9\times 10^6\cdot 0.144 = 1.42\times 10^6\,\mathrm{A/cm^2}$ — a
huge number, larger than any conduction current the device can carry,
which is exactly the regime where the contact stops limiting current
flow and acts as an ideal voltage source. This is *why* aluminium on
moderately doped n-silicon is the textbook ohmic contact.

**Ohmic Dirichlet at the M1 anode.** $N_A = 10^{17}\,\mathrm{cm^{-3}}$.
$\psi^\mathrm{eq}_L = -V_t\,\mathrm{asinh}(5\times 10^6) = -0.4167\,\mathrm{V}$.
At $V = 0$ on the anode: $\psi_L = -0.4167\,\mathrm{V}$; $\Phi_n = \Phi_p = 0$.
At $V = -0.6\,\mathrm{V}$ (M2 with anode swept, equivalent to forward
bias when the cathode is held): $\psi_L = -0.4167 - 0.6 = -1.0167\,\mathrm{V}$;
$\Phi_n = \Phi_p = -0.6\,\mathrm{V}$.

## Code map

| Concept | Equation | Code location |
|---|---|---|
| Ohmic Dirichlet on $\psi$ | (8.5) first eq | `semi/bcs.py:181-185, 252-256` |
| Ohmic Dirichlet on $\Phi_n,\Phi_p$ | (8.5) second eq | `semi/bcs.py:255-266` |
| Gate Dirichlet | (8.6) | `semi/bcs.py:186-189, 236-249` |
| Schottky raise | – | `semi/bcs.py:190-193` (currently raises) |
| Resolve facet refs | – | `semi/bcs.py:51-142` (`resolve_contacts`) |
| Schema valid kinds | – | `semi/bcs.py:34` |
| Equilibrium $\psi$ helper | (8.5) first eq with $V=0$ | `semi/physics/slotboom.py:73-83` |
| Schottky planning | (8.7) | `docs/IMPROVEMENT_GUIDE.md` §M16.5 |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §3.1, §3.2](../PHYSICS.md) — ohmic and gate BCs.
- [`docs/PHYSICS_INTRO.md` §4.1, §4.2](../PHYSICS_INTRO.md) — narrative version.
- [`docs/adr/0007-contact-bc-interface.md`](../adr/0007-contact-bc-interface.md) — ContactBC dataclass design.

## Common pitfalls

1. **Ohmic ≠ "no resistance".** The ideal ohmic contact has zero
   resistance at the contact interface, but it is still a Dirichlet
   condition: it pins the *potential and the carriers* there, and any
   real device's apparent contact resistance is the bulk resistance of
   the leads, not the contact itself. If you want to model finite
   contact resistance you need a Schottky-with-tunneling or
   transfer-length-method treatment, neither of which is shipped.
2. **Gate contact has no carrier BC.** A common error is to set
   $\Phi_n$ at a gate contact. The engine ignores it (the gate facets
   live on the oxide submesh boundary, where $\Phi_n$ doesn't exist),
   but the JSON validator may not catch the mistake.
3. **$\phi_{ms}$ sign.** The convention in `semi/bcs.py:186-188` is
   $\psi_\mathrm{gate} = V_g - \phi_{ms}$. Some textbooks define
   $\phi_{ms} = \Phi_m - \Phi_s$ as the metal-minus-semiconductor work
   function; some define it the other way. Match the convention to
   what the M6 MOSCAP benchmark uses ($\phi_{ms} = -0.95\,\mathrm{V}$
   for n+ poly on p-Si, [`semi/cv.py:69`](../../semi/cv.py)).
4. **Schottky failure is loud, not silent.** If you build a config
   with `contacts[*].type == "schottky"`, [`semi/bcs.py:190-193`](../../semi/bcs.py)
   raises `NotImplementedError`. Don't expect a silent fall-through to
   ohmic.
5. **Image-force lowering is missing.** Even if M16.5 lands, the
   first cut will not include image-force lowering; that is a separate
   correction that adds $\sim 30\,\mathrm{mV}$ to extracted barrier
   heights. Watch the M16.5 acceptance test (`schottky_1d` benchmark)
   for the residual mismatch.

## Exercises

**Exercise 8.1.** Compute $\phi_B$ for Au on n-GaAs given
$\Phi_m^\mathrm{Au} = 5.10\,\mathrm{eV}$ and $\chi^\mathrm{GaAs} =
4.07\,\mathrm{eV}$.

**Exercise 8.2.** A Schottky contact has $\phi_B = 0.85\,\mathrm{eV}$
and $A^* = 110\,\mathrm{A/cm^2/K^2}$. Compute the saturation current at
300 K and at 400 K. By what factor does $J_s$ increase?

**Exercise 8.3.** Show that the ohmic-Dirichlet construction (8.5)
together with Boltzmann statistics gives $n_\mathrm{contact} = N_D$ on
the n-side, regardless of applied bias.

**Exercise 8.4.** Read the gate-contact handling at
[`semi/bcs.py:236-249`](../../semi/bcs.py). Why does the code apply
only a single Dirichlet BC (on $\psi$) at gate contacts, even though it
is consuming the same `bcs` list as ohmic contacts?

**Exercise 8.5.** A heavily doped pn junction contact ($N = 10^{20}\,\mathrm{cm^{-3}}$)
has a depletion region $\sim 5\,\mathrm{nm}$ under the metal. Estimate
the tunneling probability across this barrier (use the WKB approximation
$T \approx \exp(-2\sqrt{2 m^*\phi_B}\cdot W/\hbar)$). Why does this make
heavy doping the standard ohmic-contact recipe?

### Solutions

**8.1.** $\phi_B = 5.10 - 4.07 = 1.03\,\mathrm{eV}$.

**8.2.** $J_s(300) = 110\cdot 300^2\cdot \exp(-0.85/0.02585)
= 9.9\times 10^6\cdot 5.42\times 10^{-15} = 5.36\times 10^{-8}\,\mathrm{A/cm^2}$.
$J_s(400) = 110\cdot 400^2\cdot \exp(-0.85/(k\cdot 400/q))
= 1.76\times 10^7\cdot \exp(-0.85/0.0345)
= 1.76\times 10^7\cdot 1.96\times 10^{-11} = 3.45\times 10^{-4}\,\mathrm{A/cm^2}$.
Ratio $\approx 6400$. Schottky reverse leakage is *very* sensitive to
temperature.

**8.3.** $\psi_\mathrm{contact} = V_t\,\mathrm{asinh}(N_D/(2n_i)) + V$.
$\Phi_n = V$. Substitute into (3.7):
$n = n_i\exp((\psi-\Phi_n)/V_t) = n_i\exp(\mathrm{asinh}(N_D/(2n_i)))
= n_i \cdot (N_D/(2n_i) + \sqrt{1 + (N_D/(2n_i))^2})
\approx n_i \cdot N_D/n_i = N_D$ for $N_D \gg n_i$. ✓

**8.4.** Carriers (Slotboom $\Phi_n, \Phi_p$) live only on the
semiconductor submesh; gate facets are on the oxide-side boundary of
that submesh, where $\Phi_n, \Phi_p$ don't exist. The single $\psi$
Dirichlet is what couples the gate's applied voltage to the silicon's
electrostatic potential through the interfacial flux-continuity
condition. This is also why no `phi_n_bc` or `phi_p_bc` is built at
the gate.

**8.5.** $\sqrt{2 m^*\phi_B} = \sqrt{2\cdot 0.26\cdot 9.11\times 10^{-31}\cdot
1.6\times 10^{-19}\cdot 0.85} \approx 4.0\times 10^{-25}\,\mathrm{kg^{1/2}\,J^{1/2}}
= 4.0\times 10^{-25}/\hbar = 4.0\times 10^{-25}/1.055\times 10^{-34}
= 3.8\times 10^9\,\mathrm{m^{-1}}$. Then
$T \approx \exp(-2\cdot 3.8\times 10^9\cdot 5\times 10^{-9}) = \exp(-38) = 3.1\times 10^{-17}$.
That is small per electron — but at $10^{20}\,\mathrm{cm^{-3}}$
electron density and saturation velocity $\sim 10^7\,\mathrm{cm/s}$,
the tunneling current is $qNvT \sim 10^{-2}\,\mathrm{A/cm^2}$, well above
practical thermionic limits. Heavy doping makes the depletion thin
enough that even highly suppressed tunneling carries a useful current,
which is what makes the contact effectively ohmic.

## Further reading

- **Sze and Ng (2007), Chapter 3** — the canonical Schottky-contact
  reference, including thermionic-emission diffusion theory.
- **Pierret (1996), §13–14** for ohmic-contact engineering and the
  doping-tunneling tradeoff.
- **Hu (2010), Chapter 6** for an undergraduate-friendly treatment with
  numerical examples on Pt, Au, and Al on Si.
- **Crowell, C. R., and Sze, S. M. (1966).** "Current transport in
  metal-semiconductor barriers." *Solid State Electron.* 9, 1035.
  Foundation paper for thermionic-emission diffusion theory.
- **Bardeen, J. (1947).** "Surface states and rectification at a
  metal-semiconductor contact." *Phys. Rev.* 71, 717. The interface-state
  picture that explains why real $\phi_B$ is often Fermi-level-pinned
  and not equal to $\Phi_m - \chi$.
