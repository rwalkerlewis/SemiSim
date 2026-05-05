# 20 — Material parameter database

## Learning objectives

- List the materials shipped in `semi/materials.py` (Si, Ge, GaAs;
  SiO₂, HfO₂, Si₃N₄) and recognize the parameters each one carries.
- Trace each numerical value to a primary source (Sze 3rd ed. tables;
  Altermatt 2003 for $n_i^\mathrm{Si}$).
- Explain why insulators carry only $\varepsilon_r$ and how the engine
  enforces the semiconductor-only fields are zero for insulators.
- Recognize the temperature dependence of the parameters and where the
  300 K assumption is baked in.
- Predict the impact of replacing one material with another in a given
  benchmark.

## Physical motivation

Every device simulation needs material parameters: relative permittivity
to compute the displacement field, bandgap and effective DOS to compute
$n_i$, mobilities to compute the currents. Hard-coding numbers into
each runner would be unmaintainable; the engine factors them out into
a single `Material` dataclass and a `MATERIALS` registry. This chapter
catalogues each entry, gives the physical origin of each number, and
flags assumptions (temperature, bandgap narrowing, doping dependence)
that future milestones will need to relax.

## Materials catalogue

### Silicon (`Si`)

```python
Material(
    name="Si",
    role="semiconductor",
    epsilon_r=11.7,
    Eg=1.12,                 # eV
    chi=4.05,                # eV
    Nc=cm3_to_m3(2.86e19),   # m^-3
    Nv=cm3_to_m3(3.10e19),
    n_i=cm3_to_m3(1.0e10),   # Altermatt 2003
    mu_n=cm2_to_m2(1400.0),  # m^2/Vs
    mu_p=cm2_to_m2(450.0),
)
```

- **$\varepsilon_r = 11.7$.** Sze 3rd ed. Appendix C; matches every
  textbook within 1%. At cryogenic temperatures it rises to 12.0
  (~270 K to 4 K).
- **$E_g = 1.12\,\mathrm{eV}$.** At 300 K. Decreases by 0.4 meV/K
  above 100 K; at 77 K, $E_g = 1.16\,\mathrm{eV}$. The engine uses the
  300 K value throughout; cryogenic device simulation would need a
  $T$-dependent $E_g$ (M16+ work).
- **$\chi = 4.05\,\mathrm{eV}$.** Electron affinity. Used in
  $\phi_{ms}$ for gate contacts (Ch. 9) and in Schottky barrier height
  (Ch. 8).
- **$N_c = 2.86\times 10^{19}\,\mathrm{cm^{-3}}$.** Sze 3rd ed. Table 7,
  derived from $m_n^*\approx 1.08\,m_0$ via (2.6).
- **$N_v = 3.10\times 10^{19}\,\mathrm{cm^{-3}}$.** Sze Table 7,
  $m_p^*\approx 1.15\,m_0$.
- **$n_i = 1.0\times 10^{10}\,\mathrm{cm^{-3}}$.** Altermatt et al.
  2003 consensus value; supersedes the older Sze 1.45×10¹⁰. Modern
  TCAD tools (Sentaurus, Atlas) match this; comparing against legacy
  texts may give 19 mV offsets in $V_{bi}$
  ([`docs/PHYSICS.md` §4.2](../PHYSICS.md) note).
- **$\mu_n = 1400\,\mathrm{cm^2/Vs}$, $\mu_p = 450\,\mathrm{cm^2/Vs}$.**
  Undoped, low-field, 300 K bulk silicon. Caughey-Thomas (M16.1, shipped)
  uses these as the reference $\mu_0$ in (5.8); Lombardi (M16.2,
  planned) refines for surface scattering.

### Germanium (`Ge`)

```python
Material("Ge", role="semiconductor",
         epsilon_r=16.0, Eg=0.66, chi=4.0,
         Nc=cm3_to_m3(1.04e19), Nv=cm3_to_m3(6.0e18),
         n_i=cm3_to_m3(2.0e13),
         mu_n=cm2_to_m2(3900.0), mu_p=cm2_to_m2(1900.0))
```

- $E_g = 0.66\,\mathrm{eV}$ — much smaller than silicon, so $n_i$ is
  $\sim 2000\times$ larger ($n_i = 2\times 10^{13}$ vs Si's $10^{10}$).
- $\mu_n, \mu_p$ are higher than silicon — Ge is a higher-mobility
  semiconductor. SiGe heterostructures exploit this.
- Sze 3rd ed. Table 7 source.

### Gallium arsenide (`GaAs`)

```python
Material("GaAs", role="semiconductor",
         epsilon_r=12.9, Eg=1.424, chi=4.07,
         Nc=cm3_to_m3(4.7e17), Nv=cm3_to_m3(7.0e18),
         n_i=cm3_to_m3(2.1e6),
         mu_n=cm2_to_m2(8500.0), mu_p=cm2_to_m2(400.0))
```

- $E_g = 1.424\,\mathrm{eV}$ — direct gap. $n_i$ is much smaller
  ($\sim 2\times 10^6$) because the gap is wider and indirect-gap
  thermal generation isn't available.
- $N_c = 4.7\times 10^{17}\,\mathrm{cm^{-3}}$ — *much* smaller than Si
  because $m_n^*\approx 0.067\,m_0$ in GaAs (light electrons). This is
  the hallmark of direct-gap semiconductors.
- $\mu_n = 8500\,\mathrm{cm^2/Vs}$ — about 6× silicon's. Drives the
  HEMT mobility advantage.

### Silicon dioxide (`SiO2`)

```python
Material("SiO2", role="insulator", epsilon_r=3.9)
```

- The standard MOS gate dielectric. $\varepsilon_r = 3.9$ at low
  frequency.
- All semiconductor-only fields ($E_g$, $\chi$, $N_c$, $N_v$, $n_i$,
  $\mu_n$, $\mu_p$) default to zero. The actual SiO₂ bandgap is ~9 eV
  but this is irrelevant to the DD model used by the engine (no
  carriers in oxide); see Ch. 1 §pitfall 4.

### Hafnium oxide (`HfO2`)

```python
Material("HfO2", role="insulator", epsilon_r=25.0)
```

- High-k dielectric replacing SiO₂ in advanced CMOS gate stacks.
- $\varepsilon_r = 25$: a 10 nm HfO₂ layer is electrically equivalent to
  about 1.56 nm of SiO₂ (Ch. 9 Exercise 9.3) — much thinner EOT for
  the same physical thickness, which suppresses gate leakage.

### Silicon nitride (`Si3N4`)

```python
Material("Si3N4", role="insulator", epsilon_r=7.5)
```

- Used as a passivation layer and in SONOS memory.

### `Material` dataclass

```python
@dataclass
class Material:
    name: str
    role: str             # "semiconductor" or "insulator"
    epsilon_r: float
    Eg: float = 0.0
    chi: float = 0.0
    Nc: float = 0.0
    Nv: float = 0.0
    n_i: float = 0.0
    mu_n: float = 0.0
    mu_p: float = 0.0

    @property
    def epsilon(self) -> float:
        return self.epsilon_r * EPS0

    def is_semiconductor(self) -> bool:
        return self.role == "semiconductor"

    def is_insulator(self) -> bool:
        return self.role == "insulator"
```

The `epsilon` property gives the absolute permittivity in F/m. The
boolean role distinction is what the BC builder, the runner, and the
multi-region Poisson form use to decide whether the material carries
carriers, doping, etc.

### `MATERIALS` registry and the HTTP endpoint

[`semi/materials.py:49-89`](../../semi/materials.py) is a module-level
`dict[str, Material]`. Lookups via `get_material(name)`. The HTTP
server's `GET /materials` endpoint
([`kronos_server/routes/health.py`](../../kronos_server/routes/health.py))
returns this dict for UI autocomplete and validation. The UI never
imports the engine; it reads the JSON dump and constructs the Region
inputs against this allowed-name list.

## Worked numerical example: replace Si with GaAs in pn_1d

The M1 benchmark uses Si. Suppose a user replaces `"material": "Si"`
with `"material": "GaAs"` while keeping doping at $10^{17}\,\mathrm{cm^{-3}}$
on each side:

- $V_t$ unchanged ($V_t = kT/q = 25.85\,\mathrm{mV}$ at 300 K).
- $n_i^\mathrm{GaAs} = 2.1\times 10^6\,\mathrm{cm^{-3}}$ instead of $10^{10}$.
- $V_{bi} = V_t\ln(N_AN_D/n_i^2) = V_t\ln(10^{34}/(2.1\times 10^6)^2)
  = V_t\ln(10^{34}/4.41\times 10^{12}) = V_t\ln(2.27\times 10^{21})
  = 0.02585\cdot 49.16 = 1.271\,\mathrm{V}$. Compared to Si's
  0.834 V.
- $\varepsilon_r^\mathrm{GaAs} = 12.9$ vs Si's 11.7. Slight increase
  in $L_D$, slight decrease in $W$ (wider $V_{bi}$ winning).
- Mobility: $\mu_n^\mathrm{GaAs} = 8500\,\mathrm{cm^2/Vs}$, 6× Si.
  Forward currents scale roughly $\propto \mu$; the GaAs diode at
  $V_F = 1.0\,\mathrm{V}$ would have ~6× the forward current of Si at
  $V_F = 1.0\,\mathrm{V}$, but it would not turn on until $V_F\sim 1.27\,\mathrm{V}$
  due to the larger $V_{bi}$. (Real GaAs LEDs are biased above 1.4 V,
  consistent with this picture.)

Replace doping or dimensions, and you can repeat the same calculation.

## Code map

| Concept | Code location |
|---|---|
| `Material` dataclass | `semi/materials.py:17-46` |
| `MATERIALS` registry | `semi/materials.py:49-89` |
| `get_material(name)` | `semi/materials.py:92-98` |
| `list_materials()` | `semi/materials.py:101-102` |
| `epsilon` property | `semi/materials.py:38-40` |
| `is_semiconductor` / `is_insulator` | `semi/materials.py:42-46` |
| Conversion helpers `cm3_to_m3`, `cm2_to_m2` | `semi/constants.py:29-46` |
| `GET /materials` HTTP endpoint | `kronos_server/routes/health.py` |
| Reference material in runner | `semi/runners/_common.py:reference_material` |

## Existing-docs cross-reference

- [`docs/PHYSICS.md` §4](../PHYSICS.md) — material parameter tables.
- [`docs/ARCHITECTURE.md`](../ARCHITECTURE.md) Layer 3 — `materials.py`
  is pure-Python core.
- [`docs/IMPROVEMENT_GUIDE.md`](../IMPROVEMENT_GUIDE.md) §M17 — heterojunctions
  will need cellwise chi(x), Eg(x) overrides per region.

## Common pitfalls

1. **Material names are case-sensitive.** `"Si"`, `"SiO2"`, `"GaAs"`,
   etc. A typo gives `KeyError`; the schema validator does not catch
   this because the field is just a string.
2. **Insulator material has zero semiconductor fields.** If you query
   `MATERIALS["SiO2"].n_i`, you get 0. Code that uses the `n_i` field
   without first checking `is_semiconductor` will divide by zero (or
   produce NaN) somewhere downstream. The Slotboom-form residual
   correctly avoids this by integrating the carrier rows only on the
   semiconductor submesh (Ch. 14).
3. **Sze tables are at 300 K.** All values are room-temperature.
   Cryogenic operation (e.g. 4 K for quantum computing) requires
   different $E_g$, $n_i$, $\mu$, and the parameter-database approach
   would need a temperature axis.
4. **No bandgap narrowing.** Heavy doping ($N \gtrsim 10^{19}\,\mathrm{cm^{-3}}$)
   reduces the bandgap by tens of meV in silicon (Slotboom-de Graaf
   model). The shipped engine assumes constant $E_g$; this is part of
   why the M16.4 Fermi-Dirac milestone is paired with high-doping
   benchmarks where bandgap narrowing also matters.
5. **No doping-dependent mobility.** The shipped engine's $\mu_n, \mu_p$
   are fixed values. Real silicon mobility falls from
   $1400\,\mathrm{cm^2/Vs}$ at low doping to $\sim 100\,\mathrm{cm^2/Vs}$
   at $10^{19}\,\mathrm{cm^{-3}}$ due to ionized-impurity scattering.
   Caughey-Thomas (M16.1) has only field dependence; doping dependence
   would need a Klaassen-style composite (post-M16+).

## Exercises

**Exercise 20.1.** Compute $n_i$ for Ge at 300 K from $N_c, N_v, E_g$
in `Material("Ge")` via (2.9). Compare to the stored value.

**Exercise 20.2.** A user replaces SiO₂ with HfO₂ at 5 nm thickness on
the M6 MOSCAP. Compute $C_{ox}$ and predict the new $V_T$ assuming
$\phi_{ms}$ unchanged.

**Exercise 20.3.** Why does silicon's $N_c$ exceed silicon's $n_i$ by
nine orders of magnitude? Connect to the band-structure picture in Ch. 2.

**Exercise 20.4.** Read [`semi/runners/_common.py:reference_material`](../../semi/runners/_common.py).
What does it return? Why does the Scaling object need a single
"reference material" rather than per-region material parameters?

**Exercise 20.5.** Sketch what additions to the `Material` dataclass
M16.1 (Caughey-Thomas) and M17 (heterojunctions) would require. Why
does the engine not already carry $v_\mathrm{sat}$ on `Material`?

### Solutions

**20.1.** $n_i^2 = N_c N_v\exp(-E_g/kT)
= 1.04\times 10^{19}\cdot 6\times 10^{18}\cdot\exp(-0.66/0.02585)
= 6.24\times 10^{37}\cdot \exp(-25.53)
= 6.24\times 10^{37}\cdot 8.18\times 10^{-12} = 5.10\times 10^{26}\,\mathrm{cm^{-6}}$.
$n_i = 2.26\times 10^{13}\,\mathrm{cm^{-3}}$. Stored: $2.0\times 10^{13}$.
Match within 13% (consistent with the slight effective-mass-value
discrepancy across Sze tables).

**20.2.** $C_{ox}^\mathrm{HfO_2}(5\,\mathrm{nm}) = 25\cdot 8.854\times 10^{-12}/5\times 10^{-9}
= 4.43\times 10^{-2}\,\mathrm{F/m^2}$. About 6.4× SiO₂'s value at the same
thickness.
Depletion-charge term in $V_T$: same $\sqrt{2\varepsilon_s qN_a\cdot 2|\phi_B|}$
divided by 6.4× larger $C_{ox}$, so 6.4× smaller. For Hu Fig. 5-18
($N_a = 5\times 10^{16}$): original term was 0.333 V; new term is
$0.333/6.4 = 0.052\,\mathrm{V}$. $V_T = V_{fb} + 2|\phi_B| + 0.052
= -0.950 + 0.798 + 0.052 = -0.100\,\mathrm{V}$ (negative — the device
turns on at *negative* gate voltage, which means the p-body MOSCAP is
in inversion at $V_g = 0$).

**20.3.** $N_c$ counts the (3D parabolic-band) states available within
$\sim kT$ of the conduction-band edge — about $10^{19}\,\mathrm{cm^{-3}}$
in silicon. $n_i$ counts electrons *thermally promoted* to those
states across an $E_g = 1.12\,\mathrm{eV}$ gap. The promotion factor
is $\exp(-E_g/2kT) \sim 10^{-9}$, so $n_i \sim N_c \cdot 10^{-9}\sim 10^{10}\,\mathrm{cm^{-3}}$.
The big gap makes thermal carriers very rare relative to available states.

**20.4.** [`semi/runners/_common.py:reference_material`](../../semi/runners/_common.py)
returns the first semiconductor material in the regions list — the
"reference" material whose $\mu_n$, $n_i$ are used as the scaling
references $\mu_0, n_i^\mathrm{ref}$. Multi-material devices (M17,
planned) would need per-region parameters in the residual; the
*scaling* still uses a single reference (so the scaled equations have
a single $\mu_0$ etc.), and per-region departures are absorbed as
$\hat\mu_n^{(\mathrm{region})}/\hat\mu_n^{(\mathrm{ref})}$ ratios.
Currently every benchmark has one semiconductor material, so this
distinction doesn't matter; M17 will surface it.

**20.5.** M16.1 Caughey-Thomas adds $v_\mathrm{sat,n}, v_\mathrm{sat,p},
\beta_n, \beta_p$ as *physics.mobility* JSON fields, not as `Material`
fields, because they are model parameters rather than material
constants. (Per-material defaults could go on `Material`, but the
shipped M16.1 puts them on the schema with global defaults.) M17
heterojunctions would add per-region $\chi$, $E_g$ overrides via JSON
`material_overrides`, leaving `Material` itself unchanged. Adding
$v_\mathrm{sat}$ to `Material` is reasonable but would couple the
material registry to the mobility model; the current factoring
(material = bulk parameters; mobility model = velocity-saturation)
keeps responsibilities separate.

## Further reading

- **Sze and Ng (2007), Appendix C.** The reference table for
  silicon, germanium, GaAs, and other III-V semiconductors. Source
  for most of `semi/materials.py`'s numbers.
- **Altermatt, P. P., et al. (2003).** "Reassessment of the intrinsic
  carrier density in crystalline silicon." *J. Appl. Phys.* 93,
  1598. The source for the modern $n_i^\mathrm{Si} = 1.0\times 10^{10}$.
- **Klaassen, D. B. M. (1992).** "A unified mobility model for device
  simulation — I." *Solid-State Electron.* 35, 953. Doping-dependent
  mobility model that future engine versions could adopt.
- **Slotboom, J. W., and de Graaff, H. C. (1976).** "Measurements of
  bandgap narrowing in Si bipolar transistors." *Solid-State
  Electron.* 19, 857. Bandgap narrowing in heavy doping; not yet in
  the engine.
