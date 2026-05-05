# Appendix A — Constants and units

## Universal physical constants

All values are the 2019 SI redefinition (exact). Loaded in
[`semi/constants.py:14-21`](../../semi/constants.py).

| Symbol | Value | Units | Notes |
|---|---|---|---|
| $q$ | $1.602176634\times 10^{-19}$ | C | Elementary charge. Exact (2019 SI). |
| $k_B$ | $1.380649\times 10^{-23}$ | J/K | Boltzmann constant. Exact (2019 SI). |
| $\varepsilon_0$ | $8.8541878128\times 10^{-12}$ | F/m | Vacuum permittivity. CODATA 2018. |
| $\hbar$ | $1.054571817\times 10^{-34}$ | J·s | Reduced Planck. CODATA 2018. |
| $m_0$ | $9.1093837015\times 10^{-31}$ | kg | Electron rest mass. CODATA 2018. |
| $N_A$ | $6.02214076\times 10^{23}$ | 1/mol | Avogadro. Exact (2019 SI). |
| $c$ | $299792458$ | m/s | Speed of light. Exact (2019 SI). |

## Derived constants used by the engine

| Symbol | Definition | Value at 300 K | Units |
|---|---|---|---|
| $V_t$ | $kT/q$ | $0.025852$ | V |
| $\beta$ | $1/(kT)$ | $38.68$ | 1/eV |
| $\varepsilon_0 V_t/q$ | $L_D^2$ scale | $1.428\times 10^{-12}/N_\mathrm{m^{-3}}$ | m² |

## Material parameters (300 K)

See Ch. 20 / [`semi/materials.py`](../../semi/materials.py) for the
full database. Quick reference:

| Material | $\varepsilon_r$ | $E_g$ (eV) | $\chi$ (eV) | $N_c$ (cm⁻³) | $N_v$ (cm⁻³) | $n_i$ (cm⁻³) | $\mu_n$ (cm²/Vs) | $\mu_p$ (cm²/Vs) |
|---|---|---|---|---|---|---|---|---|
| Si | 11.7 | 1.12 | 4.05 | 2.86×10¹⁹ | 3.10×10¹⁹ | 1.0×10¹⁰ | 1400 | 450 |
| Ge | 16.0 | 0.66 | 4.0 | 1.04×10¹⁹ | 6.0×10¹⁸ | 2.0×10¹³ | 3900 | 1900 |
| GaAs | 12.9 | 1.424 | 4.07 | 4.7×10¹⁷ | 7.0×10¹⁸ | 2.1×10⁶ | 8500 | 400 |
| SiO₂ | 3.9 | – | – | – | – | – | – | – |
| HfO₂ | 25.0 | – | – | – | – | – | – | – |
| Si₃N₄ | 7.5 | – | – | – | – | – | – | – |

Insulators have only $\varepsilon_r$; the carrier-related fields are
zero by convention (Ch. 20).

## Unit convention

| Quantity | Internal (SI) | JSON input | Conversion |
|---|---|---|---|
| Length | m | m | – |
| Potential | V | V | – |
| Density (carriers, doping) | m⁻³ | **cm⁻³** | `cm3_to_m3`: × 10⁶ |
| Mobility | m²/(V·s) | **cm²/(V·s)** | `cm2_to_m2`: × 10⁻⁴ |
| Current density | A/m² | A/m² | – |
| Time | s | s | – |
| Lifetimes $\tau_n, \tau_p$ | s | s | – |
| Saturation velocity (M16.1) | m/s | **cm/s** | × 10⁻² |
| Temperature | K | K | – |

The two **bolded** units are the device-physics legacy: every textbook
quotes doping in cm⁻³ and mobility in cm²/(V·s), and the JSON contract
preserves this. The conversion happens *exactly once* at ingest in
[`semi/constants.py:29-46`](../../semi/constants.py).

## Common conversions

- 1 nm = $10^{-9}$ m
- 1 µm = $10^{-6}$ m
- 1 cm = $10^{-2}$ m
- 1 cm⁻³ = $10^6$ m⁻³
- 1 cm²/(V·s) = $10^{-4}$ m²/(V·s)
- 1 eV = $1.602\times 10^{-19}$ J
- 1 mV = $10^{-3}$ V
- $V_t$ at 300 K = 25.85 mV = 25.85×10⁻³ V

## Useful numerical landmarks (Si at 300 K, $C_0 = 10^{17}\,\mathrm{cm^{-3}}$)

| Quantity | Value | Source |
|---|---|---|
| $V_t$ | 25.85 mV | $kT/q$ at 300 K |
| $\sqrt{\varepsilon_0 V_t/qC_0}$ (bare $L_D$) | 3.78 nm | Ch. 12 |
| $\sqrt{\varepsilon_r\varepsilon_0V_t/qC_0}$ ($L_D$ in Si) | 12.96 nm | Ch. 12 |
| $\lambda^2$ at $L_0 = 2\,\mu\mathrm{m}$ | $3.57\times 10^{-6}$ | Ch. 12 |
| $V_{bi}$ for $10^{17}/10^{17}$ | 0.834 V | Ch. 7 |
| $W$ at $V = 0$ | 147 nm | Ch. 7 |
| $|E_\mathrm{max}|$ at $V = 0$ | 113 kV/cm | Ch. 7 |
| $\tau_n = \tau_p$ | 100 ns | engine default |
| $D_n$ in Si | $3.62\times 10^{-3}$ m²/s | Einstein from $\mu_n$ |
| $L_n = \sqrt{D_n\tau_n}$ | 19 µm | Ch. 5 |
| $J_s$ for symmetric $10^{17}$ Si diode | $4.78\times 10^{-8}$ A/m² | Ch. 5 |
| $C_{ox}$ for 5 nm SiO₂ | $6.91\times 10^{-3}$ F/m² | Ch. 9 |
| $C_{ox}$ for 10 nm SiO₂ | $3.45\times 10^{-3}$ F/m² | Ch. 9 |
| Silicon breakdown field | $\sim 3\times 10^5$ V/cm | textbook |

## Mixed-unit pitfalls

1. **JSON doping in cm⁻³.** Forgetting this → engine sees $10^{17}\,\mathrm{m^{-3}}$
   instead of $10^{17}\,\mathrm{cm^{-3}} = 10^{23}\,\mathrm{m^{-3}}$ — six
   orders of magnitude under-doped.
2. **JSON mobility in cm²/(V·s).** Forgetting → engine sees $1400\,\mathrm{m^2/(V\cdot s)}$
   instead of $1400\,\mathrm{cm^2/(V\cdot s)} = 0.14\,\mathrm{m^2/(V\cdot s)}$.
3. **eV vs V.** Trap energy $E_t$ is specified in **eV** in the schema
   but enters the formula as $E_t/V_t$ where $V_t$ is in V. The conversion
   $E_t[\mathrm{eV}]/V_t[\mathrm{V}]$ is dimensionless because 1 eV / 1 V = 1.
   The engine handles this at [`semi/runners/bias_sweep.py:90`](../../semi/runners/bias_sweep.py)
   (`E_t_over_Vt = E_t_eV / sc.V0`).
4. **Current density per unit area vs per unit length in 1D.** A 1D
   simulation gives $J$ in A/m² because the integration measure is
   per-cell (no transverse area concept). A 2D simulation gives current
   *per unit z-width* — A/m. A 3D simulation gives total current — A.
   Be careful when comparing across dimensions; the M14.3 MOSFET
   verifier divides by `L_drain_line` (the contact's vertical extent in
   2D) to recover $I_D/W$ matching the Pao-Sah formula.

## Rounding conventions

The engine's tolerances and benchmarks use the following rounding
conventions:

- $V_t$ at 300 K: 25.85 mV (4 significant figures).
- $V_{bi}$ for symmetric $10^{17}$ Si: 0.83 or 0.834 V depending on
  precision (the asinh-form gives 0.8334; log-form 0.8334; agreement
  to 4 places).
- $W$, $|E_\mathrm{max}|$: 3 significant figures
  ($147\,\mathrm{nm}$, $113\,\mathrm{kV/cm}$).
- MMS rates: reported to 3 decimal places (e.g., 2.000, 1.998).
- Conservation residuals: scientific notation, e.g. $1.5\times 10^{-17}$.

## Typical scales by device class

| Device | $L_0$ | $C_0$ (cm⁻³) | $\lambda^2$ | Operating biases |
|---|---|---|---|---|
| Long pn diode | 100 µm | $10^{15}$ | $\sim 10^{-3}$ | $\pm 1$ V |
| Short pn diode (M1) | 2 µm | $10^{17}$ | $\sim 10^{-6}$ | $\pm 1$ V |
| MOSCAP (M6) | 0.5 µm | $10^{17}$ | $\sim 10^{-5}$ | $\pm 2$ V |
| MOSFET (M12) | 1 µm | $10^{17}$ | $\sim 10^{-6}$ | 0 to 2 V |
| Thin-channel MOSFET | 50 nm | $10^{18}$ | $\sim 10^{-3}$ | 0 to 1 V |

## See also

- [`semi/constants.py`](../../semi/constants.py) — module source.
- Ch. 12 — full nondimensionalization derivation.
- Ch. 20 — material parameter database.
- Appendix C — symbol/glossary.
