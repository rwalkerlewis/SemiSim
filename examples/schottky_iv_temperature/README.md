# schottky_iv_temperature

## What this is

A Pt-on-n-Si Schottky diode I-V characterization across three
temperatures (250 K, 300 K, 350 K). Geometry mirrors
[`benchmarks/schottky_1d`](../../benchmarks/schottky_1d) for parameter
consistency: 1D, uniform N_D = 1e16 cm^-3 over a 5 um device, Pt
barrier height 0.85 eV (Sze 3rd ed Table 5), Boltzmann statistics.
The new axis is temperature: the anchor JSON pins T = 300 K and the
smoke verifier reruns the configuration at T = 250 K and T = 350 K
via in-process re-runs that override `physics.temperature`. The
intended use case is "how does my Schottky diode I-V shift when the
package warms up by 50 K," and the config is meant to be cloned and
re-parametrized for your own metal-semiconductor junction.

This is an example, not a V&V gate. For the analytical thermionic-
emission match at fixed T, see
[`benchmarks/schottky_1d`](../../benchmarks/schottky_1d) and ADR 0015.

## Physics features exercised

- **M16.5 Schottky thermionic-emission Robin BC.**
  `contacts[0].type = "schottky"` with
  `barrier_height_eV = 0.85` (Pt on n-Si). The BC pair is the
  metal-Fermi-level psi Dirichlet plus a thermionic-emission Robin
  form on the electron continuity row; see
  [`docs/adr/0015-schottky-robin-bc.md`](../../docs/adr/0015-schottky-robin-bc.md)
  for the V&V scope discussion.
- **Temperature dependence.** `physics.temperature` is the parametric
  knob; the verifier sweeps over {250 K, 300 K, 350 K}. The
  temperature enters thermionic emission through the Richardson
  pre-factor `J_sat = A* T^2 exp(-q phi_B / kT)` and the thermal
  voltage `kT/q`.
- **SRH recombination** (default, mostly inactive at the moderate
  doping and low forward bias of this geometry).

## How to run

From the repository root:

```
docker compose run --rm benchmark schottky_iv_temperature
```

Output lands in `results/schottky_iv_temperature/`:

- `schottky_iv_temperature.png`: three I-V curves overlaid on a
  semilog-y axis, color-coded by temperature.
- The result JSON written by the engine is in the standard
  `results/schottky_iv_temperature/` location for the anchor
  T = 300 K run; the companion runs land in
  `_T_<temperature>K/` subdirectories.

## Expected output

The thermionic-emission saturation current scales as

    J_sat = A* T^2 exp(-q phi_B / (k T))

For Pt on n-Si (phi_B = 0.85 eV) the dominant temperature dependence
is the exponential. The order-of-magnitude expectations at V_F = 0 V
are

    T = 250 K:  J_sat ~  3 x 10^-7 A/m^2
    T = 300 K:  J_sat ~  3 x 10^-3 A/m^2
    T = 350 K:  J_sat ~  3 x 10^+0 A/m^2

(roughly a 10^4 swing across the temperature range). At V_F = 0.3 V
the diode is well into the exponential turn-on regime and the
ordering I(350 K) > I(300 K) > I(250 K) must hold; this is the
qualitative gate the smoke verifier checks. The simulation matches
those order-of-magnitude expectations to within roughly a factor of
3 (the 5 um geometry adds a finite bulk-drift contribution that the
ideal thermionic-emission formula does not capture).

If the reported ordering at V_F = 0.3 V is wrong (for example
I(250 K) > I(300 K)), check:

- The barrier height (`contacts[0].barrier_height_eV`). It is the
  dominant temperature-coefficient knob; a 0.05 eV change in barrier
  height shifts every saturation current by roughly an order of
  magnitude at room temperature.
- The temperature itself (`physics.temperature`). Below ~150 K the
  Boltzmann approximation breaks more thoroughly; the simple
  thermionic-emission formula is accurate but the FEM
  drift-diffusion path may need additional care.
- The doping (`doping[0].profile.N_D`). Heavier doping shrinks the
  depletion region and pushes the device into the field-emission
  (tunneling) regime where pure thermionic emission underpredicts.

## How to adapt

Most users will want to start by changing one of:

- **Barrier height.** `contacts[0].barrier_height_eV`. Common
  metal-on-Si barrier heights are tabulated in Sze 3rd ed Table 5
  (Au: 0.80 eV, Pt: 0.85 eV, Ni: 0.61 eV). Change this value to
  match your metal.
- **Temperature set.** Edit the verifier (or override
  `physics.temperature` in your own copy of the JSON) to whatever
  temperature points are relevant to your operating envelope.
- **Doping.** Change `doping[0].profile.N_D`. Heavier doping shifts
  the device toward field emission; for very thin barriers
  (>1e18 cm^-3) the M16.6 BBT kernel becomes relevant and you would
  add `physics.tunneling = {bbt: true}` (see
  [`benchmarks/zener_1d`](../../benchmarks/zener_1d) for the
  reverse-bias BBT case).
- **Bias range.** `contacts[0].voltage_sweep`. The default
  [0.0, 0.5] V step 0.025 V covers the exponential turn-on; widen
  to characterize a specific operating point.

## Choice of parametric mechanism

The schema does not have a built-in "sweep over `physics.temperature`"
mechanism; the available `voltage_sweep` field operates on contact
bias only. The cleanest fit for the example pattern is the
in-process re-run approach: ship a single anchor JSON at the
nominal operating point and have the smoke verifier loop the other
parameter values via `semi.run.run()` calls on deepcopies of the
config. This mirrors the pattern used by `diode_velsat_1d`,
`diode_auger_1d`, and `diode_fermi_dirac_1d` for their companion
sweeps; future work that needs broad parameter studies as
first-class config can revisit this once a `parametric` schema
block is justified by enough callers.
