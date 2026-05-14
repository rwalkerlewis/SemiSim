# power_diode_reverse_recovery

## What this is

A practical 1D long-base power-rectifier diode driven through a
canonical reverse-recovery transient: forward conduction at +0.7 V,
abrupt voltage reversal to -2.0 V over a 10 ns linear ramp, then
reverse blocking. The waveform exhibits the standard reverse-recovery
shape, namely a forward conduction current, a current dip during
stored-charge clearance after the voltage reverses, and recovery
toward the small reverse-blocking leakage. The intended use case is
"how does my rectifier look during turn-off," and the config is meant
to be cloned and re-parametrized for your own switching scenario.

This is an example, not a V&V gate. For the analytical small-signal
voltage_t correctness reference, see audit case 06
(`tests/audit/test_06_transient_fft_vs_ac_sweep.py`) and the
demonstrators [`benchmarks/pn_1d_pulse`](../../benchmarks/pn_1d_pulse)
and [`benchmarks/diode_sine_1d`](../../benchmarks/diode_sine_1d). For
the Auger correctness gate see
[`benchmarks/diode_auger_1d`](../../benchmarks/diode_auger_1d).

## Physics features exercised

- **M16.7 voltage_t (table variant).**
  `contacts[0].voltage_t.type = "table"` with a piecewise-linear
  trajectory from +0.7 V (t in [0, 50] ns), through a 10 ns linear
  ramp to -2.0 V (t in [50, 60] ns), followed by reverse blocking
  at -2.0 V (t in [60, 200] ns). The materialized table has 49
  sample points across the full trajectory.
- **M16.3 Auger recombination.**
  `physics.recombination.auger = true` with the schema default Si
  Dziewior-Schmid coefficients (C_n = 2.8e-31 cm^6/s,
  C_p = 9.9e-32 cm^6/s). Auger contributes meaningfully under
  forward conduction because the light doping (1e15 / 5e15 cm^-3)
  puts the device into high injection at +0.7 V.
- **SRH recombination** with tau_n = tau_p = 100 ns. The reverse-
  recovery time scale is set primarily by the SRH lifetime; Auger
  trims the high-injection tail of the forward conduction
  carrier density.
- **Backward-Euler / BDF2 transient.** `solver.type = "transient"`,
  `dt = 1 ns` (initial), `t_end = 200 ns`, BDF2 (`order = 2`).
  `bc_ramp_steps = 20` ramps the +0.7 V initial-condition bias
  through 20 steady-state continuation steps before the time loop
  starts.
- **M18 adaptive dt.** `solver.adaptive.enabled = true`. The
  controller in `semi/continuation.py` (the same class the bias-
  sweep ramp uses) drives dt between `dt_min = 1 ps` and
  `dt_max = 5 ns`, halving on SNES failure and growing on
  consecutive easy SNES solves. dt is also clamped at every
  `voltage_t.table.times[i]` breakpoint so the integrator never
  straddles a slope change. See
  [docs/M18_STARTER_PROMPT.md](../../docs/M18_STARTER_PROMPT.md)
  and [docs/adr/0017-adaptive-transient-dt.md](../../docs/adr/0017-adaptive-transient-dt.md)
  for the controller design and the rationale.

## How adaptive dt is used

The `voltage_t.table` waveform has 49 sample points; the largest
slope changes are at t = 50 ns (start of the +0.7 V to -2.0 V
ramp) and t = 60 ns (end of the ramp into reverse blocking). With
fixed `dt = 1 ns` (the v0.24.0 setting), the BDF2 step at the
ramp boundary mixes a forward-conduction history with a sharp
reverse-bias BC update; SNES stagnates at the second or third
post-transition step because the linearization point is far from
the new operating regime. M18 fixes this in two ways:

1. **Breakpoint clamping.** The runner pre-computes the sorted
   list of interior `times[i]` from the table and clamps each
   step's dt so it lands exactly on the next breakpoint. The
   integrator never straddles a slope change.
2. **Halve on failure.** Across the t = 50 ns to t = 60 ns
   transition the controller is expected to halve dt down to the
   ~50 ps neighborhood (the exact floor depends on the SNES rtol
   / atol setting and how aggressively the ramp turns the
   minority-carrier population around). After the settled
   reverse-blocking tail begins (roughly t > 100 ns) the
   controller grows dt back toward `dt_max = 5 ns` over a handful
   of steps.

The `solver.max_steps = 2000` cap accounts for the smaller dt
during the transition; the v0.24.0 cap (250) was sized for fixed
dt = 1 ns and would exhaust during the halving window.

## How to run

From the repository root:

```
docker compose run --rm benchmark power_diode_reverse_recovery
```

Output lands in `results/power_diode_reverse_recovery/`:

- `iv_recovery.png`: I(t) and V(t) on twin y-axes with the
  conduction, transition, recovery, and settled regions annotated.
- The result JSON written by the engine (with the per-timestep I-V
  table) is in the standard
  `results/power_diode_reverse_recovery/` location.

## Expected output

The waveform has four qualitative phases:

1. **Forward conduction** (t in [0, 50] ns at V = +0.7 V).
   Peak forward current density is order-of-magnitude 1e6 A/m^2
   (the light doping gives a low-resistance forward path; the
   exact value depends on the diffusion-length / device-length
   ratio plus the Auger limiter at high injection).
2. **Transition** (t in [50, 60] ns at V going from +0.7 V to
   -2.0 V). The current crosses zero and goes sharply negative as
   the stored charge starts to be swept out.
3. **Reverse-recovery dip** (t in [60, 100] ns at V = -2.0 V).
   The minimum (most negative) current is order-of-magnitude
   1e5 A/m^2, smaller than the forward peak by roughly 10x. The
   recovery time scale is roughly tau scaled by the depletion-
   transit time, on the order of 50-100 ns.
4. **Settled reverse blocking** (t > ~100 ns at V = -2.0 V).
   The current settles toward the small reverse-saturation
   leakage; final-time current is meaningfully closer to zero
   than the recovery minimum.

If the qualitative shape is wrong (for example no current dip
appears after t = 60 ns, or the post-transition current never
recovers toward zero), check:

- The voltage_t table (`contacts[0].voltage_t.values`) actually
  reaches the negative reverse-bias value at the right time.
- The carrier lifetime (`physics.recombination.tau_n` /
  `tau_p`). Longer tau makes the recovery slower and the dip
  larger; shorter tau makes the recovery faster and the dip
  smaller. The shipped value (100 ns) is moderate.
- The doping (`doping[0].profile`). Heavier doping shrinks the
  diffusion length and the stored charge, so the dip is smaller
  and the recovery is faster.
- The `bc_ramp_steps` (`solver.bc_ramp_steps`). The 20-step
  default in this config ramps the +0.7 V initial-condition bias
  in 20 steady-state continuation steps; the high-injection
  regime at +0.7 V on light doping is nonlinear and the SNES
  needs more than the runner default of 10 steps to converge.

## How to adapt

Most users will want to start by changing one of:

- **Voltage trajectory.** Edit `contacts[0].voltage_t.times` and
  `values`. The shipped trajectory is a single +0.7 V to -2.0 V
  reversal with a 10 ns transition; multi-pulse switching, slower
  ramps, or asymmetric forward / reverse swings are all possible
  by extending the table.
- **Doping.** Change `doping[0].profile.N_A_left` (p side) and
  `N_D_right` (n side). For a more realistic high-voltage
  rectifier, use heavier p-side and light n-side (asymmetric pn-
  diode).
- **Carrier lifetime.** `physics.recombination.tau_n` / `tau_p`
  controls the recovery time scale. Slow rectifier diodes
  intended for line-frequency rectification have tau ~ 1-10 us.
  Fast-recovery rectifiers have tau ~ 50-200 ns.
- **Time grid.** `solver.dt` and `solver.t_end`. The shipped
  configuration runs at 1 ns step out to 200 ns. If you change
  the trajectory time scales, adjust dt accordingly so the
  voltage_t table is well-resolved by the time loop.
- **Auger.** `physics.recombination.auger = false` turns Auger
  off. For lightly-doped diodes at moderate forward bias the
  effect is small; for heavy-injection regimes (very high
  forward current, deep saturation) Auger is essential.
