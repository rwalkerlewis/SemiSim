# nmos_idvgs

## What this is

A practical 2D n-channel MOSFET with a 100 nm gate length, 2 nm SiO2
gate oxide, p-type body (N_A = 1e16 cm^-3), and Gaussian n+
source/drain implants (peak N_D = 1e20 cm^-3). The anchor JSON sweeps
V_GS in [0, 0.5] V at V_DS = 0.05 V (linear regime), exercising
subthreshold conduction through the onset of weak inversion. The
intended use case is "what does an Id-Vgs measurement on a
short-channel-ish bulk Si NMOS look like in the framework," and the
config is meant to be cloned and re-parametrized for your own device.

The bias_sweep coupled SNES inner solve uses the PETSc `nleqerr` line
search (the Deuflhard 2004 natural monotonicity test; ADR 0018, M18.1)
rather than the M12-era `bt` default. The `bt` backtracking line
search cannot find a descent direction at the depletion-to-inversion
onset for this device class; `nleqerr` is opt-in via
`solver.snes.line_search` (M18.1 / schema v2.10.0).

This is an example, not a V&V gate. For the analytical Pao-Sah I-V
correctness reference at long-channel geometry, see
[`benchmarks/mosfet_2d`](../../benchmarks/mosfet_2d). For the
companion C-V analysis, see [`benchmarks/mos_2d`](../../benchmarks/mos_2d).

## Physics features exercised

- **M16.1 Caughey-Thomas bulk mobility.** `physics.mobility.bulk_model
  = "caughey_thomas"`. Velocity saturation in the channel under
  drain-side high field is moderate at the shipped V_DS = 0.05 V; the
  knob is on for users who clone this config and bump V_DS.
- **M16.2 Lombardi surface mobility.** `physics.mobility.model =
  "lombardi"` with `interface_facet_tag = 4` pointing at the gate
  facet. Surface acoustic-phonon and surface-roughness scattering in
  the inversion layer; the dominant ingredient for quantitative MOSFET
  I-V above threshold.
- **M16.4 Fermi-Dirac statistics.** `physics.statistics =
  "fermi_dirac"`. The Gaussian source/drain implants peak at
  N_D = 1e20 cm^-3, well above the Boltzmann break point
  (~1e19 cm^-3). The basic Blakemore approximation applies; see
  [`docs/IMPROVEMENT_GUIDE.md`](../../docs/IMPROVEMENT_GUIDE.md) § 4
  M16.4.
- **SRH recombination** with tau_n = tau_p = 10 ns. Auger and
  tunneling are intentionally off here so the run anchors the
  surface-mobility / FD path; they are exercised in
  [`examples/power_diode_reverse_recovery`](../power_diode_reverse_recovery)
  and the M16.6 `zener_1d` benchmark respectively.
- **M18.1 SNES nleqerr line search.** `solver.snes.line_search:
  "nleqerr"`. Required for this device under FD statistics; see
  [`docs/adr/0018-snes-line-search-nleqerr.md`](../../docs/adr/0018-snes-line-search-nleqerr.md)
  for the line-search decision and rejected alternatives.

## How to run

From the repository root:

```
docker compose run --rm benchmark nmos_idvgs
```

(The CLI in `scripts/run_benchmark.py` falls back to `examples/` when
the requested name is not found under `benchmarks/`, so the same
invocation works for both directories.)

Output lands in `results/nmos_idvgs/`:

- `id_vgs_overlay.png`: linear and semilog Id-Vgs curve.
- The result JSON written by the engine (with the per-step I-V table)
  is in the standard `results/nmos_idvgs/` location.

## Expected output

The anchor sweep at V_DS = 0.05 V shows monotonically increasing
J_drain through the subthreshold regime, with positive
transconductance dJ_drain / dV_GS > 0 over V_GS in [0.3, 0.5] V (the
inversion-onset window for this lighter body doping). The smoke
verifier gates exactly that qualitative shape; tight quantitative
matching against the Pao-Sah formula lives under
[`benchmarks/mosfet_2d`](../../benchmarks/mosfet_2d).

To extract V_T from the linear-regime curve, find the V_GS value at
which the linear extrapolation of the steepest-slope tangent crosses
zero current. With N_A = 1e16 cm^-3 body, V_T lands in roughly
[0.2, 0.4] V; the companion C-V geometry in
[`benchmarks/mos_2d`](../../benchmarks/mos_2d) gives the V_T from a
quasi-static gate sweep and is a useful cross-check.

## How to adapt

Most users will want to start by changing one of:

- **Body doping.** Change `doping[0].profile.N_A`. The shipped value
  N_A = 1e16 cm^-3 is the M18.1 reparam for SNES convergence under FD
  statistics. Heavier body doping (1e17 to 5e17) raises V_T and
  presently stalls the bias_sweep SNES at depletion-to-inversion onset
  for this short-channel geometry; that regime is tracked as a backlog
  item under "Bias-sweep SNES robustness, phase 2" (ADR 0018
  "Deferred"). For exploratory use at heavier doping, expect to halve
  the V_GS step and lower `solver.continuation.min_step` aggressively.
- **Channel length.** Edit the mesh extents in x and the Gaussian S/D
  centers proportionally. Shrinking the channel (sub-50 nm) introduces
  short-channel effects (DIBL, V_T roll-off) that the present
  Slotboom kernel captures only qualitatively.
- **Gate oxide thickness.** Edit the y extent (`extents[1][1]`) and
  the oxide region's lower y bound, then re-derive the mesh resolution
  so the oxide remains exactly one cell tall (the engine tags cells by
  centroid; if the oxide region's vertical extent contains no cell
  midpoint, the oxide effectively disappears).
- **S/D doping.** Change the Gaussian peak. Above 1e20 cm^-3 the
  Fermi-Dirac correction is essential; below 5e18 cm^-3 you can switch
  back to `physics.statistics = "boltzmann"` and the run will be
  meaningfully cheaper.
- **Temperature.** `physics.temperature` (default 300 K) is the
  thermal-budget knob. Above ~400 K you may want to add the Arora
  doping-dependent mobility branch (not yet shipped); below ~150 K the
  Boltzmann approximation breaks even more thoroughly.
- **Line search.** `solver.snes.line_search` accepts `bt`, `nleqerr`,
  `cp`, `l2`, `basic`. The shipped `nleqerr` is required for this
  device; users porting this config to a less-stiff device may
  benefit from leaving the default `bt` in place for faster solves.

## Notes on CI runtime

The shipped configuration runs in roughly 35 seconds on a typical GHA
runner: 6 V_GS samples at 0.1 V step, with the threshold-transition
step at V_GS = 0.1 V taking 22 SNES iterations under `nleqerr` (the
remaining 5 samples are 1-iter solves each). Total: 27 SNES iterations
across the full sweep.
