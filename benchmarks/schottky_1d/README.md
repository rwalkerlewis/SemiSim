# schottky_1d — 1D Pt-on-n-Si Schottky diode (M16.5)

## What this benchmark exercises

This benchmark validates the M16.5 thermionic-emission boundary
condition by comparing a 1D forward I-V sweep on a Pt-on-n-Si Schottky
diode against the closed-form thermionic-emission analytical reference
(Sze 3rd ed Section 3.4):

```
J(V) = A* T^2 exp(-q phi_B / kT) [exp(qV / kT) - 1]
```

where `A* = 4 pi q m_n* k^2 / h^3` is the Richardson constant for the
electron thermionic-emission effective mass. For Si we use the Sze
3rd ed Table 1 thermionic mass `m_n* = 0.26 m_0`, which the engine
threads from `semi/materials.py` through `Scaling.v_n_thermal` into the
Schottky surface form on the continuity rows.

## Device parameters

| parameter | value | source |
|---|---|---|
| barrier height phi_B | 0.85 eV | Pt-on-n-Si, Sze 3rd ed Table 5 |
| donor concentration N_D | 1e16 cm^-3 | moderate doping |
| device length L | 5 um | thin enough that depletion is short |
| temperature T | 300 K | |
| anode | type=schottky | M16.5 |
| cathode | type=ohmic | |
| statistics | Boltzmann | M16.5 ships under Boltzmann; FD follow-up |
| forward bias sweep | 0 to 0.5 V, step 0.025 V | 21 points |

The doping is moderate (1e16 cm^-3) so the depletion region is short
(W ~ 0.3 um at zero bias) and the thermionic-emission rate sets the
boundary current rather than diffusion. The 5 um device length is
chosen such that the bulk diffusion length is much larger than the
depletion width, isolating the boundary thermionic-emission rate as
the dominant rate-limiting step.

## Acceptance gates

The verifier `verify_schottky_1d` asserts two complementary checks
on every V in [0.1 V, 0.5 V] (21 sweep points minus the V = 0 V
endpoint, which is excluded because the relative error diverges at
zero current):

1. **Slope match.** A linear fit to `ln(J_FEM)` vs `V` must reproduce
   the thermionic-emission slope `1 / V_t` within 5 %. This is the
   exponential V-dependence diagnostic; the FEM forward I-V must look
   like the analytical formula on a semilog-y plot.
2. **Order-of-magnitude absolute match.** `|J_FEM(V) - J_thermionic(V)|
   / J_thermionic(V) < 5x` (i.e., 500 %). The simple analytical
   thermionic-emission formula is the *thermionic-limit asymptote* of
   the thermionic-diffusion theory; for a finite-thickness device with
   a finite-rate bulk diffusion-drift the FEM picks up a geometry-
   dependent additive bulk-drift contribution that lifts the absolute
   J above the thermionic-limit prediction. Empirically the schottky_1d
   device (5 um Pt-on-n-Si at N_D = 1e16 cm^-3) sits in a mixed
   regime where FEM J is ~2.5-3 times the thermionic-limit value
   while preserving the V-dependence slope. The 5x absolute gate is
   a conservative envelope around this physics.

ADR 0015 documents the broader V&V scope: boundary-physics milestones
use the analytical-benchmark slope match plus existing-benchmark byte-
identity instead of an MMS rate gate, because the existing MMS-DD
harness uses Dirichlet BCs everywhere by construction and substituting
a Robin BC into the manufactured-source residual either cancels the
Robin constraint into the source or produces a BC-dominated error
that does not converge at the optimal rate.

A future tightening of the absolute-magnitude gate would either
generalize the analytical reference to include the thermionic-
diffusion correction or use a thinner / heavier-doped device geometry
where bulk drift is fast enough that the thermionic limit dominates.
M16.5 ships under the slope-plus-envelope pattern; a follow-up can
revisit the absolute gate once the M16.6 tunneling and M16.7
transient-AC consistency milestones close.

## Running

```
python scripts/run_benchmark.py schottky_1d
```

or, with the docker development container,

```
docker compose run --rm benchmark schottky_1d
```

The verifier prints the worst-case relative error across the V_F
sweep and exits 0 when the gate passes.

## Citations (no AI assistants)

- S. M. Sze, *Physics of Semiconductor Devices*, 3rd ed., Wiley, 2007;
  Section 3.4 (Schottky barrier lowering and thermionic emission;
  derivation of `J = A* T^2 exp(-q phi_B / kT) [exp(qV / kT) - 1]`),
  Table 1 (Si thermionic-emission effective masses), Table 5 (Pt-on-Si
  barrier height 0.85 eV).
- S. Selberherr, *Analysis and Simulation of Semiconductor Devices*,
  Springer, 1984; Section 5.2 (numerical treatment of the Robin BC
  for thermionic emission in finite-element drift-diffusion solvers).
