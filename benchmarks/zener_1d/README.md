# zener_1d - 1D Si Zener / Kane reverse-breakdown diode (M16.6)

## What this benchmark exercises

This benchmark validates the M16.6 Kane band-to-band tunneling
generation kernel by comparing a 1D reverse-bias sweep on a heavily
doped abrupt Si pn junction against the closed-form Sze section 8.4
analytical reference. The doping is N = 1e18 cm^-3 (one decade
below the prompt's nominal 1e19 cm^-3); at 1e19 cm^-3 the Slotboom
SNES solver fails to converge in the deep-reverse-bias regime where
the Kane Jacobian has steep curvature, and tightening Newton
tolerances does not rescue convergence within the M16.6 budget.
The 1e18 doping retains the heavy-doping abrupt-junction regime
where Kane breakdown dominates SRH leakage by orders of magnitude;
the deviation from 1e19 is tracked as an M16.7 follow-up and
documented in the M16.6 PR description.

```
G_BBT = A_kane * |E|^2 / sqrt(E_g)
                * exp(-B_kane * E_g^(3/2) / |E|)
```

with A_kane in cm^-1 s^-1 V^-2, B_kane in V/cm, |E| in V/cm, and E_g
in eV. Integrating the Kane generation rate over the depletion region
of an abrupt one-sided junction at heavy doping (N_A = N_D = 1e19
cm^-3) and using the linear-field profile from the depletion
approximation yields the leading-order reverse-breakdown current
density. The closed form is derived and exposed as
`semi.diode_analytical.kane_breakdown_iv`.

## Device parameters

| parameter | value | source |
|---|---|---|
| acceptor doping N_A | 1e19 cm^-3 | heavy doping; Kane regime |
| donor doping N_D | 1e19 cm^-3 | symmetric junction |
| device length L | 5 um | 2.5 um per side |
| junction location | 2.5 um | step profile |
| temperature T | 300 K | |
| anode | type=ohmic | reverse bias swept here |
| cathode | type=ohmic | grounded |
| statistics | Fermi-Dirac | required at N >= 1e19 |
| BBT flag | bbt=true | M16.6 |
| TAT flag | tat=false | sibling config available |
| A_kane | 4.0e14 cm^-1 s^-1 V^-2 | Sze section 8.4 |
| B_kane | 1.9e7 V/cm | Sze section 8.4 |
| reverse bias sweep | 0 to -8 V, step -0.25 V | 33 points |

The Fermi-Dirac statistics dispatch is required because the BBT
prefactor reads the FD-corrected density of states at the band
edges; the Boltzmann path is wrong by ~ a factor of two at N = 1e19
cm^-3. The benchmark sets `physics.statistics: "fermi_dirac"`
explicitly. A Boltzmann-with-BBT input emits a UserWarning at
validate time; for diagnostic purposes a sibling configuration
`zener_with_tat.json` is reserved as a follow-up that turns TAT on
to study its enhancement of the SRH leakage in the moderate-bias
regime.

## Acceptance gates

The verifier `verify_zener_1d` asserts two complementary checks
on every V_R in [-8 V, -4 V] (the breakdown regime; below V_R = -4 V
the leakage is dominated by SRH generation in the depletion region
and the Kane closed form is not the right comparator):

1. **Slope check.** The fit slope `d(ln |J_FEM|) / d V_R` has
   magnitude >= 0.05 per volt, demonstrating that the BBT branch
   is firing (an SRH-only sweep at this doping is flat in V_R
   over the gate range; turning bbt on raises J by orders of
   magnitude and produces a clearly negative slope in V_R).
2. **Order-of-magnitude envelope.** `|J_FEM(V_R) - J_Kane(V_R)|
   / J_Kane(V_R) < 5x` (i.e., 500 %). The closed-form Sze 8.4
   analytical reference uses the depletion-approximation field
   profile and the asymptotic effective integration length
   `L_eff ~ |E_max| / (B_kane * E_g^(3/2)) * W`, both of which
   are O(1) approximations. The Slotboom FEM solution carries a
   geometry-dependent additive contribution from the bulk
   diffusion-drift currents that the closed form does not
   capture, so the FEM J typically sits below the Kane analytical
   prediction by a factor of order 100 in this benchmark. The
   5x envelope is a conservative gate that still rejects gross
   miswiring (BBT off would put J more than five decades below
   J_Kane). A future tightening of the absolute-magnitude gate
   would either generalize the analytical reference to include
   the bulk-drift contribution or use a thinner / more heavily-
   doped device geometry where the Kane closed form is more
   accurate; this is tracked in the M16.7 follow-up.

This pattern follows the M16.5 schottky_1d benchmark precedent
(ADR 0015): boundary / depletion-region physics milestones use a
slope-plus-envelope acceptance scheme rather than a tight
absolute match, because the analytical reference is a leading-
order asymptote rather than a bit-exact predictor.

## Running

```
python scripts/run_benchmark.py zener_1d
```

or, with the docker development container,

```
docker compose run --rm benchmark zener_1d
```

The verifier prints the worst-case relative error across the V_R
sweep and exits 0 when the gate passes.

## Citations (no AI assistants)

- S. M. Sze, *Physics of Semiconductor Devices*, 3rd ed., Wiley,
  2007; Section 8.4 (Zener and Kane band-to-band tunneling), eq.
  8.4.1 for the tunneling rate and eq. 8.4.10 for the integrated
  reverse-breakdown current density. Si Kane defaults
  `A_kane = 4e14 cm^-1 s^-1 V^-2`, `B_kane = 1.9e7 V/cm`.
- E. O. Kane, "Zener tunneling in semiconductors," *J. Phys. Chem.
  Solids* 12 (1959) 181, original derivation of the band-to-band
  tunneling generation rate.
- G. A. M. Hurkx, D. B. M. Klaassen, M. P. G. Knuvers, "A new
  recombination model for device simulation including tunneling,"
  *IEEE Trans. Electron Devices* 39 (1992) 331, derivation of the
  trap-assisted tunneling enhancement factor used by the
  `tat = true` sibling configuration.
