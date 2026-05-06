# 0015. Schottky contacts as Robin BCs (V&V departure from ADR 0006)

- Status: Accepted
- Date: 2026-05-05

## Context

[ADR 0006](0006-verification-and-validation-strategy.md) (V&V
strategy) requires that every new physics module ship with an MMS
verifier whose finest-pair convergence rates clear `L^2 >= 1.99`
and `H^1 >= 0.99` on every meaningful block. M16.1
(Caughey-Thomas mobility), M16.2 (Lombardi surface mobility),
M16.3 (Auger recombination), and M16.4 (Fermi-Dirac statistics)
each landed an MMS variant in `semi/verification/mms_dd.py`
(Variants D, E, F, G respectively) that clears that gate.

M16.5 adds **thermionic-emission Robin boundary conditions** at
metal-semiconductor (Schottky) contacts:

```
J_n . n = q v_n_thermal (n - n_eq exp((V_a - phi_n) / V_t))
J_p . n = -q v_p_thermal (p - p_eq exp((phi_p - V_a) / V_t))
```

This is a **boundary** physics change, not a domain residual
change. The MMS-DD harness in `semi/verification/mms_dd.py` was
designed against the assumption that every contact is a Dirichlet
boundary (the existing Variants A through G all impose Dirichlet
on `psi`, `phi_n`, and `phi_p` at every contact and then
manufacture a smooth interior solution that satisfies the chosen
boundary trace). Substituting a Robin condition into that
machinery has two failure modes:

1. The Robin constraint
   `J . n = alpha (u - u_inf)` with manufactured `u_inf` derived
   from the chosen interior solution **cancels into the
   manufactured source term**: the residual is then identical to
   the Variant A through G residual with a slightly different
   source, so the rate gate measures the bulk DD discretization,
   not the BC physics. The gate becomes uninformative.
2. Manufacturing a smooth interior solution that **does not
   satisfy** the Robin constraint (so the constraint is
   non-trivially active) produces a BC-dominated error layer at
   the Schottky facet, which converges at a sub-optimal rate
   under continuous P1 elements (the trace error is
   `O(h^{1/2})` in `H^1`, not `O(h)`). The gate will fail at
   `L^2 >= 1.99` or `H^1 >= 0.99` on every mesh, which is the
   correct physical answer but a useless test outcome: every
   future Robin-BC change would have to relax the threshold,
   and the gate would no longer detect real regressions.

There is no clean third path. Robin-BC convergence theory at the
optimal rate requires either curved-element correction at the
boundary or a custom dual-norm error estimator, both of which are
outside the scope of M16 and would substantially complicate
`semi/verification/mms_dd.py`.

## Decision

**M16.5 ships without an MMS Variant H.** Its V&V gate is the pair:

1. **Analytical-benchmark match.** The new
   `benchmarks/schottky_1d/` device (1D Pt-on-n-Si Schottky diode,
   `N_D = 1e16 cm^-3`, `V_F` sweep `[0, 0.5] V`) must match the
   closed-form thermionic-emission I-V (per Sze 3rd ed § 3.4 and
   Selberherr 1984 § 5) within 10 % from `V_F = 0.1 V` to
   `V_F = 0.5 V`. Below `V_F = 0.1 V` the relative error blows
   up at the zero-current limit; above `V_F = 0.5 V` the
   diffusion-limited regime takes over and the
   thermionic-emission analytical reference is no longer the
   right comparator.
2. **Existing-benchmark byte-identity.** Every benchmark in the
   repo today (none use Schottky contacts) must produce results
   bit-identical to v0.20.0. The four anchors that gate this in
   CI are the same anchors M16.1 through M16.4 carry forward:
   `pn_1d_bias` `J(V=0.6 V) = 1.635e+03 A/m^2`,
   `diode_velsat_1d` `56.27 %` @ `V_F = 0.9 V` and `0.19 %` @
   `V_F = 0.3 V`, `diode_auger_1d` `>20 %` SRH-vs-(SRH+Auger)
   divergence at `V_F = 0.9 V`, and `diode_fermi_dirac_1d`
   `7.37 %` Boltzmann-vs-FD V_bi divergence at
   `N_D = 1e20 cm^-3`.

This pair (analytical-benchmark plus byte-identity) is the
**boundary-physics analog** of the per-module MMS rate gate that
ADR 0006 requires for domain-physics modules. The two halves of
the pair guarantee separately:

- The analytical-benchmark match catches sign or coefficient
  errors in the new Schottky form (an MMS rate gate would also
  catch these, but the absolute-current comparison is the
  stronger physical statement).
- The byte-identity match catches accidental coupling into the
  ohmic / gate / insulating paths (an MMS rate gate cannot see
  this because each variant only exercises its own physics
  block).

Future BC-adding milestones (Schottky variants under FD
statistics, ohmic-contact recombination boundary terms, traps at
the metal-semiconductor interface) follow the same convention.
The per-module MMS rule in ADR 0006 applies to **domain physics
modules**: changes to the bulk PDE residual, source terms, or
constitutive relations. The analytical-benchmark plus
byte-identity pair is the V&V gate for **boundary physics**.

## Consequences

Easier:

- The M16.5 PR ships a real V&V gate without forcing a relaxation
  of the `L^2 >= 1.99 / H^1 >= 0.99` threshold that domain-
  physics modules clear cleanly.
- Future BC-physics work (M16.6 trap-assisted boundary
  recombination, future heterojunction interface conditions
  under M17) has a documented precedent for the V&V shape; no
  per-PR re-litigation of "do we need an MMS variant for this?".
- The `benchmarks/schottky_1d/` verifier is more physically
  meaningful as a Schottky regression check than any MMS variant
  could be, because the closed-form thermionic-emission I-V is
  the textbook result every device engineer recognizes.

Harder:

- A future regression that breaks the Schottky form at a rate
  higher than the analytical-benchmark tolerance (10 %) will
  not be caught by the schottky_1d verifier. This is acceptable
  because:
  (a) any rate-losing bug in the bulk DD discretization will
      still trip Variants A through G;
  (b) any sign or coefficient bug in the Schottky form will
      blow the 10 % gate immediately, since the thermionic-
      emission current is exponentially sensitive to barrier
      height and applied voltage;
  (c) the byte-identity gate catches accidental coupling into
      the existing contact paths.
  The residual risk is a Schottky-only regression that preserves
  the order of magnitude but loses (say) a factor of two; the
  10 % gate flags this trivially.
- ADR 0006 must be amended to reference this ADR. Phase F of
  M16.5 adds a one-line footnote to ADR 0006 § Decision item 4
  pointing here.

## Cross-references

- Sze, Physics of Semiconductor Devices, 3rd ed., § 3.4
  ("Schottky barriers and Ohmic contacts"). Derivation of the
  thermionic-emission I-V; Table 5 gives Pt-on-n-Si
  `phi_B = 0.85 eV` and Si Richardson constants (electron
  `A* = 110 A / cm^2 / K^2`, hole `32 A / cm^2 / K^2`).
- Selberherr, Analysis and Simulation of Semiconductor Devices,
  Springer 1984, § 5.2 ("Boundary conditions for the carrier
  continuity equations"). Standard numerical treatment of the
  Robin BC as a surface integral on the continuity rows.
- [ADR 0006 V&V strategy](0006-verification-and-validation-strategy.md)
  (amended in M16.5 Phase F with the footnote that points here).
- [ADR 0007 Contact BC interface](0007-contact-bc-interface.md)
  (the resolver in `semi/bcs.py` already validates the
  `schottky` kind; M16.5 adds the dedicated builder dispatch
  this ADR anticipated).
- `docs/IMPROVEMENT_GUIDE.md` § M16.5 (acceptance criteria).
- `benchmarks/schottky_1d/README.md` (benchmark device and
  verifier description).
