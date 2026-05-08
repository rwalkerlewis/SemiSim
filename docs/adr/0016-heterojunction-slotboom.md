# 0016. Heterojunction-aware Slotboom and ohmic equilibrium

- Status: Accepted
- Date: 2026-05-06

## Context

[ADR 0004](0004-slotboom-variables-for-dd.md) mandates Slotboom
variables `(psi, phi_n, phi_p)` as the primary unknowns of the
drift-diffusion block. Under single-material Boltzmann
statistics, the carrier density is recovered pointwise from

```
n = n_i exp((psi - phi_n) / V_t),
p = n_i exp((phi_p - psi) / V_t),
```

with the intrinsic density `n_i = sqrt(N_C N_V) exp(-Eg /
(2 V_t))` a scalar property of the reference material. The
continuity flux

```
J_n = -q mu_n n grad(phi_n)
```

is a pure gradient with a positive coefficient; standard
Galerkin on `phi_n` is stable without stabilization, which is
the load-bearing reason ADR 0004 is locked.

[ADR 0006](0006-verification-and-validation-strategy.md)
mandates that every new domain-physics module ship with an MMS
verifier whose finest-pair convergence rates clear `L^2 >= 1.99`
and `H^1 >= 0.99` on every meaningful block. Variants A through
H in `semi/verification/mms_dd.py` covered the residual stack
through M16.6.

[ADR 0007](0007-contact-bc-interface.md) defines the contact-BC
interface and the resolver in `semi/bcs.py`. The ohmic-
equilibrium psi calculation in that resolver today reads `chi`
implicitly from the reference material via
`make_scaling_from_config`.

M17 introduces position-dependent material parameters via
cellwise DG0 fields:

```
chi(x), Eg(x), N_C(x), N_V(x), n_i(x) = sqrt(N_C(x) N_V(x))
                                         exp(-Eg(x) / (2 V_t)).
```

At an AlGaAs/GaAs heterojunction, all of `chi`, `Eg`, `N_C`,
`N_V` jump at the interface, and `n_i(x)` therefore inherits a
step discontinuity. Three questions follow:

1. **Does Slotboom survive position-dependent `n_i(x)`?** Under
   the Boltzmann substitution `n = n_i(x) exp((psi - phi_n) /
   V_t)`, the carrier density `n` itself becomes discontinuous
   across the heterojunction. This is *physically correct*: the
   self-consistent Poisson solution at a heterojunction does
   exhibit a step in `n`, with the chemical potential `phi_n`
   continuous across the interface. The continuity flux
   `J_n = -q mu_n n grad(phi_n)` is well-defined because
   `grad(phi_n)` is the gradient of a C0 quantity and the
   discontinuity in `n` is absorbed in the algebraic
   substitution upstream of the divergence. The Slotboom
   advantage is preserved.
2. **Does ADR 0006 extend to the discontinuous-coefficient
   case?** Manufactured solutions for PDEs with discontinuous
   coefficients are a separate convergence-theory problem
   (interface-fitted vs unfitted meshes, jump conditions, dual-
   norm error estimators). The straightforward MMS construction
   that Variants A through H rely on does not extend cleanly to
   a discontinuous `n_i(x)`. The smooth-coefficient case (a
   continuous chi(x), Eg(x) ramp) does extend cleanly: it is
   Variant I.
3. **Where does the ohmic-contact equilibrium psi calculation
   go wrong if it keeps reading reference-material chi?** At a
   heterojunction with an ohmic contact in the wide-gap region
   and a second ohmic contact in the narrow-gap region, the
   Anderson rule gives a built-in band-edge offset of
   `chi_R - chi_L` between the two contacts. If the resolver
   computes `psi_eq` using the reference-material chi for both
   contacts, the offset is missed and the equilibrium psi
   profile is wrong by exactly that constant. This is silent
   under bit-identity (single-material configs collapse) but
   shows up as soon as a heterojunction is configured.

## Decision

M17 ships with the following V&V and BC structure:

1. **Generalized-Slotboom substitution.** The Slotboom helpers
   in `semi/physics/slotboom.py` accept a position-dependent
   `n_i_hat` argument (UFL Function or scalar; UFL handles
   both). The Poisson and DD form builders thread a per-cell
   DG0 `n_i_hat` field constructed by
   `semi/physics/heterojunction.py::build_dg0_material_fields`.
   The continuity-flux shape is unchanged. ADR 0004 is
   preserved. This is structurally analogous to the M16.4
   generalized-Slotboom path under Fermi-Dirac statistics
   (where the prefactor `gamma_n` enters the same substitution
   slot). Under FD-plus-heterojunction configurations the two
   extensions compose orthogonally: `gamma_n(x) * n_i(x)`.

2. **Ohmic-contact equilibrium psi reads local chi.** The
   resolver in `semi/bcs.py` looks up the material at each
   ohmic contact's region and computes

   ```
   psi_eq(R) = (chi(R) - chi_ref) / V_t
               + ln(N_doping / n_i(R))
   ```

   in dimensionless internal units. The single-material path
   (every region with the reference material, no
   `material_overrides`) collapses to the existing formula
   bit-identically because `chi(R) == chi_ref` and `n_i(R) ==
   n_i_ref` for every contact. The Schottky equilibrium psi
   helper from M16.5 already accepts a local-chi argument; the
   M17 patch generalizes the ohmic helper in the same shape.

3. **MMS Variant I covers the smooth-coefficient case.**
   `semi/verification/mms_dd.py` ships Variant I with
   `chi(x) = chi_0 + delta_chi * x/L` and
   `Eg(x) = Eg_0 + delta_Eg * x/L` chosen to give an O(0.1)
   shift in `n_i(x)` across the domain. The manufactured
   solution `psi_e, phi_n_e, phi_p_e` is unchanged from
   Variant C; only the substitution rule for `n` and `p`
   varies. Acceptance: finest-pair `L^2 >= 1.99` and `H^1 >=
   0.99` on each block. The textbook P1 rate gate applies
   because the coefficients are smooth; no decoupling caveat.
   The equilibrium-Poisson MMS in
   `semi/verification/mms_poisson.py` gains a counterpart
   smooth-ramp variant under the same gate.

4. **The discontinuous case is gated by the HEMT 2D
   benchmark.** `benchmarks/hemt_2d/` (AlGaAs/GaAs, V_GS sweep
   `[0.0, 1.0] V`) is the V&V departure for the
   discontinuous-coefficient case. The verifier compares the
   integrated 2DEG sheet density against a published
   Poisson-Schrodinger reference with a `15 %` tolerance for
   `V_GS in [0.4, 1.0] V`. The classical FEM solver matches
   the classical-electrostatic reference within `~5 %`; the
   `15 %` tolerance accommodates the gap to the fully-quantum
   Poisson-Schrodinger reference. This pair (Variant I plus
   `hemt_2d`) is the V&V structure for the M17 module:
   per-physics-module MMS for the smooth case, published-
   reference benchmark for the discontinuous case. Future
   heterojunction work (SiGe HBT, InGaAs/InP, quantum
   confinement corrections) reuses this shape.

## Consequences

Easier:

- M17 ships a real V&V gate without forcing a relaxation of
  the `L^2 >= 1.99 / H^1 >= 0.99` threshold. The
  smooth-coefficient case clears the textbook gate; the
  discontinuous case is gated by a physically-meaningful 2DEG
  match against a published reference.
- Future heterojunction-physics work has a documented
  precedent. The per-module MMS rule from ADR 0006 applies to
  smooth-coefficient extensions; the published-reference plus
  byte-identity pair is the V&V gate for discontinuous-
  coefficient extensions.
- The Slotboom helpers, Poisson form builder, DD form
  builder, and ohmic-BC resolver all extend without changing
  shape: the position-dependent path collapses to the scalar
  path bit-identically when the configuration is single-
  material. Bit-identity on every existing benchmark and PR
  #85 example is the regression gate.
- Composition with M16.4 Fermi-Dirac is orthogonal:
  `gamma_n(x) * n_i(x)` enters the same substitution slot as
  scalar `n_i` did in v0.21.0.

Harder:

- A future regression that breaks the heterojunction code at
  the discontinuous-coefficient boundary (interface jump
  conditions, integrated 2DEG sheet density) at a rate higher
  than the `15 %` `hemt_2d` tolerance will not be caught by
  the classical reference alone. This is acceptable because
  Variant I catches any rate-losing bug in the smooth-
  coefficient bulk discretization, and the `15 %` 2DEG gate
  is exponentially sensitive to the AlGaAs barrier height and
  the gate work function (any sign or coefficient bug in the
  band-edge alignment trips the gate immediately).
- ADR 0006 must be amended in spirit but not in letter:
  Variant I satisfies ADR 0006 for the smooth-coefficient
  case; the discontinuous case is documented as a deferral
  here rather than a direct ADR 0006 amendment, on the same
  precedent as ADR 0015 (Schottky as Robin BC). The
  per-physics-module MMS rule applies to **smooth-coefficient
  domain physics**; the published-reference plus byte-
  identity pair is the V&V structure for discontinuous-
  coefficient and boundary-physics modules.

## Cross-references

- Stern, "Self-consistent results for n-type Si inversion
  layers," Phys. Rev. B 5, 4891 (1972). Foundational treatment
  of self-consistent Poisson-Schrodinger at a 2DEG.
- Vurgaftman, Meyer, Ram-Mohan, "Band parameters for III-V
  compound semiconductors and their alloys," J. Appl. Phys.
  89, 5815 (2001). Source for the AlGaAs_0p3 band-edge
  parameters at 300 K.
- Sze, Physics of Semiconductor Devices, 3rd ed. § 1.5
  (carrier statistics, position-dependent `n_i(x)`) and
  § 2.10 (heterojunction band alignment, Anderson rule).
- [ADR 0004 Slotboom variables for drift-diffusion](0004-slotboom-variables-for-dd.md)
  (preserved by the generalized-Slotboom path documented
  here).
- [ADR 0006 Verification & Validation strategy](0006-verification-and-validation-strategy.md)
  (satisfied for the smooth-coefficient case via Variant I;
  partially deferred for the discontinuous case via the
  `hemt_2d` benchmark, on the same precedent as ADR 0015).
- [ADR 0007 Contact BC interface](0007-contact-bc-interface.md)
  (the ohmic resolver in `semi/bcs.py` extends to read local
  chi without changing the resolver's external interface).
- [ADR 0015 Schottky contacts as Robin BCs](0015-schottky-robin-bc.md)
  (precedent for V&V departures from the strict ADR 0006 MMS
  gate when the physics class does not admit a clean MMS
  variant).
- `docs/IMPROVEMENT_GUIDE.md` § M17 (acceptance criteria).
- `docs/PHYSICS.md` § 1.1, § 1.2, § 1.3 (electrostatics,
  carrier statistics under position-dependent material
  parameters, continuity in Slotboom current form at a
  heterojunction).
- `benchmarks/hemt_2d/README.md` (HEMT device description and
  verifier).
