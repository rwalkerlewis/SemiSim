## Context

You are working in `kronos-semi` at `v0.23.0` (post-M16.7, post-
examples-catalogue PR #85, package version `0.23.0`, schema
`2.7.0`). The repo is at
`https://github.com/rwalkerlewis/kronos-semi`; `main` is at
`3a23383` (`examples: catalogue of practical device configs (#85)`).

Your assignment is **M17: Heterojunction / position-dependent
band structure**. The maintainer has chosen M17 over M19 (3D
MOSFET capstone); M19 remains in the roadmap as a follow-up
milestone.

M17 is the largest single physics addition since M14.3. It is
substantially heavier than any single M16.x slice. Reasons:

- It promotes `n_i`, `N_C`, `N_V`, `Eg`, `chi` from scalar
  fields on `Scaling` and the reference material to *position-
  dependent cellwise DG0 fields* on the mesh. Every
  callsite that consumes these as scalars needs a position-
  aware path.
- It extends the materials model with `AlGaAs` (with a chosen
  Al fraction, ~0.3 for the HEMT benchmark) and per-region
  material overrides via the schema.
- It touches the Poisson form, the DD form, the Slotboom
  helpers, and the ohmic-contact equilibrium psi calculation.
  All four are bit-identity-critical against v0.23.0 when
  the configuration is single-material.
- It adds an MMS variant (smooth chi/Eg ramp; no discontinuity)
  and an HEMT 2D benchmark (discontinuous; the analytical
  reference is a published Poisson-Schrodinger solution within
  a 15 % quantum-confinement caveat).
- It requires a new ADR (0016) extending the M16.4 generalized-
  Slotboom argument to position-varying `n_i(x)` and the
  ohmic-contact equilibrium psi adjustment.

`PLAN.md` § "Next task" names M17 as one of two candidates.
The maintainer has selected M17.

> **M17: heterojunctions.** Position-dependent electron
> affinity chi(x) and band gap E_g(x) on multi-region meshes.
> Touches `semi/physics/poisson.py` and
> `semi/physics/slotboom.py` (the band-edge offsets are no
> longer constants). Depends on M16.4 (Fermi-Dirac, needed at
> heterojunction barriers); unblocked.

This prompt **is** the starter. Save it as
`docs/M17_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below), preserving the M-numbered convention from
M16.x.

## Branch and PR rules

- Work on a fresh branch off `main`:
  `git checkout -b dev/m17-heterojunction`. Do **not** rebase
  or commit onto `main` directly.
- One milestone, one PR. Do not bundle M17 with M19 (3D MOSFET)
  or any other physics work.
- Push every phase commit to `origin/dev/m17-heterojunction`
  immediately after it lands locally. Open the PR after Phase 0
  so reviewers can watch phases land.
- Title the PR `M17: Heterojunction / position-dependent band
  structure`. Use the PR description template at the bottom.
- **Run all seven phases consecutively without pausing for
  confirmation, status reports, or "should I continue"
  prompts.** Do not call `ask_user_input` or any equivalent
  between phases. The phases are an execution plan, not a
  checklist of separately-permissioned tasks. Run Phase 0
  through Phase F end to end. Push commits as they land.
  The only legitimate reasons to stop are: (a) a hard blocker
  that re-reading this prompt, the linked ADRs, or the codebase
  cannot resolve; (b) a `pytest` or `ruff` failure that
  persists after a reasonable repair attempt; (c) the Stop
  conditions section reports green and the PR is ready for
  review. "Should I move on to the next phase" and "do you
  want me to also add X" are not legitimate reasons. Default
  behavior: pick the most defensible interpretation, document
  inline, continue.
- **No AI-assistant credits in any shipped artifact.** The
  `includeCoAuthoredBy: false` settings fix is presumed in
  effect (it held through M16.6, M16.7, and PR #85); verify
  with `git log --format=%B dev/m17-heterojunction | grep -i
  'co-authored-by: claude'` before opening and merging the PR.
  Empty output is the gate.

## Required reading (in order; ~75 minutes; this milestone is
the heaviest M-numbered slice to date)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.23.0
   (post-M16.7 + post PR #85) and that no other PR is in
   flight on `dev/m17-*`. Read the M16.4 entry in "Completed
   work log" carefully (the generalized-Slotboom path under
   FD is the structural prior art for M17's position-dependent
   `n_i(x)` substitution).
2. `docs/IMPROVEMENT_GUIDE.md` § M16.4 (the structural
   precedent for "extend the Slotboom substitution rule
   without breaking ADR 0004"), § M17 (Why / Deliverable /
   Acceptance / Dependencies), and § 1 (Honest current state).
3. `docs/ROADMAP.md` § Capability matrix. M17 is currently
   Planned; the row moves to shipped at the end of this PR.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M17 touches
   the schema (Layer 2; new `material_overrides` and
   `heterojunction` keys), pure-Python material model
   extensions (Layer 3), Poisson and DD form builders
   (Layer 4), the BC layer (Layer 4; ohmic equilibrium psi),
   the verification harness (Layer 5; MMS Variant I), and the
   benchmark (Layer 5; HEMT 2D).
5. `docs/PHYSICS.md`:
   - § 1.1 (electrostatics; the heterojunction makes
     `eps_r(x)`, `chi(x)`, `Eg(x)` position-dependent in the
     Poisson source).
   - § 1.2 (carrier statistics; `n_i(x) = sqrt(N_C(x) N_V(x))
     exp(-Eg(x) / (2 V_t))`. Both terms inherit position-
     dependence at a heterojunction).
   - § 1.3 (continuity in Slotboom current form). The
     continuity flux `J = -q mu n grad(phi_n)` retains its
     shape at a heterojunction; only the substitution rule
     for n in terms of phi_n changes.
   - § Verification & Validation. Variants A through H live
     in `semi/verification/mms_dd.py`. M17 adds **Variant I**
     (smooth chi(x), Eg(x) ramp; no discontinuity).
6. `docs/adr/`:
   - 0001 (JSON contract), 0002 (nondimensionalization), 0004
     (Slotboom DD; preserved through M17). 0006 (V&V
     strategy; per-physics-module MMS rule applies to M17).
   - 0007 (BC interface; the ohmic equilibrium psi
     calculation in `semi/bcs.py` needs a one-line edit to
     read chi from the *local* material rather than the
     reference material; document via the new ADR 0016).
   - 0015 (Schottky as Robin BC; M17 does not interact with
     this).
7. `semi/materials.py` end to end. Today: Si, Ge, GaAs,
   SiO2, HfO2, Si3N4. AlGaAs (Al_0.3 Ga_0.7 As at room
   temperature) is added in Phase B with explicit Vurgaftman-
   2001 derived parameters.
8. `semi/scaling.py` end to end. The `Scaling` dataclass
   carries `n_i: float, N_C: float | None, N_V: float | None`
   today. Under M17 these become *reference values* that the
   Poisson and DD form builders multiply by per-region
   correction factors derived from per-cell DG0 chi, Eg, Nc,
   Nv fields. The reference material is unchanged.
9. `semi/physics/poisson.py` L70-L78 (the `rho_hat = ni_hat *
   (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn` inline). The
   substitution under M17 becomes `ni_hat(x) * (...)` where
   `ni_hat(x)` is a UFL expression evaluated against the per-
   cell n_i field built from the per-region material. The
   multi-region variant `_mr` at L100-L120 has the same
   structure.
10. `semi/physics/slotboom.py` end to end. The `n_from_slotboom`
    and `p_from_slotboom` helpers gain a position-dependent
    `n_i_hat` argument. Today the call sites pass a scalar
    `Scaling.n_i / sc.C0`; under M17 they pass a UFL expression
    that picks up the cellwise DG0 n_i field. The numpy
    counterparts get the same treatment for the verifier path.
11. `semi/physics/drift_diffusion.py` L150-L160 (the
    `n_hat = n_from_slotboom(psi, phi_n, ni_hat)` and the
    `R_SRH` / `R_Auger` blocks). The `ni_hat` argument
    becomes position-dependent.
12. `semi/bcs.py` `_ohmic_psi_eq` (or wherever the equilibrium
    psi at an ohmic contact is computed today). Today it
    reads `chi` from the reference material. Under M17 it
    reads `chi` from the *local* material at the contact
    facet's cell. Same shape applies to gate-contact psi
    where work-function compensation against chi shows up
    (M16.5 added `_schottky_psi_eq` with similar local-chi
    sensitivity; cross-reference that for the pattern).
13. `semi/verification/mms_dd.py` for the Variant H structural
    template (M16.6). Variant I follows the same shape but
    manufactures smooth chi(x), Eg(x), Nc(x), Nv(x) ramps and
    propagates them through the manufactured weak source.
14. `semi/verification/mms_poisson.py` for the equilibrium-
    Poisson MMS coverage; the chi/Eg position-dependence
    needs an MMS check on the equilibrium Poisson form too
    (the manufactured rho includes the position-varying
    intrinsic source).
15. `benchmarks/moscap_axisym_2d/` and `benchmarks/mos_2d/`
    as the closest existing two-region benchmarks (silicon +
    SiO2). Structurally close to the HEMT 2D in geometry but
    SiO2 is an insulator (no heterojunction in the
    semiconductor sense; the existing infrastructure already
    handles material discontinuity for the Poisson permittivity
    but not for chi / Eg / n_i).
16. **Stern & Sarma 1984 or any equivalent published 2DEG
    sheet-density vs gate voltage reference** for the AlGaAs/
    GaAs HEMT. The IMPROVEMENT_GUIDE acceptance specifies
    "within 15 % of a published self-consistent Poisson-
    Schrodinger reference"; the 15 % tolerance is generous
    because the FEM solver here is purely classical (no
    quantum confinement). Any reference giving n_s vs V_GS
    in the linear range works; cite Pozela &
    Reklaitis 1980, Stern 1972, or a Sentaurus / Atlas
    validation paper. Do **not** cite ChatGPT or any AI tool.

## Heterojunction interaction (read this twice before Phase B)

ADR 0004 is the load-bearing constraint for M17 the same way
it was for M16.4. The natural question "does a heterojunction
invalidate Slotboom?" has the same published answer: no, but
the substitution rule extends. Under single-material
Boltzmann statistics:

```
n = n_i * exp((psi - phi_n) / V_t)
```

with `n_i` a scalar. At a heterojunction, `n_i` becomes
position-dependent because:

```
n_i(x) = sqrt(N_C(x) * N_V(x)) * exp(-Eg(x) / (2 V_t))
```

and N_C(x), N_V(x), Eg(x) all jump at the material interface.
The generalized-Slotboom form is:

```
n = n_i(x) * exp((psi - phi_n) / V_t)
```

with `n_i(x)` a DG0 cellwise field. Across the
heterojunction, n itself has a step discontinuity. This is
*physically correct* — the carrier density does step at a
heterojunction in a self-consistent classical Poisson
solution. The Slotboom advantage holds: phi_n is the primary
unknown and is C0-continuous (the FE space ensures it); n is
derived from phi_n and chi and is allowed to be discontinuous.
The continuity flux `J = -q mu n grad(phi_n)` remains valid
because grad(phi_n) is the gradient of a continuous field
and the discontinuity in n cancels with the discontinuity in
chi when the algebraic substitution is propagated through
correctly.

**Implementation consequence.** The Slotboom helpers
(`n_from_slotboom`, `p_from_slotboom`) accept a position-
dependent `n_i_hat` argument (the same shape as a UFL
Function or a fem.Constant; UFL handles both). The Poisson
form builder builds `n_i_hat` from a per-cell DG0 field
constructed at runner-setup time from the per-region
material. The DD form builder receives the same `n_i_hat`
and threads it through the Slotboom helpers. **Do not
rewrite the residual.** The shape is unchanged; only the
substitution rule for n is.

If the configuration is single-material (every region has
the same material with no `material_overrides`), the per-
cell n_i field is a constant; the position-dependent path
collapses to the single-material path bit-identically. This
is acceptance test 2.

The ohmic-contact equilibrium psi calculation reads chi from
the *local* material at the contact facet, not the reference
material. Today the calculation uses `Scaling.chi_ref`
implicitly (the chi of the reference material was used to
compute n_i). Under M17, the equilibrium psi at an ohmic
contact in region R is:

```
psi_eq(R) = chi(R) - chi_ref + V_t * ln(N_doping / n_i(R))
```

Where chi(R), n_i(R) are the local-material values. The
existing single-material formula collapses when chi(R) ==
chi_ref and n_i(R) == n_i_ref.

If FD is also enabled (M16.4), the FD-correction prefactor
`gamma_n` from M16.4 multiplies into n_i(x) in the same
substitution rule. The two extensions compose orthogonally.
Cross-check the algebra against PHYSICS.md § 1.3 before
committing Phase B; if the cancellation does not work out
cleanly under your derivation, stop and open a new ADR
extension before continuing. ADR 0016 documents the M17 path
and is shipped in Phase 0; if you find a derivation issue,
update ADR 0016 before any code change.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** The new schema entries
  `regions[].material` (already exists; reuse), `regions[].
  material_overrides: {chi_eV, Eg_eV}` (new), and a
  `heterojunction: true` flag on a region (new) must be
  expressible in `schemas/input.v2.json` (current strict
  v2.7.0), validated by `semi/schema.py`, and exercised by
  the HEMT benchmark.
- **Schema versioning is binding.** Additive minor bump
  v2.7.0 to v2.8.0. v2.0.0 through v2.7.0 inputs must
  continue to validate and produce bit-identical results to
  v0.23.0 (no input today specifies a heterojunction).
  Update `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Five layers, enforced.**
  - `semi/materials.py` (Layer 3) gains the AlGaAs entry
    plus a per-material `n_i_at_T(T)` helper. Pure-Python.
    No dolfinx imports.
  - `semi/scaling.py` (Layer 3) gains a per-region
    `material_at(region_tag) -> Material` lookup table on
    `Scaling`. Pure-Python.
  - The DG0 cellwise field construction lives in a new
    `semi/physics/heterojunction.py` (Layer 4 with deferred
    UFL imports; the field-construction helpers do
    `import dolfinx.fem` and `import ufl` inside function
    bodies).
  - The Poisson and DD form-builder edits stay in their
    existing files (Layer 4); imports already deferred.
- **Slotboom variables for DD.** ADR 0004 is locked. The
  generalized-Slotboom path (see "Heterojunction interaction"
  above) is the supported route. Do not introduce SUPG,
  streamline diffusion, or a primary-density form. If the
  cancellation argument breaks down under your derivation,
  stop and update ADR 0016 before any code change.
- **Per-physics-module MMS rule applies.** ADR 0006 mandates
  Variant I for M17. Acceptance L2 rate >= 1.99 / H1 rate
  >= 0.99 at the finest pair on each block, gated in
  `scripts/run_verification.py mms_dd`. **The discontinuous-
  heterojunction case is not in scope for the MMS variant**
  (manufactured solutions for discontinuous-coefficient PDEs
  have their own ADR-style debate that is not entered here);
  Variant I uses smooth chi(x), Eg(x) ramps. The discontinuous
  case is exercised by the HEMT benchmark, not MMS.
- **Per-runner threading.** Every runner that constructs
  Poisson or DD form (`equilibrium`, `bias_sweep`, `transient`,
  `ac_sweep`, `mos_cv`, `mos_cap_ac`, `resistor_3d`)
  constructs the per-cell DG0 chi, Eg, Nc, Nv, n_i fields at
  setup time and passes them to the form builders. **Seven
  runners** like M16.4. The construction is centralized in
  `semi/physics/heterojunction.py::build_dg0_material_fields`;
  each runner calls it once. Bit-identity is preserved when
  the configuration is single-material because the constructed
  field is a constant.
- **Bit-identity quantitative anchors.** Every existing
  benchmark and every PR #85 example must be bit-identical to
  v0.23.0 when no heterojunction is configured. Anchors:
  - `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2`
  - `diode_velsat_1d` 56.27 % @ V_F = 0.9 V, 0.19 % @
    V_F = 0.3 V
  - `diode_auger_1d` >20 % SRH-vs-(SRH+Auger) divergence at
    V_F = 0.9 V
  - `diode_fermi_dirac_1d` 7.37 % FD-vs-Boltzmann V_bi at
    N_D = 1e20 cm^-3
  - `schottky_1d` worst-case <10 % thermionic match
  - `zener_1d` worst-case <20 % Kane match
  - PR #85 example smoke gates green
- **One milestone, one PR.** Do not bundle M17 with M19 (3D
  MOSFET) or M19.1 (MPI) or M20 (HTTP server). Do not retire
  the `mosfet_2d` `allow-failure: "true"` flag (carry-over
  from M16.1-M16.7).
- **No em dashes in prose or code comments.**
- **Coverage gate is 95 in the gated `docker-fem-tests` job.**
  M16.5 needed a coverage-recovery follow-up commit; M16.6
  needed two; M16.7 needed none; the examples PRs needed
  none. M17 form-builder edits touch ~7 files in `semi/`;
  the HEMT benchmark runs in a separate CI matrix entry
  whose coverage is not merged. **Plan unit tests before
  Phase E** (the HEMT benchmark) so the new form-builder
  branches are covered by the gated suite. Specifically: a
  `tests/fem/test_heterojunction_fields.py` exercises the
  DG0 field construction; a
  `tests/fem/test_heterojunction_assembly.py` exercises the
  Poisson and DD form builders on a tiny two-region 1D mesh.
  These tests run in the gated job.

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/`
and `pytest tests/`. Do not advance with red tests. **Do not
pause between phases for confirmation.** Phase gates are
mechanical: green tests, push, advance to the next phase.
Run all seven phases end to end.

---

### Phase 0: ship this starter prompt and ADR 0016

1. Save this entire file as `docs/M17_STARTER_PROMPT.md`.
   Strip nothing; commit verbatim.
2. Author **ADR 0016: Heterojunction-aware Slotboom and
   ohmic equilibrium** at
   `docs/adr/0016-heterojunction-slotboom.md`. Status:
   Accepted, dated the day of the PR. Content:
   - Context: ADR 0004 mandates Slotboom variables for DD.
     ADR 0006 mandates an MMS verifier per physics module.
     M17 introduces position-dependent `n_i(x), N_C(x),
     N_V(x), chi(x), Eg(x)` via cellwise DG0 fields. The
     Slotboom substitution rule extends from
     `n = n_i exp((psi-phi_n)/V_t)` to
     `n = n_i(x) exp((psi-phi_n)/V_t)` without changing the
     continuity-flux shape. Cite the M16.4 generalized-
     Slotboom precedent.
   - Decision: M17 ships with (a) the position-dependent
     n_i(x) substitution above; (b) an ohmic-contact
     equilibrium psi calculation that reads local chi
     rather than reference chi; (c) MMS Variant I exercising
     smooth chi(x), Eg(x) ramps (not discontinuous); (d) an
     HEMT 2D benchmark exercising the discontinuous case
     within a 15 % quantum-confinement caveat.
   - Consequences: ADR 0004 is preserved (continuity flux
     shape unchanged). ADR 0006 is satisfied for the
     smooth-coefficient case (Variant I) and partially
     deferred for the discontinuous case (per-physics-module
     MMS does not extend to discontinuous coefficient PDEs;
     the HEMT benchmark's published-reference acceptance
     replaces the MMS gate for the discontinuous case). The
     M16.4 FD-correction prefactor `gamma_n` composes
     orthogonally with n_i(x) under FD-plus-heterojunction
     configurations.
   - Cross-reference: Stern 1972, Vurgaftman 2001 III-V
     parameters review, Sze 3rd ed § 1.5 / § 2.10. Do
     **not** cite any AI tool.
3. Append a one-line "Author M17 starter prompt and ADR 0016"
   entry to `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under
   `[Unreleased]`. (PR #85 left an Unreleased block; appending
   here is fine.)
4. Push the branch and open the PR with the template at the
   bottom of this file.

**Commit message:** `docs: ship M17 starter prompt and ADR 0016 (M17)`

**Acceptance:** the two new files exist on
`dev/m17-heterojunction`; CI green on docs-only diff; no
AI-credit lines.

---

### Phase A: schema surface for heterojunctions and AlGaAs

Pure-Python only. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.7.0 to 2.8.0. Bump
   `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`. The
   `examples` array gains `"2.8.0"` ahead of existing entries.
2. Extend the `regions[]` entry:
   - Add `material_overrides` sub-object with optional fields
     `chi_eV` (number, electron affinity in eV), `Eg_eV`
     (number, bandgap in eV), `Nc_per_cm3` (number, CB
     effective DOS in cm^-3), `Nv_per_cm3` (number, VB
     effective DOS in cm^-3). All optional; if any field is
     unset, the value comes from `regions[].material`. Any
     field set overrides the corresponding material field
     for that region only.
   - Add `heterojunction` boolean (default `false`). When
     `true` on a region join, the form builders build a
     discontinuous chi/Eg/Nc/Nv DG0 field (no smoothing
     across the interface). The default `false` smooths
     the interface (one-cell-wide transition), which is
     sufficient for the bit-identity guarantee on existing
     benchmarks. Document the default semantics in the
     description.
   - The existing `material` field stays unchanged.
3. Add AlGaAs to the supported-materials list documented in
   `docs/schema/reference.md` (cite Vurgaftman 2001 for the
   parameters; the entry shipped in Phase B). Note that
   AlGaAs is parametric in Al fraction; the materials.py
   entry is fixed at Al_0.3 and the Al fraction is encoded
   in the material name (`AlGaAs_0p3`). This is the
   defensible default; if a user wants a different Al
   fraction, they ship `material_overrides: {chi_eV: ...,
   Eg_eV: ...}` directly.
4. Default-fill: an input with no `material_overrides` and
   no `heterojunction: true` is bit-identical to v0.23.0.
5. Pure-Python tests in a new
   `tests/test_heterojunction_schema.py`:
   - Every existing benchmark JSON validates unchanged
     against v2.8.0.
   - `regions[].material_overrides: {chi_eV: 4.07, Eg_eV:
     1.42}` validates.
   - `regions[].heterojunction: true` validates.
   - A negative `Eg_eV` is rejected.
   - An unknown sub-field under `material_overrides` is
     rejected (additionalProperties: false).
6. Update `docs/schema/reference.md` versioning table with
   the v2.8.0 row; the change-log column reads "additive:
   per-region material_overrides + heterojunction switch
   (M17)".

**Commit message:** `feat(schema): regions[].material_overrides + heterojunction (schema 2.8.0); existing inputs unchanged (M17)`

**Acceptance:** all existing tests pass; new tests pass; no
FEM behavior change.

---

### Phase B: position-dependent material fields

The structural heart of M17. Promotes `n_i, N_C, N_V, chi,
Eg` from scalars to DG0 cellwise fields. Bit-identity-
preserving when the configuration is single-material.

1. Extend `semi/materials.py`:
   - Add `AlGaAs_0p3` entry: Al_0.3 Ga_0.7 As at 300 K.
     Vurgaftman 2001 derived parameters: `Eg = 1.798 eV`
     (direct gap at x = 0.3, 1.424 + 1.247*0.3), `chi =
     3.74 eV` (4.07 - 1.1*0.3, electron affinity bowing per
     Anderson rule), `Nc = 5.5e17 cm^-3` (interpolated),
     `Nv = 1.0e19 cm^-3` (interpolated), `n_i = 1.0e2
     cm^-3` at 300 K (much smaller than GaAs because of
     the wider gap), `mu_n = 1500 cm^2/V/s` (compensated
     by alloy scattering vs GaAs's 8500), `mu_p = 100
     cm^2/V/s`, `epsilon_r = 12.0`, `m_n_star = 0.063` (per
     Vurgaftman). Cite Vurgaftman 2001 in the entry comment.
   - Add `n_i_at_T(T)` instance method on `Material`
     returning `sqrt(Nc * Nv) * exp(-Eg * EV_TO_J / (2 *
     KB_J * T))`. The existing `n_i` field remains as the
     300 K reference; `n_i_at_T(T)` is the temperature-
     dependent computation used by the heterojunction
     field construction. This is a generalization the
     materials module needed eventually anyway (M16.5 and
     PR #85 examples both use `physics.temperature` and
     compute n_i at that temperature outside the materials
     module; centralizing is cleaner).
2. New module `semi/physics/heterojunction.py`:
   - Module docstring: cite ADR 0016 and the Phase B
     conventions; explain the DG0 field construction as the
     position-dependent counterpart to the scalar fields on
     `Scaling`.
   - Public API:

     ```python
     def build_dg0_material_fields(
         mesh, cell_tags, regions_cfg, sc, T, *, fdim=None,
     ):
         """
         Build per-cell DG0 fields for chi, Eg, Nc, Nv, n_i,
         epsilon_r given a mesh, cell-region tags, and the
         regions configuration. Returns a dict of
         dolfinx.fem.Function objects indexed by field name.

         For each cell, look up the region the cell belongs
         to via cell_tags, retrieve the material from
         regions_cfg[region]["material"], apply
         material_overrides if present, and assemble the
         per-cell value into the DG0 function. The n_i field
         uses the per-cell Eg, Nc, Nv at temperature T.

         The fields are scaled to dimensionless internal
         units consumed by the form builders: chi_hat = chi
         / V_t, Eg_hat = Eg / V_t, Nc_hat = Nc / sc.C0,
         Nv_hat = Nv / sc.C0, n_i_hat = n_i / sc.C0,
         epsilon_r_hat = epsilon_r (already dimensionless).
         """
     ```

     Imports `dolfinx.fem` and `ufl` inside the function
     body. Pure-Python module-scope.
   - Private helper `_resolve_region_material(region_cfg)`
     that applies `material_overrides` to the material
     fetched from `MATERIALS`. Pure-Python.
   - Private helper `_chi_ref_offset(material_field, sc)`
     returning the per-cell `chi - chi_ref` UFL expression
     used by the ohmic equilibrium psi calculation in
     Phase D.
3. Pure-Python unit tests in a new
   `tests/test_heterojunction_fields.py`:
   - Single-material 1D mesh produces a constant DG0 chi
     field (every cell value equal).
   - Two-region 1D mesh (Si + GaAs) produces a chi field
     with a step at the region boundary (one cell on the Si
     side has chi = 4.05 eV, one cell on the GaAs side has
     chi = 4.07 eV).
   - `material_overrides: {chi_eV: 3.5}` on the GaAs region
     produces a chi field with the override value on the
     GaAs side.
   - `_resolve_region_material` returns a Material instance
     with the overrides applied; the original MATERIALS dict
     is unchanged.
   - `n_i_at_T(T)` round-trip: Si at T = 300 K matches the
     existing `MATERIALS["Si"].n_i` within 1 %.
   - `n_i_at_T(T)` AlGaAs_0p3 at T = 300 K is much smaller
     than GaAs (the wider gap suppresses n_i). Document the
     ratio in the test docstring.
4. FEM-side tests in
   `tests/fem/test_heterojunction_assembly.py`:
   - A two-region 1D mesh (Si + GaAs) with no heterojunction
     flag set produces a Poisson form that assembles
     without error and is well-conditioned.
   - The same mesh with `heterojunction: true` on the
     region join produces a Poisson form that also
     assembles without error. The chi field has a step at
     the boundary (verify by sampling the DG0 field at
     points near the interface).
   - A single-material configuration produces an identical
     residual to the v0.23.0 path within 1e-12 relative
     (proves bit-identity at the residual level on a tiny
     mesh, before the full benchmark CI matrix runs in
     Phase F).

**Commit message:** `feat(physics): heterojunction.py DG0 material fields (M17)`

**Acceptance:** new tests pass; no FEM behavior change on
existing benchmarks (no benchmark uses material_overrides or
heterojunction yet; the field construction reduces to a
constant for single-material configs).

---

### Phase C: thread the position-dependent fields through Poisson and DD

Layer 4 form-builder edits. Bit-identity-preserving when the
input field is constant.

1. `semi/physics/poisson.py`:
   - `build_equilibrium_poisson_form` and the multi-region
     variant gain a `heterojunction_fields: dict | None =
     None` keyword. When None, the existing path runs (the
     Phase 0 baseline). When set, the dict carries the DG0
     fields from Phase B.
   - Replace the `ni_hat` scalar with the position-dependent
     UFL expression from `heterojunction_fields["n_i_hat"]`
     when provided. The `eps_r_ufl` UFL expression already
     handles position-dependence via the existing multi-
     region variant; cross-reference that pattern.
   - The discontinuous-coefficient PDE is fine for the
     existing P1 Lagrange psi space because the source has
     a step but psi itself is C0; the variational formulation
     handles the discontinuity in the source naturally.
2. `semi/physics/drift_diffusion.py`:
   - `build_dd_block_residual` and `_mr` gain the same
     `heterojunction_fields: dict | None = None` keyword.
     Forward to the Slotboom helpers in the inline at L150-
     L160 (and L309-L316 in `_mr`).
   - The Auger and Hurkx and BBT inlines from M16.3 / M16.6
     all consume `n_hat` and `p_hat`; they pick up the
     position-dependent values for free via the Slotboom
     helpers.
3. `semi/physics/slotboom.py`:
   - `n_from_slotboom`, `p_from_slotboom`, and the numpy
     counterparts already accept a UFL `n_i_hat` argument.
     Document in the docstring that under M17 this argument
     can be a scalar (single-material) or a UFL Function
     (heterojunction). No signature change required if the
     existing typing is permissive enough; otherwise widen
     the type hint.
4. `semi/runners/equilibrium.py`,
   `semi/runners/bias_sweep.py`,
   `semi/runners/transient.py`,
   `semi/runners/ac_sweep.py`,
   `semi/runners/mos_cv.py`,
   `semi/runners/mos_cap_ac.py`,
   `semi/runners/resistor_3d.py`:
   - At runner setup, call
     `heterojunction.build_dg0_material_fields` and pass the
     resulting dict to the form builders as
     `heterojunction_fields=...`. The build is centralized;
     each runner calls it once.
   - Bit-identity is preserved on every existing
     benchmark because the constructed DG0 field is a
     constant for single-material configs.
5. Re-run the full benchmark CI matrix at this point (or at
   least the four anchors: pn_1d_bias, diode_velsat_1d,
   diode_auger_1d, schottky_1d). Confirm bit-identity
   numerically before advancing.

**Commit message:** `feat(runners): thread heterojunction_fields through Poisson and DD form builders (M17)`

**Acceptance:** all existing benchmarks bit-identical to
v0.23.0; the FEM smoke test from Phase B step 4 passes.

---

### Phase D: ohmic-contact equilibrium psi reads local chi

Smaller surface area than Phase C. The ohmic equilibrium psi
calculation in `semi/bcs.py` reads chi from the *local*
material rather than the reference material.

1. `semi/bcs.py`:
   - Locate the `_ohmic_psi_eq` (or equivalent name)
     helper. Today it computes equilibrium psi using
     `Scaling.chi` (or implicitly the reference material's
     chi via the n_i scalar). Rewrite to accept a
     `local_chi_eV` argument and compute the
     `chi_local - chi_ref` shift correctly. The single-
     material path collapses to the existing formula when
     `local_chi_eV == Scaling.chi_ref_eV`.
   - The build_psi_dirichlet_bcs and build_dd_dirichlet_bcs
     paths look up the contact's region from the facet
     tag, fetch the local material from
     `regions_cfg[region]`, and pass the local chi to
     `_ohmic_psi_eq`.
   - The Schottky equilibrium psi from M16.5
     (`_schottky_psi_eq`) similarly accepts a local-chi
     argument; cross-reference that pattern.
   - Gate-contact psi (work-function compensation) reads
     the local chi from the semiconductor side of the gate
     stack; today this works via the SiO2 / Si interface
     handling in the multi-region variant. M17 doesn't
     change the gate-stack behavior beyond ensuring chi is
     read from the local Si region (it already is in
     practice; document that fact in the docstring).
2. Pure-Python tests in `tests/test_bcs.py`:
   - `_ohmic_psi_eq` with `local_chi_eV == reference_chi_eV`
     returns the same value as the v0.23.0 path within
     1e-12 relative (regression).
   - `_ohmic_psi_eq` with a 0.5 eV chi shift returns a value
     shifted by 0.5 V from the reference (the chi shift
     translates directly into a psi shift in the Anderson-
     rule electrostatic alignment).
3. FEM-side test in
   `tests/fem/test_heterojunction_ohmic.py`:
   - Two-region 1D mesh (Si + AlGaAs_0p3, ohmic contacts on
     both ends) at equilibrium. Verify that the
     equilibrium psi profile has the expected band-edge
     alignment across the heterojunction (the discontinuity
     in chi shows up in psi via the Anderson rule). The
     test asserts the psi step magnitude matches the chi
     step within 1 %.

**Commit message:** `feat(bcs): ohmic equilibrium psi reads local chi (M17)`

**Acceptance:** all existing benchmarks remain bit-identical
to v0.23.0; the new heterojunction test passes.

---

### Phase E: MMS Variant I

ADR 0006 mandates an MMS verifier for the smooth-coefficient
case. M17's discontinuous case is exercised by the HEMT
benchmark in Phase F, not MMS.

1. Extend `semi/verification/mms_dd.py`:
   - `VARIANTS = ("A", "B", "C", "D", "E", "F", "G", "H",
     "I")`.
   - Variant I uses the same Variant C `psi_e, phi_n_e,
     phi_p_e` manufactured fields plus smooth chi(x) and
     Eg(x) ramps. For 1D: `chi(x) = chi_0 + delta_chi *
     x/L` and `Eg(x) = Eg_0 + delta_Eg * x/L` with
     `delta_chi`, `delta_Eg` engineered to give an O(0.1)
     shift in `n_i(x)` across the domain (so the
     position-dependence is materially exercised but not
     so extreme that the Slotboom substitution is
     numerically pathological).
   - `_build_weak_sources` substitutes the position-
     dependent `n_i(x)` into the manufactured weak source.
     The manufactured solution `psi_e, phi_n_e, phi_p_e`
     does not change; only the substitution rule for n
     does.
   - `run_one_level` constructs a one-region cfg with
     `material_overrides` set to the manufactured chi/Eg
     ramps via a small per-cell-tag mechanism (or by
     manufacturing on a fine sub-region partition; pick
     whichever is cleaner with the existing harness).
   - CLI study runs Variant I on the same Ns_1d / Ns_2d
     that Variant H (M16.6) runs on, or one level coarser
     if the finest level hits the residual floor.
   - Acceptance: finest-pair L2 rate >= 1.99 and H1 rate
     >= 0.99 on each block. The position-dependence is
     smooth so the rate gate is the textbook P1 gate; no
     decoupling caveat applies (M16.6's psi-only rate gate
     does not extend here because chi(x) couples directly
     into both n and p substitution).
2. Extend `semi/verification/mms_poisson.py` with an
   equilibrium-Poisson MMS variant exercising the
   position-dependent chi/Eg/Nc/Nv. Same smooth-ramp
   approach. Acceptance: L2 rate >= 1.99 / H1 rate >= 0.99.
   This is the equilibrium counterpart to MMS-DD Variant I.
3. New file `tests/fem/test_mms_heterojunction.py`:
   - 1D-rate test, 2D-rate test for MMS-DD Variant I.
   - 1D-rate test, 2D-rate test for the MMS-Poisson
     equilibrium variant.
   Mirror `tests/fem/test_mms_tunneling.py` (M16.6) for
   structure.
4. Update `docs/PHYSICS.md` § Verification & Validation:
   - Add a "Variant I (heterojunction smooth ramp, M17)"
     bullet under Variant H. Append the rates to the
     finest-pair-rates table once measured.
5. Update `docs/mms_dd_derivation.md` with the manufactured-
   source derivation under position-dependent n_i(x). The
   derivation factors cleanly because the Slotboom
   substitution preserves its shape; cite the M16.4
   precedent.

**Commit message:** `feat(verification): MMS-DD heterojunction Variant I + MMS-Poisson smooth-ramp (M17)`

**Acceptance test 3 (MMS gate)**: `python
scripts/run_verification.py mms_dd` includes Variant I and
reports L2 rate >= 1.99 finest-pair on each block; the
Poisson MMS rate gate is also met.

---

### Phase F: the HEMT 2D benchmark

The analytical anchor: an AlGaAs/GaAs HEMT 2D structure with
the 2DEG sheet density at the heterojunction within 15 % of a
published Poisson-Schrodinger reference. The 15 % tolerance
is generous because the FEM solver is purely classical (no
quantum confinement) and the Poisson-Schrodinger reference
includes confinement effects; the agreement on integrated
sheet density at moderate gate biases is the gate.

1. New directory `benchmarks/hemt_2d/`:
   - `hemt.json`: 2D AlGaAs/GaAs HEMT. Geometry:
     - Top gate (Schottky contact, work function 4.7 eV;
       reuse M16.5 Schottky infrastructure).
     - 30 nm AlGaAs_0p3 barrier.
     - At y = 30 nm, Si delta-doping plane (modeled as a
       1 nm-thick uniformly-doped slab at N_D = 5e18 cm^-3
       to give a ~5e12 cm^-2 areal donor density; the
       delta-function approximation is built into the
       benchmark mesh).
     - 100 nm GaAs channel (lightly doped, p-type to
       suppress parasitic conduction; N_A = 1e15 cm^-3).
     - 200 nm GaAs buffer below the channel.
     - Source and drain ohmic contacts at the GaAs channel
       edges (left and right).
     - Total 2D extent: 1 um wide x 330 nm tall.
     - Mesh: ~150 x ~80 cells with refinement near the
       heterojunction (5 nm cell size on the GaAs side, 1
       nm cell size in the AlGaAs near the heterojunction).
   - `physics.statistics: "fermi_dirac"` (heavy n+ doping;
     FD is required for the channel charge accuracy).
   - `physics.mobility: {model: "constant"}` (the 2DEG
     sheet density depends on the electrostatic charge
     balance, not the channel mobility; minimize variables
     for the V&V gate).
   - `regions[].heterojunction: true` on the AlGaAs / GaAs
     join.
   - V_GS sweep [0.0, 1.0] V step 0.1 V (11 points).
   - `README.md`: device description, the analytical
     reasoning for the 15 % tolerance (cite Stern 1972 or
     Pozela & Reklaitis 1980 for the published Poisson-
     Schrodinger reference; document the no-quantum-
     confinement caveat). Do **not** cite any AI tool.
2. Extend `semi/diode_analytical.py` (or a new
   `semi/hemt_analytical.py`):
   - Add `hemt_2deg_classical(V_GS, barrier_thickness,
     N_doping_delta, chi_AlGaAs, Eg_AlGaAs, chi_GaAs,
     Eg_GaAs, T)` returning the classical-electrostatic
     2DEG sheet density at a given gate bias. The
     Pozela-Reklaitis derivation gives a closed-form
     relation `n_s(V_GS) = (eps_AlGaAs / (q
     barrier_thickness)) * (V_GS - V_T_HEMT)` in the
     linear-regime above threshold. The classical reference
     differs from the published Poisson-Schrodinger
     reference by ~5-15 % depending on bias; the verifier
     uses the classical reference (which the FEM solver
     should match within 5 %) and the 15 % tolerance
     accommodates the gap to the quantum reference.
   - Pure-Python; no dolfinx imports.
3. New verifier `verify_hemt_2d` in
   `scripts/run_benchmark.py`:
   - Run the V_GS sweep.
   - At each V_GS, integrate `n(y) - N_D(y)` over y from
     y = 30 nm (the heterojunction) to y = 60 nm (30 nm
     into the GaAs channel; this captures the full 2DEG
     extent at the relevant biases). Convert to areal
     density n_s in cm^-2.
   - Compute the classical-electrostatic reference via
     `hemt_2deg_classical`.
   - Assert `|n_s_FEM(V_GS) - n_s_classical(V_GS)| /
     n_s_classical(V_GS) < 0.15` for V_GS in [0.4, 1.0] V
     (the linear regime above threshold). Below 0.4 V the
     2DEG is sub-threshold and the relative error blows
     up; the verifier excludes that range.
   - Mirror the structure of `verify_zener_1d` and
     `verify_schottky_1d`.
4. Wire the new benchmark into CI: add a step in
   `.github/workflows/ci.yml` matrix near
   `diode_fermi_dirac_1d` and `zener_1d`:
   ```yaml
   - name: hemt_2d
     args: hemt_2d
   ```
   **Do not** mark `allow-failure: true`. The benchmark must
   pass cleanly. If it doesn't (worst case: the 2DEG
   extraction is off because of mesh resolution near the
   heterojunction), refine the mesh or adjust the
   classical reference; do not relax the 15 % gate.
5. Plotter `plot_hemt_2d`: n_s vs V_GS on a linear axis,
   with the classical reference overlaid. Companion subplot
   showing the y-profile of n at V_GS = 1.0 V (the 2DEG
   peak near the heterojunction visible).
6. Coverage hook: confirm the FEM smoke tests from Phase B
   and Phase D cover the heterojunction code paths; the
   HEMT benchmark itself runs in a separate CI matrix entry
   whose coverage is not merged into the gate.

**Commit message:** `feat(benchmark): hemt_2d AlGaAs/GaAs heterojunction 2DEG verifier (M17)`

**Acceptance test 1**: `python scripts/run_benchmark.py hemt_2d`
exits 0 with the 15 % 2DEG match passing on every V_GS in
[0.4, 1.0] V.

---

### Phase G: closeout

1. `PLAN.md`:
   - Move M17 from "Next task" to "Completed work log" with
     PR number, deliverables, schema minor bump (2.7.0 to
     2.8.0), and acceptance-test results (cite observed
     worst-case 2DEG match from Phase F and observed
     Variant I L2 / H1 rates from Phase E).
   - Update the "Honest current state" gap line: M17 is
     shipped; the remaining gap-list items are M19 (3D MOSFET
     capstone), M19.1 (MPI parallel), M20 (HTTP server
     hardening), and the post-M19 items (impact ionization,
     etc.).
   - Set "Next task" to **M19 (3D MOSFET capstone)** with a
     pointer to the (yet-to-be-authored) M19 starter prompt.
     M19 depends on M16.1 (mobility) and is unblocked.
   - Refresh "Current state" with the new package version
     (bump 0.23.0 to 0.24.0).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M17 Done in § 4 with a one-line summary and a
     CHANGELOG anchor.
   - Update the gap list at L86-87 (currently lists M17
     among the "next-tier gaps"); rewrite to:
     ```
     - **M16 umbrella complete and M17 heterojunctions
       shipped.** Remaining roadmap items are M19 (3D MOSFET
       capstone), M19.1 (MPI parallel), and M20 (HTTP server
       hardening). Post-M19 items (impact ionization,
       advanced trap models) are tracked separately.
     ```
   - Append a § 9 changelog entry under `[0.24.0]`. Move
     the existing `[Unreleased]` block into `### Released`.
3. `docs/PHYSICS_INTRO.md`:
   - § 6: append a heterojunctions bullet:
     ```
     - **Heterojunctions.** Position-dependent electron
       affinity and bandgap with the AlGaAs/GaAs HEMT 2D
       benchmark as the reference (M17,
       `benchmarks/hemt_2d/`). Anderson-rule band alignment;
       classical electrostatic 2DEG within 15 % of a
       published Poisson-Schrodinger reference (the gap
       reflects the lack of quantum confinement in the
       classical solver).
     ```
     Insert after the BBT/TAT bullet from M16.6 and before
     the M16.7 transient bullet.
   - § 7: drop the "no heterojunctions" bullet entirely.
     The remaining § 7 bullets are impact ionization, 3D
     MOSFET, and quantum confinement (the new addition; the
     HEMT benchmark's 15 % tolerance reflects this gap).
4. `docs/ROADMAP.md`:
   - Update the M17 row in the capability matrix from
     Planned to shipped.
   - Promote M19 to the top of the Planned column.
5. `CHANGELOG.md`:
   - New `[0.24.0]` entry with the M17 line items.
   - Update the schema-banner comment at the top to mention
     v2.8.0.
6. `pyproject.toml` and `semi/__init__.py`: bump 0.23.0 to
   0.24.0.
7. Push every commit to `origin/dev/m17-heterojunction`.
   Verify no AI-credit trailers anywhere
   (`git log --format=%B dev/m17-heterojunction | grep -i
   'co-authored-by: claude'` returns empty).

**Commit message:** `docs: close out M17 (PLAN, IMPROVEMENT_GUIDE, PHYSICS_INTRO, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact.
- [ ] `git log --format=%B dev/m17-heterojunction | grep -i
      'co-authored-by: claude'` returns empty.
- [ ] Pure-Python core remains dolfinx-free
      (`tests/test_lazy_imports.py` clean).
- [ ] Single-material configs (no `material_overrides`, no
      `heterojunction: true`) are bit-identical to v0.23.0
      on every existing benchmark and every PR #85 example.
      Numerical anchors as listed in Conventions.
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004; preserved via the generalized-
      Slotboom path documented in ADR 0016).
- [ ] `make_scaling_from_config` still on every solve path;
      heterojunction.build_dg0_material_fields composes on
      top of it.
- [ ] No PETSc / UFL types leak into `kronos_server` public
      API.
- [ ] Schema bumped per minor (2.7.0 to 2.8.0); v2.0.0
      through v2.7.0 inputs still validate.
- [ ] MMS rate gate L2 >= 1.99 / H1 >= 0.99 active for
      Variants A through I (Variant I is smooth-coefficient
      only; the discontinuous heterojunction case is
      exercised by hemt_2d, not MMS, per ADR 0016).
- [ ] Coverage gate holds at 95 in the gated `docker-fem-tests`
      job.

## Anti-goals

- Do not start M19 (3D MOSFET) in this PR. M19 is its own
  PR; the M19 starter prompt is the follow-up deliverable
  for this PR's reviewer.
- Do not add quantum confinement corrections (Schrodinger
  coupling, Bohm potential, density-gradient). The 15 % HEMT
  tolerance reflects the lack of quantum confinement; a
  proper quantum-corrected HEMT is its own M-numbered
  milestone.
- Do not introduce SUPG, streamline diffusion, or a primary-
  density form. ADR 0004 is locked.
- Do not add MMS variants for the discontinuous-coefficient
  case. ADR 0016 documents that the smooth-coefficient
  Variant I plus the HEMT benchmark are the V&V gate;
  manufactured solutions for discontinuous-coefficient PDEs
  are out of scope.
- Do not add interface trap densities (D_it) at the
  heterojunction. Interface traps are a follow-up scope; the
  HEMT benchmark assumes ideal (trap-free) interfaces.
- Do not retire the `mosfet_2d` `allow-failure: "true"` flag.
- Do not bundle with M19, M19.1, or M20.
- Do not add `Co-Authored-By` trailers, generated-by footers,
  or any other AI-credit marker.

## Stop conditions

You are done when:

1. The M17 PR is opened on branch `dev/m17-heterojunction`.
2. All three acceptance tests in `docs/IMPROVEMENT_GUIDE.md`
   § M17 pass in CI:
   - A1: `benchmarks/hemt_2d` 2DEG sheet density within
     15 % of the published reference at every V_GS in
     [0.4, 1.0] V.
   - A2: every existing benchmark and every PR #85 example
     bit-identical to v0.23.0 (single-material default).
   - A3: MMS Variant I L2 >= 1.99 / H1 >= 0.99 finest-pair
     on each block.
3. ADR 0016 (heterojunction-aware Slotboom and ohmic
   equilibrium) is on `main`.
4. PLAN.md, IMPROVEMENT_GUIDE.md, PHYSICS_INTRO.md,
   ROADMAP.md, CHANGELOG.md reflect the closeout.
5. Package version bumped 0.23.0 to 0.24.0 in
   `pyproject.toml` and `semi/__init__.py`.
6. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI
   assistant. Squash-merge body is also clean.
7. Coverage gate holds at 95 in the gated `docker-fem-tests`
   job, not via a follow-up commit.
8. PR reviewed, CI green (modulo the documented `allow-
   failure` on `mosfet_2d`), merged.

## PR description template

```
## Summary

M17: Heterojunction / position-dependent band structure. The
largest single physics PR since M14.3.

Promotes n_i, N_C, N_V, chi, Eg from scalar fields on the
reference material to position-dependent cellwise DG0 fields
on the mesh. Extends the materials model with AlGaAs_0p3
(Vurgaftman 2001 derived parameters). Schema additive minor
bump v2.7.0 to v2.8.0 (regions[].material_overrides + a
heterojunction switch on the region join).

ADR 0016 documents the heterojunction-aware Slotboom path:
the substitution rule extends from `n = n_i exp((psi-phi_n)/
V_t)` to `n = n_i(x) exp((psi-phi_n)/V_t)` without changing
the continuity-flux shape. ADR 0004 preserved. Composes
orthogonally with M16.4 Fermi-Dirac (gamma_n combines
multiplicatively with n_i(x) when both are enabled).

Single-material configs (no material_overrides, no
heterojunction: true) are bit-identical to v0.23.0 on every
existing benchmark and every PR #85 example.

New benchmark: hemt_2d AlGaAs/GaAs HEMT 2D, 2DEG sheet
density within 15 % of a published Poisson-Schrodinger
reference at every V_GS in [0.4, 1.0] V. The 15 % tolerance
reflects the lack of quantum confinement in the classical
solver; per ADR 0016 this is the V&V gate for the
discontinuous case.

New MMS variants: Variant I (smooth chi/Eg ramp, finest-pair
L2 >= 1.99 / H1 >= 0.99 on each block) plus an equilibrium-
Poisson MMS smooth-ramp variant.

## Acceptance tests

(Numbered in docs/IMPROVEMENT_GUIDE.md § M17.)

- [ ] A1 (2DEG match): hemt_2d 2DEG sheet density within
      15 % of published reference at every V_GS in
      [0.4, 1.0] V (observed worst: <fill in>)
- [ ] A2 (byte-identity): every existing benchmark and PR #85
      example bit-identical to v0.23.0
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
      diode_auger_1d >20 % divergence at 0.9 V;
      diode_fermi_dirac_1d 7.37 % FD-vs-Boltzmann V_bi at
      N_D=1e20 cm^-3; schottky_1d worst-case <10 %;
      zener_1d worst-case <20 % Kane match)
- [ ] A3 (MMS gate): scripts/run_verification.py mms_dd
      Variant I L2 >= 1.99 / H1 >= 0.99 finest-pair on
      each block (observed: <fill in>)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95 (in the gated
      docker-fem-tests job)
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark hemt_2d
- [ ] docker compose run --rm benchmark pn_1d_bias
      (single-material, bit-identical)
- [ ] docker compose run --rm benchmark diode_velsat_1d
      (single-material, bit-identical)
- [ ] docker compose run --rm benchmark diode_auger_1d
      (single-material, bit-identical)
- [ ] docker compose run --rm benchmark diode_fermi_dirac_1d
      (single-material, bit-identical)
- [ ] docker compose run --rm benchmark schottky_1d
      (single-material, bit-identical)
- [ ] docker compose run --rm benchmark zener_1d
      (single-material, bit-identical)
- [ ] PR #85 example smoke verifiers green

## Notes

- The mosfet_2d CI matrix entry remains `allow-failure:
  "true"`, unchanged in scope from M16.7.
- AlGaAs_0p3 added to materials.py with Vurgaftman 2001
  derived parameters at 300 K; other Al fractions can be
  modeled via `material_overrides: {chi_eV, Eg_eV, ...}`.
- The 15 % HEMT tolerance is the gap between the classical
  electrostatic solver and a fully-quantum Poisson-
  Schrodinger reference; quantum confinement corrections
  are a future M-numbered milestone, not this PR.
- ADR 0016 is the V&V departure documentation; ADR 0006 is
  satisfied for the smooth-coefficient case (Variant I) and
  partially deferred for the discontinuous case (replaced
  by the hemt_2d benchmark gate).
- Coverage gate held at 95 in the gated docker-fem-tests
  job via tests/fem/test_heterojunction_*.py, avoiding the
  follow-up-commit pattern from M16.5 / M16.6.
```

## Hand-off

When M17 lands and is reviewed, the next pickup is **M19 (3D
MOSFET capstone)**. M19 ships:

- A 3D MOSFET on a gmsh-sourced unstructured mesh.
- Exercises the M15 GPU linear-solver path under non-trivial
  geometry.
- Exercises the M16.1 Caughey-Thomas mobility under 3D field
  conditions.
- Expected to need MPI parallel orchestration (M19.1) for
  runtime to be tolerable on a single-machine CI runner.

M19 depends on M16.1 (already shipped) and M15 (already
shipped); it is unblocked. The M17 heterojunction infrastructure
shipped in this PR is *not* a dependency of M19 (3D MOSFET is
single-material Si by default), but the position-dependent
material field machinery is reusable for any future 3D
heterojunction device benchmark (e.g., 3D HEMT, 3D HBT) that
M19 might inspire.

After M19, M19.1 (MPI parallel) and M20 (HTTP server
hardening) close the post-M16 roadmap. Post-M20 candidates
include impact ionization (post-M19), quantum confinement
corrections (would close the HEMT 15 % gap), and advanced
trap models (interface traps, deep-level transient
spectroscopy support).
