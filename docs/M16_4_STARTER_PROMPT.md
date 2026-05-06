## Context

You are working in `kronos-semi` at `v0.19.0` (post-M16.3, package
version `0.19.0`, schema `2.3.0`). The repo is at
`https://github.com/rwalkerlewis/kronos-semi`; `main` is at
`22dc145` (`M16.3: Auger recombination (#79)`).

Your assignment is **M16.4: Fermi-Dirac statistics (gated)**, the
fourth physics-completeness slice of the M16 umbrella. M16.1
(Caughey-Thomas mobility), M16.2 (Lombardi surface mobility), and
M16.3 (Auger recombination) have merged. M16.4 is structurally the
heaviest of the M16.x slices to date because it touches the carrier
statistics that the Slotboom primary unknowns are derived against,
and ADR 0004 locks Slotboom in. **Read the "Slotboom interaction"
section below before any code change.**

`PLAN.md` § "Next task" names M16.4 explicitly:

> **M16.4: Fermi-Dirac statistics (gated)** on a fresh branch
> `dev/m16.4-fermi-dirac`. To be picked up via
> `docs/M16_4_STARTER_PROMPT.md` (yet to be authored, in the same
> shape as `docs/M16_3_STARTER_PROMPT.md`); acceptance tests in
> `docs/IMPROVEMENT_GUIDE.md` § M16.4. Boltzmann breaks above
> ~1e19 cm^-3, which is the source/drain extension regime of every
> modern MOSFET.

This prompt **is** the file PLAN.md is asking you to author. Save
it as `docs/M16_4_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below), preserving the convention M16.3 inherited from
M16.2 and M16.1.

This prompt does not restate the M16.4 deliverable, the acceptance
tests, or the rationale; those live in `docs/IMPROVEMENT_GUIDE.md`
§ M16.4. Read the guide first; this prompt only tells you the order
in which to execute and which invariants must remain load-bearing.

## Branch and PR rules

- Work on a fresh branch off `main`: `git checkout -b dev/m16.4-fermi-dirac`.
  Do **not** rebase or commit onto `main` directly.
- One milestone, one PR. Do not bundle M16.4 with M16.5 (Schottky)
  or any other physics work. M16.4 is independent of M16.2
  Lombardi and M16.3 Auger (composes orthogonally with both); do
  not retouch their surface area.
- Push every phase commit to `origin/dev/m16.4-fermi-dirac`
  immediately after it lands locally. Open the PR after Phase 0 so
  reviewers can watch phases land.
- Title the PR `M16.4: Fermi-Dirac statistics (gated)`. Use the PR
  description template at the bottom of this prompt.
- **Do not credit Claude, Claude Code, or any AI assistant in
  commit messages, PR descriptions, code comments, or any other
  shipped artifact.** Standard "Co-authored-by" trailers and
  generated-by lines are off the table for this PR.

## Required reading (do not skip; ~60 minutes; this milestone is
heavier than the M16.1 / M16.2 / M16.3 reading load)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.19.0 (post-M16.3) and
   that no other PR is in flight on `dev/m16.4-*`. The "Next task"
   section names M16.4; if it does not, stop and align with the
   maintainer before continuing. Read the M16.3 entry in
   "Completed work log" (top of the log) end to end.
2. `docs/IMPROVEMENT_GUIDE.md` § M16 (umbrella context), § M16.1 /
   § M16.2 / § M16.3 (the just-shipped milestones, for context on
   the schema versioning and benchmark-verifier patterns), § M16.4
   (Why / Deliverable / Acceptance / Dependencies), and § 1
   (Honest current state).
3. `docs/ROADMAP.md` § Capability matrix and § Honest gap. M16.4
   is listed as Planned; the row moves to shipped at the end of
   this PR.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M16.4 touches
   the schema (Layer 2), pure-Python statistics helpers (Layer 3),
   the Slotboom UFL builders and Poisson form (Layer 4), and the
   benchmark / verifier (Layer 5). The pure-Python statistics
   module must remain dolfinx-free.
5. `docs/PHYSICS.md`:
   - § 1.2 (carrier statistics under Boltzmann, the derivation
     M16.4 generalizes).
   - § 1.3 (continuity in Slotboom current form; this is what
     ADR 0004 locks in and what the FD path must remain
     compatible with).
   - § 2.1 / § 2.2 / § 2.5 (scaling). Under FD the Poisson source
     `n_i (exp(-psi) - exp(psi))` becomes `N_C F_{1/2}((-psi -
     E_C/V_t)/...) - N_V F_{1/2}((psi - E_V/V_t)/...)`; the
     scaling-derivation lines in § 2.2 are touched.
   - § Verification & Validation (Variants A through F MMS-DD;
     you are adding Variant G).
6. `docs/adr/`:
   - **0004 (Slotboom variables for DD).** This is the load-
     bearing constraint. The Slotboom derivation in ADR 0004
     and PHYSICS.md § 1.3 is explicitly Boltzmann-based:
     `n = n_i exp((psi - phi_n) / V_t)` is substituted into
     `J_n = q mu_n n E + q D_n grad(n)` and the Einstein relation
     `D_n = mu_n V_t` collapses the result to `J_n = -q mu_n n
     grad(phi_n)`. Under FD the Einstein relation acquires a
     density-dependent correction factor (the "FD Einstein
     factor", `g(eta) = F_{1/2}(eta) / F_{-1/2}(eta)` where eta
     is the reduced Fermi level), and `n = n_i exp((psi -
     phi_n) / V_t)` no longer holds. The supported textbook fix
     is the **generalized Slotboom variables** (sometimes called
     Bernoulli-Slotboom or "modified Slotboom"): redefine
     `n_hat = gamma_n(...) * exp(psi - phi_n)` where `gamma_n`
     is the Fermi-Dirac correction. The continuity rows keep
     their `J = -q mu n grad(phi)` shape; the only change is
     that `n` is now built via the generalized expression.
     **This means the change to the DD form builder is a
     replacement of the `n_from_slotboom` / `p_from_slotboom`
     helpers, not a rewrite of the residual.** ADR 0004 stays
     intact. Open a new ADR (number 0015) if the implementation
     needs a fundamentally different DD form; do **not** deviate
     from ADR 0004 silently.
   - 0001 (JSON contract), 0002 (nondimensionalization), 0006
     (V&V strategy; M16.4 needs an MMS variant), 0007 (BC
     interface; do not touch).
7. `docs/mms_dd_derivation.md`. Variants D / E / F added field-,
   surface-, and density-dependent terms; Variant G adds an FD
   correction factor that is a smooth nonlinear function of
   `(psi - phi_n)` and a hard prerequisite to a meaningful
   manufactured solution at high doping.
8. `semi/physics/slotboom.py` end to end. The current module is
   ~50 lines: `n_from_slotboom`, `p_from_slotboom`, NumPy
   counterparts. Under FD these two helpers grow a `statistics`
   keyword argument or are renamed and a dispatcher is added (see
   Phase B below for the chosen layout).
9. `semi/physics/poisson.py` L73 and L112 (the
   `rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi)) + N_hat_fn`
   inlines in `build_equilibrium_poisson_form` and the
   multi-region variant). Under FD the
   `(exp(-psi) - exp(psi))` factor becomes
   `(F_{1/2}-based term for holes) - (F_{1/2}-based term for
   electrons)`; the substitution is via the same statistics
   dispatcher.
10. `semi/physics/drift_diffusion.py` L150-L151 (the
    `n_hat = n_from_slotboom(psi, phi_n, ni_hat)` and
    `p_hat = p_from_slotboom(psi, phi_p, ni_hat)` lines in
    `build_dd_block_residual`) and L309-L316 (the same in
    `build_dd_block_residual_mr`).
11. `semi/verification/mms_dd.py` L337-L338 (the manufactured
    `n_e = ni_hat * ufl.exp(psi_e)` / `p_e = ni_hat *
    ufl.exp(-psi_e)` for thermal-equilibrium probing). Under
    Variant G these become FD-aware; Variants A-F continue to
    use the existing Boltzmann manufactured forms unchanged
    (they are gated by the `statistics: "boltzmann"` config).
12. `semi/verification/mms_poisson.py` L8-L16, L175-L181,
    L199 (the manufactured Poisson source). Same FD substitution
    pattern as the production form.
13. `semi/fem/scharfetter_gummel.py` L234-L235 (the
    Scharfetter-Gummel Bernoulli function `B(x)`). The SG
    primitives are dead-on-active-path per ADR 0012 (M14.3
    deletion); the file is kept for historical reasons. M16.4
    does **not** touch this file.

## Slotboom interaction (read this twice before Phase B)

ADR 0004 is the load-bearing constraint for M16.4. The natural
question "does FD invalidate Slotboom?" has a published answer:
no, but the substitution rule changes. Under Boltzmann:

```
n = n_i * exp((psi - phi_n) / V_t)         (textbook Slotboom)
```

Under Fermi-Dirac, with reduced Fermi level
`eta_n = (E_F_n - E_C) / kT`:

```
n = N_C * F_{1/2}(eta_n)
```

where `F_{1/2}` is the order-1/2 Fermi-Dirac integral. The
generalized-Slotboom form preserved by ADR 0004 reads:

```
n = n_i * gamma_n(psi, phi_n) * exp((psi - phi_n) / V_t)
```

where `gamma_n` is the FD-correction prefactor. In closed form
with the Blakemore approximation (`F_{1/2}(eta) ~ 1 / (exp(-eta) +
0.27)`), `gamma_n` is a smooth bounded function of
`(psi - phi_n)` that approaches 1 in the non-degenerate limit.
The continuity row `div(mu n grad(phi_n)) = R` is unchanged in
shape; only the substitution rule for `n` in terms of the
primary unknowns changes.

**Implementation consequence.** The Slotboom helpers
(`n_from_slotboom`, `p_from_slotboom`) become dispatch points on
`statistics`. The DD residual builder, the Poisson residual
builder, and the MMS harnesses call these helpers; the only edits
to the residual builders are adding a `statistics_cfg` keyword and
forwarding it to the helper. The form structure does not change.
**Do not rewrite the DD residual.**

The Einstein-relation correction (`D = mu V_t g(eta)` with
`g(eta) = F_{1/2}(eta) / F_{-1/2}(eta)`) does **not** appear in
the Slotboom residual because the Slotboom-form current is
`J = -q mu n grad(phi_n)` regardless of the statistics: the
Einstein factor cancels against the FD-correction prefactor in
the substitution. This is exactly why generalized Slotboom is
the standard production-FD path; cite Schenk 1998 and the
Sentaurus device manual in the Phase B docstring. Cross-check the
algebra against PHYSICS.md § 1.3 before committing Phase B; if the
cancellation does not work out cleanly under your derivation, stop
and open a new ADR (0015) before continuing.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** The new schema entry
  `physics.statistics: "boltzmann" | "fermi_dirac"` (currently a
  one-element enum at `"boltzmann"` only) must be expressible in
  `schemas/input.v2.json` (current strict v2.3.0 from M16.3),
  validated by `semi/schema.py`, and exercised by at least one
  benchmark JSON.
- **Schema versioning is binding.** This is an additive minor
  bump v2.3.0 to v2.4.0 (the field name `statistics` already
  exists with a one-element enum; widening the enum is additive,
  not breaking). v2.0.0 / v2.1.0 / v2.2.0 / v2.3.0 inputs must
  continue to validate and produce bit-identical results to
  v0.19.0. Update `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Five layers, enforced.**
  - Pure-Python `semi/physics/statistics.py` (Layer 3) ships the
    Blakemore approximation and the closed-form `gamma_n`,
    `gamma_p` helpers as NumPy functions. **No `import dolfinx`,
    `import ufl`, `import petsc4py` at module scope.**
  - The UFL versions of the same helpers live in
    `semi/physics/slotboom.py` (Layer 4); imports inside function
    bodies. The pattern matches the existing
    `semi/physics/recombination.py` (Layer 4 with both NumPy and
    UFL helpers in one module, UFL imports deferred to function
    bodies).
  - The verification reference (full `F_{1/2}` via
    `scipy.special.fdk` or equivalent; `mpmath.polylog(1.5, ...)`
    is an acceptable backup if `scipy.special.fdk` is not on
    the docker-fem image) lives next to the Blakemore helper in
    `semi/physics/statistics.py` and is **only called from the
    benchmark verifier and the unit tests**, not from the
    production residual.
- **Slotboom variables for DD.** ADR 0004 is locked. The
  generalized-Slotboom path (see "Slotboom interaction" above) is
  the supported FD route. Do not introduce SUPG, streamline
  diffusion, or a primary-density form. If the cancellation
  argument breaks down under your derivation, stop and open
  ADR 0015 before continuing.
- **Every new physics module needs an MMS verifier.** ADR 0006.
  No exceptions. Acceptance L2 rate >= 1.99 and H1 rate >= 0.99
  at the finest pair, gated in `scripts/run_verification.py
  mms_dd`.
- **Every new analytical model needs a benchmark.**
  `diode_fermi_dirac_1d` is the analytical anchor; it must show
  >15 % FD-vs-Boltzmann divergence on V_bi at N_D = 1e20 cm^-3
  in the n+ region and the FD result must match a stand-alone
  scipy-driven reference within 1e-3 (relative).
- **Per-runner threading.** The `bias_sweep`, `transient`, and
  `ac_sweep` runners thread `statistics_cfg` through to the form
  builders the same way M16.3 threaded `recomb_cfg`. The
  `equilibrium`, `mos_cv`, `mos_cap_ac`, and `resistor_3d`
  runners need the same threading because Poisson is now FD-
  aware too. **All seven runners get the threading in this PR.**
  This is a wider surface than M16.3 because both Poisson and DD
  consume the FD substitution.
- **One milestone, one PR.** Do not bundle M16.4 with M16.5
  (Schottky) or any other physics work, and do not retire the
  `mosfet_2d` `allow-failure: "true"` flag (that is a separate
  follow-up, unchanged in scope from M16.2 / M16.3).
- **No em dashes in prose or code comments.** Use commas,
  periods, parentheses, or colons. Re-grep new diffs after each
  phase.
- **No AI-assistant credits in shipped artifacts.** Commit
  messages, PR descriptions, code, comments, docs, and benchmark
  README files do not mention Claude, Claude Code, Anthropic, or
  any AI tool. Standard "Co-authored-by" trailers for AI
  assistants are off the table.
- **Physics-style variable names are allowed.** `eta_n`,
  `gamma_n`, `F_half`, `N_C`, `N_V`, `psi_minus_phi_n`, etc.
- **Coverage gate is 95 on `semi/`.** The new
  `semi/physics/statistics.py` module must carry pure-Python
  unit tests for the closed-form expressions (Blakemore vs
  scipy reference, gamma_n / gamma_p limits, Einstein-factor
  cancellation algebra); FEM wiring is covered by the MMS
  test and the benchmark.
- **Boltzmann-default byte-identity.** Every existing benchmark
  with `physics.statistics: "boltzmann"` (the default) must
  produce bit-identical results to v0.19.0. Numerical anchors:
  `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2` (M16.1; reuse),
  `diode_velsat_1d` 56.27 % @ V_F = 0.9 V, 0.19 % @ V_F = 0.3 V
  (M16.1 / M16.2 / M16.3; reuse), and `diode_auger_1d` >20 %
  divergence at V_F = 0.9 V (M16.3; reuse).

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/`
and `pytest tests/`. Do not advance with red tests.

---

### Phase 0: ship this starter prompt

Pure docs, no code. Establishes the convention M16.3 inherited
from M16.2 and M16.1, and gives reviewers the contract for the
rest of the PR.

1. Save this entire file (the prose you are reading now,
   including the PR template at the bottom) as
   `docs/M16_4_STARTER_PROMPT.md`. Strip nothing; the file is
   committed verbatim. **Do not** add a footer crediting Claude
   or any AI tool.
2. Append a one-line "Author M16.4 starter prompt" entry to
   `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under a new
   `[Unreleased]` heading if one is not already present.
3. Push the branch and open the PR with the template at the
   bottom of this file. The PR opens with all acceptance-test
   boxes unchecked; tick them as the phases land.

**Commit message:** `docs: ship M16.4 starter prompt (M16.4)`

**Acceptance:** the file exists at
`docs/M16_4_STARTER_PROMPT.md` on `dev/m16.4-fermi-dirac`; CI
green on docs-only diff; no AI-credit lines in the diff.

---

### Phase A: schema surface for FD

Pure-Python only. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.3.0 to 2.4.0. Confirm
   the major-version gate in `semi/schema.py` still accepts
   2.4.0; bump `SCHEMA_SUPPORTED_MINOR` accordingly. The
   `examples` array at the top of the schema gains `"2.4.0"`
   ahead of the existing `"2.3.0"` / `"2.2.0"` / `"2.1.0"` /
   `"2.0.0"` entries.
2. Widen `physics.statistics`:
   - The existing one-element enum `["boltzmann"]` becomes
     `["boltzmann", "fermi_dirac"]`. Default stays `"boltzmann"`.
   - Rewrite the description: "Carrier statistics. boltzmann
     (default) is the V&V reference and is bit-identical to
     v0.19.0. fermi_dirac substitutes the Blakemore
     approximation for `F_{1/2}` in the Poisson source and the
     Slotboom carrier-density helpers; required for quantitative
     I-V at heavy doping (>~1e19 cm^-3) where Boltzmann breaks."
3. Default-fill: an input with no `physics.statistics`, or with
   `statistics == "boltzmann"`, produces an exact bit-equivalent
   solve to v0.19.0.
4. Pure-Python tests in `tests/test_schema.py` and a new
   `tests/test_statistics_schema.py`:
   - Every existing benchmark JSON validates unchanged against
     v2.4.0.
   - `physics.statistics: "fermi_dirac"` validates.
   - `physics.statistics: "foo"` is rejected.
   - Default-fill leaves `statistics` as `"boltzmann"` when
     unset.
5. Update `docs/schema/reference.md` versioning table with the
   v2.4.0 row; the change-log column reads "additive: Fermi-
   Dirac statistics gated by physics.statistics dispatch
   (M16.4)".

**Commit message:** `feat(schema): physics.statistics fermi_dirac dispatch (schema 2.4.0); boltzmann default unchanged (M16.4)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change is observable on any benchmark.

---

### Phase B: Blakemore helpers and the FD-vs-Boltzmann reference

Pure-Python statistics module plus the UFL helpers that consume
it. Layer 3 module is dolfinx-free; Layer 4 helpers in
`semi/physics/slotboom.py` defer UFL imports.

1. New module `semi/physics/statistics.py`:
   - Module docstring: cite Blakemore 1982 and the Sentaurus
     device manual for the approximation
     `F_{1/2}(eta) ~ 1 / (exp(-eta) + zeta(eta))` with
     `zeta(eta) = 3 sqrt(pi/2) (eta + 2.13 + (|eta - 2.13|^2.4
     + 9.6)^(5/12))^{-3/2}`. The simpler-and-good-enough
     form for the production path is the basic Blakemore:
     `F_{1/2}(eta) ~ 1 / (exp(-eta) + 0.27)`. Document
     accuracy: <1 % vs the full integral for `eta < 5`,
     ~5 % at `eta = 7`. Heavy doping in the M16.4 acceptance
     benchmark sits at `eta ~ 4-5`, well within the
     accurate-Blakemore regime.
   - Public API:

     ```python
     def fermi_dirac_half_blakemore(eta):
         """
         Basic Blakemore approximation:
             F_{1/2}(eta) ~ 1 / (exp(-eta) + 0.27)
         Pure NumPy. Pass scalar or array. Returns same shape.
         """

     def fermi_dirac_half_reference(eta):
         """
         Reference F_{1/2}(eta) via scipy.special.fdk (k=1/2)
         where available, falling back to mpmath.polylog(1.5, -exp(eta))
         scaled by -1 for the verification path. Used by the unit
         tests and the diode_fermi_dirac_1d benchmark verifier;
         not called from the production residual.
         """

     def gamma_n_blakemore(psi_minus_phi_n_over_Vt, eta_offset):
         """
         FD-correction prefactor in the generalized-Slotboom form:
             n = n_i * gamma_n * exp(psi - phi_n)
         where gamma_n -> 1 as the system becomes non-degenerate
         (Boltzmann limit) and < 1 in degenerate regions.
         Returns a NumPy scalar/array. The eta_offset argument
         absorbs the (E_C - E_F_intrinsic) / kT contribution from
         the band-edge convention; it is a per-material constant
         computed once at scaling time. See ADR 0004 derivation
         note in PHYSICS.md § 1.2.
         """

     def gamma_p_blakemore(phi_p_minus_psi_over_Vt, eta_offset):
         """Hole counterpart of gamma_n_blakemore."""

     def einstein_factor_blakemore(eta):
         """
         g(eta) = F_{1/2}(eta) / F_{-1/2}(eta), the FD Einstein
         correction. Used only in the unit tests and the algebraic
         cross-check that gamma_n / einstein_factor_blakemore
         cancels in the generalized-Slotboom current expression
         (see ADR 0004 derivation note).
         """
     ```
   - Private helper `_eta_offset_for_material(N_C, n_i)` computing
     the band-edge offset from `n_i` and `N_C` so callers do not
     have to derive it themselves; document the formula
     `eta_offset = ln(n_i / N_C)` (Boltzmann reduces to this).
2. Extend `semi/physics/slotboom.py`:
   - Add a `statistics_cfg: dict | None = None` keyword to
     `n_from_slotboom`, `p_from_slotboom`, `n_from_slotboom_np`,
     `p_from_slotboom_np`. Default `None` is equivalent to
     `{"statistics": "boltzmann"}` (preserves byte-identity).
   - When `statistics_cfg["statistics"] == "fermi_dirac"`:
     multiply the Boltzmann form by the UFL or NumPy version of
     `gamma_n_blakemore` / `gamma_p_blakemore` evaluated at the
     same Slotboom argument. The eta_offset is computed from
     the scaling object and the material `N_C` / `N_V`; thread
     these through `Scaling` (see step 3 below).
   - Document in the module docstring that the Boltzmann path is
     unchanged and that the FD path is the generalized-Slotboom
     form.
3. Extend `semi/scaling.py`:
   - Add `N_C` and `N_V` (effective density of states for
     conduction and valence bands) as fields on the `Scaling`
     dataclass, populated from the material at scaling-build
     time. Si defaults: `N_C = 2.86e19 cm^-3`, `N_V = 3.10e19
     cm^-3`. Read from `semi/materials.py` if a slot exists; add
     it if it does not (Si entry only; the other materials get
     `None` placeholders and are ignored when `statistics ==
     "boltzmann"`).
   - The eta_offset for each carrier is computed from
     `N_C`, `N_V`, and `n_i` and exposed as an attribute on
     `Scaling` (`Scaling.eta_offset_n`, `Scaling.eta_offset_p`).
4. Pure-Python unit tests in `tests/test_statistics.py`:
   - Blakemore-vs-scipy reference: `fermi_dirac_half_blakemore`
     matches `fermi_dirac_half_reference` within 5 % for
     `eta in [-5, 7]`.
   - `gamma_n_blakemore` non-degenerate limit: as
     `(psi - phi_n) / V_t -> -inf` (large-positive
     `eta = -psi - eta_offset_n` for the n-type case), `gamma_n
     -> 1` within 1e-4.
   - `gamma_n_blakemore` degenerate limit: at `eta = 5`,
     `gamma_n` is materially smaller than 1 (cite the value
     against the Blakemore-vs-scipy reference; this is a
     correctness test, not a regression).
   - Einstein-factor cancellation: a closed-form algebraic check
     that `gamma_n * (1/g(eta_n))` simplifies to 1 in the
     Slotboom-current expression. Implementation: pick five
     `(psi, phi_n)` pairs, compute `n_FD` two ways (production
     residual via `n_from_slotboom` and analytical via
     `N_C * F_{1/2}(eta)`), assert agreement within 1 % across
     the range that the Blakemore approximation covers.
   - Default-fill: `n_from_slotboom_np` with `statistics_cfg =
     None` returns the Boltzmann result bit-identically (5+
     sample points; `np.array_equal` not `np.allclose`).

**Commit message:** `feat(physics): Blakemore Fermi-Dirac statistics (M16.4)`

**Acceptance:** new tests pass; no FEM behavior change because
no runner is yet calling the new dispatch (Phase C wires it in).

---

### Phase C: thread the dispatch through the form builders and runners

Plumbing pass. Every form builder and every runner gets a
`statistics_cfg` keyword; the default-None path is bit-identical
to v0.19.0.

1. `semi/physics/poisson.py`:
   - Add `statistics_cfg: dict | None = None` to
     `build_equilibrium_poisson_form` and the multi-region
     variant. When None (or `{"statistics": "boltzmann"}`), the
     existing `rho_hat = ni_hat * (ufl.exp(-psi) - ufl.exp(psi))
     + N_hat_fn` line is unchanged.
   - When `statistics_cfg["statistics"] == "fermi_dirac"`,
     substitute the Blakemore-based n and p expressions
     (n at equilibrium = `n_from_slotboom(psi, 0, ni_hat,
     statistics_cfg=...)`; `phi_n = phi_p = 0` at equilibrium so
     the existing math is `ni_hat * exp(psi)` and `ni_hat *
     exp(-psi)`). Use the new helpers from `semi/physics/slotboom.py`
     directly; do not duplicate the Blakemore math here.
2. `semi/physics/drift_diffusion.py`:
   - Add `statistics_cfg: dict | None = None` to
     `build_dd_block_residual` and `_mr`. Forward to the
     Slotboom helpers at L150-L151 and L309-L316.
3. Runners, all seven, in this order
   (`bias_sweep`, `transient`, `ac_sweep`, `equilibrium`,
   `mos_cv`, `mos_cap_ac`, `resistor_3d`): read
   ```python
   stat_cfg = {"statistics": phys.get("statistics", "boltzmann")}
   ```
   from the config and pass `statistics_cfg=stat_cfg` to every
   form builder call. The `equilibrium` runner is the simplest
   wire-through (Poisson only); the `mos_cap_ac` runner has
   both Poisson and a sensitivity form; the others have full
   Poisson + DD.
4. The full benchmark CI matrix (every existing benchmark with
   the default `statistics: "boltzmann"`) must run bit-identical
   to v0.19.0. This is acceptance test 1 in IMPROVEMENT_GUIDE
   § M16.4. Verify locally with the three numerical anchors
   listed in "Conventions" before pushing.

**Commit message:** `feat(runners): thread statistics_cfg into Poisson and DD form builders (M16.4)`

**Acceptance:** all existing benchmarks run bit-identical to
v0.19.0; no `diode_fermi_dirac_1d` benchmark runs yet (Phase E
adds it).

---

### Phase D: MMS-DD Variant G

ADR 0006 mandates an MMS verifier for every new physics module.
Extend the existing MMS-DD harness with Variant G.

1. Extend `semi/verification/mms_dd.py`:
   - `VARIANTS = ("A", "B", "C", "D", "E", "F", "G")`.
   - Variant G uses the same `psi_e, phi_n_e, phi_p_e`
     manufactured fields as Variant C but evaluates the FD
     substitution (Blakemore) instead of the Boltzmann
     substitution at every callsite. Engineer the manufactured
     amplitudes / `N_C` / `N_V` so the FD-correction prefactor
     `gamma_n` deviates from 1 by ~10-20 % near the manufactured
     extremum (the Blakemore-Boltzmann difference is materially
     exercised; the MMS does not collapse to Variant C).
   - `_build_weak_sources` substitutes the FD versions of
     `n_e` and `p_e` into the manufactured Poisson source and
     the manufactured DD continuity flux.
   - `run_one_level` passes `statistics_cfg = {"statistics":
     "fermi_dirac"}` to `build_dd_block_residual` and
     `build_equilibrium_poisson_form` when `variant == "G"`.
   - CLI study runs Variant G on the same Ns_1d / Ns_2d that
     Variant E and Variant F use, or one level coarser if the
     finest level hits the residual floor (document the choice
     in the docstring).
   - Acceptance: finest-pair L2 rate >= 1.99 and H1 rate >= 0.99
     on each block (psi, phi_n, phi_p). Same gate as Variants D,
     E, F.
2. New file `tests/fem/test_mms_fermi_dirac.py`:
   - Two tests, one for the 1D rate, one for the 2D rate.
     Mirror `tests/fem/test_mms_caughey_thomas.py` (M16.1),
     `tests/fem/test_mms_lombardi.py` (M16.2), and
     `tests/fem/test_mms_auger.py` (M16.3) end to end.
3. Update `docs/PHYSICS.md` § Verification & Validation:
   - Add a "Variant G (Fermi-Dirac statistics, M16.4)" bullet
     under the existing Variant F bullet, summarizing the
     Blakemore-vs-Boltzmann shift target and citing the gate.
   - Append the Variant G rates to the finest-pair-rates table
     once they are measured.
4. Update `docs/mms_dd_derivation.md` with the manufactured-
   source derivation under FD. The derivation should explicitly
   show the generalized-Slotboom substitution and the
   Einstein-factor cancellation (closed-form, not numerical).

**Commit message:** `feat(verification): MMS-DD Fermi-Dirac variant G (M16.4)`

**Acceptance test 3 (MMS rate gate)**: `python
scripts/run_verification.py mms_dd` includes Variant G and
reports L2 rate >= 1.99 at the finest pair on each block.

---

### Phase E: the diode_fermi_dirac_1d benchmark

The analytical anchor: a 1D pn diode at `N_D = 1e20 cm^-3` in
the n+ region where Boltzmann and FD diverge by >15 % on V_bi
and the FD result matches a stand-alone scipy-driven reference
within 1e-3.

1. New directory `benchmarks/diode_fermi_dirac_1d/`:
   - `diode_fermi_dirac.json`: 1D pn diode at N_A = 1e17 cm^-3
     (p-side), N_D = 1e20 cm^-3 (n+ side; heavy enough to push
     the n+ region into degenerate FD territory),
     20 um total, `physics.statistics: "fermi_dirac"`. Two
     companion configs (or one config with a parametric flag):
     `statistics: "boltzmann"` and `statistics: "fermi_dirac"`.
   - `README.md`: device description, the analytical reasoning
     for the >15 % V_bi divergence (cite Sze 3rd ed § 1.5.4 for
     the FD V_bi correction or equivalent textbook reference;
     do **not** cite ChatGPT or any AI assistant).
2. Extend `semi/diode_analytical.py`:
   - Add `vbi_fermi_dirac(N_A, N_D, T, N_C, N_V, n_i)` returning
     the analytical built-in voltage under FD via
     `vbi = (E_F_n - E_F_p) / q` with each Fermi level computed
     from the bulk-doping FD inversion (`F_{1/2}(eta) =
     N / N_C`, solve numerically for `eta`). Pure-Python; no
     dolfinx imports.
   - Add `vbi_boltzmann(N_A, N_D, T, n_i)` returning the
     textbook `V_t * ln(N_A * N_D / n_i^2)` for the comparison.
3. New verifier in `scripts/run_benchmark.py`:
   `verify_diode_fermi_dirac_1d`:
   - Run the FD sweep, run the Boltzmann sweep.
   - Assert `(V_bi_Boltzmann - V_bi_FD) / V_bi_Boltzmann > 0.15`
     (acceptance test 2 part 1).
   - Assert `|V_bi_FD_FEM - V_bi_FD_analytical| / V_bi_FD_analytical
     < 1e-3` (acceptance test 2 part 2). The FEM V_bi is read
     off the equilibrium psi profile at the bulk extrema (same
     pattern `verify_pn_1d` uses).
   - At low doping (the Boltzmann-Boltzmann regression on the
     existing `pn_1d_bias` config), confirm the FD path collapses
     to Boltzmann within 1e-6 relative on V_bi: under
     `physics.statistics: "fermi_dirac"`, run the existing
     `pn_1d_bias` config and assert the FD result is within
     1e-6 of the v0.19.0 Boltzmann result. This is a "did the
     dispatch path actually integrate the FD code" sanity check
     and a regression guard.
   - Mirror the structure of `verify_diode_velsat_1d` (M16.1
     verifier) and `verify_diode_auger_1d` (M16.3 verifier) for
     consistency.
4. Wire the new benchmark into CI: add a step in
   `.github/workflows/ci.yml` matrix near the `diode_auger_1d`
   entry:
   ```yaml
   - name: diode_fermi_dirac_1d
     args: diode_fermi_dirac_1d
   ```
   Document the runtime: roughly 2x `pn_1d_bias` because the
   verifier runs both a Boltzmann and an FD sweep internally.
   **Do not** mark this entry `allow-failure: true`; the
   benchmark must pass cleanly.
5. Notebook (optional but recommended for parity with shipped
   benchmarks): bump number from M16.3 (`07_diode_auger_1d`
   was the optional M16.3 notebook; if it shipped, this is
   `08_diode_fermi_dirac_1d`).

**Commit message:** `feat(benchmark): diode_fermi_dirac_1d FD-vs-Boltzmann verifier (M16.4)`

**Acceptance test 2**: `python scripts/run_benchmark.py
diode_fermi_dirac_1d` exits 0 with the V_bi-divergence and
analytical-match verifier passing.

---

### Phase F: closeout

1. `PLAN.md`:
   - Move M16.4 from "Next task" to "Completed work log" with
     the PR number, deliverables, schema minor bump
     (2.3.0 to 2.4.0), and acceptance-test results (cite the
     observed L2/H1 rates from Phase D, the worst-case
     `|V_bi_FD_FEM - V_bi_FD_analytical| / V_bi_FD_analytical`
     from Phase E, and the Boltzmann-vs-FD V_bi divergence).
   - Set "Next task" to **M16.5 Schottky contacts** with a
     pointer to the (yet-to-be-authored) M16.5 starter prompt.
   - Refresh "Current state" with the new package version
     (bump 0.19.0 to 0.20.0 for the schema minor and the new
     physics dispatch).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M16.4 Done in § 4 with a one-line summary and a
     CHANGELOG anchor.
   - Append a § 9 changelog entry under `[0.20.0]`.
3. `docs/ROADMAP.md`:
   - Update the M16.4 row in the capability matrix from
     Planned to Done; cite the L2 rate, the analytical-match
     worst case, and the Boltzmann-vs-FD V_bi divergence.
4. `CHANGELOG.md`:
   - New `[0.20.0]` entry with the M16.4 line items.
5. `pyproject.toml` and `semi/__init__.py`: bump 0.19.0 to 0.20.0.
6. Push every commit to `origin/dev/m16.4-fermi-dirac`
   immediately. Include the SHA, log line, and `git status`
   snapshot in each gate report.

**Commit message:** `docs: close out M16.4 (PLAN, IMPROVEMENT_GUIDE, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by
      this PR.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact (commits, PR
      description, code, comments, docs, README files).
- [ ] Pure-Python core (`semi/constants.py`, `semi/materials.py`,
      `semi/scaling.py`, `semi/doping.py`, `semi/schema.py`,
      `semi/cv.py`, `semi/timestepping.py`, `semi/bcs.py`,
      `semi/results.py`, `semi/diode_analytical.py`,
      `semi/continuation.py`, `semi/compute.py`,
      `semi/physics/statistics.py`) remains dolfinx-free
      (`tests/test_lazy_imports.py` clean).
- [ ] Boltzmann-default path is bit-identical to v0.19.0 on every
      benchmark. Numerical anchors:
      `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2` (M16.1; reuse),
      `diode_velsat_1d` 56.27 % @ V_F = 0.9 V, 0.19 % @ V_F = 0.3
      V (M16.1 / M16.2 / M16.3; reuse), and `diode_auger_1d`
      >20 % divergence at V_F = 0.9 V (M16.3; reuse).
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004). The generalized-Slotboom FD path
      preserves ADR 0004; if your derivation breaks the
      cancellation argument, ADR 0015 ships first.
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `kronos_server` public API.
- [ ] Schema bumped per minor (2.3.0 to 2.4.0); v2.0.0,
      v2.1.0, v2.2.0, and v2.3.0 inputs still validate.
- [ ] MMS rate gate L2 >= 1.99 and H1 >= 0.99 active in
      `scripts/run_verification.py` for Variants A through G.
- [ ] Coverage gate holds at 95.

## Anti-goals

- Do not start M16.5 (Schottky) in this PR. M16.5 is its own PR
  with its own benchmark; the M16.5 starter prompt is the
  follow-up deliverable for this PR's reviewer.
- Do not change anything in `semi/bcs.py`, `semi/runners/equilibrium.py`'s
  algorithmic structure (only the `statistics_cfg` threading is
  in scope), or `semi/physics/mobility.py`. The FD dispatch
  goes through the statistics module, the Slotboom helpers, the
  Poisson form, the DD form, and all seven runners.
- Do not refactor the Slotboom helpers into a class hierarchy or
  a strategy pattern. The dispatch is a `statistics_cfg` dict
  threaded through; one `if` statement on
  `statistics_cfg.get("statistics", "boltzmann")` per helper.
  Match the M16.1 / M16.3 dispatch style exactly.
- Do not introduce SUPG, streamline diffusion, or a primary-
  density form. ADR 0004 is locked.
- Do not introduce a new ADR unless the generalized-Slotboom
  cancellation argument fails under your derivation. If it does,
  open ADR 0015 first and pause the PR.
- Do not lower the MMS rate gate to "passes on the coarse mesh".
  The gate is L2 >= 1.99 and H1 >= 0.99 at the finest pair,
  every block. If a rate degrades, the implementation has a bug.
- Do not retire the `allow-failure: "true"` flag on the
  `mosfet_2d` CI matrix entry. That carry-over from M16.1 /
  M16.2 / M16.3 is its own follow-up PR.
- Do not bundle this PR with the M14.2.x backlog or with M16.5 /
  M16.6 / M19 / M19.1 / M20. Those are separate PRs.

## Stop conditions

You are done with this prompt when:

1. The M16.4 PR is opened on branch `dev/m16.4-fermi-dirac`.
2. All three acceptance tests in `docs/IMPROVEMENT_GUIDE.md`
   § M16.4 are green in CI:
   - A1: every existing benchmark with `statistics ==
     "boltzmann"` (default) is bit-identical to v0.19.0
     (gated by the existing benchmark CI matrix; cite the
     M16.1, M16.2, and M16.3 anchors verbatim in the PR
     description).
   - A2: `benchmarks/diode_fermi_dirac_1d` shows the
     >15 % Boltzmann-vs-FD V_bi divergence at N_D = 1e20 cm^-3
     and the FD result matches the analytical scipy-driven
     reference within 1e-3.
   - A3: MMS Variant G gate `L2 >= 1.99` finest-pair on each
     block.
3. PLAN.md, IMPROVEMENT_GUIDE.md, ROADMAP.md, and CHANGELOG.md
   reflect the close-out.
4. The package version has bumped 0.19.0 to 0.20.0 in
   `pyproject.toml` and `semi/__init__.py`.
5. PR is reviewed, CI green (modulo the documented `allow-
   failure` on `mosfet_2d`, unchanged in scope from v0.19.0),
   and merged.
6. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI assistant.

## PR description template

```
## Summary

M16.4: Fermi-Dirac statistics (gated), the fourth physics-
completeness slice of the M16 umbrella. Production path uses
the Blakemore approximation `F_{1/2}(eta) ~ 1 / (exp(-eta) +
0.27)`; verification path uses the full integral via
scipy.special.fdk. Generalized-Slotboom substitution preserves
ADR 0004: the continuity rows keep their `J = -q mu n
grad(phi)` shape; only the substitution rule for n in terms of
the primary unknowns changes. No new unknowns. Independent of
M16.2 surface mobility and M16.3 Auger; composes orthogonally
with both.

Schema additive minor bump v2.3.0 to v2.4.0
(physics.statistics enum widened from `["boltzmann"]` to
`["boltzmann", "fermi_dirac"]`). v2.0.0 / v2.1.0 / v2.2.0 /
v2.3.0 inputs continue to validate; statistics=boltzmann (the
default) is bit-identical to v0.19.0.

New benchmark: diode_fermi_dirac_1d (1D pn diode, N_A = 1e17
cm^-3, N_D = 1e20 cm^-3) demonstrates >15 % Boltzmann-vs-FD V_bi
divergence and matches the analytical scipy-driven FD reference
within 1e-3.

New MMS variant: Blakemore-FD substitution manufactured
solution, Variant G. Finest-pair L2 rate >= 1.99 and H1 rate
>= 0.99 on each block.

## Acceptance tests

(All three numbered in docs/IMPROVEMENT_GUIDE.md § M16.4.)

- [ ] A1: every existing benchmark with statistics=boltzmann
      (default) is bit-identical to v0.19.0
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
      diode_auger_1d >20 % divergence at 0.9 V)
- [ ] A2 (divergence): diode_fermi_dirac_1d Boltzmann vs FD
      V_bi at N_D = 1e20 cm^-3 differs by > 15 %
      (observed: <fill in>)
- [ ] A2 (analytical match): diode_fermi_dirac_1d FD V_bi
      within 1e-3 of analytical scipy reference
      (observed worst: <fill in>)
- [ ] A3 (MMS gate): scripts/run_verification.py mms_dd
      Variant G L2 >= 1.99 finest-pair on each block
      (observed: <fill in>)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark diode_fermi_dirac_1d
- [ ] docker compose run --rm benchmark pn_1d_bias
      (boltzmann default, bit-identical)
- [ ] docker compose run --rm benchmark diode_velsat_1d
      (boltzmann default, bit-identical)
- [ ] docker compose run --rm benchmark diode_auger_1d
      (boltzmann default, bit-identical)

## Notes

- The mosfet_2d CI matrix entry remains `allow-failure: "true"`,
  unchanged in scope from M16.3. Retiring that flag is a separate
  follow-up.
- ADR 0004 is preserved via the generalized-Slotboom path; the
  Einstein-factor cancellation is verified algebraically in the
  unit tests and numerically in the MMS Variant G gate.
- Si N_C = 2.86e19 cm^-3, N_V = 3.10e19 cm^-3 added to
  semi/scaling.py and semi/materials.py; other materials carry
  None placeholders and are ignored under statistics=boltzmann.
```

## Hand-off

When M16.4 lands and is reviewed, the next pickup is M16.5
(Schottky contacts), authored at the time by the contributor in
the same shape as this prompt. M16.5 touches `semi/bcs.py`
(adding a `type: "schottky"` contact and a Robin-style boundary
form on the continuity rows) and is the first M16.x slice to
extend the BC layer; the previous M16.x slices stayed inside the
form builders and the schema. Subsequent milestones (M16.6
BBT/TAT, M19, M19.1, M20) each ship as their own PR with their
own starter prompt.
