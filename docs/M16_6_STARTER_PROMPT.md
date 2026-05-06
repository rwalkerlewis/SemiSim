## Context

You are working in `kronos-semi` at `v0.21.0` (post-M16.5, package
version `0.21.0`, schema `2.5.0`). The repo is at
`https://github.com/rwalkerlewis/kronos-semi`; `main` is at
`12e2110` (`M16.5: Schottky contacts (#82)`).

Your assignment is **M16.6: BBT and TAT tunneling**, the sixth and
largest physics-completeness slice of the M16 umbrella. M16.1
(Caughey-Thomas mobility), M16.2 (Lombardi surface mobility), M16.3
(Auger recombination), M16.4 (Fermi-Dirac statistics), and M16.5
(Schottky contacts) have merged. M16.6 is "the largest single physics
addition in M16" per IMPROVEMENT_GUIDE § M16.6: it ships two
generation/recombination kernels (Kane band-to-band tunneling and
Hurkx trap-assisted tunneling), a `zener_1d` reverse-breakdown
benchmark, and depends on M16.4 because BBT requires the FD-corrected
density of states.

`PLAN.md` § "Next task" names M16.6 explicitly:

> **M16.6: BBT and TAT tunneling** on a fresh branch
> `dev/m16.6-tunneling`. Acceptance tests in
> `docs/IMPROVEMENT_GUIDE.md` § M16.6. M16.6 is the largest single
> physics addition in M16; reverse-breakdown benchmark
> (`benchmarks/zener_1d`) at heavy doping (1e19) within 20 % of a
> Kane analytical reference. Depends on M16.4 (BBT density of
> states uses the same Blakemore-FD path).

This prompt **is** the starter. Save it as
`docs/M16_6_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below).

This prompt does not restate the M16.6 deliverable, the acceptance
tests, or the rationale; those live in `docs/IMPROVEMENT_GUIDE.md`
§ M16.6. Read the guide first; this prompt only tells you the order
in which to execute and which invariants must remain load-bearing.

## Branch and PR rules

- Work on a fresh branch off `main`: `git checkout -b dev/m16.6-tunneling`.
  Do **not** rebase or commit onto `main` directly.
- One milestone, one PR. Do not bundle M16.6 with M16.7 (transient
  FFT vs AC validation) or any other physics work.
- Push every phase commit to `origin/dev/m16.6-tunneling`
  immediately after it lands locally. Open the PR after Phase 0.
- Title the PR `M16.6: BBT and TAT tunneling`. Use the PR
  description template at the bottom of this prompt.
- **No AI-assistant credits in any shipped artifact.** Commit
  messages, PR description, code comments, docs, README files,
  the § 9 changelog entry, and the CHANGELOG `[0.22.0]` entry do
  not mention Claude, Claude Code, Anthropic, or any AI tool. If
  Claude Code is configured to auto-append a `Co-Authored-By:`
  trailer, **disable that before the first commit**. The
  per-prompt prohibition has been violated on M16.4 and M16.5
  squash bodies despite explicit prohibitions; the fix is at the
  Claude Code settings layer (`includeCoAuthoredBy: false` in
  `~/.claude/settings.json` or `.claude/settings.json` in the
  project), not at the prompt-text layer. Verify with
  `git log --format=%B dev/m16.6-tunneling` before opening the
  PR; if any trailer is present, fix the settings and re-author
  the commits before pushing.

## Required reading (in order; ~70 minutes; this is the heaviest
M16.x slice)

1. `PLAN.md` in full. Confirm `main` is at v0.21.0 (post-M16.5)
   and that no other PR is in flight on `dev/m16.6-*`. Read the
   M16.5 entry in "Completed work log" (top of the log) end to
   end so you know what surface area is freshly settled (the
   ContactBC `barrier_height_eV` field, the `schottky_facets`
   parameter on the DD form builder, and the M16.5 coverage-
   recovery follow-up that landed inside the M16.5 PR).
2. `docs/IMPROVEMENT_GUIDE.md` § M16 (umbrella context),
   § M16.1 / § M16.2 / § M16.3 / § M16.4 / § M16.5 (just-shipped
   milestones), § M16.6 (Why / Deliverable / Acceptance /
   Dependencies), and § 1 (Honest current state).
3. `docs/ROADMAP.md` § Capability matrix. M16.6 row moves from
   Planned to shipped at the end of this PR. M16.7 then becomes
   the only remaining M16 row before the umbrella closes.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M16.6 touches
   the schema (Layer 2), pure-Python kernel helpers (Layer 3 for
   the closed-form Kane and Hurkx math; the existing
   `semi/physics/recombination.py` module already mixes
   pure-Python and UFL helpers without dolfinx at module scope),
   the DD form builder (Layer 4, for the inline kernels), and
   the benchmark / verifier (Layer 5).
5. `docs/PHYSICS.md`:
   - § 1.4 (recombination kernel; SRH then Auger). The Kane and
     Hurkx kernels are added as additive generation/recombination
     terms next to Auger.
   - § Verification & Validation. Variants A through G live in
     `semi/verification/mms_dd.py`. M16.6 adds **Variant H**
     (BBT and TAT generation kernels in the manufactured weak
     source). Unlike M16.5 (boundary physics, ADR 0015 carve-
     out), M16.6 is domain physics, so the per-module MMS rule
     from ADR 0006 still binds.
6. `docs/adr/`:
   - 0001 (JSON contract), 0002 (nondimensionalization), 0004
     (Slotboom DD; preserved). 0006 (V&V strategy; the per-
     module MMS rule applies to M16.6).
   - 0007 (BC interface; not touched).
   - 0015 (Schottky as Robin BC; the boundary-physics carve-
     out does not apply to M16.6 because Kane and Hurkx are
     domain kernels, not BCs).
7. `docs/mms_dd_derivation.md`. Variants F (Auger) and G (FD
   statistics) are the closest structural analogs. Variant H
   adds two new generation-rate terms in the manufactured weak
   source; the manufactured solution does not need new fields,
   just new kernel evaluations at the manufactured Slotboom
   densities.
8. `semi/physics/recombination.py` end to end. The current
   module ships SRH (`srh_rate`, `srh_rate_np`) and Auger
   (`auger_rate`, `auger_rate_np`) as both UFL and NumPy
   helpers. M16.6 adds the same shape for Kane and Hurkx
   (UFL + NumPy + scaled-coefficient helpers).
9. `semi/physics/drift_diffusion.py` L221-L240 (the Auger inline
   block at `R = R_SRH + R_Auger` when `recomb_cfg["auger"]` is
   on; L613-L632 in `_mr`). The Kane and Hurkx terms inline
   alongside, gated by `recomb_cfg["bbt"]` and
   `recomb_cfg["tat"]`. Kane is a *generation* term (negative
   sign in R for forward operation; positive in R under reverse
   bias where it dominates); Hurkx is an enhancement of SRH in
   high-field regions. Read the Auger inline carefully; the BBT
   and TAT inlines mirror the structure but use distinct
   sign conventions (see Phase B step 1 below).
10. `semi/physics/statistics.py`. The Kane prefactor depends on
    the FD-corrected density of states at the band edges; reuse
    the M16.4 helpers `fermi_dirac_half_blakemore`,
    `gamma_n_blakemore`, `gamma_p_blakemore` rather than
    re-deriving the FD math.
11. `semi/runners/bias_sweep.py` L80-L95. The runner already
    threads `recomb_cfg` to the form builder (M16.3); adding
    BBT and TAT keys is a one-block extension.
12. `benchmarks/diode_auger_1d/` (M16.3) and
    `benchmarks/schottky_1d/` (M16.5) as structural templates
    for `benchmarks/zener_1d/`.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** The new schema entries
  `physics.tunneling: {bbt: bool, tat: bool}` plus the model
  parameters (Kane coefficients A_kane, B_kane; Hurkx
  coefficients tau_n_min, tau_p_min, F_kT, alpha) must be
  expressible in `schemas/input.v2.json` (current strict v2.5.0
  from M16.5), validated by `semi/schema.py`, and exercised by
  at least one benchmark JSON.
- **Schema versioning is binding.** Additive minor bump v2.5.0 to
  v2.6.0. v2.0.0 through v2.5.0 inputs must continue to validate
  and produce bit-identical results to v0.21.0. Update
  `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Five layers, enforced.** The Kane and Hurkx kernels live in
  `semi/physics/recombination.py` (Layer 4 with deferred UFL
  imports). The closed-form NumPy versions are pure-Python; the
  UFL versions defer `import ufl` to the function body. Match
  the M16.3 Auger pattern exactly.
- **Slotboom variables for DD.** ADR 0004 is locked. Both
  kernels enter the continuity rows as scalar source terms
  evaluated at Slotboom-derived `n_hat`, `p_hat`, exactly like
  Auger. No SUPG, no streamline diffusion, no primary-density
  form.
- **Every new physics module needs an MMS verifier.** ADR 0006.
  No exceptions for domain physics. Acceptance L2 rate >= 1.99
  and H1 rate >= 0.99 at the finest pair, gated in
  `scripts/run_verification.py mms_dd`. The boundary-physics
  carve-out from ADR 0015 does **not** apply to M16.6.
- **Every new analytical model needs a benchmark.** `zener_1d`
  is the analytical anchor; reverse-bias breakdown current
  matches the closed-form Kane reference within 20 % from
  V_R = 4 V to V_R = 8 V on a heavily-doped (1e19) abrupt
  junction.
- **Per-runner threading.** `recomb_cfg` already threads through
  `bias_sweep`, `transient`, and `ac_sweep` (per the M16.3
  Auger PR). The new `bbt` and `tat` keys land alongside
  `auger`. Equilibrium runner does not consume recomb (no
  continuity rows). `mos_cv`, `mos_cap_ac`, and `resistor_3d`
  thread recomb through but neither tunneling kernel is
  meaningful for those benchmarks (low fields, no breakdown
  regime); the keys default to False and the kernels do not
  fire.
- **One milestone, one PR.** Do not bundle with M16.7
  (transient FFT vs AC) or any other physics work. Do not
  retire the `mosfet_2d` `allow-failure: "true"` carry-over
  from M16.1-M16.5.
- **No em dashes in prose or code comments.** Use commas,
  periods, parentheses, or colons. Re-grep new diffs after each
  phase.
- **No AI-assistant credits.** Reiterated; see Branch and PR
  rules above. The Claude Code settings fix
  (`includeCoAuthoredBy: false`) is the durable solution.
- **Coverage gate is 95 on `semi/`.** M16.5 needed a coverage-
  recovery follow-up commit (the schottky pre-solve in
  `bias_sweep.py` was uncovered by the gated suite because
  `schottky_1d` runs in a separate CI job whose coverage is
  not merged). **Same risk applies to `zener_1d` if you wire
  Kane / Hurkx pre-solve hooks into bias_sweep.** Plan unit
  tests **before** Phase E that exercise the BBT / TAT
  branches via a tiny config in the gated `docker-fem-tests`
  job (mirror M16.5's `tests/fem/test_schottky_surface.py`
  approach, which is what eventually closed the M16.5
  coverage gap).
- **Recomb-off byte-identity.** Every existing benchmark with
  `physics.tunneling.bbt == false` and
  `physics.tunneling.tat == false` (the defaults) must produce
  bit-identical results to v0.21.0. Numerical anchors:
  `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2` (M16.1; reuse),
  `diode_velsat_1d` 56.27 % @ V_F = 0.9 V, 0.19 % @ V_F = 0.3 V
  (M16.1; reuse), `diode_auger_1d` >20 %
  SRH-vs-(SRH+Auger) divergence at V_F = 0.9 V (M16.3; reuse),
  `diode_fermi_dirac_1d` Boltzmann-vs-FD V_bi 7.37 % at
  N_D = 1e20 cm^-3 (M16.4; reuse), `schottky_1d` worst-case
  thermionic-emission match within 10 % (M16.5; reuse).

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/`
and `pytest tests/`. Do not advance with red tests.

---

### Phase 0: ship this starter prompt

Pure docs.

1. Save this entire file (the prose you are reading now,
   including the PR template at the bottom) as
   `docs/M16_6_STARTER_PROMPT.md`. Strip nothing; commit
   verbatim. **Do not** add a footer crediting any AI tool.
2. Append a one-line "Author M16.6 starter prompt" entry to
   `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under
   `[Unreleased]`.
3. Push the branch and open the PR with the template at the
   bottom of this file. Verify the squash-body preview has no
   AI-credit trailer before opening.

**Commit message:** `docs: ship M16.6 starter prompt (M16.6)`

**Acceptance:** the file exists at
`docs/M16_6_STARTER_PROMPT.md` on `dev/m16.6-tunneling`; CI
green on docs-only diff; no AI-credit lines anywhere.

---

### Phase A: schema surface for tunneling

Pure-Python only. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.5.0 to 2.6.0. Bump
   `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`. The `examples`
   array gains `"2.6.0"` ahead of the existing entries.
2. Add a new sub-object `physics.tunneling` next to the existing
   `physics.recombination` block:
   - `bbt` (boolean, default `false`): "Whether the Kane band-
     to-band tunneling generation kernel is added to the
     continuity rows. When false, the tunneling generation is
     bit-identical to v0.21.0 (zero contribution)."
   - `tat` (boolean, default `false`): "Whether the Hurkx
     trap-assisted tunneling enhancement of the SRH rate is
     active. When false, SRH reduces to the M14.3 textbook
     form bit-identical to v0.21.0."
   - `A_kane` (number, default `4.0e14`, units `cm^-1 s^-1
     V^-2`, minimum 0): "Kane band-to-band prefactor.
     Si default per Sze § 8.4 (4e14)."
   - `B_kane` (number, default `1.9e7`, units `V/cm`,
     minimum 0): "Kane band-to-band exponent coefficient.
     Si default 1.9e7 V/cm."
   - `tau_n_min`, `tau_p_min` (numbers, default `1.0e-9` s):
     the SRH lifetimes used by the Hurkx enhancement at the
     trap-tunneling-only limit.
   - `F_kT` (number, default `1.4e7`, units `V/cm`): the Hurkx
     characteristic field; "field above which Hurkx
     enhancement is order unity (Hurkx 1992)."
   - `alpha` (number, default `2.0`, dimensionless): the Hurkx
     exponent.
3. Conditional validation in `semi/schema.py`: when
   `physics.tunneling.bbt == true` and
   `physics.statistics == "boltzmann"`, emit a `UserWarning`
   (not an error) at validate time noting that BBT is most
   accurate under FD statistics and that the Boltzmann path
   uses a leading-order correction. The benchmark
   (`zener_1d`) ships with `statistics: "fermi_dirac"` set
   explicitly.
4. Default-fill: an input with no `physics.tunneling` block is
   bit-identical to v0.21.0.
5. Pure-Python tests in `tests/test_schema.py` and a new
   `tests/test_tunneling_schema.py`:
   - Every existing benchmark JSON validates unchanged against
     v2.6.0.
   - `physics.tunneling.bbt: true` validates with default
     parameters.
   - `physics.tunneling.tat: true` validates with default
     parameters.
   - Negative `A_kane` is rejected.
   - Default-fill leaves both flags at `false` when unset.
   - The BBT-with-Boltzmann warning fires; the BBT-with-FD path
     does not warn.
6. Update `docs/schema/reference.md` versioning table with the
   v2.6.0 row; the change-log column reads "additive: BBT and
   TAT tunneling dispatch (M16.6)".

**Commit message:** `feat(schema): physics.tunneling bbt + tat dispatch (schema 2.6.0); existing branches unchanged (M16.6)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change is observable on any benchmark.

---

### Phase B: Kane and Hurkx kernels

Pure-Python plus UFL helpers. Layer 4 module with deferred UFL
imports; closed-form math is dolfinx-free.

1. Extend `semi/physics/recombination.py`. **Sign conventions
   matter; misreading these is the most common Kane/Hurkx bug.**
   - Kane band-to-band:
     ```
     G_BBT = A_kane * |E|^2 / sqrt(E_g) * exp(-B_kane * E_g^(3/2) / |E|)
     ```
     where `|E|` is the magnitude of the local electric field
     (`E = -grad(psi)`), `E_g` is the band gap (Si: 1.12 eV).
     This is a **generation** term: it adds to both n and p
     continuity (a tunneled electron-hole pair appears). In the
     `R = G - U` convention used in `drift_diffusion.py`, BBT
     contributes `-G_BBT` to R.
   - Hurkx trap-assisted:
     ```
     R_TAT = (1 + Gamma(F)) * R_SRH
     Gamma(F) = 2 sqrt(3 pi) * (F / F_kT)^(alpha-1)
                * exp((F / F_kT)^2)
     ```
     where `F` is the local field magnitude. Hurkx **enhances**
     SRH (Gamma >= 0); when both `tat` and `srh` are off,
     R_TAT collapses to zero. When `tat` is on but `srh` is on
     too, the existing SRH inline is replaced by
     `(1 + Gamma) * R_SRH`. The composition with Auger and
     Schottky thermionic emission is additive in the same
     `R = R_SRH * (1 + Gamma_TAT) + R_Auger - G_BBT` shape.
2. Public API additions:

   ```python
   def bbt_rate(E_field_magnitude, E_g, A_kane_hat, B_kane_hat):
       """
       UFL expression for scaled BBT generation rate.
           G_hat = A_kane_hat * |E_hat|^2 / sqrt(E_g_hat)
                   * exp(-B_kane_hat * E_g_hat^(3/2) / |E_hat|)
       Returns a UFL expression scaled by C0/t0. Sign
       convention: positive (this is generation; the caller
       subtracts it from R).
       """

   def bbt_rate_np(E_field_magnitude, E_g, A_kane, B_kane):
       """NumPy counterpart."""

   def hurkx_gamma(F, F_kT_hat, alpha):
       """
       UFL expression for the Hurkx field-enhancement factor
       Gamma. Returns a dimensionless UFL expression.
       """

   def hurkx_gamma_np(F, F_kT, alpha):
       """NumPy counterpart."""

   def scaled_kane_coefficients(A_kane_si, B_kane_si, sc):
       """
       Convert the SI Kane coefficients (cm-based units in the
       JSON) into the dimensionless ratios consumed by the UFL
       builder. Returns a dict of fem.Constant-ready scalars.
       """
   ```
3. Pure-Python unit tests in
   `tests/test_recombination.py`:
   - `bbt_rate_np` zero-field limit: `G -> 0` as `|E| -> 0`.
   - `bbt_rate_np` strong-field scaling: at `|E| = 1.5e6 V/cm`
     and `E_g = 1.12 eV`, the rate matches a manual calculation
     of the Si textbook value within 1 %.
   - `bbt_rate_np` band-gap dependence: doubling `E_g` reduces
     `G_BBT` by orders of magnitude (the
     `exp(-B_kane * E_g^(3/2) / |E|)` factor dominates).
   - `hurkx_gamma_np` zero-field limit: `Gamma(0) = 0`.
   - `hurkx_gamma_np` characteristic-field crossover: at
     `F = F_kT`, `Gamma` is order unity (within a factor of
     two; the Hurkx form is not exactly 1 at F = F_kT).
   - `scaled_kane_coefficients` round-trip: SI -> scaled ->
     implied SI reproduces the input within 1e-12 relative.

**Commit message:** `feat(physics): Kane BBT and Hurkx TAT kernels (M16.6)`

**Acceptance:** new tests pass; no FEM behavior change.

---

### Phase C: wire BBT and TAT into the DD form builder

Inline both kernels next to the Auger inline. The `R_SRH * (1 +
Gamma)` structure replaces the bare `R_SRH` when `tat` is on;
the BBT generation subtracts from R when `bbt` is on.

1. `semi/physics/drift_diffusion.py`:
   - Both `build_dd_block_residual` and
     `build_dd_block_residual_mr` get the BBT and TAT inlines.
     The structure:
     ```python
     R_base = (n_hat * p_hat - ni_hat ** 2) / (
         tau_p * (n_hat + n1) + tau_n * (p_hat + p1)
     )
     # M16.6 TAT: enhance SRH with Hurkx field-enhancement
     if recomb_cfg is not None and recomb_cfg.get("tat", False):
         F_mag = ufl.sqrt(ufl.dot(ufl.grad(psi), ufl.grad(psi))) \
                 * sc.V0 / sc.L0  # convert scaled grad to V/cm
         Gamma = hurkx_gamma(F_mag, F_kT_hat, alpha)
         R_SRH = (1.0 + Gamma) * R_base
     else:
         R_SRH = R_base
     # M16.6 BBT: subtract generation
     if recomb_cfg is not None and recomb_cfg.get("bbt", False):
         G_BBT = bbt_rate(F_mag, E_g_hat, A_kane_hat, B_kane_hat)
         R = R_SRH + R_Auger - G_BBT
     else:
         R = R_SRH + R_Auger
     ```
     (the `F_mag` expression is shared between TAT and BBT;
     compute once if both are on, factor out via a UFL local
     binding.)
   - Add `E_g` to `Scaling` and to the Si entry in
     `semi/materials.py` if not already present (it should be;
     check first). Si default 1.12 eV at T = 300 K.
2. Runner threading: extend the existing `recomb_cfg` extraction
   in `bias_sweep`, `transient`, and `ac_sweep` to read the
   `physics.tunneling.bbt` / `tat` flags and the Kane / Hurkx
   parameters into the same dict that already carries `auger`,
   `C_n`, `C_p`. One-block extension per runner.
3. The `mos_cv`, `mos_cap_ac`, and `resistor_3d` runners thread
   the empty-tunneling defaults (one-line addition each).
   Tunneling does not fire in these benchmarks (low fields, no
   breakdown regime); they remain bit-identical to v0.21.0.

**Commit message:** `feat(runners): thread tunneling flags into DD form; BBT and TAT off by default (M16.6)`

**Acceptance:** all M16.1-M16.5 benchmarks remain bit-identical
to v0.21.0; no `zener_1d` benchmark runs yet (Phase E).

---

### Phase D: MMS-DD Variant H

ADR 0006 mandates an MMS verifier for every new domain-physics
module. Variant H exercises both BBT and TAT kernels in the
manufactured weak source.

1. Extend `semi/verification/mms_dd.py`:
   - `VARIANTS = ("A", "B", "C", "D", "E", "F", "G", "H")`.
   - New module constants `MMS_H_A_KANE_FOR_FORM`,
     `MMS_H_B_KANE_FOR_FORM`, `MMS_H_F_KT_FOR_FORM`,
     `MMS_H_ALPHA_FOR_FORM` engineered so each kernel shifts
     the total recombination rate by O(0.1) at the typical
     manufactured amplitudes. The pattern matches Variant F
     (Auger).
   - Manufactured solution: reuse the Variant C `psi_e,
     phi_n_e, phi_p_e` expressions. The BBT and TAT kernels
     compose additively into the manufactured weak source.
   - `_build_weak_sources` substitutes `R_SRH * (1 + Gamma_TAT)
     + R_Auger - G_BBT` evaluated at the manufactured fields
     into the manufactured weak source.
   - `run_one_level` passes `recomb_cfg = {"auger": True,
     "bbt": True, "tat": True, ...}` when `variant == "H"`.
   - CLI study runs Variant H on the same Ns_1d / Ns_2d as
     Variants F and G.
   - Acceptance: finest-pair L2 rate >= 1.99 and H1 rate >=
     0.99 on each block (psi, phi_n, phi_p).
2. New file `tests/fem/test_mms_tunneling.py`:
   - 1D-rate test, 2D-rate test. Mirror
     `tests/fem/test_mms_auger.py`.
3. Update `docs/PHYSICS.md` § Verification & Validation:
   - Add a "Variant H (BBT + TAT tunneling, M16.6)" bullet
     under Variant G. Append the rates to the finest-pair table
     once measured.
4. Update `docs/mms_dd_derivation.md` with the manufactured-
   source derivation under BBT/TAT. Variant H is purely
   additive (no field substitution change), so the derivation
   is short.

**Commit message:** `feat(verification): MMS-DD tunneling Variant H (M16.6)`

**Acceptance test 2 (MMS gate)**: `python
scripts/run_verification.py mms_dd` includes Variant H and
reports L2 rate >= 1.99 finest-pair on each block.

---

### Phase E: the zener_1d benchmark

The analytical anchor: a 1D abrupt junction at heavy doping
(N_A = N_D = 1e19 cm^-3) where reverse-bias breakdown current
is dominated by Kane band-to-band tunneling and matches the
closed-form Kane integral within 20 % from V_R = 4 V to
V_R = 8 V.

1. New directory `benchmarks/zener_1d/`:
   - `zener.json`: 1D abrupt junction, N_A = N_D = 1e19 cm^-3,
     5 um total, junction at 2.5 um. `physics.statistics:
     "fermi_dirac"` (BBT depends on the FD-corrected DOS;
     Boltzmann is wrong by a factor of ~2 at N=1e19).
     `physics.tunneling: {bbt: true, tat: false}`. V_R sweep
     [-8, 0] V step 0.25 V (33 points). The TAT-on configuration
     ships as a sibling config `zener_with_tat.json` for
     diagnostic comparison; the gating verifier uses the
     BBT-only path.
   - `README.md`: device description, the Kane integral
     reference, the analytical reasoning for the 20 %
     tolerance (the Kane closed-form is itself a leading-order
     approximation; tighter than 20 % would be over-claiming).
     Cite Sze § 8.4 and Hurkx 1992 (no AI-tool citations).
2. Extend `semi/diode_analytical.py`:
   - Add `kane_breakdown_iv(V_R, N, E_g, A_kane, B_kane, ...)`
     returning the closed-form Kane reverse-bias current.
     Pure-Python; the integral has a closed-form approximation
     under abrupt-junction assumptions (depletion-approximation
     field profile, integrate Kane G_BBT over the depletion
     region). Cite Sze § 8.4 in the docstring.
3. New verifier in `scripts/run_benchmark.py`:
   `verify_zener_1d`:
   - Run the BBT-only sweep.
   - Compute the analytical Kane reverse current via
     `kane_breakdown_iv`.
   - Assert
     `|J_FEM(V_R) - J_Kane(V_R)| / J_Kane(V_R) < 0.20` for every
     V_R in `[-8, -4]` V (the breakdown regime; below 4 V the
     leakage is dominated by SRH generation and the Kane
     reference is not the right comparator).
   - Mirror the `verify_diode_velsat_1d` (M16.1) and
     `verify_schottky_1d` (M16.5) patterns.
4. Wire into CI: add a step in `.github/workflows/ci.yml`
   matrix near the `schottky_1d` entry:
   ```yaml
   - name: zener_1d
     args: zener_1d
   ```
   **Do not** mark `allow-failure: true`.
5. Plotter: `plot_zener_1d` produces a semilog-y reverse I-V
   with the FEM and Kane-analytical curves overlaid. Mirror
   `plot_schottky_1d`.
6. **Coverage hook (M16.5 lesson learned).** The pre-solve
   hooks for tunneling in `semi/runners/bias_sweep.py` need
   coverage from the gated `docker-fem-tests` job (not just
   from the `zener_1d` benchmark CI job, which runs separately
   and whose coverage is not merged into the gate). Add a
   small end-to-end test in `tests/fem/test_tunneling_surface.py`
   on a tiny zener_1d-like config that exercises every
   tunneling-on branch of `bias_sweep.py`. Plan this test
   **before** running the full benchmark; it is the same
   structural fix that closed the M16.5 coverage gap in a
   follow-up commit. Doing it here in Phase E avoids the
   follow-up commit M16.5 needed.

**Commit message:** `feat(benchmark): zener_1d Kane reverse-breakdown verifier (M16.6)`

**Acceptance test 1**: `python scripts/run_benchmark.py
zener_1d` exits 0 with the 20 % Kane match passing on every
V_R in [-8, -4] V.

---

### Phase F: closeout

1. `PLAN.md`:
   - Move M16.6 from "Next task" to "Completed work log" with
     PR number, deliverables, schema minor bump
     (2.5.0 to 2.6.0), and acceptance-test results (cite
     observed worst-case `|J_FEM - J_Kane| / J_Kane` from Phase
     E and observed Variant H L2 / H1 rates from Phase D).
   - Update the gap list at L195-200 (M16.5 left only "tunneling
     (BBT or TAT) and transient FFT vs AC sweep validation"; now
     the only remaining M16 gap is M16.7). Replace with:
     ```
     **Physics gaps:** no transient FFT vs AC sweep validation.
     M16.7. (Field-dependent mobility shipped in M16.1; Lombardi
     surface mobility in M16.2; Auger in M16.3; Fermi-Dirac in
     M16.4; Schottky in M16.5; BBT and TAT tunneling in M16.6.)
     ```
   - Set "Next task" to **M16.7 transient FFT vs AC sweep
     validation** (the final M16 slice; no new physics, just a
     V&V cross-check between the M13.1 transient and M14 AC
     paths).
   - Refresh "Current state" with the new package version
     (bump 0.21.0 to 0.22.0).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M16.6 Done in § 4 with a one-line summary and a
     CHANGELOG anchor.
   - Update L81-82 gap line (M16.5 left it at "no tunneling
     (BBT or TAT) and no transient FFT vs AC sweep validation");
     replace with:
     ```
     - **No transient FFT vs AC sweep validation.** M16.7. (M16.1
       Caughey-Thomas mobility, M16.2 Lombardi surface mobility,
       M16.3 Auger, M16.4 Fermi-Dirac, M16.5 Schottky contacts,
       and M16.6 BBT and TAT tunneling have all shipped; see § 4
       for the per-milestone Done entries.)
     ```
   - Append a § 9 changelog entry under `[0.22.0]`. Move the
     existing `[Unreleased]` block into `### Released`.
3. `docs/PHYSICS_INTRO.md`:
   - § 6: append a tunneling bullet:
     ```
     - **Band-to-band and trap-assisted tunneling.** Kane BBT
       generation and Hurkx-enhanced SRH, both as additive
       generation/recombination kernels (M16.6,
       `benchmarks/zener_1d/`).
     ```
     Insert in alphabetical-ish order after the Schottky bullet.
   - § 7: drop the "no tunneling" bullet entirely; the only
     remaining § 7 bullets are impact ionization, hetero-
     junctions, 3D MOSFET, and transient-vs-AC validation.
4. `docs/ROADMAP.md`:
   - Update the M16.6 row in the capability matrix from Planned
     to shipped; cite the worst-case verifier match and the
     observed Variant H rates.
5. `CHANGELOG.md`:
   - New `[0.22.0]` entry with the M16.6 line items.
   - Update the schema-banner comment at the top to mention
     v2.6.0 (`additive minor; M16.6 tunneling dispatch; shipped
     with [0.22.0] below`).
6. `pyproject.toml` and `semi/__init__.py`: bump 0.21.0 to
   0.22.0.
7. Push every commit to `origin/dev/m16.6-tunneling`
   immediately. Verify no AI-credit trailers anywhere
   (`git log --format=%B dev/m16.6-tunneling | grep -i
   'co-authored-by: claude'` returns empty).

**Commit message:** `docs: close out M16.6 (PLAN, IMPROVEMENT_GUIDE, PHYSICS_INTRO, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by
      this PR.
- [ ] No mention of Claude, Claude Code, Anthropic, or any AI
      assistant in any shipped artifact.
- [ ] `git log --format=%B dev/m16.6-tunneling | grep -i
      'co-authored-by: claude'` returns empty before opening
      and merging the PR.
- [ ] Pure-Python core remains dolfinx-free
      (`tests/test_lazy_imports.py` clean).
- [ ] Both-flags-off path is bit-identical to v0.21.0 on every
      benchmark. Numerical anchors as listed in Conventions.
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004).
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `kronos_server` public API.
- [ ] Schema bumped per minor (2.5.0 to 2.6.0); v2.0.0 through
      v2.5.0 inputs still validate.
- [ ] MMS rate gate L2 >= 1.99 / H1 >= 0.99 active for Variants
      A through H.
- [ ] Coverage gate holds at 95 in **the gated docker-fem-tests
      job**, not just in the union of CI matrix coverage. The
      M16.5 follow-up landed because this distinction was
      missed; do not repeat.

## Anti-goals

- Do not start M16.7 (transient FFT vs AC) in this PR. M16.7 is
  the final M16 slice and its own PR.
- Do not change anything in `semi/bcs.py`,
  `semi/physics/mobility.py`, or `semi/physics/statistics.py`.
  The tunneling kernels live in `semi/physics/recombination.py`;
  the Kane BBT reuses the M16.4 FD helpers but does not modify
  them.
- Do not introduce SUPG, streamline diffusion, or a primary-
  density form. ADR 0004 is locked.
- Do not extend ADR 0015's boundary-physics carve-out to cover
  domain physics. M16.6 is domain physics; it gets MMS Variant H.
- Do not lower the 20 % Kane tolerance. The Kane closed form is
  a leading-order approximation; 20 % is the textbook gate, and
  the FEM result should match it within that bound.
- Do not retire the `mosfet_2d` `allow-failure: "true"` flag.
- Do not bundle with M14.2.x backlog or M16.7 / M19 / M19.1 /
  M20.
- Do not add `Co-Authored-By` trailers, generated-by footers,
  or any other AI-credit marker. If Claude Code is auto-
  appending one, fix the settings before the first commit.

## Stop conditions

You are done when:

1. The M16.6 PR is opened on branch `dev/m16.6-tunneling`.
2. Both acceptance tests in `docs/IMPROVEMENT_GUIDE.md` § M16.6
   pass in CI:
   - A1: `benchmarks/zener_1d` matches the Kane analytical
     reverse-bias breakdown current within 20 % from V_R = 4 V
     to 8 V.
   - A2: every existing benchmark with both tunneling flags
     off (default) is bit-identical to v0.21.0.
3. MMS Variant H L2 >= 1.99 / H1 >= 0.99 finest-pair on each
   block.
4. PLAN.md, IMPROVEMENT_GUIDE.md, PHYSICS_INTRO.md, ROADMAP.md,
   CHANGELOG.md reflect the closeout.
5. Package version bumped 0.21.0 to 0.22.0 in `pyproject.toml`
   and `semi/__init__.py`.
6. No commit, PR description, code comment, or doc in this PR
   mentions Claude, Claude Code, Anthropic, or any AI
   assistant. Squash-merge body is also clean.
7. Coverage gate holds at 95 in the gated `docker-fem-tests`
   job, not via a follow-up commit.
8. PR reviewed, CI green (modulo the documented `allow-failure`
   on `mosfet_2d`), merged.

## PR description template

```
## Summary

M16.6: Band-to-band and trap-assisted tunneling, the sixth
and largest physics-completeness slice of the M16 umbrella.
Closed-form additive kernels:
- Kane BBT generation:
    G_BBT = A_kane |E|^2 / sqrt(E_g)
            * exp(-B_kane E_g^(3/2) / |E|)
- Hurkx TAT enhancement:
    R_TAT = (1 + Gamma(F)) * R_SRH
both inlined alongside the Auger inline in the DD block residual
builder. No new unknowns; no change to Slotboom primary form.

Schema additive minor bump v2.5.0 to v2.6.0 (new
physics.tunneling sub-object with bbt / tat boolean flags and
Kane / Hurkx parameters). v2.0.0 through v2.5.0 inputs continue
to validate; both-flags-off path is bit-identical to v0.21.0.

New benchmark: zener_1d (1D abrupt junction, N_A = N_D = 1e19
cm^-3, V_R sweep [-8, 0] V, FD statistics) matches the
closed-form Kane reverse-breakdown current within 20 % from
V_R = 4 V to 8 V.

New MMS variant: BBT + TAT additive source manufactured
solution, Variant H. Finest-pair L2 rate >= 1.99 and H1 rate
>= 0.99 on each block.

## Acceptance tests

(Both numbered in docs/IMPROVEMENT_GUIDE.md § M16.6.)

- [ ] A1 (Kane match): zener_1d FEM reverse-bias current within
      20 % of analytical Kane from V_R = 4 V to 8 V (observed
      worst: <fill in>)
- [ ] A2 (byte-identity): every existing benchmark with both
      tunneling flags off bit-identical to v0.21.0
      (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2;
      diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V;
      diode_auger_1d >20 % divergence at 0.9 V;
      diode_fermi_dirac_1d 7.37 % FD-vs-Boltzmann V_bi at
      N_D=1e20 cm^-3; schottky_1d worst-case <10 %)
- [ ] MMS Variant H: scripts/run_verification.py mms_dd Variant H
      L2 >= 1.99 finest-pair on each block (observed: <fill in>)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95 (in the gated
      docker-fem-tests job, not via separate matrix coverage)
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark zener_1d
- [ ] docker compose run --rm benchmark pn_1d_bias
      (both flags off, bit-identical)
- [ ] docker compose run --rm benchmark diode_velsat_1d
      (both flags off, bit-identical)
- [ ] docker compose run --rm benchmark diode_auger_1d
      (both flags off, bit-identical)
- [ ] docker compose run --rm benchmark diode_fermi_dirac_1d
      (both flags off, bit-identical)
- [ ] docker compose run --rm benchmark schottky_1d
      (both flags off, bit-identical)

## Notes

- The mosfet_2d CI matrix entry remains `allow-failure: "true"`,
  unchanged in scope from M16.5.
- Si Kane defaults: A = 4e14 cm^-1 s^-1 V^-2, B = 1.9e7 V/cm
  per Sze § 8.4.
- Si Hurkx defaults: F_kT = 1.4e7 V/cm, alpha = 2.0 per
  Hurkx 1992.
- BBT depends on FD statistics; the zener_1d benchmark sets
  physics.statistics: "fermi_dirac" explicitly. A
  Boltzmann-with-BBT input emits a UserWarning at validate
  time.
- Coverage gate held at 95 in the gated docker-fem-tests job
  via tests/fem/test_tunneling_surface.py, avoiding the
  follow-up commit pattern from M16.5.
```

## Hand-off

When M16.6 lands, the next pickup is **M16.7 transient FFT vs
AC sweep validation**, the final M16 slice. M16.7 ships no new
physics; it adds a V&V cross-check between the M13.1 transient
and M14 AC paths (the FFT of a transient response at a DC
operating point should match the AC sweep at the same point).
After M16.7 the M16 umbrella closes; PLAN.md and PHYSICS_INTRO
gap lists become "no remaining M16 gaps" and the next
milestone in the roadmap is M17 (heterojunctions) or M19 (3D
MOSFET capstone), whichever the maintainer prioritizes.
