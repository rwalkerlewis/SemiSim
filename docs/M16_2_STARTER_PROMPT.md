# M16.2 starter prompt: Lombardi surface mobility

Paste-ready prompt for Claude in VS Code. Mirrors the convention of
`docs/M9_STARTER_PROMPT.md`, `docs/M15_STARTER_PROMPT.md`,
`docs/M14_3_STARTER_PROMPT.md`, and `docs/M16_1_STARTER_PROMPT.md`.
Pick this up on a fresh branch `dev/m16.2-lombardi` after M16.1 has
merged. Do not bundle with any other milestone work.

---

## Context

You are working in `SemiSim` at `v0.17.0` (post-M16.1, package
version `0.17.0`, schema `2.1.0`). The repo is at
`https://github.com/rwalkerlewis/SemiSim`; the only commit on
`main` since v0.17.0 landed is the doc-only `0f09973`
(`docs(physics-study-guide): fix equation rendering on GitHub` PR #77),
which contains no code change.

Your assignment is **M16.2: Lombardi surface mobility**, the second
physics-completeness slice of the M16 umbrella. M16.1 (Caughey-Thomas
field-dependent bulk mobility) has merged; M16.2 layers a composite
surface model on top of the bulk dispatch landed in M16.1.

`PLAN.md` § "Next task" names M16.2 explicitly:

> **M16.2: Lombardi surface mobility** on a fresh branch
> `dev/m16.2-lombardi`. To be picked up via
> `docs/M16_2_STARTER_PROMPT.md` (yet to be authored, in the same
> shape as `docs/M16_1_STARTER_PROMPT.md`); acceptance tests in
> `docs/IMPROVEMENT_GUIDE.md` § M16.2.

This prompt **is** the file PLAN.md is asking you to author. Save it
as `docs/M16_2_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below) so that subsequent contributors picking up M16.3 see
the same convention M16.1 inherited from M14.3 and M15.

This prompt does not restate the M16.2 deliverable, the acceptance
tests, or the rationale; those live in `docs/IMPROVEMENT_GUIDE.md`
§ M16.2. Read the guide first; this prompt only tells you the order
in which to execute and which invariants must remain load-bearing.

## Branch and PR rules

- Work on a fresh branch off `main`: `git checkout -b dev/m16.2-lombardi`.
  Do **not** rebase or commit onto `main` directly.
- One milestone, one PR. Do not bundle M16.2 with M16.3 (Auger),
  M16.4 (Fermi-Dirac), or any other physics work, and do not retouch
  the M16.1 Caughey-Thomas surface area beyond the additive
  composition required by Lombardi.
- Push every phase commit to `origin/dev/m16.2-lombardi` immediately
  after it lands locally. Open the PR after Phase 0 so reviewers can
  watch phases land.
- Title the PR `M16.2: Lombardi surface mobility`. Use the PR
  description template at the bottom of this prompt.

## Required reading (do not skip; ~45 minutes)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.17.0 (post-M16.1) and
   that no other PR is in flight on `dev/m16.2-*`. The "Next task"
   section names M16.2; if it does not, stop and align with the
   maintainer before continuing. Read the M16.1 entry in
   "Completed work log" (top of the log) end to end; your work
   composes on top of it.
2. `docs/IMPROVEMENT_GUIDE.md` § M16 (the umbrella context),
   § M16.1 (the just-shipped milestone, for context on the dispatch
   surface you are extending), § M16.2 (Why / Deliverable / Acceptance
   / Dependencies), and § 1 (Honest current state).
3. `docs/ROADMAP.md` § Capability matrix and § Honest gap. M16.2 is
   listed as Planned; the row moves to shipped at the end of this PR.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M16.2 touches
   physics (Layer 4) and the schema (Layer 2) and the benchmark
   (Layer 5). Layer 3 (pure-Python core) and the BC layer
   (`semi/bcs.py`) are **not** touched.
5. `docs/PHYSICS.md`:
   - § 1.5 / § 2.5 (mobility coefficient on the Slotboom flux). The
     Lombardi composite multiplies the same mu0 the M16.1
     Caughey-Thomas model multiplies; the dispatch slot is the same.
   - § Verification & Validation (Variant A through D MMS-DD; you
     are adding Variant E).
   - § 6 (MOS reference); the Pao-Sah verifier window in mosfet_2d
     is the analytical anchor for Acceptance test 2.
6. `docs/adr/` in this order: 0001 (JSON contract), 0002
   (nondimensionalization), 0004 (Slotboom DD; M16.2 must remain in
   Slotboom form), 0006 (V&V strategy; M16.2 needs an MMS variant),
   0007 (BC interface; do not touch).
7. `docs/mms_dd_derivation.md` for the MMS-DD harness. Variant D
   added a gradient-dependent mobility; Variant E adds a surface-
   distance-dependent prefactor on top of the gradient dependence.
8. `semi/physics/mobility.py` end to end. The dispatch
   `build_mobility_expressions` is your extension point. The M16.1
   ADR-style note at the top of the file (the docstring) is the
   project-wide convention for unit-handling commentary in new
   physics modules; the Lombardi block must add an analogous note
   that explains the surface-normal direction, the
   distance-to-interface field, and the unit choices for the
   Lombardi parameters.
9. `semi/verification/mms_dd.py`. The Variant D pattern
   (`MMS_D_VSAT_*_FOR_FORM`, `_build_weak_sources`, `run_one_level`'s
   variant dispatch, the `VARIANTS = ("A", "B", "C", "D")` tuple) is
   the structural template for Variant E.
10. `benchmarks/mosfet_2d/`. The Pao-Sah verifier window in
    `scripts/run_benchmark.py::verify_mosfet_2d` widens upward in
    M16.2 (from `[V_T + 0.2, V_T + 0.6]` V to extend through
    `[V_T + 0.4, V_T + 1.0]` V); the surface-mobility regime is
    where the bulk Caughey-Thomas model alone is no longer good
    enough.

## Conventions (project rules, not suggestions)

- **JSON is the contract.** The new schema entries
  `physics.mobility.model: "lombardi"` and the parameters under
  `physics.mobility.lombardi` (see Phase A) must be expressible in
  `schemas/input.v2.json` (current strict v2.1.0 from M16.1),
  validated by `semi/schema.py`, and exercised by at least one
  benchmark JSON.
- **Schema versioning is binding.** This is an additive minor bump
  v2.1.0 to v2.2.0. v2.0.0 and v2.1.0 inputs (no
  `physics.mobility.lombardi` block; `model` not set to `"lombardi"`)
  must continue to validate and produce bit-identical results to
  v0.17.0. Update `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Five layers, enforced.** Lombardi extensions land in
  `semi/physics/mobility.py` (Layer 4). Imports for `dolfinx`,
  `ufl`, `petsc4py` go inside function bodies, not at module scope.
  The Lombardi parameters live in the JSON schema (Layer 1) and the
  schema loader (Layer 2); pure Python.
- **Slotboom variables for DD.** ADR 0004 is locked. The Lombardi
  composite multiplies the Slotboom flux; do not re-derive the
  continuity rows in primary-density form. Lombardi composes with
  M16.1 Caughey-Thomas via the resistor-sum
  `1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr` per IMPROVEMENT_GUIDE
  § M16.2, where `mu_bulk` is the constant or Caughey-Thomas branch
  selected by a sub-key (see Phase A).
- **Every new physics module needs an MMS verifier.** ADR 0006.
  No exceptions. Acceptance L2 rate >= 1.99 and H1 rate >= 0.99 at
  the finest pair, gated in `scripts/run_verification.py mms_dd`.
- **Every new analytical model needs a benchmark.** mosfet_2d with
  the widened Pao-Sah window is the analytical anchor; the
  verifier threshold in `verify_mosfet_2d` tightens from 20 % to
  10 % in the new window per IMPROVEMENT_GUIDE § M16.2.
- **One milestone, one PR.** Do not bundle M16.2 with M16.1 polish
  (M16.1 has shipped), with M16.3+ (those are separate PRs), or
  with the `mosfet_2d` `allow-failure: "true"` SNES line-search
  follow-up tracked in `.github/workflows/ci.yml` L145-L152. The
  SNES fallback is its own carry-over item; if Lombardi happens to
  unblock the depletion-onset stagnation, that is a happy
  side effect and a separate PR retires the `allow-failure` flag.
- **No em dashes in prose or code comments.** Use commas, periods,
  parentheses, or colons. Re-grep new diffs after each phase.
- **Physics-style variable names are allowed.** `mu_AC`, `mu_sr`,
  `E_perp`, `B_n`, `C_n`, `delta_n`, `lambda_n`, `N_total`, etc.
- **Coverage gate is 95 on `semi/`.** The new Lombardi paths in
  `semi/physics/mobility.py` must carry pure-Python unit tests for
  the closed-form composition math; the FEM wiring is covered by
  the MMS test and the benchmark.
- **Constant-mobility byte-identity.** Every existing benchmark with
  `physics.mobility.model: "constant"` (the default) must produce
  bit-identical results to v0.17.0. This is Acceptance test 3
  (project-wide convention; not enumerated in IMPROVEMENT_GUIDE
  § M16.2 because the same gate applied in M16.1 and was given a
  numerical anchor there: `pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2`,
  byte-identical to pre-M16.1; reuse this anchor verbatim and add
  the M16.2 confirmation alongside it in the PR description).
- **Caughey-Thomas byte-identity.** With `lombardi`-mode off, the
  Caughey-Thomas branch must remain bit-identical to v0.17.0; the
  M16.1 `diode_velsat_1d` divergence/convergence verifier is the
  numerical anchor.

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/` and
`pytest tests/`. Do not advance with red tests.

---

### Phase 0: ship this starter prompt

Pure docs, no code. Establishes the convention M16.1 inherited from
M14.3 and M15 and gives reviewers the contract for the rest of the PR.

1. Save this entire file (the prose you are reading now, including
   the PR template at the bottom) as
   `docs/M16_2_STARTER_PROMPT.md`. Strip nothing; the file is
   committed verbatim.
2. Append a one-line "Author M16.2 starter prompt" entry to
   `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under a new
   `[Unreleased]` heading if one is not already present.
3. Push the branch and open the PR with the template at the bottom
   of this file. The PR opens with all acceptance-test boxes
   unchecked; tick them as the phases land.

**Commit message:** `docs: ship M16.2 starter prompt (M16.2)`

**Acceptance:** the file exists at `docs/M16_2_STARTER_PROMPT.md`
on `dev/m16.2-lombardi`; CI green on docs-only diff.

---

### Phase A: schema surface for Lombardi dispatch

Pure-Python only. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.1.0 to 2.2.0. Confirm the
   major-version gate in `semi/schema.py` still accepts 2.2.0; bump
   `SCHEMA_SUPPORTED_MINOR` accordingly. The `examples` array at the
   top of the schema gains `"2.2.0"` ahead of the existing
   `"2.1.0"` and `"2.0.0"` entries.
2. Extend the existing `physics.mobility` block:
   - Add `"lombardi"` to the `model` enum after `"caughey_thomas"`.
   - Update the `model` field description to enumerate all three
     branches (`constant`, `caughey_thomas`, `lombardi`) and to
     state that Lombardi composes the bulk branch
     (`bulk_model: "constant" | "caughey_thomas"`) with the
     Lombardi surface terms via the resistor sum
     `1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr`.
   - Add a `bulk_model` sub-property (enum
     `["constant", "caughey_thomas"]`, default `"constant"`)
     that selects the bulk branch when `model == "lombardi"`.
     When `model != "lombardi"`, `bulk_model` is ignored
     (document this in the description; do not enforce it through
     conditional schema logic).
   - Add a `lombardi` sub-object with `additionalProperties:
     false` and the following parameters (all numbers, all with
     defaults; SI Lombardi parameter conventions, derived from
     Lombardi 1988 and the Synopsys Sentaurus mobility chapter,
     cited in the Phase B docstring):

     | key        | default     | units      | meaning                                      |
     |------------|-------------|------------|----------------------------------------------|
     | `B_n`      | `4.75e7`    | cm/s       | electron acoustic-phonon prefactor           |
     | `B_p`      | `9.93e6`    | cm/s       | hole acoustic-phonon prefactor               |
     | `C_n`      | `1.74e5`    | cm^(5/3) V^(-2/3) s^(-1) | electron Coulomb-phonon coefficient |
     | `C_p`      | `8.84e5`    | cm^(5/3) V^(-2/3) s^(-1) | hole Coulomb-phonon coefficient     |
     | `lambda_n` | `0.125`     | -          | electron doping exponent                     |
     | `lambda_p` | `0.0317`    | -          | hole doping exponent                         |
     | `delta_n`  | `5.82e14`   | cm^2 V s^-1 | electron surface-roughness coefficient      |
     | `delta_p`  | `2.0546e14` | cm^2 V s^-1 | hole surface-roughness coefficient          |
     | `interface_facet_tag` | `null` (required at solve time when `model == "lombardi"`) | int | Mesh facet tag identifying the Si/SiO2 interface; the surface normal is read from this facet's outward normal at form-assembly time. |

   - Validation rule: when `model == "lombardi"` and
     `interface_facet_tag` is null, the loader raises a clear error
     (see Phase B). The schema declares `interface_facet_tag` as
     `["integer", "null"]` and lets the loader enforce the
     conditional requirement; do not encode this with
     `if/then/else` in JSON Schema (project convention; we keep the
     schema flat and let `semi/schema.py` do the conditional check).
3. Default-fill: an input with no `physics.mobility` block, or with
   `physics.mobility.model in ("constant", "caughey_thomas")`,
   produces an exact bit-equivalent solve to v0.17.0. The
   `bulk_model` and `lombardi` sub-object are silently ignored when
   `model != "lombardi"`.
4. Pure-Python tests in `tests/test_mobility_schema.py`
   (extend the existing M16.1 file; do **not** create a new file):
   - Every existing benchmark JSON validates unchanged against
     v2.2.0.
   - `physics.mobility.model: "lombardi"` with a non-null
     `interface_facet_tag` validates with defaults.
   - `physics.mobility.model: "lombardi"` with a null
     `interface_facet_tag` is rejected by the schema loader (the
     conditional check in `semi/schema.py`, not pure JSON Schema).
   - An unknown model name (`"foo"`) is still rejected (regression
     of the M16.1 test).
   - An extra property under `physics.mobility.lombardi` is rejected
     (`additionalProperties: false`).
   - A `bulk_model` value not in the enum is rejected.
5. Update `docs/schema/reference.md` versioning table with the
   v2.2.0 row; the change-log column reads "additive: Lombardi
   surface mobility (M16.2)".

**Commit message:** `feat(schema): physics.mobility lombardi dispatch (schema 2.2.0); existing branches unchanged (M16.2)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change is observable on any benchmark.

---

### Phase B: the Lombardi composite builder

The closed-form Lombardi expression as a UFL builder, composed with
the M16.1 bulk dispatch via the resistor sum. Layer 4; imports inside
function bodies.

1. Extend `semi/physics/mobility.py`:
   - Add a module-level docstring section "Lombardi surface
     mobility (M16.2)" matching the convention of the M16.1
     section. State: the composite formula (resistor sum), the
     three component formulas with units, the surface-normal
     direction convention (outward from semiconductor into oxide;
     `E_perp = grad(psi) . n_hat` with the sign chosen so that
     `E_perp > 0` in the inversion regime of an n-channel MOSFET),
     and the unit conversions from JSON cm/s and cm^2 V s^-1 into
     the scaled-form units consumed by the UFL builder.
   - Public API additions (no rename of existing helpers):

     ```python
     def lombardi_mu_AC(B, C, N_total, E_perp, lam, T):
         """
         Acoustic-phonon term:
             mu_AC = B / E_perp + C * N_total^lam / (E_perp^(1/3) * T)
         Inputs in scaled units; returns UFL or scalar.
         """

     def lombardi_mu_sr(delta, E_perp):
         """
         Surface-roughness term:
             mu_sr = delta / E_perp^2
         """

     def lombardi_compose(mu_bulk, mu_AC, mu_sr):
         """
         Resistor-sum composition:
             1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr
         """

     def lombardi_unit_conversions(lombardi_cfg, sc):
         """
         Pure-Python conversion of the JSON Lombardi parameters
         from device-physics units (cm/s, cm^2 V s^-1, cm^(5/3)
         V^(-2/3) s^-1) into the dimensionless ratios consumed by
         the UFL builders against
         `E_perp_for_form = sc-derived 1/m`.
         Returns a dict of fem.Constant-ready scalars.
         """
     ```
   - Extend `build_mobility_expressions` to dispatch on
     `model == "lombardi"`. The Lombardi branch:
     1. Reads the `bulk_model` sub-key (default `"constant"`) and
        constructs `mu_bulk_n`, `mu_bulk_p` by recursively calling
        the existing constant or Caughey-Thomas branch with the
        same arguments. **Do not duplicate code**; refactor the
        constant and Caughey-Thomas branches into `_build_constant`
        and `_build_caughey_thomas` private helpers and have the
        Lombardi branch call the right one.
     2. Constructs the perpendicular-field expression. The
        Lombardi model defines `E_perp` as the magnitude of the
        component of `grad(psi)` along the outward semiconductor
        normal at the interface, extended to the bulk via a
        signed-distance heuristic. The simplest implementation
        (and the one this PR ships) reads
        `E_perp = abs(dot(grad(psi), n_hat))` where `n_hat` is the
        FacetNormal restricted to the interface and zero-extended
        into the bulk via a `MeshTags`-conditional UFL coefficient.
        Document the limitation: the surface-mobility branch is
        active only on cells touching `interface_facet_tag`. Cells
        more than one layer away pick up the bulk branch
        unchanged. This matches the spirit of the Lombardi
        boundary-layer model and is sufficient for the inversion-
        regime mosfet_2d acceptance test; an MMS Variant E with a
        fully-3D distance field is out of scope.
     3. Builds `mu_AC_n`, `mu_AC_p`, `mu_sr_n`, `mu_sr_p` via the
        helpers above using
        `N_total = N_A + N_D` (read from the existing doping field
        already on the function space), `lambda_n` / `lambda_p` /
        `B_*` / `C_*` / `delta_*` constants, and the `E_perp`
        expression.
     4. Composes via `lombardi_compose(mu_bulk_*, mu_AC_*, mu_sr_*)`
        per carrier and returns
        `(mu_n_expr, mu_p_expr, "lombardi")`.
     5. If `interface_facet_tag` is null at this dispatch point,
        raise `ValueError("physics.mobility.model='lombardi' "
        "requires physics.mobility.interface_facet_tag")`. Catch
        this in `semi/schema.py::validate` (Phase A) so users see
        the error at `validate(cfg)` time, not deep in the FEM
        path.
2. Pure-Python unit tests in
   `tests/test_mobility_closed_form.py` (extend the existing
   M16.1 file):
   - `lombardi_mu_AC` low-doping low-field limit: `B / E_perp`
     dominates; `C * N_total^lam / (E_perp^(1/3) * T)` is small.
     Sample at three (B, C, N, E, lam, T) tuples.
   - `lombardi_mu_sr` `E_perp -> 0` produces `+inf` (or a
     numerically large value via the same regularization eps as
     M16.1 caughey_thomas; pick `_E_PERP_EPS_SQ = 1e-30`).
   - `lombardi_compose` reduces to `mu_bulk` when
     `mu_AC, mu_sr -> inf`.
   - `lombardi_compose` reduces to `min(mu_AC, mu_sr)` when
     `mu_bulk -> inf`.
   - The sample-point Lombardi-Si-electron mobility at
     `N_total = 1e17 cm^-3, E_perp = 1e6 V/m, T = 300 K, mu0 =
     1400 cm^2/(V s)` matches a reference value computed in the
     test by NumPy with the same formula (the test is a self-check
     on the algebra, not against an external reference; that
     belongs to the mosfet_2d benchmark).

**Commit message:** `feat(physics): Lombardi composite mobility UFL builder (M16.2)`

**Acceptance:** new tests pass; no FEM behavior change on existing
benchmarks (constant default).

---

### Phase C: wire the builder through the runners

Add the dispatch hook in the drift-diffusion form builder; thread
the `interface_facet_tag` through `bias_sweep` and `mos_cap_ac` (the
two runners that drive mosfet_2d-class meshes).

1. `semi/physics/drift_diffusion.py`: the existing `mobility_cfg`
   parameter on `build_dd_block_residual` and `_mr` already
   threads through; no signature change. The new Lombardi branch
   needs the parent `MeshTags` for facets so the UFL FacetNormal
   restricted to `interface_facet_tag` is in scope. Pass the
   `facet_tags` that bias_sweep already constructs into
   `build_mobility_expressions` (extend the signature from
   `(mobility_cfg, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc)`
   to `(mobility_cfg, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0,
   sc, *, facet_tags=None)`). The existing constant and
   Caughey-Thomas branches ignore `facet_tags`; the Lombardi
   branch raises if it is None and the model is Lombardi.
2. `semi/runners/bias_sweep.py`: at runner setup, pass
   `facet_tags=ft` into the DD form builder (the variable already
   exists in the runner). No other change.
3. `semi/runners/mos_cap_ac.py`: same one-line change. The MOSCAP
   AC path also resolves a 2D mesh with an oxide interface; if the
   user requests Lombardi here, the dispatch must honor it. Note
   that the M16.1 starter explicitly excluded `mos_cap_ac` from
   Caughey-Thomas wiring; M16.2 brings it into the Lombardi
   dispatch because the C-V curve in inversion is exactly where
   surface mobility matters. (The M16.1 anti-goal does not bind
   M16.2; the dispatch needed it now.)
4. The other runners (equilibrium, transient, ac_sweep,
   resistor_3d, mos_cv) intentionally do not adopt the Lombardi
   dispatch in this PR. Equilibrium has no continuity rows;
   transient and ac_sweep currently read `mob = phys.get(
   "mobility", {})` but do not thread `mobility_cfg` to the form
   builder (this is documented in PLAN.md M16.1 entry under
   "Phase C"); resistor_3d has no inversion regime; mos_cv uses
   gate-charge integration that is mu-independent. Adopting them
   is M16.2.x or M16.4 work.

**Commit message:** `feat(runners): thread facet_tags into mobility dispatch for Lombardi (M16.2)`

**Acceptance:** all M16.1 acceptance tests still pass (constant
default plus Caughey-Thomas at v0.17.0); no Lombardi-mode benchmark
runs yet (Phase E adds it).

---

### Phase D: MMS-DD Variant E

ADR 0006 mandates an MMS verifier for every new physics module.
Extend the existing MMS-DD harness with Variant E (Variants A through
D from M16.1 are intact).

1. Extend `semi/verification/mms_dd.py`:
   - `VARIANTS = ("A", "B", "C", "D", "E")`.
   - New module constants `MMS_E_*` mirroring the M16.1 M16.1
     `MMS_D_VSAT_*_FOR_FORM` pattern: engineer the dimensionless
     ratios `(B / E_perp_typical)`,
     `(C * N_total^lam / (E_perp_typical^(1/3) * T))`,
     `(delta / E_perp_typical^2)` so each surface term moves the
     composite mobility by a few percent at the typical
     manufactured gradient; the Lombardi branch is materially
     exercised but not so dominant that the MMS forcing is
     numerically pathological. Use the same `O(0.3)` reduction
     target M16.1 used for Variant D as the design anchor.
   - Manufactured solution: reuse the existing Variant C / D
     `psi_e, phi_n_e, phi_p_e` expressions. The surface-normal
     direction in the MMS domain is the unit `e_y` for 1D and
     `e_x` for 2D (the manufactured "interface" is the line
     `y == 0` in 1D and `x == 0` in 2D); the
     `interface_facet_tag` is set on the corresponding mesh
     boundary. This avoids needing a separate signed-distance
     UFL coefficient for the verifier.
   - `_build_weak_sources` substitutes `lombardi_compose(...)`
     evaluated at the manufactured Fermi gradients (and the
     manufactured `E_perp = grad(psi_e) . n_hat`) into the
     manufactured weak source so the forcing matches the
     production residual at `phi_e` exactly. Cite the Variant D
     pattern in the docstring.
   - `run_one_level` passes `mobility_cfg = {"model": "lombardi",
     "bulk_model": "constant", "lombardi": {...},
     "interface_facet_tag": <int>}` to `build_dd_block_residual`
     when `variant == "E"`.
   - CLI study runs Variant E on the same Ns_1d / Ns_2d the M16.1
     Variant D runs on (`Ns_1d = [40, 80]`, `Ns_2d = [32, 64,
     128]`), or a sequence one level coarser if the finest level
     hits the residual floor as Variant D's 1D path did. Document
     the choice in the docstring; the project convention is "as
     few levels as cleanly demonstrate the rate gate, no fewer".
   - Acceptance: finest-pair L2 rate >= 1.99 and H1 rate >= 0.99
     on each block (psi, phi_n, phi_p). Same gate as Variant D.
2. New file `tests/fem/test_mms_lombardi.py`:
   - Two tests, one for the 1D rate, one for the 2D rate. Each
     reads the Variant E results and asserts the gate. Mirror the
     M16.1 file `tests/fem/test_mms_caughey_thomas.py` end to end.
3. Update `docs/PHYSICS.md` § Verification & Validation:
   - Add a "Variant E (Lombardi surface mobility, M16.2)" bullet
     under the existing Variant D bullet, summarizing the
     `O(0.3)` design target and the
     `MMS_E_*` parameter knobs.
   - Append the Variant E rates to the finest-pair-rates table
     once they are measured.
4. Update `docs/mms_dd_derivation.md` with the manufactured-source
   derivation for the Lombardi variant. The derivation should
   cleanly factor out the surface-normal direction so a future
   3D MMS variant can extend it without re-deriving the
   composition.

**Commit message:** `feat(verification): MMS-DD Lombardi variant E (M16.2)`

**Acceptance test 1**: `python scripts/run_verification.py mms_dd`
includes the lombardi variant and reports L2 rate >= 1.99 at the
finest pair on each block.

---

### Phase E: widened mosfet_2d Pao-Sah verifier

The analytical anchor: the mosfet_2d benchmark is re-run with
Lombardi mobility, and the Pao-Sah verifier window widens from
`[V_T + 0.2, V_T + 0.6]` V to `[V_T + 0.4, V_T + 1.0]` V. The
verifier tolerance tightens from 20 % to 10 % in the new window.

1. Extend `benchmarks/mosfet_2d/mosfet_2d.json`:
   - Bump `schema_version` from `"2.0.0"` to `"2.2.0"`.
   - Add a `physics.mobility` block setting `model: "lombardi"`,
     `bulk_model: "caughey_thomas"`, `interface_facet_tag` set to
     the existing facet tag for the silicon/oxide boundary (read
     the value off the existing `mesh.facets_by_plane`; do not
     invent a new tag), and the default Lombardi parameters.
   - Bump `voltage_sweep` upper bound on `V_GS` from `1.5 V` to
     `2.0 V` so the new verifier window
     `[V_T + 0.4, V_T + 1.0]` is fully sampled. The existing
     low-window data points (`V_T - 0.1` through `V_T + 0.6`) are
     retained for regression.
2. Extend `scripts/run_benchmark.py::verify_mosfet_2d`:
   - The verifier dispatches on the `physics.mobility.model` of
     the input config:
     - `"constant"` (the M16.1 path): unchanged window
       `[V_T + 0.2, V_T + 0.6]`, unchanged 20 % tolerance.
     - `"caughey_thomas"`: same window as constant, same 20 %
       (the M16.1 PR did not move this window; it is the M16.2
       PR's job, not M16.1's).
     - `"lombardi"`: widened window `[V_T + 0.4, V_T + 1.0]`,
       tightened 10 % tolerance.
   - The check labels carry the `model` name so the CI log makes
     the dispatch obvious.
3. Update `tests/test_mosfet_2d_verifier.py` (extend, do not
   replace): add three pure-Python assertions exercising the
   `lombardi` branch dispatch (window and tolerance correctness).
4. Update `benchmarks/mosfet_2d/README.md`: a one-paragraph
   addition describing the M16.2 Lombardi configuration and the
   new verifier window. Cite IMPROVEMENT_GUIDE § M16.2.
5. Wire the Lombardi mosfet_2d run into CI: the existing
   `mosfet_2d` matrix entry in `.github/workflows/ci.yml`
   L150-L152 is already gated `allow-failure: "true"` because of
   the M16.1-era SNES line-search stagnation in the depletion
   onset. **Do not** drop the `allow-failure` flag in this PR;
   add a comment indicating M16.2 inherits the same flag and that
   retiring it is a separate follow-up. If your Lombardi run
   passes (it might; the surface mobility softens the gradient at
   the interface and may render the depletion-onset Newton step
   well-conditioned), still leave the flag as-is and note this in
   the PR description; the maintainer will retire the flag in a
   separate cleanup PR after the run is reproducibly green.

**Commit message:** `feat(benchmark): mosfet_2d Lombardi inversion-regime Pao-Sah verifier (M16.2)`

**Acceptance test 2**: `python scripts/run_benchmark.py mosfet_2d`
with `physics.mobility.model: "lombardi"` exits 0 with the
`[V_T + 0.4, V_T + 1.0]` Pao-Sah verifier passing within 10 %.

---

### Phase F: closeout

1. `PLAN.md`:
   - Move M16.2 from "Next task" to "Completed work log" with the
     PR number, deliverables, schema minor bump
     (2.1.0 to 2.2.0), and acceptance-test results (cite the
     observed L2/H1 rates from Phase D and the worst-case
     |I_D_sim - I_D_PaoSah|/I_D_PaoSah from Phase E).
   - Set "Next task" to **M16.3 Auger recombination** with a
     pointer to the (yet-to-be-authored) M16.3 starter prompt and
     a reminder that Auger does not depend on M16.2; it pulls
     directly from M14.3.
   - Refresh "Current state" with the new package version (bump
     0.17.0 to 0.18.0 for the schema minor and the new physics
     model).
2. `docs/IMPROVEMENT_GUIDE.md`:
   - Mark M16.2 Done in § 4 with a one-line summary and a CHANGELOG
     anchor.
   - Append a § 9 changelog entry under `[0.18.0]`.
3. `docs/ROADMAP.md`:
   - Update the M16.2 row in the capability matrix from Planned to
     Done; cite the L2 rate and the worst-case verifier number.
4. `CHANGELOG.md`:
   - New `[0.18.0]` entry with the M16.2 line items.
5. `pyproject.toml` and `semi/__init__.py`: bump 0.17.0 to 0.18.0.
6. Push every commit to `origin/dev/m16.2-lombardi` immediately.
   Include the SHA, log line, and `git status` snapshot in each
   gate report.

**Commit message:** `docs: close out M16.2 (PLAN, IMPROVEMENT_GUIDE, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- [ ] No em dashes in any new prose or code comment touched by this
      PR.
- [ ] Pure-Python core (`semi/constants.py`, `semi/materials.py`,
      `semi/scaling.py`, `semi/doping.py`, `semi/schema.py`,
      `semi/cv.py`, `semi/timestepping.py`, `semi/bcs.py`,
      `semi/results.py`, `semi/diode_analytical.py`,
      `semi/continuation.py`, `semi/compute.py`) remains
      dolfinx-free (`tests/test_lazy_imports.py` clean).
- [ ] Constant-mobility path is bit-identical to v0.17.0 on every
      benchmark. Numerical anchor: `pn_1d_bias J(V=0.6 V) =
      1.635e+03 A/m^2` (M16.1 M16.1 anchor; reuse).
- [ ] Caughey-Thomas path is bit-identical to v0.17.0 on
      `diode_velsat_1d`. Numerical anchor: 56.27 % divergence at
      V_F = 0.9 V, 0.19 % convergence at V_F = 0.3 V.
- [ ] Slotboom primary unknowns retained; no SUPG / streamline
      diffusion (ADR 0004).
- [ ] `make_scaling_from_config` still on every solve path.
- [ ] No PETSc / UFL types leak into `semisim_server` public API.
- [ ] Schema bumped per minor (2.1.0 to 2.2.0); v2.0.0 and v2.1.0
      inputs still validate.
- [ ] MMS rate gate L2 >= 1.99 and H1 >= 0.99 active in
      `scripts/run_verification.py` for Variants A through E.
- [ ] Coverage gate holds at 95.

## Anti-goals

- Do not start M16.3 (Auger) in this PR. Auger is its own PR with
  its own benchmark; the M16.3 starter prompt is the follow-up
  deliverable for this PR's reviewer.
- Do not change anything in `semi/bcs.py`, `semi/scaling.py`, or
  `semi/runners/equilibrium.py`. The Lombardi dispatch goes through
  the bias_sweep and mos_cap_ac runners and the DD form builder.
- Do not introduce a `mu_n` / `mu_p` field as a new primary unknown.
  Lombardi composes closed-form; it does not add unknowns.
- Do not lower the MMS rate gate to "passes on the coarse mesh".
  The gate is L2 >= 1.99 and H1 >= 0.99 at the finest pair, every
  block. If a rate degrades, the implementation has a bug.
- Do not retire the `allow-failure: "true"` flag on the
  `mosfet_2d` CI matrix entry. That is its own follow-up PR
  after the SNES line-search behavior is independently audited.
- Do not bundle this PR with the M14.2.x backlog or with M16.3
  through M16.7 / M19 / M19.1 / M20. Those are separate PRs.
- Do not adopt the Lombardi dispatch in the equilibrium, transient,
  ac_sweep, resistor_3d, or mos_cv runners. Those are M16.2.x or
  M16.4 follow-ups.

## Stop conditions

You are done with this prompt when:

1. The M16.2 PR is opened on branch `dev/m16.2-lombardi`.
2. Both acceptance tests in `docs/IMPROVEMENT_GUIDE.md` § M16.2
   are green in CI.
3. Acceptance test 3 (every existing benchmark with constant
   mobility produces bit-identical results to v0.17.0) is gated by
   the existing benchmark CI matrix; confirm in the PR description
   that the matrix passed and cite the M16.1 `pn_1d_bias` anchor
   verbatim.
4. Acceptance test 4 (Caughey-Thomas branch bit-identical to
   v0.17.0 on `diode_velsat_1d`) confirmed in the PR description.
5. PLAN.md, IMPROVEMENT_GUIDE.md, ROADMAP.md, and CHANGELOG.md
   reflect the close-out.
6. The package version has bumped 0.17.0 to 0.18.0 in
   `pyproject.toml` and `semi/__init__.py`.
7. PR is reviewed, CI green (modulo the documented `allow-failure`
   on `mosfet_2d`, which is unchanged in scope from v0.17.0), and
   merged.

## PR description template

```
## Summary

M16.2: Lombardi surface mobility, the second physics-completeness
slice of the M16 umbrella. Closed-form composite of bulk
(constant or Caughey-Thomas, dispatched via bulk_model sub-key),
acoustic-phonon, and surface-roughness terms via the resistor sum
1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr. No new unknowns; no change
to Slotboom primary form.

Schema additive minor bump v2.1.0 to v2.2.0
(physics.mobility.model gains "lombardi"; new bulk_model and
lombardi sub-objects; new interface_facet_tag conditional). v2.1.0
and v2.0.0 inputs continue to validate; constant and
Caughey-Thomas branches are bit-identical to v0.17.0.

Updated benchmark: mosfet_2d re-runs with lombardi mobility on a
widened V_GS sweep [0.0, 2.0] V. The Pao-Sah verifier window
widens from [V_T + 0.2, V_T + 0.6] V (M14.3) to
[V_T + 0.4, V_T + 1.0] V (M16.2) and the tolerance tightens from
20 % to 10 % in the new window.

New MMS variant: gradient-and-surface-distance-dependent mobility
manufactured solution, Variant E. Finest-pair L2 rate >= 1.99 and
H1 rate >= 0.99 on each block.

## Acceptance tests

(Both numbered in docs/IMPROVEMENT_GUIDE.md § M16.2.)

- [ ] A1: scripts/run_verification.py mms_dd lombardi variant E
      L2 >= 1.99 finest-pair on each block (observed: <fill in
      from CI run>)
- [ ] A2: scripts/run_benchmark.py mosfet_2d (lombardi) within
      10 % of Pao-Sah in [V_T + 0.4, V_T + 1.0] V (observed worst:
      <fill in>)
- [ ] A3: every existing benchmark with constant mobility is
      bit-identical to v0.17.0
      (anchor: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2)
- [ ] A4: caughey_thomas branch bit-identical to v0.17.0 on
      diode_velsat_1d (anchors: 56.27 % @ 0.9 V, 0.19 % @ 0.3 V)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark mosfet_2d (lombardi config)
- [ ] docker compose run --rm benchmark diode_velsat_1d
      (caughey_thomas, bit-identical)
- [ ] docker compose run --rm benchmark pn_1d_bias
      (constant mu, bit-identical)

## Notes

- The mosfet_2d CI matrix entry remains `allow-failure: "true"`;
  retiring that flag is a separate follow-up regardless of whether
  the Lombardi run passes here. Document the observed pass/fail
  status in the PR thread.
- M16.1 anti-goal "do not adopt CT in mos_cap_ac" does not bind
  M16.2; mos_cap_ac picks up the Lombardi dispatch because the
  inversion-regime C-V curve is exactly where surface mobility
  matters.
```

## Hand-off

When M16.2 lands and is reviewed, the next pickup is M16.3 (Auger
recombination), authored at the time by the contributor in the
same shape as this prompt. Subsequent milestones (M16.4 through
M16.6, M19, M19.1, M20) each ship as their own PR with their own
starter prompt; the shape is established by M14_3, M16_1, and
this file.
