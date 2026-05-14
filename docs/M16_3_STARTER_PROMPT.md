# M16.3 starter prompt: Auger recombination

Paste-ready prompt for Claude in VS Code. Mirrors the convention of
`docs/M9_STARTER_PROMPT.md`, `docs/M15_STARTER_PROMPT.md`,
`docs/M14_3_STARTER_PROMPT.md`, `docs/M16_1_STARTER_PROMPT.md`, and
`docs/M16_2_STARTER_PROMPT.md`. Pick this up on a fresh branch
`dev/m16.3-auger` after M16.2 has merged. Do not bundle with any
other milestone work.

---

## Context

You are working in `SemiSim` at `v0.18.0` (post-M16.2, package
version `0.18.0`, schema `2.2.0`). The repo is at
`https://github.com/rwalkerlewis/SemiSim`; `main` is at `67194c4`
(`M16.2: Lombardi surface mobility (#78)`).

Your assignment is **M16.3: Auger recombination**, the third
physics-completeness slice of the M16 umbrella. M16.1 (Caughey-Thomas
field-dependent bulk mobility) and M16.2 (Lombardi surface mobility)
have merged. Auger does not compose with the mobility models; it is a
parallel slice landing on the recombination kernel
(`semi/physics/recombination.py`) rather than the mobility builder
(`semi/physics/mobility.py`), and it pulls directly from the M14.3
SRH infrastructure.

`PLAN.md` § "Next task" names M16.3 explicitly:

> **M16.3: Auger recombination** on a fresh branch `dev/m16.3-auger`.
> To be picked up via `docs/M16_3_STARTER_PROMPT.md` (yet to be
> authored, in the same shape as `docs/M16_2_STARTER_PROMPT.md`);
> acceptance tests in `docs/IMPROVEMENT_GUIDE.md` § M16.3. Auger does
> not depend on M16.2 surface mobility; it pulls directly from the
> M14.3 SRH infrastructure and adds the closed-form
> `R_Auger = (C_n * n + C_p * p) * (n*p - n_i^2)` to the recombination
> kernel.

This prompt **is** the file PLAN.md is asking you to author. Save it
as `docs/M16_3_STARTER_PROMPT.md` as the first commit of the PR
(Phase 0 below), preserving the convention M16.2 inherited from M16.1
and M14.3.

This prompt does not restate the M16.3 deliverable, the acceptance
tests, or the rationale; those live in `docs/IMPROVEMENT_GUIDE.md`
§ M16.3. Read the guide first; this prompt only tells you the order
in which to execute and which invariants must remain load-bearing.

## Branch and PR rules

- Work on a fresh branch off `main`: `git checkout -b dev/m16.3-auger`.
  Do **not** rebase or commit onto `main` directly.
- One milestone, one PR. Do not bundle M16.3 with M16.4 (Fermi-Dirac),
  M16.5 (Schottky), M16.6 (BBT/TAT), or any other physics work. Auger
  is independent of M16.2 Lombardi; do not retouch the Lombardi
  surface area.
- Push every phase commit to `origin/dev/m16.3-auger` immediately
  after it lands locally. Open the PR after Phase 0 so reviewers can
  watch phases land.
- Title the PR `M16.3: Auger recombination`. Use the PR description
  template at the bottom of this prompt.

## Required reading (do not skip; ~40 minutes)

Per `CONTRIBUTING.md` "Before you start", in order:

1. `PLAN.md` in full. Confirm `main` is at v0.18.0 (post-M16.2) and
   that no other PR is in flight on `dev/m16.3-*`. The "Next task"
   section names M16.3; if it does not, stop and align with the
   maintainer before continuing. Read the M16.2 entry in
   "Completed work log" (top of the log) end to end so you know what
   surface area is freshly settled (the Lombardi dispatch and the
   mosfet_2d `allow-failure: "true"` carry-over).
2. `docs/IMPROVEMENT_GUIDE.md` § M16 (umbrella context), § M16.1 and
   § M16.2 (the just-shipped milestones, for context on the schema
   versioning and benchmark-verifier patterns), § M16.3
   (Why / Deliverable / Acceptance / Dependencies), and § 1 (Honest
   current state).
3. `docs/ROADMAP.md` § Capability matrix and § Honest gap. M16.3 is
   listed as Planned; the row moves to shipped at the end of this PR.
4. `docs/ARCHITECTURE.md` for the five-layer rule. M16.3 touches
   physics (Layer 4), the schema (Layer 2), and the benchmark
   (Layer 5). Layer 3 (pure-Python core) and the BC layer
   (`semi/bcs.py`) are **not** touched.
5. `docs/PHYSICS.md`:
   - § 1.4 (recombination kernel; SRH today). The Auger expression
     is added next to the SRH paragraph; cross-reference it.
   - § 2.5 (scaled drift-diffusion; the recombination term enters
     the continuity rows in scaled form `R_hat = R / (C0/t0)`).
   - § Verification & Validation (Variants A through E MMS-DD; you
     are adding Variant F).
6. `docs/adr/` in this order: 0001 (JSON contract), 0002
   (nondimensionalization), 0004 (Slotboom DD; M16.3 must remain in
   Slotboom form because n_hat and p_hat enter the Auger kernel via
   the same `n_from_slotboom` / `p_from_slotboom` helpers SRH already
   uses), 0006 (V&V strategy; M16.3 needs an MMS variant), 0007 (BC
   interface; do not touch).
7. `docs/mms_dd_derivation.md`. Variants D and E added gradient-
   dependent and surface-distance-dependent mobility; Variant F adds
   a recombination kernel that is cubic in carrier density (versus
   SRH's bilinear-over-linear form) at the typical manufactured
   amplitudes.
8. `semi/physics/recombination.py` end to end. The current module
   ships `srh_rate` (UFL), `srh_rate_np` (NumPy), and `scaled_tau`
   (pure Python). The same three-helper pattern extends to Auger.
9. `semi/physics/drift_diffusion.py` lines around the inlined SRH
   `R = (n_hat p_hat - ni_hat^2) / (tau_p (n_hat + n1) + tau_n
   (p_hat + p1))` (currently L155-L160 in
   `build_dd_block_residual` and L318-L323 in
   `build_dd_block_residual_mr`). Note that the SRH expression is
   **inlined** in the DD form builder rather than imported from
   `semi/physics/recombination.py`; the reason is "share the same
   `ni_hat` Constant with the Poisson block and avoid UFL-type
   surprises" (see the comment at L153). Auger has the same `ni_hat`
   requirement; inline the same way (see Phase B).
10. `semi/runners/bias_sweep.py` L84-L90 (the
    `rec = phys.get("recombination", {})` block; tau_n / tau_p / E_t
    are read here and passed through to `build_dd_block_residual` as
    scalars). `semi/runners/transient.py` L204 and
    `semi/runners/ac_sweep.py` L269 have the same
    `rec = phys.get(...)` pattern.
11. `schemas/input.v2.json` L540-L580 (the existing
    `physics.recombination` block). The `auger: bool` field is
    already present as a forward-compat placeholder (default `false`,
    description "Placeholder Auger-recombination flag; not yet
    implemented in the engine"); your PR turns the placeholder into
    a real flag and adds C_n / C_p alongside.
12. `benchmarks/pn_1d_bias/pn_junction_bias.json` (the structural
    template for `diode_auger_1d`; same 1D pn geometry, different
    doping, different sweep).

## Conventions (project rules, not suggestions)

- **JSON is the contract.** The Auger flag and parameters
  (`physics.recombination.auger`, `physics.recombination.C_n`,
  `physics.recombination.C_p`) must be expressible in
  `schemas/input.v2.json` (current strict v2.2.0 from M16.2),
  validated by `semi/schema.py`, and exercised by at least one
  benchmark JSON.
- **Schema versioning is binding.** This is an additive minor bump
  v2.2.0 to v2.3.0 (the field name `auger` already exists as a bool
  placeholder; promoting it to a real flag plus adding C_n / C_p is
  additive, not breaking). v2.0.0 / v2.1.0 / v2.2.0 inputs must
  continue to validate and produce bit-identical results to v0.18.0.
  Update `SCHEMA_SUPPORTED_MINOR` in `semi/schema.py`.
- **Five layers, enforced.** The Auger kernel ships in
  `semi/physics/recombination.py` (Layer 4 if it has any UFL helper;
  the existing module already mixes pure-Python and UFL helpers
  without a `dolfinx` import at module scope, because the UFL
  `srh_rate` does its `import ufl` inside the function body). Match
  that pattern: the new `auger_rate` does `import ufl` inside the
  function body. The `auger_rate_np` pure-Python counterpart needs
  no FEM imports.
- **Slotboom variables for DD.** ADR 0004 is locked. Auger is
  algebraically coupled into the Slotboom-form continuity rows via
  `n_hat = n_from_slotboom(...)`, `p_hat = p_from_slotboom(...)`,
  the same n_hat / p_hat the SRH inline already builds. Do not
  re-derive the continuity rows.
- **Every new physics module needs an MMS verifier.** ADR 0006. No
  exceptions. Acceptance L2 rate >= 1.99 and H1 rate >= 0.99 at the
  finest pair, gated in `scripts/run_verification.py mms_dd`.
- **Every new analytical model needs a benchmark.** The
  `diode_auger_1d` benchmark is the analytical anchor; it must show
  >20 % SRH-vs-(SRH+Auger) divergence at high injection and match an
  analytical high-injection long-diode reference within 5 %.
- **One milestone, one PR.** Do not bundle M16.3 with M16.4 (FD) or
  any other physics work, and do not retire the `mosfet_2d`
  `allow-failure: "true"` flag (that is a separate follow-up,
  unchanged in scope from M16.2).
- **No em dashes in prose or code comments.** Use commas, periods,
  parentheses, or colons. Re-grep new diffs after each phase.
- **Physics-style variable names are allowed.** `C_n`, `C_p`,
  `n_hat`, `p_hat`, `R_Auger`, `R_SRH`, etc.
- **Coverage gate is 95 on `semi/`.** The new Auger paths in
  `semi/physics/recombination.py` must carry pure-Python unit tests
  for the closed-form expressions; the FEM wiring is covered by the
  MMS test and the benchmark.
- **Auger-off byte-identity.** Every existing benchmark with
  `physics.recombination.auger == false` (the default) must produce
  bit-identical results to v0.18.0. This is Acceptance test 1 in
  IMPROVEMENT_GUIDE § M16.3 (worded against v0.15.0 there because
  the section was authored before M16.1 / M16.2; the constant-
  mobility / SRH-only path has been bit-identical from v0.15.0
  through v0.18.0, so reuse the v0.18.0 numerical anchors).
  Numerical anchors for the regression: `pn_1d_bias J(V=0.6 V) =
  1.635e+03 A/m^2` (M16.1 anchor; reuse) and `diode_velsat_1d`
  56.27 % @ V_F = 0.9 V, 0.19 % @ V_F = 0.3 V (M16.1 / M16.2
  anchor; reuse).

## Phases, one commit per phase

Do not bundle. After each phase, run `ruff check semi/ tests/` and
`pytest tests/`. Do not advance with red tests.

---

### Phase 0: ship this starter prompt

Pure docs, no code. Establishes the convention M16.2 inherited from
M16.1 and M14.3, and gives reviewers the contract for the rest of the
PR.

1. Save this entire prompt (the prose you are reading now, including
   the PR template at the bottom) as
   `docs/M16_3_STARTER_PROMPT.md`. Strip nothing; the file is
   committed verbatim.
2. Append a one-line "Author M16.3 starter prompt" entry to
   `docs/IMPROVEMENT_GUIDE.md` § 9 changelog under a new
   `[Unreleased]` heading if one is not already present.
3. Push the branch and open the PR with the template at the bottom
   of this prompt. The PR opens with all acceptance-test boxes
   unchecked; tick them as the phases land.

**Commit message:** `docs: ship M16.3 starter prompt (M16.3)`

**Acceptance:** the file exists at `docs/M16_3_STARTER_PROMPT.md` on
`dev/m16.3-auger`; CI green on docs-only diff.

---

### Phase A: schema surface for Auger

Pure-Python only. No FEM behavior change.

1. Bump `schemas/input.v2.json` minor: 2.2.0 to 2.3.0. Confirm the
   major-version gate in `semi/schema.py` still accepts 2.3.0; bump
   `SCHEMA_SUPPORTED_MINOR` accordingly. The `examples` array at the
   top of the schema gains `"2.3.0"` ahead of the existing `"2.2.0"`
   / `"2.1.0"` / `"2.0.0"` entries.
2. Promote and extend `physics.recombination`:
   - The existing `auger` field stays a `boolean` with default
     `false`; rewrite the description to "Whether the Auger
     recombination kernel `R_Auger = (C_n * n + C_p * p) * (n p -
     n_i^2)` is added to the SRH rate. When false, the kernel is
     bit-identical to v0.18.0."
   - Add `C_n` (number, default `2.8e-31`, units `cm^6/s`,
     description "Electron Auger coefficient. Si default
     (Dziewior-Schmid) is 2.8e-31 cm^6/s; pass the equivalent in
     scaled units to the form builder is handled internally.").
   - Add `C_p` (number, default `9.9e-32`, units `cm^6/s`,
     description "Hole Auger coefficient. Si default
     (Dziewior-Schmid) is 9.9e-32 cm^6/s.").
3. Default-fill: an input with no `physics.recombination.auger`, or
   with `auger == false`, produces an exact bit-equivalent solve to
   v0.18.0. The `C_n` / `C_p` fields are silently filled with
   defaults but not consumed when `auger == false`.
4. Pure-Python tests in `tests/test_recombination.py` (extend the
   existing file; do **not** create a new file):
   - Every existing benchmark JSON validates unchanged against
     v2.3.0.
   - `physics.recombination.auger: true` with default C_n / C_p
     validates.
   - A negative `C_n` is rejected (`minimum: 0` in the schema; add
     this constraint).
   - An extra property under `physics.recombination` is rejected
     (regression of the M14.3 strict-mode guard).
   - Default-fill behavior: unset `auger` lands as `false`; unset
     `C_n` / `C_p` land at the Si defaults.
5. Update `docs/schema/reference.md` versioning table with the
   v2.3.0 row; the change-log column reads "additive: Auger
   recombination kernel (M16.3)".

**Commit message:** `feat(schema): physics.recombination.auger + C_n/C_p (schema 2.3.0); existing branches unchanged (M16.3)`

**Acceptance:** all existing tests pass; new tests pass; no FEM
behavior change is observable on any benchmark.

---

### Phase B: the Auger kernel

The closed-form Auger expression as a NumPy helper, a UFL helper, and
an inline branch in the DD form builder. Layer 4 imports inside
function bodies; pure-Python helpers free of FEM imports.

1. Extend `semi/physics/recombination.py`:
   - Add a module docstring section "Auger recombination (M16.3)"
     stating: the dimensional formula
     `R_Auger = (C_n n + C_p p) (n p - n_i^2)`, the high-injection
     limit `R_Auger -> (C_n + C_p) n^3` (when `n ~ p >> n_i, N`),
     the scaled form derivation, and the units required of
     `C_n_hat` (dimensionless when computed via
     `C_hat = C_SI * C0^2 * t0`; SI `C` has units `cm^6/s`, so
     convert to `m^6/s` via `1e-12` first; cite ADR 0002 for the
     scaling convention).
   - Public API additions:
     ```python
     def auger_rate(n_hat, p_hat, n_i_hat, C_n_hat, C_p_hat):
         """UFL expression for scaled Auger recombination rate."""
     def auger_rate_np(n, p, n_i, C_n, C_p):
         """NumPy counterpart of auger_rate."""
     def scaled_auger_C(C_si, C0, t0):
         """Dimensionless Auger coefficient: C_hat = C_SI * C0^2 * t0"""
     ```
2. Pure-Python unit tests in `tests/test_recombination.py`:
   - `auger_rate_np` low-injection limit
   - `auger_rate_np` high-injection limit: rate approaches
     `(C_n + C_p) X^3` within 1 %
   - `auger_rate_np` equilibrium: when `n p == n_i^2`, rate is
     identically zero
   - `scaled_auger_C` round-trip

**Commit message:** `feat(physics): Auger recombination kernel (M16.3)`

---

### Phase C: wire Auger into the DD form builder and runners

1. `semi/physics/drift_diffusion.py`: both `build_dd_block_residual`
   and `build_dd_block_residual_mr` grow a keyword-only
   `recomb_cfg: dict | None = None` parameter. When `recomb_cfg` is
   None or `recomb_cfg.get("auger", False)` is False, Auger term is
   not added (bit-identical to v0.18.0). When True, inline:
   ```python
   np_minus_nieq = n_hat * p_hat - ni_hat ** 2
   R_SRH = np_minus_nieq / (tau_p * (n_hat + n1) + tau_n * (p_hat + p1))
   R_Auger = (C_n_hat * n_hat + C_p_hat * p_hat) * np_minus_nieq
   R = R_SRH + R_Auger
   ```
2. `semi/runners/bias_sweep.py`, `semi/runners/transient.py`,
   `semi/runners/ac_sweep.py`: read auger parameters from
   `rec = phys.get("recombination", {})` and pass `recomb_cfg`
   through.

**Commit message:** `feat(runners): thread recomb_cfg into DD form builder; Auger off by default (M16.3)`

---

### Phase D: MMS-DD Variant F

ADR 0006 mandates an MMS verifier for every new physics module.

1. Extend `semi/verification/mms_dd.py` with
   `VARIANTS = ("A", "B", "C", "D", "E", "F")` and Variant F logic.
2. New file `tests/fem/test_mms_auger.py` with 1D and 2D rate tests
   mirroring `test_mms_caughey_thomas.py` and `test_mms_lombardi.py`.
3. Update `docs/PHYSICS.md` § V&V with Variant F bullet.
4. Update `docs/mms_dd_derivation.md` with the Variant F
   manufactured-source derivation.

**Commit message:** `feat(verification): MMS-DD Auger variant F (M16.3)`

---

### Phase E: the diode_auger_1d benchmark

1. New `benchmarks/diode_auger_1d/` directory with
   `diode_auger.json` (1D pn diode, N_A = N_D = 1e15 cm^-3, 20 um,
   V_F sweep [0, 0.9] V step 0.05 V) and `README.md`.
2. Extend `semi/diode_analytical.py` with
   `shockley_iv_with_auger(...)`.
3. New verifier in `scripts/run_benchmark.py`: `verify_diode_auger_1d`
   asserting >20 % SRH-vs-Auger divergence at V_F=0.9 V and <5 %
   analytical match.
4. Wire into CI `.github/workflows/ci.yml` (no `allow-failure`).
5. Optional notebook `notebooks/06_diode_auger_1d.ipynb` (or
   `07_*` if `06_*` is taken).

**Commit message:** `feat(benchmark): diode_auger_1d SRH-vs-Auger verifier (M16.3)`

---

### Phase F: closeout

1. `PLAN.md`: move M16.3 to completed log; set next task to M16.4.
2. `docs/IMPROVEMENT_GUIDE.md`: mark M16.3 Done; append `[0.19.0]`
   changelog entry.
3. `docs/ROADMAP.md`: update M16.3 row to Done.
4. `CHANGELOG.md`: new `[0.19.0]` entry.
5. `pyproject.toml` and `semi/__init__.py`: bump 0.18.0 to 0.19.0.
6. Push all commits to `origin/dev/m16.3-auger`.

**Commit message:** `docs: close out M16.3 (PLAN, IMPROVEMENT_GUIDE, ROADMAP, CHANGELOG)`

---

## Invariants checklist (re-verify before each commit)

- No em dashes in any new prose or code comment touched by this PR.
- Pure-Python core remains dolfinx-free.
- Auger-off path is bit-identical to v0.18.0 on every benchmark.
- Slotboom primary unknowns retained; no SUPG / streamline diffusion
  (ADR 0004).
- Schema bumped per minor (2.2.0 to 2.3.0); v2.0.0, v2.1.0, and
  v2.2.0 inputs still validate.
- MMS rate gate L2 >= 1.99 and H1 >= 0.99 active for Variants A
  through F.
- Coverage gate holds at 95.

## Anti-goals

- Do not start M16.4 (Fermi-Dirac) in this PR.
- Do not change anything in `semi/bcs.py`, `semi/scaling.py`,
  `semi/runners/equilibrium.py`, or `semi/physics/mobility.py`.
- Do not refactor the SRH inline in the DD form builder into an
  imported helper.
- Do not introduce a recombination "model" enum.
- Do not lower the MMS rate gate.
- Do not retire the `allow-failure: "true"` flag on `mosfet_2d`.

## Stop conditions

Done when:
1. The M16.3 PR is opened on branch `dev/m16.3-auger`.
2. Both acceptance tests in `docs/IMPROVEMENT_GUIDE.md` § M16.3 are
   green in CI.
3. MMS Variant F gate reports L2 >= 1.99 at the finest pair on each
   block.
4. PLAN.md, IMPROVEMENT_GUIDE.md, ROADMAP.md, and CHANGELOG.md
   reflect the close-out.
5. Package version bumped 0.18.0 to 0.19.0 in `pyproject.toml` and
   `semi/__init__.py`.
6. PR reviewed, CI green, and merged.

## PR description template

```
## Summary

M16.3: Auger recombination, the third physics-completeness slice of the M16 umbrella. Closed-form additive kernel R_Auger = (C_n n + C_p p) (n p - n_i^2) inlined alongside the existing SRH expression in the DD block residual builder. No new unknowns; no change to Slotboom primary form. Independent of M16.2 surface mobility.

Schema additive minor bump v2.2.0 to v2.3.0 (physics.recombination.auger promoted from forward-compat placeholder to a real flag; new C_n / C_p parameters with Si Dziewior-Schmid defaults). v2.0.0 / v2.1.0 / v2.2.0 inputs continue to validate; auger=false branch is bit-identical to v0.18.0.

New benchmark: diode_auger_1d (1D pn diode, N_A = N_D = 1e15 cm^-3, V_F sweep [0, 0.9] V) demonstrates >20 % SRH-only underprediction at V_F = 0.9 V and SRH+Auger matches the analytical high-injection long-diode reference within 5 %.

New MMS variant: cubic-in-density Auger source manufactured solution, Variant F. Finest-pair L2 rate >= 1.99 and H1 rate >= 0.99 on each block.

## Acceptance tests

- [ ] A1: every existing benchmark with auger=false (default) is bit-identical to v0.18.0 (anchors: pn_1d_bias J(V=0.6 V) = 1.635e+03 A/m^2; diode_velsat_1d 56.27 % @ 0.9 V, 0.19 % @ 0.3 V)
- [ ] A2 (divergence): diode_auger_1d SRH-only vs SRH+Auger at V_F = 0.9 V differs by > 20 % (observed: <fill in>)
- [ ] A2 (analytical match): diode_auger_1d SRH+Auger within 5 % of analytical high-injection long-diode at V_F = 0.9 V (observed worst: <fill in>)
- [ ] MMS Variant F: scripts/run_verification.py mms_dd auger variant F L2 >= 1.99 finest-pair on each block (observed: <fill in>)

## Test plan

- [ ] ruff check semi/ tests/
- [ ] pytest tests/
- [ ] pytest --cov=semi --cov-fail-under=95
- [ ] python scripts/run_verification.py all
- [ ] docker compose run --rm benchmark diode_auger_1d
- [ ] docker compose run --rm benchmark pn_1d_bias (auger=false default, bit-identical)
- [ ] docker compose run --rm benchmark diode_velsat_1d (auger=false default, bit-identical)

## Notes

- The mosfet_2d CI matrix entry remains allow-failure: "true", unchanged in scope from M16.2.
- The Auger inline in the DD form builder follows the same pattern as the SRH inline (shared ni_hat Constant; UFL CSE on the (n_hat * p_hat - ni_hat^2) factor).
- Si default C_n = 2.8e-31 cm^6/s, C_p = 9.9e-32 cm^6/s (Dziewior-Schmid; cited in the recombination.py docstring).
```
