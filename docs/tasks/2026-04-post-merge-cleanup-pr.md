# PR: docs: post-merge cleanup, documentation reorg, and loose-ends

This file is the draft body for the PR opened from
`docs/post-merge-cleanup` against `main`. Copy into the GitHub PR
description before requesting review.

## Summary

Documentation, organization, and follow-through pass for the work
that landed in PR #64 (axisymmetric 2D MOSCAP, schema 1.3.0,
MOSCAP analytics, gmsh benchmark, notebook 05 scaffold).

No production code in `semi/` is modified. The post-merge cleanup
covers:

- A real `docs/` tree with `theory/`, `schema/`, `benchmarks/`, and `tasks/` subdirs and a navigable [`docs/index.md`](../../docs/index.md) TOC.
- README slimmed: design notes extracted into `docs/theory/`, inline schema example replaced by a link to `docs/schema/reference.md`, status badges added (CI, license, Python, two Colab notebooks), Day-5 / M14.2 marked done, Verification section refreshed with current pytest counts and headline MOSCAP numbers.
- CHANGELOG reformatted to standard Keep-a-Changelog with preamble and explicit call-out of the schema 1.3.0 bump; existing entries preserved verbatim.
- CONTRIBUTING expanded with axisymmetric, gmsh `.geo`, test split, and benchmark layout sections.
- New regression test `tests/test_moscap_axisym_cv_fem.py` that compares an FEM-extracted C-V to the analytical reference and skips cleanly when the FEM CSV is absent.
- Bidirectional links across README, docs, benchmarks, notebooks, and `semi/` source.
- Historical task prompt moved into `docs/tasks/`.
- Markdown link check verifies all 56 markdown files have valid relative links.

## Acceptance checklist

- [x] Branch `docs/post-merge-cleanup` cut from current `main`. *(CI green is the user's responsibility once they push.)*
- [x] No production code in `semi/` modified beyond docstring/comment updates needed to support the docs reorg. *(No `semi/` files touched in this PR.)*
- [x] README accurately reflects current `main`: M14.2 marked done, current test counts (237 passed / 22 skipped pure-Python; 15/15 MOSCAP analytical anchors green), both Colab badges, status badges row.
- [x] `docs/` directory exists with the prescribed structure; [`docs/index.md`](../../docs/index.md) is a navigable TOC.
- [x] Theory notes ([scaling](../theory/scaling.md), [slotboom](../theory/slotboom.md), [dolfinx_choice](../theory/dolfinx_choice.md), [axisymmetric](../theory/axisymmetric.md), [moscap_cv](../theory/moscap_cv.md)) live in `docs/theory/` and are linked from the README's slimmed-down design-philosophy paragraph.
- [x] [`docs/schema/reference.md`](../schema/reference.md) exists and is more complete than the README's prior inline example.
- [ ] Notebook 05 has been executed end-to-end on Colab; all 7 required figures live in `notebooks/figures/05_moscap_axisym/` and are embedded in the notebook and in `docs/benchmarks/moscap_axisym_2d.md`. **DEFERRED.** See [Open Questions](#open-questions) below; an environment with FEM-on-Colab is required and is not available to this branch's automation.
- [x] [`tests/test_moscap_axisym_cv_fem.py`](../../tests/test_moscap_axisym_cv_fem.py) added; passes (or skips cleanly when `fem_cv.csv` is absent). Today: 1 passed (analytical anchor self-check), 3 skipped.
- [x] [`CHANGELOG.md`](../../CHANGELOG.md) is in Keep-a-Changelog format with preamble, an `[Unreleased]` section that explicitly calls out the schema 1.3.0 bump, and prior version sections preserved verbatim.
- [x] [`CONTRIBUTING.md`](../../CONTRIBUTING.md) covers axisymmetric benchmarks, gmsh `.geo` convention, pure-Python-vs-FEM test split, and benchmark directory layout.
- [ ] Each Open Question from PR #64 has a corresponding GitHub issue with `[follow-up]` prefix and appropriate labels. **Drafted below for the user to file** (per the issues-tracking decision in the cleanup task brief).
- [x] Cross-linking pass complete: README ↔ `docs/index.md` ↔ CHANGELOG ↔ CONTRIBUTING; benchmark READMEs ↔ `docs/benchmarks/` ↔ notebooks; theory notes link to the `semi/` modules they document.
- [x] No broken markdown links anywhere in the repo. (Verified by an in-tree relative-link checker; 56 files scanned, 0 broken.)
- [x] PR description includes the acceptance checklist with each item linked to the relevant file.

## Open Questions

1. **Notebook 05 has not been executed end-to-end.** This branch was assembled in an environment without dolfinx and without Colab automation. The seven required figures (geometry sketch, mesh plot with corner inset, equilibrium $\psi$, carrier-density panels, band diagram, C-V plot vs Hu Fig. 5-18, convergence study) need a working FEM-on-Colab session. The branch-specific Colab URL for executing it is:

   ```
   https://colab.research.google.com/github/rwalkerlewis/SemiSim/blob/docs/post-merge-cleanup/notebooks/05_moscap_axisym_cv.ipynb
   ```

   Once the figures land they should be saved to
   `notebooks/figures/05_moscap_axisym/` and embedded into the notebook
   markdown and into
   [`docs/benchmarks/moscap_axisym_2d.md`](../benchmarks/moscap_axisym_2d.md).

2. **`benchmarks/moscap_axisym_2d/fem_cv.csv` is not produced by this PR.** The new regression test
   [`tests/test_moscap_axisym_cv_fem.py`](../../tests/test_moscap_axisym_cv_fem.py)
   is wired to compare against it but skips cleanly until it lands.

## Follow-up issues to file (drafts)

The user requested the issues be compiled in the PR body for them to
file manually. Each title is prefixed `[follow-up]` per the brief.
Cross-reference PR #64 (the merge that surfaced these) and this PR
(the cleanup) in each body.

### 1. `[follow-up] Rigorous AC small-signal HF method for MOSCAP C-V`
- **Labels:** `physics`, `enhancement`
- **Body:**
  > Today, `compute_hf_cv_depletion_clamp` in
  > [`semi/cv.py`](../../semi/cv.py) implements the textbook
  > depletion-approximation clamp: `C_HF(V_g) = C_LF(V_g)` until
  > `psi_s == 2|phi_B|`, then saturates at `C_min`. This matches Hu
  > Fig. 5-18 to within a few percent for the standard anchor case
  > but is not a true AC simulation.
  >
  > A rigorous HF method would solve `(J + j omega M) du = -(dF/dV) dV`
  > at high omega around the converged DC operating point and report
  > `C(omega) = -Im(Y) / (2 pi f)`. The infrastructure already exists
  > in [`semi/runners/ac_sweep.py`](../../semi/runners/ac_sweep.py)
  > for the two-terminal pn diode case (M14, ADR 0011); extending
  > it to the gate contact of an axisymmetric MOSCAP is the work.
  >
  > Acceptance: at one anchor bias (V_g = V_t + 0.5 V), the AC-extracted
  > C(1 MHz) matches the depletion-clamp `C_min` within a documented
  > tolerance, and the curve transitions smoothly from C_LF at low
  > omega to C_min at high omega.
  >
  > See [`docs/theory/moscap_cv.md`](../theory/moscap_cv.md) for the
  > derivation and the explicit honest-flag on the current
  > implementation. Surfaced by PR #64; tracked in this cleanup PR.

### 2. `[follow-up] Cartesian-2D MOSCAP variant`
- **Labels:** `enhancement`, `physics`
- **Body:**
  > The axisymmetric path in PR #64 is the meridian-half-plane
  > cylindrical formulation. A pure-Cartesian 2D MOSCAP variant
  > (planar gate, side BC for a wide device) would let us directly
  > compare the cylindrical and planar limits and provide a sanity
  > check on the axisymmetric implementation.
  >
  > Acceptance: a `benchmarks/moscap_2d_cartesian/` directory with
  > input JSON, reference CSV, and a verifier; a notebook that
  > overlays the Cartesian and axisymmetric C-V to demonstrate
  > agreement in the planar limit (large gate radius).
  >
  > Surfaced by PR #64; tracked in this cleanup PR.

### 3. `[follow-up] Notebook 05 end-to-end Colab run with figures landed`
- **Labels:** `documentation`
- **Body:**
  > The merged notebook is a scaffold; output cells are empty.
  > Reproducing Hu Fig. 5-18 visually requires running the notebook
  > top-to-bottom on Colab via the FEM-on-Colab installer pattern
  > used by notebook 01.
  >
  > Acceptance: seven figures saved to
  > `notebooks/figures/05_moscap_axisym/` (geometry, mesh + corner
  > inset, equilibrium psi, n/p panels at three biases, band
  > diagram at three biases, LF/HF C-V vs Hu Fig. 5-18 with
  > attribution, convergence vs R and refinement), embedded in the
  > notebook markdown and in
  > [`docs/benchmarks/moscap_axisym_2d.md`](../benchmarks/moscap_axisym_2d.md);
  > a `benchmarks/moscap_axisym_2d/fem_cv.csv` extracted from the
  > FEM run; the existing regression test
  > [`tests/test_moscap_axisym_cv_fem.py`](../../tests/test_moscap_axisym_cv_fem.py)
  > goes from 1-passed-3-skipped to 4-passed.
  >
  > Surfaced by PR #64; tracked in this cleanup PR.

### 4. `[follow-up] Strict input-schema additionalProperties (v2.0.0 schema bump)`
- **Labels:** `enhancement`
- **Body:**
  > The input schema lacks `"additionalProperties": false` on most
  > nested objects, so a UI typo (`"voltag"` instead of `"voltage"`)
  > validates and is silently ignored. The manifest schema is
  > already strict. Flipping the input schema is a v2.0.0 break;
  > worth scheduling alongside the next intentional schema break.
  >
  > See "Known schema caveats" in
  > [`docs/schema/reference.md`](../schema/reference.md).

### 5. `[follow-up] Per-block schema_version`
- **Labels:** `enhancement`
- **Body:**
  > Today there is one top-level `schema_version`. A bump in
  > the `physics` block forces every caller to upgrade atomically.
  > Per-block versioning would allow partial upgrades.
  >
  > Surfaced when M14.2 added `coordinate_system` and forced a
  > minor bump to 1.3.0 even though only axisymmetric users
  > consume the new field.

## Files in this PR

Documentation:
- New: [`docs/index.md`](../index.md), [`docs/theory/scaling.md`](../theory/scaling.md), [`docs/theory/slotboom.md`](../theory/slotboom.md), [`docs/theory/dolfinx_choice.md`](../theory/dolfinx_choice.md), [`docs/theory/axisymmetric.md`](../theory/axisymmetric.md), [`docs/theory/moscap_cv.md`](../theory/moscap_cv.md), [`docs/schema/reference.md`](../schema/reference.md), [`docs/benchmarks/pn_junction_1d.md`](../benchmarks/pn_junction_1d.md), [`docs/benchmarks/moscap_axisym_2d.md`](../benchmarks/moscap_axisym_2d.md), [`docs/tasks/2026-04-axisymmetric-moscap.md`](2026-04-axisymmetric-moscap.md) (renamed from `prompts/`).
- Updated: [`README.md`](../../README.md), [`CHANGELOG.md`](../../CHANGELOG.md), [`CONTRIBUTING.md`](../../CONTRIBUTING.md), [`benchmarks/moscap_axisym_2d/README.md`](../../benchmarks/moscap_axisym_2d/README.md), [`benchmarks/pn_1d/README.md`](../../benchmarks/pn_1d/README.md), [`.gitignore`](../../.gitignore).

Tests:
- New: [`tests/test_moscap_axisym_cv_fem.py`](../../tests/test_moscap_axisym_cv_fem.py).

## Verification

```text
$ python3 -m pytest -q --ignore=tests/test_moscap_axisym_cv.py
238 passed, 25 skipped, 6 deselected, 1 warning in 1.69s

$ python3 -m pytest -q tests/test_moscap_axisym_cv_fem.py
1 passed, 3 skipped in 0.14s
```

Markdown link checker: 0 broken links across 56 tracked `.md` files.

## How this PR was assembled

This PR was assembled locally per the user's chosen workflow
("local commits only, you'll push & open PR"; figures deferred per
the explicit allowance in the cleanup task brief; follow-up issues
compiled here for the user to file). The repository's two-agent
worker/reviewer pattern in
[`CLAUDE.md`](../../CLAUDE.md) was honored: no architectural
decisions were taken, no `semi/` source code modified.
