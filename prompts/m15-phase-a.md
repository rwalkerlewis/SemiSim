# Phase A — refresh docs for v0.14.x reality (M15 setup)

**Source spec:** `docs/M15_STARTER_PROMPT.md` Phase A.
**Branch:** `dev/m15-phase-a` off `main` at `799934d`.
**Type:** documentation only, no code paths changed.

## Preconditions verified by reviewer

- PR #65 is **merged** (commit `799934d` on `main`). PLAN.md has not yet
  been updated to reflect the merge; it still calls PR #65 "Active dev
  branch" and points "Next task" at it. The starter prompt's gate
  ("stop and wait if PR #65 still open") is therefore satisfied.
- `pyproject.toml` version is `0.14.1`. The starter-prompt expectation
  was `0.14.2` post-PR-#65. The version was not bumped at PR-#65 merge.
  Phase A treats current state as v0.14.1 with PR #65 already on main;
  Phase E (M15 close-out) will bump 0.14.x → 0.15.0. Do not bump in
  Phase A.
- Latest entry in `PLAN.md` "Completed work log" is M14.2 (2026-04-30).
- `docs/IMPROVEMENT_GUIDE.md` ends at §9 (line 611). No §10 exists yet.
- `README.md` contains **no** "Day N" or "Week N" residual references
  (the starter prompt's premise is stale on this point). The only
  surviving day-language is the phrase "Day-1 math helper" at
  README.md line 397 plus filename references to
  `tests/check_day1_math.py` at README.md lines 398 and 414.

## Scope corrections vs starter prompt §Phase A.6

The starter prompt asks the worker to remove "Day N" and "Week N"
references from README and to drop a "(planned for Day 2+)" Slotboom
parenthetical. Neither exists in the current README. Reduced README
scope to:

- Replace "The Day-1 math helper" wording with milestone-anchored
  language at `README.md` L397-L398.
- Update the two `tests/check_day1_math.py` filename references at
  L398 and L414 to the new `tests/check_analytical_math.py`.

The README §Status, §Design notes, and §Roadmap sections already use
M-numbering consistently. Leave them alone.

## Tasks (one commit)

### 1. PLAN.md

- §"Repository" → "Active dev branch": replace the PR #65 sentence
  with `Active branch: main (PR #65 merged via 799934d on
  2026-04-30; v0.14.2 tag pending)`.
- §"Current state" first paragraph: drop the "is being tagged
  together with the post-merge cleanup in the next release" clause;
  state plainly that PR #65 has merged and the v0.14.2 tag is
  outstanding.
- §"Next task": replace the entire current block (the "Tag v0.14.2"
  paragraph plus the two "After that" bullets) with a single
  paragraph naming **M15: GPU linear solver path** as the next task,
  with a one-line summary and a pointer to
  `docs/IMPROVEMENT_GUIDE.md` §4 M15 and `docs/M15_STARTER_PROMPT.md`.
- Add a new section `## Backlog` immediately after `## Next task` if
  one does not already exist. Move the validation-Phase-2 bullet
  (text starting "Physics validation suite, Phase 2 — external
  validation against Sze and Nicollian-Brews...") into it verbatim.
  Do not delete any factual content; the audit-suite status sentences
  inside that bullet remain accurate.
- Do **not** edit any entry in §"Completed work log". It is
  append-only.

### 2. docs/IMPROVEMENT_GUIDE.md

§1 "Honest current state":

- Replace `## 1. Honest current state (as of v0.8.0)` with
  `## 1. Honest current state (as of v0.14.1)`.
- Refresh "What exists and works" to reflect M9 through M14.2.
  Authoritative source is the §Roadmap table in `PLAN.md` and the
  §"Completed work log" entries for M14.2, M14.1, M14, M13.1, M13,
  M12, M11, M10, M9. Do not invent claims; lift wording from PLAN.md.
- Refresh "What does not exist" by removing the items that have
  shipped (artifact serialization, server surface, schema versioning,
  transient solver, AC small-signal, axisymmetric coordinate
  systems). Remaining honest gaps after refresh: GPU (M15), physics
  completeness (M16), heterojunctions (M17), Cartesian-2D MOSCAP
  variant and gate-driven HF C-V (M14.2.x in `PLAN.md`).

§4 "Milestones":

- For each of M9, M10, M11, M12, M13, M13.1 (if listed), M14, M14.1,
  M14.2: at the top of its §4 block, add a one-line bold status
  badge `**Status: Done (YYYY-MM-DD).**` followed by a single-line
  summary. Replace the deliverable / acceptance / scope / dependency
  body with the line `Full deliverable, acceptance tests, and scope
  preserved in §10 (Shipped milestone detail).` Then physically move
  the original body verbatim to a new appendix §10. Acceptance-test
  text moves with the body; it remains useful as a regression
  reference.
- Leave M15, M16, M17, M18 sections fully intact in §4.
- Note: §4 currently lacks an explicit M13.1 sub-section. If absent,
  insert a one-line stub `### M13.1 — Slotboom transient close-out
  **Status: Done (2026-04-27).** v0.14.1 closed out via Slotboom
  primary unknowns (ADR 0014); detail in §10.` under the M13 stub.

Create `## 10. Shipped milestone detail` at end of file (after §9
change log). Each subsection mirrors the §4 heading
(e.g. `### M9 — Result artifact writer + manifest`) and contains the
verbatim moved body.

§6 "UI integration checklist":

- Tick (`[x]`) the M9, M10, M11, `/materials`, `/schema` boxes.
- Leave `/capabilities`, WebSocket, CORS, OpenAPI, end-to-end
  integration boxes unchanged. Append to the `/capabilities` line
  the parenthetical `(M15 will close this; GPU availability is a
  capability)`.

§9 "Change log for this document":

- Append `- **2026-04-30** — Refresh of §1, §4, §6 to reflect v0.14.x
  reality post-M14.2; introduce §10 shipped-milestone appendix.`

### 3. README.md

- L397-L398: replace `The Day-1 math helper \`tests/check_day1_math.py\`
  covers thermal voltage,` with `The analytical-math helper
  \`tests/check_analytical_math.py\` covers thermal voltage,`. Keep
  the rest of the sentence intact.
- L414: replace `python tests/check_day1_math.py    # runs offline,
  no dolfinx required` with
  `python tests/check_analytical_math.py    # runs offline, no
  dolfinx required`.

### 4. Rename test helper

```bash
git mv tests/check_day1_math.py tests/check_analytical_math.py
```

Update remaining filename references:

- `CONTRIBUTING.md` L40 (`python tests/check_day1_math.py`).
- `.github/workflows/ci.yml` L40 (`run: python tests/check_day1_math.py`).
- `tests/test_axisym_moscap_math.py` L4 (the comment
  `mirrors the style of \`tests/check_day1_math.py\``).

## Verification

```bash
ruff check semi/ tests/
pytest tests/ -v 2>&1 | tail -30
```

Expected: all 237+ pure-Python tests still pass. No code paths were
modified. The renamed helper still executes when invoked directly.

If anything red, do not commit; flag back to the reviewer.

## Commit

Single commit. Message:

```
docs: refresh IMPROVEMENT_GUIDE for v0.14.x; remove residual day-language; rename check_day1_math
```

Do **not** push. The owner reviews locally and pushes manually.

## Invariants checklist

- [ ] No em dashes introduced in any prose line touched
      (PLAN.md invariant 8). Existing em dashes elsewhere in the
      file are out of scope for this commit.
- [ ] No source files (`semi/`, `semisim_server/`) modified.
- [ ] `tests/check_analytical_math.py` runs to completion offline
      (no dolfinx) when invoked directly.
- [ ] `PLAN.md` "Completed work log" entries are unchanged.
- [ ] Only the single rename appears in `git status` for `tests/`.
