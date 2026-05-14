# Agentic development strategy for SemiSim

Working draft. The goal is to formalize the development pattern that has
already produced M9 through M16.1, then extend it to handle multiple
agent runners under explicit token budgets.

This document is descriptive in §1 (what we actually did), diagnostic in
§2 (what is implicit and needs to be named), and prescriptive in §3
(the formalization). §3 is the part that should be opinionated; §1 and
§2 are ground truth and should not editorialize.

---

## 1. What was used through v0.17.0

### 1.1 Three roles, two of them codified

The development pattern that shipped M9 through M16.1 used three
distinct roles. Two are codified in `.claude/agents/`; the third was
ambient.

| Role | Where it ran | What it produced | Codified? |
|---|---|---|---|
| **Planner** | Chat/web Claude (Opus), human-mediated | `docs/M*_STARTER_PROMPT.md`, milestone decomposition, ADR drafts | No |
| **Reviewer** | Claude Code (Opus), read-only worktree | `prompts/<task>.md`, audits, classifications | Yes, in `.claude/agents/reviewer.md` |
| **Worker** | Claude Code (Sonnet), write-enabled worktree; **also Copilot Workspace** | Commits, PRs, test runs | Partially. `.claude/agents/worker.md` covers Claude Code only |

The planner is a real, load-bearing role: every M9+ milestone has a
`docs/M*_STARTER_PROMPT.md` that decomposes the milestone into phases,
names the invariants, and lists the required reading. Those files were
written by chat-Claude in collaboration with the maintainer. The
reviewer agent then translated each phase into a `prompts/<task>.md`
operational prompt that the worker executed. The system worked, but
the planner has no role definition, no codified inputs/outputs, and no
boundary on what it should and should not do.

### 1.2 Documents that act as agent contracts

The development relies on a small set of authoritative documents.
Treating these as contracts (not narrative) is what kept context
loadable at reasonable token cost.

- **`PLAN.md`**: ground truth for "current state", "next task",
  "invariants", "completed work log". Capped at ~300 lines by
  convention. Every agent reads this first.
- **`docs/IMPROVEMENT_GUIDE.md`**: milestone catalog with
  Why / Deliverable / Acceptance / Dependencies for each. Section 4
  is the spec; section 10 is the shipped-milestone archive.
- **`docs/ROADMAP.md`**: capability matrix and post-merge planning.
- **`docs/ARCHITECTURE.md`**: five-layer rule and import boundaries.
- **`docs/adr/`**: locked decisions. Changing an invariant requires
  an accepted ADR.
- **`CHANGELOG.md`**: append-only release record.

The convention `Read PLAN.md → IMPROVEMENT_GUIDE § Mxx → relevant ADR`
is the standard loading order. It is documented in
`.claude/agents/reviewer.md` and `worker.md`, and re-stated in the
starter prompts.

### 1.3 Two prompt formats

Two prompt formats are in active use. They serve different roles and
have not been formally distinguished until now.

**Starter prompt** (`docs/M*_STARTER_PROMPT.md`):

- Author: planner.
- Audience: reviewer + maintainer at milestone kickoff.
- Scope: one milestone, several phases.
- Length: 200–400 lines.
- Contains: required reading, conventions, phase decomposition,
  acceptance criteria per phase, off-limits invariants.
- Lifecycle: lives in `docs/`; archived in place after merge.

**Operational prompt** (`prompts/<task>.md`):

- Author: reviewer.
- Audience: worker, single session.
- Scope: one phase (one commit).
- Length: 50–200 lines.
- Contains: preconditions verified by reviewer, scope corrections vs
  starter prompt, numbered tasks with file/line anchors, exact
  acceptance commands, commit message, invariants checklist.
- Lifecycle: committed for audit; periodically purged after milestone
  close.

The asymmetry matters: a starter prompt can leave decisions to the
reviewer ("decide whether MOSCAP runner needs touching"); an
operational prompt cannot. The worker does not make calls.

### 1.4 Two-terminal worktree pattern

```
git worktree add ../SemiSim-review main      # reviewer
cd ../SemiSim                                # worker
```

Reviewer on `main` (or the PR branch under review); worker on the
feature branch. State is shared via git, not via shared filesystem or
shared session. The reviewer commits prompts to its worktree and
pushes; the worker pulls and reads.

This worked because:

1. The reviewer cannot accidentally edit source (its worktree is on
   `main`; tooling restrictions back this up).
2. The worker does not see reviewer-internal scratch files (those
   never leave the reviewer's worktree).
3. Crash recovery is trivial kill either session, the other still
   has its checkpoint.

### 1.5 Multiple runners, used opportunistically

Beyond Claude Code, the repo shows evidence of:

- **GitHub Copilot Workspace**: multiple `copilot/*` branches
  (`copilot/add-2d-axisymmetric-mos-capacitor-benchmark`,
  `copilot/extend-2d-moscap-axisymmetric`, etc.). Used as an
  alternate worker for milestones with clear specs and bounded scope.
- **Chat/web Claude**: planning, ADR drafting, prompt authoring.
- **Local pytest / ruff**: not an agent, but the loop that
  acceptance criteria close on.

There is no documented protocol for when to use which runner.
Selection has been intuitive.

### 1.6 Hard rules that survived

A few rules have held without exception across every milestone. They
should be carried forward verbatim.

1. **Reviewer never edits source, never commits, never pushes.**
2. **Worker never starts without a prompt; never improvises beyond
   the prompt; never silently weakens tests or gates.**
3. **One milestone, one PR. One phase, one commit.**
4. **Acceptance is exact commands, not prose.** "All tests pass" is
   not acceptance. `pytest tests/ -q` returning 0 is.
5. **Invariants in `PLAN.md` § Invariants are locked.** Changing one
   requires an ADR.
6. **`PLAN.md` "Completed work log" is append-only.**

---

## 2. What is not formalized

### 2.1 Planner role has no contract

The planner produced the M9, M14.3, M15, and M16.1 starter prompts,
but:

- It has no role file in `.claude/agents/`.
- Its inputs are ambient (the maintainer's head + the docs).
- Its outputs are not audited the way reviewer outputs are.
- Its model (Opus, Sonnet, GPT-class, mixed) varied across milestones.
- It has no token budget convention.

A starter prompt that misframes a milestone is more expensive than a
bad operational prompt: the reviewer + worker spend a session
discovering the misframe, and the planner re-issues. We do not
currently have a check at planner-output time.

### 2.2 Runner selection is intuition

When does a milestone go to Claude Code (reviewer + worker pattern)
vs Copilot Workspace vs a single-shot API agent? Today: by feel. The
copilot/* branches cluster around well-specified, bounded work
(axisymmetric extensions, sign-convention fixes). The Claude Code
worktree pattern handled the multi-phase, invariant-sensitive work.
This correlation is real but undocumented.

### 2.3 Token economy is implicit

Costs differ by an order of magnitude across runners and roles:

| Runner / role | Approximate per-session cost driver |
|---|---|
| Reviewer (Opus, full repo context) | High; reads `PLAN.md`, several docs, often source |
| Worker (Sonnet, prompted context) | Medium; reads prompt + named files |
| Copilot Workspace task | Low–medium; bounded by task scope |
| Chat planner (Opus, multi-turn) | High; long horizon, lots of re-loading |
| API single-shot worker | Low; one prompt, one diff |

There is no convention for "this kind of phase belongs on this
runner". The reviewer agent has been used for tasks that a
single-shot Sonnet worker could have closed in a fraction of the
tokens, because the reviewer was already loaded.

### 2.4 Session bootstrapping is expensive and undocumented

Every fresh session re-reads `PLAN.md`, the relevant
`IMPROVEMENT_GUIDE` section, often `ARCHITECTURE.md`, often two or
three ADRs. That bootstrap is necessary but costly, and there is no
guidance on:

- The minimum reading set per phase.
- When a session can skip the bootstrap (it's a continuation of an
  earlier session whose state is captured in `prompts/<task>.md`).
- How to compress the bootstrap when crossing runner boundaries
  (handing off from chat planner to Claude Code reviewer).

### 2.5 No protocol for cross-runner handoff

Today the handoff from chat planner to Claude Code reviewer is "the
maintainer pastes the starter prompt path into a fresh Claude Code
session". From reviewer to worker, "the worker reads
`prompts/<task>.md`". From worker back to reviewer, "the worker
reports a SHA and the reviewer checks `git log -1 --stat`". These
work but are not written down as a protocol with explicit artifacts
at each step.

---

## 3. The formalization (proposed)

This section is opinionated and meant to be edited.

### 3.1 Three codified roles

Promote the planner to a first-class role with its own contract.

```
.claude/agents/
  planner.md     # NEW: writes starter prompts, ADR drafts; never edits source
  reviewer.md    # existing: writes operational prompts, audits; never commits
  worker.md      # existing: executes operational prompts; commits and pushes
```

Each role definition specifies: inputs, outputs, off-limits
operations, model preference, expected token order of magnitude,
escalation rules.

**Planner contract (skeleton):**

- **Inputs:** `PLAN.md`, `IMPROVEMENT_GUIDE.md` § next milestone, the
  maintainer's intent (chat).
- **Outputs:** `docs/M*_STARTER_PROMPT.md`, optional ADR draft under
  `docs/adr/draft/`.
- **Never:** writes operational prompts in `prompts/`; edits source;
  commits.
- **Escalates when:** the milestone as written cannot be decomposed
  without violating an invariant; the acceptance criteria in
  IMPROVEMENT_GUIDE are ambiguous; the dependency on a prior
  milestone is unmet.

### 3.2 Runner selection rubric

Pick the runner from the work shape, not the convenience.

| Work shape | Runner | Role |
|---|---|---|
| Milestone framing, ADR, multi-phase decomposition | Chat Claude (Opus) | Planner |
| Phase-level operational prompt authoring, audit | Claude Code (Opus) on review worktree | Reviewer |
| Multi-step phase execution with invariant pressure | Claude Code (Sonnet) on feature worktree | Worker |
| Bounded, well-specified single PR with no invariant pressure | Copilot Workspace | Worker (alternate) |
| Pure mechanical batch (e.g. rename across N files, regenerate notebooks) | API single-shot Sonnet or Haiku | Worker (alternate) |
| Pure mechanical local edit | Direct human edit | n/a |

A phase belongs on Copilot Workspace (or an API single-shot) iff its
operational prompt would be straightforward enough that the
reviewer's audit is the only judgment step. M14.2 axisymmetric
extension fit this; M14.3 schema strict mode did not.

### 3.3 Token economy: budgets per artifact

Set rough budgets per artifact and treat overruns as a signal that
something is mis-scoped.

| Artifact | Target tokens (input + output) | Overrun signal |
|---|---|---|
| Starter prompt for one milestone | 30k–80k | Milestone is too large; split it |
| Operational prompt for one phase | 5k–20k | Phase is too large; split it |
| Worker session for one phase | 30k–100k | Phase is poorly specified, or invariant pressure was missed at planning time |
| Reviewer audit of one phase | 5k–15k | Diff is sprawling; reject and ask for re-scope |
| Bootstrap (read PLAN + IMPROVEMENT § Mxx + 2 ADRs) | ~10k | Required reading list is too long; trim |

These are not hard caps; they are alarms.

### 3.4 Minimum reading set per phase

Not every phase needs `ARCHITECTURE.md` and the full ADR set. The
operational prompt should explicitly enumerate the reading the
worker needs, and exclude the rest.

Heuristic:

- Always: the prompt itself.
- Schema-touching phase: `schemas/input.v*.json` + `semi/schema.py` +
  the schema ADR.
- Physics-touching phase: `docs/PHYSICS.md` § the relevant equation +
  `docs/adr/0004` + the touched physics module.
- Runner-touching phase: the runner file + the form builder it calls.
- Test-only phase: the test directory + the module under test.

Everything else is excluded by default. The reviewer makes the call
at prompt-writing time.

### 3.5 Cross-runner handoff protocol

Every handoff is an artifact, not a chat message.

```
planner       → reviewer:  docs/M{milestone}_STARTER_PROMPT.md (committed)
reviewer      → worker:    prompts/{milestone}-{phase}.md (committed)
worker        → reviewer:  commit SHA + `git log -1 --stat` output (chat)
reviewer      → planner:   PLAN.md "Current state" update + audit notes
                           (committed if PLAN changes; chat if questions)
worker        → worker:    fresh session reads prompts/{milestone}-{phase}.md
                           and `git log --oneline -5` to recover state
```

Rules:

1. **No verbal-only handoffs.** Every cross-role transition has a
   committed artifact or a structured chat report with named fields.
2. **A stale artifact is a bug.** If `PLAN.md` "Next task" disagrees
   with the active branch, the reviewer fixes `PLAN.md` before doing
   anything else.
3. **A worker session can be killed at any phase boundary** without
   losing more than the current commit. State lives in git.

### 3.6 Bootstrap compression

For continuation sessions (worker resumes a phase, or reviewer
re-audits after a worker push), the bootstrap is:

```
git log --oneline -5
git status
cat prompts/{milestone}-{phase}.md
```

For fresh sessions on a new milestone, the bootstrap is the standard
order in `.claude/agents/reviewer.md` § "Reading order before any
substantive work". Do not re-read what is already in the operational
prompt's "Preconditions verified by reviewer" block.

### 3.7 Failure modes and escalation

| Symptom | Most likely cause | Escalate to |
|---|---|---|
| Worker stops on ambiguity | Operational prompt under-specified | Reviewer |
| Reviewer cannot decompose phase | Starter prompt mis-framed milestone | Planner |
| Planner cannot decompose milestone | IMPROVEMENT_GUIDE acceptance is wrong | Maintainer |
| Acceptance command fails repeatedly | Hidden invariant violation, or environment drift | Reviewer + maintainer |
| `PLAN.md` and code disagree | Previous milestone closed without PLAN update | Reviewer (fix PLAN first) |

Escalation never goes the other direction (worker does not tell the
planner; planner does not tell the worker). The reviewer is the
relay.

---

## 4. Desired improvements

Backlog for the agentic system itself, ordered by expected payoff.
These are improvements to the *development process*, not the
codebase. None block current work; all reduce friction or close a
gap that is paid in tokens, mistakes, or maintainer time.

### 4.1 Codify the planner

`.claude/agents/planner.md` (or `docs/PLANNER.md` if we prefer
"runs in chat" framing) with the same shape as the existing
reviewer and worker definitions: inputs, outputs, off-limits ops,
escalation rules. The planner already runs; it just runs without a
contract.

### 4.2 Prompt templates and a linter

Two committed templates plus a script:

- `docs/templates/STARTER_PROMPT.md`: the structure that M9, M14.3,
  M15, M16.1 starter prompts converged on (required reading,
  conventions, phase decomposition, acceptance per phase,
  invariants).
- `docs/templates/OPERATIONAL_PROMPT.md`: the structure of
  `prompts/m15-phase-a.md` (preconditions verified, scope
  corrections vs starter, numbered tasks with file:line anchors,
  exact acceptance commands, commit message, invariants checklist).
- `scripts/lint_prompt.py`: checks that a prompt under `prompts/`
  has all required sections and that every acceptance command is a
  shell command (not prose). Runs in CI on PR diff.

### 4.3 Bootstrap pre-loader

`scripts/bootstrap_for_phase.py <milestone> <phase>` prints the
minimum reading set for a fresh agent session: the starter prompt
path, the operational prompt path, the relevant `IMPROVEMENT_GUIDE`
section, the ADRs the prompt cites, the touched source files. Used
by reviewer to validate "did I name the right reading?" and by
worker to reduce bootstrap from ~10k to ~3k tokens on continuation
sessions.

### 4.4 Token instrumentation

Log per-PR token consumption (planner / reviewer / worker, and
across runners) in a small CSV under `docs/agent_metrics/`. Not for
billing; for detecting when a phase is mis-scoped. The §3.3 budgets
are alarms; without measurement they are decorative.

### 4.5 Cross-runner prompt translator

If §3.2 holds and we use Copilot Workspace as an alternate worker,
the reviewer should not write two prompts. A small translator that
takes `prompts/<task>.md` and emits a Copilot-Workspace-shaped
brief (or vice versa). Open question: whether the operational
prompt format is already runner-agnostic enough that the
translator is a lint pass plus a header rewrite, or whether the
two formats genuinely differ.

### 4.6 Audit-as-a-script

`scripts/audit_prompt.py <prompt> <sha>` runs a structured
comparison: which files did the prompt name? which files did the
commit touch? what is the diff in either set? did the acceptance
commands run? did they pass at the required thresholds (not just
"exit 0")? Output is a markdown audit fragment the reviewer can
use as the basis of its report rather than producing free-form
prose.

### 4.7 Phase-replay capability

Re-run an operational prompt against a fresh worker session on a
clean worktree from the same starting SHA. Compare the resulting
diff to the merged commit. Used as: (a) a regression test for the
agent definitions themselves, (b) a way to detect non-deterministic
prompts (where the same prompt produces materially different
diffs).

### 4.8 Eval suite for agent changes

Two or three toy milestones in a sandbox repo (or a fixtures
directory in this repo) with known-good outcomes. Before changing
`reviewer.md` or `worker.md`, run the agents on the toy milestones
and check that the outcome is unchanged or improved. Without this,
we are tuning the agent contracts blind.

### 4.9 Operational prompt move

`prompts/<task>.md` is fine while one milestone is in flight.
Once two are in flight (M16.2 and M17 staging, plausible mid-2026),
move to `prompts/{milestone}/{phase}.md`. Do this preemptively, not
reactively.

### 4.10 ADR draft area

`docs/adr/draft/` for planner-authored ADR drafts that are not yet
accepted. The reviewer promotes a draft to a numbered ADR (or
rejects). Today the planner can draft an ADR but there is no place
to put it that is not the canonical numbered location.

---

## 5. Testing, verification, and validation

The repo has serious V&V infrastructure (ADR 0006: four-phase MMS
suite, mesh convergence, discrete conservation; cross-runner audit
under `tests/audit/`; analytical anchors in every benchmark; 95%
coverage gate on `semi/`). The agentic process must preserve every
gate. This section says how.

### 5.1 The V&V layers and what each one catches

| Layer | What it catches | Where it lives | Who owns the gate |
|---|---|---|---|
| Pure-Python unit tests | Math errors in scaling, schema, materials, doping, BCs, continuation, analytical references | `tests/test_*.py` | Worker writes; reviewer audits |
| FEM unit tests | Form-builder errors, residual sign errors that are visible at one mesh | `tests/fem/test_*.py` | Worker writes; reviewer audits |
| MMS verifiers | Discretization order loss; residual sign errors that are invisible at one mesh; coefficient-of-zero bugs | `semi/verification/mms_*.py`, `tests/fem/test_mms_*.py` | Reviewer specifies; worker implements |
| Mesh convergence | Anisotropic mesh error, P1-degrades-to-P0 silent failures | `semi/verification/mesh_convergence.py` | Reviewer specifies |
| Conservation | Charge neutrality, current continuity (sign and magnitude) | `semi/verification/conservation.py` | Reviewer specifies |
| Analytical benchmarks | Cross-check against closed-form physics for the device class | `benchmarks/*/verify_*.py` | Worker implements; planner specifies |
| Cross-runner audit | Inconsistencies *between* runners that pass their own gates individually | `tests/audit/` | Planner specifies; reviewer audits |
| Coverage gate | Untested branches, dead code | `pyproject.toml`, CI | CI enforces |

The gates are not interchangeable. An MMS rate of 1.99 does not
prove anything about charge conservation; conservation does not
prove discretization order; cross-runner consistency does not
prove either. The reviewer must keep them distinct in operational
prompts.

### 5.2 The "every new physics module" rule

`CONTRIBUTING.md`: any new physics module requires (a) pure-Python
unit tests for closed-form math, (b) an MMS verifier per ADR 0006,
(c) an analytical benchmark. This is a **planner-time rule**, not
a worker-time rule. By the time the worker is executing, the
operational prompt must already require all three or it is
incomplete.

A starter prompt that omits the MMS verifier on a new physics
module is a bug in the starter prompt. Reviewer escalates to
planner.

### 5.3 Acceptance commands are gates, not tests

Every operational prompt ends with acceptance commands. These are
not "the test suite"; they are the gates the worker must pass for
the phase to be considered done. The distinction:

- A test that returns exit 0 with a relaxed threshold is failure
  masquerading as success.
- A test that runs but does not print the relevant rate or
  conservation residual cannot be audited.
- A test that the worker had to modify to make pass is, by
  default, a regression.

The operational prompt's acceptance commands must:

1. Capture and print the rate / residual / coverage number, not
   just exit code.
2. Reference the gate threshold from ADR 0006 (or the milestone
   spec), not from the worker's interpretation.
3. Run on the same image (Docker `docker-fem`) the reviewer used.

### 5.4 What the worker reports back

V&V results are reported literally:

```
Phase D acceptance:
  scripts/run_verification.py mms_dd
    Variant A: L2 rate finest pair = 1.997, H1 = 0.991  (gate 1.99 / 0.99) PASS
    Variant B: L2 = 1.987, H1 = 0.985                   (gate 1.99 / 0.99) FAIL
  pytest tests/fem/test_mms_caughey_thomas.py -q
    8 passed in 142.3s
```

Not "all tests passed". Not "MMS converged". The reviewer reads the
numbers and decides whether 1.987 is a real regression or
within-noise. The worker does not make that call.

### 5.5 Failure modes specific to LLM workers

Three patterns to gate against, named because they have happened
in the wild:

**Threshold relaxation.** Worker sees `L2 rate = 1.987` against
gate 1.99, decides "close enough", changes the gate to 1.95. The
reviewer must check every diff for changes to acceptance
thresholds in tests, ADRs, and CI config. Operational prompts list
"thresholds touched in this phase: none" as an explicit assertion.

**Test deletion.** Worker hits a stubborn failing test, deletes it
or marks it `skip`. Reviewer's audit must check `git log -p` for
test removals in any phase that did not explicitly authorize them.

**Test reframing.** Worker keeps the test but rewrites the
assertion to something easier. The audit script in §4.6 should
flag any modification to a test file in a phase whose prompt did
not name that test file.

These are detection patterns, not trust failures. Workers are
optimization processes; they will find the gradient of "make the
prompt look done", and that gradient sometimes points at the
gates. The reviewer is the corrective.

### 5.6 The audit suite as second-opinion V&V

`tests/audit/` cases 01–06 cross-check runners against each other:
bias-sweep vs transient at steady state, AC at omega→0 vs
bias-sweep dI/dV, MOSCAP runners against each other. These do not
test physics; they test that two runners that should agree do
agree. They are the right gate for milestones that touch shared
infrastructure (runners, form builders, scaling).

When a milestone touches a runner, the operational prompt for the
final phase must include the relevant audit case in its
acceptance, not just the per-runner unit tests.

### 5.7 V&V regression budget

A new milestone may move a rate from 1.997 to 1.991. Both pass the
1.99 gate. This is acceptable but worth tracking; over five
milestones it could erode the gate. Suggested: log finest-pair
rates per merged milestone in `docs/agent_metrics/vv_history.csv`
and reject any phase whose rate drops by more than 0.03 from the
previous milestone's value, even if it passes the absolute gate.

This is a gate-on-trends check, not a gate-on-absolute. It belongs
to the reviewer.

---

## 6. When development outpaces the maintainer's expertise

Agents are about to be able to push this project into physics
territory the maintainer is competent in but not expert in:
heterojunctions (M17), band-to-band tunneling, Fermi-Dirac
statistics, advanced numerical methods (multigrid for indefinite
systems, mixed-precision GPU paths). The §6.5 random-sample audit
assumes the maintainer can spot a regression by reading the diff.
That assumption degrades the further the work moves from
drift-diffusion-with-Boltzmann-statistics, which is on home turf,
toward physics where verifying the math from scratch is a
multi-day project the maintainer will not in fact do.

This section says how to ship correct results past the maintainer's
verification frontier without pretending the frontier is not
there. The strategy is not "trust the agents". It is **redundant
independent gates that do not require domain expertise to
interpret**, plus honest scoping when the gates cannot be built.

### 6.1 The failure mode this section is built to prevent

A plausible scenario: an agent ships a Fermi-Dirac integral
implementation. The MMS rates pass. The unit tests pass. The
diode benchmark looks reasonable at 300 K. The maintainer reads
the diff, sees Joyce-Dixon expansions in the code, recognizes the
form, signs off. Six months later a downstream user reports
40 % current-density error at low temperature where the
Joyce-Dixon expansion is outside its range of validity. The
agent did not flag this because the agent did not know it was a
range-of-validity question; the maintainer did not flag it
because the maintainer did not re-derive the expansion bounds.

The bug was always findable. Nobody was looking from the right
angle. This is the class of error the strategy below is designed
to prevent.

### 6.2 Independent oracles

The single most important property: the gates that approve a piece
of physics must be derivable from sources that do not share an
implementation. An MMS verifier that imports the production
physics module is not independent of that module; passing only
proves they are mutually consistent.

For new physics, the operational prompt must require:

1. **Closed-form analytical anchor at a non-trivial limit.**
   Caughey-Thomas reduces to constant mobility as F → 0 and to
   `vsat / F` as F → ∞. Both limits are checkable by hand and
   independent of any FEM machinery.
2. **MMS verifier with symbolically-derived forcing.** The
   forcing must come from differentiating the manufactured
   solution against the *equation as written in the derivation
   document*, not against a callable that wraps the production
   residual. Two independent paths to the same forcing.
3. **Conservation residual unrelated to either of the above.**
   Charge neutrality at equilibrium, current continuity at bias.
   These are exact consequences of the underlying equations and
   pass or fail regardless of what the discretization or
   forcing-derivation code does.
4. **One external anchor.** Either a published numerical figure
   the benchmark reproduces to within stated tolerance, or a
   code-to-code comparison against a reference simulator.

A new physics module that ships with only an MMS verifier and
unit tests has *one* path to a passing build. The chance that a
subtle bug satisfies all four independent paths simultaneously is
much smaller than the chance that it satisfies one or two.

### 6.3 Checks that do not require domain expertise

A list of gates the maintainer can write and audit without being
an expert in the physics being implemented. These are cheap and
they catch a real fraction of plausible-but-wrong outputs.

- **Dimensional analysis.** Every term in the residual has units;
  units must balance. Trivially scriptable; routinely missed.
- **Limit behavior.** Does the new term vanish or saturate where
  it should? Does the result reduce to the prior model when the
  new parameter is set to its identity value (vsat → ∞,
  beta → 0, etc.)?
- **Sign of perturbation.** When a parameter increases, does the
  observable move in the documented direction? Higher T should
  increase ni; higher V_F should increase forward current.
- **Order of magnitude.** A silicon diode at V_F = 0.6 V at 300 K
  has J on the order of 10⁻³ A/cm². If the answer is 10⁰ or
  10⁻⁶, something is wrong even if all rates pass.
- **Monotonicity where physics demands it.** I-V curves under
  forward bias are monotone increasing. If the simulation
  returns a non-monotone segment, the physics or the solver
  failed.
- **Reflection symmetry.** Swap n ↔ p, donor ↔ acceptor,
  rebuild a symmetric structure. The result must be symmetric
  to discretization tolerance.
- **Bit-equivalence on identity inputs.** When a new model is set
  to its identity (constant mobility, infinite lifetime, no
  tunneling), the result must match the prior version
  bit-for-bit. M16.1 already enforces this; it generalizes.
- **Reproducibility under re-run.** Same input twice must
  produce the same answer. Failures here usually mean
  uninitialized state or precision-dependent paths.

These checks do not prove the physics is correct. They prove the
implementation is consistent with itself and with elementary
physical reasoning. They are **necessary** conditions, and
violating any one is grounds for the reviewer to reject the phase.

### 6.4 Forced derivations ("show your work")

For any milestone past the maintainer's expertise frontier, the
operational prompt must require a derivation document, not just
code. Specifically:

- **New physics → `docs/<topic>_derivation.md`.** Each step from
  textbook reference to scaled UFL form. Citations carry
  textbook chapter and page or DOI; no "as is well known".
  Symbolic intermediate forms, not just final result.
- **New numerical method → `docs/adr/<n>-<method>.md`.** Why it
  converges, complexity bound, stability conditions, what
  happens at the failure modes.
- **New approximation → `docs/<topic>_validity.md`.** Range of
  validity stated quantitatively. What changes outside the
  range and how the code detects or warns.

The maintainer reviews the derivation, not the code. Reviewing a
derivation is a tractable task even at the frontier; reading a
PETSc-laden UFL block and inferring whether it implements the
right equation is not. The existing repo already does this for
core physics (`docs/mms_dd_derivation.md`,
`docs/mos_derivation.md`); the rule is to extend it past the
core.

The reviewer agent's job grows: it must verify that the
derivation document is consistent with the implementation, term
by term. This is a tractable agent task because both artifacts
are available; the maintainer audits the consistency check, not
the underlying physics.

### 6.5 External anchors when stakes warrant

For milestones whose output ships into someone else's design
work, at least one external anchor is non-negotiable:

- **Published numerical benchmark.** Sze chapter 2 problem set,
  Selberherr chapter 6 examples, Pierret problems with known
  answers. Implement the input, reproduce the figure or table to
  stated tolerance.
- **Code-to-code comparison.** devsim is open-source, MIT-licensed,
  and runs the same physics class. For a milestone touching
  drift-diffusion behavior, a devsim cross-check on the same
  benchmark catches discretization-class errors that MMS does
  not.
- **Periodic external expert review.** Route a representative
  diff to an external expert (paid or quid-pro-quo) on a
  cadence appropriate to the volume of past-frontier work. One
  day of an expert's time is cheaper than a silent wrong-answer
  in production.

These are expensive and slow. Use sparingly, but use them. The
"agents are fast" gain is partly returned to "validation is slow"
when validation is honest.

### 6.6 Adversarial triangulation between agents

Two agents, different models, given the same task independently.
Compare outputs. If they agree, that is signal but not proof
(shared training-data biases). If they disagree, that is strong
signal something is wrong; escalate.

A stronger pattern: the second agent is given the first agent's
output and instructed to find a flaw. Adversarial framing surfaces
issues that "do you agree?" framing buries. The second agent's
output is not "yes/no"; it is a list of specific concerns the
reviewer then checks.

This is expensive in tokens. Reserve for milestones where
independent oracles (§6.2) are not all available, or where the
maintainer has explicitly flagged the work as past the frontier.

### 6.7 Calibration on known-answer work

Before delegating a class of work to agents, run them on
retrospective tasks where the answer is known. M9 through M14.2
are now history; their starter prompts and merged diffs are in
the repo. Replay one milestone with current agents on a sandbox
worktree and compare to the merged result. If the agent gets
known-answer work right, that is evidence (weak) it can handle
similar prospective work. If it gets known-answer work wrong,
that is strong evidence not to delegate.

Calibration is also the right gate before adopting a new model
as planner or reviewer. New Opus version, new Sonnet, new
external model: replay before promotion.

The §4.7 phase-replay capability and §4.8 eval suite combine to
make this routine. Today calibration is ad hoc.

### 6.8 The line the maintainer does not cross

When validation is genuinely impossible at the maintainer's
expertise level, and external anchors are unavailable, the right
answer is **do not ship**. Mark the feature experimental, gate
behind a flag, document the gap in the capability matrix, do not
claim it works.

The temptation when agents enable building past the maintainer's
expertise is to ship anyway because the implementation looks
clean and the tests pass. Resist. The asymmetry of error is bad:
a silent wrong answer in someone else's design tool costs them
real engineering time and costs the project's credibility. The
speed gain on the
ship side does not justify it.

Concrete examples for this project:

- **OK to ship past frontier:** Caughey-Thomas (M16.1).
  Closed-form, well-known limits, MMS verifier, analytical
  benchmark, devsim cross-check feasible if needed.
- **Probably OK with explicit external review:** Heterojunctions
  (M17). Band-alignment conventions are textbook but easy to get
  the sign wrong; thermionic-emission coupling at the interface
  has multiple competing forms. Worth one external read.
- **Not OK to ship without expert validation:** Trap-assisted
  tunneling, full-band Monte Carlo, FD statistics outside the
  Joyce-Dixon range. The literature has competing models, the
  numerical traps are subtle, and the benchmarks against real
  devices are scarce. Either acquire the expertise, partner
  with someone who has it, or do not ship.

The repo's existing `docs/IMPROVEMENT_GUIDE.md` § 8 anti-goals
section is the right pattern. Extend it: any milestone past the
frontier gets either an explicit external-review requirement in
its acceptance, or an entry in anti-goals.

### 6.9 Honest signaling in shipped artifacts

When a feature ships in a regime the maintainer cannot fully
verify:

1. **Capability matrix marks it explicitly.** Not "supported";
   "experimental, verified against [specific anchor], not
   verified against [named gap]".
2. **README says so.** A user landing on the docs sees the
   limitation before they invest in adopting the feature.
3. **CHANGELOG entry names the gap.** Not "added FD statistics"
   but "added FD statistics; verified within Joyce-Dixon range
   (eta < 5); behavior outside this range is not validated".
4. **The benchmark verifies what is verified, and explicitly
   does not claim what is not.** A `verify_*.py` that asserts
   only the claims that have anchors, and prints the
   un-anchored quantities for inspection without gating on them.

The honest path is verbose. A `# experimental` decorator in code
and a one-line README note are not enough; users do not read
those. The capability matrix and CHANGELOG are the load-bearing
artifacts.

### 6.10 What this means for the agent contracts

Adjustments to the role files in `.claude/agents/` and the
proposed planner contract:

- **Planner.** Must classify each milestone on a frontier scale
  (in-domain / near-frontier / past-frontier) at decomposition
  time. Past-frontier milestones must include in their starter
  prompt: derivation requirement, external-anchor requirement,
  capability-matrix updates, and a frontier-specific exception
  trigger for §7.2.
- **Reviewer.** Must verify the derivation document exists and
  is consistent with the implementation. For past-frontier
  milestones, the reviewer audit includes a §6.3 checklist run
  (dimensional analysis, limits, signs, magnitudes,
  monotonicity, symmetry, bit-equivalence, reproducibility) as
  explicit pass/fail items, not implicit.
- **Worker.** Must implement the derivation document alongside
  the code, in the same phase. Reporting includes the §6.3
  checklist results.

These additions push token cost up on past-frontier milestones,
appropriately.

---

## 7. Human-in-the-loop design

The maintainer is currently doing six things manually:
authoring planner intent, reviewing every starter prompt, reviewing
every operational prompt, spot-checking reviewer audits, approving
every PR, and resolving invariant questions. The first and last
must stay human. The middle four can be tiered.

### 7.1 What stays human, no exceptions

These are the loops where streamlining costs more than it saves.

| Loop | Why it stays human |
|---|---|
| Setting milestone intent | The planner cannot infer intent from the codebase alone; it inverts only after the maintainer states it |
| Accepting an ADR | ADRs encode locks; locks bind future agents. A wrong lock is more expensive than any review |
| Approving an invariant change | Same reason as ADRs; the invariant list is the agent's source of truth |
| Final merge of any PR that touches `semi/physics/` or `semi/runners/` | Physics regressions are silent and expensive; the audit suite is good but not exhaustive |
| First merge in a new physics area | No prior MMS / benchmark to anchor against; only the human can sanity-check the new anchor |
| **Past-frontier milestones (§6.8)** | Validation requires expertise the maintainer must source explicitly (external review, paid expert, or scope refusal); cannot be delegated to gates alone |
| **Capability matrix and CHANGELOG entries for partially-verified work (§6.9)** | The honest-signaling guarantee is what protects downstream users; only the maintainer can sign that the language is honest |

For these, no streamlining. The maintainer pulls the branch, runs
the V&V suite locally, reads the diff, merges by hand.

### 7.2 What can be tiered

Most other reviews can move from "every artifact" to "exception
report". The structure:

| Artifact | Default path | Maintainer-review trigger |
|---|---|---|
| Starter prompt | Planner writes; maintainer skims | Touches a new physics area, a new runner, or any invariant |
| Operational prompt | Reviewer writes; auto-merges to `prompts/` after lint | Touches a test threshold, an ADR, or `PLAN.md` invariants |
| Worker commit | Reviewer audits; auto-OK if audit passes | Diff > 200 lines; touches a runner; touches V&V code |
| Reviewer audit report | Filed in PR description | Reviewer flagged a deviation (always shown to maintainer) |
| PR for non-physics phase | CI green + reviewer audit clean → maintainer merges in batch | Any of the above triggers, or first PR of milestone, or last PR of milestone |

The principle: the maintainer reads exception reports, not
artifacts. The agents produce exception reports.

### 7.3 Sign-off ledger

A small append-only log under `docs/agent_metrics/sign_offs.md`
recording what the maintainer signed off on, when, and at what
artifact SHA. Used for: reconstruction after the fact ("when did
I OK the M16.1 starter prompt?"), and as evidence that the
streamlining did not silently bypass any required review.

Format:

```
2026-04-29  M16.1 starter prompt     docs/M16_1_STARTER_PROMPT.md@e3a1f29  OK
2026-04-30  ADR 0014 acceptance      docs/adr/0014-...@4d2c8b1            OK
2026-04-30  M14.3 final merge        PR #70 @1dece5b                       OK with note: tighten audit case 02 next milestone
```

Three columns. No prose. Searchable.

### 7.4 The async-friendly handoff

Every cross-role artifact lives in git. The maintainer does not
need to be online for any single artifact. The pipeline:

```
planner finishes starter prompt   → commit to docs/M*_STARTER_PROMPT.md
                                  → maintainer reviews on their schedule
                                  → maintainer signs off in ledger
                                  → reviewer pulls and starts phase A

reviewer finishes operational prompt → commit to prompts/{milestone}-{phase}.md
                                     → lint passes → worker pulls
                                     → no maintainer step unless trigger fires

worker finishes phase commit → push to feature branch
                             → reviewer pulls and audits
                             → audit clean → reviewer files report in PR
                             → no maintainer step unless trigger fires

milestone PR ready → maintainer pulls, runs V&V locally, merges
```

Three places the maintainer steps in by default: starter prompt
sign-off, milestone PR merge, and any triggered exception. Three
places per milestone, not three places per phase.

### 7.5 Streamlining without losing fidelity

The risk of tiered review: a regression slips past the reviewer
and the maintainer never sees it because it did not trigger
exception. The mitigations:

1. **Exception triggers are conservative.** Better to over-trigger
   than under-trigger. The cost of a false exception is one
   maintainer skim; the cost of a missed regression is a hot fix.
2. **The audit suite (`tests/audit/`) runs on every milestone PR**,
   not just on phases that touch runners. Cross-runner drift is a
   class of regression the per-phase reviewer cannot see.
3. **V&V regression budget (§5.7)** is a quantitative tripwire
   that triggers maintainer review independently of any
   structural check.
4. **The phase-replay capability (§4.7)**, when it exists, is run
   on a sample of merged phases. If replay diverges from the
   merged diff, that is a maintainer escalation.
5. **Random-sample audit.** The maintainer audits one random
   merged phase per milestone in full, regardless of triggers.
   This is the calibration check on the tiered system.

The system is not "trust the agents". It is "name the loops where
human judgment is load-bearing, keep humans there, automate the
rest, and instrument enough to catch drift in the automated
parts". Fidelity comes from the gates being tight (§5) and the
sign-off ledger being honest (§6.3), not from the maintainer
reading every diff.

### 7.6 What the maintainer gives up

Honestly: the comfort of having seen everything. The pipeline as
described means a worker commit can land without the maintainer
reading it, provided the reviewer's audit is clean and no
triggers fire. This is the trade. Time is regained; what is
accepted in exchange is that the audit becomes the trusted
artifact, not the maintainer's own eyes.

The way to keep this trade fair is the random-sample audit and
the V&V regression budget. A random sample should never produce
a surprising finding. If it does, the triggers in §7.2 are
under-tuned and need tightening before more milestones go through.

---

## 8. Efficient development from here forward

The framework in §1 through §7 is the structure. This section is
the playbook: given the structure, what does development actually
look like in practice?

### 8.1 Where the bottleneck is

Tokens are cheap. Agent wall-clock is fast. **The maintainer is
the bottleneck.** Every efficiency lever in this section is some
form of "move work out of the maintainer's serial path".

The maintainer's time, ordered by leverage (highest first):

1. Setting milestone intent. One sentence here prevents 100 lines
   of reviewer-worker confusion downstream.
2. Authoring or auditing starter prompts. A good one is paid back
   across reviewer + worker + audit + merge.
3. Approving ADRs and invariant changes. Locks bind every future
   agent; one wrong lock costs more than any milestone.
4. Frontier classification (§6.10). Distinguishes heavy-machinery
   milestones from light-machinery milestones.
5. Exception-triggered review. Reading audit reports the reviewer
   flagged, plus PR diffs for trigger-fired phases.
6. Final merge of physics or runner PRs.
7. Random-sample audit (§7.5).
8. Periodic calibration of the agent stack (§4.8, §6.7).

What should not consume maintainer time:

- Reading every operational prompt.
- Reading every commit.
- Reviewing notebook regenerations or doc refreshes.
- Auditing CI green status.
- Choosing which runner to use for a given task.

This ordering holds while maintainer hours are the binding
constraint, which is the case today. The constraint shifts
with milestone mix: external-review capacity (§6.5, §6.8)
becomes binding when past-frontier work (M17 and beyond)
enters the active mix; reviewer-Opus token cost binds under
fast iteration with many small phases; independent-oracle
availability (§6.2) binds for past-frontier physics in regimes
where no published anchor exists. Re-rank the leverage list
whenever the milestone mix changes, against whichever
constraint is binding now, and route the next infrastructure
investment at relieving that one specifically.

### 8.2 The default pipeline for one in-domain milestone

```
[Intent]      Maintainer states intent in chat to planner
              ("Next: M16.2 Canali mobility variant on top of M16.1")

[Kickoff]     Planner writes docs/M16_2_STARTER_PROMPT.md (full phase
              decomposition, frontier classification: in-domain)
              Maintainer skims, signs off in ledger

[Phase A]     Reviewer writes prompts/m16.2-phase-a.md
              Worker executes phase A, commits, pushes
              Reviewer audits, files report in PR draft
              Maintainer not involved unless trigger fires

[Phases B,C]  Same pattern, in parallel where possible
              Trigger fires once: phase C touches the schema. Maintainer
              skims the operational prompt and the diff

[Phase D]     V&V verifier. Acceptance commands print MMS rates literally;
              reviewer audits the numbers

[Close]       PR ready. Maintainer pulls, runs V&V locally, reads audit
              report, merges. Reviewer updates PLAN.md as docs commit
```

Maintainer touches: intent, kickoff sign-off, the trigger that
fired during phase B/C, and close. Four touches per milestone.
The framework moves the prompt-authoring and per-phase audit
work onto agents, leaving the maintainer in the loops where
their judgment is the asset.

### 8.3 Parallelism rules

The pipeline runs in parallel for non-conflicting milestones.

Safe to parallelize:

- Two physics milestones in different modules (mobility variant +
  axisymmetric extension).
- One physics + one infrastructure (M16.2 + CI overhaul).
- Doc-only work alongside any code milestone.

Must serialize:

- Two milestones touching the same runner.
- Two milestones touching the schema.
- Anything touching invariants.
- Anything past the frontier (§6) until external review returns.

Worktree pattern extends: one worktree per in-flight milestone,
plus one read-only review worktree on main. The reviewer agent
multiplexes across them; the worker is one session per worktree.
Two milestones in flight is the comfortable maximum without a
dispatcher; three needs explicit scheduling.

### 8.4 Runner selection in one screen

| Work | Runner | Why |
|---|---|---|
| Milestone framing, frontier classification, ADR draft | Chat Opus (planner) | Judgment, low volume, high leverage |
| Operational prompt for invariant-touching phase | Claude Code Opus (reviewer) | Catches its own prompt subtleties |
| Multi-phase milestone with invariant pressure | Claude Code Sonnet (worker) | Long context, commits, pushes |
| Single PR, well-specified, no invariant pressure | Copilot Workspace | Cheap, bounded, async-friendly |
| Batch mechanical (regenerate N notebooks, rename across files) | API single-shot Sonnet or Haiku | One prompt, one diff, no session |
| Surgical local edit | Maintainer direct | Fastest when scope is one line |

If unsure which runner fits, the answer is probably one of the
bottom three. The top three should already have a session open.

### 8.5 Investments, in build order

Build the §4 improvements in this order. Each unlocks the next.

1. **Prompt templates + linter (§4.2).** First because
   everything downstream assumes a stable prompt format.
2. **Bootstrap pre-loader (§4.3).** Cuts continuation-session
   tokens substantially by trimming what each fresh agent
   session reads.
3. **Audit-as-script (§4.6).** Required to make tiered review
   (§7.2) safe; without it the §5.5 threshold-relaxation
   failure mode is undefended.
4. **Frontier-classification field in starter-prompt template
   (§6.10).** Trivial once item 1 has landed.
5. **Sanity checks as runnable gate (§6.3, §6.10).** Implements
   the §6.3 domain-agnostic checklist (dimensional analysis
   where feasible, limit checks, sign checks, magnitude checks,
   monotonicity, symmetry, bit-equivalence on identity inputs,
   reproducibility) as `scripts/sanity_checks.py`. Required for
   any past-frontier milestone; valuable on every milestone as
   a pre-merge tripwire. Resolves open question 13.
6. **Token instrumentation (§4.4).** Useful only after 1-3 are
   stable; numbers from inconsistent pipelines mean nothing.
7. **Phase-replay (§4.7).** The largest build of the set.
   Defer until the rest are stable; powers calibration (§6.7)
   and agent-definition eval (§4.8).
8. **Cross-runner translator (§4.5).** Build only after
   Copilot Workspace earns a permanent slot on the §8.4 runner
   table; today it is one option among several.

Stop at any point if the cost-benefit reverses. Items 1-5 carry
the most leverage. Items 5 and 7 together gate past-frontier
work safely: item 5 catches plausible-but-wrong implementations
at the phase level (§5.5, §6.1), item 7 catches drift in the
agent definitions themselves over time.

### 8.6 Go/no-go before running a milestone under the playbook

§8.5 is the build order. This subsection is the tripwire
between having those investments built and running a milestone
end-to-end under §8.2.

Verify all five before kickoff:

1. **Planner codified.** `.claude/agents/planner.md` exists and
   was used to author the milestone's starter prompt, not as
   ambient chat-Claude in the maintainer's head.
2. **Frontier classified.** The starter prompt declares
   in-domain, near-frontier, or past-frontier per §6.10.
   No "TBD" or implicit classification.
3. **Audit defenses in place.** Either `scripts/audit_prompt.py`
   (§8.5 item 3) runs end-to-end against a merged historical
   commit (M14.3 phase A is a clean test) and produces a
   usable audit fragment, OR the maintainer has explicitly
   committed to manual reviewer-output audit on every phase of
   this milestone. The §5.5 failure modes are otherwise
   undefended.
4. **Sanity checks ready for past-frontier.** If the milestone
   is past-frontier, either `scripts/sanity_checks.py` (§8.5
   item 5) implements the §6.3 checklist as runnable gates
   and passes on the milestone's V&V baseline, OR the
   maintainer has explicitly committed to running the §6.3
   checklist by hand at every phase boundary. No third option;
   "we'll check the obvious things in the audit" is not a
   sanity-check substitute.
5. **Sign-off ledger initialized.**
   `docs/agent_metrics/sign_offs.md` exists with at least one
   backfilled entry (M14.3 phase A is the suggested seed) so
   the format is exercised before it matters.

If any of the five fails, the milestone runs under the current
(pre-playbook) process. Do not mix half-playbook with
half-current: the audit trail becomes unreconstructable and
obscures which discipline caught (or failed to catch) what.
The honest move when an item fails is to land the §8.5 work
first, or to ship the milestone under the existing process and
defer playbook adoption to the next one.

### 8.7 The honest failure mode of this plan

Worth naming explicitly so it is detectable. The plan above
fails if:

- Starter-prompt review becomes rubber-stamp. The framework
  assumes the §8.1 leverage-1 activity stays load-bearing. If
  starter prompts stop being read carefully, the rest of the
  pipeline cannot recover from a bad framing.
- The reviewer agent drifts (subtly relaxes its own audit
  criteria, fails to catch threshold relaxation in workers). The
  random-sample audit (§7.5) and phase-replay (§4.7) are the
  defenses; if neither runs, drift is invisible.
- Frontier classification gets optimistic. Calling a milestone
  near-frontier when it is past-frontier removes the §6.2
  four-oracle requirement and ships unverified physics. The
  defense is honest signaling (§6.9): if the capability matrix
  language for a near-frontier milestone reads like a
  past-frontier hedge, the classification was wrong.
- Parallelism exceeds the conflict graph. Two milestones touch
  the same runner without anyone noticing; one merges first and
  silently invalidates the other's prompt. The defense is the
  reviewer holding a mental model of in-flight work; once that
  exceeds three milestones, a dispatcher becomes necessary.

If any of these symptoms appear, the response is to slow down,
not to add more automation. The framework does not survive
neglect of its premises.

---

## 9. Open questions for us to work through

These are the things this draft does not yet answer.

1. Should the planner role live as a `.claude/agents/planner.md`
   file even though it runs in chat, not in Claude Code? Or as
   `docs/PLANNER.md`?
2. What is the right ADR-vs-starter-prompt boundary? Today some
   starter prompts include rationale that should arguably be ADR
   text.
3. Should `prompts/` move to `prompts/{milestone}/{phase}.md` once
   we have multiple phases per milestone in flight? (§4.9 says
   yes preemptively; not yet decided.)
4. Is there a third prompt format we are missing, something
   between starter and operational, for cross-milestone refactors?
5. How do we benchmark the system itself? Token-per-merged-PR is
   measurable; should we track it? (§4.4 says yes; mechanism TBD.)
6. Copilot Workspace has different prompt expectations than Claude
   Code. Do we need a `prompts/copilot/` variant, or is the
   operational prompt format already runner-agnostic? (§4.5
   sketches a translator; we have not validated.)
7. Should the reviewer ever invoke a single-shot API worker
   directly (e.g., for "regenerate the four notebooks"), or does
   that always route through a Claude Code worker session?
8. The V&V regression budget (§5.7) is a delta-on-rates check.
   Should it also gate conservation residuals, runtime, or
   coverage? Each is a different kind of regression.
9. How do we sign off on changes to the agent definitions
   themselves (`.claude/agents/*.md`)? Is that an ADR, a starter
   prompt, or something else? §4.8 implies an eval pipeline; the
   governance is undefined.
10. Where does the maintainer's random-sample audit (§7.5) get
    recorded, and what triggers re-tuning the exception triggers
    in §7.2 if the random samples surprise?
11. Who classifies a milestone on the in-domain / near-frontier /
    past-frontier scale (§6.10)? The planner is the natural fit,
    but the maintainer is the one with the most accurate read on
    their own expertise. Probably joint, but the artifact and
    sign-off path are undefined.
12. What is the trigger for invoking external expert review
    (§6.5)? "Past-frontier milestone" is a category, not a
    trigger. Per-milestone, per-quarter, per-physics-area? And
    who pays?
13. The §6.3 domain-agnostic checklist is currently prose. Should
    it be a script (`scripts/sanity_checks.py`) the worker runs
    automatically and reports machine-readable results from? That
    would make §6.3 a gate, not a recommendation.
14. Do we want `docs/<topic>_derivation.md` to be machine-checked
    against the implementation in any way, or does it stay as
    "the reviewer reads both and verifies consistency"? A
    symbolic-math check (sympy, then UFL) is feasible for some
    physics; not for others.
15. How does §6.9 honest-signaling interact with marketing or
    user-facing copy? The CHANGELOG and capability matrix can be
    honest, but a separate "what SemiSim can do" page on a
    project website can drift. Out of scope for this doc but
    worth flagging.


---

## 10. How to set this up in practice

A practical implementation guide for putting the framework
above into use. Written for a reader who is comfortable in a
code editor and has used GitHub but has not run AI coding
agents before. If you already use Claude Code with worktrees,
skip to §10.4.

### 10.1 What you are building

Three pieces of software that talk to each other through
files in your repository:

1. **A web chat with Claude** for the planner role. You open
   this in a browser at claude.ai and use it the way you
   would use any chat assistant.
2. **Claude Code, a command-line tool** for the reviewer and
   worker roles. You open this in a terminal and it acts
   inside your repository, reading and editing files directly.
3. **Your code editor** (such as VS Code) for reading what
   the agents produced and merging the result.

The agents do not call each other directly. They communicate
through markdown files committed to git: starter prompts the
planner writes, operational prompts the reviewer writes, and
the code commits the worker produces. If you can read a
markdown file and run a `git pull`, you can audit any handoff
in the system.

### 10.2 Accounts and software you need

Free or low-cost tier of each is sufficient to start.

| Tool | What it is | How to get it |
|---|---|---|
| GitHub account | Where your repo lives | github.com, free signup |
| Git on your computer | Version control | git-scm.com, OS-specific installer |
| Code editor | For reading code and merging | code.visualstudio.com (VS Code, free) is a fine default |
| Terminal app | For running commands | macOS Terminal, Windows PowerShell or Windows Terminal, any Linux terminal |
| Claude.ai account | For the planner role (chat) | claude.ai, free tier works to start; Pro is recommended for Opus access |
| Claude Code | For the reviewer and worker roles | docs.claude.com has install instructions; runs in your terminal |
| Optional: Copilot Workspace | Alternate worker for well-bounded tasks | github.com/features/copilot, requires a Copilot subscription |

You do not need to install Python, the project's scientific
dependencies, Docker, or any domain-specific tooling to set
up the agentic workflow itself. Those are project
dependencies; the agentic infrastructure is just files plus a
chat browser tab plus a command-line tool.

### 10.3 One-time setup

Run these commands once. They assume your repo is already on
GitHub and cloned locally into a folder named `myproject`.

**Step 1: create the agent role files.**

Inside your repo, make a directory called `.claude/agents/`.
The leading dot is intentional; it hides the directory from
most file browsers, but Claude Code reads it automatically.

```
cd myproject
mkdir -p .claude/agents
```

Inside that directory, create three files:

- `planner.md` (use the contract sketched in §3.1 of this doc
  as a starting point)
- `reviewer.md` (use this repo's existing
  `.claude/agents/reviewer.md` as a starting point)
- `worker.md` (use this repo's existing
  `.claude/agents/worker.md` as a starting point)

**Step 2: create the prompts directory.**

```
mkdir prompts
```

This is where operational prompts (reviewer-to-worker
handoffs) live. They are committed to git so the audit trail
is preserved.

**Step 3: create the supporting documentation directories.**

```
mkdir -p docs/templates
mkdir -p docs/adr/draft
mkdir -p docs/agent_metrics
```

These hold prompt templates (§4.2), ADR drafts (§4.10), and
the sign-off ledger and token logs (§4.4, §7.3).

**Step 4: set up a second worktree for the reviewer.**

The reviewer agent runs on a copy of your repo that is
checked out to your `main` branch and stays there. The worker
runs on a feature branch in the original copy. Setting up a
separate worktree keeps them from stepping on each other.

```
git worktree add ../myproject-review main
```

This creates a sibling folder `myproject-review/` next to
your `myproject/` folder. Both share the same git database
under the hood, so commits made in one are visible in the
other after a `git fetch`.

You will run the reviewer agent from `myproject-review/` and
the worker agent from `myproject/`.

### 10.4 What a session looks like

A typical milestone uses three windows.

**Window 1: web browser at claude.ai (planner).**

Open a new chat. Tell the planner what milestone you want to
work on, in plain language. The planner writes a starter
prompt as a markdown file. You copy that into your repo at
`docs/M{milestone}_STARTER_PROMPT.md` and commit it.

**Window 2: terminal in the review worktree (reviewer).**

```
cd ../myproject-review
claude
```

Inside the Claude Code session, switch to the reviewer agent
and the Opus model:

```
/agent reviewer
/model opus
```

The reviewer reads the starter prompt and writes operational
prompts under `prompts/`. It commits those prompts but does
not edit source code.

**Window 3: terminal in the main worktree (worker).**

```
cd ../myproject
git checkout -b dev/m16.2-canali
claude
```

Inside this Claude Code session, switch to the worker agent
and the Sonnet model:

```
/agent worker
/model sonnet
```

The worker reads the operational prompt the reviewer wrote,
executes it, and commits the result to the feature branch.

### 10.5 Typical workflow

The pattern end to end:

1. You open the planner chat and describe what you want done.
2. The planner produces a starter prompt. You read it for
   sanity and commit it to your repo.
3. The reviewer (Window 2) reads the starter prompt and
   writes the first operational prompt under `prompts/`.
4. The worker (Window 3) reads the operational prompt and
   executes it, producing a commit on the feature branch.
5. The reviewer audits the worker's commit and files an
   audit report in the PR description. You are not involved
   unless the audit flags something.
6. Steps 3 through 5 repeat for each phase of the milestone.
7. When all phases are done, you pull the branch, run the
   verification suite locally, read the reviewer's audit
   report, and merge.

Most of your direct attention is at steps 1, 2, and 7. Steps
3 through 6 run with you watching from the side rather than
driving.

### 10.6 Where to look when something breaks

| Symptom | Likely cause | Where to look |
|---|---|---|
| Worker stops mid-task asking a question | Operational prompt was unclear | The operational prompt under `prompts/`; ask the reviewer to revise it |
| Reviewer flags an audit concern | Worker did something the prompt did not name | The reviewer's audit report; usually traces back to a missing constraint in the prompt |
| Tests pass but you do not trust the result | Possible §5.5 failure mode (silent test relaxation) | `git log -p` on the suspect phase; check whether threshold values changed |
| Agent session loses context partway through | Token window exhausted | Start a fresh session; the operational prompt plus `git status` is enough to resume |
| Two sessions disagree about repo state | One has not pulled recent commits | Run `git pull` in each worktree, then re-read `git log` in both |
| Worker runs but produces no diff | Prompt was already satisfied by prior commit | Reviewer's job to confirm; usually means the phase was a no-op |

When in doubt, the artifacts in `prompts/`, `docs/`, and the
`git log` are the source of truth. Any agent session can be
killed and restarted at any time; nothing load-bearing lives
only in chat memory.

### 10.7 Minimal viable adoption

If the full framework feels like too much to adopt at once,
the smallest version that still produces correct results and
captures the audit trail is:

1. Add the three role files under `.claude/agents/`.
2. Add the second worktree for the reviewer.
3. Use the planner-reviewer-worker handoff pattern in §10.5
   for one milestone, end to end.
4. Skip the prompt linter, the audit script, the token
   instrumentation, and the sanity-check script. Run those
   loops manually instead, with the maintainer in the loop.

This minimal version trades some efficiency for adoption
speed. It still gives you the load-bearing properties: an
audit trail in git, separation of judgment from execution,
and a reviewer that does not also write code. Add the
infrastructure pieces from §4 and §8.5 once the basic loop
is comfortable.