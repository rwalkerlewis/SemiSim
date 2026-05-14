# Claude Code workflow for SemiSim

This repo uses a **two-agent pattern** for AI-assisted development:

- A **reviewer** agent reads state, classifies issues, and writes prompts.
  It does not edit files, commit, or push.
- A **worker** agent executes specific tasks from prompts. It does not
  make architectural decisions.

Both agents are configured in `.claude/agents/`. The pattern exists to
separate "what should we do and why" from "how do we execute it
mechanically," so that judgment-heavy work gets the attention it needs
and mechanical work gets done efficiently.

## Working pattern

Run two Claude Code sessions in separate terminals, ideally on two
git worktrees of this repo so state stays in sync via git:

```bash
# From repo root, one-time setup
git worktree add ../SemiSim-review main

# Terminal 1 — reviewer (on main or PR branch under review)
cd ../SemiSim-review
claude
/agent reviewer
/model opus

# Terminal 2 — worker (on feature branch being built)
cd ../SemiSim
claude
/agent worker
/model sonnet
```

The reviewer inspects state and writes prompts to `prompts/<task>.md`.
The worker reads those prompts and executes them. When the worker pushes,
the reviewer fetches and audits.

## Model selection

- **Reviewer → Opus.** Architecture, classification, pushback, catching
  subtle issues. Judgment-heavy work.
- **Worker → Sonnet.** File edits, committed changes, following
  explicit instructions. Mechanical work.

## Ground rules both agents honor

1. **Read `PLAN.md` first.** The "Current state" and "Next task"
   sections are the source of truth for where the project is and what
   to work on.
2. **Read `docs/IMPROVEMENT_GUIDE.md`** before starting any M9+
   milestone. Acceptance tests per milestone are defined there.
3. **Respect the invariants in `PLAN.md`.** The five-layer architecture,
   the JSON-as-contract rule, the "no PETSc types leak across the
   engine API" rule. Changes to any invariant require an ADR first.
4. **Do not lower CI gates.** 95% coverage, MMS rates, analytical
   verification thresholds. If they fail, fix the cause.
5. **When unsure, stop.** The worker does not improvise. The reviewer
   does not hand off ambiguous prompts.

## Directory conventions

- `.claude/agents/` — agent role definitions (committed).
- `.claude/CLAUDE.md` — extended context that Claude Code loads
  automatically. Currently unused; this file serves that role via the
  repo root.
- `prompts/` — generated prompts from reviewer to worker. Committed
  so handoffs are visible in git history. Safe to archive or purge
  periodically.
- `docs/session-log.md` — optional running log of what each session
  did. Useful for cross-session continuity.
