# Architecture Decision Records

This directory holds the Architecture Decision Records (ADRs) for
kronos-semi. An ADR captures a single significant decision: its context,
the decision itself, and its consequences. Once accepted, an ADR is
immutable. To change an accepted decision, write a new ADR that
supersedes the old one and mark the old one as Superseded.

## When to add an ADR

Add an ADR for any change that:

- Modifies an item listed under "Invariants" in `PLAN.md`.
- Introduces a new external dependency or removes one.
- Changes the target version of a core tool (dolfinx, PETSc, Python).
- Alters the input format or breaks backward compatibility of the JSON
  schema.
- Picks between multiple defensible technical approaches where the
  choice is visible to users or future maintainers.

Do **not** add ADRs for routine code changes, bug fixes, or internal
refactors that do not cross an invariant.

## Format

Each ADR is a separate markdown file named
`NNNN-short-slug.md`, numbered sequentially from 0001. Each file follows
this skeleton:

```markdown
# NNNN. Short title

- Status: Proposed | Accepted | Deprecated | Superseded by ADR-XXXX
- Date: YYYY-MM-DD

## Context

What problem are we solving? What forces are at play?

## Decision

What did we decide to do? State it plainly.

## Consequences

What becomes easier, harder, or impossible as a result? Note any
follow-up ADRs that might be needed.
```

## Status values

- **Proposed:** under discussion; not binding.
- **Accepted:** in effect; code must conform.
- **Deprecated:** no longer applies but kept for historical context; no
  replacement needed.
- **Superseded by ADR-XXXX:** replaced; see the linked ADR for the
  current decision.

## Reference

Template and background: <https://github.com/joelparkerhenderson/architecture-decision-record>.

## Index

- [0001 JSON as input format](0001-json-as-input-format.md)
- [0002 Nondimensionalization is mandatory](0002-nondimensionalization-mandatory.md)
- [0003 Target dolfinx 0.10 API only](0003-dolfinx-0-10-api.md)
- [0004 Slotboom variables for drift-diffusion](0004-slotboom-variables-for-dd.md)
- [0005 Docker for dev, conda and FEM-on-Colab for users](0005-docker-for-dev-conda-fem-on-colab-for-users.md)
- [0006 Verification & Validation strategy](0006-verification-and-validation-strategy.md)
- [0007 Contact BC interface and `semi/bcs.py` layer](0007-contact-bc-interface.md)
- [0008 SNES tolerances](0008-snes-tolerances.md)
- [0009 Transient formulation](0009-transient-formulation.md)
- [0010 BDF time integration](0010-bdf-time-integration.md)
- [0011 AC small-signal](0011-ac-small-signal.md)
- [0012 Scharfetter-Gummel](0012-scharfetter-gummel.md)
- [0014 Slotboom transient](0014-slotboom-transient.md)
- [0015 Schottky as Robin BC](0015-schottky-robin-bc.md)
- [0016 Heterojunction-aware Slotboom and ohmic equilibrium](0016-heterojunction-slotboom.md)
