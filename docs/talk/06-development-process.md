# §6 — Development Process

**Suggested slide title:** "How It Was Built: Milestones, ADRs, and a V&V Gate"
**Target time:** 4–5 minutes

---

## Slide 6.1 — Milestone-by-Milestone Delivery

**[28:00]**

The project was built incrementally in 20 milestones, each with a single
acceptance test that gated the merge. The early milestones built the
physics foundation; later milestones added capabilities on top of a
stable base.

| Milestone | What shipped | Gate |
|---|---|---|
| M1–M2 | Equilibrium Poisson 1D; multi-region | Depletion approximation <1% |
| M3 | Adaptive continuation; bias sweep | Shockley diode within 10% |
| M4 | V&V suite: MMS, mesh convergence, conservation | 62/62 PASS |
| M5–M7 | 2D, 3D, MOS capacitor, MOSFET, unstructured gmsh | Per-device verifiers |
| M9–M10 | Artifact writer; HTTP API | CLI + integration tests |
| M11–M12 | Schema versioning; XDMF mesh ingest | Strict schema v2.0.0 |
| M13 | Transient BDF1/BDF2 | Steady-state-limit + BDF rate |
| M14 | AC small-signal; MOSCAP axisymmetric | RC within 0.4%; LF/HF C-V |
| M15 | GPU linear-solver path | ≥5× speedup at 500k DOFs |
| M16+ | Physics completeness (planned) | Per-physics analytical gate |

The key discipline: never merge without a passing quantitative gate. The
coverage gate is 95% on every push; CI rejects the PR if it falls below
that.

**Key points**
- 20 milestones shipped; each has an acceptance test in `docs/IMPROVEMENT_GUIDE.md`.
- Order is driven by physical dependencies: you can't do AC without bias sweep,
  can't do GPU without a working CPU path.
- Coverage gate: 95% minimum, enforced by CI on every push.

---

## Slide 6.2 — Architecture Decision Records

**[29:30]**

When we made a choice that would affect the whole codebase — something a
future maintainer would need to understand — we wrote an ADR. An ADR
captures the context, the decision, and the consequences in 200–300
words. Once accepted, it is immutable. To change it you write a new one.

The most important ADRs for SemiSim:

- **ADR 0001** — JSON as input format. Why not Python DSL, TOML, or YAML.
  Short answer: schema validation at the boundary, round-trip through
  any tool, AI compatibility.

- **ADR 0002** — Nondimensionalization is mandatory. The 10³⁰-condition
  number argument. The asymmetric scaling choice (fields scaled,
  coordinates in meters) and the coefficient bug it caused in M1.

- **ADR 0004** — Slotboom variables for drift-diffusion. Why not SUPG
  or Scharfetter-Gummel. Short answer: Galerkin on pure-gradient form
  is stable and gives second-order convergence with no tuning.

- **ADR 0006** — V&V strategy. Why MMS over code-to-code comparison.
  Short answer: MMS catches sign bugs and coefficient bugs that
  single-point comparison against an analytical formula does not.

ADRs live in `docs/adr/`. The rule is: any change to an item listed under
"Invariants" in `PLAN.md` requires an ADR first.

**Key points**
- ADRs document "why this choice, not the alternatives" in a permanent record.
- Locked once accepted; overriding requires a new ADR, not an edit to the old one.
- The five invariants in PLAN.md (five-layer architecture, JSON contract,
  nondimensionalization, no PETSc types in the server API, no dolfinx in Layer 3)
  each have an ADR backing them.

---

## Slide 6.3 — Verification and Validation Strategy

**[31:00]**

The V&V suite (M4) is the engineering work I am most proud of in this
project. Here is why it matters and how it is different from "running
a benchmark."

A **benchmark** asks: does the code reproduce this analytical formula?
That is validation (are we solving the right equations?). It does not
tell you whether the FEM discretization converges at the correct rate.

**Method of Manufactured Solutions** (MMS) asks: if we inject a known
exact solution as a forcing term, does the error decrease at the
theoretical rate as we refine the mesh? MMS catches:
- Sign errors in a residual (the convergence rate collapses)
- Coefficient errors (the error magnitude is wrong at all resolutions)
- Stabilization bugs (order degrades on high-Péclet elements)
- Coupling errors between the three PDE blocks

The MMS suite in SemiSim covers: Poisson 1D and 2D, multi-region
Poisson 2D, and coupled drift-diffusion in three variants (Poisson-only,
full three-block with negligible recombination, and full three-block with
realistic SRH). The finest-pair L² rates land at theoretical 2.000.

Beyond MMS, **conservation checks** at each bias step verify that the
total current J_n + J_p is constant in space (as it must be from the
continuity equations) and that total charge integrates to zero at
equilibrium. These are cheap, exact checks that catch whole classes
of residual bugs without a mesh sweep.

**Key points**
- MMS vs. benchmark: MMS verifies convergence order; a benchmark verifies physics.
- Three MMS variants exercise progressively more of the DD residual.
- Conservation checks are cheap, exact, and run at every bias step in CI.
- Result: 62/62 V&V PASS; CI gate; any convergence regression fails immediately.

**Transition:** Let me close with a look at where the project is going
and what I'd do differently if I started over.
