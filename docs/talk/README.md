# Talk: SemiSim Speaker Notes

This folder contains a complete set of speaker notes for a 25–30 minute
talk on SemiSim: what the code does, how it was built, and why it was
designed the way it was.

## Files

| File | Contents | Suggested time |
|---|---|---|
| [01-opening.md](01-opening.md) | Hook, problem statement, what a device simulator is | 2–3 min |
| [02-physics.md](02-physics.md) | The three governing equations in plain English | 5–6 min |
| [03-numerical-challenges.md](03-numerical-challenges.md) | Why naïve FEM fails; nondimensionalization, Slotboom, continuation | 5–6 min |
| [04-architecture.md](04-architecture.md) | Five-layer design; JSON contract; HTTP API | 4–5 min |
| [05-capabilities-and-demos.md](05-capabilities-and-demos.md) | Live capabilities tour: devices, solvers, meshes, GPU path | 5–6 min |
| [06-development-process.md](06-development-process.md) | Milestones, ADRs, V&V strategy, CI, coverage gate | 4–5 min |
| [07-closing.md](07-closing.md) | Roadmap, honest gaps, takeaways, Q&A seed questions | 2–3 min |
| [08-comsol-comparison.md](08-comsol-comparison.md) | Comprehensive comparison with COMSOL Semiconductor Module | 8–12 min |
| [glossary.md](glossary.md) | Glossary of all physics, numerical, device, and software terms/acronyms | reference |

## Total target: 27–34 min (trim or expand §5 to hit your window)

For a talk that includes the COMSOL comparison, add 8–12 min or replace
§7 with a condensed version and use §8 as the closing discussion.
The glossary is intended as a leave-behind or pre-read, not a spoken slide.

## How to use these notes

Each file is structured as:

```
## Slide: <suggested slide title>
**[time]** Transition cue or delivery note.

<paragraph-form speaker notes>

**Key points** (bullet list for easy scanning at the podium)

**Transition:** bridging sentence to the next slide.
```

The notes are written to be spoken, not read verbatim. Treat them as a
script draft: rehearse once, then speak from the key-points bullets.
