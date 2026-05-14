# 0001. JSON as input format

- Status: Accepted
- Date: 2026-04-20

## Context

SemiSim needs a way for users to specify a device: mesh, regions,
doping profiles, contacts, physics options, solver options, and output
preferences. Candidate formats are:

- JSON
- TOML
- YAML
- A Python DSL (the user writes `Device(...)` calls in Python)

The audience is mixed: TCAD engineers who may not know Python well, AI
coding assistants that benefit from a schema-validatable format, and
automated test harnesses.

## Decision

Use JSON as the only supported input format. Validate inputs against a
jsonschema Draft-07 schema defined in `semi/schema.py`.

## Consequences

Easier:

- Schema validation catches malformed inputs at the boundary with
  precise error messages, before any FEM code runs.
- JSON round-trips cleanly through web APIs, UIs, JavaScript tooling,
  and test fixtures.
- AI agents and human reviewers can inspect a device spec statically
  without executing any code.
- Cases are reproducible: a JSON file plus a package version fully
  determines the simulation.

Harder:

- Deeply nested configurations are verbose compared to a Python DSL.
- No comments in JSON proper; we rely on descriptive `name` and
  `description` fields to document intent.
- Default values must be maintained in two places (the schema and any
  downstream code that reads the config); schema tests guard against
  drift.

Rejected alternatives:

- **TOML:** elegant for shallow configs, but the schema has nested
  arrays of objects (doping entries inside regions, facets-by-plane
  inside mesh) where TOML's array-of-tables syntax becomes awkward.
- **YAML:** whitespace-significant parsing is fragile, and YAML's
  implicit type coercion has known foot-guns (for example, `1e17`
  parsing as a string in some parsers, or `no` parsing as boolean).
- **Python DSL:** running arbitrary Python to define a case is a
  security and reproducibility problem (a case description should not
  require executing code), prevents static analysis by AI assistants,
  and cannot be produced by a JS frontend.
