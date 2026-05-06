# JSON input schema reference

The input JSON contract is defined by two coexisting schemas
(both Draft-07): the strict
[`schemas/input.v2.json`](../../schemas/input.v2.json) (current)
and the legacy
[`schemas/input.v1.json`](../../schemas/input.v1.json) (deprecated,
kept for one minor cycle). Every object node carries a UI-facing
`description`. The engine loader is
[`semi/schema.py`](../../semi/schema.py); call
`semi.schema.get_schema(major)` for the schema for a specific major.

This page summarises the schema. For the authoritative definition,
the full annotated source of truth is the JSON file.

## Versioning

- `schema_version` is required at the top level. Format `MAJOR.MINOR.PATCH`.
- The engine accepts majors listed in `ENGINE_SUPPORTED_SCHEMA_MAJORS`
  in [`semi/schema.py`](../../semi/schema.py); currently `(1, 2)`.
- v1 inputs continue to validate against `input.v1.json` (loose: extra
  fields are silently ignored) but emit a `DeprecationWarning`. They
  will be removed after one minor release cycle.
- v2 inputs validate against `input.v2.json` (strict:
  `additionalProperties: false` on every object node, so typos in
  contact / mesh / solver fields fail validation rather than being
  silently dropped).
- Minor/patch skew is accepted silently within a major.
- Current schema version: **2.6.0** (M16.6 BBT and TAT tunneling
  dispatch); v2.0.0, v2.1.0, v2.2.0, v2.3.0, v2.4.0, and v2.5.0
  inputs continue to validate (additive minors).

Schemas are published to
`https://rwalkerlewis.github.io/kronos-semi/schemas/` on every release
tag; see [`.github/workflows/publish-schemas.yml`](../../.github/workflows/publish-schemas.yml).

### v1 -> v2 migration

For the M14.3 transition, every benchmark JSON in the tree is at
`schema_version: "2.0.0"`. To migrate a downstream input:

1. Bump `schema_version` from `"1.x.y"` to `"2.0.0"`.
2. Run validation. v2 will reject any field name that does not match a
   defined property; fix typos surfaced by the validator.
3. The data shape is otherwise unchanged; v2 is a pure-strictness bump.

History:
- **1.0.0** (M11) — extracted from `semi.schema.SCHEMA`.
- **1.1.0** (M13) — `solver.type = "transient"` and the BDF fields.
- **1.2.0** (M14) — `solver.type = "ac_sweep"` and the `dc_bias` / `ac` blocks.
- **1.3.0** (M14.2) — `coordinate_system` (`"cartesian"` / `"axisymmetric"`).
- **1.4.0** (M15) — `solver.backend` and `solver.compute`.
- **2.0.0** (M14.3) — strict mode (`additionalProperties: false` on every
  object node); v1 deprecated.
- **2.1.0** (M16.1) — added `physics.mobility.model: caughey_thomas` plus
  `vsat_n`, `vsat_p`, `beta_n`, `beta_p`. v2.0.0 inputs continue to validate.
- **2.2.0** (M16.2): additive Lombardi surface mobility. Adds
  `physics.mobility.model: lombardi`, the `bulk_model` selector
  (`constant` / `caughey_thomas`), the `lombardi` sub-object
  (`B_n`, `B_p`, `C_n`, `C_p`, `lambda_n`, `lambda_p`, `delta_n`,
  `delta_p`), and `interface_facet_tag` (declared `["integer", "null"]`
  in JSON Schema; the loader enforces non-null when
  `model == "lombardi"`). v2.0.0 and v2.1.0 inputs continue to validate.
- **2.3.0** (M16.3): additive Auger recombination kernel. Promotes
  `physics.recombination.auger` from a forward-compat placeholder to
  a real flag and adds `physics.recombination.C_n` (Si default
  2.8e-31 cm^6/s) and `physics.recombination.C_p` (Si default
  9.9e-32 cm^6/s; both Dziewior-Schmid). When `auger == false` (the
  default), the recombination kernel is bit-identical to v0.18.0.
  v2.0.0, v2.1.0, and v2.2.0 inputs continue to validate.
- **2.4.0** (M16.4): additive Fermi-Dirac statistics gated by
  `physics.statistics` dispatch. Widens the `physics.statistics`
  enum from `["boltzmann"]` to `["boltzmann", "fermi_dirac"]`. The
  `boltzmann` branch (default) is bit-identical to v0.19.0; the
  `fermi_dirac` branch substitutes the Blakemore approximation for
  F_{1/2} in the generalized-Slotboom helpers (ADR 0004) and the
  Poisson source. v2.0.0, v2.1.0, v2.2.0, and v2.3.0 inputs continue
  to validate.
- **2.5.0** (M16.5): additive Schottky contact type and
  `barrier_height_eV` (M16.5). Widens the `contacts[].type` enum from
  `["ohmic", "gate", "insulating"]` to add `"schottky"`, and adds
  `contacts[].barrier_height_eV` (declared `["number", "null"]` in
  JSON Schema; the loader enforces non-null and non-negative when
  `type == "schottky"`). The Schottky path applies a metal-Fermi-level
  Dirichlet on psi and a thermionic-emission Robin BC on the
  continuity rows; ohmic, gate, and insulating contacts remain
  bit-identical to v0.20.0. v2.0.0, v2.1.0, v2.2.0, v2.3.0, and
  v2.4.0 inputs continue to validate.
- **2.6.0** (M16.6): additive BBT and TAT tunneling dispatch
  (M16.6). Adds the `physics.tunneling` sub-object with `bbt` and
  `tat` boolean flags plus the Kane (`A_kane`, `B_kane`) and Hurkx
  (`tau_n_min`, `tau_p_min`, `F_kT`, `alpha`) parameters. Both flags
  default to false; the kernel is bit-identical to v0.21.0 when
  both are off. A UserWarning fires at validate time when `bbt`
  is true and `statistics` is `"boltzmann"` (BBT is most accurate
  under Fermi-Dirac at heavy doping). v2.0.0 through v2.5.0 inputs
  continue to validate.

## Top-level fields

| Field | Required | Type | Notes |
|---|---|---|---|
| `schema_version` | ✓ | string | Semver. Currently `"1.3.0"`. |
| `coordinate_system` | – | enum | `"cartesian"` (default) or `"axisymmetric"`. |
| `name` | ✓ | string | Case name; used in output paths. |
| `description` | – | string | Free-form. |
| `dimension` | ✓ | int | 1, 2, or 3. Must equal 2 when `coordinate_system="axisymmetric"`. |
| `mesh` | ✓ | object | Builtin box or external gmsh `.msh`. |
| `regions` | ✓ | object | name → material, role. |
| `doping` | ✓ | array | Per-region doping profiles. |
| `contacts` | ✓ | array | Boundary conditions. |
| `physics` | – | object | Recombination model parameters. |
| `solver` | – | object | Solver type and tolerances. |

## `mesh` block

Two variants chosen by `source`:

- `source: "builtin"` — structured box mesh.
  - `extents`: `[[xmin, xmax], ...]` per dimension, in meters.
  - `resolution`: cells per dimension.
  - `regions_by_box`: optional list of axis-aligned sub-boxes that
    tag cells by region (required for multi-region devices).
  - `facets_by_plane`: list of axis-aligned planes that tag boundary
    facets, referenced by name in `contacts[*].facet`.

- `source: "file"` — external gmsh `.msh` or XDMF.
  - `path`: relative to the JSON file's directory.
  - `format`: `"gmsh"` or `"xdmf"`. The XDMF branch reads the mesh via
    `dolfinx.io.XDMFFile.read_mesh` and then optional cell and facet
    tag sets via `read_meshtags`. Default grid names are `mesh`,
    `cell_tags`, and `facet_tags`; override with
    `xdmf_mesh_name`, `xdmf_cell_tags_name`, and
    `xdmf_facet_tags_name` if the file uses different ones.
  - Physical groups in the `.msh` (or named meshtag sets in the
    `.xdmf`) supply region and facet tags.
  - For axisymmetric benchmarks, the meridian mesh is built from a
    `.geo` template (e.g.
    [`benchmarks/moscap_axisym_2d/moscap_axisym.geo`](../../benchmarks/moscap_axisym_2d/moscap_axisym.geo))
    and a generated `.msh` is committed alongside it.

## `coordinate_system` and axisymmetric constraints

When `coordinate_system == "axisymmetric"`, cross-field validation
in [`semi/schema.py`](../../semi/schema.py) enforces:

1. `dimension == 2`.
2. The radial extent is non-negative; the first mesh coordinate is `r ≥ 0`.
3. No Dirichlet contact lives on the symmetry axis `r = 0`. The axis is a natural (no-flux) boundary.
4. The outer radial wall `r = R` is treated as homogeneous Neumann; choose `R` large enough that the solution is insensitive to the cutoff (`R ≳ 5 W_dmax` for MOSCAP).

See [`docs/theory/axisymmetric.md`](../theory/axisymmetric.md) for the derivation.

## `doping` profiles

Each entry targets one region with one profile. Densities are in
cm⁻³ (device-physics convention) and are converted to m⁻³ internally.

- `uniform`: `N_D`, `N_A`.
- `step`: `axis`, `location`, `N_{D,A}_{left,right}`.
- `gaussian`: `axis`, `center`, `sigma`, `peak`, `donor` (bool).

## `contacts`

| Field | Notes |
|---|---|
| `name` | Symbolic contact name, used in output IV files. |
| `facet` | Either an int tag or a string referencing `mesh.facets_by_plane[*].name`. |
| `type` | `"ohmic"`, `"gate"`, `"insulating"`. |
| `voltage` | DC value in V (default 0). |
| `voltage_sweep` | `{start, stop, step}` for bias-sweep runners. |
| `workfunction` | Required for `type: "gate"` (V). |
| `oxide_region` | Required for gate; the insulator region key. |

## `solver`

| `solver.type` | Runner | Notes |
|---|---|---|
| `"equilibrium"` | `run_equilibrium` | Equilibrium Poisson only. |
| `"drift_diffusion"` | `run` | Coupled Slotboom DD at fixed bias. |
| `"bias_sweep"` | `run_bias_sweep` | Adaptive continuation along a contact's `voltage_sweep`. |
| `"mos_cv"` | `run_mos_cv` | 2D MOS C–V via `numpy.gradient(Q, V)`. |
| `"mos_cap_ac"` | `run_mos_cap_ac` | C–V via analytic dQ/dV (M14.1). |
| `"transient"` | `run_transient` | BDF1/BDF2 (M13). Fields: `t_end`, `dt`, `order`, `max_steps`, `output_every`. |
| `"ac_sweep"` | `run_ac_sweep` | Small-signal AC at a DC operating point (M14). |

Common `solver` sub-keys: `tolerances`, `continuation`, `jacobian_shift`.

## Minimal 1D pn-junction example

```json
{
  "schema_version": "1.0.0",
  "name": "pn_junction_1d",
  "dimension": 1,
  "mesh": {
    "source": "builtin",
    "extents": [[0.0, 2.0e-6]],
    "resolution": [400],
    "facets_by_plane": [
      {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
      {"name": "cathode", "tag": 2, "axis": 0, "value": 2.0e-6}
    ]
  },
  "regions": {
    "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"}
  },
  "doping": [
    {
      "region": "silicon",
      "profile": {
        "type": "step", "axis": 0, "location": 1.0e-6,
        "N_D_left": 0.0,  "N_A_left": 1.0e17,
        "N_D_right": 1.0e17, "N_A_right": 0.0
      }
    }
  ],
  "contacts": [
    {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
    {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0}
  ]
}
```

## Axisymmetric example

[`benchmarks/moscap_axisym_2d/moscap_axisym.json`](../../benchmarks/moscap_axisym_2d/moscap_axisym.json)
is the live reference for the axisymmetric path.

## Known schema caveats (carried forward from M11)

- **Permissive about unknown keys.** Top-level and most nested
  objects lack `additionalProperties: false`, so a UI typo in a
  field name validates and is silently ignored. The manifest
  schema is strict; flipping the input schema is a v2.0.0 change.
- **No per-block `schema_version`.** A bump in the `physics` block
  forces every caller to upgrade atomically.

## See also

- [`schemas/manifest.v1.json`](../../schemas/manifest.v1.json) — output manifest schema (M9).
- [`semi/schema.py`](../../semi/schema.py) — loader and version gate.
- [`semi/io/`](../../semi/io/) — manifest writer and reader.
