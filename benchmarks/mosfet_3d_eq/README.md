# 3D planar n-MOSFET, equilibrium (M19 precursor)

3D planar n-MOSFET on a gmsh-sourced unstructured mesh, equilibrium
Poisson solve at V_GS = V_DS = V_BS = 0 V. No bias sweep, no Pao-Sah
verifier; those are M19 proper deliverables. See
`docs/M19_PRECURSOR_STARTER_PROMPT.md` and ADR 0019 for the full
scope statement.

## Physics scope

- Single semiconductor material (Si) with a Gaussian n+ implant
  profile in the source / drain regions on a p-type body.
  Doping-only S/D differentiation per ADR 0019 (one Si region
  with tag 1; the JSON `doping` block carries the profiles).
- SiO2 gate oxide as a second region (tag 4).
- Four contacts: source (tag 10), drain (tag 11), gate (tag 12),
  body (tag 13), all at V = 0 V.
- Boltzmann statistics, SRH recombination only.
- M19 proper adds the V_GS bias sweep, the M16.4 FD statistics
  path for the heavy n+ implants, and the Pao-Sah analytical
  comparison verifier.

## Verification target

The equilibrium verifier `semi/verification/mosfet_3d_eq.py`
gates five geometric / shape-based properties at V = 0:

1. `n_max` in the source / drain implant volumes within a factor
   of 3 of N_D peak (1e20 cm^-3).
2. `n_min` in the channel (50 nm below the oxide, channel
   x window) below 1e12 cm^-3.
3. `p_max` in the body bulk (>500 nm below the surface) within
   10 % of N_A (1e17 cm^-3).
4. Charge neutrality |n - p - N_D + N_A| < 1e15 cm^-3 in the
   body bulk.
5. Built-in junction potential: psi at the lateral edges of the
   S/D implants is at least 3 * V_t below psi_body.

## Geometry summary

| Dimension      | Value     |
|----------------|-----------|
| L_x (total)    | 1.5 um    |
| L_y (total)    | 0.5 um    |
| L_z (Si depth) | 2.0 um    |
| L_g (channel)  | 250 nm    |
| T_ox (oxide)   | 5 nm      |
| x_g_lo .. x_g_hi | 0.625 .. 0.875 um (channel x window) |

L_y is 0.5 um (half the prompt's nominal 1 um) so the cell
count stays inside the 200k-500k budget at the prompt's oxide
thickness. ADR 0019 records the deviation; widening L_y back
to 1 um is part of M19 proper.

## Mesh quality summary (as shipped)

| Region    | Max edge (as shipped) | Prompt nominal |
|-----------|-----------------------|----------------|
| Gate oxide | 3 nm                  | 1 nm           |
| Channel Si | 15 nm                 | 10 nm          |
| Bulk Si    | 150 nm                | 100 nm         |

Isotropic tet meshing at 1 nm in the oxide volume produces
millions of cells in the lateral footprint; the precursor ships
at the values above. ADR 0019 records the deviation. The mesh
quality module shipped in Phase B gates the as-shipped thresholds.

Observed mesh stats on the committed `mosfet_3d.msh`:

| Quantity           | Value    |
|--------------------|----------|
| Total cells (tets) | 210883   |
| Si body cells      | 70005    |
| Oxide cells        | 140878   |
| Vertices           | 47386    |
| Source facets      | 2558     |
| Drain facets       | 2641     |
| Gate facets        | 32152    |
| Body facets        | 90       |

gmsh wall time for regeneration: ~9 s in the dev image.

## Regeneration

The benchmark ships the generated `.msh` under `fixtures/` so CI
does not need gmsh at runtime (the `fixtures/` subdirectory is
the gitignore whitelist exception per the resistor_3d precedent;
the bare `benchmarks/<name>/*.msh` location is git-ignored).
Regenerate from the committed `.geo` with:

```
docker compose run --rm dev gmsh -3 -format msh4 \
  benchmarks/mosfet_3d_eq/fixtures/mosfet_3d.geo \
  -o benchmarks/mosfet_3d_eq/fixtures/mosfet_3d.msh
```

gmsh version: tested against the dev image's gmsh 4.14.0 with
OCC 7.6.3. The `.geo` works internally in micrometres and uses
`Mesh.ScalingFactor = 1e-6` to write meters at output (OCC's
hardcoded `Precision::Confusion` is 1e-7, which would collapse
the 5 nm oxide if the geometry were authored in metres
directly). Mesh node coordinates in the `.msh` file are in
meters, matching `semi/mesh.py`'s ingest convention.

## Status

- Phase A: gmsh `.geo` + generated `.msh`. Done.
- Phase B: mesh quality module + ingest hook. Pending.
- Phase C: schema additions for `mesh.quality_gate`. Pending.
- Phase D: this benchmark's JSON config. Pending.
- Phase E: equilibrium verifier + CI matrix entry. Pending.
- Phase F: PLAN / README / ROADMAP / CHANGELOG sweep, ADR 0019
  Accepted. Pending.
