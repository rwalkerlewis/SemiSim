# Axisymmetric 2D MOSCAP with LF/HF C-V Curves

Source: user-provided GitHub Copilot task, 2026-04-30. Authoritative.

## Context

- Branch off `main` at HEAD `6176903` (M14 audit cases 02/05 closed).
- A prior agent attempt exists on `origin/copilot/extend-2d-moscap-axisymmetric-again`. Treat that as discardable history; start a fresh branch `feature/axisymmetric-moscap-cv` from `main`.
- A planar 2D MOSCAP already exists at `benchmarks/mos_2d/mos_cap.json`. Extend it; do not duplicate.
- Drift-diffusion + Slotboom modules already exist under `semi/physics/` (`drift_diffusion.py`, `slotboom.py`). Reuse them; the user prompt's note about DD being unmerged is stale.
- Repo conventions in `/memories/repo/SemiSim.md` (no em dashes, mesh in physical meters, scaled Poisson coefficient `lambda2 * L0^2 * eps_r`, etc.) apply.

## Goal

Reproduce Hu Fig. 5-18 (LF vs. HF C-V split for an MOS capacitor) using a cylindrical-axisymmetric 2D model in the meridian half-plane.

## Task (verbatim user spec)

> Extend the existing 2D MOSCAP example in `rwalkerlewis/SemiSim` to an
> axisymmetric (cylindrical) 2D model, and reproduce the low-frequency
> (quasi-static) vs. high-frequency C-V curve split shown in Fig. 5-18 of
> Chenming Hu, Modern Semiconductor Devices for Integrated Circuits,
> Chapter 5 (https://www.chu.berkeley.edu/wp-content/uploads/2020/01/Chenming-Hu_ch5-1.pdf).

### Physics: axisymmetric formulation

MOSCAP is a circular gate stack with rotational symmetry about z. Solve
in the meridian half-plane $(r, z) \in [0, R] \times [-T_{Si}, T_{ox}]$:
silicon body $z \in [-T_{Si}, 0]$, oxide $z \in [0, T_{ox}]$, gate at
$z = T_{ox}$ for $r \le R_g$.

Volume measure: $dV = 2\pi r\, dr\, dz$. Multiply every volume integrand
by `r` (the $2\pi$ cancels). This applies to Poisson AND continuity.

BCs:

| Boundary | Where | Condition |
|---|---|---|
| Symmetry axis | $r = 0$ | Natural (no Dirichlet). |
| Outer radial | $r = R$ | Homogeneous Neumann; choose $R \ge 5 W_{dmax}$. |
| Bulk contact | $z = -T_{Si}$ | Ohmic: $\psi = \psi_{bulk}$, $\Phi_n = \Phi_p = 0$. |
| Gate | $z = T_{ox}, r \le R_g$ | Dirichlet $\psi = V_g - V_{fb}$. |
| Field oxide top | $z = T_{ox}, R_g < r \le R$ | Insulating Neumann. |

Use gmsh via `dolfinx.io.gmshio` for tagged subdomains and facets, with
graded refinement near the Si/SiO2 interface and the gate edge.

### LF vs. HF C-V

For each $V_g$ in the sweep:
1. Solve coupled Poisson-DD (Slotboom) DC.
2. LF: $C_{LF} = -dQ_s/dV_g$ via centered finite difference of total
   semiconductor charge. (Hu Eq. 5.6.1.)
3. HF: small-signal solve of Poisson with minority delta-charge frozen,
   OR depletion-approximation clamp $W_{dep} = W_{dmax}$ once $\phi_s
   \ge 2\phi_B$ (state which method is used in docstring + notebook).
4. Plot $C/C_{ox}$ vs $V_g$, both curves on one axis.

### Reference parameters

- p-type Si, $N_a = 5 \times 10^{16}$ cm^-3
- $T_{ox} = 10$ nm SiO2
- N+ poly-Si gate
- $T = 300$ K
- $R_g = 50$ um, $R = 200$ um, $T_{Si} = 5$ um
- $V_g$ sweep: -2.0 V to +2.0 V, 81 points

### Files to add or modify

- `semi/physics/poisson_axisym.py` (or `coordinate_system: "axisymmetric"`
  switch in existing `poisson.py`)
- `semi/physics/drift_diffusion_axisym.py` (or same switch pattern)
- `semi/mesh.py`: extend with axisymmetric 2D, OR add gmsh `.geo` under
  `benchmarks/moscap_axisym_2d/`
- `semi/schema.py` + `schemas/input.v1.json`: accept
  `coordinate_system: "axisymmetric"` (default `"cartesian"`); validate
  `dimension: 2` and meridian extent
- `semi/cv.py`: `compute_lf_cv(...)` and `compute_hf_cv(...)`
- `benchmarks/moscap_axisym_2d/moscap_axisym.json`
- `benchmarks/moscap_axisym_2d/reference_cv.csv` (analytical Cox, Cmin,
  Vfb, Vt)
- `notebooks/05_moscap_axisym_cv.ipynb` mirroring
  `notebooks/01_pn_junction_1d.ipynb`
- `tests/check_axisym_moscap_math.py` (no dolfinx)
- `tests/test_moscap_axisym_cv.py` (pytest)
- `CHANGELOG.md` entry

### Required figures (commit to `notebooks/figures/` or repo equivalent)

1. Geometry sketch of meridian half-plane.
2. Mesh plot with refinement zoom at gate corner.
3. Equilibrium $\psi(r,z)$ at $V_g = 0$.
4. $n(r,z)$ and $p(r,z)$ at deep accumulation, threshold, strong inversion.
5. Band diagram along $r=0$ at the same three biases.
6. C-V plot with $V_{fb}, V_t$ marked, side-by-side with Hu Fig 5-18 thumbnail.
7. Convergence study: $C_{HF,min}$ vs $R$ and vs mesh refinement.

At least one figure must show the field revolved around the z-axis
(pyvista `extrude_rotate` or equivalent).

## Branching and PR workflow

1. `git checkout -b feature/axisymmetric-moscap-cv` from `main`.
2. Conventional Commits, small reviewable commits per concern (schema,
   mesh, physics, solver, benchmark, notebook, tests, docs).
3. Push, open PR titled `feat: axisymmetric 2D MOSCAP with LF/HF C-V curves`.
4. PR description: physics summary, embedded images, acceptance
   checklist tied to files/tests, deviations note, "Open questions"
   section if anything is ambiguous.
5. Run `pytest`; all green.
6. Do NOT merge.

## Acceptance criteria (verbatim)

- [ ] All work on `feature/axisymmetric-moscap-cv`; PR open against `main`; no commits to `main`.
- [ ] Existing pytest passes; new tests pass.
- [ ] CI green on PR.
- [ ] Axisymmetric weak forms r-weighted; documented in module docstrings.
- [ ] No Dirichlet at $r=0$.
- [ ] JSON schema accepts `coordinate_system` with validation.
- [ ] LF returns to $C_{ox}$ in strong inversion (within 2%).
- [ ] HF stays at $C_{min}$ in strong inversion (within 2%).
- [ ] LF and HF coincide for $V_g < V_t$.
- [ ] C-V plot visually consistent with Hu Fig 5-18.
- [ ] All seven figures committed and embedded in the notebook.
- [ ] Notebook runs top-to-bottom; README Colab badge section extended.
- [ ] PR description includes figures, checklist, physics writeup.

## Constraints

- Honor invariants in `PLAN.md` (five-layer architecture, JSON-as-contract,
  no PETSc types across engine API).
- No em dashes in any prose, code comments, commit messages, or PR body.
- Do not lower CI gates.
- Reuse existing Slotboom DD path; do not re-derive sign conventions.
- Mesh coordinates in physical meters; scaled Poisson coefficient is
  `L_D^2 * eps_r = lambda2 * L0^2 * eps_r`, not bare `lambda2 * eps_r`.

## Reporting format

End with:
1. Branch and PR URL.
2. List of new/modified files.
3. `pytest` summary line.
4. Any deviations or open questions for review.
