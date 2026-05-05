# 19 — Linear solvers and GPU backends

## Learning objectives

- Explain what MUMPS does (sparse direct LU factorization) and why it
  is the default per-Newton-step linear solver in the engine.
- Recognize the GPU-accelerated alternatives — AMGX, hypre BoomerAMG —
  as algebraic-multigrid preconditioners for Krylov methods.
- Read the `solver.backend` and `solver.compute` schema fields and trace
  how they thread to PETSc options.
- State the engine's bit-identity guarantee: the GPU path produces the
  same answer as the CPU path to floating-point precision *by default*.
- Identify the `GET /capabilities` endpoint in the HTTP server and
  understand how the UI gates a backend dropdown on runtime availability.

## Physical motivation

Every Newton iteration solves a sparse linear system $J\delta\mathbf{u}
= -F$. For sub-100k DOF problems, MUMPS LU factorization is robust and
gives tight error control; for ~1M DOF problems (3D MOSFETs, FinFETs),
direct LU memory and time both blow up super-linearly and we need
iterative methods with strong preconditioners. The M15 milestone added
PETSc-CUDA / PETSc-HIP with AMGX or hypre BoomerAMG; the goal is a 5×
linear-solve speedup on a 500 kDOF benchmark, which the
`benchmarks/poisson_3d_gpu/` test gates on a nightly self-hosted runner.

This chapter is about *numerical* methods, not physics, but the
right-tool-for-the-right-size choice matters because the linear solver
typically dominates wall time once the problem is non-trivial. The
shipped engine keeps the CPU path as the reference; the GPU is opt-in
per ADR's M15 acceptance criterion.

## Derivation from first principles

### Sparse direct LU (MUMPS)

For a sparse matrix $A \in \mathbb{R}^{N\times N}$ with $\mathrm{nnz}(A)$
non-zeros, an LU factorization $A = LU$ writes $A$ as a lower-triangular
times an upper-triangular matrix. Forward and backward substitution
solve $L\mathbf{y} = \mathbf{b}$ and $U\mathbf{x} = \mathbf{y}$ in
$O(\mathrm{nnz}(L) + \mathrm{nnz}(U))$ work each. The asymptotic memory
cost of LU on a 3D problem is $O(N^{4/3})$ in the worst case (vs
$O(N)$ for the original $A$); the time cost is $O(N^2)$. Beyond
~100 kDOF on a single workstation, this scaling is prohibitive.

**MUMPS** (Multifrontal Massively Parallel sparse direct Solver) is a
mature multi-frontal sparse LU library. The dolfinx-PETSc default is
MUMPS-LU; kronos-semi inherits this:

```python
DEFAULT_PETSC_OPTIONS = {
    ...
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}
```

`ksp_type: preonly` means "apply the preconditioner once, no Krylov
iteration" — i.e. use the LU factorization as a direct solve, not as a
preconditioner. MUMPS handles ill-conditioning robustly; the M14
transient runner bumps `mat_mumps_icntl_14: 200` (200% extra workspace)
to absorb delayed pivots from the Jacobian shift (Ch. 17).

### Iterative solvers and AMG

For large sparse problems, iterative Krylov methods (CG, GMRES, BiCGSTAB)
solve $A\mathbf{x} = \mathbf{b}$ by building a sequence of approximations
$\mathbf{x}^{(k)}$ that converge to $\mathbf{x}$. The *number* of iterations
to reach tolerance scales with $\sqrt{\kappa(A)}$ for SPD systems and
$\kappa(A)$ for non-SPD; for the discrete DD Jacobian, $\kappa$ can be
in the millions. A **preconditioner** $M^{-1}$ approximates $A^{-1}$ so
that $M^{-1}A$ has $\kappa = O(1)$ and the iteration count is bounded.

**Algebraic multigrid (AMG)** preconditioners exploit the
elliptic-like character of the Poisson stiffness: errors at coarse
length scales are cheap to resolve on a coarse grid. AMG builds a
hierarchy of grids and smoothers from $A$ alone (no geometric mesh
information needed) and is highly effective for diffusion-dominated
problems.

- **AMGX** (NVIDIA): GPU-native AMG, written in CUDA. Best on H100/A100/V100.
- **hypre BoomerAMG**: classical AMG, originally CPU but with a CUDA
  port (`pc_hypre_use_gpu = true`). Works on both NVIDIA and AMD with
  the right build.

The engine's `solver.backend` field selects between `cpu-mumps` (default),
`gpu-amgx`, `gpu-hypre`, and `auto` (best available, falls back to CPU
if no device).

### Backend resolution and options threading

[`semi/compute.py`](../../semi/compute.py) (M15 deliverable) probes the
PETSc build at runtime and returns the available backends. The runtime
probe checks:
- Is PETSc compiled with CUDA / HIP?
- Is AMGX available?
- Is hypre's GPU path available?

[`semi/solver.py:152-183`](../../semi/solver.py) (`_resolve_backend_options`)
merges the user's `solver.backend` choice with the defaults and the
`solver.compute` overrides:

```json
"solver": {
  "type": "equilibrium",
  "backend": "gpu-amgx",
  "compute": {
    "device": "cuda",
    "linear_solver": "gmres",
    "preconditioner": "amgx"
  }
}
```

Backend → PETSc options translation:
- `cpu-mumps`: `ksp_type=preonly, pc_type=lu, pc_factor_mat_solver_type=mumps`.
- `gpu-amgx`: `ksp_type=gmres, pc_type=amgx, mat_type=aijcusparse, vec_type=cuda`.
- `gpu-hypre`: `ksp_type=gmres, pc_type=hypre, pc_hypre_use_gpu=true`.

The CPU path stays bit-identical: setting `backend: "cpu-mumps"` (or
omitting the field) gives the *exact same* PETSc options as before M15,
and the resulting answers match v0.14.1 to floating-point precision.

### Bit-identity guarantee

The GPU path is **opt-in**. The default behaviour is `cpu-mumps`, which
is the CI-tested, MMS-verified, V&V-confirmed truth. The GPU path's
acceptance test ([`docs/IMPROVEMENT_GUIDE.md` §M15](../IMPROVEMENT_GUIDE.md))
demands that the GPU output is within $10^{-8}$ relative L2 of the
CPU-MUMPS reference on every existing benchmark. The GPU is therefore
a *speed* upgrade, not a *physics* upgrade — a deliberate ADR design
choice to keep the verification truth on a single, well-understood path.

### `GET /capabilities` endpoint

The HTTP server ([`kronos_server/routes/health.py`](../../kronos_server/routes/health.py))
exposes `/capabilities` returning `{"backends": [...], "device": ...,
"compute": {...}}`. The UI uses this to gate the backend dropdown:
greyed-out if `gpu-amgx` is not available on this engine build,
selectable if it is. The UI never has to introspect the PETSc build
itself.

### Manifest fields

For traceability, the M15 manifest (`schemas/manifest.v1.json` v1.1.0)
adds `solver.ksp_iters` (Krylov iterations per Newton step) and
`solver.linear_solve_wall_s` (wall time of the linear solve). These
are populated from PETSc's KSP iteration counter and a `time.monotonic()`
bracket around the solve call ([`semi/solver.py:230-247`](../../semi/solver.py)).

## Key results

- MUMPS LU: robust direct solver, dominates wall time at ~10kDOF+.
- AMG preconditioning: $\kappa(M^{-1}A) = O(1)$ for elliptic problems.
- `solver.backend` schema field: routes between CPU and GPU paths.
- `solver.compute`: fine-grained PETSc options.
- Bit-identity: CPU path is the reference; GPU within $10^{-8}$ rel L2.
- `GET /capabilities`: runtime backend availability for UI integration.

## Worked numerical example

**M15 acceptance benchmark** (`benchmarks/poisson_3d_gpu/`):
- Mesh: $80\times 80\times 80$ uniform 3D box → $531441$ vertices.
- Equation: pure Poisson with cellwise $\varepsilon_r$ and a Gaussian
  source.
- Hardware: NVIDIA A100 (80 GB) or AMD MI250.

Acceptance criterion: $T_\mathrm{cpu}/T_\mathrm{gpu} \ge 5$ on the
linear-solve portion alone (excluding mesh build, IO).

Typical wall times reported in the M15 acceptance run:
- CPU MUMPS LU on a 32-core EPYC: ~120 s linear solve.
- GPU AMGX on A100: ~22 s linear solve.
- Speedup: 5.5×, clears the gate.

Output match: max-norm relative difference $\|\psi_\mathrm{gpu} - \psi_\mathrm{cpu}\|_\infty/\|\psi_\mathrm{cpu}\|_\infty
\sim 10^{-9}$, comfortably below the $10^{-8}$ rel-L2 acceptance gate.

## Code map

| Concept | Code location |
|---|---|
| Default PETSc options | `semi/solver.py:20-30` (`DEFAULT_PETSC_OPTIONS`) |
| Backend resolution | `semi/solver.py:152-183` (`_resolve_backend_options`) |
| Backend probe | `semi/compute.py` (M15) |
| GPU AMGX options | (per backend dict in `semi/compute.py::backend_settings_from_cfg`) |
| `KSP iters` and `linear_solve_wall_s` | `semi/solver.py:230-247` |
| `GET /capabilities` | `kronos_server/routes/health.py` |
| Schema fields `solver.backend, solver.compute` | `schemas/input.v2.json` |
| Manifest 1.1.0 KSP fields | `schemas/manifest.v1.json` |
| `KRONOS_BACKEND` env override | `semi/compute.py` |

## Existing-docs cross-reference

- [`docs/gpu.md`](../gpu.md) — definitive M15 reference with install recipes.
- [`docs/IMPROVEMENT_GUIDE.md` §M15, §5](../IMPROVEMENT_GUIDE.md) — original deliverable and GPU strategy.
- [`docs/adr/0008-snes-tolerances.md`](../adr/0008-snes-tolerances.md) — SNES tolerances that drive the linear-solve precision.
- [`docs/ROADMAP.md`](../ROADMAP.md) — capability matrix entry for M15.

## Common pitfalls

1. **GPU isn't free.** The first frequency / Newton step pays a
   factorization-setup cost on the device; subsequent reuse benefits.
   For very small problems (sub-1k DOF), the launch overhead exceeds
   the speedup; always profile your specific workload.
2. **Real PETSc only.** The GPU path is in real PETSc only;
   `solver.type: "ac_sweep"` (which uses a complex linear system) falls
   back to `cpu-mumps` for the AC linear-solve portion regardless of
   `backend`. The DC operating-point sub-solve still respects the
   requested backend ([`docs/gpu.md`](../gpu.md) Known limits).
3. **Double precision only.** The `solver.compute.precision` field is
   reserved for fp32 paths; any value other than `"double"` is
   rejected at schema validation today.
4. **Strict GPU = explicit failure.** Asking for `backend: "gpu-amgx"`
   on a CPU-only PETSc build is a hard error (no silent fallback). Use
   `backend: "auto"` to opt into the fallback. Acceptance test A3 in
   ADR M15 documents this.
5. **Assembly stays on the host.** dolfinx 0.10's UFL assembly is
   CPU-bound; only the linear solve runs on the device. For problems
   where assembly is the bottleneck (very nonlinear, high-quadrature),
   the GPU path saves only the linear-solve fraction. dolfinx-cuda is
   experimental and not yet wired in (M15+ deferred work).

## Exercises

**Exercise 19.1.** A 1D 100k-DOF problem takes 12 s of MUMPS LU. A 3D
1M-DOF problem with the same per-DOF density takes... how long, by
the $O(N^2)$ scaling estimate? Why is this a useful sanity check for
when GPU AMG becomes worth it?

**Exercise 19.2.** Read [`semi/compute.py::backend_settings_from_cfg`](../../semi/compute.py).
What PETSc option dict does `gpu-amgx` produce, and which key tells
PETSc to allocate the matrix on the GPU?

**Exercise 19.3.** A user submits a JSON with `backend: "gpu-amgx"` on
a CPU-only build. What happens at runtime, per the M15 acceptance
test A3? Compare with `backend: "auto"`.

**Exercise 19.4.** Check the manifest of a v0.16.0+ run on a benchmark
of your choice. Look for `solver.ksp_iters` and `solver.linear_solve_wall_s`.
Why are these useful diagnostics?

**Exercise 19.5.** A heterojunction simulation (M17, planned) will
double the number of unknowns in regions with chi(x) variation. How
does this affect the choice of CPU-MUMPS vs GPU-AMG?

### Solutions

**19.1.** $O(N^2)$: $T_\mathrm{3D}/T_\mathrm{1D} = (10^6/10^5)^2 = 100$.
So 1200 s = 20 min of CPU MUMPS for the 3D problem. With the M15
5× speedup target, GPU AMG would do this in 4 min — qualitatively a
different user experience. The 100k-DOF crossover is the rough
threshold.

**19.2.** Reading the engine source, `gpu-amgx` produces approximately
`{"ksp_type": "gmres", "pc_type": "amgx", "mat_type": "aijcusparse",
"vec_type": "cuda", "pc_factor_mat_solver_type": None}`. The
`mat_type: aijcusparse` key tells PETSc the matrix lives in a CUDA-aware
sparse format on the device. The `None` for `pc_factor_mat_solver_type`
unsets the CPU-MUMPS factorization key.

**19.3.** Strict request, no fallback: the engine raises an error at
solve time, citing that `gpu-amgx` is not available on this build.
The user must either install AMGX or change to `cpu-mumps` / `auto`.
With `auto`: the engine probes available backends, picks the best
present (e.g. `cpu-mumps` if no GPU), and prints a warning to the
manifest's `warnings` list.

**19.4.** Look in `runs/<run_id>/manifest.json` under `solver`. The
fields are populated for every run. Useful diagnostics:
- `ksp_iters > 100` flags either ill-conditioning or AMG setup
  failure on the GPU path.
- `linear_solve_wall_s` proportional to total wall time → linear
  solve is the bottleneck → GPU is worth attempting.
- Across a bias sweep, comparing per-step values shows where Newton
  conditioning degrades.

**19.5.** Heterojunctions add chi(x) and Eg(x) as cellwise DG0 fields
(M17 deliverable). The DOF count itself doesn't double — these are
*coefficients* in existing UFL forms, not new unknowns. What does
change is the *coupling structure*: at a heterojunction interface, the
band-edge discontinuity introduces a thin transition layer where the
solution varies rapidly. AMG may converge slower on this layer (the
interpolation operators within AMG typically assume smoother solutions);
the M17 acceptance test will surface any AMG-tuning needed.

## Further reading

- **Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed.
  (2003).** The Krylov-method reference. Free online.
- **Trottenberg, Oosterlee, and Schuller, *Multigrid* (2001).** The
  classic AMG textbook.
- **MUMPS user guide** (Amestoy et al., 2001-current). Distributed
  with the MUMPS package.
- **PETSc manual:** https://petsc.org/release/manual/. Specifically
  the SNES, KSP, PC, and Mat sections.
- **AMGX docs:** https://github.com/NVIDIA/AMGX
- **hypre docs:** https://hypre.readthedocs.io/
- **`docs/gpu.md`** in this repo — installation recipes and the
  acceptance numbers.
