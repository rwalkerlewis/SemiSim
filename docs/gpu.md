# GPU linear-solver path (M15)

SemiSim can run the linear solve at each Newton step on the GPU
through PETSc's CUDA or HIP build. The CPU-MUMPS path remains the
default and is bit-identical to v0.14.1.

## Selecting a backend

Set `solver.backend` in the JSON input:

```json
{
  "solver": {
    "type": "equilibrium",
    "backend": "gpu-amgx",
    "compute": {
      "device": "cuda",
      "linear_solver": "gmres",
      "preconditioner": "amgx"
    }
  }
}
```

| `solver.backend` | What happens |
|---|---|
| `cpu-mumps` (default) | Direct LU on CPU via MUMPS. Bit-identical to pre-M15. |
| `gpu-amgx` | KSP on GPU with the AMGX algebraic-multigrid PC. |
| `gpu-hypre` | KSP on GPU with hypre BoomerAMG (`pc_hypre_use_gpu=true`). |
| `auto` | Best available; falls back to `cpu-mumps` on CPU-only PETSc. |

Querying availability programmatically:

```python
from semi import compute
print(compute.available_backends())   # e.g. ["cpu-mumps", "gpu-amgx"]
print(compute.device_info())          # PETSc build info
```

The HTTP server exposes the same data at `GET /capabilities` so the
UI can gate the backend dropdown.

### Override at run time

The environment variable `SEMISIM_BACKEND` overrides
`solver.backend` if set. Useful for switching the entire benchmark
suite onto the GPU for one CI run without editing JSON files.

## Install recipes

### PETSc + CUDA + AMGX (recommended)

```bash
# conda
conda install -c conda-forge petsc=3.21 petsc4py "petsc=*=cuda*" amgx
```

```dockerfile
# In the SemiSim Dockerfile, replace the petsc apt line with
RUN pip install --no-binary petsc4py "petsc[cuda]" amgx-python
```

### PETSc + HIP + hypre (AMD)

```bash
# Build PETSc from source against ROCm
./configure --with-hip=1 --with-hypre=1 --download-hypre \
            --with-hypre-extra-configure-arguments="--with-cuda=no --with-hip=yes"
make all
```

### hypre + CUDA

```bash
./configure --with-cuda=1 --download-hypre \
            --with-hypre-extra-configure-arguments="--with-cuda=yes --enable-unified-memory"
```

## Acceptance numbers (A2)

The 3D Poisson box at `benchmarks/poisson_3d_gpu/` (uniform doping,
80x80x80 = 531441 nodes) is required to deliver

    cpu_wall_s / gpu_wall_s >= 5.0

on a single A100 / MI250 with PETSc 3.21+ and AMGX 2.4+ or hypre 2.32+.
The nightly self-hosted CI job at `.github/workflows/gpu-nightly.yml`
gates this number.

## Known limits

- **Assembly stays on the host.** Only the linear solve runs on the
  device. Form assembly is still CPU-side per the dolfinx 0.10 API.
- **Real scalars only.** Complex PETSc scalars are not yet wired
  through the GPU path (used by the small-signal AC sweep), so
  `solver.type = "ac_sweep"` falls back to `cpu-mumps` for the
  AC-side LinearProblem; the DC operating-point sub-solve still
  picks up the requested backend.
- **Double precision only.** The `solver.compute.precision` field
  is reserved for a future fp32 path; for now any value other than
  `"double"` is rejected at schema validation.
- **No silent fallback (acceptance test A3).** Asking for an
  unavailable explicit GPU backend is a hard error. Use
  `backend: "auto"` to opt in to the fallback behaviour.
