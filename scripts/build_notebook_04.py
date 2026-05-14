"""
Build the Colab notebook 04 that walks through the 3D doped resistor
V-I benchmark. M7 content: 3D ohmic V-I linearity on a builtin box
mesh and on a gmsh `.msh` fixture, compared to R = L / (q N_D mu_n A).
"""
from pathlib import Path

import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(r"""# SemiSim · Notebook 04: 3D Doped Resistor V-I

First 3D benchmark in SemiSim. A 1 µm × 200 nm × 200 nm uniformly n-doped silicon bar ($N_D = 10^{18}$ cm⁻³, $\mu_n = 1400$ cm$^2$/(V·s)) with ohmic contacts on the $x = 0$ and $x = L$ faces, swept over $V \in [-0.01, +0.01]$ V at 0.005 V steps. The notebook:

1. Runs the bipolar sweep on the **builtin** box mesh (`[64, 16, 16]` hex-as-tet cells) and overlays the simulated $I(V)$ against the ohmic line $I = V / R$ with $R = L / (q\, N_D\, \mu_n\, A)$.
2. Runs the same sweep on the **gmsh `.msh` fixture** (`fixtures/box.msh`, an unstructured tetrahedral mesh of the same geometry) to confirm the `dolfinx.io.gmsh.read_from_msh` path in `semi/mesh.py` produces the same $R_{sim}$ within 1%.

The two runs share a JSON schema, a physics stack, and a verifier. What differs is *where the mesh comes from*: the builtin variant is generated in-process by `semi.mesh._build_from_builtin`; the gmsh variant is parsed from the committed `.msh` file by `semi.mesh._build_from_file`. Both are 3D tetrahedral meshes by the time the runner sees them.

## Notebook set

This is the fourth and final Colab walkthrough that exercises the [end-of-M7 capability matrix](https://github.com/rwalkerlewis/SemiSim#status) from the README:

| Notebook | Covers |
|----------|--------|
| [01 — Equilibrium Poisson on a 1D pn junction](./01_pn_junction_1d.ipynb) | Nonlinear Poisson with Boltzmann carriers; depletion-approximation comparison |
| [02 — Bias sweep on the 1D pn junction](./02_pn_junction_bias.ipynb) | Coupled Slotboom drift-diffusion with SRH; Shockley forward + SNS reverse |
| [03 — MOS capacitor C-V (2D multi-region)](./03_mos_cv.ipynb) | Gate contact BCs with $\phi_{ms}$; Si/SiO$_2$ multi-region; depletion-approximation C-V |
| **04 — 3D doped resistor V-I** (this one) | Ohmic V-I linearity; builtin box mesh vs. gmsh `.msh` fixture |

**What this notebook does:**
1. Installs `dolfinx` on Colab via [FEM on Colab](https://fem-on-colab.github.io/) and ensures the `gmsh` Python module is present so the `.msh` loader works (see §1 for the Option-A pin).
2. Clones the [SemiSim](https://github.com/rwalkerlewis/SemiSim) repo.
3. Loads `benchmarks/resistor_3d/resistor.json` (builtin mesh), runs the 5-point bipolar sweep, and plots $I(V)$ vs. $V/R$.
4. Loads `benchmarks/resistor_3d/resistor_gmsh.json` (gmsh fixture), runs the same sweep, overlays on the same axes.
5. Compares $R_{sim}$ across the two meshes; confirms both agree with $R_{theory}$ to within 1%.

**Wall time on Colab.** Budget 5-10 minutes. 3D assemblies are roughly an order of magnitude slower than the 2D MOS sweep per bias point, but the sweep is only 5 points per variant (10 solves total).
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 1. Install dolfinx and gmsh

Two installs in this cell:

- **`dolfinx`** via FEM-on-Colab's `release-real` build (~30 s on first run, identical to Notebooks 01-03).
- **`gmsh`** via `pip` (**Option A pin**). The `.msh` loader in `semi/mesh.py` imports `from dolfinx.io import gmsh as _gmsh_io`, and `dolfinx.io.gmsh.read_from_msh` calls into the `gmsh` Python module at parse time. FEM-on-Colab's release image does not always ship the `gmsh` package bundled, and when it is missing the import surfaces as a `ModuleNotFoundError` only at the moment the first gmsh benchmark runs — by which point the notebook has already walked the reviewer through the builtin variant and it is too late to gracefully recover. Pinning `!pip install gmsh` up-front eliminates that failure mode; the install is idempotent and cheap (~5 s) so there is no cost to running it even when gmsh is already present."""))

cells.append(nbf.v4.new_code_cell(r"""try:
    import dolfinx
    print(f"dolfinx already present: version {dolfinx.__version__}")
except ImportError:
    !wget -q "https://fem-on-colab.github.io/releases/fenicsx-install-release-real.sh" -O "/tmp/fenicsx-install.sh" && bash "/tmp/fenicsx-install.sh"
    import dolfinx
    print(f"dolfinx installed: version {dolfinx.__version__}")

major, minor = (int(x) for x in dolfinx.__version__.split('.')[:2])
assert (major, minor) >= (0, 10), f"Need dolfinx >= 0.10; got {dolfinx.__version__}"

# Option A: ensure the gmsh Python module is available for the .msh loader.
# gmsh's Python bindings link against libGLU even when used headlessly.
# Colab's FEM-on-Colab image does not always ship libGLU, so install it
# explicitly before the pip install. This is idempotent if already present.
!apt-get install -y -q libglu1-mesa
!pip install -q gmsh
import gmsh as _gmsh_probe
print(f"gmsh Python module present: version {_gmsh_probe.__version__}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 2. Clone the repo and install the package

Same pattern as Notebooks 01-03 — `pip install -e` so any source edits you make persist in the Colab session."""))

cells.append(nbf.v4.new_code_cell(r"""import os
if not os.path.exists('SemiSim'):
    !git clone -q https://github.com/rwalkerlewis/SemiSim.git

%cd SemiSim
!pip install -q -e .
print("Package installed.")
"""))

cells.append(nbf.v4.new_code_cell(r"""import numpy as np
import matplotlib.pyplot as plt

from semi import schema, run
from semi.constants import Q, cm3_to_m3
from semi.materials import get_material

print(f"semi version: {__import__('semi').__version__}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 3. Builtin-mesh variant: `benchmarks/resistor_3d/resistor.json`

The JSON drives the simulator: device geometry, doping, contacts, solver continuation. For this benchmark the 3D box is constructed in-process by `dolfinx.mesh.create_box` with the extents and resolution given under `mesh.extents` and `mesh.resolution`. The `contacts` block pins `V_{left} = 0` and sweeps `V_{right}` from $-0.01$ V to $+0.01$ V in $0.005$ V steps.

The bipolar continuation in `semi.runners.bias_sweep` walks $0 \to -V_{max} \to +V_{max}$ (each intermediate bias visited twice because the path crosses zero twice); the verifier dedupes by $V$ after the sweep, and we do the same below."""))

cells.append(nbf.v4.new_code_cell(r"""cfg_builtin = schema.load("benchmarks/resistor_3d/resistor.json")

sweep_right = next(c for c in cfg_builtin["contacts"] if c.get("voltage_sweep"))
sweep = sweep_right["voltage_sweep"]
ext = cfg_builtin["mesh"]["extents"]
res = cfg_builtin["mesh"]["resolution"]
print(f"Case:       {cfg_builtin['name']}")
print(f"Dimension:  {cfg_builtin['dimension']}D")
print(f"Mesh:       {cfg_builtin['mesh']['source']}, "
      f"{res[0]} x {res[1]} x {res[2]} cells "
      f"over {ext[0][1]*1e6:.2f} um x {ext[1][1]*1e9:.0f} nm x {ext[2][1]*1e9:.0f} nm")
print(f"Doping:     uniform N_D = {cfg_builtin['doping'][0]['profile']['N_D']:.0e} cm^-3")
print(f"Sweep:      V_{sweep_right['name']} from {sweep['start']} V to {sweep['stop']} V "
      f"in {sweep['step']} V steps")
print(f"Solver:     {cfg_builtin['solver']['type']}")
"""))

cells.append(nbf.v4.new_code_cell(r"""result_builtin = run.run(cfg_builtin)

info = result_builtin.solver_info
print(f"Sweep complete: {len(result_builtin.iv)} iv rows recorded "
      f"(bipolar path visits intermediate points twice)")
print(f"Final SNES: iterations = {info.get('iterations', 'n/a')}, "
      f"converged = {info.get('converged', 'n/a')}, "
      f"reason = {info.get('reason', 'n/a')}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""### $R_{theory}$ and $R_{sim}$ extraction

Analytical resistance for the uniformly n-doped ohmic bar (see `docs/resistor_derivation.md` §2):

$$R_{theory} \;=\; \frac{L}{q\, N_D\, \mu_n\, A}$$

with $L = 1$ µm, $A = (200\,\text{nm})^2 = 4 \times 10^{-14}$ m$^2$, $N_D = 10^{24}$ m⁻³, $\mu_n = 1400$ cm$^2$/(V·s) $= 0.14$ m$^2$/(V·s), $q = 1.6022 \times 10^{-19}$ C. That gives $R_{theory} \approx 1.116$ kΩ.

The runner records `{"V", "J"}` pairs into `result.iv` where `J` is the facet-averaged current density on the swept contact (A/m$^2$). Converting to current: $I = J \cdot A$. The geometry is read from `mesh.extents` for the builtin variant and from the loaded mesh's bounding box for the gmsh variant — mirroring `scripts/run_benchmark.py::_resistor_geometry` so the notebook's comparison matches the verifier's line-for-line."""))

cells.append(nbf.v4.new_code_cell(r"""# Geometry helper: mirrors scripts.run_benchmark._resistor_geometry so the
# notebook's R_theory uses the same L, A the verifier uses.
def _resistor_geometry(cfg, mesh_obj):
    mesh_cfg = cfg["mesh"]
    if mesh_cfg.get("source") == "builtin":
        ext = mesh_cfg["extents"]
        Lx = float(ext[0][1] - ext[0][0])
        Wy = float(ext[1][1] - ext[1][0])
        Wz = float(ext[2][1] - ext[2][0])
        return Lx, Wy * Wz
    x = mesh_obj.geometry.x
    Lx = float(x[:, 0].max() - x[:, 0].min())
    Wy = float(x[:, 1].max() - x[:, 1].min())
    Wz = float(x[:, 2].max() - x[:, 2].min())
    return Lx, Wy * Wz


def _doping_and_mobility(cfg):
    prof = cfg["doping"][0]["profile"]
    N_D = cm3_to_m3(float(prof["N_D"]))
    mu_n_SI = float(cfg["physics"]["mobility"].get("mu_n", 1400.0)) * 1.0e-4
    return N_D, mu_n_SI


# Dedup iv by V (bipolar path re-visits intermediate points) and convert
# J (A/m^2) to I (A) by multiplying by the swept-contact area.
def _extract_IV(result, cfg):
    L_x, A_contact = _resistor_geometry(cfg, result.mesh)
    N_D, mu_n_SI = _doping_and_mobility(cfg)
    R_theory = L_x / (Q * N_D * mu_n_SI * A_contact)

    by_V = {}
    for row in result.iv:
        by_V[round(float(row["V"]), 9)] = float(row["J"])
    V = np.array(sorted(by_V.keys()))
    J = np.array([by_V[v] for v in V])
    I = J * A_contact
    return V, I, R_theory, L_x, A_contact


V_b, I_b, R_th_b, L_b, A_b = _extract_IV(result_builtin, cfg_builtin)
nonzero = np.abs(V_b) > 1.0e-12
R_sim_b = V_b[nonzero] / I_b[nonzero]
err_b = np.abs(R_sim_b - R_th_b) / R_th_b
worst_b = int(np.argmax(err_b))

print(f"  L = {L_b*1e6:.3f} um,  A = {A_b*1e18:.3f} (nm)^2 "
      f"(W^2 = {(A_b**0.5)*1e9:.1f} nm square)")
print(f"  R_theory = {R_th_b:.3e} Ohm  ({R_th_b*1e-3:.3f} kOhm)")
print(f"  I(V=0) = {float(I_b[~nonzero][0]):+.3e} A  (should be numerical noise)")
for v, i, r in zip(V_b[nonzero], I_b[nonzero], R_sim_b, strict=True):
    print(f"    V = {v:+.4f} V   I = {i:+.3e} A   R_sim = {r:.3e} Ohm "
          f"(rel err {abs(r-R_th_b)/R_th_b*100:.3f}%)")
print(f"  Worst V-I linearity error (builtin): {err_b[worst_b]*100:.3f}% "
      f"at V = {V_b[nonzero][worst_b]:+.4f} V  (verifier tolerance: 1.000%)")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 4. gmsh-fixture variant: `benchmarks/resistor_3d/resistor_gmsh.json`

Same device, same sweep, but the mesh now comes from `benchmarks/resistor_3d/fixtures/box.msh` via `semi.mesh._build_from_file`. That function calls `dolfinx.io.gmsh.read_from_msh`, which requires the `gmsh` Python module installed in step 1. Physical groups in the `.msh` carry the region tag (1 = silicon) and facet tags (10 = contact_left, 20 = contact_right), so the JSON references them by integer rather than by axis-plane like the builtin variant does.

The fixture is an unstructured tetrahedral mesh of the same rectangular bar, produced offline by `fixtures/box.geo` and committed to the repo so the test is reproducible.

**Why this matters.** It proves that the benchmark is not a builtin-mesh artifact: the verifier passes at 1% whether the mesh is the structured hex-as-tet grid or an unstructured Delaunay tet mesh with heterogeneous element sizes. It also regression-tests the `.msh` parse path, which is the on-ramp for arbitrary external geometries (the next realistic step past a rectangular bar)."""))

cells.append(nbf.v4.new_code_cell(r"""cfg_gmsh = schema.load("benchmarks/resistor_3d/resistor_gmsh.json")

sweep_g = next(c for c in cfg_gmsh["contacts"] if c.get("voltage_sweep"))["voltage_sweep"]
print(f"Case:       {cfg_gmsh['name']}")
print(f"Mesh:       {cfg_gmsh['mesh']['source']} ({cfg_gmsh['mesh']['format']}), "
      f"path={cfg_gmsh['mesh']['path']}")
print(f"Doping:     uniform N_D = {cfg_gmsh['doping'][0]['profile']['N_D']:.0e} cm^-3")
print(f"Sweep:      V from {sweep_g['start']} V to {sweep_g['stop']} V "
      f"in {sweep_g['step']} V steps (5 points)")
"""))

cells.append(nbf.v4.new_code_cell(r"""result_gmsh = run.run(cfg_gmsh)

info = result_gmsh.solver_info
print(f"Sweep complete: {len(result_gmsh.iv)} iv rows recorded")
print(f"Final SNES: iterations = {info.get('iterations', 'n/a')}, "
      f"converged = {info.get('converged', 'n/a')}, "
      f"reason = {info.get('reason', 'n/a')}")

msh = result_gmsh.mesh
x_coords = msh.geometry.x
n_cells_local = msh.topology.index_map(msh.topology.dim).size_local
print(f"gmsh mesh: {x_coords.shape[0]} vertices, {n_cells_local} local cells "
      f"(unstructured tetrahedra)")

V_g, I_g, R_th_g, L_g, A_g = _extract_IV(result_gmsh, cfg_gmsh)
nonzero_g = np.abs(V_g) > 1.0e-12
R_sim_g = V_g[nonzero_g] / I_g[nonzero_g]
err_g = np.abs(R_sim_g - R_th_g) / R_th_g
worst_g = int(np.argmax(err_g))

print(f"\n  gmsh bounding-box geometry: L = {L_g*1e6:.3f} um, "
      f"A = {A_g*1e18:.3f} (nm)^2")
print(f"  R_theory (gmsh geom)    = {R_th_g:.3e} Ohm  ({R_th_g*1e-3:.3f} kOhm)")
for v, i, r in zip(V_g[nonzero_g], I_g[nonzero_g], R_sim_g, strict=True):
    print(f"    V = {v:+.4f} V   I = {i:+.3e} A   R_sim = {r:.3e} Ohm "
          f"(rel err {abs(r-R_th_g)/R_th_g*100:.3f}%)")
print(f"  Worst V-I linearity error (gmsh): {err_g[worst_g]*100:.3f}% "
      f"at V = {V_g[nonzero_g][worst_g]:+.4f} V  (verifier tolerance: 1.000%)")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 5. V-I linearity comparison

Both meshes are overlaid on a single axes against the ohmic line $I = V / R_{theory}$. With $R_{theory} \approx 1.116$ kΩ and the sweep confined to $|V| \le 10$ mV, the bar never leaves the low-field linear regime and the three curves should lie on top of each other to within 1%."""))

cells.append(nbf.v4.new_code_cell(r"""R_theory = R_th_b  # identical device geometry; builtin extents are exact
V_line = np.linspace(V_b.min() * 1.05, V_b.max() * 1.05, 50)
I_line = V_line / R_theory

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(V_line * 1e3, I_line * 1e6, "k-", lw=1.2,
        label=rf"ohmic line $I = V / R$  ($R = {R_theory*1e-3:.3f}$ k$\Omega$)")
ax.plot(V_b * 1e3, I_b * 1e6, "o", color="C0", ms=8,
        label="builtin box mesh (64x16x16)")
ax.plot(V_g * 1e3, I_g * 1e6, "s", color="C3", ms=6, fillstyle="none", mew=1.8,
        label=r"gmsh fixture (`fixtures/box.msh`)")
ax.axhline(0.0, color="grey", lw=0.6, alpha=0.5)
ax.axvline(0.0, color="grey", lw=0.6, alpha=0.5)
ax.set_xlabel(r"$V$ (mV)")
ax.set_ylabel(r"$I$ ($\mu$A)")
ax.set_title(r"3D resistor V-I: simulation vs. $I = V / R$")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", fontsize=9)
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(7, 3.8))
ax.plot(V_b[nonzero] * 1e3, err_b * 100, "o-", color="C0", ms=7, lw=1.5,
        label="builtin")
ax.plot(V_g[nonzero_g] * 1e3, err_g * 100, "s-", color="C3", ms=6, lw=1.5,
        fillstyle="none", mew=1.8, label="gmsh fixture")
ax.axhline(1.0, color="k", ls="--", lw=1.0,
           label="1% verifier tolerance")
ax.set_xlabel(r"$V$ (mV)")
ax.set_ylabel(r"$|R_{sim} - R_{theory}| / R_{theory}$ (%)")
ax.set_title(r"V-I linearity error vs. $R_{theory}$")
ax.grid(True, alpha=0.3)
ax.legend(loc="best", fontsize=9)
fig.tight_layout()
plt.show()

delta_R = (R_sim_g.mean() - R_sim_b.mean()) / R_sim_b.mean() * 100
print(f"  Mean R_sim (builtin): {R_sim_b.mean()*1e-3:.4f} kOhm")
print(f"  Mean R_sim (gmsh):    {R_sim_g.mean()*1e-3:.4f} kOhm")
print(f"  gmsh vs builtin:      {delta_R:+.3f}%  (expected < 1% by the benchmark design)")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## Summary

- **3D ohmic V-I** on a uniformly n-doped silicon bar passes the V-I linearity verifier at 1% on both the builtin `create_box` mesh and the `gmsh` `.msh` fixture, confirming the simulator extends cleanly from 1D/2D to 3D tetrahedra.
- **Mesh-source equivalence** — `R_{sim}$ from the unstructured gmsh mesh matches $R_{sim}$ from the structured builtin mesh within a fraction of a percent, well inside the 1% tolerance. The `.msh` loader path in `semi/mesh.py::_build_from_file` is therefore a drop-in alternative to the builtin path: a user who needs arbitrary external geometry can switch by editing the JSON and committing an `.msh` fixture, with no physics-code change required.
- **Option-A install pin** — the explicit `pip install gmsh` in §1 is what makes `schema.load("resistor_gmsh.json")` work on Colab. `dolfinx.io.gmsh.read_from_msh` requires the `gmsh` Python module at parse time, and FEM-on-Colab's `release-real` image does not always bundle it; pinning `gmsh` up-front avoids a late `ModuleNotFoundError` that would otherwise only surface in §4.
- **Where to find the derivation and the verifier**:
  - Analytical resistance formula, device-parameter table, and the rationale for the $|V| \le 10$ mV sweep window: `docs/resistor_derivation.md` §2.
  - V-I linearity verifier specification (5-point sweep, zero-bias floor, sign check, 1% tolerance): `docs/resistor_derivation.md` §3.
  - Verifier implementation: `verify_resistor_3d()` in `scripts/run_benchmark.py`, registered under the key `"resistor_3d"`.

**Notebook set complete.** Notebooks 01-04 together walk a new reviewer through all four device classes in the end-of-M7 capability matrix: 1D equilibrium Poisson, 1D coupled drift-diffusion under bias, 2D multi-region MOS C-V, and 3D ohmic V-I with two mesh sources. The rest of the simulator's capability — MMS-verified convergence rates, charge-neutrality conservation, SRH reverse-bias generation, adaptive continuation diagnostics — is exercised by the test suite (206 tests, 95.58% coverage) and the V&V battery (62/62 PASS). See the [repo](https://github.com/rwalkerlewis/SemiSim) and its [capability matrix](https://github.com/rwalkerlewis/SemiSim#status)."""))

nb.cells = cells

out_path = Path(__file__).resolve().parent.parent / "notebooks" / "04_resistor_3d.ipynb"
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as f:
    nbf.write(nb, f)
print(f"Wrote {out_path} with {len(cells)} cells")
