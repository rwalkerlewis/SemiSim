"""
Build the Colab notebook 03 that walks through the 2D MOS capacitor C-V
benchmark. M6 content: multi-region (Si/SiO2) equilibrium Poisson
with gate contact BCs and a quasi-static C-V sweep against the
depletion-approximation theory curve.
"""
from pathlib import Path

import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(r"""# SemiSim · Notebook 03: MOS Capacitor C-V (2D, Si/SiO$_2$)

First 2D and first multi-region benchmark in SemiSim. A p-type silicon substrate (500 nm, $N_A = 10^{17}$ cm⁻³) carries a 5 nm SiO$_2$ gate oxide; the gate contact sweeps from deep accumulation through depletion and into strong inversion. The notebook:

1. Runs multi-region equilibrium Poisson at $V_{gate} = 0$ and renders the 2D $\psi(x, y)$ contour.
2. Runs the full $V_{gate}$ sweep and compares the simulated $C(V_{gate}) = dQ_{gate}/dV_{gate}$ against the depletion-approximation theory curve.

The multi-region plumbing is a single dolfinx mesh with cellwise $\varepsilon_r$ (Si = 11.7, SiO$_2$ = 3.9) on a DG0 function, space-charge restricted to silicon cells via `dx(subdomain_id=semi_tag)`, and the oxide carrying only the Laplacian term. Gate contacts enter as voltage-shifted Dirichlet BCs on $\psi$ (offset by $\phi_{ms}$), with natural zero-flux at the Si/SiO$_2$ interface.

## Notebook set

This is the third of four Colab walkthroughs that exercise the [end-of-M7 capability matrix](https://github.com/rwalkerlewis/SemiSim#status) from the README:

| Notebook | Covers |
|----------|--------|
| [01 — Equilibrium Poisson on a 1D pn junction](./01_pn_junction_1d.ipynb) | Nonlinear Poisson with Boltzmann carriers; depletion-approximation comparison |
| [02 — Bias sweep on the 1D pn junction](./02_pn_junction_bias.ipynb) | Coupled Slotboom drift-diffusion with SRH; Shockley forward + SNS reverse |
| **03 — MOS capacitor C-V** (this one) | Gate contact BCs with $\phi_{ms}$; Si/SiO$_2$ multi-region; depletion-approximation C-V |
| [04 — 3D doped resistor V-I](./04_resistor_3d.ipynb) | Ohmic V-I linearity; builtin box mesh vs. gmsh `.msh` fixture |

**What this notebook does:**
1. Installs `dolfinx` on Colab via [FEM on Colab](https://fem-on-colab.github.io/) (one wget, ~30 s)
2. Clones the [SemiSim](https://github.com/rwalkerlewis/SemiSim) repo
3. Loads `benchmarks/mos_2d/mos_cap.json`
4. Solves equilibrium at $V_{gate} = 0$ and plots the 2D $\psi$ field with `matplotlib.tri.Triangulation` + `tricontourf`
5. Runs the full sweep $V_{gate} \in [-0.9, +1.2]$ V in 0.05 V steps (43 bias points)
6. Extracts $C_{sim}(V_{gate}) = dQ_{gate}/dV_{gate}$ by centered finite difference and overlays the depletion-approximation theory

Wall time on Colab: ~10-20 minutes. The sweep is the heavy part (43 multi-region 2D solves on a 4 × 505 triangulated mesh, ~2.5 k DOFs each).
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 1. Install dolfinx

~30 s on first run. Cached across cell re-runs within the same Colab session."""))

cells.append(nbf.v4.new_code_cell(r"""try:
    import dolfinx
    print(f"dolfinx already present: version {dolfinx.__version__}")
except ImportError:
    !wget -q "https://fem-on-colab.github.io/releases/fenicsx-install-release-real.sh" -O "/tmp/fenicsx-install.sh" && bash "/tmp/fenicsx-install.sh"
    import dolfinx
    print(f"dolfinx installed: version {dolfinx.__version__}")

major, minor = (int(x) for x in dolfinx.__version__.split('.')[:2])
assert (major, minor) >= (0, 10), f"Need dolfinx >= 0.10; got {dolfinx.__version__}"
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 2. Clone the repo and install the package

Same pattern as Notebooks 01 and 02 — `pip install -e` so any source edits you make persist in the Colab session."""))

cells.append(nbf.v4.new_code_cell(r"""import os
if not os.path.exists('SemiSim'):
    !git clone -q https://github.com/rwalkerlewis/SemiSim.git

%cd SemiSim
!pip install -q -e .
print("Package installed.")
"""))

cells.append(nbf.v4.new_code_cell(r"""import copy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from semi import schema, run, materials, constants

print(f"semi version: {__import__('semi').__version__}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 3. Load the MOS benchmark

The JSON drives the simulator: device geometry, doping, contacts, solver type. For C-V sweeps `solver.type == "mos_cv"` dispatches to `semi.runners.mos_cv.run_mos_cv`, which loops over `contacts[gate].voltage_sweep` and records `(V_gate, Q_gate)` into `result.iv` at each bias."""))

cells.append(nbf.v4.new_code_cell(r"""cfg = schema.load("benchmarks/mos_2d/mos_cap.json")

sweep = next(c for c in cfg["contacts"] if c["type"] == "gate")["voltage_sweep"]
print(f"Case:       {cfg['name']}")
print(f"Dimension:  {cfg['dimension']}D")
print(f"Mesh:       {cfg['mesh']['source']}, "
      f"{cfg['mesh']['resolution'][0]} x {cfg['mesh']['resolution'][1]} cells")
print(f"Regions:    {[(r, cfg['regions'][r]['material'], cfg['regions'][r]['role']) for r in cfg['regions']]}")
print(f"Doping:     uniform N_A={cfg['doping'][0]['profile']['N_A']:.0e} cm^-3 in silicon")
print(f"Gate sweep: V from {sweep['start']} V to {sweep['stop']} V in {sweep['step']} V steps "
      f"({int(round((sweep['stop']-sweep['start'])/sweep['step'])) + 1} points)")
print(f"Solver:     {cfg['solver']['type']}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 4. Equilibrium at $V_{gate} = 0$ and the 2D $\psi$ contour

Run a single-point "sweep" at $V_{gate} = 0$ so the runner returns the equilibrium potential field. At this bias the device sits in weak depletion (the body's ohmic BC pins $\psi_{body} = -\phi_F \approx -0.417$ V and the ideal gate with $\phi_{ms} = 0$ pulls the surface up toward 0 V).

**Visualization choice.** This notebook renders the 2D $\psi$ field with `matplotlib.tri.Triangulation` + `tricontourf`, identical to what `scripts/run_benchmark.py` produces for the benchmark. This choice is deliberate:

- The dolfinx mesh is already a 2D triangulation, and `tricontourf` accepts the dof coordinates directly — no resampling, no export step.
- Matplotlib is already a Colab dependency. The alternative (VTX export + `pyvista`) would pull in VTK + trame on top of FEM-on-Colab's fenicsx install, adding 30-60 s of install time and a second rendering backend to maintain for a plot that `tricontourf` handles in one call.
- `pyvista` is the right tool for 3D surfaces / isosurfaces / interactive rotation, none of which apply to a 2D contour plot on a rectangular domain."""))

cells.append(nbf.v4.new_code_cell(r"""# Clone the config with a single-point sweep at V_gate = 0 so we can
# render the 2D psi field at equilibrium before running the full sweep.
cfg_eq = copy.deepcopy(cfg)
for c in cfg_eq["contacts"]:
    if c["type"] == "gate":
        c["voltage_sweep"] = {"start": 0.0, "stop": 0.0, "step": 0.05}

result_eq = run.run(cfg_eq)
info = result_eq.solver_info
print(f"SNES at V_gate = 0 V: iterations = {info.get('iterations', 'n/a')}, "
      f"converged = {info.get('converged', 'n/a')}, reason = {info.get('reason', 'n/a')}")
print(f"Q_gate at V_gate = 0 V: {result_eq.iv[-1]['Q_gate']*1e4:+.3f} uC/cm^2")
"""))

cells.append(nbf.v4.new_code_cell(r"""x = result_eq.x_dof[:, 0]
y = result_eq.x_dof[:, 1]
psi = np.asarray(result_eq.psi_phys)

# Si/SiO2 interface y-coordinate (top of the silicon region's bounds).
si_region = cfg_eq["regions"]["silicon"]
ox_region = cfg_eq["regions"]["oxide"]
si_bounds_y = next(rb["bounds"][1] for rb in cfg_eq["mesh"]["regions_by_box"]
                    if int(rb["tag"]) == int(si_region["tag"]))
y_interface_nm = float(si_bounds_y[1]) * 1e9
x_extent_nm = [float(cfg_eq["mesh"]["extents"][0][0]) * 1e9,
               float(cfg_eq["mesh"]["extents"][0][1]) * 1e9]

tri = Triangulation(x * 1e9, y * 1e9)

fig, ax = plt.subplots(figsize=(7, 5))
tcf = ax.tricontourf(tri, psi, levels=20, cmap="viridis")
cbar = fig.colorbar(tcf, ax=ax)
cbar.set_label(r"$\psi$ (V)")
ax.plot(x_extent_nm, [y_interface_nm, y_interface_nm],
        color="white", lw=1.0, ls="--", label="Si/SiO$_2$ interface")
ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_title(r"MOS cap: $\psi(x, y)$ at $V_{gate}$ = 0 V")
ax.set_aspect("auto")
ax.legend(loc="upper left", fontsize=8)
ax.grid(False)
fig.tight_layout()
plt.show()

# Central-column vertical slice for sanity.
xmid = 0.5 * sum(cfg_eq["mesh"]["extents"][0])
near_mid = np.abs(x - xmid) < 1.0e-9
if not near_mid.any():
    near_mid = np.abs(x - xmid) < (np.abs(x - xmid).min() + 1.0e-15)
ys = y[near_mid]
psi_mid = psi[near_mid]
order = np.argsort(ys)
ys_nm = ys[order] * 1e9
psi_mid_s = psi_mid[order]

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(ys_nm, psi_mid_s, "-", lw=1.5)
ax.axvline(y_interface_nm, color="red", ls=":", label="Si/SiO$_2$ interface")
ax.set_xlabel("y (nm)")
ax.set_ylabel(r"$\psi$ (V)")
ax.set_title(r"Central-column $\psi(y)$ at $V_{gate}$ = 0 V")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
plt.show()

print(f"  psi at body (y=0):          {psi_mid_s[0]:+.4f} V  (expect ~-phi_F = -0.417 V)")
print(f"  psi at Si surface (y=500nm): "
      f"{psi_mid_s[np.argmin(np.abs(ys_nm - y_interface_nm))]:+.4f} V")
print(f"  psi at gate (y=505nm):       {psi_mid_s[-1]:+.4f} V")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 5. C-V sweep and theory overlay

Now the full sweep: $V_{gate}$ from $-0.9$ V through $+1.2$ V in 0.05 V steps (43 bias points). At each step the runner assembles the multi-region equilibrium Poisson form, solves with PETSc SNES, and integrates the silicon space charge over $\Omega_{Si}$ to recover the gate charge per unit area:

$$Q_{gate}(V_{gate}) \;=\; -\frac{q}{W_{\text{lat}}}\,\int_{\Omega_{Si}} \rho\,dA,\qquad C_{sim}(V_{gate}) = \frac{dQ_{gate}}{dV_{gate}}.$$

The theory overlay is the textbook depletion-approximation C-V:

$$\frac{1}{C_{\text{theory}}(V_{gate})} \;=\; \frac{1}{C_{ox}} + \frac{1}{C_{dep}(\psi_s)},\qquad C_{dep}(\psi_s) = \sqrt{\frac{\varepsilon_s q N_A}{2\psi_s}}$$

with $\psi_s(V_{gate})$ the surface potential that solves the control equation $V_{gate} - V_{FB} = \psi_s + \sqrt{2 \varepsilon_s q N_A \psi_s} / C_{ox}$."""))

cells.append(nbf.v4.new_code_cell(r"""result_sweep = run.run(cfg)
info = result_sweep.solver_info
print(f"Sweep complete: {len(result_sweep.iv)} bias points accepted")
print(f"Final SNES (at V_gate={result_sweep.iv[-1]['V']:+.3f} V): "
      f"iterations = {info.get('iterations', 'n/a')}, "
      f"converged = {info.get('converged', 'n/a')}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""### Verifier-window disclosure: $V_{FB}+0.1 \to V_{FB}+0.2$

The benchmark's verifier window is $V_{gate} \in [V_{FB}+0.2,\ V_T-0.1]$ V, **not** $[V_{FB}+0.1,\ V_T-0.1]$ V as originally planned on M6.

**What happened.** The first M6 implementation used $V_{FB}+0.1$ as the low edge. Under our $\psi=0$-at-intrinsic BC convention the worst relative error in that window was **10.06%** at $V_{gate} = -0.25$ V, just past the 10% tolerance. This is a depletion-approximation modeling limit: at $\psi_s \lesssim 2 V_t$ the free-carrier tail extends across a significant fraction of the depletion width and the sharp-edge approximation breaks down toward 10%. It is not a solver accuracy problem — a finer mesh does not move the number.

**The fix.** Per the reviewer's standing guidance ("if the verifier's tolerance is held at 10%, the right debugging action on a near-miss is to shrink the window, not loosen the tolerance"), the low edge moved from $V_{FB}+0.1$ to $V_{FB}+0.2$. The new window's worst error is **9.25% at $V_{gate} = -0.200$ V** (the low edge itself). See `docs/mos_derivation.md` §6.9 and `docs/PHYSICS.md` §6 for the full derivation and disclosure.

**Why surfacing this here matters.** A SemiSim reviewer reading the Colab sees the design choice where it applies, instead of having to cross-reference the PR body or CHANGELOG. The simulated C-V is correct outside the verifier window too; what is gated is the narrow regime where the *analytical* reference curve is trustworthy."""))

cells.append(nbf.v4.new_code_cell(r"""from semi.materials import get_material
from semi.constants import Q, EPS0, cm3_to_m3

# Device parameters for the theory curve. Mirrors
# scripts.run_benchmark._mos_device_params / _mos_C_theory so the
# notebook comparison is identical to the benchmark's verifier curve.
mat = get_material(cfg["regions"]["silicon"]["material"])
ox_mat = get_material(cfg["regions"]["oxide"]["material"])
N_A_SI = cm3_to_m3(cfg["doping"][0]["profile"]["N_A"])
eps_s = mat.epsilon
eps_ox = ox_mat.epsilon
ox_bounds = next(
    rb for rb in cfg["mesh"]["regions_by_box"]
    if int(rb["tag"]) == int(cfg["regions"]["oxide"]["tag"])
)
t_ox = float(ox_bounds["bounds"][1][1] - ox_bounds["bounds"][1][0])
gate = next(c for c in cfg["contacts"] if c["type"] == "gate")
phi_ms = float(gate.get("workfunction", 0.0) or 0.0)
V_t = result_sweep.scaling.V0
phi_F = V_t * float(np.log(N_A_SI / mat.n_i))
C_ox = eps_ox / t_ox
V_FB = phi_ms - phi_F
Q_B_thresh = float(np.sqrt(4.0 * eps_s * Q * N_A_SI * phi_F))
V_T = V_FB + 2.0 * phi_F + Q_B_thresh / C_ox

def psi_s_of_V(V_gate):
    V_ov = V_gate - V_FB
    if V_ov <= 0.0:
        return float("nan")
    cap = 2.0 * phi_F + Q_B_thresh / C_ox + 1.0e-9
    if V_ov >= cap:
        return float("nan")
    a = float(np.sqrt(2.0 * eps_s * Q * N_A_SI)) / C_ox
    u = (-a + float(np.sqrt(a * a + 4.0 * V_ov))) / 2.0
    return float(u * u)

def C_theory(V_arr):
    out = np.full_like(V_arr, np.nan, dtype=float)
    for i, Vg in enumerate(V_arr):
        psi_s = psi_s_of_V(float(Vg))
        if not np.isfinite(psi_s) or psi_s <= 0.0:
            continue
        C_dep = float(np.sqrt(eps_s * Q * N_A_SI / (2.0 * psi_s)))
        out[i] = C_ox * C_dep / (C_ox + C_dep)
    return out

# Extract C_sim from the iv table by centered finite difference in V_gate.
V = np.array([r["V"] for r in result_sweep.iv])
Q_gate = np.array([r["Q_gate"] for r in result_sweep.iv])
order = np.argsort(V)
V = V[order]
Q_gate = Q_gate[order]
C_sim = np.gradient(Q_gate, V)
C_th = C_theory(V)

window_lo = V_FB + 0.2
window_hi = V_T - 0.1
print(f"  phi_F  = {phi_F:.3f} V, C_ox = {C_ox*1e2:.3f} uF/cm^2")
print(f"  V_FB   = {V_FB:+.3f} V, V_T = {V_T:+.3f} V")
print(f"  window = [{window_lo:+.3f}, {window_hi:+.3f}] V (= [V_FB+0.2, V_T-0.1])")

# Worst error inside the window (excluding array endpoints where FD is one-sided).
endpoint = np.ones_like(V, dtype=bool); endpoint[0] = False; endpoint[-1] = False
mask = (V >= window_lo) & (V <= window_hi) & np.isfinite(C_th) & endpoint
if mask.any():
    rel = np.abs(C_sim[mask] - C_th[mask]) / C_th[mask]
    w = int(np.argmax(rel))
    print(f"  worst window error: {rel[w]*100:.2f}% at V_gate = {V[mask][w]:+.3f} V "
          f"(expect ~9.25% at V_gate=-0.200 V per the benchmark)")
"""))

cells.append(nbf.v4.new_code_cell(r"""fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(V, C_sim * 1e2, "o-", color="C0", ms=4, lw=1.8, label="simulation (SemiSim)")
finite = np.isfinite(C_th)
ax.plot(V[finite], C_th[finite] * 1e2, "r--", lw=1.5, label="depletion-approx theory")
ax.axhline(C_ox * 1e2, color="k", ls=":", lw=1.0, label=r"$C_{ox}$")
ax.axvspan(window_lo, window_hi, alpha=0.15, color="green",
           label=r"verifier window $[V_{FB}+0.2,\,V_T-0.1]$")
ax.axvline(V_FB, color="grey", ls=":", alpha=0.6)
ax.axvline(V_T, color="grey", ls=":", alpha=0.6)
ax.set_xlabel(r"$V_{gate}$ (V)")
ax.set_ylabel(r"$C$ ($\mu$F/cm$^2$)")
ax.set_title(r"MOS C-V vs depletion-approximation theory")
ax.grid(True, alpha=0.3)
ax.legend(loc="best", fontsize=8)
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(V, Q_gate * 1e2, "o-", ms=3, lw=1.5)
ax.axvspan(window_lo, window_hi, alpha=0.15, color="green")
ax.axvline(V_FB, color="grey", ls=":", alpha=0.6, label=rf"$V_{{FB}} = {V_FB:+.3f}$ V")
ax.axvline(V_T, color="grey", ls=":", alpha=0.6, label=rf"$V_T = {V_T:+.3f}$ V")
ax.set_xlabel(r"$V_{gate}$ (V)")
ax.set_ylabel(r"$Q_{gate}$ ($\mu$C/cm$^2$)")
ax.set_title(r"$Q_{gate}(V_{gate})$ (integrated silicon space charge)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## Summary

- **2D equilibrium** at $V_{gate} = 0$ V renders cleanly via `tricontourf` on the dolfinx mesh's native triangulation. The body contact pins $\psi_{body} = -\phi_F$ and the ideal gate pulls the surface toward 0 V (weak depletion at zero gate bias under our BC convention).
- **C-V sweep** (43 bias points, $[-0.9, +1.2]$ V) reproduces the depletion branch between $V_{FB}$ and $V_T$ and drops monotonically from near $C_{ox}$ in accumulation through a minimum at strong-inversion onset.
- **Depletion-approximation agreement**: worst error 9.25% at $V_{gate} = -0.200$ V, inside the benchmark's 10% tolerance.
- **Verifier-window shift $V_{FB}+0.1 \to V_{FB}+0.2$** (see §5 above): a 10.06% near-miss at the original low edge was resolved by shrinking the window per the "hold the tolerance, shrink the window" reviewer rule, not by loosening the 10% bound.
- **Regimes outside the verifier window** — accumulation ($V_{gate} < V_{FB}$, hole pile-up at the interface) and strong inversion ($V_{gate} > V_T$, minority-carrier inversion layer) — are rendered on the plot for reference but are not gated, because the depletion approximation does not model them and quasi-static inversion recovery needs a transient solver.

**Next up.** Notebook 04 moves from 2D quasi-static C-V to a 3D doped resistor V-I sweep, compares the builtin box mesh to a gmsh `.msh` fixture, and anchors both against the ohmic line $I = V/R$.

See [the repo](https://github.com/rwalkerlewis/SemiSim) and its [capability matrix](https://github.com/rwalkerlewis/SemiSim#status) for the full scope."""))

nb.cells = cells

out_path = Path(__file__).resolve().parent.parent / "notebooks" / "03_mos_cv.ipynb"
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as f:
    nbf.write(nb, f)
print(f"Wrote {out_path} with {len(cells)} cells")
