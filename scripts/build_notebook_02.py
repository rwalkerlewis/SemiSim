"""
Build the Colab notebook 02 that walks through the 1D pn-junction bias
sweep. M2 forward branch (Shockley / SNS reference) followed by M3
reverse branch (SRH-generation reference).
"""
from pathlib import Path

import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(r"""# SemiSim · Notebook 02: Forward and Reverse Bias on a 1D pn Junction

Couples Poisson to drift-diffusion continuity equations in Slotboom $(\psi, \Phi_n, \Phi_p)$ variables, sweeps the anode bias, and compares the simulated $J(V)$ curve to analytical references on both branches:

- **Forward bias** (0 to 0.6 V): Sah-Noyce-Shockley $J_\text{diff} + J_\text{rec}$, with the ideal long-diode Shockley diffusion curve as a second overlay. Current is dominated by depletion-region SRH recombination at low bias and by diffusion at high bias.
- **Reverse bias** (0 to -2 V): net SRH generation in the expanding depletion region, $J_\text{gen,net}(V) = (q n_i / 2 \tau_\text{eff})\,(W(V) - W(0))$, plus the ideal Shockley saturation floor $J_s$.

Both sweeps use adaptive bias continuation (Newton step halving on divergence, re-growth on easy iterations) and start from the M1 equilibrium iterate.

## Notebook set

This is the second of four Colab walkthroughs that exercise the [end-of-M7 capability matrix](https://github.com/rwalkerlewis/SemiSim#status) from the README:

| Notebook | Covers |
|----------|--------|
| [01 — Equilibrium Poisson on a 1D pn junction](./01_pn_junction_1d.ipynb) | Nonlinear Poisson with Boltzmann carriers; depletion-approximation comparison |
| **02 — Bias sweep on the 1D pn junction** (this one) | Coupled Slotboom drift-diffusion with SRH; Shockley forward + SNS reverse |
| [03 — MOS capacitor C-V (2D multi-region)](./03_mos_cv.ipynb) | Gate contact BCs with $\phi_{ms}$; Si/SiO₂ multi-region; depletion-approximation C-V |
| [04 — 3D doped resistor V-I](./04_resistor_3d.ipynb) | Ohmic V-I linearity; builtin box mesh vs. gmsh `.msh` fixture |

**What this notebook does:**
1. Installs `dolfinx` on Colab via [FEM on Colab](https://fem-on-colab.github.io/) (one wget, ~30 s)
2. Clones the [SemiSim](https://github.com/rwalkerlewis/SemiSim) repo
3. Loads `benchmarks/pn_1d_bias/pn_junction_bias.json` and runs the forward sweep
4. Plots $|J(V)|$ against Shockley diffusion and SNS references on a log-linear axis
5. Loads `benchmarks/pn_1d_bias_reverse/pn_junction_bias_reverse.json` and runs the reverse sweep
6. Plots $|J(V)|$ against the SRH-generation reference and the $J_s$ floor

If you want to modify the physics (mesh resolution, doping, $\tau$, sweep range/step), edit the JSON and re-run. The Python code lives in the package — this notebook is a thin driver.

**Wall time on Colab.** Roughly 5-15 minutes total. The reverse sweep is the slower of the two because it has more voltage points and a tighter `min_step`.
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

Same pattern as Notebook 01 — `pip install -e` so any source edits you make persist in the Colab session."""))

cells.append(nbf.v4.new_code_cell(r"""import os
if not os.path.exists('SemiSim'):
    !git clone -q https://github.com/rwalkerlewis/SemiSim.git

%cd SemiSim
!pip install -q -e .
print("Package installed.")
"""))

cells.append(nbf.v4.new_code_cell(r"""import numpy as np
import matplotlib.pyplot as plt

from semi import schema, run, materials, constants
from semi.diode_analytical import (
    shockley_long_diode_saturation,
    sns_total_reference,
    srh_generation_reference,
)

print(f"semi version: {__import__('semi').__version__}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 3. Forward-bias sweep: 0 to 0.6 V

Symmetric 1D pn junction, 20 µm total, junction at 10 µm, $N_A = N_D = 10^{17}$ cm⁻³, $\tau_n = \tau_p = 10^{-8}$ s. Anode ramped in 0.05 V steps with cathode grounded. At each bias the runner assembles a 3-block $(\psi, \Phi_n, \Phi_p)$ Jacobian, solves with PETSc SNES, then integrates the current density at the anode facet to populate `result.iv`."""))

cells.append(nbf.v4.new_code_cell(r"""cfg_fwd = schema.load("benchmarks/pn_1d_bias/pn_junction_bias.json")

sweep = cfg_fwd["contacts"][0]["voltage_sweep"]
print(f"Case:       {cfg_fwd['name']}")
print(f"Dimension:  {cfg_fwd['dimension']}D, {cfg_fwd['mesh']['resolution'][0]} elements")
print(f"Doping:     N_A={cfg_fwd['doping'][0]['profile']['N_A_left']:.0e} cm^-3 | "
      f"N_D={cfg_fwd['doping'][0]['profile']['N_D_right']:.0e} cm^-3")
print(f"Sweep:      V_{cfg_fwd['contacts'][0]['name']} from {sweep['start']} V "
      f"to {sweep['stop']} V in {sweep['step']} V steps")
print(f"Solver:     {cfg_fwd['solver']['type']} with adaptive continuation")
"""))

cells.append(nbf.v4.new_code_cell(r"""result_fwd = run.run(cfg_fwd)

info = result_fwd.solver_info
print(f"\nBias sweep: {len(result_fwd.iv)} accepted points; "
      f"final V = {result_fwd.iv[-1]['V']:.3f} V")
print(f"Final SNES: iterations = {info.get('iterations', 'n/a')}, "
      f"converged = {info.get('converged', 'n/a')}, "
      f"reason = {info.get('reason', 'n/a')}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""### Forward-bias $J(V)$ vs analytical references

Two reference curves overlay the simulated data:
- **Shockley diffusion** (black dashed): ideal long-diode $J_s (e^{V/V_t} - 1)$. Matches the simulation near $V = 0.6$ V where diffusion dominates.
- **SNS (Sah-Noyce-Shockley)** (red dotted): $J_\text{diff} + J_\text{rec}$ including the depletion-region recombination term. Matches the simulation in the low-bias regime where SRH recombination dominates."""))

cells.append(nbf.v4.new_code_cell(r"""sc = result_fwd.scaling
mat = materials.get_material(cfg_fwd["regions"]["silicon"]["material"])
V_t = sc.V0

dop = cfg_fwd["doping"][0]["profile"]
N_A = dop["N_A_left"] * 1e6     # cm^-3 -> m^-3
N_D = dop["N_D_right"] * 1e6
mu_n_SI = cfg_fwd["physics"]["mobility"]["mu_n"] * 1.0e-4   # cm^2/Vs -> m^2/Vs
mu_p_SI = cfg_fwd["physics"]["mobility"]["mu_p"] * 1.0e-4
tau_n = cfg_fwd["physics"]["recombination"]["tau_n"]
tau_p = cfg_fwd["physics"]["recombination"]["tau_p"]
n_i = mat.n_i
eps = mat.epsilon

V_fwd = np.array([r["V"] for r in result_fwd.iv])
J_fwd = np.abs(np.array([r["J"] for r in result_fwd.iv]))

J_s, L_n, L_p = shockley_long_diode_saturation(
    N_A, N_D, n_i, mu_n_SI, mu_p_SI, tau_n, tau_p, V_t=V_t,
)
J_shockley = J_s * (np.exp(V_fwd / V_t) - 1.0)
J_sns, _Jd, _Jr, _ = sns_total_reference(
    N_A, N_D, n_i, eps, mu_n_SI, mu_p_SI, tau_n, tau_p, V_t=V_t, V=V_fwd,
)

fig, ax = plt.subplots(figsize=(7, 4.5))
fwd_mask = V_fwd > 0.05
ax.semilogy(V_fwd[fwd_mask], J_fwd[fwd_mask], "o-", color="C0",
            lw=2, ms=5, label="simulation (SemiSim)")
ax.semilogy(V_fwd[fwd_mask], np.maximum(J_shockley[fwd_mask], 1e-30),
            "k--", lw=1.5, label="Shockley diffusion $J_s(e^{V/V_t}-1)$")
ax.semilogy(V_fwd[fwd_mask], np.maximum(J_sns[fwd_mask], 1e-30),
            "r:", lw=2, label="SNS $J_\\mathrm{diff}+J_\\mathrm{rec}$")
ax.set_xlabel("V (V)")
ax.set_ylabel("|J| (A/m$^2$)")
ax.set_title("Forward-bias J(V) on the 1D pn junction")
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="lower right")
fig.tight_layout()
plt.show()

print(f"\n  J_s (Shockley long-diode): {J_s:.3e} A/m^2")
print(f"  L_n = {L_n*1e6:.2f} um, L_p = {L_p*1e6:.2f} um")
print(f"  Forward-sweep points: {len(V_fwd)}; V in [{V_fwd.min():.3f}, {V_fwd.max():.3f}] V")

row06 = min(result_fwd.iv, key=lambda r: abs(r["V"] - 0.6))
if abs(row06["V"] - 0.6) < 0.03:
    J_sim_06 = abs(float(row06["J"]))
    J_sh_06 = float(J_s * (np.exp(row06["V"] / V_t) - 1.0))
    err06 = abs(J_sim_06 - J_sh_06) / J_sh_06
    print(f"\n  At V={row06['V']:.2f} V: |J_sim| = {J_sim_06:.3e}, "
          f"Shockley = {J_sh_06:.3e}, rel err = {err06*100:.1f}%")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## 4. Reverse-bias sweep: 0 to -2 V

Same device, now sweeping the anode negative. The relevant physics shifts from diffusion/recombination to depletion-region **generation**. At $V = 0$ generation and recombination balance; under reverse bias the depletion region grows from $W(0)$ to $W(V)$ and the extra volume contributes a net generation current

$$J_\text{gen,net}(V) \;=\; \frac{q\,n_i}{2\,\tau_\text{eff}}\,\bigl(W(V) - W(0)\bigr),\qquad \tau_\text{eff} = \sqrt{\tau_n \tau_p}.$$

The ideal Shockley diffusion saturation $J_s$ is orders of magnitude smaller for $\tau \sim 10^{-8}$ s and shows up as the dashed floor."""))

cells.append(nbf.v4.new_code_cell(r"""cfg_rev = schema.load("benchmarks/pn_1d_bias_reverse/pn_junction_bias_reverse.json")

sweep = cfg_rev["contacts"][0]["voltage_sweep"]
print(f"Case:       {cfg_rev['name']}")
print(f"Sweep:      V_{cfg_rev['contacts'][0]['name']} from {sweep['start']} V "
      f"to {sweep['stop']} V in {sweep['step']} V steps")
print(f"Continuation: min_step={cfg_rev['solver']['continuation']['min_step']}, "
      f"max_halvings={cfg_rev['solver']['continuation']['max_halvings']}")
"""))

cells.append(nbf.v4.new_code_cell(r"""result_rev = run.run(cfg_rev)

info = result_rev.solver_info
print(f"\nBias sweep: {len(result_rev.iv)} accepted points; "
      f"final V = {result_rev.iv[-1]['V']:.3f} V")
print(f"Final SNES: iterations = {info.get('iterations', 'n/a')}, "
      f"converged = {info.get('converged', 'n/a')}, "
      f"reason = {info.get('reason', 'n/a')}")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""### Reverse-bias $|J(V)|$ vs SRH-generation reference

The reference overlay is $J_s + |J_\text{gen,net}(V)|$. The $J_s$ floor (black dashed) is the ideal diffusion saturation, the dotted curve adds the depletion-generation growth. The simulation should sit on the dotted curve across most of the saturation window."""))

cells.append(nbf.v4.new_code_cell(r"""V_rev = np.array([r["V"] for r in result_rev.iv])
J_rev = np.abs(np.array([r["J"] for r in result_rev.iv]))

J_gen_net, W0, V_bi = srh_generation_reference(
    N_A, N_D, n_i, eps, tau_n, tau_p, V_t=V_t, V=V_rev,
)
J_total_ref = np.abs(J_gen_net) + J_s

fig, ax = plt.subplots(figsize=(7, 4.5))
rev_mask = V_rev < -0.05
ax.semilogy(V_rev[rev_mask], np.maximum(J_rev[rev_mask], 1e-30), "o-",
            color="C3", lw=2, ms=5, label="simulation (SemiSim)")
ax.semilogy(V_rev[rev_mask], np.maximum(J_total_ref[rev_mask], 1e-30), "r:",
            lw=2, label=r"$J_s + |J_{\mathrm{gen,net}}(V)|$")
ax.axhline(J_s, color="k", ls="--", lw=1.2, label=f"$J_s$ = {J_s:.2e} A/m$^2$")
ax.set_xlabel("V (V)")
ax.set_ylabel("|J| (A/m$^2$)")
ax.set_title("Reverse-bias J(V) with SRH generation")
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="lower left")
fig.tight_layout()
plt.show()

print(f"\n  V_bi = {V_bi:.4f} V, W(0) = {W0*1e9:.1f} nm")
print(f"  tau_eff = sqrt(tau_n tau_p) = {(tau_n*tau_p)**0.5:.2e} s")

sat = [(r["V"], abs(r["J"])) for r in result_rev.iv if -2.0 <= r["V"] <= -0.5]
if sat:
    V_sat = np.array([v for v, _ in sat])
    J_sat = np.array([j for _, j in sat])
    J_ref_sat, _, _ = srh_generation_reference(
        N_A, N_D, n_i, eps, tau_n, tau_p, V_t=V_t, V=V_sat,
    )
    J_ref_sat = np.abs(J_ref_sat) + J_s
    rel_err = np.abs(J_sat - J_ref_sat) / J_ref_sat
    worst = int(np.argmax(rel_err))
    print(f"  Saturation window [-2, -0.5] V: {len(V_sat)} points")
    print(f"  Worst rel err vs reference: {rel_err[worst]*100:.1f}% at V={V_sat[worst]:.2f} V")
"""))

cells.append(nbf.v4.new_markdown_cell(r"""## Summary

- Forward sweep (0 to 0.6 V) reproduces both the SNS recombination branch at low bias and the Shockley diffusion branch at high bias. The verifier (`scripts/run_benchmark.py pn_1d_bias`) asserts ≤15% vs SNS on [0.15, 0.55] V and ≤10% vs Shockley diffusion at $V = 0.6$ V.
- Reverse sweep (0 to -2 V) reproduces SRH-generation saturation with the correct $W(V)$ scaling. The verifier (`pn_1d_bias_reverse`) asserts ≤20% vs the generation reference on [-2, -0.5] V.
- Both sweeps solve with adaptive Newton continuation — if the solver struggles at a bias step it halves and re-tries; if it converges easily it grows back up.
- The current-continuity check (interior sampling of $J_\text{total}$ inside the device) holds to within 5% on the forward branch and 15% on the reverse branch at the final bias; these tolerances are built into the benchmark V&V.

**Next up.** Notebook 03 moves from 1D pn to a 2D MOS capacitor, adds a gate contact with metal-semiconductor work-function difference $\phi_{ms}$, and runs a C-V sweep against the depletion-approximation theory curve.

See [the repo](https://github.com/rwalkerlewis/SemiSim) and its [capability matrix](https://github.com/rwalkerlewis/SemiSim#status) for the full scope."""))

nb.cells = cells

out_path = Path(__file__).resolve().parent.parent / "notebooks" / "02_pn_junction_bias.ipynb"
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as f:
    nbf.write(nb, f)
print(f"Wrote {out_path} with {len(cells)} cells")
