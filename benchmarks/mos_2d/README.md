# 2D MOS capacitor benchmark (M6, M14.1)

First 2D benchmark and first multi-region device in SemiSim.

As of M14.1 (2026-04-26) the C-V curve is extracted using the
analytic-AC differential capacitance runner
(`solver.type == "mos_cap_ac"`) rather than the original
finite-difference dQ/dV path (`mos_cv`). See "C-V extraction" below.

## Device

- Silicon substrate: 500 nm x 500 nm, uniform N_A = 1e17 cm^-3 (p-type)
- SiO2 gate oxide: 500 nm x 5 nm on top
- Body contact: ohmic, y = 0
- Gate contact: ideal-gate Dirichlet (phi_ms = 0), y = 505 nm
- Left / right edges: natural (Neumann) boundaries

## Mesh

Uniform 2D rectangle, triangles: 4 cells laterally, 505 cells
vertically (1 nm per cell). Interface at y = 500 nm is a grid line;
the oxide gets exactly 5 cells.

Total DOFs: 2530 (P1 Lagrange on 4 x 505 rectangle -> 5 x 506 nodes).

## Physics

Multi-region equilibrium Poisson (`solver.type == "mos_cap_ac"`, M14.1):

- Full-mesh stiffness `L_D^2 * eps_r(x) * grad psi . grad v` with
  cellwise DG0 `eps_r` (Si: 11.7, SiO2: 3.9)
- Boltzmann space-charge `n_i (exp(-psi) - exp(psi)) + N_net`
  restricted to silicon via `dx(subdomain_id=semi_tag)`
- Oxide: Laplacian only (no space charge, no Slotboom variables)

At each V_gate the gate charge per area is

    Q_gate(V_gate) = -(q / W_lat) * integral_{Omega_Si} rho dA

and the differential capacitance is extracted analytically per bias
point (see "C-V extraction" below).

## C-V extraction (M14.1)

The legacy `mos_cv` runner records `Q_gate(V_gate)` and the verifier
extracts `C_sim` by `numpy.gradient(Q, V)`. That centered finite
difference is the omega -> 0 numerical limit of admittance and works,
but it picks up two avoidable error sources:

* O(h^2) finite-difference error from the bias step h plus residual
  jitter from per-step SNES termination.
* One-sided differences at sweep endpoints, which forced the verifier
  to drop the first and last bias point from its window.

The M14.1 `mos_cap_ac` runner solves one extra linear system per bias
point to get the same number directly. After convergence at V_gate it
solves

    K(psi_0) * delta_psi_hat = 0    with delta_psi_hat(gate) = 1/V_t

where K is the Jacobian of the multi-region Poisson residual at the
converged psi_0 and the gate-Dirichlet perturbation 1/V_t encodes the
unit V_gate sensitivity. The differential capacitance is then

    C(V_gate) = (q C_0 / W_lat) integral_{Omega_Si}
                ni_hat (exp(-psi_0) + exp(psi_0)) * delta_psi_hat dx

This is exact to discretisation, free of FD step error, and uses every
bias point in the verifier window (no endpoint loss). It matches the
spirit of the M14 `ac_sweep` runner at omega = 0 but is built around
the multi-region equilibrium Poisson stack, which `ac_sweep` does not
support today (it targets ohmic two-terminal pn diodes; see
docs/adr/0011-ac-small-signal.md for the locked M14 scope).

The legacy `mos_cv` solver type is retained for backwards compatibility
but deprecated. Both paths agree on `Q_gate(V_gate)` byte-identically;
the verifier and plotter auto-detect `iv[i]["C_ac"]` and prefer it
when present.

### Headline numbers vs the FD-based path

|              | worst rel err | window samples | endpoint usable |
| ------------ | ------------- | -------------- | --------------- |
| `mos_cv` FD  | 9.25 %        | 14 / 43        | no              |
| `mos_cap_ac` | 6.79 %        | 16 / 43        | yes             |

Both at V_gate = -0.200 V (depletion onset, where the depletion
approximation itself starts losing fidelity); both inside the fixed
10 % tolerance.

## Verification target

Depletion-approximation C-V curve within 10% in the window
[V_FB + 0.2, V_T - 0.1] V. For the ideal-gate M6 device under our
psi=0-at-intrinsic BC convention:

    V_FB = phi_ms - phi_F = 0 - 0.417 V = -0.417 V
    V_T  = V_FB + 2 phi_F + sqrt(4 eps_s q N_A phi_F) / C_ox = +0.658 V
    window = [-0.217, +0.558] V

Worst observed error in the window with the M14.1 analytic-AC method:
6.79% at V_gate = -0.200 V (down from 9.25% with the legacy FD path).

Monotone non-increasing C_sim across the window is also checked (C
drops from near C_ox in accumulation through a minimum at strong
inversion onset).

## Excluded regimes

- **Accumulation** (V_gate < V_FB): holes pile up at the interface;
  the depletion approximation does not model the accumulation
  surface charge. C_sim recovers toward C_ox but is not gated here.
- **Strong inversion** (V_gate > V_T): minority-carrier inversion
  layer forms. Quasi-static C-V recovery toward C_ox is not well
  modelled without a transient solver. Not gated.

The full V_gate sweep [-0.9, +1.2] V is rendered on `cv.png` for
reference (simulation line crosses through both excluded regions);
the verifier window is shaded.

## Running

```
docker compose run --rm benchmark mos_2d
```

Produces four PNGs in `results/mos_2d/`:
- `psi_2d.png`: tricontourf of psi(x, y) at the last V_gate
- `potentials_1d.png`: central-column psi(y) slice
- `cv.png`: C_sim vs depletion-approximation theory
- `qv.png`: Q_gate vs V_gate

## References

- `docs/mos_derivation.md`: full derivation gate (device, interface
  conditions, submesh formulation, gate BC, C-V theory, MMS
  construction).
- `docs/PHYSICS.md` Section 6: condensed physics reference including
  the V_FB shift from the BC convention.
