# ADR 0012: Scharfetter-Gummel edge-flux discretisation for transient continuity

**Date:** 2026-04-26
**Status:** Accepted (supersedes the (n, p) Galerkin choice in ADR 0009); primary-unknown assumption amended by ADR 0014; the SG block-assembly module (`semi/fem/sg_assembly.py`, ~792 LOC) was deleted in M14.3 (PR on `dev/m14.3-housekeeping`) since it has no production caller; the per-edge primitives (`semi/fem/scharfetter_gummel.py`: Bernoulli function and edge fluxes) are retained for future 2D high-Peclet work. If M13.2 ever revives the SG flux path, the deleted block-assembly is reachable in git history at the parent of the M14.3 deletion commit.
**Milestone:** M13.1 — Transient solver close-out

---

## Context

ADR 0009 selected the (psi, n_hat, p_hat) primary-density formulation for the transient continuity equations and discretised the convection-diffusion blocks with **standard P1 Galerkin** (the same form builders the steady-state Slotboom runner uses, with `n_hat` substituted for `n_i exp(psi - phi_n)`).

That combination ships in M13 (PR #30) and fails. The merged transient runner reproduces issue #34's "Bug #2": at tight SNES atol the Newton iterates blow up to ~1e24 m^-3 carrier densities, the SRH denominator collapses, and the iterate locks at a non-physical fixed point. Mesh refinement does not change the failure mode; loosening atol does, by stopping Newton before the blow-up. The M13 close-out (PR #35) xfail'd `tests/fem/test_transient_steady_state.py` with a pointer to this ADR.

The structural problem is the same one Scharfetter and Gummel diagnosed in 1969: standard Galerkin (or central-difference FD) on the convection-diffusion drift term `mu_n n grad psi` is unstable when the cell Peclet number `|grad psi * h| / (2 V_t)` exceeds 1, which is almost everywhere in a real device. The 1D Read-diode paper (Scharfetter & Gummel 1969, IEEE TED ED-16(1):64-77) absorbs the exponential carrier-density variation along an edge into the basis function itself. The result is unconditionally stable in the cell Peclet number and recovers the correct upwind limit at high field.

The 2D extension is non-trivial: the natural construction lives on the edges of simplicial elements with Nedelec H(curl) edge basis functions, and the per-edge SG flux scatters into two stiffness-matrix entries per edge. Brezzi, Marini, and Pietra (1989, SIAM J. Numer. Anal. 26:1342-1355) gave the first rigorous treatment; Bochev et al. (Sandia 2011-3865) wrote the open implementation note that this ADR draws from for the multi-dimensional construction.

The earlier sibling-branch attempt (commit `323d30a` on `dev/m13-mesh-and-slotboom-test`) tried to fix M13 by switching the unknowns to Slotboom (psi, phi_n, phi_p) with a chain-rule time term. That formulation guarantees positivity of `n` and `p` by construction (they are `n_i exp(...)`) but introduces a Jacobian that is ill-conditioned by ~25 orders of magnitude across the device, and MUMPS LU fails with `DIVERGED_LINEAR_SOLVE`. The attempt was abandoned with a note that the textbook fix is SG. This ADR is that fix.

## Decision

Discretise the transient continuity convection-diffusion operator in `semi/runners/transient.py` with **Scharfetter-Gummel edge fluxes**, while keeping the (psi, n_hat, p_hat) primary unknowns from ADR 0009 unchanged. Carrier densities remain the primary unknowns (no chain-rule time term, no Slotboom transform). Newton's positivity-of-(n, p) failure mode is addressed by the SG flux, which is the upwind scheme that Galerkin failed to be: it does not amplify negative excursions into a divergent fixed point because the edge flux is a convex combination weighted by `B(+- dpsi)` rather than a centred difference. SRH source/sink terms remain as quadrature on cell vertices (lumped-mass) and are unchanged.

### Bernoulli function

```
B(x) = x / (exp(x) - 1)        for x != 0
B(0) = 1                       (Taylor series, exact)
```

Identity: `B(x) - B(-x) = -x` (algebraic; numerator `2 cosh(x) - 2`, denominator `-(2 cosh(x) - 2)`, ratio `-1`). The earlier draft of this prompt stated `B(x) + B(-x) = -x`, which is wrong; the sum is `x coth(x/2)`. The corrected identity is the one used in the unit tests.

Stable evaluation in four regimes:

- `x == 0`: return 1.0 exactly.
- `|x| < eps_taylor` (default 1e-3): degree-4 Taylor, `1 - x/2 + x*x/12 - x**4/720`.
- `x >> 0` (default `x > 30`): rewrite as `x exp(-x) / (1 - exp(-x))`, where `1 - exp(-x)` is bounded away from zero and `x exp(-x)` underflows benignly.
- `x << 0` (default `x < -30`): rewrite as `x / (exp(x) - 1)` directly is fine because `exp(x) -> 0` and the denominator `-> -1`, so `B(x) -> -x` as `x -> -inf`.
- Mid-range `1e-3 <= |x| <= 30`: closed form `x / (exp(x) - 1)`.

### 1D edge flux

For an edge connecting nodes `i, j` with edge length `h_ij` and scaled potential difference `dpsi = psi_hat_j - psi_hat_i`, the SG electron flux on the edge, in the **Sandia / Farrell-et-al convention** (the same convention the existing `evaluate_partial_currents` resolves to once the SemiSim `phi_n` sign convention is unwound):

```
F_n_ij = (mu_n * V_t / h_ij) * ( n_j * B(+dpsi) - n_i * B(-dpsi) )
```

`F_n_ij` is the electron contribution to the conventional current density on the edge from i to j. At `dpsi = 0` with `n_j > n_i`, `F > 0` (diffusion in the `+grad(n)` direction). At uniform `n` with `dpsi > 0`, `F < 0` (electrons drift toward j; conventional current opposite).

The hole flux in the same convention:

```
F_p_ij = (mu_p * V_t / h_ij) * ( p_i * B(+dpsi) - p_j * B(-dpsi) )
```

The (i, j) positions in the hole expression are swapped relative to the electron expression. Starting from the electron form `[n_j B(+dpsi) - n_i B(-dpsi)]`, the device-equation symmetry `(n <-> p, psi -> -psi)` yields, after re-grouping by the i, j positions of the swapped flux, the hole spelling above. Cross-checked against Spevak Section 3.2.6 (which uses the opposite *electron-particle-flux* convention `B_S(x) = -B(-x)`): Spevak's `phi_p` equals `-F_p_ij` here, sign-flip exactly as expected. Guard 2b verifies the symmetry numerically to 1e-10 relative.

### Sign convention (Sandia / Farrell-et-al)

The SG flux defined above represents the **electron contribution to the conventional current density** along the edge from i to j (similarly for holes). This is the same convention used in:

- arXiv 1911.00377 Eq (30) (Farrell, Doan, Hopf, Koprucki, Glitzky 2019).
- Sandia OSTI 2011-3865 Eq (18) (Bochev et al. 2011), after substituting `a_ij = -(psi_j - psi_i) / (2 beta)`.

Spevak's TU Wien thesis Section 3.2.6 uses the **opposite** convention (electron *particle* flux, in the direction of electron motion). The relation `B_Spevak(x) = -B(-x)` maps between the two; the underlying physics is identical, just the sign of the reported flux flips. Future readers should not interpret the Spevak sign as a disagreement with Sandia/arxiv on physics.

### Consistency with `evaluate_partial_currents`

The existing `semi/postprocess.py::evaluate_partial_currents` computes
`J_n = q mu_n n grad(phi_n) . n_outward` (Slotboom form). Unwinding the
SemiSim `phi_n` sign convention (the BC sets `phi_n = +V_applied` at
ohmic contacts, which is sign-flipped from the textbook
`phi_n = -V_applied`):

- `phi_n_semisim = -phi_n_textbook`
- `grad(phi_n_semisim) = -grad(phi_n_textbook)`
- `J_n_semisim = +q mu_n n grad(phi_n_semisim) = -q mu_n n grad(phi_n_textbook) = J_n_textbook`

So `evaluate_partial_currents` returns the textbook conventional electron current density. At zero field with `n_j > n_i` along the outward direction, this is `+q D_n grad(n) > 0`. The SG flux above also gives `> 0` at the same configuration. They agree on sign. M13.1 changes the IV recorder in `transient.py` to use the SG edge fluxes directly instead of recomputing from `(n, psi)` via the Slotboom intermediate; the terminal-current sign at ohmic contacts remains the M12 / ADR 0008 sign convention end-to-end.

### 2D simplicial edge-flux assembly (lumped-mass approximation)

For each simplicial element, iterate over its edges. For each edge `e_ij` of length `h_ij` connecting nodes `i, j`:

1. Compute `dpsi = psi_hat_j - psi_hat_i` from the per-element vertex DOFs.
2. Compute the SG flux scalar `F_n_ij` per the 1D formula above.
3. Scale by the dual-side length `|partial C_ij|`. For a P1 simplicial mesh in the lumped-mass approximation this equals `area(K) / 3` in 2D for each of the three element edges (area-share of the dual control volume around vertex i contributed by element K, attributed equally to each of the two edges incident at i within K). In 3D the analogous fraction is `volume(K) / 4` per edge per tet, but SemiSim does not target 3D transient before M16 so this is documented for future work, not implemented.
4. Scatter `F_n_ij * |partial C_ij|` into the residual entry for vertex `i` with sign `-` (efflux from i) and into vertex `j` with sign `+` (influx into j). Hole flux mirror this.

**This is an approximation.** The full Brezzi-Marini-Pietra / Bochev construction weights the per-edge SG flux by the lowest-order Nedelec H(curl) edge basis function projected onto the dual control-volume side; for P1 simplicial elements with the lumped-mass / nearest-neighbour control volume the weight collapses to the dual-side length `|partial C_ij|` exactly only for isotropic, well-shaped elements. For higher-order elements, anisotropic meshes, or curved control volumes the full Nedelec weighting (Sandia 2011-3865 Eqs 22-29) is required. M13.1 ships the lumped approximation; if a future benchmark exposes an anisotropic-mesh failure the fix is to swap in the full edge-element weighting without changing the per-edge flux scalar. This deferred work is captured in the **Forward notes** section below.

### Issue #34 "Bug #3" sign-convention fix

Issue #34's "Bug #3" notes that the existing `evaluate_partial_currents` reconstructs `phi_n` from `(n, psi)` and feeds the Slotboom flux formula, which has opposite sign from the natural `(n, p)` flux. The M13.1 fix is to switch the transient runner's IV recorder to call directly into the SG edge-flux assembler using the `(n, p)` primary unknowns, with the Sandia / Farrell convention above. Because the SG convention agrees with `evaluate_partial_currents`'s sign (see "Consistency with `evaluate_partial_currents`" above), the post-merge artefact is one terminal-current sign convention used end-to-end in the transient runner and the IV recorder; the steady-state runner (`bias_sweep.py`) keeps the existing Slotboom path and same sign.

## Consequences

**Positive:**

- The (n, p) blow-up failure mode of M13 / ADR 0009's Galerkin choice is structurally eliminated. The SG flux is upwind by construction; carrier densities cannot diverge to 1e24 m^-3 because the per-edge flux is bounded by `min(n_i, n_j) * exp(|dpsi|) * mu_n V_t / h` regardless of central-difference oscillations.
- (psi, n, p) primary unknowns retained, so the lumped mass matrix from M13 / ADR 0010 carries over unchanged, the BDF1/BDF2 time-stepper is unchanged, and the Newton outer loop is unchanged.
- The Slotboom-chain-rule attempt (commit 323d30a)'s 25-OOM Jacobian conditioning failure is avoided because (n, p) primary unknowns make the mass matrix linear and well-conditioned.
- The reuse path for steady-state DD is open: the SG edge-flux assembler can replace the Galerkin convection blocks in `bias_sweep.py` in a future PR. M13.1 only does it for the transient runner to keep the change scoped, but no architectural barrier remains.

**Negative:**

- Per-edge assembly is more code than UFL stiffness assembly. The SG module is pure-Python NumPy / dolfinx low-level, not UFL. Coverage on the new module needs to be >= 95% in its own right (per the project's Invariant 4), and the load-bearing piece is the Bernoulli function (~10 unit tests).
- The lumped-mass weight collapse is a real approximation and is documented as such. Anisotropic meshes will eventually need the full Nedelec weighting.

**Neutral:**

- Recombination R, Poisson assembly, BCs, and IV evaluation are all unchanged in expression; only the convection-diffusion blocks of the residual change.
- Steady-state runner (`run_bias_sweep`) is **unchanged** in this PR; it stays on Slotboom + Galerkin. SG can be ported there later without any ADR change.

## Test plan (commitments)

These are the gates `tests/fem/test_transient_steady_state.py` and `tests/mms/test_transient_convergence.py` will need to clear, in addition to the new SG-specific tests under `tests/fem/test_scharfetter_gummel.py`:

1. **Bernoulli mpmath cross-check.** 1000 log-uniform-in-`|x|` points across `[-50, 50]`, mixed signs, 1e-12 relative tolerance vs `mpmath.mp.mpf` evaluation of `x / (exp(x) - 1)` at 50 decimal digits (with the `B(0) = 1` Taylor branch for `|x| < 1e-300`). The mpmath result is the independent reference; the `numpy` direct evaluation is excluded for `|x| > 30` because of overflow.

2. **Two-direction sign test, electrons.** Edge with `n_i = 1e16`, `n_j = 1e17`, `psi_i = 0`, `psi_j = +0.1 V`. SG flux must agree with the **midpoint-Galerkin reference** `mu_n V_t (dn/h - n_mid * dpsi/(V_t h))`, `n_mid = (n_i + n_j)/2`, to 5% relative at this small cell Peclet number (`|dpsi|/(2) ~ 0.05 << 1` in scaled units). Sign error: orders of magnitude. Factor-of-2 error: 2x. Agreement to 5%: sign and scaling correct. Then swap `i` and `j`, recompute, confirm flux flips sign with the same magnitude. The Slotboom 2-node closed form is **not** used here because it derives from the same `B(x)` and only catches typos.

3. **Two-direction sign test, holes (guard 2b).** Same edge but with `p_i = 1e17`, `p_j = 1e16`, `psi_i = 0`, `psi_j = +0.1 V` (field drives holes from `j` toward `i`, opposite of electrons). The hole flux must satisfy the device-equation symmetry under simultaneous `(n <-> p, psi -> -psi, anode <-> cathode)`: it must equal the negative of the corresponding electron flux configuration to 1e-10 relative. The hole sign is the most error-prone part of SG and this guard is non-negotiable.

4. **1D steady-state agreement gate, 1e-4 relative, hard.** Promoted from a "soft signal" to the only hard merge gate. `transient` ramped to steady state at V = 0.5 V on the `pn_1d_bias` config must agree with `bias_sweep` terminal current to 1e-4 relative at every mesh refinement N in `{100, 200, 400, 1000}`. Coarse-mesh trap detection: error must be **monotone decreasing** in N. If error at N = 400 is lower than at N = 1000 the PR does not merge regardless of MMS rates.

5. **MMS rates, soft signals.** BDF1 rate floor 0.95, BDF2 rate floor 1.90 on the existing pairwise-difference test. These are sanity checks; the steady-state agreement gate is the real one.

6. **xfail flip.** `tests/fem/test_transient_steady_state.py` flips from xfail to passing without modification.

7. **Coverage.** The new `semi/fem/scharfetter_gummel.py` module is gated at 95% in its own right; in particular every Bernoulli regime has a unit test.

## Forward notes

- **Anisotropic meshes / higher-order elements.** The lumped-mass weight `|partial C_ij|` collapse from the full Brezzi-Marini-Pietra / Bochev construction is an approximation that is exact for isotropic, well-shaped P1 simplicial elements only. If a future benchmark on an anisotropic mesh fails the steady-state agreement gate, the fix is to swap in the full Nedelec H(curl) edge-element weighting from Sandia 2011-3865 Eqs 22-29 without changing the per-edge SG flux scalar. The interface is: replace the `dual_side_length(K, e)` helper with a call into a UFL-form `nedelec_edge_weight(K, e, ds)` that integrates the Nedelec basis against the dual-side normal.

- **3D extension.** The 1D and 2D paths land in M13.1. 3D simplicial transient is deferred to whenever 3D transient becomes a target (currently not on the roadmap; M15 is GPU linear-solver scaling, M16 is physics completeness). The 3D dual-side weight is `volume(K) / 4` per edge per tet under the lumped approximation.

- **SG in steady-state runner.** `bias_sweep.py` keeps Slotboom + Galerkin in M13.1. The SG edge-flux assembler is general enough to drop into the steady-state convection-diffusion blocks if a future benchmark exposes a Galerkin failure there too. No ADR change required for that move.

## Sources actually used

### Caught-by-cross-check note

The M13.1 prompt's first-draft SG formula was

```
F_n_ij_prompt = (mu_n * V_t / h_ij) * ( B(-dpsi) * n_j - B(+dpsi) * n_i )
```

with `dpsi = psi_j - psi_i`. The B-arguments were swapped relative to the standard convention. The midpoint-Galerkin cross-check guard (test plan #2) caught this before any 1D residual code was wired: at small Peclet (`dpsi=0.005V`, `Pe=0.097`) the user-prompt formula disagreed with the midpoint-Galerkin reference by 27%, all on a sign flip in the leading-order drift-term Taylor expansion. The correct formula

```
F_n_ij = (mu_n * V_t / h_ij) * ( n_j * B(+dpsi) - n_i * B(-dpsi) )
```

agrees with the midpoint-Galerkin reference to 0.35% at the same Peclet. Two open sources (arXiv 1911.00377 Eq 30 and Sandia OSTI 2011-3865 Eq 18) confirm the corrected form. Future-you will appreciate knowing this guard fired; the unit test `test_sg_edge_flux_n_agrees_with_midpoint_galerkin_at_small_peclet` is the same guard in regression form.

### Open sources read

Open sources read end-to-end during the M13.1 sourcing pass:

- **Farrell, Doan, Hopf, Koprucki, Glitzky, "Generalized Scharfetter-Gummel schemes for electro-thermal transport in degenerate semiconductors using the Kelvin formula for the Seebeck coefficient", arXiv:1911.00377v3, https://arxiv.org/pdf/1911.00377.** Source for the standard Bernoulli function `B(x) = x / (exp(x) - 1)` (Eq 30) and the 1D edge flux in conventional-current convention `J_{n,K,L} = M_n k_B T_{K,L} [n_L B(X) - n_K B(-X)]` with `X = q (phi_L - phi_K) / (k_B T_{K,L})`. Open arXiv preprint, full PDF text-extractable.

- **Spevak, Stefan, "Numerical Aspects of the Simulation of Electromigration", TU Wien doctoral thesis Section 3.2.6, https://www.iue.tuwien.ac.at/phd/spevak/node66.html.** Source for the alternate-convention electron flux `phi_n = U_th mu_n (n_i B_Spevak(-dpsi) - n_j B_Spevak(dpsi))` with `B_Spevak(x) = x / (exp(-x) - 1)` (sign-flipped from the standard `B`). Confirms the 1D form against the standard convention; under `B_Spevak(x) = -B(-x)` the two reduce to the same flux up to overall sign. Open thesis from Selberherr's TU Wien institute.

- **Bochev, Lehoucq et al., "Control Volume Finite Element Method with Multidimensional Edge Element Scharfetter-Gummel upwinding. Part 1. Formulation", Sandia Report SAND2011-3865, OSTI 1020517, https://www.osti.gov/servlets/purl/1020517/.** Source for the multi-dimensional simplicial edge-flux construction. Eq (18) gives the per-edge SG flux `J_ij = (D_n / |e_ij|) [n_j B(-2 a_ij) - n_i B(2 a_ij)]` with `a_ij = -(psi_j - psi_i) / (2 beta)`. Eqs (19), (22), (24)-(29) give the per-edge to per-vertex assembly and the Nedelec H(curl) edge-element weighting. Open Sandia release.

Cited but not directly read; their content is reproduced verbatim in the open sources above:

- **Scharfetter, D. L. and Gummel, H. K., "Large-signal analysis of a silicon Read diode oscillator", IEEE Trans. Electron Devices ED-16(1):64-77, 1969, doi:10.1109/T-ED.1969.16566.** Original 1D SG paper. The discretisation is reproduced exactly in arXiv 1911.00377 Eq 30 and Sandia 2011-3865 Eq 18.

- **Bank, R. E. and Rose, D. J., "Some error estimates for the box method", SIAM J. Numer. Anal. 24(4):777-787, 1987, JSTOR 2157588.** Box-method error analysis, paywalled. The SG flux it analyses is identical to Sandia 2011-3865 Eq 18 and a Springer chapter "On the Scharfetter-Gummel Box-Method" (located but not fully read) reproduces it.

- **Brezzi, F., Marini, L. D., and Pietra, P., "Two-dimensional exponential fitting and applications to drift-diffusion models", SIAM J. Numer. Anal. 26(6):1342-1355, 1989, doi:10.1137/0726078.** The 2D simplicial extension. Paywalled; the construction is reproduced and extended to 3D in Sandia 2011-3865.

- **Selberherr, S., "Analysis and Simulation of Semiconductor Devices", Springer, 1984, ch. 4.4.** Textbook reference. The box-method discretisation with `B(x) = x/(exp(x) - 1)` from this textbook is reproduced verbatim in the Spevak TU Wien thesis (same institution) and in dozens of open lecture notes.
