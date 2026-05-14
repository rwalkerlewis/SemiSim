# benchmarks/moscap_axisym_2d

Axisymmetric (cylindrical) 2D MOS capacitor reproducing Hu chapter 5
Fig. 5-18: low-frequency (quasi-static) and high-frequency C-V curves
on the same axes.

See also:
[`docs/benchmarks/moscap_axisym_2d.md`](../../docs/benchmarks/moscap_axisym_2d.md)
(landing page),
[`docs/theory/axisymmetric.md`](../../docs/theory/axisymmetric.md)
(weak-form derivation),
[`docs/theory/moscap_cv.md`](../../docs/theory/moscap_cv.md)
(LF/HF C-V derivation), and
[`notebooks/05_moscap_axisym_cv.ipynb`](../../notebooks/05_moscap_axisym_cv.ipynb)
(end-to-end demonstration).

## Files

- `moscap_axisym.geo`: gmsh geometry for the meridian half-plane
  `(r, z) in [0, R] x [-T_Si, T_ox]`. Builds a graded mesh with fine
  triangles at the Si/SiO2 interface and the gate edge `r = R_g`.
- `moscap_axisym.json`: SemiSim input config (schema 1.3.0,
  `coordinate_system: "axisymmetric"`).
- `reference_cv.csv`: analytical reference C-V from Hu chapter 5
  (closed form), used by the regression test as the ground truth.

The compiled `.msh` is `.gitignore`d. Regenerate it from the `.geo`
before running the benchmark:

```bash
gmsh -2 -format msh22 \
  -o benchmarks/moscap_axisym_2d/moscap_axisym.msh \
  benchmarks/moscap_axisym_2d/moscap_axisym.geo
```

## Physics

Cylindrical-symmetry weak forms multiply every volume integrand by `r`
(the radial coordinate). The implementation lives in
`semi/physics/axisymmetric.py`. The symmetry axis `r = 0` is a natural
no-flux boundary (no Dirichlet allowed; the schema validator rejects
configs that try). The outer wall `r = R = 200 um` is homogeneous
Neumann. R is chosen `>> 5 W_dmax` so the solution near the axis is
insensitive to the cutoff; the convergence study in
`notebooks/05_moscap_axisym_cv.ipynb` quantifies this.

## C-V extraction

Two curves are computed by `semi/cv.py`:

- `compute_lf_cv_fem`: `C_LF = -dQ_s/dV_g` via centered finite
  difference of the total semiconductor charge over the gate area.
  Hu Eq. 5.6.1.
- `compute_hf_cv_depletion_clamp`: high-frequency C from the FEM
  surface potential, with the depletion width `W_dep` clamped at
  `W_dmax` once `phi_s >= 2 phi_B`. Documented as the chosen HF
  method for this benchmark; a fully principled small-signal solve
  with frozen minority charge is a future extension.

## Reference parameters

| Parameter | Value |
|---|---|
| `N_a`     | 5e16 cm^-3 (p-type Si) |
| `T_ox`    | 10 nm SiO2 |
| `phi_ms`  | -0.95 V (N+ poly gate) |
| `Q_f`     | 0 |
| `T`       | 300 K |
| `R_g`     | 50 um |
| `R`       | 200 um |
| `T_Si`    | 5 um |
| Gate sweep | -2.0 V to +2.0 V, 81 points |

Analytical anchors (also encoded in the CSV header lines):

| Quantity | Value |
|---|---|
| `V_fb`     | -0.950 V |
| `V_t`      | +0.181 V |
| `phi_B`    | 0.399 V |
| `W_dmax`   | 144 nm |
| `C_ox`     | 3.453e-3 F/m^2 |
| `C_min`    | 5.97e-4 F/m^2 (= 0.173 C_ox) |
