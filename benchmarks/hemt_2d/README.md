# hemt_2d: AlGaAs/GaAs HEMT 2DEG sheet density

Acceptance benchmark for M17 (heterojunction / position-dependent
band structure). Validates the cellwise DG0 chi/Eg/n_i fields built
by `semi/physics/heterojunction.py` and the local-chi ohmic
equilibrium psi from `semi/bcs.py` (Phase D) against a published
classical-electrostatic 2DEG reference within 15 % at every
`V_GS in [0.4, 1.0] V`.

## Device

Layered 2D AlGaAs / GaAs HEMT, top to bottom:

| Layer    | Material   | y range   | Doping           |
| -------- | ---------- | --------- | ---------------- |
| Barrier  | AlGaAs_0p3 | 300-330 nm | uniform N_D = 5e18 cm^-3 |
| Channel  | GaAs       | 200-300 nm | uniform N_A = 1e15 cm^-3 |
| Buffer   | GaAs       | 0-200 nm  | uniform N_A = 1e15 cm^-3 |

Lateral extent: 200 nm. Total thickness: 330 nm. Mesh: 4 x 132 cells
(uniform vertical resolution 2.5 nm; the heterojunction at y=300nm
sits on a grid line).

The barrier carries uniform donor doping at 5e18 cm^-3 (representing
a delta-doping plane spread across 30 nm of barrier; the integrated
areal density 5e18 cm^-3 * 30 nm = 1.5e13 cm^-2 lies in the standard
HEMT modulation-doping range). Future refinement may resolve the
delta plane explicitly via a 1 nm sub-region.

The `regions[barrier].heterojunction: true` flag activates the M17
position-dependent material parameter path: chi, Eg, Nc, Nv, n_i are
built per cell as DG0 fields, the Slotboom substitution becomes
`n = n_i(x) exp((psi - phi_n) / V_t)` (ADR 0016), and the back-
contact ohmic equilibrium psi reads chi from the local material
(Anderson-rule band alignment).

## Contacts

| Contact | Facet | Kind     | Notes                                 |
| ------- | ----- | -------- | ------------------------------------- |
| back    | y=0   | ohmic    | reference body contact                |
| gate    | y=330nm | Schottky | Pt-on-AlGaAs; phi_B = 0.95 eV (M16.5) |

The gate is swept V_GS in [0.0, 1.0] V (11 points). The Schottky path
is the M16.5 thermionic-emission Robin BC; ADR 0015 documents the
V&V scope.

## Physics

- `physics.statistics: "fermi_dirac"` (M16.4): the heavy modulation
  doping in the AlGaAs barrier requires Fermi-Dirac statistics for
  the channel-charge balance.
- `physics.mobility: {"model": "constant"}` (M0): the 2DEG sheet
  density depends on the electrostatic balance, not channel mobility.
  Future `caughey_thomas` (M16.1) variants are out of scope.

## Verifier and acceptance

`scripts/run_benchmark.py hemt_2d` invokes `verify_hemt_2d`, which:

1. Runs the V_GS sweep via `run_bias_sweep`.
2. At each V_GS, integrates `n(y) - N_D(y)` over y in
   [300 nm, 330 nm] (the barrier-side accumulation region) and over
   [200 nm, 300 nm] (the channel-side 2DEG). The signed sheet
   density `n_s(V_GS)` is the result.
3. Compares to the classical-electrostatic closed form
   `n_s(V_GS) = (eps_AlGaAs * eps_0 / (q * d_barrier)) * (V_GS - V_T)`
   in the linear regime above threshold (Pozela & Reklaitis 1980;
   Sze 3rd ed Section 7.4).
4. Asserts `|n_s_FEM - n_s_classical| / n_s_classical < 0.15` at
   every `V_GS in [0.4, 1.0] V`. The `V_GS = 0.0 - 0.3 V` range is
   excluded because the 2DEG is sub-threshold and the relative error
   blows up.

## Why 15 %

The 15 % tolerance reflects the gap between the *classical
electrostatic* solver (which kronos-semi is) and the *fully quantum*
Poisson-Schrodinger references usually cited for HEMT 2DEG sheet
densities (Stern 1972, Stern & Sarma 1984, Pozela & Reklaitis 1980).
The classical solver matches the classical-electrostatic reference
within ~5 %; the quantum reference is ~5-15 % below the classical
because quantum confinement penalizes the band-edge accumulation.
The 15 % budget covers both.

A proper quantum-confined HEMT verifier (Schrodinger coupling,
density-gradient correction, or Bohm potential) is a future
M-numbered milestone, not M17. ADR 0016 documents this V&V
departure: the smooth-coefficient case is gated by MMS Variant I
(`semi/verification/mms_dd.py`); the discontinuous-coefficient case
is gated by this benchmark's classical-electrostatic comparison.

## References

- Pozela & Reklaitis, "Electron transport properties in GaAs at high
  electric fields," Solid-State Electronics 23, 927 (1980).
- Stern, "Self-consistent results for n-type Si inversion layers,"
  Phys. Rev. B 5, 4891 (1972). Foundational self-consistent
  Poisson-Schrodinger framework reused for AlGaAs/GaAs.
- Vurgaftman, Meyer, Ram-Mohan, "Band parameters for III-V compound
  semiconductors and their alloys," J. Appl. Phys. 89, 5815 (2001).
  Source for the AlGaAs_0p3 band-edge parameters at 300 K
  (`semi/materials.py::AlGaAs_0p3`).
- Sze, *Physics of Semiconductor Devices*, 3rd ed., Section 7.4
  ("Heterojunction field-effect transistors"). Classical
  electrostatic 2DEG derivation.
- ADR 0016 (`docs/adr/0016-heterojunction-slotboom.md`). M17 V&V
  scope and the algebraic preservation of the Slotboom continuity-
  flux shape under position-dependent `n_i(x)`.
