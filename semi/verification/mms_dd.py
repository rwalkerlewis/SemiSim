"""
Method of Manufactured Solutions (MMS) for the coupled drift-diffusion
block residual.

The production residual (`semi/physics/drift_diffusion.py`) assembles
three coupled blocks in Slotboom form:

    Poisson:   -div( L_D^2 eps_r grad psi_hat )            - (p_hat - n_hat + N_hat) = 0
    Electron:  -div( L_0^2 mu_n_hat n_hat grad phi_n_hat ) - R_hat                   = 0
    Hole:      -div( L_0^2 mu_p_hat p_hat grad phi_p_hat ) + R_hat                   = 0

with Slotboom `n_hat = ni_hat exp(psi - phi_n)`, `p_hat = ni_hat exp(phi_p - psi)`
and SRH `R_hat`. For MMS we set `N_hat = 0` (Variants A through D) or a
constant value (Variant E) and subtract a weak-form manufactured source
from each block so that a chosen smooth exact triple
`(psi_e, phi_n_e, phi_p_e)` is the solution of the modified coupled
system. Each block's discretization error then decays at the FE rate.

Seven variants exercise progressively more of the residual:

    A  `phi_n_e = phi_p_e = 0`, `tau_hat = infinity`:   Poisson block only
    B  full triple, `tau_hat = infinity`:               full coupling, R ~ 0
    C  full triple, Si lifetimes:                       full coupling with SRH
    D  Variant C plus Caughey-Thomas field-dependent mobility (M16.1).
    E  Variant C plus Lombardi composite surface mobility (M16.2). The
       composite is the resistor sum
           1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr
       with mu_bulk = mu_n_over_mu0 (constant), mu_AC the Lombardi
       acoustic-phonon term, and mu_sr the surface-roughness term.
       The perpendicular-field expression is
           E_perp_for_form = abs(grad(psi_e) . n_hat)
       with n_hat the unit vector along the configured
       `interface_normal_axis` (axis 0 for both 1D and 2D MMS so the
       manufactured "interface" is at x = 0 on the boundary; the
       direction matches the prompt's MMS Variant E geometry).
       N_hat is held at a constant `MMS_E_N_HAT_CONST` so the
       Lombardi C-term sees a non-zero N_total^lambda; the manufactured
       Poisson source absorbs the constant N contribution.
    F  Variant C plus the Auger recombination kernel (M16.3). The
       additive Auger term
           R_Auger = (C_n_hat n_hat + C_p_hat p_hat)
                     * (n_hat p_hat - n_i_hat^2)
       is appended to the SRH rate. C_n_hat and C_p_hat are engineered
       so the Auger contribution is comparable in magnitude to SRH at
       the typical manufactured amplitudes (cubic-in-density vs SRH's
       bilinear-over-linear; see derivation Section 3.5). The
       reverse-engineered JSON Auger coefficients are passed via
       `recomb_cfg = {"auger": True, "C_n": ..., "C_p": ...}` so
       `build_dd_block_residual` produces the same expression the
       manufactured weak source uses.
    G  Variant C plus the Fermi-Dirac generalized-Slotboom
       substitution (M16.4). The Slotboom helpers gain the Blakemore
       prefactor `gamma(eta) = 1 / (1 + 0.27 * exp(eta))` so the
       densities become `n = ni * gamma_n * exp(psi - phi_n)` and the
       hole counterpart. ADR 0004 (Slotboom variables) is preserved
       because the FD Einstein factor cancels against gamma in the
       continuity flux (the closed identity `g(eta) * gamma(eta) = 1`
       holds exactly under the basic Blakemore form). The
       reverse-engineered `Scaling.N_C` / `Scaling.N_V` values are
       written by `run_one_level` so `sc.eta_offset_n` /
       `sc.eta_offset_p` resolve to the engineered constants the
       manufactured weak source uses.

Variant D substitutes the closed-form
    mu(F) = mu0 / (1 + (mu0 * F_par / vsat)^beta)^(1/beta)
for the constant-mobility `fem.Constant` in both the production form
and the manufactured weak source. F_par is the carrier-specific
|grad(phi_e)|, taken from the same Slotboom-quasi-Fermi gradient that
the production form sees (ADR 0004). The MMS uses engineered
`vsat_for_form` and `beta` values so the CT factor is materially
different from 1 at the typical gradient of the manufactured solution
(~30 % mu reduction near the midplane); see
`MMS_D_VSAT_N_FOR_FORM` and `MMS_D_BETA_N` below.

Variant E substitutes the Lombardi composite (closed-form, see
`semi/physics/mobility.py::lombardi_compose`) for `mu_n_eff` and
`mu_p_eff` in both the production form and the manufactured weak
source. The MMS_E_* constants in this module are engineered so the
acoustic-phonon and surface-roughness terms each shift the composite
mobility by a few percent at the typical manufactured perpendicular
gradient, matching the O(0.3) reduction target M16.1 used for
Variant D. The JSON Lombardi parameters fed to
`build_mobility_expressions` are reverse-engineered from these
for-form constants so the production form sees the identical closed
form (matches the M16.1 Variant D vsat reverse-engineering pattern).

See `docs/mms_dd_derivation.md` for the approved derivation; every weak
source below tracks the sections of that document.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .mms_poisson import (
    EPS_R_DEFAULT,
    L_0_REF,
    _all_boundary_facets,
    _build_mesh,
    build_mms_scaling,
)

# ---------------------------------------------------------------------------
# Module-level constants. See derivation Section 2.
# ---------------------------------------------------------------------------
VARIANTS: tuple[str, ...] = ("A", "B", "C", "D", "E", "F", "G", "H", "I")

#: Variant I (M17): smooth chi(x), Eg(x) ramps producing position-
#: dependent n_i(x). The manufactured solution reuses the Variant C
#: full sin-product triple; the active branch exercises the
#: heterojunction substitution rule
#:     n = n_i(x) * exp((psi - phi_n) / V_t)
#: documented in ADR 0016. The discontinuous-coefficient case (HEMT 2D)
#: is gated by the published-reference `benchmarks/hemt_2d/` verifier
#: rather than MMS, per ADR 0016. Variant I rate gates: L^2 >= 1.99 /
#: H^1 >= 0.99 finest-pair on each block, the textbook P1 gate. Until
#: the run_one_level Variant I branch is fully wired (Phase E follow-
#: up: build per-cell DG0 chi_hat(x) / Eg_hat(x) ramps, pass via
#: `heterojunction_fields` to build_dd_block_residual, mirror the
#: position-dependent n_i in `_build_weak_sources`), Variant I is
#: registered here for forward-compatibility and run_one_level raises
#: NotImplementedError when invoked. The infrastructure (Phases B-D)
#: is shipped; the MMS-side wiring is the missing piece. M17.
MMS_I_DELTA_CHI_EV: float = 0.05
MMS_I_DELTA_EG_EV: float = 0.05

#: Default amplitude triple `(A_psi, A_n, A_p)` used by the pytest gate.
#: `A_p = -A_n` ensures Variant C's SRH numerator
#: `ni_hat^2 * (exp(phi_p - phi_n) - 1)` is pointwise nonzero (derivation 1).
DEFAULT_AMPS: tuple[float, float, float] = (0.5, 0.3, -0.3)

#: Larger-amplitude triple reserved for the CLI sweep. Keeps |exp(...)|
#: bounded by exp(1.5) ~ 4.48, comfortably inside PETSc double precision.
NONLINEAR_AMPS: tuple[float, float, float] = (1.0, 0.5, -0.5)

#: Variant C SRH lifetime (Si mid-gap default, seconds).
TAU_SI_S: float = 1.0e-7

#: Variant A/B "off" lifetime in scaled units. Keeps the SRH code path
#: live but drives `R_e` below machine precision.
TAU_OFF_HAT: float = 1.0e+20

#: Mobility ratios pinned across the mesh sweep so convergence measures
#: discretization error only. `mu_n_over_mu0 = 1.0` because mu_0 is the
#: electron mobility itself; `mu_p / mu_n = 0.45 / 1.4` from Si defaults
#: in `semi/materials.py`.
MU_N_OVER_MU0: float = 1.0
MU_P_OVER_MU0: float = 0.45 / 1.4

#: Variant D Caughey-Thomas parameters. Engineered to push the
#: dimensionless ratio (mu0_hat * F_par_e / vsat_for_form)^beta into
#: the O(0.3) regime at the typical manufactured gradient
#: (|grad(phi_n_e)| ~ A_n * pi / L ~ 4.7e5 m^-1 with the default
#: amplitudes and L = L_0_REF = 2 um), so mu(F)/mu0 ~ 0.93 and the CT
#: nonlinearity is materially exercised. Values are not physical Si
#: numbers; the goal is to verify discretization rate, not match a
#: device. The mobility module's vsat_for_form has units 1/m so we set
#: it directly here rather than going through the cm/s -> 1/m converter.
MMS_D_VSAT_N_FOR_FORM: float = 1.5e6
MMS_D_VSAT_P_FOR_FORM: float = 1.5e6
MMS_D_BETA_N: float = 2.0
MMS_D_BETA_P: float = 1.0

#: Variant E Lombardi parameters in for-form (scaled) units.
#: Engineered so each surface term shifts the composite mu by a few
#: percent at the typical manufactured perpendicular gradient. With
#: amplitudes (A_psi=0.5) and L = L_0_REF = 2 um, the typical
#: |grad(psi_e) . e_x| at x close to the interface is
#: A_psi * pi / L ~ 7.85e5 m^-1; mu_AC ~ B/E ~ 5 (so 1/mu_AC ~ 0.2
#: contributes ~20% of 1/mu_bulk), mu_sr ~ delta/E^2 ~ 5 likewise.
#: The Lombardi composite then reduces mu by ~30% vs the bulk branch
#: (matching the M16.1 Variant D design anchor).
MMS_E_B_FOR_FORM_N: float = 3.93e6
MMS_E_B_FOR_FORM_P: float = 3.93e6
MMS_E_C_EFF_FOR_FORM_N: float = 92.0
MMS_E_C_EFF_FOR_FORM_P: float = 92.0
MMS_E_LAMBDA_N: float = 0.125
MMS_E_LAMBDA_P: float = 0.0317
MMS_E_DELTA_FOR_FORM_N: float = 3.08e12
MMS_E_DELTA_FOR_FORM_P: float = 3.08e12
MMS_E_T_K: float = 300.0

#: Variant E constant N_hat. The Lombardi C-term reads
#: `N_total_hat^lambda`; we hold N_hat at this constant value so the
#: per-cell N_total^lambda is well-defined and non-trivial. The
#: manufactured Poisson source absorbs the resulting (p - n + N_hat)
#: with N_hat = MMS_E_N_HAT_CONST.
MMS_E_N_HAT_CONST: float = 1.0

#: Variant E interface normal axis. Axis 0 corresponds to the line
#: x = 0 (the boundary) for both 1D and 2D MMS meshes; n_hat = e_x.
MMS_E_INTERFACE_NORMAL_AXIS: int = 0

#: Variant F Auger coefficients in for-form (scaled) units. With
#: C_0_REF = 1e22 m^-3, t_0 ~ 1.1e-9 s, and the default amplitudes
#: (psi 0.5, phi_n 0.3, phi_p -0.3), n_hat / p_hat at the manufactured
#: peak are O(1e-6) (= ni_hat = 1e-6 inflated by exp(O(0.5))). The
#: ratio R_Auger / R_SRH at the peak is approximately
#: 4 * C_n_hat * n_hat^2 * tau_n_hat ~ 4 * C_n_hat * 1e-12 * 91, so
#: C_n_hat ~ 8e8 brings the Auger contribution to ~30 % of SRH (the
#: same O(0.3) reduction target M16.1 used for Variant D and M16.2
#: used for Variant E). Values are engineered for materially-exercised
#: discretization, not physical Si; the JSON Auger coefficients passed
#: to build_dd_block_residual are reverse-engineered from these
#: for-form constants in run_one_level so the production form sees
#: the identical closed form.
MMS_F_C_N_HAT_FOR_FORM: float = 8.0e8
MMS_F_C_P_HAT_FOR_FORM: float = 8.0e8

#: Variant G Fermi-Dirac eta-offsets. The reduced Fermi level under
#: the generalized-Slotboom substitution is
#:     eta_n = (psi_e - phi_n_e) + eta_offset_n
#: With the default amplitudes (A_psi=0.5, A_n=0.3) the (psi - phi_n)
#: drive ranges over roughly [-0.8, +0.8]. Setting eta_offset_n = -1.0
#: pushes the typical eta to roughly [-1.8, -0.2], so the Blakemore
#: prefactor gamma_n = 1 / (1 + 0.27 * exp(eta)) ranges from ~0.96
#: (~4 % FD-vs-Boltzmann shift in the deep-non-degenerate region) to
#: ~0.82 (~18 % shift in the manufactured-degenerate peak), materially
#: exercising the Blakemore branch without driving the closed form
#: outside its sub-5 % accuracy window. The hole side mirrors the
#: electron side; eta_offset_p = -1.0 is the symmetric choice.
#: M16.4.
MMS_G_ETA_OFFSET_N: float = -1.0
MMS_G_ETA_OFFSET_P: float = -1.0

#: Variant H Kane and Hurkx for-form constants (M16.6). Engineered so
#: each tunneling kernel shifts the total recombination rate by O(0.1)
#: of R_SRH at the typical manufactured amplitudes, mirroring the
#: Variant F O(0.3) reduction target. R_SRH at the manufactured peak
#: is ~ |np - n_i^2| / denom ~ -1.0e-9 in scaled units (ni_hat=1e-6
#: with the default A_psi=0.5 / A_n=0.3 / A_p=-0.3 amplitudes); the
#: BBT and TAT magnitudes are tuned so each contributes ~ 0.1 of that.
#:
#: The peak L_0-scaled field magnitude is
#:     F_hat_peak = L_0 * |grad(psi_e)|_peak = A_psi * 2*pi = pi,
#: so:
#:   - Kane prefactor MMS_H_A_KANE_FOR_FORM = 3.5e-11 with
#:     MMS_H_B_KANE_FOR_FORM = 0.5, MMS_H_E_G_HAT_FOR_FORM = 1.0 gives
#:     G_BBT_hat = 3.5e-11 * pi^2 * exp(-0.5 / pi) ~ 3.0e-10 at peak;
#:   - Hurkx characteristic field MMS_H_F_KT_FOR_FORM = 20.0 with
#:     MMS_H_ALPHA = 2.0 gives Gamma(F_hat=pi) ~ 0.16 (16 % SRH
#:     enhancement, ~ 0.1 |R_SRH| absolute shift).
#:
#: Reverse-engineered JSON values feed into build_dd_block_residual via
#: run_one_level so the production form sees the identical closed forms
#: at the manufactured triple. The conversions in
#: semi/physics/recombination.py are
#:     A_kane_hat = A_kane_cm * 100 * V_0^(3/2) * t_0 / (L_0^2 * C_0)
#:     B_kane_hat = B_kane_cm * 100 * L_0 * sqrt(V_0)
#:     F_kT_hat   = F_kT_cm   * 100 * L_0 / V_0
#: so the inverse formulas in run_one_level produce JSON-side numbers
#: that hit the for-form constants exactly.
MMS_H_A_KANE_FOR_FORM: float = 3.5e-11
MMS_H_B_KANE_FOR_FORM: float = 0.5
MMS_H_F_KT_FOR_FORM: float = 20.0
MMS_H_ALPHA: float = 2.0
MMS_H_E_G_HAT_FOR_FORM: float = 1.0


@dataclass(frozen=True)
class MMSDDCase:
    """Specification for one mesh-level MMS-DD run."""

    dim: int                                 # 1 or 2
    N: int                                   # cells per side
    variant: str                             # "A", "B", or "C"
    L: float = L_0_REF                       # device length, m
    A_psi: float = DEFAULT_AMPS[0]
    A_n: float = DEFAULT_AMPS[1]
    A_p: float = DEFAULT_AMPS[2]
    eps_r: float = EPS_R_DEFAULT
    cell_kind: str = "triangle"              # 2D only


@dataclass(frozen=True)
class MMSDDResult:
    """Per-level numerical result with per-block error norms."""

    dim: int
    N: int
    variant: str
    h: float
    n_dofs: int                              # identical across the three blocks
    e_L2_psi: float
    e_H1_psi: float
    e_L2_phi_n: float
    e_H1_phi_n: float
    e_L2_phi_p: float
    e_H1_phi_p: float
    snes_iters: int
    solve_time_s: float
    cell_kind: str = "triangle"
    A_psi: float = DEFAULT_AMPS[0]
    A_n: float = DEFAULT_AMPS[1]
    A_p: float = DEFAULT_AMPS[2]


# ---------------------------------------------------------------------------
# Exact triple and manufactured weak sources (derivation Sections 1, 3)
# ---------------------------------------------------------------------------


def _exact_triple_ufl(mesh, dim: int, L: float,
                      A_psi: float, A_n: float, A_p: float,
                      variant: str):
    """
    Return `(psi_e, phi_n_e, phi_p_e)` as UFL expressions.

    For Variants B and C we use the full sin-product triple; for Variant
    A the Fermi levels are identically zero. This choice matters for the
    weak-source construction in `_build_weak_sources`: Variant A only
    produces `f_psi_weak` (the other two manufactured sources are zero
    and are not subtracted from the production forms, which avoids
    compiling zero integrands).
    """
    import ufl

    # Variants D, E, F, G, and H share the Variant B/C smooth-sin
    # triple; the active branch is what differs (CT mobility for D,
    # Lombardi mobility for E, Auger recombination for F, Fermi-Dirac
    # statistics for G, BBT and TAT tunneling for H), evaluated in
    # `_build_weak_sources` and dispatched into the production form
    # via `run_one_level`.
    full_triple = variant in ("B", "C", "D", "E", "F", "G", "H")
    x = ufl.SpatialCoordinate(mesh)
    if dim == 1:
        psi_e = A_psi * ufl.sin(2.0 * ufl.pi * x[0] / L)
        if not full_triple:
            phi_n_e = None
            phi_p_e = None
        else:
            phi_n_e = A_n * ufl.sin(ufl.pi * x[0] / L)
            phi_p_e = A_p * ufl.sin(ufl.pi * x[0] / L)
    elif dim == 2:
        psi_e = (
            A_psi
            * ufl.sin(2.0 * ufl.pi * x[0] / L)
            * ufl.sin(3.0 * ufl.pi * x[1] / L)
        )
        if not full_triple:
            phi_n_e = None
            phi_p_e = None
        else:
            phi_n_e = (
                A_n
                * ufl.sin(ufl.pi * x[0] / L)
                * ufl.sin(ufl.pi * x[1] / L)
            )
            phi_p_e = (
                A_p
                * ufl.sin(ufl.pi * x[0] / L)
                * ufl.sin(ufl.pi * x[1] / L)
            )
    else:
        raise ValueError(f"Unsupported dim {dim}")
    return psi_e, phi_n_e, phi_p_e


def _build_weak_sources(
    mesh,
    spaces,
    sc,
    *,
    variant: str,
    eps_r_value: float,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    psi_e,
    phi_n_e,
    phi_p_e,
    N_hat_const: float = 0.0,
):
    """
    Build the three weak-form manufactured sources.

    For Variants B and C this returns `(f_psi_weak, f_n_weak, f_p_weak)`.
    For Variant A only `f_psi_weak` is built and the continuity-block
    entries are `None` -- the caller must leave the production continuity
    forms unmodified (they are already satisfied by `(psi_e, 0, 0)` because
    `grad(phi_n_e) = 0` and `R_e = 0` there).

    The weak-form construction follows Section 3 of
    `docs/mms_dd_derivation.md`. We never form `ufl.div(ufl.grad(...))`
    directly; integrating by parts against the test function uses the
    same gradient discretization as the production residual, so MMS
    verifies the exact discretization in use (compare
    `mms_poisson.manufactured_source_weak`).
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from semi.physics.slotboom import n_from_slotboom, p_from_slotboom

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(mesh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(mesh, PETSc.ScalarType(sc.L0 ** 2))
    eps_r = fem.Constant(mesh, PETSc.ScalarType(eps_r_value))
    ni_hat = fem.Constant(mesh, PETSc.ScalarType(sc.n_i / sc.C0))
    mu_n = fem.Constant(mesh, PETSc.ScalarType(mu_n_over_mu0))
    mu_p = fem.Constant(mesh, PETSc.ScalarType(mu_p_over_mu0))
    tau_n = fem.Constant(mesh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(mesh, PETSc.ScalarType(tau_p_hat))

    if variant == "A":
        # Variant A: phi_n_e = phi_p_e = 0, so n_e * p_e = ni_hat^2 exactly
        # and R_e vanishes pointwise. Slotboom densities reduce to
        # ni_hat * exp(+/- psi_e), matching the Poisson MMS carrier term.
        n_e = ni_hat * ufl.exp(psi_e)
        p_e = ni_hat * ufl.exp(-psi_e)
        f_psi_weak = (
            L_D2 * eps_r * ufl.inner(ufl.grad(psi_e), ufl.grad(v_psi))
            - (p_e - n_e) * v_psi
        ) * ufl.dx
        return f_psi_weak, None, None

    # Variants B, C, D, E, F, G: full triple. Variant G activates the
    # Fermi-Dirac dispatch on the Slotboom helpers (the generalized-
    # Slotboom substitution under Blakemore basic). The eta_offset
    # values match what `run_one_level` writes to `sc.N_C` / `sc.N_V`,
    # so the production form sees the identical closed form.
    if variant == "G":
        stat_cfg_local = {"statistics": "fermi_dirac"}
        eta_off_n_local = MMS_G_ETA_OFFSET_N
        eta_off_p_local = MMS_G_ETA_OFFSET_P
    else:
        stat_cfg_local = None
        eta_off_n_local = None
        eta_off_p_local = None
    n_e = n_from_slotboom(
        psi_e, phi_n_e, ni_hat,
        statistics_cfg=stat_cfg_local, eta_offset_n=eta_off_n_local,
    )
    p_e = p_from_slotboom(
        psi_e, phi_p_e, ni_hat,
        statistics_cfg=stat_cfg_local, eta_offset_p=eta_off_p_local,
    )

    # SRH rate, inlined to share the same ni_hat Constant with the Poisson
    # block (matches production code in `build_dd_block_residual`). E_t = 0
    # is enforced module-wide, so n1 = p1 = ni_hat.
    n1 = ni_hat * math.exp(0.0)
    p1 = ni_hat * math.exp(-0.0)
    np_minus_nieq_e = n_e * p_e - ni_hat * ni_hat
    R_base_e = np_minus_nieq_e / (
        tau_p * (n_e + n1) + tau_n * (p_e + p1)
    )
    # Variant H: Hurkx field-enhancement Gamma(F) multiplies the SRH
    # rate. Built from the L_0-scaled gradient of psi_e to match the
    # production form's field expression.
    if variant == "H":
        eps_F_const_h = fem.Constant(mesh, PETSc.ScalarType(1.0e-30))
        L0_const_h = fem.Constant(mesh, PETSc.ScalarType(sc.L0))
        F_hat_e = L0_const_h * ufl.sqrt(
            ufl.dot(ufl.grad(psi_e), ufl.grad(psi_e)) + eps_F_const_h
        )
        F_kT_hat_c = fem.Constant(
            mesh, PETSc.ScalarType(MMS_H_F_KT_FOR_FORM)
        )
        alpha_c = fem.Constant(mesh, PETSc.ScalarType(MMS_H_ALPHA))
        from semi.physics.recombination import hurkx_gamma
        Gamma_e = hurkx_gamma(F_hat_e, F_kT_hat_c, alpha_c)
        R_e = (1.0 + Gamma_e) * R_base_e
    else:
        R_e = R_base_e
    if variant == "F":
        # Auger contribution at the manufactured triple, using the
        # for-form C_n_hat / C_p_hat. The reverse-engineered JSON
        # values fed into build_dd_block_residual via run_one_level
        # produce the identical closed form when their cm^6/s ->
        # m^6/s -> dimensionless conversion lands on these constants.
        C_n_hat_c = fem.Constant(mesh, PETSc.ScalarType(MMS_F_C_N_HAT_FOR_FORM))
        C_p_hat_c = fem.Constant(mesh, PETSc.ScalarType(MMS_F_C_P_HAT_FOR_FORM))
        R_e = R_e + (C_n_hat_c * n_e + C_p_hat_c * p_e) * np_minus_nieq_e
    if variant == "H":
        # Subtract the Kane BBT generation. The same L0-scaled field
        # used by Gamma above feeds bbt_rate; constants match the
        # for-form values reverse-engineered into JSON cm-based units
        # in run_one_level.
        from semi.physics.recombination import bbt_rate
        A_kane_hat_c = fem.Constant(
            mesh, PETSc.ScalarType(MMS_H_A_KANE_FOR_FORM)
        )
        B_kane_hat_c = fem.Constant(
            mesh, PETSc.ScalarType(MMS_H_B_KANE_FOR_FORM)
        )
        E_g_hat_c = fem.Constant(
            mesh, PETSc.ScalarType(MMS_H_E_G_HAT_FOR_FORM)
        )
        G_BBT_e = bbt_rate(F_hat_e, E_g_hat_c, A_kane_hat_c, B_kane_hat_c)
        R_e = R_e - G_BBT_e

    if variant == "D":
        # Caughey-Thomas: substitute the closed-form mobility evaluated
        # at the manufactured gradient. The same vsat_for_form / beta
        # values feed the production form via run_one_level so the
        # weak source matches the production residual at phi_e exactly.
        from semi.physics.mobility import caughey_thomas_mu

        vsat_n_c = fem.Constant(mesh, PETSc.ScalarType(MMS_D_VSAT_N_FOR_FORM))
        vsat_p_c = fem.Constant(mesh, PETSc.ScalarType(MMS_D_VSAT_P_FOR_FORM))
        beta_n_c = fem.Constant(mesh, PETSc.ScalarType(MMS_D_BETA_N))
        beta_p_c = fem.Constant(mesh, PETSc.ScalarType(MMS_D_BETA_P))
        eps_n_c = fem.Constant(mesh, PETSc.ScalarType(1.0e-30))
        eps_p_c = fem.Constant(mesh, PETSc.ScalarType(1.0e-30))
        F_par_n_e = ufl.sqrt(
            ufl.dot(ufl.grad(phi_n_e), ufl.grad(phi_n_e)) + eps_n_c
        )
        F_par_p_e = ufl.sqrt(
            ufl.dot(ufl.grad(phi_p_e), ufl.grad(phi_p_e)) + eps_p_c
        )
        mu_n_eff = caughey_thomas_mu(mu_n, F_par_n_e, vsat_n_c, beta_n_c)
        mu_p_eff = caughey_thomas_mu(mu_p, F_par_p_e, vsat_p_c, beta_p_c)
    elif variant == "E":
        # Lombardi composite: substitute
        #   1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr
        # evaluated at the manufactured E_perp = abs(grad(psi_e) . e_x)
        # and at the constant manufactured N_total = MMS_E_N_HAT_CONST.
        # The reverse-engineered Lombardi JSON parameters in
        # run_one_level make build_mobility_expressions produce the
        # identical closed form, so the weak source matches the
        # production residual at (psi_e, phi_n_e, phi_p_e) exactly.
        from semi.physics.mobility import (
            lombardi_compose,
            lombardi_mu_AC,
            lombardi_mu_sr,
        )

        gdim = mesh.geometry.dim
        normal_components = [0.0] * gdim
        normal_components[MMS_E_INTERFACE_NORMAL_AXIS] = 1.0
        n_hat_vec = ufl.as_vector([
            fem.Constant(mesh, PETSc.ScalarType(c)) for c in normal_components
        ])
        eps_perp = fem.Constant(mesh, PETSc.ScalarType(1.0e-30))
        E_perp_signed = ufl.dot(ufl.grad(psi_e), n_hat_vec)
        E_perp_e = ufl.sqrt(E_perp_signed * E_perp_signed + eps_perp)

        T_const = fem.Constant(mesh, PETSc.ScalarType(MMS_E_T_K))
        N_total_const = fem.Constant(mesh, PETSc.ScalarType(N_hat_const))
        B_n_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_B_FOR_FORM_N))
        B_p_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_B_FOR_FORM_P))
        C_n_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_C_EFF_FOR_FORM_N))
        C_p_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_C_EFF_FOR_FORM_P))
        lam_n_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_LAMBDA_N))
        lam_p_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_LAMBDA_P))
        delta_n_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_DELTA_FOR_FORM_N))
        delta_p_c = fem.Constant(mesh, PETSc.ScalarType(MMS_E_DELTA_FOR_FORM_P))

        mu_AC_n = lombardi_mu_AC(B_n_c, C_n_c, N_total_const, E_perp_e, lam_n_c, T_const)
        mu_AC_p = lombardi_mu_AC(B_p_c, C_p_c, N_total_const, E_perp_e, lam_p_c, T_const)
        mu_sr_n = lombardi_mu_sr(delta_n_c, E_perp_e)
        mu_sr_p = lombardi_mu_sr(delta_p_c, E_perp_e)
        mu_n_eff = lombardi_compose(mu_n, mu_AC_n, mu_sr_n)
        mu_p_eff = lombardi_compose(mu_p, mu_AC_p, mu_sr_p)
    else:
        mu_n_eff = mu_n
        mu_p_eff = mu_p

    # Constant-N contribution to the manufactured Poisson source: the
    # production residual carries `+ N_hat_fn * v_psi` so the weak
    # forcing must mirror it. Variants A through D set N_hat_fn = 0
    # and `N_hat_const` defaults to 0; Variant E sets a constant.
    N_hat_const_c = fem.Constant(mesh, PETSc.ScalarType(N_hat_const))
    f_psi_weak = (
        L_D2 * eps_r * ufl.inner(ufl.grad(psi_e), ufl.grad(v_psi))
        - (p_e - n_e + N_hat_const_c) * v_psi
    ) * ufl.dx
    f_n_weak = (
        L0_sq * mu_n_eff * n_e * ufl.inner(ufl.grad(phi_n_e), ufl.grad(v_n))
        - R_e * v_n
    ) * ufl.dx
    f_p_weak = (
        L0_sq * mu_p_eff * p_e * ufl.inner(ufl.grad(phi_p_e), ufl.grad(v_p))
        + R_e * v_p
    ) * ufl.dx
    return f_psi_weak, f_n_weak, f_p_weak


# ---------------------------------------------------------------------------
# Per-level solve (derivation Sections 4, 5)
# ---------------------------------------------------------------------------


def run_one_level(case: MMSDDCase, *, sc=None) -> MMSDDResult:
    """
    Solve the MMS-DD problem on a single mesh level and return per-block
    L^2 and H^1 seminorm errors together with SNES diagnostics.

    Reuses `build_dd_block_residual` so future edits to the production
    residual are automatically exercised. BCs are homogeneous Dirichlet
    on all three fields (exact, not interpolated, because every exact
    factor is `sin(k*pi*x/L)` with `k` a positive integer and so vanishes
    on the boundary; see derivation Section 4).
    """
    from dolfinx import fem
    from petsc4py import PETSc

    from semi.physics.drift_diffusion import (
        build_dd_block_residual,
        make_dd_block_spaces,
    )
    from semi.solver import solve_nonlinear_block

    from ._norms import h1_seminorm_error_squared, l2_error_squared

    if case.variant not in VARIANTS:
        raise ValueError(f"variant {case.variant!r} not in {VARIANTS}")
    if case.variant == "I":
        # M17 Phase E follow-up: the smooth-coefficient heterojunction
        # rate gate is registered for forward-compatibility but the
        # MMS-side wiring (per-cell DG0 chi(x) / Eg(x) ramps,
        # position-dependent n_i(x) in `_build_weak_sources`) is
        # tracked as a deferred follow-up. The Phase B / C / D
        # infrastructure (heterojunction.build_dg0_material_fields,
        # form-builder threading, ohmic local-chi shift) is shipped;
        # the MMS verifier extension lands in a follow-up PR. The
        # discontinuous-coefficient HEMT case has its own V&V gate via
        # `benchmarks/hemt_2d/` per ADR 0016.
        raise NotImplementedError(
            "MMS Variant I (M17 heterojunction smooth ramp) is registered "
            "in VARIANTS for forward-compatibility; full rate-gate "
            "implementation is a Phase E follow-up. See ADR 0016 and "
            "docs/M17_STARTER_PROMPT.md Phase E for the wiring plan."
        )
    if sc is None:
        sc = build_mms_scaling(L=case.L)

    mesh = _build_mesh(case.dim, case.N, case.L, case.cell_kind)
    spaces = make_dd_block_spaces(mesh)

    # N_hat = 0 for MMS Variants A-D; Variant E uses a constant
    # MMS_E_N_HAT_CONST so the Lombardi C-term sees a non-trivial
    # N_total^lambda. The manufactured Poisson source is adjusted
    # accordingly inside _build_weak_sources.
    N_hat_fn = fem.Function(spaces.V_psi, name="N_hat_const")
    N_hat_const_value = (
        MMS_E_N_HAT_CONST if case.variant == "E" else 0.0
    )
    N_hat_fn.x.array[:] = N_hat_const_value
    N_hat_fn.x.scatter_forward()

    # Initial guess: zero everywhere. Trivially satisfies the
    # homogeneous Dirichlet BC and the charge-neutrality identity
    # (n_init = p_init = ni_hat, rho_init = 0 since N_hat = 0).
    # Variant H also starts at zero; Newton converges to the BBT-
    # active basin via the engineered mild kernel parameters and
    # the dimension-aware SNES atol below.
    for fn in (spaces.psi, spaces.phi_n, spaces.phi_p):
        fn.x.array[:] = 0.0
        fn.x.scatter_forward()

    # Per-variant lifetime choice (derivation Section 2, SRH row).
    # Variants C, D, E, F, G, and H use Si lifetimes (D adds CT
    # mobility, E adds the Lombardi composite, F adds Auger
    # recombination, G adds Fermi-Dirac statistics, H adds Kane and
    # Hurkx tunneling, all on top of the C electronics).
    if case.variant in ("C", "D", "E", "F", "G", "H"):
        tau_n_hat = tau_p_hat = TAU_SI_S / sc.t0
    else:
        tau_n_hat = tau_p_hat = TAU_OFF_HAT

    # Variant G writes engineered N_C / N_V values to the Scaling
    # object so `sc.eta_offset_n` and `sc.eta_offset_p` resolve to the
    # MMS_G_ETA_OFFSET_* constants the manufactured weak source uses.
    # eta_offset = ln(n_i / N_C) so N_C = n_i * exp(-eta_offset). The
    # boltzmann-default branches on every other variant skip this
    # write; sc.eta_offset_n is never read on those branches.
    if case.variant == "G":
        sc.N_C = sc.n_i * math.exp(-MMS_G_ETA_OFFSET_N)
        sc.N_V = sc.n_i * math.exp(-MMS_G_ETA_OFFSET_P)

    # Per-variant mobility dispatch. Variant D activates Caughey-Thomas
    # via the production-form mobility_cfg; the JSON-side vsat_*_cm_per_s
    # values are reverse-engineered from the module-level
    # MMS_D_VSAT_*_FOR_FORM so build_mobility_expressions produces the
    # exact closed form the manufactured weak source uses
    # (caughey_thomas_vsat_for_form: vsat_for_form = vsat_SI / (sc.mu0 *
    # sc.V0); vsat_SI = vsat_cm_per_s * 1e-2).
    if case.variant == "D":
        scale = sc.mu0 * sc.V0 * 100.0
        mobility_cfg = {
            "model": "caughey_thomas",
            "vsat_n": MMS_D_VSAT_N_FOR_FORM * scale,
            "vsat_p": MMS_D_VSAT_P_FOR_FORM * scale,
            "beta_n": MMS_D_BETA_N,
            "beta_p": MMS_D_BETA_P,
        }
    elif case.variant == "E":
        # Reverse-engineer the JSON Lombardi parameters so
        # `lombardi_unit_conversions` returns the same for-form values
        # the manufactured weak source uses
        # (MMS_E_*_FOR_FORM). The conversions in
        # semi/physics/mobility.py::lombardi_unit_conversions are:
        #   B_for_form    = B_cm * 1e-2 / (V_t * mu0_ref)
        #   delta_for_form = delta_cm * 1e-4 / (V_t**2 * mu0_ref)
        #   C_for_form_geom = C_cm * 1e-(10/3) / (V_t**(1/3) * mu0_ref)
        #   C_eff = C_for_form_geom * sc.C0**lambda
        # Inverting:
        cm_to_m = 1.0e-2
        cm5_3 = cm_to_m ** (5.0 / 3.0)
        cm2 = cm_to_m ** 2

        def _to_cm_B(B_for_form: float) -> float:
            return B_for_form * sc.V0 * sc.mu0 / cm_to_m

        def _to_cm_delta(delta_for_form: float) -> float:
            return delta_for_form * (sc.V0 ** 2) * sc.mu0 / cm2

        def _to_cm_C(C_eff_for_form: float, lam: float) -> float:
            C_geom = C_eff_for_form / (sc.C0 ** lam)
            return C_geom * (sc.V0 ** (1.0 / 3.0)) * sc.mu0 / cm5_3

        mobility_cfg = {
            "model": "lombardi",
            "bulk_model": "constant",
            "interface_facet_tag": 1,
            "interface_normal_axis": MMS_E_INTERFACE_NORMAL_AXIS,
            "temperature_K": MMS_E_T_K,
            "lombardi": {
                "B_n": _to_cm_B(MMS_E_B_FOR_FORM_N),
                "B_p": _to_cm_B(MMS_E_B_FOR_FORM_P),
                "C_n": _to_cm_C(MMS_E_C_EFF_FOR_FORM_N, MMS_E_LAMBDA_N),
                "C_p": _to_cm_C(MMS_E_C_EFF_FOR_FORM_P, MMS_E_LAMBDA_P),
                "lambda_n": MMS_E_LAMBDA_N,
                "lambda_p": MMS_E_LAMBDA_P,
                "delta_n": _to_cm_delta(MMS_E_DELTA_FOR_FORM_N),
                "delta_p": _to_cm_delta(MMS_E_DELTA_FOR_FORM_P),
            },
        }
    else:
        mobility_cfg = None

    # Variant F passes recomb_cfg to build_dd_block_residual so the
    # production form sees the identical Auger closed form the
    # manufactured weak source uses. The runner-side conversion is
    #   C_n_hat = C_n_cm * 1e-12 * sc.C0**2 * sc.t0
    # so we invert to obtain the cm^6/s value the JSON contract
    # accepts.
    recomb_cfg = None
    if case.variant == "F":
        cm_to_si = 1.0e-12
        scale_aug = cm_to_si * (sc.C0 ** 2) * sc.t0
        recomb_cfg = {
            "auger": True,
            "C_n": MMS_F_C_N_HAT_FOR_FORM / scale_aug,
            "C_p": MMS_F_C_P_HAT_FOR_FORM / scale_aug,
        }
    elif case.variant == "H":
        # Reverse-engineer the JSON cm-based values so
        # scaled_kane_coefficients / scaled_hurkx_F_kT lands on the
        # MMS_H_*_FOR_FORM constants used by the manufactured weak
        # source. Also overwrite sc.E_g so the production form's
        # scaled_E_g(sc.E_g, sc) = MMS_H_E_G_HAT_FOR_FORM.
        A_kane_cm_h = (
            MMS_H_A_KANE_FOR_FORM
            * sc.L0 ** 2 * sc.C0
            / (sc.V0 ** 1.5 * sc.t0)
            / 100.0
        )
        B_kane_cm_h = (
            MMS_H_B_KANE_FOR_FORM
            / (sc.L0 * math.sqrt(sc.V0))
            / 100.0
        )
        F_kT_cm_h = (
            MMS_H_F_KT_FOR_FORM * sc.V0 / sc.L0 / 100.0
        )
        recomb_cfg = {
            "bbt": True,
            "tat": True,
            "A_kane": A_kane_cm_h,
            "B_kane": B_kane_cm_h,
            "F_kT": F_kT_cm_h,
            "alpha": MMS_H_ALPHA,
        }
        sc.E_g = MMS_H_E_G_HAT_FOR_FORM * sc.V0

    # Variant E needs a `facet_tags` MeshTags so the lombardi UFL
    # dispatch passes the runner-wiring guard. The MMS doesn't have a
    # geometrically meaningful Si/SiO2 interface; `_build_lombardi`
    # only consults `interface_normal_axis` for the normal direction,
    # so we tag every boundary facet with value 1 as a sentinel.
    variant_facet_tags = None
    if case.variant == "E":
        from dolfinx import mesh as _mesh_mod

        fdim_local = mesh.topology.dim - 1
        mesh.topology.create_connectivity(fdim_local, mesh.topology.dim)
        boundary_facets_local = _all_boundary_facets(mesh, fdim_local)
        indices = np.asarray(boundary_facets_local, dtype=np.int32)
        values = np.ones(indices.shape, dtype=np.int32)
        order = np.argsort(indices)
        variant_facet_tags = _mesh_mod.meshtags(
            mesh, fdim_local, indices[order], values[order]
        )

    # Variant G activates the Fermi-Dirac dispatch on the production
    # residual; eta_offset_n / eta_offset_p come from the engineered
    # sc.N_C / sc.N_V written above.
    statistics_cfg = (
        {"statistics": "fermi_dirac"} if case.variant == "G" else None
    )

    # Production residual, per-block.
    F_prod = build_dd_block_residual(
        spaces, N_hat_fn, sc,
        case.eps_r, MU_N_OVER_MU0, MU_P_OVER_MU0,
        tau_n_hat, tau_p_hat, 0.0,   # E_t / V_t = 0 (mid-gap traps)
        mobility_cfg=mobility_cfg,
        facet_tags=variant_facet_tags,
        recomb_cfg=recomb_cfg,
        statistics_cfg=statistics_cfg,
    )

    # Manufactured weak sources.
    psi_e, phi_n_e, phi_p_e = _exact_triple_ufl(
        mesh, case.dim, case.L,
        case.A_psi, case.A_n, case.A_p, case.variant,
    )
    f_psi_weak, f_n_weak, f_p_weak = _build_weak_sources(
        mesh, spaces, sc,
        variant=case.variant,
        eps_r_value=case.eps_r,
        mu_n_over_mu0=MU_N_OVER_MU0,
        mu_p_over_mu0=MU_P_OVER_MU0,
        tau_n_hat=tau_n_hat,
        tau_p_hat=tau_p_hat,
        psi_e=psi_e,
        phi_n_e=phi_n_e,
        phi_p_e=phi_p_e,
        N_hat_const=N_hat_const_value,
    )
    F_psi = F_prod[0] - f_psi_weak
    if case.variant == "A":
        # f_n_weak = f_p_weak = 0. Leave continuity forms unmodified; they
        # are already satisfied by (psi_e, 0, 0) (derivation Section 3.4).
        F_phi_n = F_prod[1]
        F_phi_p = F_prod[2]
    else:
        F_phi_n = F_prod[1] - f_n_weak
        F_phi_p = F_prod[2] - f_p_weak

    # Homogeneous Dirichlet BCs on all three fields, exact.
    fdim = mesh.topology.dim - 1
    mesh.topology.create_connectivity(fdim, mesh.topology.dim)
    boundary_facets = _all_boundary_facets(mesh, fdim)
    zero = fem.Constant(mesh, PETSc.ScalarType(0.0))
    bdofs_psi = fem.locate_dofs_topological(spaces.V_psi, fdim, boundary_facets)
    bdofs_phi_n = fem.locate_dofs_topological(spaces.V_phi_n, fdim, boundary_facets)
    bdofs_phi_p = fem.locate_dofs_topological(spaces.V_phi_p, fdim, boundary_facets)
    bcs = [
        fem.dirichletbc(zero, bdofs_psi, spaces.V_psi),
        fem.dirichletbc(zero, bdofs_phi_n, spaces.V_phi_n),
        fem.dirichletbc(zero, bdofs_phi_p, spaces.V_phi_p),
    ]

    # SNES tolerances: rtol and max_it match the derivation (Section 5),
    # but atol is driven below the continuity-block initial residual
    # scale. Section 7.3 identifies the block-residual scale disparity
    # (L_D^2 * ...  vs  L_0^2 * n_i_hat * ...) as the reason atol must
    # reach below every block's native floor; the derivation's nominal
    # 1e-16 is tuned to the psi block alone and prematurely terminates
    # SNES on the continuity blocks in 2D (where the integration-area
    # L^2 ~ 4e-12 shrinks the phi_n / phi_p residuals to ~1e-18 at
    # initial). Driving atol to ~0 with stol controlling termination
    # lets Newton iterate until every block has actually been solved.
    petsc_options = {
        "snes_rtol": 1.0e-14,
        "snes_atol": 0.0,
        "snes_stol": 1.0e-12,
        "snes_max_it": 80,
    }
    # Variant H: the Kane / Hurkx Jacobian is locally ill-conditioned
    # in the BBT-active basin; Newton converges in O(30) iterations
    # but tends to overshoot near the discrete solution where the
    # exp(-B/F) factor saturates. A modest snes_atol above the
    # discretization floor and a higher snes_max_it let SNES exit
    # cleanly without driving Newton into the F-domain failure mode
    # (reason=-8, NaN in exp). Discretization-rate gates remain at
    # the L2 >= 1.99 / H1 >= 0.99 acceptance thresholds.
    if case.variant == "H":
        # Variant H absolute tolerance is dimension-aware: in 1D the
        # manufactured continuity-row residual is O(1e-13) at the
        # typical amplitudes, but in 2D the integration area
        # (L_0^2 ~ 4e-12) shrinks the per-block residual to O(1e-15)
        # and SNES with atol < 1e-14 stagnates near machine noise.
        # The dimension-keyed atol keeps Newton iterating just enough
        # to resolve the discretization-error envelope without
        # pushing into floating-point round-off.
        if case.dim == 1:
            petsc_options["snes_atol"] = 1.0e-15
        else:
            petsc_options["snes_atol"] = 1.0e-13
        petsc_options["snes_max_it"] = 100
    prefix = (
        f"mms_dd_dim{case.dim}_N{case.N}_var{case.variant}"
        f"_amp{case.A_psi:g}_"
    )
    t0 = time.perf_counter()
    # Variant H: small Jacobian shift covers the null-pivot risk that
    # MUMPS reports when the BBT-driven Jacobian is locally
    # ill-conditioned (the exp(-B/F) factor's derivative is
    # super-sensitive near F = 0). Other variants pass shift=0.0 so
    # the Variant A through G byte-identity invariant holds.
    jac_shift = 1.0e-14 if case.variant == "H" else 0.0
    info = solve_nonlinear_block(
        [F_psi, F_phi_n, F_phi_p],
        [spaces.psi, spaces.phi_n, spaces.phi_p],
        bcs, prefix=prefix, petsc_options=petsc_options,
        jacobian_shift=jac_shift,
    )
    dt = time.perf_counter() - t0
    if not info.get("converged", False):
        raise RuntimeError(
            f"MMS-DD did not converge at dim={case.dim}, N={case.N}, "
            f"variant={case.variant}: reason={info.get('reason')}, "
            f"iters={info.get('iterations')}"
        )

    # Per-block error norms. For Variant A the exact Fermi fields are 0,
    # so `phi_n_e_err` below is a UFL Zero and the error reduces to the
    # discrete function's own L^2 and H^1 seminorm.
    psi_e_err = psi_e
    phi_n_e_err = phi_n_e if phi_n_e is not None else 0.0
    phi_p_e_err = phi_p_e if phi_p_e is not None else 0.0

    e_L2_psi = float(np.sqrt(max(l2_error_squared(spaces.psi, psi_e_err), 0.0)))
    e_H1_psi = float(np.sqrt(max(h1_seminorm_error_squared(spaces.psi, psi_e_err), 0.0)))
    e_L2_phi_n = float(np.sqrt(max(l2_error_squared(spaces.phi_n, phi_n_e_err), 0.0)))
    e_H1_phi_n = float(np.sqrt(max(h1_seminorm_error_squared(spaces.phi_n, phi_n_e_err), 0.0)))
    e_L2_phi_p = float(np.sqrt(max(l2_error_squared(spaces.phi_p, phi_p_e_err), 0.0)))
    e_H1_phi_p = float(np.sqrt(max(h1_seminorm_error_squared(spaces.phi_p, phi_p_e_err), 0.0)))

    n_dofs = int(spaces.V_psi.dofmap.index_map.size_global)
    return MMSDDResult(
        dim=case.dim,
        N=case.N,
        variant=case.variant,
        h=case.L / case.N,
        n_dofs=n_dofs,
        e_L2_psi=e_L2_psi, e_H1_psi=e_H1_psi,
        e_L2_phi_n=e_L2_phi_n, e_H1_phi_n=e_H1_phi_n,
        e_L2_phi_p=e_L2_phi_p, e_H1_phi_p=e_H1_phi_p,
        snes_iters=int(info.get("iterations", -1)),
        solve_time_s=dt,
        cell_kind=case.cell_kind,
        A_psi=case.A_psi, A_n=case.A_n, A_p=case.A_p,
    )


# ---------------------------------------------------------------------------
# Convergence sweeps (derivation Section 6)
# ---------------------------------------------------------------------------


def run_convergence_study(
    *,
    dim: int,
    variant: str,
    Ns: list[int],
    L: float = L_0_REF,
    A_psi: float = DEFAULT_AMPS[0],
    A_n: float = DEFAULT_AMPS[1],
    A_p: float = DEFAULT_AMPS[2],
    eps_r: float = EPS_R_DEFAULT,
    cell_kind: str = "triangle",
    sc=None,
) -> list[MMSDDResult]:
    """Run MMS-DD on a sequence of refinements. Returns one result per N."""
    if sc is None:
        sc = build_mms_scaling(L=L)
    results: list[MMSDDResult] = []
    for N in Ns:
        case = MMSDDCase(
            dim=dim, N=N, variant=variant, L=L,
            A_psi=A_psi, A_n=A_n, A_p=A_p,
            eps_r=eps_r, cell_kind=cell_kind,
        )
        results.append(run_one_level(case, sc=sc))
    return results


# Column order used by the tables / CSVs.
_BLOCKS: tuple[str, ...] = ("psi", "phi_n", "phi_p")
_ROW_COLUMNS: tuple[str, ...] = (
    "N", "h", "N_dofs",
    "e_L2_psi", "rate_L2_psi", "e_H1_psi", "rate_H1_psi",
    "e_L2_phi_n", "rate_L2_phi_n", "e_H1_phi_n", "rate_H1_phi_n",
    "e_L2_phi_p", "rate_L2_phi_p", "e_H1_phi_p", "rate_H1_phi_p",
    "snes_iters", "solve_time_s",
)


def to_table_rows(results: list[MMSDDResult]) -> list[dict]:  # pragma: no cover
    """Build per-level rows with observed L^2/H^1 rates for every block."""
    from ._convergence import observed_rates

    hs = [r.h for r in results]
    rates: dict[str, list[float]] = {}
    for b in _BLOCKS:
        rates[f"L2_{b}"] = observed_rates(hs, [getattr(r, f"e_L2_{b}") for r in results])
        rates[f"H1_{b}"] = observed_rates(hs, [getattr(r, f"e_H1_{b}") for r in results])
    rows: list[dict] = []
    for i, r in enumerate(results):
        row = {
            "N": r.N,
            "h": r.h,
            "N_dofs": r.n_dofs,
            "e_L2_psi": r.e_L2_psi,
            "rate_L2_psi": rates["L2_psi"][i],
            "e_H1_psi": r.e_H1_psi,
            "rate_H1_psi": rates["H1_psi"][i],
            "e_L2_phi_n": r.e_L2_phi_n,
            "rate_L2_phi_n": rates["L2_phi_n"][i],
            "e_H1_phi_n": r.e_H1_phi_n,
            "rate_H1_phi_n": rates["H1_phi_n"][i],
            "e_L2_phi_p": r.e_L2_phi_p,
            "rate_L2_phi_p": rates["L2_phi_p"][i],
            "e_H1_phi_p": r.e_H1_phi_p,
            "rate_H1_phi_p": rates["H1_phi_p"][i],
            "snes_iters": r.snes_iters,
            "solve_time_s": r.solve_time_s,
        }
        rows.append(row)
    return rows


def report_table(rows: list[dict], header: str = "") -> str:  # pragma: no cover
    """Pretty multi-line stdout table showing psi L^2 + rate for each block."""
    from ._convergence import format_table

    cols = [
        "N", "h", "N_dofs",
        "e_L2_psi", "rate_L2_psi",
        "e_L2_phi_n", "rate_L2_phi_n",
        "e_L2_phi_p", "rate_L2_phi_p",
    ]
    return format_table(rows, cols, header=header)


def write_artifacts(  # pragma: no cover
    rows: list[dict],
    out_dir: Path,
    *,
    title: str,
    csv_name: str = "convergence.csv",
    plot_name: str = "convergence.png",
) -> None:
    """Write `convergence.csv` and `convergence.png` to `out_dir`."""
    from ._convergence import write_convergence_csv, write_loglog_plot

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_convergence_csv(out_dir / csv_name, rows, list(_ROW_COLUMNS))
    hs = [r["h"] for r in rows]
    series = {
        "L2 psi":    [r["e_L2_psi"] for r in rows],
        "L2 phi_n":  [r["e_L2_phi_n"] for r in rows],
        "L2 phi_p":  [r["e_L2_phi_p"] for r in rows],
    }
    write_loglog_plot(
        out_dir / plot_name, hs, series,
        title=title,
        theoretical_rates={
            "L2 psi": 2.0,
            "L2 phi_n": 2.0,
            "L2 phi_p": 2.0,
        },
    )


# ---------------------------------------------------------------------------
# Default mesh sequences
# ---------------------------------------------------------------------------
PYTEST_NS_1D: list[int] = [40, 80, 160]
PYTEST_NS_2D: list[int] = [16, 32, 64]
CLI_NS_1D: list[int] = [40, 80, 160, 320]
CLI_NS_2D: list[int] = [16, 32, 64]
# Variant D 2D sequence. The [16, 32, 64] sequence used by A/B/C
# bottoms out at finest-pair rate 1.990 from triangle-mesh
# boundary-layer effects, just under the M16.1 acceptance floor of
# 1.99. One extra refinement level (N=128) pulls the finest-pair
# rate cleanly to ~1.997.
CLI_NS_2D_D: list[int] = [32, 64, 128]


def run_cli_study(out_dir: Path) -> dict[str, list[dict]]:  # pragma: no cover
    """
    Artifact-production sweep used by `scripts/run_verification.py`.

    Runs (per variant, seven variants total) three studies:
      - `1d_<variant>_linear`     (default amps, Ns=CLI_NS_1D)
      - `1d_<variant>_nonlinear`  (NONLINEAR_AMPS, Ns=CLI_NS_1D)
      - `2d_<variant>`            (default amps, Ns=CLI_NS_2D)

    Variant D activates the M16.1 Caughey-Thomas mobility dispatch on
    top of the Variant C electronics; Variant E activates the M16.2
    Lombardi composite; Variant F activates the M16.3 Auger
    recombination kernel; Variant G activates the M16.4 Fermi-Dirac
    statistics dispatch. D, E, F, and G share the residual-floor
    sensitivity and boundary-layer behavior so each reuses the
    truncated-1D / extended-2D N sequences D pioneered.

    Each study writes `convergence.csv` and `convergence.png` into
    `out_dir/<label>/`, mirroring the Poisson MMS artifact layout.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sc = build_mms_scaling()

    studies: dict[str, list[dict]] = {}

    for variant in VARIANTS:
        # M17 Phase E: Variant I is registered for forward-compatibility
        # but the MMS-side wiring is a deferred follow-up; skip it in
        # the CLI sweep so existing variants A-H continue to ship
        # rate-gated artifacts. ADR 0016 documents the V&V gate
        # structure (Variant I MMS for the smooth-coefficient case,
        # `benchmarks/hemt_2d/` for the discontinuous case).
        if variant == "I":
            continue
        # 1D default amplitudes. Variant D (Caughey-Thomas) on the
        # finest 1D mesh (N=320) reaches the double-precision residual
        # floor (~1e-22) before SNES_rtol trips at 1e-14*||F_initial||,
        # so reason -8 (DIVERGED_LINE_SEARCH) is reported even though
        # the residual is effectively zero. Use one less refinement
        # level so SNES converges cleanly; the discretization rate is
        # already demonstrated at N=160 with rate ~ 2.000 (rates 1.999
        # at N=80 and 2.000 at N=160 in the linear-amp run). M16.3
        # Variant F shares the same nonlinear-residual behavior as D
        # and E.
        ns_1d = CLI_NS_1D[:-1] if variant in ("D", "E", "F", "G") else CLI_NS_1D
        ns_1d_d_nonlinear = (
            CLI_NS_1D[:-1] if variant in ("D", "E", "F", "G") else CLI_NS_1D
        )
        res = run_convergence_study(
            dim=1, variant=variant, Ns=ns_1d,
            A_psi=DEFAULT_AMPS[0], A_n=DEFAULT_AMPS[1], A_p=DEFAULT_AMPS[2],
            sc=sc,
        )
        rows = to_table_rows(res)
        label = f"1d_{variant}_linear"
        write_artifacts(
            rows, out_dir / label,
            title=(
                f"MMS-DD 1D variant {variant} (linear, "
                f"A_psi={DEFAULT_AMPS[0]}, A_n={DEFAULT_AMPS[1]}, "
                f"A_p={DEFAULT_AMPS[2]})"
            ),
        )
        studies[label] = rows

        # 1D nonlinear amplitudes
        res = run_convergence_study(
            dim=1, variant=variant, Ns=ns_1d_d_nonlinear,
            A_psi=NONLINEAR_AMPS[0], A_n=NONLINEAR_AMPS[1], A_p=NONLINEAR_AMPS[2],
            sc=sc,
        )
        rows = to_table_rows(res)
        label = f"1d_{variant}_nonlinear"
        write_artifacts(
            rows, out_dir / label,
            title=(
                f"MMS-DD 1D variant {variant} (nonlinear, "
                f"A_psi={NONLINEAR_AMPS[0]}, A_n={NONLINEAR_AMPS[1]}, "
                f"A_p={NONLINEAR_AMPS[2]})"
            ),
        )
        studies[label] = rows

        # 2D default amplitudes, triangles. Variant D needs one extra
        # refinement level to clear the M16.1 acceptance gate
        # (rate >= 1.99); the [16, 32, 64] sequence used by A/B/C
        # bottoms out at ~1.99 from triangle-mesh boundary-layer
        # effects. See CLI_NS_2D_D for the rationale.
        ns_2d = CLI_NS_2D_D if variant in ("D", "E", "F", "G") else CLI_NS_2D
        res = run_convergence_study(
            dim=2, variant=variant, Ns=ns_2d,
            A_psi=DEFAULT_AMPS[0], A_n=DEFAULT_AMPS[1], A_p=DEFAULT_AMPS[2],
            cell_kind="triangle", sc=sc,
        )
        rows = to_table_rows(res)
        label = f"2d_{variant}"
        write_artifacts(
            rows, out_dir / label,
            title=(
                f"MMS-DD 2D variant {variant} (right-diagonal triangles, "
                f"A_psi={DEFAULT_AMPS[0]}, A_n={DEFAULT_AMPS[1]}, "
                f"A_p={DEFAULT_AMPS[2]})"
            ),
        )
        studies[label] = rows

    return studies
