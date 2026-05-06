"""
Coupled drift-diffusion block residual (Slotboom form).

Primary unknowns: (psi_hat, phi_n_hat, phi_p_hat) in scaled units.
Densities are recovered pointwise via the Slotboom relations in
:mod:`semi.physics.slotboom`.

Scaled equations (see docs/PHYSICS.md section 2):

    Poisson:
        -div( L_D^2 eps_r grad psi_hat ) = p_hat - n_hat + N_hat

    Electron continuity:
        -div( L_0^2 mu_n_hat n_hat grad phi_n_hat ) = R_hat

    Hole continuity:
        -div( L_0^2 mu_p_hat p_hat grad phi_p_hat ) = -R_hat

where L_D^2 = lambda2 * L_0^2, mu_n_hat = mu_n / mu_0, and R_hat is the
scaled SRH rate (see :mod:`semi.physics.recombination`). The mesh stays
in physical meters (Invariant 3, see PLAN.md), so L_0^2 and L_D^2 appear
explicitly as UFL constants in the forms.

The residuals returned here are weak (Galerkin) forms with homogeneous
Neumann natural boundaries. Dirichlet contact data on psi, phi_n, phi_p
are applied via `dolfinx.fem.DirichletBC` by the caller.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass
class DDBlockSpaces:
    """Container for the three P1 subspaces and their unknowns."""
    V_psi: Any
    V_phi_n: Any
    V_phi_p: Any
    psi: Any
    phi_n: Any
    phi_p: Any


def make_dd_block_spaces(msh) -> DDBlockSpaces:
    """
    Create three independent P1 Lagrange spaces on `msh` and their
    unknown Functions.

    Using three separate scalar spaces (rather than a vector MixedElement
    on the same mesh) keeps the block residual assembly explicit and
    lets us re-use the M1 equilibrium Poisson solver for the initial
    guess.
    """
    from dolfinx import fem

    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    V_phi_n = fem.functionspace(msh, ("Lagrange", 1))
    V_phi_p = fem.functionspace(msh, ("Lagrange", 1))
    psi = fem.Function(V_psi, name="psi_hat")
    phi_n = fem.Function(V_phi_n, name="phi_n_hat")
    phi_p = fem.Function(V_phi_p, name="phi_p_hat")
    return DDBlockSpaces(V_psi, V_phi_n, V_phi_p, psi, phi_n, phi_p)


def build_dd_block_residual(
    spaces: DDBlockSpaces,
    N_hat_fn,
    sc,
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    E_t_over_Vt: float = 0.0,
    mobility_cfg: dict | None = None,
    *,
    facet_tags=None,
    recomb_cfg: dict | None = None,
    statistics_cfg: dict | None = None,
    schottky_facets: list | None = None,
    ref_mat=None,
):
    """
    Build the three-block residual for the coupled drift-diffusion system.

    Parameters
    ----------
    spaces : DDBlockSpaces
        Function spaces and unknowns (psi, phi_n, phi_p in scaled units).
    N_hat_fn : dolfinx.fem.Function
        Scaled net doping interpolated in V_psi or a compatible space.
    sc : semi.scaling.Scaling
        Scaling object; provides lambda2, L0, n_i/C0 ratio.
    eps_r : float | dolfinx.fem.Function
        Relative permittivity. Scalar (single-region fast path,
        byte-identical with M2-M5) or a DG0 cellwise Function on the
        parent mesh for the multi-region Poisson coefficient jump.
    mu_n_over_mu0, mu_p_over_mu0 : float
        Low-field mobility ratios (mu / sc.mu0). Under the
        caughey_thomas branch these are the `mu0` arguments to the
        closed form (the low-field reference).
    tau_n_hat, tau_p_hat : float
        Scaled lifetimes, tau / t0.
    E_t_over_Vt : float
        Trap level relative to the intrinsic level, divided by V_t.
    mobility_cfg : dict, optional
        The `cfg["physics"]["mobility"]` sub-dict. `None` (default) or
        `{"model": "constant"}` is bit-identical to pre-M16.1.
        `{"model": "caughey_thomas", ...}` substitutes a closed-form
        velocity-saturation expression (carrier-specific
        |grad(phi_n)| / |grad(phi_p)| under ADR 0004 Slotboom flux
        form; see docs/PHYSICS.md section 1.3 for the flux derivation
        and `semi.physics.mobility` for the closed form).
    recomb_cfg : dict, optional
        The merged `cfg["physics"]["recombination"]` and
        `cfg["physics"]["tunneling"]` sub-dicts. `None` (default) or
        `{"auger": False, "bbt": False, "tat": False, ...}` is
        bit-identical to pre-M16.3 (SRH-only recombination kernel,
        no tunneling). When `recomb_cfg.get("auger", False)` is True,
        the Auger kernel
        `R_Auger = (C_n_hat n_hat + C_p_hat p_hat) (n_hat p_hat
        - n_i_hat^2)` is added to the SRH rate; `C_n_hat` and
        `C_p_hat` are read from `recomb_cfg["C_n"]` / `["C_p"]` (in
        cm^6/s; converted to scaled units inline). M16.3.
        When `recomb_cfg.get("tat", False)` is True, the SRH rate is
        replaced by `(1 + Gamma_TAT(F)) * R_SRH` with the Hurkx
        field-enhancement factor `Gamma_TAT` evaluated at the
        scaled field magnitude `|grad(psi_hat)|`. M16.6.
        When `recomb_cfg.get("bbt", False)` is True, the Kane
        band-to-band generation rate `G_BBT` is subtracted from the
        net recombination (Kane is a generation, contributing
        `-G_BBT` to R). M16.6. See `semi.physics.recombination` for
        the closed forms and `scaled_kane_coefficients` /
        `scaled_hurkx_F_kT` for the unit conversions.
    statistics_cfg : dict, optional
        The `cfg["physics"]` sub-slice for carrier statistics. `None`
        (default) or `{"statistics": "boltzmann"}` is bit-identical
        to pre-M16.4 (textbook Slotboom). When
        `statistics_cfg["statistics"] == "fermi_dirac"`, the
        electron and hole densities are built via the generalized-
        Slotboom helpers in `semi.physics.slotboom` (Blakemore
        prefactor); the eta_offset values are read from
        `sc.eta_offset_n` / `sc.eta_offset_p`. ADR 0004 is preserved
        because the Einstein-factor cancellation in the
        generalized-Slotboom current means the continuity-row shape
        `J = -q mu n grad(phi)` is unchanged. M16.4.
    schottky_facets : list, optional
        List of ``(facet_tag, params)`` entries describing each
        Schottky contact whose continuity rows should carry the
        thermionic-emission Robin surface form. ``params`` must
        provide ``barrier_height_eV`` (the metal-semiconductor barrier
        in eV; Sze 3rd ed Table 5) and ``V_applied`` (the applied
        voltage in volts at this contact). ``ref_mat`` is consulted
        for ``Nc``, ``Nv``, and ``n_i`` so the thermionic equilibrium
        densities ``n_eq = N_C exp(-phi_B / V_t)`` and
        ``p_eq = N_V exp(-(E_g - phi_B) / V_t) = n_i^2 / n_eq`` can be
        formed in scaled units. ``None`` (default) or empty list is
        bit-identical to v0.20.0. When non-empty, ``facet_tags`` and
        ``ref_mat`` must also be supplied. M16.5.
    ref_mat : Material, optional
        Reference material; required when ``schottky_facets`` is
        non-empty. Used for the thermionic equilibrium-density
        formula. M16.5.

    Returns
    -------
    list of ufl.Form
        [F_psi, F_phi_n, F_phi_p] -- pass to the blocked NonlinearProblem.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .mobility import build_mobility_expressions
    from .slotboom import n_from_slotboom, p_from_slotboom

    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p
    msh = spaces.V_psi.mesh

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(msh, PETSc.ScalarType(sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    # The lombardi branch reads psi for the perpendicular field and
    # the absolute net doping field as N_total_hat. Other branches
    # ignore these kwargs.
    mu_n_hat, mu_p_hat, _mob_model = build_mobility_expressions(
        mobility_cfg, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc,
        psi=psi, facet_tags=facet_tags, N_total_hat=ufl.algebra.Abs(N_hat_fn),
    )
    tau_n = fem.Constant(msh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(msh, PETSc.ScalarType(tau_p_hat))

    eta_offset_n_val = (
        sc.eta_offset_n
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    eta_offset_p_val = (
        sc.eta_offset_p
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    n_hat = n_from_slotboom(
        psi, phi_n, ni_hat,
        statistics_cfg=statistics_cfg, eta_offset_n=eta_offset_n_val,
    )
    p_hat = p_from_slotboom(
        psi, phi_p, ni_hat,
        statistics_cfg=statistics_cfg, eta_offset_p=eta_offset_p_val,
    )

    # SRH (and optional Auger, M16.3) rate inlined here so we can share
    # the same ni_hat Constant with the Poisson block and avoid UFL-type
    # surprises.
    import math as _math
    n1 = ni_hat * _math.exp(E_t_over_Vt)
    p1 = ni_hat * _math.exp(-E_t_over_Vt)
    np_minus_nieq = n_hat * p_hat - ni_hat * ni_hat
    R_base = np_minus_nieq / (
        tau_p * (n_hat + n1) + tau_n * (p_hat + p1)
    )
    # M16.6: shared dimensionless field magnitude. Kane BBT and Hurkx
    # TAT both read the L_0-scaled magnitude `L_0 * |grad(psi_hat)|`,
    # which is the natural dimensionless ratio because
    # `|E_phys| = (V_0 / L_0) * (L_0 * |grad(psi_hat)|)`. The eps_F
    # term is a tiny additive guard so the ufl.sqrt derivative
    # remains finite at the rare interior points where grad(psi_hat)
    # is exactly zero (typical solves never see this in practice but
    # the SNES Jacobian assembly evaluates derivatives at the initial
    # guess where the field can be zero).
    tun_cfg = recomb_cfg or {}
    bbt_on = bool(tun_cfg.get("bbt", False))
    tat_on = bool(tun_cfg.get("tat", False))
    if bbt_on or tat_on:
        eps_F_const = fem.Constant(msh, PETSc.ScalarType(1.0e-30))
        L0_const = fem.Constant(msh, PETSc.ScalarType(sc.L0))
        F_hat_mag = L0_const * ufl.sqrt(
            ufl.dot(ufl.grad(psi), ufl.grad(psi)) + eps_F_const
        )
    else:
        F_hat_mag = None

    if tat_on:
        # M16.6 Hurkx TAT enhancement. The SRH rate is multiplied by
        # (1 + Gamma(F)); Gamma vanishes as F -> 0 so the bulk is
        # unchanged and the depletion-region peak field gets the
        # super-exponential enhancement.
        from .recombination import hurkx_gamma, scaled_hurkx_F_kT
        F_kT_cm = float(tun_cfg.get("F_kT", 1.4e7))
        alpha = float(tun_cfg.get("alpha", 2.0))
        F_kT_hat_val = scaled_hurkx_F_kT(F_kT_cm, sc)
        F_kT_hat = fem.Constant(msh, PETSc.ScalarType(F_kT_hat_val))
        Gamma_TAT = hurkx_gamma(F_hat_mag, F_kT_hat, alpha)
        R_SRH = (1.0 + Gamma_TAT) * R_base
    else:
        R_SRH = R_base

    R = R_SRH
    if recomb_cfg is not None and recomb_cfg.get("auger", False):
        # JSON contract: C_n / C_p in cm^6/s. Convert to m^6/s (1e-12)
        # then to dimensionless C_hat = C_SI * C0^2 * t0; see ADR 0002
        # and semi/physics/recombination.py module docstring (M16.3).
        C_n_hat_val = float(recomb_cfg.get("C_n", 2.8e-31)) * 1.0e-12 \
            * sc.C0 ** 2 * sc.t0
        C_p_hat_val = float(recomb_cfg.get("C_p", 9.9e-32)) * 1.0e-12 \
            * sc.C0 ** 2 * sc.t0
        C_n_hat = fem.Constant(msh, PETSc.ScalarType(C_n_hat_val))
        C_p_hat = fem.Constant(msh, PETSc.ScalarType(C_p_hat_val))
        R = R + (C_n_hat * n_hat + C_p_hat * p_hat) * np_minus_nieq

    if bbt_on:
        # M16.6 Kane BBT generation. G_BBT is positive; it subtracts
        # from R because the R = U - G convention here treats
        # generation as a negative contribution to net recombination.
        from .recombination import (
            bbt_rate,
            scaled_E_g,
            scaled_kane_coefficients,
        )
        A_kane_cm = float(tun_cfg.get("A_kane", 4.0e14))
        B_kane_cm = float(tun_cfg.get("B_kane", 1.9e7))
        kane = scaled_kane_coefficients(A_kane_cm, B_kane_cm, sc)
        E_g_eV_val = sc.E_g if sc.E_g > 0.0 else 1.12
        E_g_hat_val = scaled_E_g(E_g_eV_val, sc)
        A_kane_hat = fem.Constant(
            msh, PETSc.ScalarType(kane["A_kane_hat"])
        )
        B_kane_hat = fem.Constant(
            msh, PETSc.ScalarType(kane["B_kane_hat"])
        )
        E_g_hat = fem.Constant(msh, PETSc.ScalarType(E_g_hat_val))
        G_BBT = bbt_rate(F_hat_mag, E_g_hat, A_kane_hat, B_kane_hat)
        R = R - G_BBT

    rho_hat = p_hat - n_hat + N_hat_fn

    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * ufl.dx
        - rho_hat * v_psi * ufl.dx
    )
    F_phi_n = (
        L0_sq * mu_n_hat * n_hat * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * ufl.dx
        - R * v_n * ufl.dx
    )
    F_phi_p = (
        L0_sq * mu_p_hat * p_hat * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * ufl.dx
        + R * v_p * ufl.dx
    )

    if schottky_facets:
        F_phi_n_extra, F_phi_p_extra = _build_schottky_surface_forms(
            msh, facet_tags, schottky_facets, sc, ref_mat,
            psi=psi, phi_n=phi_n, phi_p=phi_p, v_n=v_n, v_p=v_p,
            statistics_cfg=statistics_cfg,
            eta_offset_n=eta_offset_n_val,
            eta_offset_p=eta_offset_p_val,
        )
        F_phi_n = F_phi_n + F_phi_n_extra
        F_phi_p = F_phi_p + F_phi_p_extra
    return [F_psi, F_phi_n, F_phi_p]


def _build_schottky_surface_forms(
    msh, facet_tags, schottky_facets, sc, ref_mat,
    *, psi, phi_n, phi_p, v_n, v_p,
    statistics_cfg=None, eta_offset_n=None, eta_offset_p=None,
    submesh_kwargs: dict | None = None,
):
    """
    Assemble the M16.5 Schottky thermionic-emission surface forms on the
    electron and hole continuity rows.

    The Robin condition at a Schottky contact uses the Sentaurus-style
    formulation (Selberherr 1984 Section 5.2; see Sze 3rd ed Section
    3.4 for the underlying thermionic-emission derivation):

        J_n . n_face = q v_n_th (n - n_eq)
        J_p . n_face = -q v_p_th (p - p_eq)

    where ``v_n_th = sqrt(kT / (2 pi m_n*))`` is the electron Richardson
    velocity, ``n_eq = N_C exp(-phi_B / V_t)`` is the metal-side
    equilibrium electron density (a constant, fixed by the barrier
    height), and ``p_eq = n_i^2 / n_eq`` is its hole counterpart.
    Mass-action law preserves ``n_eq * p_eq = n_i^2``.

    The applied bias enters through the ``psi`` Dirichlet at the
    Schottky boundary (``semi/bcs.py``); under Boltzmann statistics with
    Slotboom variables, the local boundary density satisfies
    ``n_local = n_eq exp((V_a - phi_n) / V_t)``, so the Robin term
    reduces to ``J_n . n = q v_n_th n_eq (exp((V_a - phi_n) / V_t) - 1)``,
    which under low injection (phi_n ~= 0 at the boundary) recovers the
    standard analytical thermionic-emission I-V
    ``J = A* T^2 exp(-q phi_B / kT) (exp(qV_a / kT) - 1)`` from Sze
    eq 3.4.10. The (V_a - phi_n) inside the exponential lives in
    ``n_local`` via the Slotboom relation, NOT in the Robin reference
    density n_eq.

    Substituting the Robin condition into the weak form's IBP boundary
    integral gives a residual prefactor of ``sc.t0 * v_n_th`` on the
    carrier-density / thermionic-flux difference, integrated against
    the test function on the Schottky facet.

    See ADR 0015 for the V&V scope (analytical-benchmark plus byte-
    identity gates instead of an MMS rate gate, since the existing MMS
    harness uses Dirichlet BCs everywhere by construction).
    """
    import math as _math

    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .slotboom import n_from_slotboom, p_from_slotboom

    if facet_tags is None:
        raise RuntimeError(
            "Schottky surface form requested but facet_tags is None; "
            "the form builder cannot resolve the Schottky facet measure."
        )
    if ref_mat is None:
        raise RuntimeError(
            "Schottky surface form requested but ref_mat is None; "
            "the thermionic equilibrium density formula needs Nc, Nv, "
            "and n_i from the reference material."
        )

    submesh_kwargs = submesh_kwargs or {}
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags,
                     **submesh_kwargs)
    V_t = sc.V0
    C0 = sc.C0
    n_i = ref_mat.n_i
    if n_i is None or n_i <= 0.0:
        raise RuntimeError(
            f"Schottky surface form requires the reference material to "
            f"expose a positive intrinsic density n_i; got {n_i}."
        )

    # Precompute sc.v_n_thermal / sc.v_p_thermal once so the property
    # accessor (which checks m_*_star) runs at form-build time, not on
    # the SNES jacobian assembly path.
    v_n_th = sc.v_n_thermal
    v_p_th = sc.v_p_thermal
    prefactor_n = sc.t0 * v_n_th
    prefactor_p = sc.t0 * v_p_th

    # Self-consistent CB / VB density-of-states from the thermionic
    # effective masses on `Scaling`. Computing N_C / N_V here (rather
    # than reading the DOS values on the material) ensures the Robin
    # form satisfies the textbook identity `A* T^2 = q v_R N_C` exactly
    # for the same m* used in `v_n_th`. Without this, the thermionic-
    # emission analytical match misses by a factor proportional to
    # (m_DOS / m_thermionic)^{3/2}.
    from ..constants import HBAR, KB, M0
    h_planck = 2.0 * _math.pi * HBAR
    kT = KB * sc.T
    Nc_TE = 2.0 * (
        2.0 * _math.pi * sc.m_n_star * M0 * kT / (h_planck ** 2)
    ) ** 1.5
    Nv_TE = 2.0 * (
        2.0 * _math.pi * sc.m_p_star * M0 * kT / (h_planck ** 2)
    ) ** 1.5

    # The Slotboom n_from_slotboom helper expects ni_hat as a
    # UFL-compatible Constant; build one on the same mesh that owns the
    # surface measure.
    ni_hat = fem.Constant(msh, PETSc.ScalarType(n_i / C0))

    F_phi_n_extra = 0
    F_phi_p_extra = 0
    for facet_tag, params in schottky_facets:
        phi_b = float(params["barrier_height_eV"])
        # V_applied is unused inside the Robin form: the Sentaurus-style
        # BC J_n . n = q v_n_th (n - n_eq) only references the fixed
        # metal-side equilibrium density. The applied bias enters the
        # solution exclusively through the psi Dirichlet at the
        # Schottky facet (set in semi/bcs.py); the Slotboom relation
        # then gives n_local = n_eq exp((V_a - phi_n) / V_t), and the
        # Robin term recovers the analytical thermionic-emission I-V
        # at low injection (phi_n ~= 0 at the boundary).
        # Use the self-consistent thermionic-emission DOS Nc_TE / Nv_TE
        # so q v_R N_C = A* T^2 holds with the same m* on both sides.
        n_eq_phys = Nc_TE * _math.exp(-phi_b / V_t)
        # p_eq_TE: hole-side metal-equilibrium density, computed from
        # the thermionic-mass-derived N_V and the gap-minus-barrier
        # exponent. This is identical in spirit to (n_i^2) / n_eq under
        # mass-action, but uses the thermionic effective masses for
        # both bands, keeping the Robin form internally consistent.
        Eg = ref_mat.Eg if ref_mat.Eg > 0.0 else 0.0
        p_eq_phys = Nv_TE * _math.exp(-(Eg - phi_b) / V_t) if Eg > 0.0 \
            else (n_i * n_i) / n_eq_phys
        n_eq_hat = fem.Constant(msh, PETSc.ScalarType(n_eq_phys / C0))
        p_eq_hat = fem.Constant(msh, PETSc.ScalarType(p_eq_phys / C0))
        pref_n = fem.Constant(msh, PETSc.ScalarType(prefactor_n))
        pref_p = fem.Constant(msh, PETSc.ScalarType(prefactor_p))

        n_hat_local = n_from_slotboom(
            psi, phi_n, ni_hat,
            statistics_cfg=statistics_cfg, eta_offset_n=eta_offset_n,
        )
        p_hat_local = p_from_slotboom(
            psi, phi_p, ni_hat,
            statistics_cfg=statistics_cfg, eta_offset_p=eta_offset_p,
        )
        # Electron Robin BC sign: kronos-semi tests
        # -div(mu n grad(phi_n)) = R, so the IBP boundary integral
        # -int (mu n grad(phi_n)).n_face v ds picks up a +/- sign
        # depending on whether n_face_outward points into or out of
        # the metal. The carrier-thermionic-emission relation in the
        # conventional-current sign (Sze 3rd ed eq 3.4) gives a
        # negative-of-(n - n_eq) contribution to F_phi_n in this
        # convention; bench-tuning against the closed-form I-V on
        # benchmarks/schottky_1d picks out this sign.
        #
        # Hole Robin BC: deliberately omitted. For an n-type Schottky
        # with phi_B in the upper half of the bandgap (Pt-on-n-Si:
        # phi_Bn = 0.85, phi_Bp = Eg - phi_Bn = 0.27), modeling the
        # metal as an infinite hole reservoir at p_eq would inject
        # unrealistically large minority hole currents that the bulk
        # minority diffusion does not actually supply. Standard
        # textbook Schottky analytical I-V (Sze 3rd ed Section 3.4)
        # is electron-thermionic-only; hole minority injection is a
        # second-order effect that requires a separate model. The
        # hole continuity row at the Schottky facet keeps the
        # natural homogeneous-Neumann condition (`J_p . n_face = 0`),
        # which models the contact as a hole-blocking boundary.
        # M16.6 (tunneling) and follow-ups can revisit this.
        ds_facet = ds(int(facet_tag))
        F_phi_n_extra = F_phi_n_extra - pref_n * (
            n_hat_local - n_eq_hat
        ) * v_n * ds_facet
        # Suppress unused-variable lint while documenting that the
        # hole-side equilibrium density is computed (above) for
        # completeness even though the hole Robin row is omitted.
        _ = (p_hat_local, p_eq_hat, pref_p)
    return F_phi_n_extra, F_phi_p_extra


@dataclass
class DDBlockSpacesMR:
    """
    Multi-region block spaces for the MOS capacitor and similar devices.

    V_psi lives on the parent mesh (spans silicon + oxide); V_phi_n and
    V_phi_p live on the semiconductor submesh only (Slotboom variables
    are ill-defined in an ideal insulator; see docs/mos_derivation.md
    section 2.2). `entity_map` is the parent<->submesh cell mapping
    returned by `dolfinx.mesh.create_submesh`; it threads through the
    mixed-domain form compilation (`fem.form(..., entity_maps=[em])`)
    and into the SNES solver (`NonlinearProblem(..., entity_maps=[em])`).
    """
    V_psi: Any
    V_phi_n: Any
    V_phi_p: Any
    psi: Any
    phi_n: Any
    phi_p: Any
    parent_mesh: Any
    submesh: Any
    entity_map: Any


def make_dd_block_spaces_mr(msh, submesh, entity_map) -> DDBlockSpacesMR:
    """Create the multi-region P1 spaces.

    V_psi lives on the parent mesh; V_phi_n and V_phi_p live on the
    semiconductor submesh. Unknown Functions are instantiated on their
    respective spaces.
    """
    from dolfinx import fem

    V_psi = fem.functionspace(msh, ("Lagrange", 1))
    V_phi_n = fem.functionspace(submesh, ("Lagrange", 1))
    V_phi_p = fem.functionspace(submesh, ("Lagrange", 1))
    psi = fem.Function(V_psi, name="psi_hat")
    phi_n = fem.Function(V_phi_n, name="phi_n_hat")
    phi_p = fem.Function(V_phi_p, name="phi_p_hat")
    return DDBlockSpacesMR(
        V_psi=V_psi, V_phi_n=V_phi_n, V_phi_p=V_phi_p,
        psi=psi, phi_n=phi_n, phi_p=phi_p,
        parent_mesh=msh, submesh=submesh, entity_map=entity_map,
    )


def build_dd_block_residual_mr(
    spaces: DDBlockSpacesMR,
    N_hat_fn,
    sc,
    eps_r,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    tau_n_hat: float,
    tau_p_hat: float,
    cell_tags,
    semi_tag: int,
    E_t_over_Vt: float = 0.0,
    mobility_cfg: dict | None = None,
    *,
    facet_tags=None,
    recomb_cfg: dict | None = None,
    statistics_cfg: dict | None = None,
    schottky_facets: list | None = None,
    ref_mat=None,
):
    """Multi-region (submesh-based) block residual for the coupled DD system.

    Structural differences from the single-region path:

    - The Poisson stiffness term integrates over the full parent mesh
      (oxide included), using the cellwise DG0 eps_r(x) Function. The
      space-charge term integrates only over silicon via the parent
      mesh's `dx(subdomain_data=cell_tags, subdomain_id=semi_tag)`
      measure, which is restricted to semiconductor cells where
      (phi_n, phi_p) live on the submesh.
    - The continuity-equation residuals integrate on the submesh
      (`Measure('dx', domain=submesh)`), reading psi from the parent
      mesh.

    The mixed-domain mapping threads through `NonlinearProblem`'s
    `entity_maps` argument; callers must pass `entity_maps=[spaces.entity_map]`
    to `solve_nonlinear_block`. The derivation is in
    docs/mos_derivation.md section 4.3.

    Returns
    -------
    list of ufl.Form
        [F_psi, F_phi_n, F_phi_p] -- pass to `solve_nonlinear_block`
        with `entity_maps=[spaces.entity_map]`.
    """
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    from .mobility import build_mobility_expressions
    from .slotboom import n_from_slotboom, p_from_slotboom

    psi = spaces.psi
    phi_n = spaces.phi_n
    phi_p = spaces.phi_p
    msh = spaces.parent_mesh
    submesh = spaces.submesh

    v_psi = ufl.TestFunction(spaces.V_psi)
    v_n = ufl.TestFunction(spaces.V_phi_n)
    v_p = ufl.TestFunction(spaces.V_phi_p)

    L_D2 = fem.Constant(msh, PETSc.ScalarType(sc.lambda2 * sc.L0 ** 2))
    L0_sq = fem.Constant(submesh, PETSc.ScalarType(sc.L0 ** 2))
    if isinstance(eps_r, (int, float)):
        eps_r_ufl = fem.Constant(msh, PETSc.ScalarType(float(eps_r)))
    else:
        eps_r_ufl = eps_r
    ni_hat_parent = fem.Constant(msh, PETSc.ScalarType(sc.n_i / sc.C0))
    ni_hat_sub = fem.Constant(submesh, PETSc.ScalarType(sc.n_i / sc.C0))
    # The lombardi branch reads psi (parent-mesh electrostatic
    # potential) and the absolute net doping. Other branches ignore
    # the optional kwargs.
    mu_n_hat, mu_p_hat, _mob_model = build_mobility_expressions(
        mobility_cfg, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc,
        psi=psi, facet_tags=facet_tags, N_total_hat=ufl.algebra.Abs(N_hat_fn),
    )
    tau_n = fem.Constant(submesh, PETSc.ScalarType(tau_n_hat))
    tau_p = fem.Constant(submesh, PETSc.ScalarType(tau_p_hat))

    # Measures
    dx_parent = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    dx_semi_parent = dx_parent(int(semi_tag))
    dx_sub = ufl.Measure("dx", domain=submesh)

    eta_offset_n_val = (
        sc.eta_offset_n
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )
    eta_offset_p_val = (
        sc.eta_offset_p
        if statistics_cfg is not None
        and statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"
        else None
    )

    # Parent-mesh Slotboom (for the Poisson source term, which is integrated
    # over silicon cells of the parent mesh). phi_n, phi_p live on the
    # submesh; entity_maps threads them through at form compilation.
    n_hat_parent = n_from_slotboom(
        psi, phi_n, ni_hat_parent,
        statistics_cfg=statistics_cfg, eta_offset_n=eta_offset_n_val,
    )
    p_hat_parent = p_from_slotboom(
        psi, phi_p, ni_hat_parent,
        statistics_cfg=statistics_cfg, eta_offset_p=eta_offset_p_val,
    )
    rho_hat = p_hat_parent - n_hat_parent + N_hat_fn

    # Submesh-mesh Slotboom (for the continuity blocks). Here psi is the
    # parent-mesh Function and entity_maps pulls it onto submesh cells.
    n_hat_sub = n_from_slotboom(
        psi, phi_n, ni_hat_sub,
        statistics_cfg=statistics_cfg, eta_offset_n=eta_offset_n_val,
    )
    p_hat_sub = p_from_slotboom(
        psi, phi_p, ni_hat_sub,
        statistics_cfg=statistics_cfg, eta_offset_p=eta_offset_p_val,
    )

    import math as _math
    n1 = ni_hat_sub * _math.exp(E_t_over_Vt)
    p1 = ni_hat_sub * _math.exp(-E_t_over_Vt)
    np_minus_nieq_sub = n_hat_sub * p_hat_sub - ni_hat_sub * ni_hat_sub
    R_base_mr = np_minus_nieq_sub / (
        tau_p * (n_hat_sub + n1) + tau_n * (p_hat_sub + p1)
    )
    # M16.6: shared L_0-scaled field magnitude on the parent mesh
    # (psi lives there). See the docstring of the single-region path
    # above for the dimensionless-ratio rationale.
    tun_cfg_mr = recomb_cfg or {}
    bbt_on_mr = bool(tun_cfg_mr.get("bbt", False))
    tat_on_mr = bool(tun_cfg_mr.get("tat", False))
    if bbt_on_mr or tat_on_mr:
        eps_F_const_mr = fem.Constant(msh, PETSc.ScalarType(1.0e-30))
        L0_const_mr = fem.Constant(msh, PETSc.ScalarType(sc.L0))
        F_hat_mag_mr = L0_const_mr * ufl.sqrt(
            ufl.dot(ufl.grad(psi), ufl.grad(psi)) + eps_F_const_mr
        )
    else:
        F_hat_mag_mr = None

    if tat_on_mr:
        from .recombination import hurkx_gamma, scaled_hurkx_F_kT
        F_kT_cm_mr = float(tun_cfg_mr.get("F_kT", 1.4e7))
        alpha_mr = float(tun_cfg_mr.get("alpha", 2.0))
        F_kT_hat_val_mr = scaled_hurkx_F_kT(F_kT_cm_mr, sc)
        F_kT_hat_mr = fem.Constant(
            submesh, PETSc.ScalarType(F_kT_hat_val_mr)
        )
        Gamma_TAT_mr = hurkx_gamma(F_hat_mag_mr, F_kT_hat_mr, alpha_mr)
        R_SRH_mr = (1.0 + Gamma_TAT_mr) * R_base_mr
    else:
        R_SRH_mr = R_base_mr

    R = R_SRH_mr
    if recomb_cfg is not None and recomb_cfg.get("auger", False):
        # M16.3 Auger inline. Submesh Constants because the continuity
        # blocks integrate on the semiconductor submesh; the Auger
        # term lives only there (no continuity rows in the oxide).
        C_n_hat_val_mr = float(recomb_cfg.get("C_n", 2.8e-31)) * 1.0e-12 \
            * sc.C0 ** 2 * sc.t0
        C_p_hat_val_mr = float(recomb_cfg.get("C_p", 9.9e-32)) * 1.0e-12 \
            * sc.C0 ** 2 * sc.t0
        C_n_hat_sub = fem.Constant(submesh, PETSc.ScalarType(C_n_hat_val_mr))
        C_p_hat_sub = fem.Constant(submesh, PETSc.ScalarType(C_p_hat_val_mr))
        R = R + (
            C_n_hat_sub * n_hat_sub + C_p_hat_sub * p_hat_sub
        ) * np_minus_nieq_sub

    if bbt_on_mr:
        from .recombination import (
            bbt_rate,
            scaled_E_g,
            scaled_kane_coefficients,
        )
        A_kane_cm_mr = float(tun_cfg_mr.get("A_kane", 4.0e14))
        B_kane_cm_mr = float(tun_cfg_mr.get("B_kane", 1.9e7))
        kane_mr = scaled_kane_coefficients(
            A_kane_cm_mr, B_kane_cm_mr, sc
        )
        E_g_eV_val_mr = sc.E_g if sc.E_g > 0.0 else 1.12
        E_g_hat_val_mr = scaled_E_g(E_g_eV_val_mr, sc)
        A_kane_hat_mr = fem.Constant(
            submesh, PETSc.ScalarType(kane_mr["A_kane_hat"])
        )
        B_kane_hat_mr = fem.Constant(
            submesh, PETSc.ScalarType(kane_mr["B_kane_hat"])
        )
        E_g_hat_mr = fem.Constant(
            submesh, PETSc.ScalarType(E_g_hat_val_mr)
        )
        G_BBT_mr = bbt_rate(
            F_hat_mag_mr, E_g_hat_mr, A_kane_hat_mr, B_kane_hat_mr
        )
        R = R - G_BBT_mr

    F_psi = (
        L_D2 * eps_r_ufl * ufl.inner(ufl.grad(psi), ufl.grad(v_psi)) * dx_parent
        - rho_hat * v_psi * dx_semi_parent
    )
    F_phi_n = (
        L0_sq * mu_n_hat * n_hat_sub
        * ufl.inner(ufl.grad(phi_n), ufl.grad(v_n)) * dx_sub
        - R * v_n * dx_sub
    )
    F_phi_p = (
        L0_sq * mu_p_hat * p_hat_sub
        * ufl.inner(ufl.grad(phi_p), ufl.grad(v_p)) * dx_sub
        + R * v_p * dx_sub
    )

    if schottky_facets:
        # Schottky contacts on the multi-region path are unusual (the
        # MOSCAP / MOSFET 2D devices have no metal-semiconductor
        # contacts on the silicon submesh), but the API accepts them
        # for completeness. The thermionic-emission Robin form is
        # applied on the parent mesh's facet measure, with phi_n and
        # phi_p restricted to the semiconductor submesh via the same
        # entity_maps mechanism the volume continuity rows already use.
        F_phi_n_extra, F_phi_p_extra = _build_schottky_surface_forms(
            msh, facet_tags, schottky_facets, sc, ref_mat,
            psi=psi, phi_n=phi_n, phi_p=phi_p, v_n=v_n, v_p=v_p,
            statistics_cfg=statistics_cfg,
            eta_offset_n=eta_offset_n_val,
            eta_offset_p=eta_offset_p_val,
        )
        F_phi_n = F_phi_n + F_phi_n_extra
        F_phi_p = F_phi_p + F_phi_p_extra
    return [F_psi, F_phi_n, F_phi_p]
