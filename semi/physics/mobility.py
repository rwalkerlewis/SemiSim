"""
Field-dependent carrier mobility models (M16.1, M16.2).

Layer 4 (FEM). Imports of dolfinx, ufl, and petsc4py are deferred to
function bodies so the pure-Python core stays dolfinx-free per
PLAN.md invariant 4 (covered by tests/test_lazy_imports.py).

Currently supported:
    constant       mu(F) = mu0 (no field dependence). The V&V reference;
                   bit-identical to pre-M16.1.
    caughey_thomas mu(F) = mu0 / (1 + (mu0 * F_par / vsat)^beta)^(1/beta)
                   Closed-form velocity saturation. F_par is the
                   magnitude of the gradient of the carrier-specific
                   quasi-Fermi potential; mu0 is the low-field
                   electron / hole mobility.
    lombardi       Composite surface mobility (M16.2). Resistor sum
                       1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr
                   where mu_bulk is the constant or caughey_thomas
                   bulk branch (selected by `bulk_model`), mu_AC is
                   the Lombardi acoustic-phonon term, and mu_sr is
                   the surface-roughness term. The closed-form math
                   is exposed via `lombardi_mu_AC`, `lombardi_mu_sr`,
                   and `lombardi_compose` (pure-Python; testable
                   against literal restatements). The UFL dispatch in
                   `build_mobility_expressions` requires the runner
                   to thread the parent mesh facet_tags so the
                   interface normal can be located; Phase C wires
                   that through bias_sweep and mos_cap_ac.

Asymptotics for caughey_thomas:
    F_par -> 0:    mu(F) -> mu0
    F_par -> inf:  mu(F) * F_par -> vsat / mu0 * mu0 = vsat (drift
                                                               saturation)

Asymptotics for lombardi (resistor sum 1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr):
    E_perp -> 0:        mu_AC, mu_sr -> +inf, mu -> mu_bulk (the
                        composite naturally relaxes to the bulk away
                        from the interface).
    mu_bulk -> +inf:    mu -> 1 / (1/mu_AC + 1/mu_sr) (parallel
                        combination of the surface terms).

Scaling. The DD form integrates `mu_hat * carrier_density * grad(phi)`
against `grad(test)` in scaled units (Slotboom form, see
docs/PHYSICS.md section 2.5 and ADR 0004). The mesh stays in physical
meters (PLAN.md invariant 3), so `ufl.grad(phi_n_hat)` has units of
1/m; the field magnitude is `F = V_t * grad(phi_hat)` with units V/m,
and the SI velocity argument of the saturation formula is
`mu0 * F = mu0 * V_t * grad(phi_hat)`. Dividing by `vsat` gives a
dimensionless ratio when

    vsat_for_form = vsat_SI / (sc.mu0 * sc.V0)        [units 1/m]

so the closed form

    mu_hat(F) = mu0_hat / (1 + (mu0_hat * F_par / vsat_for_form)^beta)^(1/beta)

is dimensionless. `F_par` here is `sqrt(grad(phi).grad(phi) + eps)`
straight from UFL (units 1/m). See `caughey_thomas_vsat_for_form` for
the conversion helper.

The regularization epsilon is tiny (1e-30 in scaled-form units) so
that the `(F/vsat)^beta` term is below the floor of float precision
and shifts mu only at machine precision; near the equilibrium initial
guess where `grad(phi) = 0` it averts a `0^0` UFL nan.

Lombardi surface mobility (M16.2)
=================================

The Lombardi 1988 composite (a.k.a. Synopsys Sentaurus "EnormalDependence")
adds two surface-scattering terms to a bulk mobility branch via the
resistor sum

    1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr.

Component formulas (Lombardi 1988 / Sentaurus convention; SI for the
resulting mobility, device-physics units for the JSON parameters):

    mu_AC(E_perp, N_total) = B / E_perp
                              + C * N_total^lambda / (E_perp^(1/3) * T)
    mu_sr(E_perp)          = delta / E_perp^2

with

    B       acoustic-phonon prefactor          (cm/s)
    C       Coulomb-phonon coefficient         (cm^(5/3) V^(-2/3) s^-1)
    lambda  doping exponent                    (dimensionless)
    delta   surface-roughness coefficient      (cm^2 V s^-1)
    T       device temperature                 (K)
    E_perp  perpendicular field at interface   (V/cm in JSON;
                                                converted to V/m
                                                internally)

`E_perp` is the magnitude of `grad(psi) . n_hat` where `n_hat` is the
outward semiconductor normal at the Si/SiO2 interface; the sign is
chosen so that `E_perp > 0` in the inversion regime of an n-channel
MOSFET. The current Phase B implementation builds the closed form in
pure Python (UFL-compatible). The UFL dispatch reads
`interface_facet_tag` from the validated config; the FEM-side wiring
of `n_hat` and the per-cell indicator that disables the surface
terms in the bulk lands in Phase C alongside the runner pass-through
of `facet_tags`.

Unit conversions: the JSON parameters live in device-physics units
(cm/s for B; cm^(5/3) V^(-2/3) s^-1 for C; cm^2 V s^-1 for delta).
`lombardi_unit_conversions` returns a dict of dimensionless constants
ready to be wrapped in `fem.Constant` against the scaled `E_perp`
expression `grad(psi_hat) . n_hat` whose units are 1/m. The
derivation parallels `caughey_thomas_vsat_for_form` and is documented
inline in `lombardi_unit_conversions`.
"""
from __future__ import annotations

from typing import Any

# Regularization for sqrt(grad(phi).grad(phi)) at quadrature points where
# the gradient is identically zero (equilibrium initial guess). 1e-30 in
# (1/m)^2 units leaves the dimensionless ratio (mu0 * F / vsat)^beta
# below double-precision underflow.
_F_PAR_EPS_SQ: float = 1.0e-30

# Regularization for E_perp -> 0 in the Lombardi surface terms. Same
# floor as _F_PAR_EPS_SQ; applied as `E_perp**2 + _E_PERP_EPS_SQ`
# inside `lombardi_mu_sr` and as `(E_perp + sqrt(_E_PERP_EPS_SQ))`
# inside `lombardi_mu_AC` so the divergence at zero perpendicular
# field returns a numerically large (not NaN) mu_AC / mu_sr; the
# resistor-sum composition then drives mu -> mu_bulk in the bulk.
_E_PERP_EPS_SQ: float = 1.0e-30


def caughey_thomas_mu(mu0, F_par, vsat, beta):
    """
    Caughey-Thomas field-dependent mobility, dimensionless or UFL form.

        mu(F) = mu0 / (1 + (mu0 * F_par / vsat)^beta)^(1/beta)

    Parameters
    ----------
    mu0 : float | ufl expression | fem.Constant
        Low-field mobility. Pass either a Python scalar (for the unit
        tests) or a UFL-compatible quantity in scaled units (for the
        DD form builder).
    F_par : float | ufl expression
        Magnitude of the parallel field. In SI use V/m and a SI vsat;
        in scaled units use the gradient of the scaled quasi-Fermi
        potential (1/m) and pass the corresponding scaled vsat (also
        1/m, see `caughey_thomas_vsat_for_form`).
    vsat : float | fem.Constant
        Saturation velocity. Units must match F_par * mu0.
    beta : float | fem.Constant
        Caughey-Thomas exponent (Si defaults: beta_n = 2, beta_p = 1).

    Returns
    -------
    Same kind as the inputs (Python scalar or UFL expression).
    """
    ratio = (mu0 * F_par) / vsat
    return mu0 / (1.0 + ratio ** beta) ** (1.0 / beta)


def constant_mu(mu0):
    """Identity wrapper. Keeps `build_mobility()` returning the same
    shape (UFL-compatible object) for either model branch. The constant
    branch reuses the input `fem.Constant` so byte-equivalence with
    pre-M16.1 is preserved at the form-compilation layer.
    """
    return mu0


def lombardi_mu_AC(B, C, N_total, E_perp, lam, T):
    """
    Lombardi acoustic-phonon (and Coulomb-phonon) surface-mobility term.

        mu_AC = B / E_perp + C * N_total^lambda / (E_perp^(1/3) * T)

    Inputs may be Python scalars (for unit tests) or UFL-compatible
    expressions in scaled-form units (for the FEM builder). Returns
    the same kind. The two summands carry the "acoustic phonon" and
    "Coulomb phonon" contributions of the Lombardi 1988 model; the
    combined mu_AC reduces to B / E_perp at low doping and approaches
    the C * N_total^lambda term as the doping increases.
    """
    return B / E_perp + C * N_total ** lam / (E_perp ** (1.0 / 3.0) * T)


def lombardi_mu_sr(delta, E_perp):
    """
    Lombardi surface-roughness mobility term.

        mu_sr = delta / E_perp^2

    Diverges as E_perp -> 0; in the resistor-sum composition this
    sends 1/mu_sr -> 0 so the composite reduces to the parallel
    combination of mu_bulk and mu_AC away from the interface.
    Inputs may be Python scalars or UFL expressions.
    """
    return delta / (E_perp ** 2)


def lombardi_compose(mu_bulk, mu_AC, mu_sr):
    """
    Resistor-sum composition of the Lombardi composite mobility.

        1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr

    Limits:
        mu_AC, mu_sr -> +inf:  mu -> mu_bulk
        mu_bulk -> +inf:       mu -> 1 / (1/mu_AC + 1/mu_sr)
                                  = mu_AC * mu_sr / (mu_AC + mu_sr)
                                  ~ min(mu_AC, mu_sr) when one term
                                    dominates
    """
    return 1.0 / (1.0 / mu_bulk + 1.0 / mu_AC + 1.0 / mu_sr)


def lombardi_unit_conversions(lombardi_cfg: dict[str, Any], sc) -> dict[str, float]:
    """
    Convert the JSON Lombardi parameters from device-physics units
    (cm/s for B; cm^(5/3) V^(-2/3) s^-1 for C; cm^2 V s^-1 for delta)
    into the dimensionless ratios consumed by the UFL builder against
    a scaled `E_perp_for_form = grad(psi_hat) . n_hat` (units 1/m).

    Derivation (parallels `caughey_thomas_vsat_for_form`):

        Physical perpendicular field
            E_perp_SI = V_t * E_perp_for_form              [V/m]

        Acoustic-phonon term, SI:
            mu_AC_SI = B_SI / E_perp_SI
                     + C_SI * N_total_SI^lam / (E_perp_SI^(1/3) * T)
            with B_SI = B_cm * 1e-2  [m/s]
                 C_SI = C_cm * (1e-2)^(5/3) [m^(5/3) V^(-2/3) s^-1]
                 N_total_SI = N_total_cm * 1e6  [m^-3]

        Dimensionless mu_AC_hat = mu_AC_SI / sc.mu0:

            B_for_form = B_SI / (V_t * sc.mu0)              [1/m]
                       (so B_for_form / E_perp_for_form is
                        dimensionless and equals mu_AC_hat
                        contribution from the B term)
            C_for_form = C_SI * N_total_SI^lam
                          / (V_t^(1/3) * T * sc.mu0)        [1/m^(1/3)]
                       (so C_for_form / E_perp_for_form^(1/3) is
                        dimensionless)

        Surface-roughness term, SI:
            mu_sr_SI = delta_SI / E_perp_SI^2
            with delta_SI = delta_cm * (1e-2)^2  [m^2 V s^-1]

            delta_for_form = delta_SI
                              / (V_t^2 * sc.mu0)            [1/m^2]

    `lombardi_cfg` is the validated `physics.mobility.lombardi`
    sub-dict; `sc` is the `semi.scaling.Scaling` object. `N_total` is
    not folded in here (it is a per-cell field on the function space
    and the UFL builder substitutes it directly); the C_for_form
    helper below carries only the geometric/scaling factor and the
    per-cell N_total contribution stays in the UFL expression.

    Returns a dict with keys `B_n_for_form`, `B_p_for_form`,
    `C_n_for_form`, `C_p_for_form`, `lambda_n`, `lambda_p`,
    `delta_n_for_form`, `delta_p_for_form`. The `_for_form` values
    are intended to be wrapped in `fem.Constant`.
    """
    cm_to_m = 1.0e-2
    cm5_3_to_m5_3 = cm_to_m ** (5.0 / 3.0)
    cm2_to_m2 = cm_to_m ** 2
    V_t = sc.V0
    mu0_ref = sc.mu0

    def _B(B_cm: float) -> float:
        B_SI = float(B_cm) * cm_to_m
        return B_SI / (V_t * mu0_ref)

    def _C_geom(C_cm: float) -> float:
        # geometric / scaling factor on the Coulomb-phonon term;
        # the per-cell N_total^lambda contribution lives in UFL.
        C_SI = float(C_cm) * cm5_3_to_m5_3
        return C_SI / (V_t ** (1.0 / 3.0) * mu0_ref)

    def _delta(delta_cm: float) -> float:
        delta_SI = float(delta_cm) * cm2_to_m2
        return delta_SI / (V_t ** 2 * mu0_ref)

    return {
        "B_n_for_form": _B(lombardi_cfg.get("B_n", 4.75e7)),
        "B_p_for_form": _B(lombardi_cfg.get("B_p", 9.93e6)),
        "C_n_for_form": _C_geom(lombardi_cfg.get("C_n", 1.74e5)),
        "C_p_for_form": _C_geom(lombardi_cfg.get("C_p", 8.84e5)),
        "lambda_n": float(lombardi_cfg.get("lambda_n", 0.125)),
        "lambda_p": float(lombardi_cfg.get("lambda_p", 0.0317)),
        "delta_n_for_form": _delta(lombardi_cfg.get("delta_n", 5.82e14)),
        "delta_p_for_form": _delta(lombardi_cfg.get("delta_p", 2.0546e14)),
    }


def caughey_thomas_vsat_for_form(vsat_cm_per_s: float, sc) -> float:
    """
    Convert a saturation velocity in cm/s to the dimensionless form
    consumed by `caughey_thomas_mu` against `F_par = sqrt(grad(phi_hat)
    . grad(phi_hat))` in 1/m (the gradient of the scaled quasi-Fermi
    potential under the SemiSim mesh-in-meters convention).

    Derivation (cite docs/PHYSICS.md section 2.5):
        physical field    F = V_t * grad(phi_hat)          [V/m]
        physical velocity v = mu0_SI * F                   [m/s]
        dimensionless arg (mu0_SI * F) / vsat_SI
                      = (mu0_hat * mu0_ref * V_t * grad(phi_hat)) / vsat_SI
                      = mu0_hat * grad(phi_hat) / [vsat_SI / (mu0_ref * V_t)]
        so the scaled vsat that closes the dimensionless ratio is
            vsat_for_form = vsat_SI / (sc.mu0 * sc.V0)     [1/m]

    `vsat_cm_per_s` is the JSON value in cm/s; convert to m/s via 1e-2
    before dividing.
    """
    vsat_SI = float(vsat_cm_per_s) * 1.0e-2
    return vsat_SI / (sc.mu0 * sc.V0)


def _build_constant(phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0):
    """Constant-mobility branch. Bit-identical to pre-M16.1: returns
    `fem.Constant` of the low-field ratios on each carrier's mesh."""
    from dolfinx import fem
    from petsc4py import PETSc

    msh_n = phi_n.function_space.mesh
    msh_p = phi_p.function_space.mesh
    mu_n_expr = fem.Constant(msh_n, PETSc.ScalarType(mu_n_over_mu0))
    mu_p_expr = fem.Constant(msh_p, PETSc.ScalarType(mu_p_over_mu0))
    return mu_n_expr, mu_p_expr, "constant"


def _build_caughey_thomas(
    mob: dict[str, Any], phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc
):
    """Caughey-Thomas closed-form velocity-saturation branch.
    Bit-identical to v0.17.0 on diode_velsat_1d (M16.1 anchor)."""
    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    msh_n = phi_n.function_space.mesh
    msh_p = phi_p.function_space.mesh

    vsat_n = caughey_thomas_vsat_for_form(mob.get("vsat_n", 1.0e7), sc)
    vsat_p = caughey_thomas_vsat_for_form(mob.get("vsat_p", 8.0e6), sc)
    beta_n = float(mob.get("beta_n", 2.0))
    beta_p = float(mob.get("beta_p", 1.0))

    mu0_n = fem.Constant(msh_n, PETSc.ScalarType(mu_n_over_mu0))
    mu0_p = fem.Constant(msh_p, PETSc.ScalarType(mu_p_over_mu0))
    vsat_n_c = fem.Constant(msh_n, PETSc.ScalarType(vsat_n))
    vsat_p_c = fem.Constant(msh_p, PETSc.ScalarType(vsat_p))
    beta_n_c = fem.Constant(msh_n, PETSc.ScalarType(beta_n))
    beta_p_c = fem.Constant(msh_p, PETSc.ScalarType(beta_p))
    eps_n = fem.Constant(msh_n, PETSc.ScalarType(_F_PAR_EPS_SQ))
    eps_p = fem.Constant(msh_p, PETSc.ScalarType(_F_PAR_EPS_SQ))

    grad_phi_n = ufl.grad(phi_n)
    grad_phi_p = ufl.grad(phi_p)
    F_par_n = ufl.sqrt(ufl.dot(grad_phi_n, grad_phi_n) + eps_n)
    F_par_p = ufl.sqrt(ufl.dot(grad_phi_p, grad_phi_p) + eps_p)

    mu_n_expr = caughey_thomas_mu(mu0_n, F_par_n, vsat_n_c, beta_n_c)
    mu_p_expr = caughey_thomas_mu(mu0_p, F_par_p, vsat_p_c, beta_p_c)
    return mu_n_expr, mu_p_expr, "caughey_thomas"


def _build_lombardi(
    mob: dict[str, Any],
    psi,
    phi_n,
    phi_p,
    mu_n_over_mu0,
    mu_p_over_mu0,
    sc,
    *,
    facet_tags,
    N_total_hat,
):
    """
    Lombardi composite-mobility branch (M16.2).

    Builds the resistor-sum composite
        1/mu = 1/mu_bulk + 1/mu_AC + 1/mu_sr
    where mu_bulk is the constant or caughey_thomas branch (selected
    by `bulk_model`), mu_AC is the Lombardi acoustic-phonon term, and
    mu_sr is the surface-roughness term. The closed-form math comes
    from the public `lombardi_mu_AC`, `lombardi_mu_sr`, and
    `lombardi_compose` helpers; this function only wraps them in UFL
    against scaled fields.

    The perpendicular field is approximated as
        E_perp_for_form = sqrt((grad(psi) . n_hat)^2 + eps)
    where n_hat is a unit vector along the last geometric axis (axis
    `dim - 1` of the parent mesh). The Lombardi reference geometry has
    the Si/SiO2 interface at the top of the silicon body (the highest
    coordinate of the depth axis), so the last-axis convention covers
    every benchmark currently in the tree (mosfet_2d 2D, diode_1d in
    the unused-Lombardi default). Cells far from the interface where
    grad(psi) . n_hat -> 0 see mu_AC, mu_sr -> +inf via the
    regularized division, and the resistor sum collapses to mu_bulk
    in the bulk; the Lombardi correction is naturally localized to
    the surface region without an explicit cell indicator.
    """
    if facet_tags is None:
        raise ValueError(
            "physics.mobility.model='lombardi' requires the runner to "
            "pass facet_tags=ft so the FEM builder can locate the "
            "interface_facet_tag on the mesh; bias_sweep and the MMS "
            "Variant E harness do this from Phase C of M16.2 onward"
        )
    if N_total_hat is None:
        raise ValueError(
            "physics.mobility.model='lombardi' requires the runner to "
            "pass N_total_hat (scaled total ionized impurity field "
            "abs(N_net_hat) for uncompensated profiles)"
        )
    if psi is None:
        raise ValueError(
            "physics.mobility.model='lombardi' requires the runner to "
            "pass psi=spaces.psi so the perpendicular electrostatic "
            "field can be evaluated on the parent mesh"
        )

    interface_tag = mob.get("interface_facet_tag")
    if interface_tag is None:
        raise ValueError(
            "physics.mobility.model='lombardi' requires "
            "physics.mobility.interface_facet_tag (an integer mesh "
            "facet tag identifying the Si/SiO2 interface)"
        )

    import ufl
    from dolfinx import fem
    from petsc4py import PETSc

    bulk_model = mob.get("bulk_model", "constant")
    if bulk_model == "constant":
        mu_bulk_n, mu_bulk_p, _ = _build_constant(
            phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0
        )
    elif bulk_model == "caughey_thomas":
        mu_bulk_n, mu_bulk_p, _ = _build_caughey_thomas(
            mob, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc
        )
    else:
        raise ValueError(
            f"physics.mobility.bulk_model={bulk_model!r} not supported "
            "under the lombardi composite; expected "
            "'constant' or 'caughey_thomas'"
        )

    msh = psi.function_space.mesh
    gdim = msh.geometry.dim
    # Unit vector along the last geometric axis (the depth axis of the
    # mosfet_2d benchmark; matches the manufactured-interface
    # orientation in MMS Variant E).
    normal_axis = int(mob.get("interface_normal_axis", gdim - 1))
    if not (0 <= normal_axis < gdim):
        raise ValueError(
            f"interface_normal_axis={normal_axis} is out of bounds for "
            f"a mesh with geometric dimension {gdim}"
        )
    n_hat_components = [0.0] * gdim
    n_hat_components[normal_axis] = 1.0
    n_hat = ufl.as_vector([
        fem.Constant(msh, PETSc.ScalarType(c)) for c in n_hat_components
    ])

    cfg_lombardi = mob.get("lombardi", {}) or {}
    conv = lombardi_unit_conversions(cfg_lombardi, sc)

    T_K = float(mob.get("temperature_K", 300.0))
    T_const = fem.Constant(msh, PETSc.ScalarType(T_K))
    eps_perp = fem.Constant(msh, PETSc.ScalarType(_E_PERP_EPS_SQ))

    grad_psi = ufl.grad(psi)
    E_perp_signed = ufl.dot(grad_psi, n_hat)
    E_perp = ufl.sqrt(E_perp_signed * E_perp_signed + eps_perp)

    def _surface_terms(carrier: str):
        if carrier == "n":
            B = conv["B_n_for_form"]
            C_geom = conv["C_n_for_form"]
            lam = conv["lambda_n"]
            delta = conv["delta_n_for_form"]
        else:
            B = conv["B_p_for_form"]
            C_geom = conv["C_p_for_form"]
            lam = conv["lambda_p"]
            delta = conv["delta_p_for_form"]
        # Fold sc.C0**lam into C_geom so the UFL expression sees the
        # bare dimensionless N_total_hat^lam (no extra material
        # constants in the form).
        C_eff = C_geom * (sc.C0 ** lam)
        B_c = fem.Constant(msh, PETSc.ScalarType(B))
        C_c = fem.Constant(msh, PETSc.ScalarType(C_eff))
        delta_c = fem.Constant(msh, PETSc.ScalarType(delta))
        lam_c = fem.Constant(msh, PETSc.ScalarType(lam))
        mu_AC = lombardi_mu_AC(B_c, C_c, N_total_hat, E_perp, lam_c, T_const)
        mu_sr = lombardi_mu_sr(delta_c, E_perp)
        return mu_AC, mu_sr

    mu_AC_n, mu_sr_n = _surface_terms("n")
    mu_AC_p, mu_sr_p = _surface_terms("p")

    mu_n_expr = lombardi_compose(mu_bulk_n, mu_AC_n, mu_sr_n)
    mu_p_expr = lombardi_compose(mu_bulk_p, mu_AC_p, mu_sr_p)
    return mu_n_expr, mu_p_expr, "lombardi"


def build_mobility_expressions(
    mobility_cfg: dict[str, Any] | None,
    phi_n,
    phi_p,
    mu_n_over_mu0: float,
    mu_p_over_mu0: float,
    sc,
    *,
    psi=None,
    facet_tags=None,
    N_total_hat=None,
) -> tuple[Any, Any, str]:
    """
    Dispatch on `mobility_cfg["model"]` and return UFL-compatible
    mobility expressions substitutable for the legacy
    `fem.Constant(mu_n_over_mu0)` / `fem.Constant(mu_p_over_mu0)` in
    the DD residual builders.

    Parameters
    ----------
    mobility_cfg : dict or None
        The validated `cfg["physics"]["mobility"]` sub-dict, or `None`
        which is equivalent to `{"model": "constant"}` (preserves
        pre-M16.1 byte-identity).
    phi_n, phi_p : dolfinx.fem.Function
        The carrier-specific scaled quasi-Fermi potentials. The
        Caughey-Thomas branch reads `ufl.grad(phi_n)` for the electron
        continuity row and `ufl.grad(phi_p)` for the hole continuity
        row (ADR 0004 Slotboom flux form).
    mu_n_over_mu0, mu_p_over_mu0 : float
        Low-field mobility ratios (mu / sc.mu0) already computed by
        the caller.
    sc : semi.scaling.Scaling
        Scaling object; used by `caughey_thomas_vsat_for_form` and
        `lombardi_unit_conversions`.
    psi : dolfinx.fem.Function, optional
        Scaled electrostatic potential on the parent mesh; required
        by the lombardi branch (the perpendicular field at the
        interface is the gradient of psi). Ignored by the constant
        and caughey_thomas branches.
    facet_tags : dolfinx.mesh.MeshTags, optional
        Parent-mesh facet tags carrying the interface label; required
        by the lombardi branch and ignored by the others. The runner
        already constructs this for BC resolution and forwards it
        through `build_dd_block_residual`.
    N_total_hat : ufl expression, optional
        Scaled total ionized impurity concentration field
        N_total / sc.C0. For uncompensated profiles (every shipped
        benchmark) `abs(N_net_hat)` is the right value; the runner
        hands this in via the same Function used for the Poisson RHS.

    Returns
    -------
    (mu_n_expr, mu_p_expr, model_name) : tuple
        UFL or fem.Constant objects substitutable for the legacy
        `fem.Constant(mu_n_over_mu0)` / `fem.Constant(mu_p_over_mu0)`
        in the DD residual. `model_name` echoes the dispatch for
        diagnostics.

    Branches
    --------
    constant (default):
        Returns `(fem.Constant(msh, mu_n_over_mu0),
                  fem.Constant(msh, mu_p_over_mu0))`. Bit-identical to
        pre-M16.1 since the existing DD form already wraps the inputs
        in fem.Constants of the same scalar values.
    caughey_thomas:
        Builds `caughey_thomas_mu(mu0_const, F_par, vsat_for_form,
        beta)` for each carrier where `F_par` is
        `sqrt(grad(phi) . grad(phi) + eps)` of the carrier-specific
        scaled quasi-Fermi potential. Bit-identical to v0.17.0.
    lombardi:
        Resistor-sum composite of the bulk branch (selected by
        `bulk_model`) with the Lombardi acoustic-phonon and
        surface-roughness terms. Requires `psi`, `facet_tags`, and
        `N_total_hat` to be passed; raises ValueError otherwise.
    """
    mob = mobility_cfg or {}
    model = mob.get("model", "constant")

    if model == "constant":
        return _build_constant(phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0)

    if model == "caughey_thomas":
        return _build_caughey_thomas(
            mob, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc
        )

    if model == "lombardi":
        return _build_lombardi(
            mob, psi, phi_n, phi_p, mu_n_over_mu0, mu_p_over_mu0, sc,
            facet_tags=facet_tags, N_total_hat=N_total_hat,
        )

    raise ValueError(
        f"physics.mobility.model={model!r} not supported; "
        "expected 'constant', 'caughey_thomas', or 'lombardi'"
    )
