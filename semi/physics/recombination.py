"""
Net recombination kernels.

Shockley-Read-Hall recombination (always available).

Dimensional form (PHYSICS.md 1.4):

    R_SRH = (n p - n_i^2) / (tau_p (n + n1) + tau_n (p + p1))

with n1 = n_i exp(E_t / V_t), p1 = n_i exp(-E_t / V_t), and E_t the
trap energy measured from the intrinsic level (zero for a mid-gap
trap).

Scaled form. Substituting n = C0 n_hat, p = C0 p_hat, n_i = C0 n_i_hat,
and dividing numerator and denominator by C0:

    R / (C0 / t0)  =  (n_hat p_hat - n_i_hat^2)
                      / ( tau_p_hat (n_hat + n1_hat)
                        + tau_n_hat (p_hat + p1_hat) )

where the scaled lifetimes are tau_hat = tau / t0. The return value is
scaled by the density rate C0/t0, which is what the scaled continuity
equation expects.

The pure-Python helper :func:`srh_rate_np` evaluates the same formula
on numpy arrays and is used by the unit tests. The UFL helper
:func:`srh_rate` consumes UFL expressions (densities built via
Slotboom) and returns a UFL expression suitable for a form builder.

Auger recombination (M16.3).

Dimensional form (PHYSICS.md 1.4, Auger paragraph):

    R_Auger = (C_n n + C_p p) (n p - n_i^2)

Two-carrier band-to-band Auger; same drive (n p - n_i^2) as SRH but a
density-weighted prefactor (C_n n + C_p p) instead of the SRH lifetime
denominator. Units: SI `C_n`, `C_p` in m^6/s; the JSON contract uses
cm^6/s (Si Dziewior-Schmid: C_n = 2.8e-31, C_p = 9.9e-32). At high
injection (n ~ p >> n_i, N) the kernel limits to

    R_Auger -> (C_n + C_p) n^3,

cubic in carrier density and dominant over the bilinear-over-linear
SRH form. Auger is additive to SRH: the form builder evaluates
R = R_SRH + R_Auger when the auger flag is on (see ADR 0001 for the
JSON contract; see ADR 0002 for the scaling convention).

Scaled form. Substituting n = C0 n_hat into the Auger expression and
dividing by C0/t0 (the density rate the scaled continuity equation
expects):

    R / (C0/t0) = C0^2 t0 (C_n n_hat + C_p p_hat) (n_hat p_hat
                                                    - n_i_hat^2)

so the dimensionless Auger coefficient is

    C_hat = C_SI * C0^2 * t0

with C_SI in m^6/s. SI conversion from cm^6/s to m^6/s is 1e-12
(cm^6 = 1e-12 m^6). The pure-Python helper :func:`scaled_auger_C`
performs the dimensional reduction and returns the dimensionless
ratio; runners convert cm^6/s -> m^6/s before calling.
"""
from __future__ import annotations

import math

import numpy as np


def srh_rate(n_hat, p_hat, n_i_hat, tau_n_hat, tau_p_hat, E_t_over_Vt=0.0):
    """
    UFL expression for scaled SRH recombination rate.

    Parameters
    ----------
    n_hat, p_hat : UFL expressions
        Scaled carrier densities (for example built with
        :mod:`semi.physics.slotboom`).
    n_i_hat : UFL expression or float
        Scaled intrinsic density.
    tau_n_hat, tau_p_hat : UFL expression or float
        Scaled electron and hole lifetimes (tau / t0).
    E_t_over_Vt : float
        Trap level referenced to the intrinsic Fermi level, divided by
        V_t. Zero for a mid-gap trap.

    Returns
    -------
    UFL expression
        Scaled recombination rate (per C0/t0).
    """
    import ufl
    n1_hat = n_i_hat * math.exp(E_t_over_Vt)
    p1_hat = n_i_hat * math.exp(-E_t_over_Vt)
    num = n_hat * p_hat - n_i_hat * n_i_hat
    den = tau_p_hat * (n_hat + n1_hat) + tau_n_hat * (p_hat + p1_hat)
    # Explicit UFL call clarifies evaluation when n_i_hat is a fem.Constant
    return ufl.as_ufl(1.0) * num / den


def srh_rate_np(n, p, n_i, tau_n, tau_p, E_t_over_Vt=0.0):
    """
    NumPy counterpart of :func:`srh_rate`.

    Accepts any units (dimensional or scaled) so long as all six
    arguments are consistent. Returns the rate in the matching unit
    system (for example m^-3 s^-1 if the inputs are SI).
    """
    n1 = n_i * np.exp(E_t_over_Vt)
    p1 = n_i * np.exp(-E_t_over_Vt)
    num = n * p - n_i ** 2
    den = tau_p * (n + n1) + tau_n * (p + p1)
    return num / den


def scaled_tau(tau_seconds, t0_seconds):
    """Dimensionless lifetime, tau_hat = tau / t0."""
    return float(tau_seconds) / float(t0_seconds)


# ---------------------------------------------------------------------------
# M16.3 Auger recombination.
# ---------------------------------------------------------------------------


def auger_rate(n_hat, p_hat, n_i_hat, C_n_hat, C_p_hat):
    """
    UFL expression for the scaled Auger recombination rate

        R_hat_Auger = (C_n_hat n_hat + C_p_hat p_hat)
                      * (n_hat p_hat - n_i_hat^2).

    The factor `(n_hat p_hat - n_i_hat^2)` is shared with the SRH
    kernel; an inline assembly in the form builder reuses it via UFL's
    automatic common subexpression elimination (see
    `semi/physics/drift_diffusion.py`).

    Parameters
    ----------
    n_hat, p_hat : UFL expressions
        Scaled carrier densities (Slotboom).
    n_i_hat : UFL expression or float
        Scaled intrinsic density.
    C_n_hat, C_p_hat : UFL expression or float
        Dimensionless Auger coefficients (C_SI * C0^2 * t0); see
        :func:`scaled_auger_C`.

    Returns
    -------
    UFL expression
        Scaled Auger recombination rate (per C0/t0).
    """
    import ufl
    np_minus_nieq = n_hat * p_hat - n_i_hat * n_i_hat
    prefactor = C_n_hat * n_hat + C_p_hat * p_hat
    return ufl.as_ufl(1.0) * prefactor * np_minus_nieq


def auger_rate_np(n, p, n_i, C_n, C_p):
    """
    NumPy counterpart of :func:`auger_rate`.

    Accepts any consistent unit system (dimensional or scaled). With
    SI inputs (n, p, n_i in m^-3 and C_n, C_p in m^6/s), the return
    value is in m^-3 s^-1.
    """
    return (C_n * n + C_p * p) * (n * p - n_i ** 2)


def scaled_auger_C(C_si, C0, t0):
    """
    Convert a dimensional Auger coefficient (m^6/s) to the
    dimensionless ratio consumed by the scaled form:

        C_hat = C_SI * C0^2 * t0

    where C0 has units m^-3 and t0 has units s. Callers reading from
    the JSON contract must first convert cm^6/s to m^6/s (factor 1e-12)
    before invoking this helper; see ADR 0002.
    """
    return float(C_si) * float(C0) ** 2 * float(t0)


# ---------------------------------------------------------------------------
# M16.6 Kane band-to-band tunneling and Hurkx trap-assisted tunneling.
# ---------------------------------------------------------------------------
#
# Kane band-to-band tunneling (Kane 1959; Sze 3rd ed section 8.4).
#
# Dimensional form:
#
#     G_BBT = A_kane * |E|^2 / sqrt(E_g)
#                  * exp(-B_kane * E_g^(3/2) / |E|)
#
# where A_kane has units cm^-1 s^-1 V^-2 (Si default 4e14), B_kane
# has units V/cm (Si default 1.9e7), |E| is the magnitude of the
# local electric field in V/cm, and E_g is the band gap in eV. This
# is a pure generation term: under the R = U - G convention used in
# semi/physics/drift_diffusion.py (R contributes positively to the
# electron continuity row and negatively to the hole continuity row),
# BBT contributes -G_BBT to R so that the same scalar appears with
# opposite signs on the two continuity rows (one tunneled
# electron-hole pair appears).
#
# Hurkx trap-assisted tunneling (Hurkx 1992).
#
# Dimensional form:
#
#     R_TAT = (1 + Gamma(F)) * R_SRH
#     Gamma(F) = 2 sqrt(3 pi) * (F / F_kT)^(alpha - 1)
#                              * exp((F / F_kT)^2)
#
# where F is the local field magnitude, F_kT the Hurkx characteristic
# field (Si default 1.4e7 V/cm), and alpha the Hurkx exponent (Si
# default 2.0). Gamma vanishes as F -> 0 (alpha > 1) so the TAT
# branch is bit-identical to v0.21.0 in the low-field bulk; it grows
# super-exponentially in the high-field depletion region where
# trap-assisted tunneling enhances SRH.
#
# Scaled form. The form builder in semi/physics/drift_diffusion.py
# evaluates the field magnitude as |grad(psi_hat)|, which is
# dimensionless (psi_hat is V_t-scaled and the gradient is taken in
# the L_0-scaled coordinate). The dimensional field is then
#
#     |E_phys|_V/m = (V_0 / L_0) * |grad(psi_hat)|
#
# To keep the UFL helpers free of unit baggage we expose
# `scaled_kane_coefficients` / `scaled_hurkx_F_kT` that produce the
# dimensionless ratios consumed by the UFL helpers. The closed-form
# closures below match the manufactured-source evaluations in
# semi/verification/mms_dd.py Variant H exactly.


def bbt_rate(E_field_magnitude, E_g_hat, A_kane_hat, B_kane_hat):
    """
    UFL expression for the scaled Kane band-to-band tunneling
    generation rate

        G_hat = A_kane_hat * |E_hat|^2 / sqrt(E_g_hat)
                         * exp(-B_kane_hat * E_g_hat^(3/2) / |E_hat|)

    Sign convention: positive (this is a generation term; the caller
    subtracts it from the recombination R when assembling the
    continuity rows).

    Parameters
    ----------
    E_field_magnitude : UFL expression
        Scaled field magnitude (|grad(psi_hat)|).
    E_g_hat : UFL expression or float
        Band gap in scaled units; see :func:`scaled_E_g`.
    A_kane_hat, B_kane_hat : UFL expression or float
        Dimensionless Kane coefficients; see
        :func:`scaled_kane_coefficients`.

    Returns
    -------
    UFL expression
        Scaled BBT generation rate (per C0/t0).
    """
    import ufl
    return (
        ufl.as_ufl(1.0)
        * A_kane_hat
        * E_field_magnitude * E_field_magnitude
        / ufl.sqrt(E_g_hat)
        * ufl.exp(
            -B_kane_hat * (E_g_hat ** 1.5) / E_field_magnitude
        )
    )


def bbt_rate_np(E_field_magnitude, E_g, A_kane, B_kane):
    """
    NumPy counterpart of :func:`bbt_rate`.

    Accepts any consistent unit system. With JSON-contract inputs
    (`|E|` in V/cm, `E_g` in eV, `A_kane` in cm^-1 s^-1 V^-2,
    `B_kane` in V/cm), the return value is in cm^-3 s^-1.
    """
    return (
        A_kane
        * E_field_magnitude ** 2
        / np.sqrt(E_g)
        * np.exp(-B_kane * (E_g ** 1.5) / E_field_magnitude)
    )


def hurkx_gamma(F, F_kT_hat, alpha):
    """
    UFL expression for the Hurkx field-enhancement factor Gamma(F)

        Gamma = 2 sqrt(3 pi) * (F / F_kT_hat)^(alpha - 1)
                             * exp((F / F_kT_hat)^2)

    Returns a dimensionless UFL expression; the caller multiplies the
    SRH rate by `(1 + Gamma)` to enable trap-assisted tunneling.

    Parameters
    ----------
    F : UFL expression
        Scaled field magnitude (|grad(psi_hat)|).
    F_kT_hat : UFL expression or float
        Hurkx characteristic field in scaled units; see
        :func:`scaled_hurkx_F_kT`.
    alpha : float
        Hurkx exponent (typically 2.0).
    """
    import ufl
    PREF = 2.0 * math.sqrt(3.0 * math.pi)
    ratio = F / F_kT_hat
    return (
        ufl.as_ufl(PREF)
        * (ratio ** (alpha - 1.0))
        * ufl.exp(ratio * ratio)
    )


def hurkx_gamma_np(F, F_kT, alpha):
    """
    NumPy counterpart of :func:`hurkx_gamma`. Accepts any consistent
    unit system; with JSON-contract inputs (`F` and `F_kT` both in
    V/cm) the return value is dimensionless.
    """
    pref = 2.0 * math.sqrt(3.0 * math.pi)
    ratio = F / F_kT
    return pref * (ratio ** (alpha - 1.0)) * np.exp(ratio * ratio)


def scaled_kane_coefficients(A_kane_cm, B_kane_cm, sc):
    """
    Convert the JSON Kane coefficients (cm-based units) into the
    dimensionless ratios consumed by :func:`bbt_rate`.

    Substituting `|E_phys| = (V_0 / L_0) * |grad(psi_hat)|` and
    `E_g_phys = V_0 * E_g_hat` into the dimensional Kane formula and
    dividing by `C_0 / t_0` yields

        A_kane_hat = A_SI * V_0^(3/2) * t_0 / (L_0^2 * C_0)
        B_kane_hat = B_SI * L_0 * sqrt(V_0)

    with SI conversions `A_SI [m^-1 s^-1 V^-2] = A_cm * 100` and
    `B_SI [V/m] = B_cm * 100`.

    Parameters
    ----------
    A_kane_cm : float
        Kane prefactor in cm^-1 s^-1 V^-2 (Si default 4.0e14).
    B_kane_cm : float
        Kane exponent coefficient in V/cm (Si default 1.9e7).
    sc : semi.scaling.Scaling
        Scaling object; provides V0, L0, C0, t0.

    Returns
    -------
    dict
        ``{"A_kane_hat": ..., "B_kane_hat": ...}``.
    """
    A_SI = float(A_kane_cm) * 100.0
    B_SI = float(B_kane_cm) * 100.0
    A_hat = A_SI * (sc.V0 ** 1.5) * sc.t0 / (sc.L0 ** 2 * sc.C0)
    B_hat = B_SI * sc.L0 * math.sqrt(sc.V0)
    return {"A_kane_hat": A_hat, "B_kane_hat": B_hat}


def scaled_hurkx_F_kT(F_kT_cm, sc):
    """
    Convert the JSON Hurkx characteristic field (V/cm) into the
    dimensionless ratio consumed by :func:`hurkx_gamma`. The scaled
    field magnitude is `|grad(psi_hat)|` (dimensionless), which
    corresponds to a physical field `(V_0 / L_0) * |grad(psi_hat)|`
    in V/m. Setting `F_kT_hat = F_kT_SI * L_0 / V_0` makes the
    `(F / F_kT_hat)` ratio in the Hurkx form dimensionless.
    """
    F_SI = float(F_kT_cm) * 100.0
    return F_SI * sc.L0 / sc.V0


def scaled_E_g(E_g_eV, sc):
    """
    Convert the band gap (eV) into the scaled-energy units consumed
    by :func:`bbt_rate`: `E_g_hat = E_g_eV / V_0` (V_0 is the
    thermal voltage in V; the eV/V identity makes this a pure
    division).
    """
    return float(E_g_eV) / sc.V0
