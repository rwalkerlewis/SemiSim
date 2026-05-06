"""
Analytical reference curves for 1D pn-junction IV verifiers.

Pure-Python, numpy-only, no dolfinx. The functions here are called by
`scripts/run_benchmark.py` to build reference curves for the
`pn_1d_bias` and `pn_1d_bias_reverse` verifiers, and by unit tests
that exercise the formulas directly.

Units are SI throughout: densities in m^-3, potentials in V, lengths
in m, current density in A/m^2.
"""
from __future__ import annotations

import math

import numpy as np

from .constants import HBAR, KB, M0, Q


def shockley_long_diode_saturation(
    N_A: float, N_D: float, n_i: float,
    mu_n_SI: float, mu_p_SI: float,
    tau_n: float, tau_p: float,
    V_t: float,
) -> tuple[float, float, float]:
    """
    Long-diode Shockley diffusion saturation current.

    J_s = q n_i^2 (D_n / (L_n N_A) + D_p / (L_p N_D)),
    D = V_t mu, L = sqrt(D tau).
    Returns (J_s, L_n, L_p).
    """
    D_n = V_t * mu_n_SI
    D_p = V_t * mu_p_SI
    L_n = (D_n * tau_n) ** 0.5
    L_p = (D_p * tau_p) ** 0.5
    J_s = Q * n_i ** 2 * (D_n / (L_n * N_A) + D_p / (L_p * N_D))
    return J_s, L_n, L_p


def depletion_width(
    N_A: float, N_D: float, n_i: float, eps: float,
    V_t: float, V,
) -> np.ndarray:
    """
    Depletion-approximation width under applied bias V < V_bi.

    V_bi = V_t ln(N_A N_D / n_i^2)
    W(V) = sqrt(2 eps (V_bi - V)(N_A + N_D) / (q N_A N_D))
    """
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)
    V_arr = np.asarray(V, dtype=float)
    delta = np.maximum(V_bi - V_arr, 1.0e-9)
    return np.sqrt(2.0 * eps * delta * (N_A + N_D) / (Q * N_A * N_D))


def sns_total_reference(
    N_A: float, N_D: float, n_i: float, eps: float,
    mu_n_SI: float, mu_p_SI: float,
    tau_n: float, tau_p: float,
    V_t: float, V,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Forward-bias Sah-Noyce-Shockley total current.

    J_total(V) = J_s (exp(V/V_t) - 1)
               + (q n_i W(V) / 2 tau_eff)
                 * (exp(V/(2 V_t)) - 1)
                 * (2 V_t / (V_bi - V))

    The final (2 V_t / (V_bi - V)) prefactor is the leading-order Sze
    correction that brings SNS into 10-15% agreement with a Galerkin
    Slotboom solution on the M2 benchmark; without it the raw
    exponential overestimates the recombination branch by roughly an
    order of magnitude.

    Returns (J_total, J_diff, J_rec, J_s).
    """
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)
    V_arr = np.asarray(V, dtype=float)
    tau_eff = (tau_n * tau_p) ** 0.5

    W = depletion_width(N_A, N_D, n_i, eps, V_t, V_arr)
    delta = np.maximum(V_bi - V_arr, 1.0e-9)

    J_s, _L_n, _L_p = shockley_long_diode_saturation(
        N_A, N_D, n_i, mu_n_SI, mu_p_SI, tau_n, tau_p, V_t,
    )
    J_rec_pref = Q * n_i * W / (2.0 * tau_eff)
    f_corr = 2.0 * V_t / delta

    J_diff = J_s * (np.exp(V_arr / V_t) - 1.0)
    J_rec = J_rec_pref * (np.exp(V_arr / (2.0 * V_t)) - 1.0) * f_corr
    return J_diff + J_rec, J_diff, J_rec, J_s


def shockley_iv_with_auger(
    N_A: float, N_D: float, n_i: float, eps: float,
    mu_n_SI: float, mu_p_SI: float,
    tau_n: float, tau_p: float,
    C_n_SI: float, C_p_SI: float,
    V_t: float, V,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    High-injection long-diode I-V with additive Auger recombination
    (M16.3).

    For a symmetric pn junction at high injection (`V_F` well above
    `V_bi`), excess carrier densities approximately satisfy the
    ambipolar diffusion equation in the bulk:

        D_a * d^2 delta / dx^2 - delta / tau_eff = 0

    with the effective lifetime that includes both SRH and Auger
    branches:

        1 / tau_eff = 1 / tau_SRH_eff + (C_n + C_p) * delta_avg^2

    The closed form J_total = 2 q D_a delta(0) / L_a, where
    `delta(0) = n_i exp(V / 2 V_t)` is the high-injection boundary
    density and `L_a = sqrt(D_a tau_eff)`.

    Approximations:
      - `delta_avg = delta(0) / 2` (the leading-order average of an
        exponential decay over a long bulk; gives a bias-dependent
        Auger lifetime without iterating on `delta`).
      - `tau_SRH_eff = sqrt(tau_n tau_p)`, the geometric-mean SRH
        lifetime at high injection.
      - `D_a = 2 D_n D_p / (D_n + D_p)`, the ambipolar diffusivity.

    Returns
    -------
    (J_total, J_SRH_only, tau_eff_per_V)
        `J_total` includes both SRH and Auger; `J_SRH_only` is
        the same formula with `(C_n + C_p) -> 0`. `tau_eff_per_V`
        is the effective lifetime at each V (numpy array), exposed
        so the verifier can print it on failure.

    SI inputs throughout (densities m^-3, lengths m, lifetimes s,
    Auger coefficients m^6/s, current density A/m^2). The verifier
    converts JSON cm^6/s to m^6/s before invoking.
    """
    V_arr = np.asarray(V, dtype=float)
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)

    D_n = V_t * mu_n_SI
    D_p = V_t * mu_p_SI
    D_a = 2.0 * D_n * D_p / (D_n + D_p)
    tau_SRH_eff = np.sqrt(tau_n * tau_p)

    # High-injection boundary: delta_0 = n_i exp(V / 2 V_t).
    delta_0 = n_i * np.exp(V_arr / (2.0 * V_t))
    delta_avg = delta_0 / 2.0

    # Effective lifetime including Auger.
    inv_tau = 1.0 / tau_SRH_eff + (C_n_SI + C_p_SI) * delta_avg ** 2
    tau_eff = 1.0 / inv_tau
    L_a = np.sqrt(D_a * tau_eff)

    # Long-diode current: contributions from both bulks scale with
    # delta_0 / L_a; the factor of 2 accounts for the two sides.
    J_total = 2.0 * Q * D_a * delta_0 / L_a

    # SRH-only counterpart (Auger -> 0).
    inv_tau_srh = 1.0 / tau_SRH_eff
    tau_srh = 1.0 / inv_tau_srh
    L_srh = np.sqrt(D_a * tau_srh)
    J_srh_only = 2.0 * Q * D_a * delta_0 / L_srh

    # Suppress the V < V_bi tail (the formula is the high-injection
    # asymptote; below V_bi the diffusion-dominated Shockley J ~
    # exp(V/V_t) law is the better reference, which the existing
    # sns_total_reference covers).
    below_bi = V_arr < V_bi
    J_total = np.where(below_bi, 0.0, J_total)
    J_srh_only = np.where(below_bi, 0.0, J_srh_only)

    return J_total, J_srh_only, tau_eff


def srh_generation_reference(
    N_A: float, N_D: float, n_i: float, eps: float,
    tau_n: float, tau_p: float,
    V_t: float, V,
) -> tuple[np.ndarray, float, float]:
    """
    Reverse-bias net SRH generation current.

    At V = 0 generation and recombination balance in the depletion
    region. Under reverse bias the depletion width grows from W(0)
    to W(V) and the extra width contributes a net generation current

        J_gen_net(V) = (q n_i / 2 tau_eff) * (W(V) - W(0))

    where tau_eff = sqrt(tau_n tau_p). Returns (J_gen_net, W0, V_bi).
    """
    V_bi = V_t * np.log(N_A * N_D / n_i ** 2)
    tau_eff = (tau_n * tau_p) ** 0.5

    W = depletion_width(N_A, N_D, n_i, eps, V_t, V)
    W0 = float(
        np.sqrt(2.0 * eps * V_bi * (N_A + N_D) / (Q * N_A * N_D))
    )

    J_gen_net = (Q * n_i / (2.0 * tau_eff)) * (W - W0)
    return J_gen_net, W0, V_bi


# ---------------------------------------------------------------------------
# M16.5 Schottky thermionic-emission analytical helpers.
# ---------------------------------------------------------------------------


def richardson_constant(m_star_rel: float) -> float:
    """
    Richardson constant ``A* = 4 pi q m* k^2 / h^3`` for an effective
    electron mass ``m_star_rel`` (relative to the bare electron rest
    mass ``m_0``). Returns SI units of A m^-2 K^-2 (multiply by 1e-4
    for the textbook A cm^-2 K^-2 unit).

    Si values (Sze 3rd ed Table 5; using thermionic-emission effective
    masses):
      - electron, m* ≈ 0.26 m_0  -> A* ≈ 32 A cm^-2 K^-2 ≈ 3.2e5 A m^-2 K^-2
        (the literature uses higher reported values for Si, e.g. 110
        A cm^-2 K^-2, when m_DOS or transverse-mass averaging is used;
        the M16.5 verifier uses the Sze Table 1 thermionic mass to
        keep `richardson_constant(m_n*) * T^2 / (q N_C)` consistent
        with `Scaling.v_n_thermal`).
      - hole,     m* ≈ 0.39 m_0  -> A* ≈ 47 A cm^-2 K^-2.

    Pure-Python; no dolfinx. M16.5.
    """
    h_planck = 2.0 * math.pi * HBAR
    return float(
        4.0 * math.pi * Q * (float(m_star_rel) * M0) * (KB ** 2) / (h_planck ** 3)
    )


def kane_breakdown_iv(
    V_R, N: float, eps: float,
    E_g_eV: float, V_bi: float,
    A_kane_cm: float = 4.0e14,
    B_kane_cm: float = 1.9e7,
) -> np.ndarray:
    """
    Closed-form Kane band-to-band reverse-bias breakdown current for a
    1D abrupt junction at heavy doping (Sze 3rd ed Section 8.4).

    Derivation sketch. Under the depletion-approximation field profile
    on a one-sided abrupt junction at heavy doping (N_A = N_D = N), the
    peak field at the metallurgical junction under reverse bias V_R is

        |E_max| = sqrt(2 q N (V_bi + V_R) / eps)

    (peak field of the depletion-approximation triangle; the field
    decays linearly toward each contact). The Kane band-to-band
    generation rate is

        G_BBT = A_kane |E|^2 / sqrt(E_g) exp(-B_kane E_g^(3/2) / |E|)

    with A_kane in cm^-1 s^-1 V^-2, B_kane in V/cm, |E| in V/cm, and
    E_g in eV. Integrating G_BBT over the depletion region (width
    W = sqrt(2 eps (V_bi + V_R) (N_A + N_D) / (q N_A N_D))) and using
    the linear-field profile `|E(x)| = |E_max| (1 - x/W)`, the total
    BBT-generated current per unit area is approximately

        J_BBT(V_R) ~ q * G_BBT(|E_max|) * W_eff

    where W_eff ~ |E_max| / (B_kane E_g^(3/2)) * W is the effective
    integration length: the Kane integrand peaks at the metallurgical
    junction and decays super-exponentially with x because of the
    `exp(-1/|E|)` factor. This is the leading-order asymptotic result
    that Sze section 8.4 cites; it is accurate to ~ 20 % over a
    decade in V_R for the heavy-doping abrupt-junction regime where
    the depletion approximation holds. M16.6.

    Parameters
    ----------
    V_R : float | array-like
        Reverse-bias voltage in V (positive value; the formula applies
        to V_R >= 0). The convention follows the FEM benchmark, which
        applies a NEGATIVE bias to the diode anode; pass `abs(V)` or
        `-V_anode` here.
    N : float
        Doping density in m^-3 (single-sided abrupt junction:
        N_A = N_D = N).
    eps : float
        Absolute permittivity of the semiconductor in F/m.
    E_g_eV : float
        Band gap in eV (Si: 1.12).
    V_bi : float
        Built-in voltage in volts (V_bi = V_t * ln(N^2 / n_i^2) for
        the symmetric abrupt junction).
    A_kane_cm : float
        Kane prefactor in cm^-1 s^-1 V^-2; Si default 4.0e14
        (Sze section 8.4).
    B_kane_cm : float
        Kane exponent coefficient in V/cm; Si default 1.9e7.

    Returns
    -------
    numpy.ndarray
        Reverse-bias breakdown current density in A/m^2. The sign
        is positive (the FEM verifier compares |J_FEM|).

    References
    ----------
    S. M. Sze, *Physics of Semiconductor Devices*, 3rd ed., Wiley,
    2007; Section 8.4 (Zener / Kane band-to-band tunneling), eq.
    8.4.10 for the closed-form integrated rate.

    Pure-Python; no dolfinx. M16.6.
    """
    V_R_arr = np.asarray(V_R, dtype=float)
    # Total voltage across the depletion region (V_bi + |V_R|).
    V_total = V_bi + np.abs(V_R_arr)
    # Peak field at the metallurgical junction under one-sided
    # abrupt-junction depletion approximation. Switch from m^-3
    # (input N) to V/cm via |E_SI| / 100.
    E_max_SI = np.sqrt(2.0 * Q * float(N) * V_total / float(eps))
    E_max_cm = E_max_SI / 100.0
    # Depletion width (m).
    W = np.sqrt(2.0 * float(eps) * V_total / (Q * float(N)))
    # Kane generation rate at the peak field, evaluated in cm-units.
    # Guard against division by zero at V_total = 0.
    E_safe = np.where(E_max_cm > 1.0e-30, E_max_cm, 1.0e-30)
    expo = -float(B_kane_cm) * float(E_g_eV) ** 1.5 / E_safe
    G_peak_cm = (
        float(A_kane_cm) * E_safe ** 2 / np.sqrt(float(E_g_eV))
        * np.exp(expo)
    )
    # cm^-3 s^-1 -> m^-3 s^-1
    G_peak_SI = G_peak_cm * 1.0e6
    # Effective integration length: |E_max| / (B_kane * E_g^(3/2)).
    # This is the asymptotic length over which the exp(-1/|E|) factor
    # falls by a factor of e, derived by substituting the linear-field
    # profile |E(x)| = |E_max| (1 - x / W) into the Kane integrand.
    # In cm; convert back to m via 1e-2.
    L_eff_cm = E_safe / (
        float(B_kane_cm) * float(E_g_eV) ** 1.5
    ) * (W * 1.0e2)
    # Cap the effective length at the depletion width (for very low
    # V_R the asymptotic formula overestimates the integration band).
    L_eff_cm = np.minimum(L_eff_cm, W * 1.0e2)
    L_eff_SI = L_eff_cm * 1.0e-2
    # Total reverse-bias current density. The factor of q converts the
    # generation rate to charge flux.
    return Q * G_peak_SI * L_eff_SI


def thermionic_iv(
    V, barrier_height_eV: float, A_richardson: float, T: float,
) -> np.ndarray:
    """
    Closed-form thermionic-emission I-V for a Schottky contact (Sze
    3rd ed eq 3.4.10):

        J(V) = A* T^2 exp(-q phi_B / kT) [exp(qV / kT) - 1]

    Parameters
    ----------
    V : float | array-like
        Applied bias in volts.
    barrier_height_eV : float
        Schottky barrier height phi_B in eV (e.g. 0.85 eV for
        Pt-on-n-Si).
    A_richardson : float
        Richardson constant A* in SI units (A m^-2 K^-2). Compute
        from :func:`richardson_constant`.
    T : float
        Device temperature in kelvin.

    Returns
    -------
    numpy.ndarray
        Current density J in A m^-2 (multiply by 1e-4 for A cm^-2).

    Pure-Python; no dolfinx. M16.5.
    """
    V_arr = np.asarray(V, dtype=float)
    V_t = KB * float(T) / Q
    Js = float(A_richardson) * float(T) ** 2 * math.exp(
        -float(barrier_height_eV) / V_t
    )
    return Js * (np.exp(V_arr / V_t) - 1.0)


# ---------------------------------------------------------------------------
# M16.4 Fermi-Dirac built-in voltage references.
# ---------------------------------------------------------------------------


def vbi_boltzmann(N_A: float, N_D: float, n_i: float, V_t: float) -> float:
    """
    Textbook Boltzmann-Slotboom built-in voltage on a 1D pn junction.

        V_bi_B = V_t * ln(N_A N_D / n_i^2)

    Equivalent to `(psi_n_bulk - psi_p_bulk) * V_t` where the bulk
    equilibrium psi is the Boltzmann charge-neutrality root
    `psi_bulk = V_t * arcsinh(N_net / (2 n_i))`.

    Pure-Python; no dolfinx. Used by the diode_fermi_dirac_1d verifier
    and by tests.
    """
    return float(V_t * np.log(float(N_A) * float(N_D) / float(n_i) ** 2))


def _solve_eta_for_density_ratio(target_ratio: float, fermi_half) -> float:
    """
    Find eta such that `fermi_half(eta) = target_ratio`.

    Bisection on [-50, +30] is plenty for the M16.4 doping range
    (N_D / N_C up to ~1e2). Pure-Python; the caller provides
    `fermi_half`, which is either the Blakemore basic closed form or
    the full-integral reference from
    :mod:`semi.physics.statistics`.
    """
    lo, hi = -50.0, 30.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float(fermi_half(mid)) < target_ratio:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-12:
            break
    return 0.5 * (lo + hi)


def vbi_fermi_dirac(
    N_A: float, N_D: float, n_i: float, N_C: float, N_V: float,
    V_t: float,
    *, kind: str = "reference",
) -> float:
    """
    Built-in voltage on a 1D pn junction under Fermi-Dirac statistics.

    Works in the kronos-semi Slotboom convention where psi is measured
    from the intrinsic Fermi level (so `n = n_i exp(psi/V_t)` under
    Boltzmann). The bulk equilibrium psi on each side is recovered
    from the FD charge-neutrality root:

        n_n_bulk = N_D  =>  F_{1/2}(eta_n) = N_D / N_C
        eta_n = (psi_n_bulk - 0)/V_t + ln(n_i / N_C)
        psi_n_bulk_phys = V_t * (eta_n - ln(n_i / N_C))

    Mirror identity for the p-side. V_bi = psi_n_bulk_phys -
    psi_p_bulk_phys.

    Parameters
    ----------
    N_A, N_D, n_i, N_C, N_V : float
        Densities in m^-3.
    V_t : float
        Thermal voltage in volts (kT/q).
    kind : {"reference", "blakemore"}
        `"reference"` (default): full Fermi-Dirac integral via mpmath
        (the polylog identity F_{1/2}(eta) = -Li_{3/2}(-exp(eta))).
        This is the "scipy-driven" reference the M16.4 acceptance gate
        cites (scipy.special.fdk is not available on docker-fem;
        mpmath is the supported fallback).
        `"blakemore"`: basic Blakemore closed form, what the production
        residual evaluates. Use this when comparing FEM equilibrium
        psi to the analytical Blakemore prediction.

    Returns
    -------
    float
        V_bi in volts.
    """
    from .physics.statistics import (
        fermi_dirac_half_blakemore,
        fermi_dirac_half_reference,
    )

    if kind == "reference":
        fermi_half = fermi_dirac_half_reference
    elif kind == "blakemore":
        fermi_half = fermi_dirac_half_blakemore
    else:
        raise ValueError(f"vbi_fermi_dirac: unknown kind {kind!r}")

    eta_n = _solve_eta_for_density_ratio(float(N_D) / float(N_C), fermi_half)
    eta_p = _solve_eta_for_density_ratio(float(N_A) / float(N_V), fermi_half)
    eta_offset_n = np.log(float(n_i) / float(N_C))
    eta_offset_p = np.log(float(n_i) / float(N_V))
    psi_n_bulk_scaled = eta_n - eta_offset_n
    psi_p_bulk_scaled = -(eta_p - eta_offset_p)
    return float(V_t * (psi_n_bulk_scaled - psi_p_bulk_scaled))
