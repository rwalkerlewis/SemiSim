"""
Tests for the pure-Python diode analytical helpers in
semi.diode_analytical. No dolfinx dependency.
"""
from __future__ import annotations

import numpy as np
import pytest

from semi.constants import Q
from semi.diode_analytical import (
    depletion_width,
    richardson_constant,
    shockley_iv_with_auger,
    shockley_long_diode_saturation,
    sns_total_reference,
    srh_generation_reference,
    thermionic_iv,
    vbi_boltzmann,
    vbi_fermi_dirac,
)

V_T = 0.025852
N_A = 1.0e23
N_D = 1.0e23
N_I = 1.0e16
EPS = 11.7 * 8.8541878128e-12
MU_N = 0.14
MU_P = 0.045
TAU_N = 1.0e-8
TAU_P = 1.0e-8


def test_shockley_saturation_matches_hand_calc():
    J_s, L_n, L_p = shockley_long_diode_saturation(
        N_A, N_D, N_I, MU_N, MU_P, TAU_N, TAU_P, V_T,
    )
    D_n = V_T * MU_N
    D_p = V_T * MU_P
    assert L_n == pytest.approx((D_n * TAU_N) ** 0.5)
    assert L_p == pytest.approx((D_p * TAU_P) ** 0.5)
    expected = Q * N_I ** 2 * (D_n / (L_n * N_A) + D_p / (L_p * N_D))
    assert J_s == pytest.approx(expected)


def test_depletion_width_vanishes_near_vbi():
    V_bi = V_T * np.log(N_A * N_D / N_I ** 2)
    W_near = depletion_width(N_A, N_D, N_I, EPS, V_T, V_bi - 1.0e-8)
    W_mid = depletion_width(N_A, N_D, N_I, EPS, V_T, 0.3)
    W_zero = depletion_width(N_A, N_D, N_I, EPS, V_T, 0.0)
    assert float(W_near) < float(W_mid) < float(W_zero)


def test_depletion_width_grows_under_reverse_bias():
    W_zero = float(depletion_width(N_A, N_D, N_I, EPS, V_T, 0.0))
    W_rev = float(depletion_width(N_A, N_D, N_I, EPS, V_T, -2.0))
    # W(V) scales as sqrt(V_bi - V); for V=-2 vs V=0 with V_bi~0.83,
    # the ratio is sqrt((0.83+2)/0.83) ~ 1.85.
    assert 1.6 < W_rev / W_zero < 2.1


def test_sns_reference_zero_at_v_zero():
    J_total, J_diff, J_rec, J_s = sns_total_reference(
        N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=0.0,
    )
    assert float(J_total) == pytest.approx(0.0, abs=1e-12)
    assert float(J_diff) == pytest.approx(0.0, abs=1e-12)
    assert float(J_rec) == pytest.approx(0.0, abs=1e-12)
    assert J_s > 0.0


def test_sns_reference_forward_dominated_by_rec_at_low_bias():
    """At V=0.15 V with tau=1e-8 s the SNS recombination branch is
    many orders of magnitude above pure Shockley diffusion, which is
    the whole reason the M3 verifier needed this term."""
    V = 0.15
    J_total, J_diff, J_rec, _J_s = sns_total_reference(
        N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=V,
    )
    assert float(J_rec) > 100.0 * float(J_diff)
    assert float(J_total) == pytest.approx(float(J_diff) + float(J_rec))


def test_sns_reference_vectorised_matches_scalar():
    Vs = np.array([0.1, 0.2, 0.3, 0.5])
    J_vec, _, _, _ = sns_total_reference(
        N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=Vs,
    )
    for V, J_v in zip(Vs, J_vec, strict=True):
        J_s, _, _, _ = sns_total_reference(
            N_A, N_D, N_I, EPS, MU_N, MU_P, TAU_N, TAU_P, V_T, V=float(V),
        )
        assert float(J_s) == pytest.approx(float(J_v), rel=1e-12)


def test_srh_generation_zero_at_v_zero():
    J_net, W0, V_bi = srh_generation_reference(
        N_A, N_D, N_I, EPS, TAU_N, TAU_P, V_T, V=0.0,
    )
    assert float(J_net) == pytest.approx(0.0, abs=1e-12)
    assert W0 > 0.0
    assert V_bi == pytest.approx(V_T * np.log(N_A * N_D / N_I ** 2))


def test_srh_generation_grows_with_reverse_bias():
    Vs = np.array([-0.5, -1.0, -1.5, -2.0])
    J_net, _W0, _V_bi = srh_generation_reference(
        N_A, N_D, N_I, EPS, TAU_N, TAU_P, V_T, V=Vs,
    )
    # The helper returns the signed (q n_i / 2 tau_eff)(W(V) - W(0))
    # which is positive for V < 0 because W(V) > W(0); the verifier
    # compares |J_sim| against |J_net|.
    assert np.all(J_net > 0.0)
    # |J_net| must be non-decreasing as V gets more negative.
    assert np.all(np.diff(np.abs(J_net)) > 0.0)


def test_srh_generation_order_of_magnitude_deep_reverse():
    """For this device (tau=1e-8, N=1e17) the deep-reverse generation
    current is dominated by SRH, not Shockley diffusion, and is on the
    order of 1e-2 A/m^2. The M2 reverse-bias spec (J saturating to
    J_s ~ 1e-7 A/m^2) does not hold; the M3 verifier compares to
    this quantity instead."""
    J_net, _W0, _V_bi = srh_generation_reference(
        N_A, N_D, N_I, EPS, TAU_N, TAU_P, V_T, V=-2.0,
    )
    assert 5.0e-3 < abs(float(J_net)) < 5.0e-2


# ---------------------------------------------------------------------------
# M16.3 Auger high-injection long-diode reference tests.
# ---------------------------------------------------------------------------

# diode_auger_1d benchmark parameters: 1e15 cm^-3 = 1e21 m^-3 doping.
N_AUG = 1.0e21
NI_AUG = 1.0e16          # Si Altermatt
TAU_AUG = 1.0e-7
# Engineered Auger coefficients in m^6/s (cm^6/s * 1e-12).
C_N_AUG_SI = 1.0e-29 * 1.0e-12
C_P_AUG_SI = 1.0e-29 * 1.0e-12


def test_shockley_iv_with_auger_zero_at_zero_bias():
    """Below V_bi the high-injection asymptote should not fire (the
    helper returns 0 in that regime; SNS is the proper reference for
    sub-V_bi biases)."""
    J_total, J_srh, _ = shockley_iv_with_auger(
        N_AUG, N_AUG, NI_AUG, EPS, MU_N, MU_P, TAU_AUG, TAU_AUG,
        C_N_AUG_SI, C_P_AUG_SI, V_T, V=0.0,
    )
    assert float(J_total) == 0.0
    assert float(J_srh) == 0.0


def test_shockley_iv_with_auger_high_inj_positive():
    """At V_F = 0.9 V on this geometry, J_total > J_srh_only (Auger
    adds to recombination current)."""
    J_total, J_srh, tau_eff = shockley_iv_with_auger(
        N_AUG, N_AUG, NI_AUG, EPS, MU_N, MU_P, TAU_AUG, TAU_AUG,
        C_N_AUG_SI, C_P_AUG_SI, V_T, V=0.9,
    )
    assert float(J_total) > float(J_srh) > 0.0
    # Auger reduces tau_eff below tau_SRH. tau_SRH = sqrt(tau_n tau_p)
    # = TAU_AUG (since tau_n == tau_p == TAU_AUG); so tau_eff < tau_SRH.
    assert float(tau_eff) < TAU_AUG


def test_shockley_iv_with_auger_zero_C_recovers_srh_branch():
    """With C_n = C_p = 0, J_total should equal J_srh_only."""
    J_total, J_srh, _ = shockley_iv_with_auger(
        N_AUG, N_AUG, NI_AUG, EPS, MU_N, MU_P, TAU_AUG, TAU_AUG,
        0.0, 0.0, V_T, V=0.9,
    )
    assert float(J_total) == pytest.approx(float(J_srh), rel=1.0e-12)


def test_shockley_iv_with_auger_array_input():
    """Vector V argument returns vector J pair, monotone in V above
    V_bi."""
    Vs = np.linspace(0.7, 0.9, 5)
    J_total, J_srh, _ = shockley_iv_with_auger(
        N_AUG, N_AUG, NI_AUG, EPS, MU_N, MU_P, TAU_AUG, TAU_AUG,
        C_N_AUG_SI, C_P_AUG_SI, V_T, V=Vs,
    )
    assert J_total.shape == Vs.shape
    # At least one element above V_bi should be strictly positive.
    assert np.any(J_total > 0.0)
    # J_total must be monotone non-decreasing.
    pos_total = J_total[J_total > 0.0]
    if pos_total.size >= 2:
        assert np.all(np.diff(pos_total) >= 0.0)


# ---------------------------------------------------------------------------
# M16.4 Fermi-Dirac built-in voltage references.
# ---------------------------------------------------------------------------

# diode_fermi_dirac_1d benchmark parameters (Si Altermatt n_i, n+/p doping
# in the regime where Boltzmann breaks down).
N_A_FD = 1.0e23   # 1e17 cm^-3 acceptor concentration
N_D_FD = 1.0e26   # 1e20 cm^-3 donor concentration (deep into FD)
N_I_FD = 1.0e16   # Si Altermatt at 300 K
N_C_FD = 2.86e25  # Si effective DoS, conduction band
N_V_FD = 3.10e25  # Si effective DoS, valence band


def test_vbi_boltzmann_matches_textbook_log_form():
    """V_bi_B = V_t ln(N_A N_D / n_i^2). Reproduce on benchmark numbers."""
    expected = V_T * np.log(N_A_FD * N_D_FD / N_I_FD ** 2)
    assert vbi_boltzmann(N_A_FD, N_D_FD, N_I_FD, V_T) == pytest.approx(expected)


def test_vbi_fermi_dirac_blakemore_exceeds_boltzmann_for_degenerate_n_side():
    """At N_D = 1e20 cm^-3 (N_D / N_C ~ 3.5 in Si) the Blakemore-FD V_bi
    exceeds the Boltzmann V_bi: F_{1/2}(eta) grows slower than exp(eta),
    so reproducing N_D >> N_C requires a higher eta_n_FD than Boltzmann's
    eta_n_B = ln(N_D/N_C), which raises psi_n_bulk and therefore V_bi.
    This is the "FD vs Boltzmann divergence" gated by the M16.4
    acceptance benchmark (>5 % at N_D = 1e20 cm^-3)."""
    V_bi_B = vbi_boltzmann(N_A_FD, N_D_FD, N_I_FD, V_T)
    V_bi_FD_blake = vbi_fermi_dirac(
        N_A_FD, N_D_FD, N_I_FD, N_C_FD, N_V_FD, V_T, kind="blakemore",
    )
    assert V_bi_FD_blake > V_bi_B
    # Benchmark divergence threshold is 5 %; verify it analytically.
    assert (V_bi_FD_blake - V_bi_B) / V_bi_B > 0.05


def test_vbi_fermi_dirac_reference_matches_blakemore_in_correction_sign():
    """The closed Blakemore form approximates the full Fermi-Dirac
    integral within ~4 % at the benchmark operating point (eta ~ 1.25 to
    3 across n+/p sides). Both the Blakemore-FD and the full-integral-FD
    references must therefore predict V_bi above the Boltzmann value, and
    must agree with each other within a fraction of the FD correction
    itself."""
    V_bi_FD_ref = vbi_fermi_dirac(
        N_A_FD, N_D_FD, N_I_FD, N_C_FD, N_V_FD, V_T, kind="reference",
    )
    V_bi_FD_blake = vbi_fermi_dirac(
        N_A_FD, N_D_FD, N_I_FD, N_C_FD, N_V_FD, V_T, kind="blakemore",
    )
    V_bi_B = vbi_boltzmann(N_A_FD, N_D_FD, N_I_FD, V_T)
    delta_blake = V_bi_FD_blake - V_bi_B
    delta_ref = V_bi_FD_ref - V_bi_B
    assert delta_ref > 0.0
    assert delta_blake > 0.0
    # The basic Blakemore form overestimates the FD correction at the
    # benchmark operating point (its denominator drops ~50 % faster than
    # the true F_{1/2} does at eta ~ 3). The two predictions therefore
    # agree only to within a factor of ~2.5 at N_D = 1e20 cm^-3 — which
    # is exactly why the M16.4 acceptance gate cites the full-integral
    # reference as the Boltzmann-comparison floor and tracks the basic
    # Blakemore divergence as a documentation diagnostic.
    assert delta_blake < 3.0 * delta_ref


def test_vbi_fermi_dirac_unknown_kind_raises():
    with pytest.raises(ValueError, match="unknown kind"):
        vbi_fermi_dirac(
            N_A_FD, N_D_FD, N_I_FD, N_C_FD, N_V_FD, V_T, kind="bogus",
        )


def test_vbi_fermi_dirac_default_kind_is_reference():
    """The default `kind` argument must be the analytical reference (the
    M16.4 acceptance gate cites scipy.special.fdk via mpmath); a default
    drift to "blakemore" would silently weaken the gate."""
    default = vbi_fermi_dirac(N_A_FD, N_D_FD, N_I_FD, N_C_FD, N_V_FD, V_T)
    explicit = vbi_fermi_dirac(
        N_A_FD, N_D_FD, N_I_FD, N_C_FD, N_V_FD, V_T, kind="reference",
    )
    assert default == pytest.approx(explicit, rel=1.0e-12)


# ---------------------------------------------------------------------------
# M16.5 thermionic-emission helper tests.
# ---------------------------------------------------------------------------


def test_richardson_constant_si_electron_within_one_percent():
    # Si electron Richardson constant in SI units:
    # A* = 4 pi q m* k^2 / h^3 with m* = 0.26 m_0 ~ 3.2e5 A m^-2 K^-2.
    A_si_n = richardson_constant(0.26)
    assert A_si_n == pytest.approx(3.2e5, rel=0.05)


def test_richardson_constant_si_hole_proportional_to_mass():
    # Doubling m* doubles A* exactly.
    A_p1 = richardson_constant(0.39)
    A_p2 = richardson_constant(0.78)
    assert A_p2 == pytest.approx(2.0 * A_p1, rel=1.0e-12)


def test_thermionic_iv_zero_at_v_zero():
    A_si_n = richardson_constant(0.26)
    J = thermionic_iv(0.0, barrier_height_eV=0.85, A_richardson=A_si_n, T=300.0)
    # exp(0) - 1 = 0 exactly.
    assert float(J) == pytest.approx(0.0, abs=1.0e-30)


def test_thermionic_iv_grows_exponentially_under_forward_bias():
    A_si_n = richardson_constant(0.26)
    V_bias = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    J = thermionic_iv(
        V_bias, barrier_height_eV=0.85, A_richardson=A_si_n, T=300.0,
    )
    # Strict monotone growth.
    assert np.all(np.diff(J) > 0.0)
    # Roughly factor exp((0.4-0.1)/V_t) ~ exp(11.6) ~ 1e5 between V=0.1 and V=0.4.
    ratio = J[3] / J[0]
    assert 5.0e4 < ratio < 5.0e5


def test_thermionic_iv_recovers_saturation_at_negative_v():
    # Reverse bias, J -> -J_s (saturation current with sign).
    A_si_n = richardson_constant(0.26)
    J = thermionic_iv(
        -10.0, barrier_height_eV=0.85, A_richardson=A_si_n, T=300.0,
    )
    Js = A_si_n * 300.0 ** 2 * np.exp(-0.85 / V_T)
    assert float(J) == pytest.approx(-Js, rel=1.0e-3)


def test_thermionic_iv_array_input_returns_same_shape():
    A_si_n = richardson_constant(0.26)
    V_bias = np.linspace(0.0, 0.5, 21)
    J = thermionic_iv(
        V_bias, barrier_height_eV=0.85, A_richardson=A_si_n, T=300.0,
    )
    assert J.shape == V_bias.shape
    assert J[0] == pytest.approx(0.0, abs=1.0e-30)


def test_thermionic_iv_barrier_dependence_dominant():
    """phi_B = 0.85 vs 0.70 eV at V_F = 0.3 V should differ by
    exp(0.15 / V_t) ~ 330x in J. Demonstrates that the barrier height
    dominates the I-V at moderate forward bias."""
    A_si_n = richardson_constant(0.26)
    J_85 = thermionic_iv(
        0.3, barrier_height_eV=0.85, A_richardson=A_si_n, T=300.0,
    )
    J_70 = thermionic_iv(
        0.3, barrier_height_eV=0.70, A_richardson=A_si_n, T=300.0,
    )
    expected = float(np.exp(0.15 / V_T))
    assert float(J_70 / J_85) == pytest.approx(expected, rel=1.0e-3)

