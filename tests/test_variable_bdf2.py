"""
Tests for the M18 variable-step BDF2 coefficient classmethod
``BDFCoefficients.variable_bdf2(omega)``. Pure-Python; no dolfinx.

The variable-step BDF2 stencil uses
``omega = dt_n / dt_{n-1}`` as the only free parameter; at ``omega == 1.0``
the triplet must collapse to the uniform-step values stored in
``BDFCoefficients._SUPPORTED_COEFFS[2]``, which is the bit-identity hinge
for the adaptive-dt path in ``semi/runners/transient.py``.
"""
from __future__ import annotations

import pytest

from semi.timestepping import BDFCoefficients


def test_variable_bdf2_uniform_collapses_to_supported():
    """omega = 1.0 reduces to the uniform-step (3/2, -2, 1/2) triplet."""
    alpha = BDFCoefficients.variable_bdf2(1.0)
    assert alpha == (1.5, -2.0, 0.5)
    assert alpha == BDFCoefficients._SUPPORTED_COEFFS[2]


def test_variable_bdf2_omega_two():
    """omega = 2.0 returns (5/3, -3, 4/3)."""
    a0, a1, a2 = BDFCoefficients.variable_bdf2(2.0)
    assert a0 == pytest.approx(5.0 / 3.0)
    assert a1 == pytest.approx(-3.0)
    assert a2 == pytest.approx(4.0 / 3.0)


def test_variable_bdf2_omega_half():
    """omega = 0.5 returns (4/3, -1.5, 1/6)."""
    a0, a1, a2 = BDFCoefficients.variable_bdf2(0.5)
    assert a0 == pytest.approx(4.0 / 3.0)
    assert a1 == pytest.approx(-1.5)
    assert a2 == pytest.approx(1.0 / 6.0)


@pytest.mark.parametrize("omega", [0.1, 0.5, 0.99, 1.0, 1.01, 1.5, 2.0, 3.7, 10.0])
def test_variable_bdf2_constant_exactness(omega):
    """Constant-in-time exactness: alpha_0 + alpha_1 + alpha_2 == 0."""
    a0, a1, a2 = BDFCoefficients.variable_bdf2(omega)
    assert a0 + a1 + a2 == pytest.approx(0.0, abs=1.0e-13)


@pytest.mark.parametrize("omega", [0.1, 0.5, 0.99, 1.0, 1.01, 1.5, 2.0, 3.7, 10.0])
def test_variable_bdf2_linear_exactness(omega):
    """
    Linear-in-time exactness: alpha_0 - alpha_2 / omega == 1.

    This is the standard non-uniform-BDF2 derivation check. For
    u(t) = a + b*t the discrete derivative
    (alpha_0 * u_{n+1} + alpha_1 * u_n + alpha_2 * u_{n-1}) / dt_n
    must equal b at every step, which (after using the constant-in-time
    identity to eliminate the offset) reduces to this single condition.
    """
    a0, _a1, a2 = BDFCoefficients.variable_bdf2(omega)
    assert a0 - a2 / omega == pytest.approx(1.0, abs=1.0e-13)


def test_variable_bdf2_zero_omega_raises():
    with pytest.raises(ValueError, match="omega must be strictly positive"):
        BDFCoefficients.variable_bdf2(0.0)


def test_variable_bdf2_negative_omega_raises():
    with pytest.raises(ValueError, match="omega must be strictly positive"):
        BDFCoefficients.variable_bdf2(-1.0)


def test_uniform_bdf2_unchanged_after_m18():
    """The _SUPPORTED_COEFFS table is the authoritative uniform-step lookup
    and must remain (1.5, -2.0, 0.5) for BDF2 after M18 (bit-identity gate)."""
    assert BDFCoefficients._SUPPORTED_COEFFS[2] == (1.5, -2.0, 0.5)
    assert BDFCoefficients._SUPPORTED_COEFFS[1] == (1.0, -1.0)
