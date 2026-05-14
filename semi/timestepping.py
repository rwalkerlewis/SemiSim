"""
BDF (backward differentiation formula) time-stepping coefficients.

Pure-Python, no dolfinx at module scope. Provides:

    BDFCoefficients(order)            -- uniform-step BDF1 / BDF2
    BDFCoefficients.variable_bdf2(omega) -- variable-step BDF2 (M18)

Design:
    order=1  (backward Euler / BDF1): coeffs = (1, -1)
    order=2  (BDF2):                  coeffs = (3/2, -2, 1/2)

The general BDF-k formula for the discrete time derivative at step n+1:

    du/dt ~ (1/dt) * sum_{k=0}^{order} alpha_k * u^{n+1-k}

where u^{n+1} is the current (unknown) value and u^{n}, u^{n-1}, ...
are the previous (known) values stored in `u_history`.

The :meth:`apply` method evaluates this sum given a list of PAST values
(not including the current unknown). The caller is responsible for adding
the alpha_0/dt contribution of the current unknown to the UFL residual.

For variable-step BDF2 (M18), the coefficients depend on the dt ratio
omega = dt_n / dt_{n-1}. The classmethod :meth:`variable_bdf2` returns
the (alpha_0, alpha_1, alpha_2) triplet for a given omega; at omega = 1.0
it reduces to the uniform-step values (3/2, -2, 1/2) bit-identically.
The adaptive-dt path in `semi/runners/transient.py` consumes these
coefficients per step.
"""
from __future__ import annotations

import numpy as np


class BDFCoefficients:
    """
    Coefficients for backward differentiation formula time integration.

    Parameters
    ----------
    order : int
        BDF order. 1 = backward Euler, 2 = BDF2. Only 1 and 2 are
        supported for M13.

    Attributes
    ----------
    order : int
    coeffs : tuple of float
        (alpha_0, alpha_1, ..., alpha_order). alpha_0 multiplies u^{n+1}
        (the current unknown); alpha_k for k >= 1 multiply the k-th
        previous solution u^{n+1-k}.

    Methods
    -------
    apply(u_history, dt) -> numpy.ndarray
        Compute the full discrete time-derivative approximation
        sum_k(alpha_k * u_history[-k-1]) / dt, where u_history[-1] is
        the most recent value and u_history[-2] is one step back.
    """

    _SUPPORTED_COEFFS = {
        1: (1.0, -1.0),
        2: (1.5, -2.0, 0.5),
    }

    def __init__(self, order: int) -> None:
        if order not in self._SUPPORTED_COEFFS:
            raise ValueError(
                f"BDF order {order!r} is not supported; use 1 or 2."
            )
        self.order: int = order
        self.coeffs: tuple[float, ...] = self._SUPPORTED_COEFFS[order]

    @classmethod
    def variable_bdf2(cls, omega: float) -> tuple[float, float, float]:
        """
        Variable-step BDF2 coefficients at dt ratio ``omega = dt_n / dt_{n-1}``.

        Parameters
        ----------
        omega : float
            Ratio of the current step to the previous step. Must be
            strictly positive.

        Returns
        -------
        (alpha_0, alpha_1, alpha_2) : tuple of float
            The three coefficients of the BDF2 stencil at non-uniform
            spacing:

                alpha_0 = (1 + 2*omega) / (1 + omega)
                alpha_1 = -(1 + omega)
                alpha_2 = omega**2 / (1 + omega)

            The discrete time derivative at t_{n+1} is
            ``(alpha_0 * u_{n+1} + alpha_1 * u_n + alpha_2 * u_{n-1}) / dt_n``.

        Notes
        -----
        - At ``omega == 1.0`` the triplet collapses to ``(1.5, -2.0, 0.5)``,
          bit-identical to the uniform-step BDF2 values stored in
          ``_SUPPORTED_COEFFS[2]``. This is the bit-identity hinge for the
          adaptive-dt path: when the controller does not change dt across
          two consecutive steps, the variable-step branch produces the same
          coefficients the uniform-step branch always produced.
        - ``alpha_0 + alpha_1 + alpha_2 == 0`` (constant-in-time exactness).
        - ``alpha_0 - alpha_2 / omega == 1`` (linear-in-time exactness;
          this is the standard non-uniform-BDF2 derivation check).
        """
        if omega <= 0.0:
            raise ValueError(f"omega must be strictly positive; got {omega!r}")
        denom = 1.0 + omega
        alpha_0 = (1.0 + 2.0 * omega) / denom
        alpha_1 = -(1.0 + omega)
        alpha_2 = (omega * omega) / denom
        return (alpha_0, alpha_1, alpha_2)

    def apply(self, u_history: list, dt: float) -> np.ndarray:
        """
        Evaluate the BDF discrete time derivative approximation.

        Parameters
        ----------
        u_history : list of numpy.ndarray
            Time history of the unknown, ordered from oldest to newest.
            Must contain at least `order + 1` entries.
            u_history[-1] is u^{n+1} (the newest value).
            u_history[-2] is u^n, etc.
        dt : float
            Time step size.

        Returns
        -------
        numpy.ndarray
            The discrete time derivative sum_k(alpha_k * u_history[-k-1]) / dt.
            For BDF1: (u^{n+1} - u^n) / dt
            For BDF2: (3/2 u^{n+1} - 2 u^n + 1/2 u^{n-1}) / dt
        """
        if len(u_history) < self.order + 1:
            raise ValueError(
                f"BDF{self.order} requires at least {self.order + 1} "
                f"history values; got {len(u_history)}."
            )
        result = np.zeros_like(np.asarray(u_history[-1], dtype=float))
        for k, alpha_k in enumerate(self.coeffs):
            result += alpha_k * np.asarray(u_history[-k - 1], dtype=float)
        return result / dt
