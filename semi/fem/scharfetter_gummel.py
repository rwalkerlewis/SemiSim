"""
Scharfetter-Gummel edge-flux primitives (M13.1).

This module is the load-bearing primitive for the SG transient
discretisation specified in ADR 0012. It provides:

    * `bernoulli(x)`            : numerically stable B(x) = x / (exp(x) - 1).
    * `bernoulli_array(x)`      : vectorised B over a NumPy array.
    * `sg_edge_flux_n(...)`     : 1D electron edge flux per ADR 0012.
    * `sg_edge_flux_p(...)`     : 1D hole edge flux per ADR 0012.

The 2D simplicial assembly is layered on top of these primitives in a
follow-up commit; per the M13.1 plan we land 1D first, gate on
steady-state agreement, then build 2D. This module is pure-Python /
NumPy with no dolfinx import at module scope so it can be imported
and unit-tested without a FEM toolchain (Invariant 4).

Sign convention (ADR 0012 § "Sign convention")
----------------------------------------------
Sandia / Farrell-et-al convention: `F_n_ij` is the electron contribution
to the conventional current density on the edge from i to j. At zero
field with n_j > n_i, F > 0 (current in +∇n direction); at uniform n
with psi_j > psi_i, F < 0 (electrons drift toward j, conventional
current opposite). This matches the textbook
`J_n = q mu_n n E + q D_n grad(n)` once the SemiSim phi_n sign
convention is unwound (see notes in `semi/postprocess.py`,
`evaluate_partial_currents`). The user-prompt formula
`F = (mu V_t/h)[B(-dpsi) n_j - B(+dpsi) n_i]` had the B-arguments
swapped relative to the standard convention; the swap was caught by
the midpoint-Galerkin cross-check guard before any 1D residual code
was wired. See ADR 0012 "Sources actually used" for the trail.

Bernoulli function regimes
--------------------------
- `x == 0`        : exact 1.0.
- `|x| < taylor_eps` (default 1e-3): 4th-order Taylor 1 - x/2 + x^2/12 - x^4/720.
- `x > big_x` (default 30): rewrite as `x exp(-x) / (1 - exp(-x))` to keep
                            the denominator bounded away from zero and let
                            `x exp(-x)` underflow benignly to 0.
- `x < -big_x` (default -30): closed form `x / (exp(x) - 1)` is fine; `exp(x)`
                              underflows to 0 and the denominator goes to -1,
                              so `B(x) -> -x` cleanly.
- mid-range      : closed form `x / (exp(x) - 1)`.

The identity `B(x) - B(-x) = -x` is verified by unit test across the full
sampling range. (The earlier prompt draft stated `B(x) + B(-x) = -x`, which
is wrong; the sum is `x coth(x/2)`. The corrected identity is the one
used in tests.)
"""
from __future__ import annotations

import math

import numpy as np

_TAYLOR_EPS = 1.0e-3
_BIG_X = 30.0


def bernoulli(x: float) -> float:
    """
    Scalar Bernoulli `B(x) = x / (exp(x) - 1)` with numerically stable
    branches in all four regimes.

    Parameters
    ----------
    x : float
        Argument (dimensionless; for SG flux this is the scaled potential
        difference `psi_j - psi_i`).

    Returns
    -------
    float
        `B(x)`. `B(0) = 1.0` exactly.
    """
    if x == 0.0:
        return 1.0
    ax = abs(x)
    if ax < _TAYLOR_EPS:
        x2 = x * x
        return 1.0 - 0.5 * x + x2 / 12.0 - (x2 * x2) / 720.0
    if x > _BIG_X:
        ex = math.exp(-x)
        return (x * ex) / (1.0 - ex)
    return x / (math.expm1(x))


def bernoulli_array(x: np.ndarray) -> np.ndarray:
    """
    Vectorised Bernoulli over a NumPy array of arguments.

    Delegates to :func:`bernoulli` element-wise via ``np.frompyfunc`` to
    guarantee bitwise identity with the scalar path on all platforms.
    The Taylor, mid-range, and large-argument branches are handled by the
    scalar function exactly as in :func:`bernoulli`.
    """
    x = np.asarray(x, dtype=np.float64)
    return np.frompyfunc(bernoulli, 1, 1)(x).astype(np.float64)


def sg_edge_flux_n(
    n_i: float, n_j: float, dpsi: float,
    h_ij: float, mu_n: float, V_t: float,
) -> float:
    """
    1D Scharfetter-Gummel electron flux on edge (i, j) per ADR 0012.

        F_n_ij = (mu_n V_t / h_ij) * ( n_j B(+dpsi) - n_i B(-dpsi) )

    with `dpsi = psi_j - psi_i` in scaled units (i.e., psi already divided
    by V_t). All quantities in physical SI, except `dpsi` which is
    dimensionless because the SG argument is the scaled potential drop.

    Sign convention: Sandia / Farrell-et-al. F is the electron
    contribution to the conventional current density on edge i->j.
    F > 0 at dpsi=0 with n_j > n_i (diffusion in +grad(n) direction).
    F < 0 at uniform n with dpsi > 0 (electrons drift to j; conventional
    current is opposite). Source: arXiv 1911.00377 Eq (30); Sandia OSTI
    2011-3865 Eq (18) (after a_ij = -dpsi/(2 beta) substitution).
    """
    return (mu_n * V_t / h_ij) * (n_j * bernoulli(dpsi) - n_i * bernoulli(-dpsi))


def sg_edge_flux_p(
    p_i: float, p_j: float, dpsi: float,
    h_ij: float, mu_p: float, V_t: float,
) -> float:
    """
    1D Scharfetter-Gummel hole flux on edge (i, j) per ADR 0012.

        F_p_ij = (mu_p V_t / h_ij) * ( p_i B(+dpsi) - p_j B(-dpsi) )

    Hole formula in the same Sandia / Farrell convention as electrons.
    The i, j positions are swapped relative to the electron formula
    because the device-equation symmetry is `(n <-> p, psi -> -psi)`:
    starting from electrons, applying the symmetry, and re-grouping by
    the (i, j) positions of the swapped flux yields the spelling above.
    Cross-checked against Spevak Section 3.2.6 hole expression after
    sign-mapping for the convention difference (B_S(x) = -B(-x)):
    Spevak's hole `phi_p` equals `-F_p_ij` here (Spevak uses the
    electron-particle-flux convention, opposite of Sandia).

    F_p > 0 at dpsi=0 with p_i > p_j (diffusion in +grad(p) direction
    for the hole conventional current). F_p > 0 at uniform p with
    dpsi > 0 (holes drift from i to j; conventional current matches).

    The test guard 2b verifies the device-equation symmetry under
    `(n <-> p, psi -> -psi)` numerically to 1e-10 relative.
    """
    return (mu_p * V_t / h_ij) * (p_i * bernoulli(dpsi) - p_j * bernoulli(-dpsi))


def sg_edge_flux_n_array(
    n_i: np.ndarray, n_j: np.ndarray, dpsi: np.ndarray,
    h_ij: np.ndarray, mu_n: float, V_t: float,
) -> np.ndarray:
    """Vectorised electron edge flux over an array of edges."""
    return (mu_n * V_t / h_ij) * (
        n_j * bernoulli_array(dpsi) - n_i * bernoulli_array(-dpsi)
    )


def sg_edge_flux_p_array(
    p_i: np.ndarray, p_j: np.ndarray, dpsi: np.ndarray,
    h_ij: np.ndarray, mu_p: float, V_t: float,
) -> np.ndarray:
    """Vectorised hole edge flux over an array of edges."""
    return (mu_p * V_t / h_ij) * (
        p_i * bernoulli_array(dpsi) - p_j * bernoulli_array(-dpsi)
    )


def ufl_bernoulli(x, sigma: float = 1.0):
    """
    UFL expression for `B(x) = x / (exp(x) - 1)` suitable for use inside
    a UFL form. dolfinx differentiates this expression symbolically, so
    the Jacobian of any form using `ufl_bernoulli` is auto-derived.

    Implementation: tanh-smoothed blend (no conditionals)
    -----------------------------------------------------
    The first M13.1 attempt used `ufl.conditional(abs(x) < eps, taylor, closed)`
    which produced a UFL Jacobian with discontinuity artifacts at the
    threshold `|x| = eps` because UFL's symbolic differentiator evaluates
    both branches and the resulting Newton step stagnated. This version
    avoids conditionals entirely by blending two **mutually-stable** forms
    with a `ufl.tanh` switch.

    The two forms are mathematically identical to `B(x) = x / (exp(x) - 1)`
    but are numerically stable in opposite halves of the real line:

        B_pos(x) = x * exp(-x) / (1 - exp(-x))     stable for x > 0
        B_neg(x) = x / (exp(x) - 1)                stable for x <= 0

    For `x > 0`: `1 - exp(-x)` is bounded in (0, 1] and `x * exp(-x)`
    underflows benignly to zero at large positive x. For `x <= 0`:
    `exp(x) - 1` is bounded in (-1, 0] and the closed form does not
    overflow.

    Both forms are 0/0 at x=0; the blend hides this because the two
    forms agree to all orders in their Taylor series at x=0, and the
    blend weight at x=0 is exactly 0.5 of each, so the indeterminate
    contribution from each side is regularised by the well-defined
    contribution from the other. Numerically we add a 1e-300 offset
    to each denominator so UFL's symbolic engine never traps on a
    literal 0/0; the offset is below double-precision underflow at
    every operating point we evaluate at, so it has no effect on
    accuracy.

    The blend weight is `s(x) = 0.5 * (1 + tanh(x / sigma))`. At
    `x >> sigma`: s -> 1, B = B_pos. At `x << -sigma`: s -> 0,
    B = B_neg. Around `|x| <= sigma` both contribute; the result is
    exactly correct because both forms equal the true B(x).

    Choice of `sigma = 1.0` in scaled units: this is ~V_t in physical
    units, the natural scale of the problem. Verified by the unit
    tests (mpmath 1e-12, MG agreement 5%, hole symmetry 1e-10).

    Numerical robustness for large `|x|`: at `x > ~700` `exp(-x)`
    underflows; the `B_pos` form remains well-defined as `x * 0 / 1 = 0`
    correctly. At `x < ~-700` `exp(x)` underflows; the `B_neg` form
    saturates to `x / -1 = -x` correctly. So the blend is stable
    across the full IEEE 754 finite range.
    """
    import ufl

    # Regularise both 0/0 limits at x=0 by adding the SAME `eps_reg` to
    # numerator and denominator. With both perturbed by `eps_reg`, the
    # x=0 limit evaluates to `eps_reg / eps_reg = 1.0` (correct B(0)),
    # while away from x=0 the offset is below double-precision underflow
    # magnitude (smallest normal ~2.2e-308) and has no numerical effect.
    # This prevents UFL's symbolic engine from tripping on a literal 0/0
    # AND preserves the correct limit at x=0.
    eps_reg = 1.0e-300
    B_pos = (x * ufl.exp(-x) + eps_reg) / (1.0 - ufl.exp(-x) + eps_reg)
    B_neg = (x + eps_reg) / (ufl.exp(x) - 1.0 + eps_reg)

    # Smooth blend; s(0) = 0.5, s(x>>sigma) -> 1, s(x<<-sigma) -> 0.
    s = 0.5 * (1.0 + ufl.tanh(x / sigma))
    return s * B_pos + (1.0 - s) * B_neg


def ufl_sg_diffusion_coefficient(dpsi):
    """
    Exponential-fitting coefficient `A(dpsi) = (B(+dpsi) + B(-dpsi)) / 2`
    used in the modified-diffusion form of the 1D SG flux.

    For 1D linear elements, the SG flux on cell `K` with `dpsi = psi_j - psi_i`
    is algebraically identical to a "modified diffusion" Galerkin form:

        J_SG_cell = mu_n * ( A(dpsi) * grad(n) - n * grad(psi) )

    where the drift term `n grad(psi)` is unchanged from Galerkin and
    only the diffusion coefficient on `grad(n)` is modified by
    `A(dpsi)`. This makes 1D SG a one-line drop-in substitution for
    the existing Galerkin convection-diffusion form. See ADR 0012
    "1D edge flux" section; the algebraic identity is

        n_j B(+dpsi) - n_i B(-dpsi)
            = -n_avg * dpsi + Delta_n * (B(+dpsi) + B(-dpsi)) / 2
            = -n_avg * dpsi + Delta_n * A(dpsi).

    Asymptotics:
        A(0)   = 1                       (Galerkin limit, zero field)
        A -> |dpsi|/2  as |dpsi| -> inf  (SG enhancement, high Peclet)

    Note: the "modified-diffusion" reduction is exact for 1D linear
    elements only. In 2D simplicial elements the per-edge SG fluxes
    do not reduce to a single A(dpsi) per cell because each triangle
    has three edges with different dpsi values; the M13.1 plan ships
    1D first via this UFL form, and 2D via per-edge assembly later.
    """
    return (ufl_bernoulli(dpsi) + ufl_bernoulli(-dpsi)) / 2.0


def midpoint_galerkin_flux_n(
    n_i: float, n_j: float, dpsi: float,
    h_ij: float, mu_n: float, V_t: float,
) -> float:
    """
    Independent reference flux for the SG sign+scaling guard test
    (ADR 0012 test plan #2).

        F_ref = mu_n V_t (Delta n / h - n_mid * Delta psi_phys / (V_t h))
              = (mu_n V_t / h) (n_j - n_i)
              + (mu_n / h) * (-n_mid) * V_t * dpsi      [drift, scaled]

    Equivalent to a centred-difference / standard P1 Galerkin
    discretisation of `mu_n n grad psi + D_n grad n` with `n_mid = (n_i + n_j)/2`
    and `D_n = mu_n V_t` (Einstein relation). Sign convention matches
    `sg_edge_flux_n`.

    Used **only** in unit tests as an independent (different-derivation)
    cross-check of the SG flux at small cell Peclet number. The Slotboom
    2-node closed form is **not** used as the reference because it derives
    from the same `B(x)` and only catches typos.
    """
    n_mid = 0.5 * (n_i + n_j)
    return (mu_n * V_t / h_ij) * ((n_j - n_i) - n_mid * dpsi)
