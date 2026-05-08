"""
Slotboom variable helpers.

Under Boltzmann statistics the carrier densities are

    n = n_i exp((psi - phi_n) / V_t)
    p = n_i exp((phi_p - psi) / V_t)

We solve for the triple (psi, phi_n, phi_p) instead of (psi, n, p) so the
continuity equations become coercive on the quasi-Fermi potentials (see
ADR 0004). Inside `semi/physics/` everything is in scaled units where
psi, phi_n, phi_p are divided by V_t and densities are divided by C_0;
`n_i_hat = n_i / C_0` is the only material-specific constant that shows
up here.

Under Fermi-Dirac (M16.4, gated by `physics.statistics: "fermi_dirac"`),
the substitution rule keeps the Slotboom shape and grows a smooth
bounded prefactor (see `semi.physics.statistics` module docstring for
the derivation):

    n = n_i * gamma_n(eta_n) * exp(psi - phi_n)         (scaled)
    p = n_i * gamma_p(eta_p) * exp(phi_p - psi)         (scaled)

with `gamma_n(eta) = 1 / (1 + 0.27 * exp(eta))` (Blakemore basic) and
`eta_n = (psi - phi_n)/V_t + eta_offset_n`. Boltzmann is recovered as
gamma -> 1, which the closed Blakemore form approaches smoothly as
`eta -> -inf`. The continuity-row shape (`J = -q mu n grad(phi)`) is
unchanged: the Einstein-factor cancellation in the generalized-Slotboom
current expression (see ADR 0004 derivation note) means no FD-specific
correction enters the residual beyond the `n` / `p` substitution.

UFL expressions in this module work with either dolfinx Functions or
fem.Constants for `psi`, `phi_n`, `phi_p`. Under M17 the `n_i_hat`
argument may also be a position-dependent `dolfinx.fem.Function`
(typically a DG0 cellwise field built by
`semi.physics.heterojunction.build_dg0_material_fields`); the
generalized-Slotboom substitution `n = n_i(x) * exp((psi - phi_n) /
V_t)` extends without changing the continuity-flux shape (ADR 0016
preserves ADR 0004). The NumPy helpers are pure Python and are
imported for tests and post-processing. This module must remain
independent of the FE space construction; it builds expressions, not
Forms.
"""
from __future__ import annotations

import numpy as np

from .statistics import _BLAKEMORE_BASIC_OFFSET, gamma_n_blakemore, gamma_p_blakemore


def _is_fermi_dirac(statistics_cfg):
    """True iff statistics_cfg dispatches to the Fermi-Dirac branch."""
    if statistics_cfg is None:
        return False
    return statistics_cfg.get("statistics", "boltzmann") == "fermi_dirac"


def n_from_slotboom(
    psi_hat, phi_n_hat, n_i_hat,
    *, statistics_cfg=None, eta_offset_n=None,
):
    """
    Scaled electron density.

    Boltzmann (default, statistics_cfg=None or "boltzmann"):

        n_hat = n_i_hat * exp(psi_hat - phi_n_hat)

    Fermi-Dirac (statistics_cfg["statistics"] == "fermi_dirac"):

        n_hat = n_i_hat * gamma_n(eta_n) * exp(psi_hat - phi_n_hat)
        eta_n = (psi_hat - phi_n_hat) + eta_offset_n
        gamma_n_blakemore(eta) = 1 / (1 + 0.27 * exp(eta))

    All arguments are UFL-compatible (Functions, Constants, expressions)
    or plain floats / numpy arrays. Returns an expression of the same
    flavor.

    Parameters
    ----------
    psi_hat, phi_n_hat, n_i_hat
        Scaled primary unknowns and intrinsic density.
    statistics_cfg : dict, optional
        Carrier-statistics dispatch; `None` (default) is bit-identical
        to pre-M16.4 (Boltzmann).
    eta_offset_n : float, optional
        Per-material reduced-Fermi offset `ln(n_i / N_C)`. Required
        when `statistics_cfg["statistics"] == "fermi_dirac"`; ignored
        on the Boltzmann path.
    """
    import ufl
    boltzmann = n_i_hat * ufl.exp(psi_hat - phi_n_hat)
    if not _is_fermi_dirac(statistics_cfg):
        return boltzmann
    if eta_offset_n is None:
        raise ValueError(
            "n_from_slotboom: statistics_cfg requests fermi_dirac but "
            "eta_offset_n is None. Pass Scaling.eta_offset_n through "
            "the form builder."
        )
    eta = (psi_hat - phi_n_hat) + ufl.as_ufl(float(eta_offset_n))
    gamma = ufl.as_ufl(1.0) / (
        ufl.as_ufl(1.0)
        + ufl.as_ufl(_BLAKEMORE_BASIC_OFFSET) * ufl.exp(eta)
    )
    return gamma * boltzmann


def p_from_slotboom(
    psi_hat, phi_p_hat, n_i_hat,
    *, statistics_cfg=None, eta_offset_p=None,
):
    """
    Scaled hole density. Boltzmann default; Fermi-Dirac dispatch by
    analogy with :func:`n_from_slotboom` using the hole-side reduced
    Fermi level `eta_p = (phi_p_hat - psi_hat) + eta_offset_p`.
    """
    import ufl
    boltzmann = n_i_hat * ufl.exp(phi_p_hat - psi_hat)
    if not _is_fermi_dirac(statistics_cfg):
        return boltzmann
    if eta_offset_p is None:
        raise ValueError(
            "p_from_slotboom: statistics_cfg requests fermi_dirac but "
            "eta_offset_p is None. Pass Scaling.eta_offset_p through "
            "the form builder."
        )
    eta = (phi_p_hat - psi_hat) + ufl.as_ufl(float(eta_offset_p))
    gamma = ufl.as_ufl(1.0) / (
        ufl.as_ufl(1.0)
        + ufl.as_ufl(_BLAKEMORE_BASIC_OFFSET) * ufl.exp(eta)
    )
    return gamma * boltzmann


def n_from_slotboom_np(
    psi_hat, phi_n_hat, n_i_hat,
    *, statistics_cfg=None, eta_offset_n=None,
):
    """NumPy counterpart of :func:`n_from_slotboom` for post-processing."""
    boltzmann = n_i_hat * np.exp(psi_hat - phi_n_hat)
    if not _is_fermi_dirac(statistics_cfg):
        return boltzmann
    if eta_offset_n is None:
        raise ValueError(
            "n_from_slotboom_np: statistics_cfg requests fermi_dirac "
            "but eta_offset_n is None."
        )
    drive = np.asarray(psi_hat, dtype=float) - np.asarray(phi_n_hat, dtype=float)
    return gamma_n_blakemore(drive, eta_offset_n) * boltzmann


def p_from_slotboom_np(
    psi_hat, phi_p_hat, n_i_hat,
    *, statistics_cfg=None, eta_offset_p=None,
):
    """NumPy counterpart of :func:`p_from_slotboom` for post-processing."""
    boltzmann = n_i_hat * np.exp(phi_p_hat - psi_hat)
    if not _is_fermi_dirac(statistics_cfg):
        return boltzmann
    if eta_offset_p is None:
        raise ValueError(
            "p_from_slotboom_np: statistics_cfg requests fermi_dirac "
            "but eta_offset_p is None."
        )
    drive = np.asarray(phi_p_hat, dtype=float) - np.asarray(psi_hat, dtype=float)
    return gamma_p_blakemore(drive, eta_offset_p) * boltzmann


def phi_n_from_np(psi_hat, n_hat, n_i_hat):
    """
    Recover the electron quasi-Fermi potential from (psi, n).

        phi_n = psi - log(n / n_i)  (all scaled by V_t and C_0).

    Guards against non-positive densities by flooring at a tiny value.
    """
    n_safe = np.maximum(n_hat, 1.0e-300)
    return psi_hat - np.log(n_safe / n_i_hat)


def phi_p_from_np(psi_hat, p_hat, n_i_hat):
    """Recover the hole quasi-Fermi potential from (psi, p)."""
    p_safe = np.maximum(p_hat, 1.0e-300)
    return psi_hat + np.log(p_safe / n_i_hat)


def equilibrium_psi_hat(N_net, n_i):
    """
    Scaled equilibrium potential for local charge neutrality, psi/V_t.

        psi_eq = V_t * asinh(N_net / (2 n_i))

    Works on scalars or numpy arrays. `N_net` and `n_i` must be in the
    same (physical) units; the returned quantity is dimensionless
    (already scaled by V_t).
    """
    return np.arcsinh(N_net / (2.0 * n_i))


def contact_phi_hat(V_applied_volts, V_t_volts):
    """
    Scaled quasi-Fermi potential at an ohmic contact under applied bias.

    Both phi_n and phi_p equal the applied contact bias (PHYSICS.md 3.1).
    In scaled units this is simply V_applied / V_t.
    """
    return float(V_applied_volts) / float(V_t_volts)
