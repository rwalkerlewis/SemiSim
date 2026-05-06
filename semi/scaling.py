"""
Nondimensional scaling for Poisson and drift-diffusion systems.

Why it matters. The raw semiconductor device equations have wildly
different scales: potentials ~ 1 V, carrier densities ~ 10^22 m^-3,
permittivity ~ 10^-11 F/m, charge ~ 10^-19 C. Feeding those directly
to Newton gives a Jacobian with condition number > 10^30 and the
solver fails. Scaling pulls everything to O(1).

Scales used:
    L0 = device characteristic length (m)
    V0 = thermal voltage kT/q (V)
    C0 = characteristic density, typically max|doping| (m^-3)
    mu0 = reference mobility (m^2/(V s))
    D0 = V0 * mu0 (Einstein diffusivity)
    t0 = L0^2 / D0
    J0 = q D0 C0 / L0

The dimensionless Poisson equation becomes

    -lambda^2 eps_r nabla^2 psi_hat = p_hat - n_hat + N_D_hat - N_A_hat

with lambda^2 = eps_0 V0 / (q C0 L0^2).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import EPS0, KB, M0, Q, thermal_voltage


@dataclass
class Scaling:
    """Scale factors for nondimensionalizing the device equations."""

    L0: float                    # characteristic length, m
    C0: float                    # characteristic density, m^-3
    T: float                     # temperature, K
    mu0: float                   # reference mobility, m^2/(V s)
    n_i: float                   # intrinsic density of reference material, m^-3
    # M16.4: effective conduction- and valence-band density of states
    # for the reference material, in m^-3. Optional because pre-M16.4
    # callers (and many tests) instantiate Scaling directly without an
    # FD path. The Boltzmann branch never reads them; the Fermi-Dirac
    # branch reads `eta_offset_n` / `eta_offset_p`, both of which
    # raise a clear error when N_C or N_V is None.
    N_C: float | None = None
    N_V: float | None = None
    # M16.5: thermionic-emission effective masses (relative to the
    # electron rest mass m_0). None on materials that have not been
    # characterized for Schottky contacts; the Schottky surface form
    # raises a clear error if the reference material exposes None
    # while a Schottky contact is requested. The Boltzmann ohmic /
    # gate paths never read these.
    m_n_star: float | None = None
    m_p_star: float | None = None

    @property
    def V0(self) -> float:
        """Thermal voltage, V."""
        return thermal_voltage(self.T)

    @property
    def D0(self) -> float:
        """Einstein diffusivity for reference mobility, m^2/s."""
        return self.V0 * self.mu0

    @property
    def t0(self) -> float:
        """Diffusion time over L0, s."""
        return self.L0 ** 2 / self.D0

    @property
    def J0(self) -> float:
        """Current density scale, A/m^2."""
        return Q * self.D0 * self.C0 / self.L0

    @property
    def lambda2(self) -> float:
        """
        Bare (eps_0-only) squared Debye-length-to-device-length ratio.
        Multiply by eps_r in the UFL form for the physical coefficient.
        Small values (<<1) mean the Debye length is much shorter than
        the device, which is the space-charge-dominant regime.
        """
        return EPS0 * self.V0 / (Q * self.C0 * self.L0 ** 2)

    @property
    def eta_offset_n(self) -> float:
        """
        Reduced Fermi-level offset for electrons, ln(n_i / N_C).

        Used by the M16.4 Fermi-Dirac dispatch in
        :mod:`semi.physics.statistics` to turn the Slotboom drive
        `(psi - phi_n)/V_t` into the absolute reduced Fermi level
        `eta_n`. Raises ValueError if N_C is unset (the Boltzmann
        branch never asks for this, so existing call sites are
        unaffected).
        """
        if self.N_C is None:
            raise ValueError(
                "Scaling.N_C is not populated; the Fermi-Dirac dispatch "
                "(physics.statistics='fermi_dirac') requires the reference "
                "material to expose a non-zero conduction-band density of "
                "states (Material.Nc). Boltzmann-default solves do not "
                "read this property and remain unaffected."
            )
        from .physics.statistics import _eta_offset_for_material
        return _eta_offset_for_material(self.N_C, self.n_i)

    @property
    def v_n_thermal(self) -> float:
        """
        Electron thermionic emission velocity in m/s,
        ``sqrt(kT / (2 pi m_n*))``. Si at 300 K with m_n* = 0.26 m_0
        evaluates to ~2.05e5 m/s (~2.05e7 cm/s). Used by the Schottky
        Robin surface form (M16.5; ADR 0015).
        """
        if self.m_n_star is None:
            raise ValueError(
                "Scaling.m_n_star is not populated; the Schottky "
                "thermionic-emission Robin BC requires the reference "
                "material to expose a thermionic-emission electron "
                "effective mass (Material.m_n_star). Ohmic / gate "
                "paths do not read this property."
            )
        return float(math.sqrt(KB * self.T / (2.0 * math.pi * self.m_n_star * M0)))

    @property
    def v_p_thermal(self) -> float:
        """Hole counterpart of :attr:`v_n_thermal`. Si at 300 K with
        m_p* = 0.39 m_0 evaluates to ~1.68e5 m/s.
        """
        if self.m_p_star is None:
            raise ValueError(
                "Scaling.m_p_star is not populated; the Schottky "
                "thermionic-emission Robin BC requires the reference "
                "material to expose a thermionic-emission hole "
                "effective mass (Material.m_p_star). Ohmic / gate "
                "paths do not read this property."
            )
        return float(math.sqrt(KB * self.T / (2.0 * math.pi * self.m_p_star * M0)))

    @property
    def eta_offset_p(self) -> float:
        """Hole counterpart of :attr:`eta_offset_n`, ln(n_i / N_V)."""
        if self.N_V is None:
            raise ValueError(
                "Scaling.N_V is not populated; the Fermi-Dirac dispatch "
                "(physics.statistics='fermi_dirac') requires the reference "
                "material to expose a non-zero valence-band density of "
                "states (Material.Nv). Boltzmann-default solves do not "
                "read this property and remain unaffected."
            )
        from .physics.statistics import _eta_offset_for_material
        return _eta_offset_for_material(self.N_V, self.n_i)

    def debye_length(self, n_ref: float | None = None) -> float:
        """
        Debye length in meters for a reference carrier density n_ref (m^-3).
        Uses C0 if n_ref not given.
        """
        if n_ref is None:
            n_ref = self.C0
        import math
        return math.sqrt(EPS0 * self.V0 / (Q * n_ref))

    def __repr__(self) -> str:
        return (
            f"Scaling(L0={self.L0:.2e} m, C0={self.C0:.2e} m^-3, T={self.T} K, "
            f"V0={self.V0:.4f} V, lambda2={self.lambda2:.3e})"
        )


def make_scaling_from_config(cfg: dict, reference_material) -> Scaling:
    """
    Construct a Scaling from a validated config dict and reference material.

    L0 is inferred from mesh extents (builtin) or defaults to 1 um (file-based).
    C0 is the max absolute doping across all profiles, floored at 1e16 cm^-3.
    """
    L0 = _infer_length(cfg)
    C0 = _infer_density(cfg)
    T = cfg["physics"]["temperature"]
    # M16.4: pull N_C / N_V from the reference material if the slot
    # carries a positive value. Insulators and pre-M16.4 materials
    # default to 0.0 in `semi.materials`; map zero to None so the
    # FD-dispatch property surfaces a clear error rather than a silent
    # log(0) on the Fermi-Dirac path.
    N_C = reference_material.Nc if reference_material.Nc > 0.0 else None
    N_V = reference_material.Nv if reference_material.Nv > 0.0 else None
    return Scaling(
        L0=L0, C0=C0, T=T,
        mu0=reference_material.mu_n,
        n_i=reference_material.n_i,
        N_C=N_C, N_V=N_V,
        m_n_star=reference_material.m_n_star,
        m_p_star=reference_material.m_p_star,
    )


def _infer_length(cfg: dict) -> float:
    mesh = cfg["mesh"]
    if mesh["source"] == "builtin":
        extents = mesh["extents"]
        return max(b[1] - b[0] for b in extents)
    return 1.0e-6  # fallback for file-based meshes


def _infer_density(cfg: dict) -> float:
    from .constants import cm3_to_m3
    C_max = cm3_to_m3(1.0e16)  # floor
    for entry in cfg["doping"]:
        p = entry["profile"]
        if p["type"] == "uniform":
            C_max = max(C_max, cm3_to_m3(p["N_D"]), cm3_to_m3(p["N_A"]))
        elif p["type"] == "step":
            for k in ("N_D_left", "N_A_left", "N_D_right", "N_A_right"):
                C_max = max(C_max, cm3_to_m3(p[k]))
        elif p["type"] == "gaussian":
            C_max = max(C_max, cm3_to_m3(p["peak"]))
    return C_max
