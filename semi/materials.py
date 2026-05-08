"""
Material parameter database for semiconductors and insulators.

All numbers are room-temperature unless noted. Sources cited inline.

n_i at 300 K for Si: 1.0e10 cm^-3 (Altermatt 2003 consensus; the older
textbook value 1.45e10 is superseded and disagrees with modern TCAD defaults
such as Sentaurus/Synopsys).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import EPS0, cm2_to_m2, cm3_to_m3, thermal_voltage


@dataclass
class Material:
    """
    Material parameter set.

    Attributes use SI internal units. For insulators, semiconductor-only
    fields (Eg, chi, Nc, Nv, n_i, mobilities) are left at zero.
    """

    name: str
    role: str                     # "semiconductor" or "insulator"
    epsilon_r: float              # relative permittivity
    Eg: float = 0.0               # bandgap, eV
    chi: float = 0.0              # electron affinity, eV
    Nc: float = 0.0               # CB effective DOS, m^-3
    Nv: float = 0.0               # VB effective DOS, m^-3
    n_i: float = 0.0              # intrinsic carrier density at 300 K, m^-3
    mu_n: float = 0.0             # electron mobility, m^2/(V s)
    mu_p: float = 0.0             # hole mobility, m^2/(V s)
    # M16.5 thermionic-emission effective masses (relative to m_0).
    # Used by the Schottky surface form to compute the carrier thermal
    # velocity v_th = sqrt(kT / (2 pi m*)). These are the Sze 3rd ed
    # Table 1 thermionic-emission masses, which differ from the
    # conductivity effective masses. None on materials that have not
    # been characterized for Schottky contacts.
    m_n_star: float | None = None
    m_p_star: float | None = None

    @property
    def epsilon(self) -> float:
        """Absolute permittivity, F/m."""
        return self.epsilon_r * EPS0

    def is_semiconductor(self) -> bool:
        return self.role == "semiconductor"

    def is_insulator(self) -> bool:
        return self.role == "insulator"

    def n_i_at_T(self, T: float) -> float:
        """Closed-form intrinsic density in m^-3 at temperature T,

            n_i(T) = sqrt(Nc * Nv) * exp(-Eg / (2 V_t))

        with `Eg` in eV and `V_t = kT/q` in volts (so `Eg / V_t` is the
        textbook reduced gap). Returns 0.0 for materials whose
        semiconductor parameters are not populated (insulators set Nc,
        Nv, Eg to zero by default; the formula is undefined there and
        we want the consumer to skip the cell rather than divide by
        zero downstream).

        For Si at 300 K with the Sze 3rd ed Nc / Nv (2.86e19 /
        3.10e19 cm^-3) and Eg = 1.12 eV the formula returns
        ~1.16e10 cm^-3, which differs from the Altermatt 2003
        consensus value 1.0e10 cm^-3 stored on `n_i` by ~16 %; the
        gap reflects that the empirical n_i at 300 K is measured
        rather than computed, and the textbook Nc / Nv at 300 K is
        a secondary fit. Consumers that need the empirical value at
        T = 300 K read `n_i` directly; consumers that need
        temperature scaling (or that need an n_i consistent with
        per-region material_overrides on Eg / Nc / Nv) read
        `n_i_at_T(T)`. M17.
        """
        if self.Nc <= 0.0 or self.Nv <= 0.0 or self.Eg <= 0.0:
            return 0.0
        V_t = thermal_voltage(T)
        return math.sqrt(self.Nc * self.Nv) * math.exp(-self.Eg / (2.0 * V_t))


MATERIALS: dict[str, Material] = {
    "Si": Material(
        name="Si",
        role="semiconductor",
        epsilon_r=11.7,
        Eg=1.12,
        chi=4.05,
        Nc=cm3_to_m3(2.86e19),       # Sze, 3rd ed., Table 7
        Nv=cm3_to_m3(3.10e19),
        n_i=cm3_to_m3(1.0e10),       # Altermatt 2003
        mu_n=cm2_to_m2(1400.0),      # undoped, 300 K
        mu_p=cm2_to_m2(450.0),
        m_n_star=0.26,               # Sze 3rd ed Table 1 thermionic mass
        m_p_star=0.39,
    ),
    "Ge": Material(
        name="Ge",
        role="semiconductor",
        epsilon_r=16.0,
        Eg=0.66,
        chi=4.0,
        Nc=cm3_to_m3(1.04e19),
        Nv=cm3_to_m3(6.0e18),
        n_i=cm3_to_m3(2.0e13),
        mu_n=cm2_to_m2(3900.0),
        mu_p=cm2_to_m2(1900.0),
    ),
    "GaAs": Material(
        name="GaAs",
        role="semiconductor",
        epsilon_r=12.9,
        Eg=1.424,
        chi=4.07,
        Nc=cm3_to_m3(4.7e17),
        Nv=cm3_to_m3(7.0e18),
        n_i=cm3_to_m3(2.1e6),
        mu_n=cm2_to_m2(8500.0),
        mu_p=cm2_to_m2(400.0),
    ),
    # M17: AlGaAs at Al fraction x = 0.3 (the Stern HEMT default), 300 K.
    # Vurgaftman, Meyer, Ram-Mohan 2001 ("Band parameters for III-V
    # compound semiconductors and their alloys", J. Appl. Phys. 89,
    # 5815) gives the Al_x Ga_{1-x} As bandgap at the Gamma valley as
    # Eg(x) = 1.424 + 1.247 x for x <= 0.45 (direct gap regime), so at
    # x = 0.3, Eg = 1.798 eV. The electron affinity follows the
    # Anderson rule with a published delta-chi of 1.1 eV per unit Al
    # fraction, so chi = 4.07 - 1.1 * 0.3 = 3.74 eV. The conduction-
    # band effective DOS Nc scales with the effective mass m_n* = 0.063
    # + 0.083 x = 0.0879 (Vurgaftman 2001), giving Nc(300 K) =
    # 5.5e17 cm^-3 to two significant figures (the GaAs reference is
    # 4.7e17 cm^-3 at m_n* = 0.063). Nv is dominated by the heavy-hole
    # band; the alloy-corrected Nv ~= 1.0e19 cm^-3 (interpolated from
    # the GaAs 7.0e18 and AlAs 1.5e19 endpoints). n_i at 300 K is
    # computed from the formula sqrt(Nc Nv) exp(-Eg / (2 V_t)) and
    # rounds to ~1.85e3 cm^-3 (much smaller than GaAs's 2.1e6 cm^-3
    # because the wider gap suppresses the intrinsic density
    # exponentially). Mobilities are dominated by alloy scattering at
    # Al fractions above 0.1: mu_n drops from GaAs's 8500 cm^2/V/s to
    # ~1500 cm^2/V/s at x = 0.3, mu_p similarly drops to ~100. The
    # epsilon_r interpolates linearly between GaAs (12.9) and AlAs
    # (10.06) per Adachi to ~12.0. m_n_star = 0.0879 m_0.
    "AlGaAs_0p3": Material(
        name="AlGaAs_0p3",
        role="semiconductor",
        epsilon_r=12.0,
        Eg=1.798,
        chi=3.74,
        Nc=cm3_to_m3(5.5e17),
        Nv=cm3_to_m3(1.0e19),
        n_i=cm3_to_m3(1.85e3),
        mu_n=cm2_to_m2(1500.0),
        mu_p=cm2_to_m2(100.0),
        m_n_star=0.0879,
        m_p_star=0.51,
    ),
    "SiO2": Material(name="SiO2", role="insulator", epsilon_r=3.9),
    "HfO2": Material(name="HfO2", role="insulator", epsilon_r=25.0),
    "Si3N4": Material(name="Si3N4", role="insulator", epsilon_r=7.5),
}


def get_material(name: str) -> Material:
    """Look up a material by name. Raises KeyError if unknown."""
    if name not in MATERIALS:
        raise KeyError(
            f"Unknown material {name!r}. Known: {sorted(MATERIALS.keys())}"
        )
    return MATERIALS[name]


def list_materials() -> list[str]:
    return sorted(MATERIALS.keys())
