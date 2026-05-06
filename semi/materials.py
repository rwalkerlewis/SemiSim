"""
Material parameter database for semiconductors and insulators.

All numbers are room-temperature unless noted. Sources cited inline.

n_i at 300 K for Si: 1.0e10 cm^-3 (Altermatt 2003 consensus; the older
textbook value 1.45e10 is superseded and disagrees with modern TCAD defaults
such as Sentaurus/Synopsys).
"""
from __future__ import annotations

from dataclasses import dataclass

from .constants import EPS0, cm2_to_m2, cm3_to_m3


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
