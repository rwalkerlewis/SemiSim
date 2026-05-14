"""
Physical constants (SI, 2019 redefinition values) and unit-conversion helpers.

Unit convention throughout SemiSim:
    - Lengths: meters
    - Potentials: volts
    - Densities: m^-3 internally (JSON input is cm^-3 per device-physics tradition)
    - Temperature: Kelvin
    - Mobilities: m^2/(V s) internally (JSON input is cm^2/(V s))
    - Current density: A/m^2
"""
from __future__ import annotations

# 2019 SI redefinition values
Q = 1.602176634e-19         # elementary charge, C
KB = 1.380649e-23           # Boltzmann constant, J/K
EPS0 = 8.8541878128e-12     # vacuum permittivity, F/m
HBAR = 1.054571817e-34      # reduced Planck constant, J s
M0 = 9.1093837015e-31       # electron rest mass, kg
NA = 6.02214076e23          # Avogadro constant, 1/mol
C_LIGHT = 299792458.0       # speed of light, m/s


def thermal_voltage(T: float) -> float:
    """kT/q in volts. At 300 K: ~0.02585 V."""
    return KB * T / Q


def cm3_to_m3(n_cm3: float) -> float:
    """Convert cm^-3 to m^-3."""
    return n_cm3 * 1.0e6


def m3_to_cm3(n_m3: float) -> float:
    """Convert m^-3 to cm^-3."""
    return n_m3 * 1.0e-6


def cm2_to_m2(mu_cm2: float) -> float:
    """Convert cm^2/(V s) to m^2/(V s)."""
    return mu_cm2 * 1.0e-4


def m2_to_cm2(mu_m2: float) -> float:
    """Convert m^2/(V s) to cm^2/(V s)."""
    return mu_m2 * 1.0e4
