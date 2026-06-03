"""
Equilibrium smoke verifier for the M19 precursor benchmark
`benchmarks/mosfet_3d_eq` (ADR 0019).

Five geometric / shape-based gates at V_GS = V_DS = V_BS = 0 V:

1. `n_max` in the source / drain implant volumes is within a
   factor of 3 of the N_D peak (1e20 cm^-3 = 1e26 m^-3); the
   factor accounts for Gaussian-profile averaging across the
   DOFs that fall inside the implant bounding box.
2. `n_min` in the channel (DOFs within 50 nm below the oxide,
   centred on the channel x-window) is below 1e12 cm^-3 (1e18
   m^-3) at V_GS = 0 V.
3. `p_max` in the body bulk (DOFs deeper than 500 nm below the
   surface and away from the depletion regions) is within 10 %
   of N_A (1e17 cm^-3 = 1e23 m^-3).
4. Charge neutrality |n - p - N_D + N_A| per DOF is below 1e15
   cm^-3 (1e21 m^-3) in absolute terms across the body bulk.
5. The built-in junction potential is set up: psi at the lateral
   edges of the source / drain n+ regions sits at least 3 * V_t
   below psi_body (i.e. the depletion regions under the n+/p
   junctions have formed).

ADR 0019 records the V&V scope departure: no MMS variant, no
analytical I-V comparison. These gates are precursor smoke
checks; the full Pao-Sah analytical comparison is M19 proper.
"""
from __future__ import annotations

from typing import Any

# Geometry constants matching the Phase A .geo (in meters).
L_X       = 1.5e-6
L_Y       = 5.0e-7
L_Z       = 2.0e-6
X_G_LO    = 6.25e-7   # gate / oxide left x
X_G_HI    = 8.75e-7   # gate / oxide right x
X_SRC_HI  = 5.0e-7    # source implant right x
X_DRN_LO  = 1.0e-6    # drain implant left x

# Doping reference values (Phase D config; cm^-3).
N_D_PEAK_CM3 = 1.0e20
N_A_BODY_CM3 = 1.0e17

# Physical-mask depths (meters).
IMPLANT_Z_DEPTH = 1.0e-7    # implant peak within 100 nm of top
CHANNEL_Z_DEPTH = 5.0e-8    # channel zone within 50 nm of oxide
BODY_DEEP_Z_TOP = 1.5e-6    # body bulk: DOFs with z < 1.5 um (so
                            # ~500 nm below the surface)

# Gate values.
N_MAX_FACTOR_OF_PEAK = 3.0
N_CHAN_MAX_CM3       = 1.0e12
P_BODY_REL_TOL       = 0.10
NEUTRALITY_ABS_CM3   = 1.0e15
PSI_DEPLETION_FACTOR = 3.0   # in units of V_t


def verify_mosfet_3d_eq(result: Any) -> list[tuple[str, bool, str]]:
    """Run the five equilibrium gates on the run result."""
    import numpy as np

    from semi.constants import cm3_to_m3
    from semi.materials import get_material

    cfg = result.cfg
    sc = result.scaling
    V_t = sc.V0

    x_dof = result.x_dof
    psi = result.psi_phys
    n = result.n_phys
    p = result.p_phys

    if (x_dof is None or psi is None or n is None or p is None):
        return [(
            "mosfet_3d_eq: result populates psi / n / p / x_dof",
            False,
            "missing one of result.x_dof, psi_phys, n_phys, p_phys",
        )]

    if x_dof.shape[1] != 3:
        return [(
            "mosfet_3d_eq: 3D mesh",
            False,
            f"x_dof.shape[1] = {x_dof.shape[1]}, expected 3",
        )]

    region_si = cfg["regions"]["si_body"]
    mat_si = get_material(region_si["material"])
    n_i = mat_si.n_i
    N_D_peak = cm3_to_m3(N_D_PEAK_CM3)
    N_A      = cm3_to_m3(N_A_BODY_CM3)
    n_chan_max = cm3_to_m3(N_CHAN_MAX_CM3)
    neutrality_abs = cm3_to_m3(NEUTRALITY_ABS_CM3)

    x = x_dof[:, 0]
    y = x_dof[:, 1]
    z = x_dof[:, 2]

    # Build per-DOF doping arrays from the JSON profile (sum of
    # uniform body + two Gaussian implants per Phase D config).
    N_D_per_dof = np.zeros_like(x)
    N_A_per_dof = np.zeros_like(x)
    for entry in cfg["doping"]:
        prof = entry["profile"]
        kind = prof["type"]
        if kind == "uniform":
            N_D_per_dof += cm3_to_m3(prof["N_D"])
            N_A_per_dof += cm3_to_m3(prof["N_A"])
        elif kind == "gaussian":
            cx, cy, cz = prof["center"]
            sx, sy, sz = prof["sigma"]
            peak_m3 = cm3_to_m3(prof["peak"])
            r2 = ((x - cx) / sx) ** 2 + ((y - cy) / sy) ** 2 + ((z - cz) / sz) ** 2
            contrib = peak_m3 * np.exp(-0.5 * r2)
            if prof["dopant"] == "donor":
                N_D_per_dof += contrib
            else:
                N_A_per_dof += contrib
        else:
            return [(
                "mosfet_3d_eq: known doping profile types",
                False,
                f"unsupported profile type {kind!r}",
            )]

    # Region masks (DOF-level; centroid-free, suitable for nodal
    # P1 unknowns).
    src_mask = (x <= X_SRC_HI) & (z >= L_Z - IMPLANT_Z_DEPTH)
    drn_mask = (x >= X_DRN_LO) & (z >= L_Z - IMPLANT_Z_DEPTH)
    chan_mask = ((x >= X_G_LO) & (x <= X_G_HI)
                 & (z >= L_Z - CHANNEL_Z_DEPTH) & (z <= L_Z))
    body_mask = z <= BODY_DEEP_Z_TOP

    checks: list[tuple[str, bool, str]] = []

    # Gate 1: n_max in S/D implant volumes within factor 3 of N_D peak.
    sd_mask = src_mask | drn_mask
    n_sd_max = float(n[sd_mask].max(initial=0.0)) if sd_mask.any() else 0.0
    lo = N_D_peak / N_MAX_FACTOR_OF_PEAK
    hi = N_D_peak * N_MAX_FACTOR_OF_PEAK
    checks.append((
        f"mosfet_3d_eq: n_max in S/D implants within factor "
        f"{N_MAX_FACTOR_OF_PEAK:.0f} of N_D peak ({N_D_PEAK_CM3:.1e} cm^-3)",
        lo <= n_sd_max <= hi,
        f"n_sd_max = {n_sd_max:.3e} m^-3 vs N_D peak = {N_D_peak:.3e} m^-3 "
        f"(window [{lo:.3e}, {hi:.3e}])",
    ))

    # Gate 2: n_min in the channel below N_CHAN_MAX_CM3.
    n_chan_min = (float(n[chan_mask].min(initial=np.inf))
                  if chan_mask.any() else float("inf"))
    chan_count = int(chan_mask.sum())
    checks.append((
        f"mosfet_3d_eq: n in channel zone has DOFs below "
        f"{N_CHAN_MAX_CM3:.1e} cm^-3 at V_GS = 0 V",
        chan_count > 0 and n_chan_min < n_chan_max,
        f"channel DOF count = {chan_count}; min n there = "
        f"{n_chan_min:.3e} m^-3 (threshold {n_chan_max:.3e} m^-3)",
    ))

    # Gate 3: p_max in body bulk within 10% of N_A.
    p_body_max = (float(p[body_mask].max(initial=0.0))
                  if body_mask.any() else 0.0)
    rel_err_p = abs(p_body_max - N_A) / N_A
    checks.append((
        f"mosfet_3d_eq: p_max in body bulk within "
        f"{int(P_BODY_REL_TOL * 100)} % of N_A ({N_A_BODY_CM3:.1e} cm^-3)",
        rel_err_p <= P_BODY_REL_TOL,
        f"p_body_max = {p_body_max:.3e} m^-3 vs N_A = {N_A:.3e} m^-3 "
        f"(rel_err {rel_err_p:.2%})",
    ))

    # Gate 4: charge neutrality |n - p - N_D + N_A| < neutrality_abs
    # in the body bulk. Compute the absolute residual per body DOF.
    rho_body = (n[body_mask] - p[body_mask]
                - N_D_per_dof[body_mask] + N_A_per_dof[body_mask])
    rho_max = float(np.abs(rho_body).max(initial=0.0))
    checks.append((
        f"mosfet_3d_eq: charge neutrality |n - p - N_D + N_A| "
        f"< {NEUTRALITY_ABS_CM3:.1e} cm^-3 in body bulk",
        rho_max < neutrality_abs,
        f"worst |rho| = {rho_max:.3e} m^-3 across {int(body_mask.sum())} "
        f"body-bulk DOFs (threshold {neutrality_abs:.3e} m^-3)",
    ))

    # Gate 5: built-in junction. The peak psi inside the source /
    # drain n+ implants sits at least PSI_DEPLETION_FACTOR * V_t
    # above psi_body, meaning the n+/p built-in potential is set
    # up (n+ side at +V_t ln(N_D/n_i); p body at -V_t ln(N_A/n_i);
    # difference ~ V_t ln(N_D N_A / n_i^2)). psi_body is the
    # median psi over the body bulk.
    psi_body = float(np.median(psi[body_mask]))

    psi_src_max = (float(psi[src_mask].max(initial=-np.inf))
                   if src_mask.any() else float("-inf"))
    psi_drn_max = (float(psi[drn_mask].max(initial=-np.inf))
                   if drn_mask.any() else float("-inf"))
    delta_threshold = PSI_DEPLETION_FACTOR * V_t

    src_ok = (psi_src_max - psi_body) >= delta_threshold
    drn_ok = (psi_drn_max - psi_body) >= delta_threshold

    checks.append((
        f"mosfet_3d_eq: built-in junction sets up psi in S/D "
        f"implants by >= {PSI_DEPLETION_FACTOR:.0f} * V_t above psi_body",
        src_ok and drn_ok,
        f"psi_body = {psi_body:.3f} V; "
        f"psi_src_max = {psi_src_max:.3f} V "
        f"(delta {psi_src_max - psi_body:.3f} V); "
        f"psi_drn_max = {psi_drn_max:.3f} V "
        f"(delta {psi_drn_max - psi_body:.3f} V); "
        f"threshold = {delta_threshold:.3f} V "
        f"(V_t = {V_t:.4f} V, n_i = {n_i:.3e} m^-3)",
    ))

    return checks
