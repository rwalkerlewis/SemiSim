"""
Position-dependent material parameters for heterojunctions (M17).

This module is the layer-4 counterpart to the scalar fields on
`semi.scaling.Scaling`. It builds per-cell DG0 fields for chi, Eg, Nc,
Nv, n_i, and epsilon_r given a mesh, cell-region tags, and the
regions configuration. The Poisson and DD form builders consume these
fields via the `heterojunction_fields` keyword (M17 Phase C). The
ohmic-contact equilibrium psi calculation in `semi.bcs` consumes the
local-material chi from the same field set (M17 Phase D).

Per ADR 0016 (Heterojunction-aware Slotboom and ohmic equilibrium),
the substitution rule extends from the single-material Boltzmann
relation `n = n_i exp((psi - phi_n) / V_t)` to its position-
dependent counterpart `n = n_i(x) exp((psi - phi_n) / V_t)` without
changing the continuity-flux shape (ADR 0004 preserved). Composes
orthogonally with the M16.4 Fermi-Dirac prefactor `gamma_n`: the two
extensions enter the same substitution slot, so a heterojunction
under FD statistics multiplies `gamma_n(x) * n_i(x)` per cell.

The field construction reduces to a constant for single-material
configurations (every region with the same material and no
material_overrides), so existing benchmarks and PR #85 examples are
bit-identical to v0.23.0 when the runner threads
`heterojunction_fields` through the form builders.
"""
from __future__ import annotations

import dataclasses
import math
from typing import TYPE_CHECKING, Any

from ..constants import cm3_to_m3, thermal_voltage
from ..materials import MATERIALS, Material, get_material

if TYPE_CHECKING:  # pragma: no cover
    from ..scaling import Scaling


_OVERRIDE_TO_FIELD: dict[str, str] = {
    "chi_eV": "chi",
    "Eg_eV": "Eg",
    "Nc_per_cm3": "Nc",
    "Nv_per_cm3": "Nv",
}


def cfg_uses_heterojunction(regions_cfg: dict[str, Any]) -> bool:
    """Return True iff any region opts into the M17 heterojunction path.

    A region opts in by setting `material_overrides` to a non-empty
    dict or `heterojunction: true`. Configurations that do neither are
    bit-identical to v0.23.0 because the runner skips
    `build_dg0_material_fields` and the form builders run the existing
    scalar `ni_hat` Constant path. M17 schema 2.8.0.
    """
    if not isinstance(regions_cfg, dict):
        return False
    for region_cfg in regions_cfg.values():
        if not isinstance(region_cfg, dict):
            continue
        overrides = region_cfg.get("material_overrides")
        if overrides:
            return True
        if region_cfg.get("heterojunction", False):
            return True
    return False


def _resolve_region_material(region_cfg: dict[str, Any]) -> Material:
    """Apply `material_overrides` to the database material for a region.

    Returns a new Material instance with override values substituted in;
    the global `MATERIALS` dict is unchanged. When no overrides are
    set, returns the database material directly so single-material
    configs that do not exercise any override path are byte-identical
    to v0.23.0 (the same Material reference, not a copy with the same
    fields).
    """
    name = region_cfg["material"]
    if name not in MATERIALS:
        raise KeyError(
            f"Unknown material {name!r} on region; known materials are "
            f"{sorted(MATERIALS.keys())}"
        )
    base = MATERIALS[name]
    overrides = region_cfg.get("material_overrides")
    if not overrides:
        return base

    fields = {f.name: getattr(base, f.name) for f in dataclasses.fields(base)}
    nci_relevant = False
    for ov_key, value in overrides.items():
        attr = _OVERRIDE_TO_FIELD.get(ov_key)
        if attr is None:
            continue
        if ov_key in ("Nc_per_cm3", "Nv_per_cm3"):
            fields[attr] = cm3_to_m3(float(value))
            nci_relevant = True
        elif ov_key == "Eg_eV":
            fields["Eg"] = float(value)
            nci_relevant = True
        else:
            fields[attr] = float(value)

    resolved = Material(**fields)
    if nci_relevant and resolved.Nc > 0.0 and resolved.Nv > 0.0 and resolved.Eg > 0.0:
        # When the user overrides any of Nc / Nv / Eg, the empirical
        # `n_i` field (a tabulated 300 K database value) is no longer
        # consistent with the band-edge parameters. Recompute n_i via
        # the closed-form `sqrt(Nc Nv) exp(-Eg / (2 V_t_300))` so the
        # heterojunction field builder sees a self-consistent material.
        # Pure-Python overrides at 300 K; the field builder
        # subsequently rescales to the runtime temperature via
        # `Material.n_i_at_T(T)` if the runner's temperature differs.
        # The `n_i` slot is recomputed at 300 K only as a fallback for
        # consumers that read `Material.n_i` directly; the form
        # builders read the per-cell `n_i_hat` field which is always
        # built from `n_i_at_T(T_runner)`.
        V_t_300 = thermal_voltage(300.0)
        resolved = dataclasses.replace(
            resolved,
            n_i=math.sqrt(resolved.Nc * resolved.Nv)
            * math.exp(-resolved.Eg / (2.0 * V_t_300)),
        )
    return resolved


def build_dg0_material_fields(
    mesh,
    cell_tags,
    regions_cfg: dict[str, Any],
    sc: Scaling,
    T: float,
    *,
    fdim: int | None = None,
) -> dict[str, Any]:
    """Build per-cell DG0 material fields for the heterojunction kernel.

    Parameters
    ----------
    mesh : dolfinx.mesh.Mesh
        The mesh the form builders run on.
    cell_tags : dolfinx.mesh.MeshTags or None
        Per-cell region tags. When None (single-region builtin meshes),
        the fields collapse to the reference-material values; this is
        the v0.23.0 byte-identity branch.
    regions_cfg : dict
        The validated `cfg["regions"]` map. Each region carries
        `material`, optional `material_overrides`, and optional
        `heterojunction` flag (M17 schema 2.8.0).
    sc : semi.scaling.Scaling
        The reference scaling object. The dimensionless internal units
        consumed by the form builders are derived from `sc`: chi_hat
        and Eg_hat divide by `sc.V0`, Nc_hat / Nv_hat / n_i_hat divide
        by `sc.C0`, and `epsilon_r_hat` is already dimensionless.
    T : float
        Device temperature in K. Threaded into the per-cell n_i field
        via `Material.n_i_at_T(T)` when `material_overrides` change
        the band-edge parameters; otherwise falls back to the
        empirical `Material.n_i` for byte-identity with v0.23.0.
    fdim : int, optional
        Reserved; not used in the present implementation. Future
        facet-restricted variants may consume it.

    Returns
    -------
    dict
        Keys: "chi_hat", "Eg_hat", "Nc_hat", "Nv_hat", "n_i_hat",
        "epsilon_r", "chi_ref_hat", "n_i_ref_hat". Each value other
        than the two `_ref_*` scalars is a `dolfinx.fem.Function` on
        the cellwise DG0 space. The `_ref_*` values are scalars
        (Python floats) carrying the reference-material values for
        the BC layer's `chi_local - chi_ref` shift (Phase D).

    Notes
    -----
    Per ADR 0016, the substitution rule for the carrier density is
    `n = n_i(x) * exp((psi - phi_n) / V_t)`. The Poisson and DD form
    builders read `n_i_hat` from this dict and thread it through the
    Slotboom helpers. For configurations without `material_overrides`
    on any region, every per-cell value collapses to a single
    constant (one for each material in the regions map); for a
    single-material config, the cellwise field equals the scalar
    `Material.n_i / sc.C0` everywhere, which is bit-identical to the
    pre-M17 path.
    """
    import dolfinx.fem as dfem
    from petsc4py import PETSc

    V_DG0 = dfem.functionspace(mesh, ("DG", 0))

    # Resolve every region's material once. tag -> resolved Material.
    tag_to_mat: dict[int, Material] = {}
    for region_cfg in regions_cfg.values():
        if not isinstance(region_cfg, dict) or "tag" not in region_cfg:
            continue
        tag_to_mat[int(region_cfg["tag"])] = _resolve_region_material(region_cfg)

    # The reference material is read off `Scaling`. The reference chi /
    # n_i are the v0.23.0 single-material values; the BC layer subtracts
    # them from the local-material values to compute psi shifts under
    # heterojunction band alignment (Phase D).
    chi_ref_eV = _resolve_reference_chi_eV(regions_cfg, sc)
    n_i_ref_m3 = float(sc.n_i)
    chi_ref_hat = chi_ref_eV / sc.V0
    n_i_ref_hat = n_i_ref_m3 / sc.C0

    # n_i fallback: when a region has no material_overrides, use the
    # empirical `Material.n_i` (database 300 K value) so single-material
    # configs are byte-identical to v0.23.0. When overrides changed
    # Nc / Nv / Eg, the resolved material's `n_i` was already
    # recomputed in `_resolve_region_material`; we still want the
    # temperature-scaling path for runners at T != 300 K, so always
    # call `n_i_at_T(T)` and only fall back to the empirical `n_i`
    # when the formula is undefined (insulators) or when the runner
    # temperature is exactly 300 K and no override has changed the
    # band-edge parameters (the byte-identity case).
    def _n_i_for_material(mat: Material, region_cfg: dict[str, Any]) -> float:
        has_overrides = bool(region_cfg.get("material_overrides"))
        if not has_overrides and abs(T - 300.0) < 1.0e-9:
            return float(mat.n_i)
        formula = mat.n_i_at_T(T)
        if formula > 0.0:
            return formula
        return float(mat.n_i)

    # Match each tag back to its region cfg so the n_i fallback can
    # detect overrides per region.
    tag_to_region_cfg: dict[int, dict[str, Any]] = {}
    for region_cfg in regions_cfg.values():
        if not isinstance(region_cfg, dict) or "tag" not in region_cfg:
            continue
        tag_to_region_cfg[int(region_cfg["tag"])] = region_cfg

    def _make_field(value_for_tag, default: float = 0.0):
        fn = dfem.Function(V_DG0)
        fn.x.array[:] = PETSc.ScalarType(default)
        for tag, mat in tag_to_mat.items():
            cells = cell_tags.find(tag)
            if cells.size == 0:
                continue
            val = value_for_tag(tag, mat)
            scalar_val = PETSc.ScalarType(float(val))
            for c in cells:
                dof = V_DG0.dofmap.cell_dofs(int(c))[0]
                fn.x.array[dof] = scalar_val
        fn.x.scatter_forward()
        return fn

    # When cell_tags is None (single-region builtin), the entire mesh
    # uses the reference material's parameters. This collapses every
    # field to a constant matching the scalar single-region path.
    if cell_tags is None:
        ref_mat: Material | None = None
        for region_cfg in regions_cfg.values():
            if isinstance(region_cfg, dict) and "material" in region_cfg:
                ref_mat = _resolve_region_material(region_cfg)
                break
        if ref_mat is None:
            ref_mat = MATERIALS["Si"]

        def _const_field(value: float):
            fn = dfem.Function(V_DG0)
            fn.x.array[:] = PETSc.ScalarType(float(value))
            return fn

        return {
            "chi_hat": _const_field(ref_mat.chi / sc.V0),
            "Eg_hat": _const_field(ref_mat.Eg / sc.V0),
            "Nc_hat": _const_field(ref_mat.Nc / sc.C0 if ref_mat.Nc > 0.0 else 0.0),
            "Nv_hat": _const_field(ref_mat.Nv / sc.C0 if ref_mat.Nv > 0.0 else 0.0),
            "n_i_hat": _const_field(float(sc.n_i) / sc.C0),
            "epsilon_r": _const_field(float(ref_mat.epsilon_r)),
            "chi_ref_hat": chi_ref_hat,
            "n_i_ref_hat": n_i_ref_hat,
        }

    chi_hat = _make_field(lambda tag, mat: mat.chi / sc.V0)
    Eg_hat = _make_field(lambda tag, mat: mat.Eg / sc.V0)
    Nc_hat = _make_field(
        lambda tag, mat: (mat.Nc / sc.C0 if mat.Nc > 0.0 else 0.0)
    )
    Nv_hat = _make_field(
        lambda tag, mat: (mat.Nv / sc.C0 if mat.Nv > 0.0 else 0.0)
    )
    n_i_hat = _make_field(
        lambda tag, mat: _n_i_for_material(mat, tag_to_region_cfg[tag]) / sc.C0
    )
    eps_r = _make_field(lambda tag, mat: float(mat.epsilon_r), default=1.0)

    return {
        "chi_hat": chi_hat,
        "Eg_hat": Eg_hat,
        "Nc_hat": Nc_hat,
        "Nv_hat": Nv_hat,
        "n_i_hat": n_i_hat,
        "epsilon_r": eps_r,
        "chi_ref_hat": chi_ref_hat,
        "n_i_ref_hat": n_i_ref_hat,
    }


def _resolve_reference_chi_eV(
    regions_cfg: dict[str, Any], sc: Scaling,
) -> float:
    """Return the chi (eV) of the reference material used to build `sc`.

    Mirrors `make_scaling_from_config`: the reference material is the
    first semiconductor region in the cfg's region map. Insulators are
    skipped because the reference for the scaling is the carrier-
    bearing material. Falls back to silicon when no semiconductor
    region is found (the same fallback the rest of the engine uses).
    """
    for region_cfg in regions_cfg.values():
        if not isinstance(region_cfg, dict):
            continue
        role = region_cfg.get("role", "semiconductor")
        if role != "semiconductor":
            continue
        mat = get_material(region_cfg["material"])
        if mat.chi > 0.0:
            return float(mat.chi)
    return float(MATERIALS["Si"].chi)
