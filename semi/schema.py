"""
JSON schema for semiconductor simulation input.

The schema itself lives in ``schemas/input.v1.json`` (Draft-07); this module
loads it once at import time and exposes it as ``SCHEMA`` for engine and
server code.

Design principles:
    - Flat where possible, nested where it aids clarity
    - Units: lengths in meters, potentials in volts, densities in cm^-3
      (device-physics tradition; converted to m^-3 internally)
    - Mesh and external files are referenced by relative path from the JSON
    - Schema validates before any physics runs; missing required fields
      or wrong types fail fast with a clear error

Top-level structure:
    {
      "schema_version": "1.0.0",     MAJOR.MINOR.PATCH; major gated by engine
      "name": "...",                 case name, used in output paths
      "dimension": 1 | 2 | 3,
      "mesh": { ... },               how to get the mesh (builtin or file)
      "regions": { name: {...} },    material per region
      "doping": [ ... ],             list of doping profile specs
      "contacts": [ ... ],           contact BCs (ohmic, gate, insulating)
      "physics": { ... },            model selections
      "solver": { ... },             Newton/continuation parameters
      "sweep": { ... },              optional bias sweep
      "output": { ... }              output controls
    }
"""
from __future__ import annotations

import json
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Any

_SCHEMAS_DIR = Path(__file__).parent.parent / "schemas"
_SCHEMA_PATHS = {
    1: _SCHEMAS_DIR / "input.v1.json",
    2: _SCHEMAS_DIR / "input.v2.json",
}
# Default path for the legacy `SCHEMA` re-export. v1 stays the documented
# default for one minor cycle so anything that does
# `from semi.schema import SCHEMA` still gets the v1 surface; v2 callers
# should use `get_schema(2)` or pass `schema_version: "2.0.0"` and let
# `validate` dispatch.
_SCHEMA_PATH = _SCHEMA_PATHS[1]

# Major versions of the input schema this engine build accepts. Inputs whose
# `schema_version` major is not in this set are rejected by `validate`.
# v1 inputs log a `DeprecationWarning` and continue to validate against
# input.v1.json (the loose schema); v2 inputs validate against input.v2.json
# (strict, additionalProperties: false on every object node). Schemas are
# published to https://rwalkerlewis.github.io/kronos-semi/schemas/ on every
# release tag; see .github/workflows/publish-schemas.yml.
ENGINE_SUPPORTED_SCHEMA_MAJORS = (1, 2)
# Backwards-compatibility re-export: external callers (kronos_server, the
# M11 schema-versioning tests) read `ENGINE_SUPPORTED_SCHEMA_MAJOR` and
# expect a single int. Keep it at the highest supported major so server
# `/schema` and `/capabilities` endpoints advertise the strict schema.
ENGINE_SUPPORTED_SCHEMA_MAJOR = max(ENGINE_SUPPORTED_SCHEMA_MAJORS)

# Minor version tracking. Bumped when a backward-compatible schema addition
# is made.
#   M13 (1.1.0): added the `transient` solver type.
#   M14 (1.2.0): added the `ac_sweep` solver type plus solver.dc_bias and
#                solver.ac sub-objects (frequency sweep specification).
#   M14.2 (1.3.0): added the top-level `coordinate_system` field with the
#                  `axisymmetric` option for cylindrical 2D MOSCAP and
#                  similar rotationally-symmetric devices.
#   M15 (1.4.0): added solver.backend (cpu-mumps | gpu-amgx | gpu-hypre |
#                auto) and solver.compute (device, precision, preconditioner,
#                linear_solver). Defaults preserve byte-equivalent CPU-MUMPS
#                behavior; resolution of `auto` happens at solve time.
#   M14.3 (2.0.0): strict mode. additionalProperties:false on every object
#                  node so input typos fail validation rather than being
#                  silently dropped. v1 inputs continue to load with a
#                  DeprecationWarning for one minor cycle (1.x -> 2.x).
#   M16.1 (2.1.0): added physics.mobility.model: "caughey_thomas" plus
#                  vsat_n, vsat_p, beta_n, beta_p parameters for closed-
#                  form velocity-saturation mobility. Additive bump:
#                  v2.0.0 inputs continue to validate (the new keys are
#                  optional with defaults that ignore the saturation
#                  branch).
#   M16.2 (2.2.0): added physics.mobility.model: "lombardi" plus the
#                  bulk_model selector, the lombardi sub-object (B_n,
#                  B_p, C_n, C_p, lambda_n, lambda_p, delta_n, delta_p),
#                  and the interface_facet_tag (required at solve time
#                  by the loader when model=='lombardi'; the JSON
#                  schema declares it as ['integer', 'null'] and the
#                  conditional requirement is enforced in
#                  _validate_mobility_lombardi). Additive bump:
#                  v2.0.0 and v2.1.0 inputs continue to validate.
#   M16.3 (2.3.0): promoted physics.recombination.auger from a
#                  forward-compat placeholder to a real flag, and
#                  added physics.recombination.C_n and
#                  physics.recombination.C_p (Si Dziewior-Schmid
#                  defaults 2.8e-31 cm^6/s and 9.9e-32 cm^6/s). When
#                  auger==false (default), the recombination kernel is
#                  bit-identical to v0.18.0; C_n and C_p are filled to
#                  defaults but unused. Additive bump: v2.0.0, v2.1.0,
#                  and v2.2.0 inputs continue to validate.
#   M16.4 (2.4.0): widened physics.statistics enum from ["boltzmann"]
#                  to ["boltzmann", "fermi_dirac"]. Default stays
#                  "boltzmann"; the boltzmann branch is bit-identical
#                  to v0.19.0 on every benchmark. The fermi_dirac
#                  branch substitutes the Blakemore approximation in
#                  the generalized-Slotboom helpers (ADR 0004) and the
#                  Poisson source. Additive bump: v2.0.0, v2.1.0,
#                  v2.2.0, and v2.3.0 inputs continue to validate.
#   M16.5 (2.5.0): widened contacts[].type enum from
#                  ["ohmic", "gate", "insulating"] to add "schottky"
#                  for metal-semiconductor contacts; added
#                  contacts[].barrier_height_eV (number | null;
#                  required when type=='schottky'). The Schottky path
#                  applies a metal-Fermi-level Dirichlet on psi and a
#                  thermionic-emission Robin BC on the continuity
#                  rows; ohmic, gate, and insulating contacts remain
#                  bit-identical to v0.20.0. Additive bump: v2.0.0,
#                  v2.1.0, v2.2.0, v2.3.0, and v2.4.0 inputs continue
#                  to validate.
#   M16.6 (2.6.0): added the physics.tunneling sub-object with bbt
#                  and tat boolean flags plus the Kane (A_kane,
#                  B_kane) and Hurkx (tau_n_min, tau_p_min, F_kT,
#                  alpha) parameters. Both flags default to false; the
#                  kernel is bit-identical to v0.21.0 when both are
#                  off. A UserWarning fires at validate time when
#                  bbt is true and statistics is "boltzmann" (BBT is
#                  most accurate under Fermi-Dirac at heavy doping).
#                  Additive bump: v2.0.0 through v2.5.0 inputs
#                  continue to validate.
#   M16.7 (2.7.0): added contacts[].voltage_t for time-varying
#                  contact voltage in the transient runner. Two
#                  variants: {type: "table", times, values} (linear
#                  interpolation between sample points with endpoint
#                  clamping) and {type: "step", t0, v0, v1} (one
#                  transition at t0). Mutually exclusive with
#                  voltage_sweep. The bias_sweep, ac_sweep, equilibrium,
#                  mos_cv, mos_cap_ac, and resistor_3d runners reject
#                  voltage_t at validate time so users see the failure
#                  before the FEM path. v2.0.0 through v2.6.0 inputs
#                  continue to validate; configs without voltage_t are
#                  bit-identical to v0.22.0.
#   M17 (2.8.0): added regions[*].material_overrides (optional sub-
#                object with chi_eV, Eg_eV, Nc_per_cm3, Nv_per_cm3
#                fields; any subset may be set to override the value
#                from regions[*].material) and regions[*].heterojunction
#                (optional bool, default false). When set, the position-
#                dependent DG0 chi/Eg/Nc/Nv/n_i fields built by
#                semi/physics/heterojunction.build_dg0_material_fields
#                pick up the override values for that region; the
#                ohmic-contact equilibrium psi calculation in
#                semi/bcs.py reads chi from the local region's material
#                instead of the reference material. v2.0.0 through
#                v2.7.0 inputs continue to validate; configs without
#                material_overrides and without heterojunction:true are
#                bit-identical to v0.23.0 (the DG0 fields collapse to
#                the scalar single-material values).
SCHEMA_SUPPORTED_MINOR = 8


@lru_cache(maxsize=8)
def _load_schema_for_major(major: int) -> dict[str, Any]:
    path = _SCHEMA_PATHS.get(int(major))
    if path is None:
        raise ValueError(
            f"No schema bundled for major version {major!r}; "
            f"engine supports {ENGINE_SUPPORTED_SCHEMA_MAJORS}"
        )
    with open(path) as f:
        return json.load(f)


def _load_schema() -> dict[str, Any]:
    """Legacy single-schema loader; preserved for backwards compatibility."""
    return _load_schema_for_major(1)


def get_schema(major: int = ENGINE_SUPPORTED_SCHEMA_MAJOR) -> dict[str, Any]:
    """Return the JSON-schema dict for the requested input major version."""
    return _load_schema_for_major(int(major))


SCHEMA: dict[str, Any] = _load_schema()


class SchemaError(Exception):
    """Raised when a config fails JSON schema validation."""
    pass


def validate(cfg: dict[str, Any]) -> dict[str, Any]:
    """
    Validate a config dict against the schema for its declared major and
    fill defaults.

    Major dispatch:
      - `schema_version: "1.x.y"` validates against `schemas/input.v1.json`
        (loose schema; unknown keys silently ignored). A
        `DeprecationWarning` is emitted once per call so v1 callers can
        track migration to v2.
      - `schema_version: "2.x.y"` validates against `schemas/input.v2.json`
        (strict; `additionalProperties: false` on every object node, so
        typos raise instead of being silently dropped).

    Uses jsonschema if available; otherwise falls back to a minimal
    required-key check (the v2 strict gate cannot be enforced without
    jsonschema).
    """
    requested_version = cfg.get("schema_version")
    if requested_version is None:
        raise SchemaError(
            "Missing required field 'schema_version' (semver string, e.g. '2.0.0')"
        )
    try:
        requested_major = int(str(requested_version).split(".")[0])
    except (ValueError, AttributeError):
        raise SchemaError(
            f"schema_version {requested_version!r} is not a valid semver string"
        ) from None
    if requested_major not in ENGINE_SUPPORTED_SCHEMA_MAJORS:
        raise SchemaError(
            f"schema_version {requested_version!r} has major {requested_major}; "
            f"engine supports majors {ENGINE_SUPPORTED_SCHEMA_MAJORS}"
        )
    if requested_major == 1:
        warnings.warn(
            "Schema v1 is deprecated; migrate to v2.0.0 "
            "(strict mode, additionalProperties: false). v1 will be "
            "removed after one minor release cycle. "
            "See docs/schema/reference.md for the migration guide.",
            DeprecationWarning,
            stacklevel=2,
        )

    schema = _load_schema_for_major(requested_major)
    try:
        from jsonschema import Draft7Validator
        validator = Draft7Validator(schema)
        errors = sorted(validator.iter_errors(cfg), key=lambda e: list(e.path))
        if errors:
            messages = []
            for e in errors:
                path = ".".join(str(p) for p in e.absolute_path) or "<root>"
                messages.append(f"  at {path}: {e.message}")
            raise SchemaError("JSON schema validation failed:\n" + "\n".join(messages))
    except ImportError:
        missing = [k for k in schema["required"] if k not in cfg]
        if missing:
            raise SchemaError(f"Missing required fields: {missing}") from None

    return _fill_defaults(cfg)


def _validate_coordinate_system(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for the `coordinate_system` setting.

    `axisymmetric` requires:
      - dimension == 2 (the meridian half-plane is a 2D domain),
      - all radial coordinates non-negative (mesh extents r >= 0 if builtin),
      - no Dirichlet (ohmic or gate) contact pinned at r = 0; the symmetry
        axis is a natural no-flux boundary.

    These checks happen after JSON-schema structural validation so the
    error messages can reference resolved values. Failures raise
    SchemaError with a clear message.
    """
    cs = cfg.get("coordinate_system", "cartesian")
    if cs == "cartesian":
        return
    if cs != "axisymmetric":
        raise SchemaError(
            f"coordinate_system={cs!r} is not supported; "
            f"expected 'cartesian' or 'axisymmetric'"
        )

    dim = int(cfg["dimension"])
    if dim != 2:
        raise SchemaError(
            f"coordinate_system='axisymmetric' requires dimension=2 "
            f"(meridian half-plane); got dimension={dim}"
        )

    mesh = cfg.get("mesh", {})
    if mesh.get("source") == "builtin":
        extents = mesh.get("extents", [])
        if extents and float(extents[0][0]) < 0.0:
            raise SchemaError(
                "coordinate_system='axisymmetric': radial extent extents[0] "
                f"must satisfy r_min >= 0; got r_min={extents[0][0]}"
            )

    # Forbid Dirichlet contacts on the symmetry axis. We detect axis facets
    # by looking for facets_by_plane entries with axis=0 and value=0 in
    # builtin meshes; for file meshes the user is responsible for not
    # tagging the axis as a contact (documented in the schema description).
    axis_facet_names: set[str] = set()
    for plane in mesh.get("facets_by_plane", []) if mesh.get("source") == "builtin" else []:
        if int(plane.get("axis", -1)) == 0 and abs(float(plane.get("value", 1.0))) < 1.0e-15:
            axis_facet_names.add(plane["name"])
    for c in cfg.get("contacts", []):
        if c.get("type") in ("ohmic", "gate") and c.get("facet") in axis_facet_names:
            raise SchemaError(
                f"coordinate_system='axisymmetric': contact {c.get('name')!r} "
                f"is pinned at r=0 on facet {c.get('facet')!r}; the symmetry "
                f"axis must remain a natural (Neumann) boundary"
            )


def _validate_compute(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for `solver.backend` against `solver.compute.device`.

    Mirrors the cross-field pattern used for `coordinate_system`. Rules:
      - `backend == "cpu-mumps"` requires `compute.device in {"auto", "cpu"}`.
      - `backend in {"gpu-amgx", "gpu-hypre"}` requires
        `compute.device in {"auto", "cuda", "hip"}`.
      - `backend == "auto"` is always accepted at validation time; the
        actual device is resolved at solve time by `semi.compute`
        (introduced in Phase C of M15).

    No-op when neither `backend` nor `compute` is present (default-fill
    sets `backend="cpu-mumps"` and leaves `compute` absent, which is the
    byte-equivalent legacy path).
    """
    solver = cfg.get("solver", {})
    backend = solver.get("backend")
    compute = solver.get("compute")
    if backend is None and compute is None:
        return
    backend = backend or "cpu-mumps"
    device = (compute or {}).get("device", "cpu" if backend == "cpu-mumps" else "auto")

    if backend == "cpu-mumps":
        if device not in ("auto", "cpu"):
            raise SchemaError(
                f"solver.backend='cpu-mumps' requires solver.compute.device "
                f"in {{'auto', 'cpu'}}; got {device!r}"
            )
    elif backend in ("gpu-amgx", "gpu-hypre"):
        if device not in ("auto", "cuda", "hip"):
            raise SchemaError(
                f"solver.backend={backend!r} requires solver.compute.device "
                f"in {{'auto', 'cuda', 'hip'}}; got {device!r}"
            )
    elif backend == "auto":
        return
    else:
        raise SchemaError(
            f"solver.backend={backend!r} is not a recognised backend; "
            f"expected one of {{'cpu-mumps', 'gpu-amgx', 'gpu-hypre', 'auto'}}"
        )


def _validate_mobility_lombardi(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for `physics.mobility.model='lombardi'`.

    The JSON schema declares `interface_facet_tag` as `['integer',
    'null']` so v2.0.0 and v2.1.0 inputs (no Lombardi block) keep
    validating. The schema cannot express "required only when
    model=='lombardi'" without conditional logic; we enforce that
    here so users see the failure at `validate(cfg)` time rather than
    deep inside the FEM path. M16.2.
    """
    mob = cfg.get("physics", {}).get("mobility", {}) or {}
    if mob.get("model") != "lombardi":
        return
    if mob.get("interface_facet_tag") is None:
        raise SchemaError(
            "physics.mobility.model='lombardi' requires "
            "physics.mobility.interface_facet_tag (an integer mesh "
            "facet tag identifying the Si/SiO2 interface)"
        )


def _validate_schottky_contacts(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for `contacts[].type='schottky'`.

    The JSON schema declares `barrier_height_eV` as `['number', 'null']`
    so v2.0.0 through v2.4.0 inputs (no Schottky contact) keep
    validating without supplying the field. The schema cannot express
    "required only when type=='schottky'" without conditional logic;
    we enforce that here so users see the failure at `validate(cfg)`
    time rather than deep inside the FEM path. M16.5.
    """
    for contact in cfg.get("contacts", []):
        if contact.get("type") != "schottky":
            continue
        phi_b = contact.get("barrier_height_eV")
        name = contact.get("name", "<unnamed>")
        if phi_b is None:
            raise SchemaError(
                f"contact {name!r} has type='schottky' but no "
                f"barrier_height_eV; the metal-semiconductor barrier "
                f"is required for the thermionic-emission Robin BC "
                f"(M16.5; see Sze 3rd ed Table 5 for material-specific "
                f"barrier heights)"
            )
        if float(phi_b) < 0.0:
            raise SchemaError(
                f"contact {name!r}: barrier_height_eV must be "
                f"non-negative; got {phi_b}"
            )


def _validate_voltage_t(cfg: dict[str, Any]) -> None:
    """
    Cross-field validation for `contacts[].voltage_t` (M16.7).

    The JSON schema declares the variant shape and the enum on
    `voltage_t.type`. Constraints that the schema cannot express:

      - `voltage_t` and `voltage_sweep` on the same contact are
        mutually exclusive.
      - For `type='table'`: `times` and `values` must have equal
        length and `times` must be strictly monotonically increasing.
      - For `type='step'`: `t0`, `v0`, `v1` are all required.
      - `voltage_t` is consumed only by the transient runner; setting
        it on a contact when `solver.type` is anything other than
        `"transient"` raises so users see the failure at validate
        time rather than as a silent no-op.
    """
    solver_type = cfg.get("solver", {}).get("type", "equilibrium")
    for contact in cfg.get("contacts", []):
        vt = contact.get("voltage_t")
        if vt is None:
            continue
        name = contact.get("name", "<unnamed>")
        if "voltage_sweep" in contact:
            raise SchemaError(
                f"contact {name!r}: voltage_t and voltage_sweep are "
                f"mutually exclusive; pick one"
            )
        if solver_type != "transient":
            raise SchemaError(
                f"contact {name!r}: voltage_t is only consumed by the "
                f"transient runner (solver.type='transient'); got "
                f"solver.type={solver_type!r}. Use `voltage` for a "
                f"fixed bias or `voltage_sweep` for a continuation ramp "
                f"(M16.7)."
            )
        kind = vt.get("type")
        if kind == "table":
            times = vt.get("times")
            values = vt.get("values")
            if times is None or values is None:
                raise SchemaError(
                    f"contact {name!r}: voltage_t type='table' requires "
                    f"both `times` and `values`"
                )
            if len(times) != len(values):
                raise SchemaError(
                    f"contact {name!r}: voltage_t.times (len {len(times)}) "
                    f"and voltage_t.values (len {len(values)}) must have "
                    f"equal length"
                )
            if len(times) < 2:
                raise SchemaError(
                    f"contact {name!r}: voltage_t.times must have at "
                    f"least 2 entries; got {len(times)}"
                )
            for i in range(1, len(times)):
                if not (float(times[i]) > float(times[i - 1])):
                    raise SchemaError(
                        f"contact {name!r}: voltage_t.times must be "
                        f"strictly monotonically increasing; "
                        f"times[{i - 1}]={times[i - 1]} >= "
                        f"times[{i}]={times[i]}"
                    )
        elif kind == "step":
            for key in ("t0", "v0", "v1"):
                if key not in vt:
                    raise SchemaError(
                        f"contact {name!r}: voltage_t type='step' "
                        f"requires `{key}`"
                    )
        else:
            raise SchemaError(
                f"contact {name!r}: voltage_t.type must be 'table' or "
                f"'step'; got {kind!r}"
            )


_TUNNELING_DEFAULTS: dict[str, Any] = {
    "bbt": False,
    "tat": False,
    "A_kane": 4.0e14,
    "B_kane": 1.9e7,
    "tau_n_min": 1.0e-9,
    "tau_p_min": 1.0e-9,
    "F_kT": 1.4e7,
    "alpha": 2.0,
}


def _warn_bbt_with_boltzmann(cfg: dict[str, Any]) -> None:
    """
    Emit a UserWarning when physics.tunneling.bbt is enabled but the
    statistics dispatch is the Boltzmann default. The Kane prefactor
    depends on the FD-corrected density of states at the band edges, so
    the Boltzmann path uses a leading-order correction. M16.6.
    """
    phys = cfg.get("physics", {})
    tun = phys.get("tunneling", {}) or {}
    if not tun.get("bbt", False):
        return
    if phys.get("statistics", "boltzmann") == "boltzmann":
        warnings.warn(
            "physics.tunneling.bbt is enabled but physics.statistics is "
            "'boltzmann'; Kane band-to-band tunneling is most accurate "
            "under Fermi-Dirac statistics at heavy doping (>~1e19 cm^-3). "
            "The Boltzmann branch uses a leading-order correction. Set "
            "physics.statistics='fermi_dirac' for quantitative breakdown "
            "current (M16.6; see docs/PHYSICS.md section 1.4).",
            UserWarning,
            stacklevel=3,
        )


_LOMBARDI_DEFAULTS: dict[str, float] = {
    "B_n": 4.75e7,
    "B_p": 9.93e6,
    "C_n": 1.74e5,
    "C_p": 8.84e5,
    "lambda_n": 0.125,
    "lambda_p": 0.0317,
    "delta_n": 5.82e14,
    "delta_p": 2.0546e14,
}


def _fill_defaults(cfg: dict[str, Any]) -> dict[str, Any]:
    """Fill in sensible defaults for optional sections."""
    cfg.setdefault("coordinate_system", "cartesian")
    _validate_coordinate_system(cfg)
    phys = cfg.setdefault("physics", {})
    phys.setdefault("temperature", 300.0)
    rec = phys.setdefault("recombination", {})
    rec.setdefault("srh", True)
    rec.setdefault("tau_n", 1.0e-7)
    rec.setdefault("tau_p", 1.0e-7)
    rec.setdefault("E_t", 0.0)
    rec.setdefault("auger", False)
    # M16.3 Si Dziewior-Schmid Auger coefficients in cm^6/s. The
    # defaults are filled even when auger is false; the form builder
    # only consumes them when the flag is true, so v0.18.0 byte-
    # identity is preserved on the auger-off branch.
    rec.setdefault("C_n", 2.8e-31)
    rec.setdefault("C_p", 9.9e-32)
    # M16.6 tunneling defaults. Both flags default to false; the
    # parameters are filled to Si textbook defaults but unused when the
    # flags are off, preserving v0.21.0 byte-identity on every existing
    # benchmark.
    tun = phys.setdefault("tunneling", {})
    for key, default in _TUNNELING_DEFAULTS.items():
        tun.setdefault(key, default)
    mob = phys.setdefault("mobility", {})
    mob.setdefault("model", "constant")
    mob.setdefault("mu_n", 1400.0)
    mob.setdefault("mu_p", 450.0)
    if mob.get("model") == "lombardi":
        mob.setdefault("bulk_model", "constant")
        lombardi = mob.setdefault("lombardi", {})
        for key, default in _LOMBARDI_DEFAULTS.items():
            lombardi.setdefault(key, default)
    _validate_mobility_lombardi(cfg)
    phys.setdefault("statistics", "boltzmann")
    _warn_bbt_with_boltzmann(cfg)

    # M16.5: every contact carries an explicit barrier_height_eV slot
    # (default null). Schottky contacts then have a non-null value
    # validated by _validate_schottky_contacts; non-Schottky contacts
    # ignore it.
    for contact in cfg.get("contacts", []):
        contact.setdefault("barrier_height_eV", None)
    _validate_schottky_contacts(cfg)

    # M17: every region carries an explicit heterojunction flag (default
    # false). Regions without material_overrides leave that field absent
    # so consumers (semi/physics/heterojunction.py) can distinguish "no
    # override" from "override set to default-shaped values"; the field
    # is optional with no default.
    regions = cfg.get("regions", {})
    if isinstance(regions, dict):
        for region_cfg in regions.values():
            if isinstance(region_cfg, dict):
                region_cfg.setdefault("heterojunction", False)

    solver = cfg.setdefault("solver", {})
    solver.setdefault("type", "equilibrium")
    solver.setdefault("backend", "cpu-mumps")
    solver.setdefault("max_iterations", 50)
    solver.setdefault("atol", 1.0e-10)
    solver.setdefault("rtol", 1.0e-8)
    solver.setdefault("damping", 1.0)
    ls = solver.setdefault("linear_solver", {})
    ls.setdefault("ksp_type", "preonly")
    ls.setdefault("pc_type", "lu")
    ls.setdefault("factor_mat_solver_type", "mumps")
    cont = solver.setdefault("continuation", {})
    cont.setdefault("steps", 10)
    cont.setdefault("adaptive", True)
    cont.setdefault("min_step", 1.0e-4)
    cont.setdefault("max_halvings", 6)
    cont.setdefault("easy_iter_threshold", 4)
    cont.setdefault("grow_factor", 1.5)

    _validate_compute(cfg)
    _validate_voltage_t(cfg)

    out = cfg.setdefault("output", {})
    out.setdefault("directory", "./results")
    out.setdefault("fields", ["potential", "n", "p"])
    out.setdefault("write_xdmf", True)
    out.setdefault("write_iv", True)

    return cfg


def load(path: str | Path) -> dict[str, Any]:
    """Load and validate a JSON config file; record its source path for relative lookups."""
    path = Path(path)
    with path.open() as f:
        cfg = json.load(f)
    cfg = validate(cfg)
    cfg["_source_path"] = str(path.resolve())
    cfg["_source_dir"] = str(path.parent.resolve())
    return cfg


def dumps(cfg: dict[str, Any], indent: int = 2) -> str:
    """Serialize a config back to JSON, stripping internal bookkeeping fields."""
    clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    return json.dumps(clean, indent=indent)
