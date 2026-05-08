"""
Tests for the M17 schema additions: regions[*].material_overrides and
regions[*].heterojunction. Pure-Python; no dolfinx required.

Mirrors the M16.5 / M16.6 / M16.7 schema-test pattern: every existing
benchmark JSON keeps validating, the new fields are optional, the
default-fill leaves material_overrides absent and heterojunction:false,
and a config without either field is bit-identical to v0.23.0.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def two_region_cfg():
    """Minimal 1D two-region cfg for exercising the M17 schema additions.

    Si on the left half, GaAs on the right. No heterojunction flag and
    no material_overrides set; this is the v0.23.0-shape baseline.
    """
    return {
        "schema_version": "2.8.0",
        "name": "heterojunction_schema_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [10],
            "regions_by_box": [
                {"name": "left",  "tag": 1, "bounds": [[0.0, 5.0e-7]]},
                {"name": "right", "tag": 2, "bounds": [[5.0e-7, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {
            "left":  {"material": "Si",   "tag": 1, "role": "semiconductor"},
            "right": {"material": "GaAs", "tag": 2, "role": "semiconductor"},
        },
        "doping": [
            {"region": "left",  "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0}},
            {"region": "right", "profile": {"type": "uniform", "N_D": 1.0e17, "N_A": 0.0}},
        ],
        "contacts": [
            {"name": "anode",   "facet": "anode",   "type": "ohmic", "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic", "voltage": 0.0},
        ],
        "solver": {"type": "equilibrium"},
    }


def test_supported_minor_advanced_to_8():
    assert schema.SCHEMA_SUPPORTED_MINOR >= 8


def test_baseline_two_region_cfg_validates(two_region_cfg):
    """No material_overrides, no heterojunction: this is the v0.23.0
    shape and must validate cleanly."""
    result = schema.validate(two_region_cfg)
    assert result["regions"]["left"]["material"] == "Si"
    assert result["regions"]["right"]["material"] == "GaAs"


def test_default_fill_sets_heterojunction_false(two_region_cfg):
    """The default-fill must mark heterojunction:false on every region
    so downstream consumers do not have to .get(default) explicitly."""
    result = schema.validate(two_region_cfg)
    assert result["regions"]["left"]["heterojunction"] is False
    assert result["regions"]["right"]["heterojunction"] is False


def test_default_fill_leaves_material_overrides_absent(two_region_cfg):
    """material_overrides is optional and intentionally not defaulted:
    consumers distinguish 'not set' (use material as-is) from 'present
    with explicit values' (apply the overrides). v0.23.0 byte-identity
    is preserved when the field is absent."""
    result = schema.validate(two_region_cfg)
    assert "material_overrides" not in result["regions"]["left"]
    assert "material_overrides" not in result["regions"]["right"]


def test_material_overrides_chi_eg_validates(two_region_cfg):
    two_region_cfg["regions"]["right"]["material_overrides"] = {
        "chi_eV": 4.07,
        "Eg_eV": 1.42,
    }
    result = schema.validate(two_region_cfg)
    overrides = result["regions"]["right"]["material_overrides"]
    assert overrides["chi_eV"] == pytest.approx(4.07)
    assert overrides["Eg_eV"] == pytest.approx(1.42)


def test_material_overrides_full_set_validates(two_region_cfg):
    """All four override fields together (chi, Eg, Nc, Nv)."""
    two_region_cfg["regions"]["right"]["material_overrides"] = {
        "chi_eV": 3.74,
        "Eg_eV": 1.798,
        "Nc_per_cm3": 5.5e17,
        "Nv_per_cm3": 1.0e19,
    }
    result = schema.validate(two_region_cfg)
    overrides = result["regions"]["right"]["material_overrides"]
    assert overrides["chi_eV"] == pytest.approx(3.74)
    assert overrides["Eg_eV"] == pytest.approx(1.798)
    assert overrides["Nc_per_cm3"] == pytest.approx(5.5e17)
    assert overrides["Nv_per_cm3"] == pytest.approx(1.0e19)


def test_material_overrides_partial_subset_validates(two_region_cfg):
    """Setting only chi_eV leaves Eg/Nc/Nv to fall back to material."""
    two_region_cfg["regions"]["right"]["material_overrides"] = {
        "chi_eV": 3.5,
    }
    result = schema.validate(two_region_cfg)
    overrides = result["regions"]["right"]["material_overrides"]
    assert overrides == {"chi_eV": 3.5}


def test_heterojunction_true_validates(two_region_cfg):
    two_region_cfg["regions"]["right"]["heterojunction"] = True
    result = schema.validate(two_region_cfg)
    assert result["regions"]["right"]["heterojunction"] is True
    assert result["regions"]["left"]["heterojunction"] is False


def test_negative_eg_rejected(two_region_cfg):
    two_region_cfg["regions"]["right"]["material_overrides"] = {"Eg_eV": -1.0}
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_zero_eg_rejected(two_region_cfg):
    two_region_cfg["regions"]["right"]["material_overrides"] = {"Eg_eV": 0.0}
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_negative_nc_rejected(two_region_cfg):
    two_region_cfg["regions"]["right"]["material_overrides"] = {"Nc_per_cm3": -1.0}
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_negative_nv_rejected(two_region_cfg):
    two_region_cfg["regions"]["right"]["material_overrides"] = {"Nv_per_cm3": -1.0}
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_unknown_material_overrides_field_rejected(two_region_cfg):
    """additionalProperties:false on material_overrides means typos
    fail validation rather than being silently dropped."""
    two_region_cfg["regions"]["right"]["material_overrides"] = {
        "chi_eV": 4.0,
        "BogusField": 1.0,
    }
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_unknown_region_field_still_rejected(two_region_cfg):
    """The existing additionalProperties:false on regions[*] catches
    typos; the M17 additions must not relax that gate."""
    two_region_cfg["regions"]["left"]["typo"] = 1.0
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_heterojunction_non_boolean_rejected(two_region_cfg):
    two_region_cfg["regions"]["right"]["heterojunction"] = "yes"
    with pytest.raises(schema.SchemaError):
        schema.validate(two_region_cfg)


def test_existing_benchmarks_validate_under_v2_8_0():
    """Every shipped benchmark JSON must continue to validate against
    the v2.8.0 strict schema (no benchmark uses material_overrides or
    heterojunction today, so every prior config is bit-identical to
    v0.23.0)."""
    benchmarks_dir = Path(__file__).parent.parent / "benchmarks"
    for json_path in sorted(benchmarks_dir.rglob("*.json")):
        if json_path.parent.name == "results":
            continue
        with json_path.open() as f:
            cfg = json.load(f)
        if not str(cfg.get("schema_version", "")).startswith("2."):
            continue
        schema.validate(copy.deepcopy(cfg))


def test_existing_examples_validate_under_v2_8_0():
    """The PR #85 examples directory must also continue to validate."""
    examples_dir = Path(__file__).parent.parent / "examples"
    if not examples_dir.exists():
        pytest.skip("examples/ directory not present")
    for json_path in sorted(examples_dir.rglob("*.json")):
        if json_path.parent.name == "results":
            continue
        with json_path.open() as f:
            cfg = json.load(f)
        if not str(cfg.get("schema_version", "")).startswith("2."):
            continue
        schema.validate(copy.deepcopy(cfg))


def test_v2_7_0_input_still_validates(two_region_cfg):
    """Backward-compatibility: a v2.7.0-shaped input (no heterojunction
    additions) keeps validating against the v2.8.0 schema."""
    two_region_cfg["schema_version"] = "2.7.0"
    schema.validate(two_region_cfg)


def test_v2_0_0_input_still_validates(two_region_cfg):
    """Backward-compatibility: a v2.0.0-shaped input keeps validating."""
    two_region_cfg["schema_version"] = "2.0.0"
    schema.validate(two_region_cfg)
