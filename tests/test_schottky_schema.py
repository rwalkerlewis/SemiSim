"""
Tests for the M16.5 schema additions: contacts[].type='schottky' and
contacts[].barrier_height_eV. Pure-Python; no dolfinx required.

Mirrors the M16.2 / M16.3 / M16.4 schema-test pattern: every existing
benchmark JSON keeps validating, the new field is optional except on
Schottky contacts where its non-null non-negative numeric value is
enforced by `_validate_schottky_contacts`.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def minimal_cfg():
    return {
        "schema_version": "2.5.0",
        "name": "schottky_test",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [100],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "anode",   "tag": 1, "axis": 0, "value": 0.0},
                {"name": "cathode", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {
            "silicon": {"material": "Si", "tag": 1, "role": "semiconductor"},
        },
        "doping": [
            {
                "region": "silicon",
                "profile": {"type": "uniform", "N_D": 1.0e16, "N_A": 0.0},
            }
        ],
        "contacts": [
            {
                "name": "anode", "facet": "anode", "type": "schottky",
                "barrier_height_eV": 0.85, "voltage": 0.0,
            },
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
    }


def test_schottky_contact_with_barrier_validates(minimal_cfg):
    result = schema.validate(minimal_cfg)
    anode = result["contacts"][0]
    assert anode["type"] == "schottky"
    assert anode["barrier_height_eV"] == pytest.approx(0.85)


def test_schottky_contact_without_barrier_rejected(minimal_cfg):
    minimal_cfg["contacts"][0]["barrier_height_eV"] = None
    with pytest.raises(schema.SchemaError, match="barrier_height_eV"):
        schema.validate(minimal_cfg)


def test_schottky_contact_missing_barrier_field_rejected(minimal_cfg):
    del minimal_cfg["contacts"][0]["barrier_height_eV"]
    with pytest.raises(schema.SchemaError, match="barrier_height_eV"):
        schema.validate(minimal_cfg)


def test_schottky_negative_barrier_rejected(minimal_cfg):
    minimal_cfg["contacts"][0]["barrier_height_eV"] = -0.1
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_unknown_contact_type_still_rejected(minimal_cfg):
    minimal_cfg["contacts"][0]["type"] = "magic_unicorn"
    minimal_cfg["contacts"][0].pop("barrier_height_eV", None)
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_non_schottky_contact_ignores_barrier(minimal_cfg):
    """Ohmic contact may carry a barrier_height_eV value but the loader
    leaves it untouched and does not reject it (default-fill is null)."""
    minimal_cfg["contacts"][0]["type"] = "ohmic"
    minimal_cfg["contacts"][0]["barrier_height_eV"] = None
    result = schema.validate(minimal_cfg)
    assert result["contacts"][0]["type"] == "ohmic"


def test_default_fill_leaves_barrier_null_when_unset():
    """Non-Schottky contacts that omit barrier_height_eV get a null
    default; existing benchmark JSONs that predate M16.5 are bit-
    identical to v0.20.0 after the additive bump."""
    cfg = {
        "schema_version": "2.4.0",
        "name": "ohmic_only",
        "dimension": 1,
        "mesh": {
            "source": "builtin",
            "extents": [[0.0, 1.0e-6]],
            "resolution": [10],
            "regions_by_box": [
                {"name": "silicon", "tag": 1, "bounds": [[0.0, 1.0e-6]]},
            ],
            "facets_by_plane": [
                {"name": "L", "tag": 1, "axis": 0, "value": 0.0},
                {"name": "R", "tag": 2, "axis": 0, "value": 1.0e-6},
            ],
        },
        "regions": {"silicon": {"material": "Si", "tag": 1,
                                "role": "semiconductor"}},
        "doping": [
            {"region": "silicon",
             "profile": {"type": "uniform", "N_D": 1e17, "N_A": 0.0}}
        ],
        "contacts": [
            {"name": "L", "facet": "L", "type": "ohmic", "voltage": 0.0},
            {"name": "R", "facet": "R", "type": "ohmic", "voltage": 0.0},
        ],
    }
    result = schema.validate(cfg)
    for contact in result["contacts"]:
        assert contact["barrier_height_eV"] is None


def test_existing_benchmarks_validate_under_v2_5_0():
    """Every shipped benchmark JSON must continue to validate against
    the v2.5.0 strict schema."""
    benchmarks_dir = Path(__file__).parent.parent / "benchmarks"
    for json_path in sorted(benchmarks_dir.rglob("*.json")):
        if json_path.parent.name == "results":
            continue
        with json_path.open() as f:
            cfg = json.load(f)
        if not str(cfg.get("schema_version", "")).startswith("2."):
            continue
        # Each benchmark may rely on relative paths; we only need to
        # confirm the JSON itself validates against the current schema.
        schema.validate(cfg)


def test_supported_minor_advanced_to_5():
    assert schema.SCHEMA_SUPPORTED_MINOR >= 5
