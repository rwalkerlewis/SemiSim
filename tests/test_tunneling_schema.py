"""
Tests for the M16.6 schema additions: physics.tunneling sub-object
with bbt and tat boolean flags plus the Kane (A_kane, B_kane) and
Hurkx (tau_n_min, tau_p_min, F_kT, alpha) parameters.

Pure-Python; no dolfinx required. Mirrors the M16.5 schema-test
pattern: every existing benchmark JSON keeps validating, the new
fields default to false / textbook Si values, and the conditional
warning fires at validate time when the BBT path runs under
Boltzmann statistics.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

from semi import schema


@pytest.fixture
def minimal_cfg():
    return {
        "schema_version": "2.6.0",
        "name": "tunneling_test",
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
                "profile": {"type": "uniform", "N_D": 1.0e19, "N_A": 0.0},
            }
        ],
        "contacts": [
            {"name": "anode", "facet": "anode", "type": "ohmic",
             "voltage": 0.0},
            {"name": "cathode", "facet": "cathode", "type": "ohmic",
             "voltage": 0.0},
        ],
    }


def test_default_fill_leaves_tunneling_flags_off(minimal_cfg):
    result = schema.validate(minimal_cfg)
    tun = result["physics"]["tunneling"]
    assert tun["bbt"] is False
    assert tun["tat"] is False
    assert tun["A_kane"] == pytest.approx(4.0e14)
    assert tun["B_kane"] == pytest.approx(1.9e7)
    assert tun["tau_n_min"] == pytest.approx(1.0e-9)
    assert tun["tau_p_min"] == pytest.approx(1.0e-9)
    assert tun["F_kT"] == pytest.approx(1.4e7)
    assert tun["alpha"] == pytest.approx(2.0)


def test_bbt_with_fermi_dirac_validates_silently(minimal_cfg):
    minimal_cfg["physics"] = {
        "tunneling": {"bbt": True},
        "statistics": "fermi_dirac",
    }
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        result = schema.validate(minimal_cfg)
    assert result["physics"]["tunneling"]["bbt"] is True
    assert result["physics"]["statistics"] == "fermi_dirac"


def test_bbt_with_boltzmann_warns(minimal_cfg):
    minimal_cfg["physics"] = {
        "tunneling": {"bbt": True},
        "statistics": "boltzmann",
    }
    with pytest.warns(UserWarning, match="fermi_dirac"):
        schema.validate(minimal_cfg)


def test_bbt_off_under_boltzmann_does_not_warn(minimal_cfg):
    minimal_cfg["physics"] = {
        "tunneling": {"bbt": False, "tat": True},
        "statistics": "boltzmann",
    }
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        schema.validate(minimal_cfg)


def test_tat_only_validates(minimal_cfg):
    minimal_cfg["physics"] = {"tunneling": {"tat": True}}
    result = schema.validate(minimal_cfg)
    assert result["physics"]["tunneling"]["tat"] is True
    assert result["physics"]["tunneling"]["bbt"] is False


def test_negative_A_kane_rejected(minimal_cfg):
    minimal_cfg["physics"] = {"tunneling": {"bbt": True, "A_kane": -1.0}}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_negative_B_kane_rejected(minimal_cfg):
    minimal_cfg["physics"] = {"tunneling": {"bbt": True, "B_kane": -1.0}}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_negative_F_kT_rejected(minimal_cfg):
    minimal_cfg["physics"] = {"tunneling": {"tat": True, "F_kT": -1.0}}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_unknown_tunneling_key_rejected(minimal_cfg):
    minimal_cfg["physics"] = {"tunneling": {"banana": 1.0}}
    with pytest.raises(schema.SchemaError):
        schema.validate(minimal_cfg)


def test_existing_benchmarks_validate_under_v2_6_0():
    """Every shipped benchmark JSON must continue to validate against
    the v2.6.0 strict schema."""
    benchmarks_dir = Path(__file__).parent.parent / "benchmarks"
    for json_path in sorted(benchmarks_dir.rglob("*.json")):
        if json_path.parent.name == "results":
            continue
        with json_path.open() as f:
            cfg = json.load(f)
        if not str(cfg.get("schema_version", "")).startswith("2."):
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            schema.validate(cfg)


def test_supported_minor_advanced_to_6():
    assert schema.SCHEMA_SUPPORTED_MINOR >= 6


def test_v2_5_0_input_still_validates(minimal_cfg):
    """A pre-M16.6 input (no tunneling block) must still validate."""
    minimal_cfg["schema_version"] = "2.5.0"
    result = schema.validate(minimal_cfg)
    # Default-fill writes the tunneling defaults so downstream consumers
    # always see a populated block.
    assert result["physics"]["tunneling"]["bbt"] is False
