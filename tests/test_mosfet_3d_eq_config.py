"""Phase D test: the mosfet_3d_eq.json config validates as v2.11.0."""
from __future__ import annotations

import json
from pathlib import Path

from semi.schema import validate

CFG_PATH = (Path(__file__).parent.parent
            / "benchmarks" / "mosfet_3d_eq" / "mosfet_3d_eq.json")


def test_config_exists():
    assert CFG_PATH.exists(), f"missing {CFG_PATH}"


def test_config_validates_as_v211():
    with CFG_PATH.open() as f:
        cfg = json.load(f)
    out = validate(cfg)
    assert out["mesh"]["quality_gate"] is True
    assert out["solver"]["type"] == "equilibrium"
    assert out["dimension"] == 3


def test_config_contact_voltages_all_zero():
    with CFG_PATH.open() as f:
        cfg = json.load(f)
    for contact in cfg["contacts"]:
        assert contact.get("voltage", 0.0) == 0.0, (
            f"contact {contact['name']} must be V=0 at equilibrium "
            f"(M19 precursor; ADR 0019 anti-goal: no bias sweep)"
        )
        # No sweep block on any contact.
        assert "voltage_sweep" not in contact


def test_config_implant_centres_inside_mesh():
    """The Gaussian implant centres must lie inside the mesh bbox."""
    with CFG_PATH.open() as f:
        cfg = json.load(f)
    for entry in cfg["doping"]:
        profile = entry["profile"]
        if profile["type"] != "gaussian":
            continue
        cx, cy, cz = profile["center"]
        assert 0.0 <= cx <= 1.5e-6
        assert 0.0 <= cy <= 5.0e-7
        assert 0.0 <= cz <= 2.005e-6


def test_config_regions_match_geo_tags():
    """Region tags 1 (si_body) and 4 (sio2_gate) must match the
    .geo's Physical Volume tags."""
    with CFG_PATH.open() as f:
        cfg = json.load(f)
    assert cfg["regions"]["si_body"]["tag"] == 1
    assert cfg["regions"]["sio2_gate"]["tag"] == 4


def test_config_contact_facet_tags_match_geo():
    """Contact facet tags 10/11/12/13 must match the .geo's
    Physical Surface tags."""
    with CFG_PATH.open() as f:
        cfg = json.load(f)
    by_name = {c["name"]: c["facet"] for c in cfg["contacts"]}
    assert by_name["source"] == 10
    assert by_name["drain"]  == 11
    assert by_name["gate"]   == 12
    assert by_name["body"]   == 13
