from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from models_v13.fixed_pwm_plant import FixedPWM60Plant
from models_v15.fixed60_high_resolution_reference import FROZEN_DYNAMIC
from models_v16 import (
    FiveZoneParameters,
    PhysicsConstrainedFiveZoneROM,
    ProjectedTwoMassFiveZoneROM,
    conservative_five_zone_average,
)
from models_v16.physics_constrained_five_zone_rom import TOTAL_HEAT_CAPACITY_J_K


ROOT = Path(__file__).resolve().parents[1]


def test_conservative_five_zone_remap_preserves_mean() -> None:
    field = np.arange(2 * 40 * 6 * 10, dtype=float).reshape(2, 40, 6, 10)
    zones = conservative_five_zone_average(field)
    assert zones.shape == (2, 5)
    np.testing.assert_allclose(zones.mean(axis=1), field.mean(axis=(1, 2, 3)))


def test_single_mass_rom_is_stable_positive_and_energy_conserving() -> None:
    parameters = FiveZoneParameters(
        np.full(5, TOTAL_HEAT_CAPACITY_J_K / 5),
        np.full(4, 5.0),
        np.full(5, 10.0),
        np.zeros(2),
        np.zeros(2),
    )
    model = PhysicsConstrainedFiveZoneROM(parameters)
    model.assert_physical()
    audit = model.instantaneous_energy_audit(
        np.asarray((30.0, 31.0, 32.0, 33.0, 34.0)), 20.0, 27.0
    )
    assert audit["energy_residual_percent"] < 1e-10


def test_projected_two_mass_rom_preserves_parent_physics() -> None:
    parent = FixedPWM60Plant(FROZEN_DYNAMIC, flow_rows=15, columns=5)
    model = ProjectedTwoMassFiveZoneROM.project_from_v15(parent)
    model.assert_physical()
    assert model.state_count == 10
    assert np.isclose(model.parameters.capacity_J_K.sum(), TOTAL_HEAT_CAPACITY_J_K)
    audit = model.instantaneous_energy_audit(np.full(10, 35.0), 20.0, 27.0)
    assert audit["energy_residual_percent"] < 1e-10


def test_two_mass_freeze_uses_current_protocol_hash() -> None:
    protocol = ROOT / "outputs_v16/preregistered_two_mass_rom_protocol.json"
    frozen = json.loads(
        (ROOT / "outputs_v16/frozen_projected_two_mass_rom.json").read_text(
            encoding="utf-8"
        )
    )
    digest = hashlib.sha256(protocol.read_bytes()).hexdigest()
    assert frozen["protocol_sha256"] == digest


def test_v16_does_not_authorize_control_or_formal_rom_claim() -> None:
    summary = json.loads(
        (ROOT / "outputs_v16/final_summary.json").read_text(encoding="utf-8")
    )
    assert not summary["step_4_high_resolution_confirmation"][
        "strict_experimental_gate_passed"
    ]
    assert not summary["step_5_five_zone_ROM"]["formal_claim_authorized"]
    assert not summary["step_7_control_decision"]["control_research_authorized"]
