import json
from pathlib import Path

import pytest

from models_v15 import Fixed60HighResolutionReference


ROOT = Path(__file__).resolve().parents[1]


def test_v15_delivered_mesh_and_scope_guard():
    plant = Fixed60HighResolutionReference()
    assert plant.state_count == 4800
    plant.assert_in_scope(ambient_C=23.0, duty=0.60, current_A=30.0)
    with pytest.raises(ValueError):
        plant.assert_in_scope(ambient_C=23.0, duty=0.80, current_A=30.0)


def test_v15_comsol_consistency_passes_without_overriding_experimental_gate():
    result = json.loads((ROOT / "outputs_v15/validation_metrics.json").read_text())
    assert result["COMSOL"]["implementation_consistency_passed"] is True
    assert result["decisions"]["strict_experimental_high_fidelity_gate_passed"] is False
    assert result["decisions"]["ROM_MM_EKF_oracle_MPC_authorized"] is False


def test_v15_55C_is_reference_only():
    result = json.loads((ROOT / "outputs_v15/final_summary.json").read_text())
    assert result["reference_temperature_C"] == 55.0
    assert result["reference_temperature_is_not_a_safety_limit"] is True
