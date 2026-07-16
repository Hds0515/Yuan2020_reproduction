import json
from pathlib import Path

import numpy as np

from comsol.v8_runtime.cooling_structure_models import PHYSICS, solve_structure


ROOT = Path(__file__).resolve().parents[1]


def test_v8_mass_and_energy_identities() -> None:
    result = solve_structure("M2", 8.0, np.array([91.0974242118531, 0.6381724843102996]), nodes=40)
    assert result["mdot_identity_relative_error"] < 1e-12
    assert result["energy_residual_relative"] < 0.005
    expected = PHYSICS.air_rho_kg_m3 * 8.0 * PHYSICS.stack_inlet_area_m2
    assert abs(result["mdot_cell_kg_s"] - expected) < 1e-12


def test_v8_digitization_gate_preserves_threshold() -> None:
    audit = json.loads(
        (ROOT / "outputs_v8/digitization_reaudit_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert audit["digitization_reliable"] is True
    assert audit["redigitization_triggered"] is False
    assert audit["validation_threshold_changed"] is False
    assert audit["overall_manual_vs_automatic_region_RMSE_K"] <= 1.0


def test_v8_stop_rule_blocks_downstream_work() -> None:
    summary = json.loads(
        (ROOT / "outputs_v8/final_summary.json").read_text(encoding="utf-8")
    )
    decision = summary["model_decision"]
    assert decision["M1_passed"] is False
    assert decision["M2_passed"] is False
    assert decision["stop_rule_triggered"] is True
    assert decision["downstream_ROM_or_Oracle_authorized"] is False
    assert decision["selected_model_sha256"] is None
