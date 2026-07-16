import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_validation_boundary_and_metrics() -> None:
    split = json.loads((ROOT / "identification_v3" / "calibration_validation_split.json").read_text())
    assert split["calibration"][1]["time_s"] == [0, 700]
    assert split["independent_validation"][0]["time_s"] == [700, 1200]
    assert split["validation_used_in_parameter_optimization"] is False
    metrics = json.loads((ROOT / "identification_v3" / "identification_metrics.json").read_text())
    assert metrics["fig14_calibration_rmse_C"] <= 0.35
    assert metrics["fig16_independent_validation_rmse_C"] <= 0.40


def test_heat_model_was_compared_after_calibration() -> None:
    comparison = pd.read_csv(ROOT / "identification_v3" / "heat_model_comparison.csv")
    assert set(comparison["heat_model"]) == {"A_thermoneutral_proxy", "B_quadratic_current"}
    chosen = comparison.sort_values("fig16_independent_validation_rmse_C").iloc[0]
    assert chosen["heat_model"] == "A_thermoneutral_proxy"


def test_frozen_parameter_hash() -> None:
    path = ROOT / "identification_v3" / "frozen_parameters_v3.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    expected = payload.pop("sha256_without_hash_field")
    actual = hashlib.sha256(json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")).hexdigest()
    assert actual == expected
