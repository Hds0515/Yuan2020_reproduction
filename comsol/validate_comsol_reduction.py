"""Compare converged COMSOL regional temperatures with the frozen 3-node steady model.

The script refuses to invent data: it exits with status 2 until the COMSOL CSV
exists and contains finite results for every predeclared inlet speed.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from models.three_node_model import load_frozen_parameters  # noqa: E402


def three_node_steady_temperature(heat_total_W: float, duty: float, direction: int = 1) -> np.ndarray:
    parameters = load_frozen_parameters(ROOT / "identification_v3" / "frozen_parameters_v3.json")
    coolant = np.array([24.8, 30.05, 35.3])
    if direction == -1:
        coolant = coolant[::-1]
    # At the present operating points all nodes remain warmer than coolant,
    # so the steady equations are linear and can be solved exactly.
    kc = parameters.K_cool_W_K * duty
    kn = parameters.K_node_W_K
    matrix = np.array([[kc + kn, -kn, 0.0], [-kn, kc + 2 * kn, -kn], [0.0, -kn, kc + kn]])
    rhs = np.full(3, heat_total_W / 3.0) + kc * coolant
    return np.linalg.solve(matrix, rhs)


def main() -> int:
    input_path = HERE / "comsol_zone_temperatures.csv"
    status_path = HERE / "status.json"
    if not input_path.exists():
        print(json.dumps({"validated": False, "reason": "COMSOL result CSV is absent", "status": str(status_path)}, indent=2))
        return 2
    data = pd.read_csv(input_path)
    required = ["uin_m_s", "T_zone1_C", "T_zone2_C", "T_zone3_C", "Tmax_C"]
    missing = [column for column in required if column not in data]
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if len(data) != 6 or not np.isfinite(data[required].to_numpy()).all():
        raise ValueError("Expected six finite predeclared operating points")
    data["normalized_fan_command"] = np.clip(data["uin_m_s"] / 12.0, 0.0, 1.0)
    # Java heat flux times the explicitly declared stack footprint.
    heat_total_W = 2425.5 * 0.30 * 0.082333333
    reduced = np.vstack([three_node_steady_temperature(heat_total_W, duty) for duty in data.normalized_fan_command])
    for i in range(3):
        data[f"three_node_T{i + 1}_C"] = reduced[:, i]
        data[f"error_T{i + 1}_C"] = reduced[:, i] - data[f"T_zone{i + 1}_C"]
    error = reduced - data[["T_zone1_C", "T_zone2_C", "T_zone3_C"]].to_numpy()
    summary = {
        "validated": True,
        "number_of_operating_points": int(len(data)),
        "regional_temperature_rmse_C": float(np.sqrt(np.mean(error**2))),
        "maximum_regional_absolute_error_C": float(np.max(np.abs(error))),
        "mapping_assumption": "normalized fan command = inlet speed / 12 m/s",
    }
    data.to_csv(HERE / "results" / "comsol_three_node_comparison.csv", index=False)
    (ROOT / "outputs" / "comsol_reduction_validation.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
