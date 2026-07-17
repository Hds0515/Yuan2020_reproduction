"""Fit PWM40 fan exponent and confirm the frozen plant at PWM60/80."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models_v14.multi_pwm_plant import MultiPWMPlant  # noqa: E402


OUT = ROOT / "outputs_v14"
FROZEN = OUT / "frozen_multi_pwm_high_fidelity_model.json"
MARKER = OUT / "final_confirmation_executed.json"
LOWER, UPPER = 0.3, 2.5


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def target_for(pwm: int) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.read_csv(OUT / "digitization/guo2024_multi_pwm_average_temperature.csv")
    selected = frame.loc[frame.PWM_percent == pwm].sort_values("time_s")
    return selected.time_s.to_numpy(float), selected.average_temperature_C.to_numpy(float)


def prediction(exponent: float, pwm: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    times, target = target_for(pwm)
    event_times = np.arange(0.0, 1101.0, 100.0)
    event_current = np.minimum(4.0 * (np.arange(len(event_times)) + 1), 40.0)
    grid = np.unique(np.concatenate((event_times[event_times >= times[0]], times)))
    indices = np.maximum(np.searchsorted(event_times, grid, side="right") - 1, 0)
    current = event_current[indices]
    voltage = np.clip(33.5 - 0.34 * current, 19.5, 33.5)
    plant = MultiPWMPlant(exponent)
    state = plant.simulate_from_regional_initial_state(
        grid,
        current,
        voltage,
        np.repeat(target[0], 3),
        ambient_C=23.0,
        duty=pwm / 100.0,
    )
    sensors = plant.sensor_field(plant.surface_from_state(state))
    predicted = sensors[np.searchsorted(grid, times)].mean(axis=(1, 2, 3))
    return times, target, predicted


def metrics(block: tuple[np.ndarray, np.ndarray, np.ndarray]) -> dict[str, float]:
    times, target, predicted = block
    error = predicted[1:] - target[1:]
    return {
        "RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
        "maximum_absolute_error_C": float(np.max(np.abs(error))),
        "mean_bias_C": float(np.mean(error)),
        "temperature_trend_spearman": float(spearmanr(target[1:], predicted[1:]).statistic),
        "scored_start_s": float(times[1]),
        "scored_end_s": float(times[-1]),
    }


def main() -> None:
    if FROZEN.exists():
        payload = json.loads(FROZEN.read_text(encoding="utf-8"))
        exponent = float(payload["fan_flow_exponent_relative_to_60pct"])
        action = "loaded existing freeze; no optimization"
    else:
        _, target, _ = prediction(1.2, 40)

        def residual(value: np.ndarray) -> np.ndarray:
            return prediction(float(value[0]), 40)[2][1:] - target[1:]

        fit = least_squares(
            residual,
            np.asarray((1.2,)),
            bounds=(np.asarray((LOWER + 1e-6,)), np.asarray((UPPER - 1e-6,))),
            xtol=1e-9,
            ftol=1e-9,
            gtol=1e-9,
            max_nfev=40,
        )
        exponent = float(fit.x[0])
        payload = {
            "schema_version": 1,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "protocol_sha256": sha256(OUT / "preregistered_multi_pwm_protocol.json"),
            "model_source_sha256": sha256(ROOT / "models_v14/multi_pwm_plant.py"),
            "digitized_data_sha256": sha256(
                OUT / "digitization/guo2024_multi_pwm_average_temperature.csv"
            ),
            "fan_flow_exponent_relative_to_60pct": exponent,
            "plant_record": MultiPWMPlant(exponent, flow_rows=6, columns=10).record(),
            "calibration": "Guo2024 Fig.11 PWM40, 150-750 s",
            "validation": "Guo2024 Fig.11 PWM60, 150-950 s",
            "final_confirmation": "Guo2024 Fig.11 PWM80, 150-950 s",
            "no_retuning_after_freeze": True,
        }
        FROZEN.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        action = "PWM40 exponent optimized then frozen"

    frozen_hash = sha256(FROZEN)
    calibration = prediction(exponent, 40)
    calibration_metrics = metrics(calibration)
    validation = prediction(exponent, 60)
    validation_metrics = metrics(validation)
    first_confirmation = not MARKER.exists()
    confirmation = prediction(exponent, 80)
    confirmation_metrics = metrics(confirmation)

    medium = MultiPWMPlant(exponent, flow_rows=6, columns=10)
    fine = MultiPWMPlant(exponent, flow_rows=12, columns=20)
    voltage = float(np.clip(33.5 - 0.34 * 40.0, 19.5, 33.5))
    medium_surface, medium_core = medium.steady_state(
        40.0, voltage, duty=0.80, ambient_C=23.0
    )
    fine_surface, _ = fine.steady_state(40.0, voltage, duty=0.80, ambient_C=23.0)
    mesh_difference = float(
        np.max(np.abs(medium.sensor_field(medium_surface) - fine.sensor_field(fine_surface)))
    )
    energy = medium.energy_audit(
        medium_surface, medium_core, 40.0, voltage, 0.80, 23.0
    )

    protocol = json.loads((OUT / "preregistered_multi_pwm_protocol.json").read_text())
    limit = protocol["acceptance_thresholds"]
    checks = {
        "calibration_RMSE": calibration_metrics["RMSE_C"] <= limit["calibration_RMSE_C"],
        "validation_RMSE": validation_metrics["RMSE_C"] <= limit["validation_RMSE_C"],
        "final_confirmation_RMSE": confirmation_metrics["RMSE_C"] <= limit["final_confirmation_RMSE_C"],
        "final_confirmation_maximum": confirmation_metrics["maximum_absolute_error_C"] <= limit["final_confirmation_maximum_error_C"],
        "temperature_trend": confirmation_metrics["temperature_trend_spearman"] >= limit["temperature_trend_spearman"],
        "v12_spatial": 1.1972718732972596 <= limit["v12_spatial_RMSE_C"],
        "v13_regional_dynamic": 1.4847182706525122 <= limit["v13_regional_dynamic_RMSE_C"],
        "mesh": mesh_difference <= limit["mesh_difference_C"],
        "energy": energy["energy_residual_percent"] <= limit["energy_residual_percent"],
        "fan_exponent_strictly_inside_bounds": LOWER < exponent < UPPER,
    }
    python_pass = bool(all(checks.values()))

    rows = []
    for pwm, purpose, block in (
        (40, "calibration", calibration),
        (60, "validation", validation),
        (80, "final_confirmation", confirmation),
    ):
        times, target, predicted = block
        for index, time_s in enumerate(times):
            rows.append(
                {
                    "PWM_percent": pwm,
                    "purpose": purpose,
                    "time_s": time_s,
                    "measured_C": target[index],
                    "predicted_C": predicted[index],
                    "error_C": predicted[index] - target[index],
                    "scored": index > 0,
                }
            )
    pd.DataFrame(rows).to_csv(OUT / "multi_pwm_validation.csv", index=False)

    result = {
        "calibration_action": action,
        "frozen_model_sha256": frozen_hash,
        "first_final_confirmation_run": first_confirmation,
        "fan_flow_exponent_relative_to_60pct": exponent,
        "delivered_state_count": medium.state_count,
        "fine_state_count": fine.state_count,
        "calibration_PWM40": calibration_metrics,
        "validation_PWM60": validation_metrics,
        "final_confirmation_PWM80": confirmation_metrics,
        "independent_spatial_RMSE_C": 1.1972718732972596,
        "independent_regional_dynamic_RMSE_C": 1.4847182706525122,
        "medium_fine_max_sensor_difference_C": mesh_difference,
        "energy_audit": energy,
        "acceptance_checks_before_comsol": checks,
        "python_high_fidelity_gate_passed": python_pass,
        "comsol_confirmation_pending": True,
        "validated_scope": protocol["validated_scope_if_all_pass"] if python_pass else "none",
    }
    (OUT / "python_high_fidelity_metrics.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    if first_confirmation:
        MARKER.write_text(
            json.dumps(
                {
                    "executed_utc": datetime.now(timezone.utc).isoformat(),
                    "frozen_model_sha256": frozen_hash,
                    "python_result": "pass" if python_pass else "fail",
                    "no_retuning_after_this_marker": True,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
