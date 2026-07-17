"""Freeze V13 on the 2 A ramp and confirm once on the 4 A ramp."""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models_v13.fixed_pwm_plant import (  # noqa: E402
    LOWER,
    NAMES,
    PRIOR,
    SCALE,
    UPPER,
    FixedPWM60Plant,
    vector_to_dynamic,
)


OUT = ROOT / "outputs_v13"
FROZEN = OUT / "frozen_fixed_pwm_high_fidelity_model.json"
MARKER = OUT / "final_confirmation_executed.json"
REGIONS = ("inlet", "intermediate", "outlet")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def data_for(scenario: str) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.read_csv(OUT / "digitization/guo2024_fixed_pwm_regional_temperatures.csv")
    frame = frame.loc[frame.scenario == scenario]
    times = np.sort(frame.time_s.unique()).astype(float)
    target = np.stack(
        [
            frame.loc[frame.flow_region == region]
            .sort_values("time_s")
            .temperature_C.to_numpy(float)
            for region in REGIONS
        ],
        axis=1,
    )
    return times, target


def current_events(scenario: str) -> tuple[np.ndarray, np.ndarray]:
    increment, interval = (2.0, 50.0) if scenario == "ramp_2A" else (4.0, 100.0)
    events = np.arange(0.0, 2000.0 + interval, interval)
    peak_start = 1000.0 - interval
    current = np.asarray(
        [
            min(increment * (index + 1), 40.0)
            if time < 1000.0
            else max(40.0 - increment * int((time - peak_start) / interval), 0.0)
            for index, time in enumerate(events)
        ]
    )
    return events, current


def step_values(event_times: np.ndarray, values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    indices = np.maximum(np.searchsorted(event_times, grid, side="right") - 1, 0)
    return values[indices]


def prediction(
    values: np.ndarray,
    scenario: str,
    *,
    flow_rows: int = 3,
    columns: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    target_times, target = data_for(scenario)
    event_times, event_current = current_events(scenario)
    grid = np.unique(np.concatenate((event_times[event_times >= target_times[0]], target_times)))
    current = step_values(event_times, event_current, grid)
    voltage = np.clip(33.5 - 0.34 * current, 19.5, 33.5)
    plant = FixedPWM60Plant(vector_to_dynamic(values), flow_rows=flow_rows, columns=columns)
    state = plant.simulate_from_regional_initial_state(
        grid, current, voltage, target[0], ambient_C=23.0, duty=0.60
    )
    sensors = plant.sensor_field(plant.surface_from_state(state))
    predicted_regions = sensors[np.searchsorted(grid, target_times)].mean(axis=(1, 3))
    return target_times, target, predicted_regions


def fit() -> np.ndarray:
    _, target, _ = prediction(PRIOR, "ramp_2A")
    scoring = np.arange(1, len(target))

    def residual(values: np.ndarray) -> np.ndarray:
        _, _, predicted = prediction(values, "ramp_2A")
        data = ((predicted[scoring] - target[scoring]) / 0.75).ravel()
        prior = 0.06 * (values - PRIOR) / SCALE
        return np.concatenate((data, prior))

    result = least_squares(
        residual,
        PRIOR,
        bounds=(LOWER + 1e-5, UPPER - 1e-5),
        max_nfev=80,
        xtol=1e-7,
        ftol=1e-7,
        gtol=1e-7,
        verbose=1,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.x


def freeze(values: np.ndarray) -> None:
    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": sha256(OUT / "preregistered_fixed_pwm_protocol.json"),
        "model_source_sha256": sha256(ROOT / "models_v13/fixed_pwm_plant.py"),
        "digitized_data_sha256": sha256(
            OUT / "digitization/guo2024_fixed_pwm_regional_temperatures.csv"
        ),
        "dynamic_parameters": asdict(vector_to_dynamic(values)),
        "plant_record": FixedPWM60Plant(vector_to_dynamic(values), flow_rows=6, columns=10).v13_record(),
        "calibration": "Guo2024 Fig.5(e), 2 A ramp",
        "final_confirmation": "Guo2024 Fig.5(f), 4 A ramp",
        "no_retuning_after_freeze": True,
    }
    FROZEN.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def load() -> np.ndarray:
    payload = json.loads(FROZEN.read_text(encoding="utf-8"))
    return np.asarray([payload["dynamic_parameters"][name] for name in NAMES])


def metrics(target_times: np.ndarray, target: np.ndarray, predicted: np.ndarray) -> dict[str, object]:
    score = np.arange(1, len(target))
    error = predicted[score] - target[score]
    region_rmse = np.sqrt(np.mean(np.square(error), axis=0))
    predicted_peak = target_times[score][np.argmax(np.max(predicted[score], axis=1))]
    target_peak = target_times[score][np.argmax(np.max(target[score], axis=1))]
    target_order = np.argsort(target[score], axis=1)
    predicted_order = np.argsort(predicted[score], axis=1)
    return {
        "regional_RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
        "each_region_RMSE_C": dict(zip(REGIONS, map(float, region_rmse), strict=True)),
        "maximum_absolute_error_C": float(np.max(np.abs(error))),
        "mean_bias_C": float(np.mean(error)),
        "peak_time_error_s": float(abs(predicted_peak - target_peak)),
        "region_order_accuracy_percent": float(np.mean(np.all(target_order == predicted_order, axis=1)) * 100.0),
    }


def main() -> None:
    if FROZEN.exists():
        values = load()
        action = "loaded existing freeze; no optimization"
    else:
        values = fit()
        freeze(values)
        action = "2 A calibration optimized then frozen"
    frozen_hash = sha256(FROZEN)

    calibration = prediction(values, "ramp_2A", flow_rows=6, columns=10)
    calibration_metrics = metrics(*calibration)
    first_confirmation = not MARKER.exists()
    confirmation = prediction(values, "ramp_4A", flow_rows=6, columns=10)
    confirmation_metrics = metrics(*confirmation)

    medium = FixedPWM60Plant(vector_to_dynamic(values), flow_rows=6, columns=10)
    fine = FixedPWM60Plant(vector_to_dynamic(values), flow_rows=12, columns=20)
    voltage_40A = float(np.clip(33.5 - 0.34 * 40.0, 19.5, 33.5))
    medium_surface, medium_core = medium.steady_state(
        40.0, voltage_40A, duty=0.60, ambient_C=23.0
    )
    fine_surface, _ = fine.steady_state(40.0, voltage_40A, duty=0.60, ambient_C=23.0)
    mesh_difference = float(
        np.max(np.abs(medium.sensor_field(medium_surface) - fine.sensor_field(fine_surface)))
    )
    energy = medium.energy_audit(
        medium_surface, medium_core, 40.0, voltage_40A, 0.60, 23.0
    )
    protocol = json.loads((OUT / "preregistered_fixed_pwm_protocol.json").read_text())
    threshold = protocol["acceptance_thresholds"]
    checks = {
        "calibration_regional_RMSE": calibration_metrics["regional_RMSE_C"] <= threshold["calibration_regional_RMSE_C"],
        "confirmation_regional_RMSE": confirmation_metrics["regional_RMSE_C"] <= threshold["confirmation_regional_RMSE_C"],
        "confirmation_each_region_RMSE": max(confirmation_metrics["each_region_RMSE_C"].values()) <= threshold["confirmation_each_region_RMSE_C"],
        "confirmation_maximum_error": confirmation_metrics["maximum_absolute_error_C"] <= threshold["confirmation_maximum_error_C"],
        "confirmation_peak_time": confirmation_metrics["peak_time_error_s"] <= threshold["confirmation_peak_time_error_s"],
        "region_order": confirmation_metrics["region_order_accuracy_percent"] >= threshold["region_order_accuracy_percent"],
        "mesh": mesh_difference <= threshold["medium_fine_sensor_difference_C"],
        "energy": energy["energy_residual_percent"] <= threshold["energy_residual_percent"],
        "parameters_strictly_inside_bounds": bool(np.all(values > LOWER) and np.all(values < UPPER)),
    }
    python_pass = bool(all(checks.values()))

    rows = []
    for scenario, purpose, block in (
        ("ramp_2A", "calibration", calibration),
        ("ramp_4A", "final_confirmation", confirmation),
    ):
        times, target, predicted = block
        for ti, time_s in enumerate(times):
            for ri, region in enumerate(REGIONS):
                rows.append(
                    {
                        "scenario": scenario,
                        "purpose": purpose,
                        "time_s": time_s,
                        "flow_region": region,
                        "measured_C": target[ti, ri],
                        "predicted_C": predicted[ti, ri],
                        "error_C": predicted[ti, ri] - target[ti, ri],
                        "scored": ti > 0,
                    }
                )
    pd.DataFrame(rows).to_csv(OUT / "fixed_pwm_dynamic_validation.csv", index=False)

    result = {
        "calibration_action": action,
        "frozen_model_sha256": frozen_hash,
        "first_final_confirmation_run": first_confirmation,
        "dynamic_parameters": dict(zip(NAMES, map(float, values), strict=True)),
        "delivered_state_count": medium.state_count,
        "fine_state_count": fine.state_count,
        "calibration": calibration_metrics,
        "final_confirmation": confirmation_metrics,
        "medium_fine_max_sensor_difference_C": mesh_difference,
        "energy_audit": energy,
        "acceptance_checks_before_comsol": checks,
        "python_high_fidelity_gate_passed": python_pass,
        "comsol_confirmation_pending": True,
        "validated_scope": protocol["validated_scope_if_passed"] if python_pass else "none",
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
