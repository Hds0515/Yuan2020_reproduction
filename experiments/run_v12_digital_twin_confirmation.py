"""Identify, freeze, and confirm the V12 40-cell thermal digital twin."""

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

from models_v12.gong40cell_digital_twin import (  # noqa: E402
    LOWER,
    PARAMETER_NAMES,
    PRIOR,
    PRIOR_SCALE,
    UPPER,
    Gong40CellDigitalTwin,
    fan_duty_with_release_lag,
    vector_to_parameters,
)


OUT = ROOT / "outputs_v12"
FROZEN = OUT / "frozen_40cell_digital_twin.json"
MARKER = OUT / "final_confirmation_executed.json"
STATIC_IDS = np.asarray((3, 5, 6, 7, 8, 9))
DYNAMIC_IDS = np.asarray((0, 1, 2, 4, 10))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def static_target() -> np.ndarray:
    frame = pd.read_csv(OUT.parent / "outputs_v11/digitization/gong2024_fig8_75_temperatures.csv")
    frame = frame.sort_values(["cell_group", "flow_row", "in_plane_column"])
    return frame.temperature_C.to_numpy(float).reshape(5, 3, 5)


def scenario_data(scenario: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    events = pd.read_csv(ROOT / "digitization_v12/gong2024_fig9_input_events.csv")
    targets = pd.read_csv(ROOT / "digitization_v12/gong2024_fig9_temperature_targets.csv")
    return (
        events.loc[events.scenario == scenario].sort_values("time_s").reset_index(drop=True),
        targets.loc[targets.scenario == scenario].sort_values("time_s").reset_index(drop=True),
    )


def piecewise(events_t: np.ndarray, values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    indices = np.maximum(np.searchsorted(events_t, grid, side="right") - 1, 0)
    return np.asarray(values)[indices]


def dynamic_prediction(
    values: np.ndarray,
    scenario: str,
    *,
    flow_rows: int = 3,
    columns: int = 5,
) -> tuple[np.ndarray, pd.DataFrame, float]:
    events, targets = scenario_data(scenario)
    grid = np.unique(np.concatenate((events.time_s.to_numpy(float), targets.time_s.to_numpy(float))))
    current = piecewise(events.time_s.to_numpy(float), events.current_A.to_numpy(float), grid)
    voltage = piecewise(events.time_s.to_numpy(float), events.stack_voltage_V.to_numpy(float), grid)
    duty = fan_duty_with_release_lag(grid, current, values[10])
    model = Gong40CellDigitalTwin(
        vector_to_parameters(values), flow_rows=flow_rows, columns=columns
    )
    states = model.simulate_anchors(
        grid,
        current,
        voltage,
        ambient_C=float(events.ambient_C.iloc[0]),
        initial_temperature_C=float(targets.center_temperature_C.iloc[0]),
        duties=duty,
    )
    surface = model.surface_from_state(states)
    sensors = model.sensor_field(surface)
    target_indices = np.searchsorted(grid, targets.time_s.to_numpy(float))
    prediction = sensors[target_indices, 2, 1, 2]
    return prediction, targets, model.last_runtime_s


def fit_parameters(target: np.ndarray) -> np.ndarray:
    calibration_groups = np.asarray((0, 1, 3, 4))

    def static_residual(candidate: np.ndarray) -> np.ndarray:
        values = PRIOR.copy()
        values[STATIC_IDS] = candidate
        model = Gong40CellDigitalTwin(vector_to_parameters(values))
        surface, _ = model.steady_state(30.0, 26.0, duty=0.57, ambient_C=25.0)
        prediction = model.sensor_field(surface)
        data = (prediction[calibration_groups] - target[calibration_groups]).ravel() / 0.75
        prior = 0.10 * (values[STATIC_IDS] - PRIOR[STATIC_IDS]) / PRIOR_SCALE[STATIC_IDS]
        return np.concatenate((data, prior))

    static_result = least_squares(
        static_residual,
        PRIOR[STATIC_IDS],
        bounds=(LOWER[STATIC_IDS] + 1e-6, UPPER[STATIC_IDS] - 1e-6),
        max_nfev=180,
        xtol=2e-8,
        ftol=2e-8,
        gtol=2e-8,
        verbose=1,
    )
    if not static_result.success:
        raise RuntimeError(static_result.message)
    values = PRIOR.copy()
    values[STATIC_IDS] = static_result.x

    _, dynamic_target = scenario_data("uav_like")

    def dynamic_residual(candidate: np.ndarray) -> np.ndarray:
        trial = values.copy()
        trial[DYNAMIC_IDS] = candidate
        prediction, _, _ = dynamic_prediction(trial, "uav_like")
        data = (prediction - dynamic_target.center_temperature_C.to_numpy(float)) / 0.75
        prior = 0.10 * (trial[DYNAMIC_IDS] - PRIOR[DYNAMIC_IDS]) / PRIOR_SCALE[DYNAMIC_IDS]
        return np.concatenate((data, prior))

    dynamic_result = least_squares(
        dynamic_residual,
        PRIOR[DYNAMIC_IDS],
        bounds=(LOWER[DYNAMIC_IDS] + 1e-5, UPPER[DYNAMIC_IDS] - 1e-5),
        max_nfev=100,
        xtol=2e-7,
        ftol=2e-7,
        gtol=2e-7,
        verbose=1,
    )
    if not dynamic_result.success:
        raise RuntimeError(dynamic_result.message)
    values[DYNAMIC_IDS] = dynamic_result.x
    return values


def freeze(values: np.ndarray) -> None:
    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": sha256(OUT / "preregistered_digital_twin_protocol.json"),
        "model_source_sha256": sha256(ROOT / "models_v12/gong40cell_digital_twin.py"),
        "data_sha256": {
            "fig8": sha256(OUT.parent / "outputs_v11/digitization/gong2024_fig8_75_temperatures.csv"),
            "fig9_inputs": sha256(ROOT / "digitization_v12/gong2024_fig9_input_events.csv"),
            "fig9_temperatures": sha256(ROOT / "digitization_v12/gong2024_fig9_temperature_targets.csv"),
        },
        "parameters": asdict(vector_to_parameters(values)),
        "identification_mesh": [40, 3, 5, 2],
        "delivered_mesh": [40, 6, 10, 2],
        "fine_confirmation_mesh": [40, 12, 20, 2],
        "calibration_partition": ["Fig.8 groups 1/2/4/5", "Fig.9(a)"],
        "final_confirmation_partition": ["Fig.8 group 3", "Fig.9(b)"],
        "no_retuning_after_freeze": True,
    }
    FROZEN.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def load_frozen() -> np.ndarray:
    payload = json.loads(FROZEN.read_text(encoding="utf-8"))
    return np.asarray([payload["parameters"][name] for name in PARAMETER_NAMES])


def spatial_metrics(prediction: np.ndarray, target: np.ndarray) -> dict[str, float]:
    error = prediction - target
    predicted_rows = np.argmax(np.max(prediction, axis=2), axis=1)
    target_rows = np.argmax(np.max(target, axis=2), axis=1)
    return {
        "RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
        "mean_bias_C": float(np.mean(error)),
        "maximum_absolute_point_error_C": float(np.max(np.abs(error))),
        "Tmax_absolute_error_C": float(abs(np.max(prediction) - np.max(target))),
        "DeltaT_absolute_error_C": float(abs(np.ptp(prediction) - np.ptp(target))),
        "hotspot_row_accuracy_percent": float(np.mean(predicted_rows == target_rows) * 100.0),
    }


def main() -> None:
    target = static_target()
    if FROZEN.exists():
        values = load_frozen()
        calibration_action = "loaded_frozen_parameters; no optimization"
    else:
        values = fit_parameters(target)
        freeze(values)
        calibration_action = "calibration partitions optimized then frozen"

    parameter_file_hash = sha256(FROZEN)
    parameters = vector_to_parameters(values)
    medium = Gong40CellDigitalTwin(parameters, flow_rows=6, columns=10)
    fine = Gong40CellDigitalTwin(parameters, flow_rows=12, columns=20)
    medium_surface, medium_core = medium.steady_state(30.0, 26.0, duty=0.57, ambient_C=25.0)
    fine_surface, _ = fine.steady_state(30.0, 26.0, duty=0.57, ambient_C=25.0)
    medium_sensor = medium.sensor_field(medium_surface)
    fine_sensor = fine.sensor_field(fine_surface)
    mesh_max = float(np.max(np.abs(medium_sensor - fine_sensor)))
    calibration_spatial = spatial_metrics(medium_sensor[[0, 1, 3, 4]], target[[0, 1, 3, 4]])

    uav_prediction, uav_target, uav_runtime = dynamic_prediction(
        values, "uav_like", flow_rows=6, columns=10
    )
    uav_error = uav_prediction - uav_target.center_temperature_C.to_numpy(float)
    calibration_dynamic = {
        "RMSE_C": float(np.sqrt(np.mean(np.square(uav_error)))),
        "maximum_absolute_error_C": float(np.max(np.abs(uav_error))),
        "runtime_s": float(uav_runtime),
    }

    # First access to both V12 confirmation subsets occurs after the parameter
    # file exists.  The immutable marker records whether this was the first run.
    first_confirmation = not MARKER.exists()
    final_spatial = spatial_metrics(medium_sensor[[2]], target[[2]])
    random_prediction, random_target, random_runtime = dynamic_prediction(
        values, "random_load", flow_rows=6, columns=10
    )
    random_error = random_prediction - random_target.center_temperature_C.to_numpy(float)
    final_dynamic = {
        "RMSE_C": float(np.sqrt(np.mean(np.square(random_error)))),
        "maximum_absolute_error_C": float(np.max(np.abs(random_error))),
        "mean_bias_C": float(np.mean(random_error)),
        "runtime_s": float(random_runtime),
    }

    energy = medium.energy_audit(
        medium_surface, medium_core, 30.0, 26.0, 0.57, 25.0
    )
    protocol = json.loads((OUT / "preregistered_digital_twin_protocol.json").read_text(encoding="utf-8"))
    threshold = protocol["acceptance_thresholds"]
    checks = {
        "calibration_spatial": calibration_spatial["RMSE_C"] <= threshold["calibration_spatial_RMSE_C"],
        "calibration_dynamic": calibration_dynamic["RMSE_C"] <= threshold["calibration_dynamic_RMSE_C"],
        "final_spatial": final_spatial["RMSE_C"] <= threshold["final_spatial_RMSE_C"],
        "final_Tmax": final_spatial["Tmax_absolute_error_C"] <= threshold["final_Tmax_absolute_error_C"],
        "final_DeltaT": final_spatial["DeltaT_absolute_error_C"] <= threshold["final_DeltaT_absolute_error_C"],
        "final_dynamic_RMSE": final_dynamic["RMSE_C"] <= threshold["final_dynamic_RMSE_C"],
        "final_dynamic_maximum": final_dynamic["maximum_absolute_error_C"] <= threshold["final_dynamic_maximum_error_C"],
        "hotspot_row": final_spatial["hotspot_row_accuracy_percent"] >= threshold["hotspot_row_accuracy_percent"],
        "mesh": mesh_max <= threshold["medium_fine_sensor_difference_C"],
        "steady_energy": energy["energy_residual_percent"] <= threshold["steady_energy_residual_percent"],
        "parameters_strictly_inside_bounds": bool(np.all(values > LOWER) and np.all(values < UPPER)),
    }
    python_accepted = bool(all(checks.values()))

    spatial_rows = []
    for group in range(5):
        purpose = "final_confirmation" if group == 2 else "calibration"
        for row in range(3):
            for column in range(5):
                spatial_rows.append(
                    {
                        "cell_group": group + 1,
                        "flow_row": row + 1,
                        "in_plane_column": column + 1,
                        "purpose": purpose,
                        "measured_C": target[group, row, column],
                        "predicted_C": medium_sensor[group, row, column],
                        "error_C": medium_sensor[group, row, column] - target[group, row, column],
                    }
                )
    pd.DataFrame(spatial_rows).to_csv(OUT / "spatial_validation.csv", index=False)
    dynamic_rows = []
    for scenario, targets, prediction in (
        ("uav_like", uav_target, uav_prediction),
        ("random_load", random_target, random_prediction),
    ):
        for row, predicted in zip(targets.to_dict("records"), prediction, strict=True):
            dynamic_rows.append(
                {
                    **row,
                    "predicted_C": float(predicted),
                    "error_C": float(predicted - row["center_temperature_C"]),
                }
            )
    pd.DataFrame(dynamic_rows).to_csv(OUT / "dynamic_validation.csv", index=False)

    results = {
        "calibration_action": calibration_action,
        "frozen_parameter_file_sha256": parameter_file_hash,
        "first_final_confirmation_run": first_confirmation,
        "parameters": dict(zip(PARAMETER_NAMES, map(float, values), strict=True)),
        "delivered_state_count": medium.state_count,
        "fine_state_count": fine.state_count,
        "calibration_spatial": calibration_spatial,
        "calibration_dynamic": calibration_dynamic,
        "final_spatial": final_spatial,
        "final_dynamic": final_dynamic,
        "medium_fine_max_sensor_difference_C": mesh_max,
        "steady_energy_audit": energy,
        "acceptance_checks_before_comsol": checks,
        "python_high_fidelity_gate_passed": python_accepted,
        "comsol_independent_confirmation_pending": True,
        "control_use_authorized": False,
    }
    (OUT / "python_validation_metrics.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )
    if first_confirmation:
        MARKER.write_text(
            json.dumps(
                {
                    "executed_utc": datetime.now(timezone.utc).isoformat(),
                    "frozen_parameter_file_sha256": parameter_file_hash,
                    "python_result": "pass" if python_accepted else "fail",
                    "no_retuning_after_this_marker": True,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
