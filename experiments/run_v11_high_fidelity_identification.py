"""Fit, freeze, and independently confirm the V11 40-cell reference model.

The data partition and thresholds are read from the preregistration.  Fig. 8
groups 2/4 and Fig. 9(a) are the only calibration data.  Once the frozen JSON
exists, this program never optimizes again; it only reproduces the confirmation
calculations with the recorded parameters.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models_v11.gong2024_40cell_reference import (
    LOWER_BOUNDS,
    PARAMETER_NAMES,
    PRIOR,
    PRIOR_SCALE,
    UPPER_BOUNDS,
    Gong40CellReference,
    vector_to_parameters,
)


OUT = ROOT / "outputs_v11"
DIGITIZED = OUT / "digitization"
FROZEN = OUT / "frozen_high_fidelity_model.json"
CONFIRMATION_MARKER = OUT / "final_confirmation_executed.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_static_target() -> np.ndarray:
    frame = pd.read_csv(DIGITIZED / "gong2024_fig8_75_temperatures.csv")
    frame = frame.sort_values(["cell_group", "flow_row", "in_plane_column"])
    return frame.temperature_C.to_numpy(dtype=float).reshape(5, 3, 5)


def load_dynamic(scenario: str) -> pd.DataFrame:
    frame = pd.read_csv(DIGITIZED / "gong2024_fig9_dynamic_anchors.csv")
    return frame.loc[frame.scenario == scenario].sort_values("time_s").reset_index(drop=True)


def static_calibration(target: np.ndarray) -> np.ndarray:
    calibration_groups = np.asarray((1, 3))

    def residual(five_values: np.ndarray) -> np.ndarray:
        values = PRIOR.copy()
        values[1:] = five_values
        model = Gong40CellReference(vector_to_parameters(values))
        prediction = model.sensor_field(model.steady_state(30.0, duty=0.57))
        data_residual = (prediction[calibration_groups] - target[calibration_groups]).ravel() / 0.75
        # Weak preregistered-prior stabilization prevents the two convection
        # scale parameters from wandering along their near-equivalent ridge.
        prior_residual = 0.10 * (values[1:] - PRIOR[1:]) / PRIOR_SCALE[1:]
        return np.concatenate((data_residual, prior_residual))

    result = least_squares(
        residual,
        PRIOR[1:],
        bounds=(LOWER_BOUNDS[1:] + 1e-6, UPPER_BOUNDS[1:] - 1e-6),
        xtol=2e-8,
        ftol=2e-8,
        gtol=2e-8,
        max_nfev=220,
        verbose=1,
    )
    if not result.success:
        raise RuntimeError(f"static calibration failed: {result.message}")
    values = PRIOR.copy()
    values[1:] = result.x
    return values


def dynamic_center_prediction(values: np.ndarray, frame: pd.DataFrame) -> np.ndarray:
    model = Gong40CellReference(vector_to_parameters(values))
    times = frame.time_s.to_numpy(dtype=float)
    sample_times, states = model.simulate_piecewise(
        times,
        frame.load_current_A.to_numpy(dtype=float),
        ambient_C=float(frame.center_temperature_C.iloc[0]),
        initial_temperature_C=float(frame.center_temperature_C.iloc[0]),
        sample_step_s=1.0,
    )
    indices = np.searchsorted(sample_times, times)
    # Representative "centre temperature": middle stack group, middle
    # through-flow region, and middle in-plane column.
    return states[indices, 19, 1, 2]


def thermal_mass_calibration(values: np.ndarray, frame: pd.DataFrame) -> np.ndarray:
    def objective(mass_fraction: float) -> float:
        candidate = values.copy()
        candidate[0] = mass_fraction
        error = dynamic_center_prediction(candidate, frame) - frame.center_temperature_C.to_numpy(dtype=float)
        return float(np.mean(np.square(error / 1.25)) + 0.04 * ((mass_fraction - PRIOR[0]) / PRIOR_SCALE[0]) ** 2)

    result = minimize_scalar(
        objective,
        bounds=(LOWER_BOUNDS[0] + 1e-5, UPPER_BOUNDS[0] - 1e-5),
        method="bounded",
        options={"xatol": 2e-5, "maxiter": 80},
    )
    if not result.success:
        raise RuntimeError(f"thermal-mass calibration failed: {result.message}")
    values = values.copy()
    values[0] = float(result.x)
    return values


def metric_block(prediction: np.ndarray, target: np.ndarray) -> dict[str, float]:
    error = prediction - target
    predicted_rows = np.argmax(np.max(prediction, axis=2), axis=1)
    target_rows = np.argmax(np.max(target, axis=2), axis=1)
    return {
        "spatial_RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
        "mean_bias_C": float(np.mean(error)),
        "maximum_absolute_point_error_C": float(np.max(np.abs(error))),
        "Tmax_absolute_error_C": float(abs(np.max(prediction) - np.max(target))),
        "DeltaT_absolute_error_C": float(abs(np.ptp(prediction) - np.ptp(target))),
        "hotspot_row_accuracy_percent": float(np.mean(predicted_rows == target_rows) * 100.0),
    }


def freeze(values: np.ndarray) -> dict[str, object]:
    parameters = vector_to_parameters(values)
    payload: dict[str, object] = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "preregistration_sha256": sha256(OUT / "preregistered_high_fidelity_protocol.json"),
        "calibration_data_sha256": {
            "fig8": sha256(DIGITIZED / "gong2024_fig8_75_temperatures.csv"),
            "fig9": sha256(DIGITIZED / "gong2024_fig9_dynamic_anchors.csv"),
        },
        "parameters": asdict(parameters),
        "bounds": {
            name: [float(lower), float(upper)]
            for name, lower, upper in zip(PARAMETER_NAMES, LOWER_BOUNDS, UPPER_BOUNDS, strict=True)
        },
        "calibration_only": ["Fig.8 groups 2 and 4", "Fig.9(a) UAV-like load"],
        "forbidden_for_retuning": ["Fig.8 groups 1, 3 and 5", "Fig.9(b) random load"],
    }
    FROZEN.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    payload["file_sha256"] = sha256(FROZEN)
    return payload


def load_frozen() -> tuple[np.ndarray, dict[str, object]]:
    payload = json.loads(FROZEN.read_text(encoding="utf-8"))
    values = np.asarray([payload["parameters"][name] for name in PARAMETER_NAMES], dtype=float)
    return values, payload


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    target = load_static_target()
    if FROZEN.exists():
        values, frozen_payload = load_frozen()
        calibration_action = "loaded_existing_freeze; no optimization executed"
    else:
        values = static_calibration(target)
        values = thermal_mass_calibration(values, load_dynamic("uav_like"))
        frozen_payload = freeze(values)
        calibration_action = "optimized_calibration_partition_then_frozen"

    model = Gong40CellReference(vector_to_parameters(values))
    medium_temperature = model.steady_state(30.0, duty=0.57)
    medium_sensor = model.sensor_field(medium_temperature)
    fine_model = Gong40CellReference(vector_to_parameters(values), flow_rows=6, in_plane_columns=10)
    fine_temperature = fine_model.steady_state(30.0, duty=0.57)
    fine_sensor = fine_model.sensor_field(fine_temperature)
    mesh_difference = float(np.max(np.abs(medium_sensor - fine_sensor)))

    calibration_metrics = metric_block(medium_sensor[[1, 3]], target[[1, 3]])
    validation_metrics = metric_block(medium_sensor[[0, 4]], target[[0, 4]])

    # The following partition is evaluated only after the parameter file is
    # frozen.  A marker prevents an accidental second "first look" run from
    # being misrepresented as a fresh confirmation.
    first_confirmation_run = not CONFIRMATION_MARKER.exists()
    final_spatial_metrics = metric_block(medium_sensor[[2]], target[[2]])
    random_frame = load_dynamic("random_load")
    random_prediction = dynamic_center_prediction(values, random_frame)
    random_target = random_frame.center_temperature_C.to_numpy(dtype=float)
    random_error = random_prediction - random_target
    final_dynamic_metrics = {
        "RMSE_C": float(np.sqrt(np.mean(np.square(random_error)))),
        "maximum_absolute_error_C": float(np.max(np.abs(random_error))),
        "mean_bias_C": float(np.mean(random_error)),
    }
    uav_frame = load_dynamic("uav_like")
    uav_prediction = dynamic_center_prediction(values, uav_frame)
    uav_error = uav_prediction - uav_frame.center_temperature_C.to_numpy(dtype=float)
    calibration_dynamic_metrics = {
        "RMSE_C": float(np.sqrt(np.mean(np.square(uav_error)))),
        "maximum_absolute_error_C": float(np.max(np.abs(uav_error))),
    }

    energy = model.energy_audit(medium_temperature, 30.0, 0.57, 25.0)
    protocol = json.loads((OUT / "preregistered_high_fidelity_protocol.json").read_text(encoding="utf-8"))
    limits = protocol["acceptance_thresholds"]
    checks = {
        "validation_spatial": validation_metrics["spatial_RMSE_C"] <= limits["validation_spatial_RMSE_C"],
        "validation_Tmax": validation_metrics["Tmax_absolute_error_C"] <= limits["validation_Tmax_absolute_error_C"],
        "validation_DeltaT": validation_metrics["DeltaT_absolute_error_C"] <= limits["validation_DeltaT_absolute_error_C"],
        "final_spatial": final_spatial_metrics["spatial_RMSE_C"] <= limits["final_spatial_RMSE_C"],
        "final_dynamic_RMSE": final_dynamic_metrics["RMSE_C"] <= limits["final_dynamic_RMSE_C"],
        "final_dynamic_maximum": final_dynamic_metrics["maximum_absolute_error_C"] <= limits["final_dynamic_maximum_error_C"],
        "hotspot_row": min(validation_metrics["hotspot_row_accuracy_percent"], final_spatial_metrics["hotspot_row_accuracy_percent"]) >= limits["hotspot_row_accuracy_percent"],
        "mesh": mesh_difference <= limits["medium_fine_aggregated_sensor_difference_C"],
        "energy": energy["energy_residual_percent"] <= limits["energy_residual_percent"],
        "parameters_strictly_inside_bounds": bool(np.all(values > LOWER_BOUNDS) and np.all(values < UPPER_BOUNDS)),
    }
    accepted = bool(all(checks.values()))

    prediction_rows: list[dict[str, object]] = []
    for group in range(5):
        purpose = "calibration" if group in (1, 3) else ("validation" if group in (0, 4) else "final_confirmation")
        for row in range(3):
            for column in range(5):
                prediction_rows.append(
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
    pd.DataFrame(prediction_rows).to_csv(OUT / "spatial_confirmation.csv", index=False)

    dynamic_rows: list[dict[str, object]] = []
    for frame, prediction in ((uav_frame, uav_prediction), (random_frame, random_prediction)):
        for record, predicted in zip(frame.to_dict("records"), prediction, strict=True):
            dynamic_rows.append({**record, "predicted_C": predicted, "error_C": predicted - record["center_temperature_C"]})
    pd.DataFrame(dynamic_rows).to_csv(OUT / "dynamic_confirmation.csv", index=False)

    results = {
        "calibration_action": calibration_action,
        "frozen_model_sha256": sha256(FROZEN),
        "first_final_confirmation_run": first_confirmation_run,
        "parameters": dict(zip(PARAMETER_NAMES, map(float, values), strict=True)),
        "calibration_spatial": calibration_metrics,
        "calibration_dynamic": calibration_dynamic_metrics,
        "validation_spatial": validation_metrics,
        "final_confirmation_spatial": final_spatial_metrics,
        "final_confirmation_dynamic": final_dynamic_metrics,
        "medium_fine_max_sensor_difference_C": mesh_difference,
        "energy_audit": energy,
        "acceptance_checks": checks,
        "high_fidelity_model_accepted": accepted,
        "scope": "validated 40-cell thermal reference only; no electrochemical state validation",
    }
    (OUT / "high_fidelity_validation_metrics.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )
    if first_confirmation_run:
        CONFIRMATION_MARKER.write_text(
            json.dumps(
                {
                    "executed_utc": datetime.now(timezone.utc).isoformat(),
                    "frozen_model_sha256": sha256(FROZEN),
                    "result": "pass" if accepted else "fail",
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
