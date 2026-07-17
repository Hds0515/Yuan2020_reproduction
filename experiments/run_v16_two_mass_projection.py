"""Run the preregistered V16 five-zone/two-mass projection experiment."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models_v13.fixed_pwm_plant import FixedPWM60Plant  # noqa: E402
from models_v15.fixed60_high_resolution_reference import FROZEN_DYNAMIC  # noqa: E402
from models_v16.projected_two_mass_five_zone_rom import (  # noqa: E402
    ProjectedParameters,
    ProjectedTwoMassFiveZoneROM,
    project_full_state_to_ten_zones,
)


OUT = ROOT / "outputs_v16"
DATA = OUT / "two_mass_projection_dataset"
PROTOCOL = OUT / "preregistered_two_mass_rom_protocol.json"
FROZEN = OUT / "frozen_projected_two_mass_rom.json"
MARKER = OUT / "two_mass_synthetic_final_executed.json"
PARENT_SHA = "05566ac5903af63fdda10c1dd941a3bfbd6b7f607409007d341e320debf81ae7"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def definitions() -> list[dict[str, object]]:
    return [
        {
            "scenario_id": "two_mass_calibration_staircase",
            "purpose": "calibration",
            "seed": None,
        },
        {
            "scenario_id": "two_mass_calibration_triangular",
            "purpose": "calibration",
            "seed": None,
        },
        {
            "scenario_id": "two_mass_validation_multisine",
            "purpose": "validation",
            "seed": None,
        },
        {
            "scenario_id": "two_mass_final_prbs_seed2601",
            "purpose": "synthetic_final_confirmation",
            "seed": 2601,
        },
    ]


def inputs(name: str, seed: int | None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    times = np.arange(0.0, 800.0 + 20.0, 20.0)
    if name == "two_mass_calibration_staircase":
        current = np.select(
            [times < 100, times < 260, times < 440, times < 620],
            [5.0, 15.0, 30.0, 40.0],
            default=10.0,
        )
    elif name == "two_mass_calibration_triangular":
        current = 8.0 + 30.0 * (1.0 - np.abs(times - 400.0) / 400.0)
    elif name == "two_mass_validation_multisine":
        current = np.clip(
            23.0
            + 10.0 * np.sin(2 * np.pi * times / 320.0)
            + 5.0 * np.sin(2 * np.pi * times / 140.0),
            5.0,
            40.0,
        )
    elif name == "two_mass_final_prbs_seed2601":
        rng = np.random.default_rng(seed)
        blocks = rng.choice(np.asarray((8.0, 16.0, 24.0, 32.0, 40.0)), size=10)
        current = blocks[np.minimum((times / 80.0).astype(int), len(blocks) - 1)]
    else:
        raise KeyError(name)
    voltage = np.clip(33.5 - 0.34 * current, 19.5, 33.5)
    return times, current.astype(float), voltage.astype(float)


def generate(
    plant: FixedPWM60Plant, definition: dict[str, object]
) -> pd.DataFrame:
    name = str(definition["scenario_id"])
    times, current, voltage = inputs(name, definition["seed"])
    state = plant.simulate_anchors(
        times,
        current,
        voltage,
        ambient_C=23.0,
        initial_temperature_C=23.0,
        duties=np.full_like(times, 0.60),
    )
    reduced = project_full_state_to_ten_zones(plant, state)
    surface = plant.surface_from_state(state)
    flat_index = np.argmax(surface.reshape(len(times), -1), axis=1)
    hot_row = np.asarray(
        np.unravel_index(flat_index, (plant.cells, plant.flow_rows, plant.columns))
    ).T[:, 1]
    hotspot_zone = hot_row // (plant.flow_rows // 5)
    frame = pd.DataFrame(
        {
            "scenario_id": name,
            "purpose": str(definition["purpose"]),
            "time_s": times,
            "current_A": current,
            "stack_voltage_V": voltage,
            "ambient_C": 23.0,
            "fan_duty": 0.60,
            **{f"surface_zone_{i + 1}_C": reduced[:, i] for i in range(5)},
            **{f"core_zone_{i + 1}_C": reduced[:, 5 + i] for i in range(5)},
            "Tmax_C": surface.max(axis=(1, 2, 3)),
            "Tmin_C": surface.min(axis=(1, 2, 3)),
            "hotspot_zone_zero_based": hotspot_zone,
        }
    )
    frame["DeltaT_C"] = frame.Tmax_C - frame.Tmin_C
    return frame


def fit_offsets(frames: list[pd.DataFrame]) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.concat(frames, ignore_index=True)
    surface = frame[[f"surface_zone_{i}_C" for i in range(1, 6)]].to_numpy(float)
    basis = np.column_stack((np.ones(len(frame)), frame.current_A.to_numpy(float) / 40.0))
    maximum = np.maximum(frame.Tmax_C.to_numpy(float) - surface.max(axis=1), 0.0)
    minimum = np.maximum(surface.min(axis=1) - frame.Tmin_C.to_numpy(float), 0.0)
    fit_max = least_squares(
        lambda values: basis @ values - maximum,
        np.asarray((0.1, 0.5)),
        bounds=(0.0, np.inf),
    ).x
    fit_min = least_squares(
        lambda values: basis @ values - minimum,
        np.asarray((0.1, 0.5)),
        bounds=(0.0, np.inf),
    ).x
    return fit_max, fit_min


def load_parameters() -> ProjectedParameters:
    payload = json.loads(FROZEN.read_text(encoding="utf-8"))["parameters"]
    return ProjectedParameters(
        np.asarray(payload["capacity_J_K"], dtype=float),
        np.asarray(payload["thermal_operator_W_K"], dtype=float),
        np.asarray(payload["tmax_offset_C_intercept_and_current40"], dtype=float),
        np.asarray(payload["tmin_offset_C_intercept_and_current40"], dtype=float),
    )


def evaluate(
    model: ProjectedTwoMassFiveZoneROM, frames: list[pd.DataFrame]
) -> tuple[dict[str, object], pd.DataFrame]:
    rows: list[pd.DataFrame] = []
    runtime = 0.0
    for frame in frames:
        initial = frame[
            [
                *[f"surface_zone_{i}_C" for i in range(1, 6)],
                *[f"core_zone_{i}_C" for i in range(1, 6)],
            ]
        ].iloc[0].to_numpy(float)
        start = perf_counter()
        state = model.simulate(
            frame.time_s.to_numpy(float),
            frame.current_A.to_numpy(float),
            frame.stack_voltage_V.to_numpy(float),
            initial,
        )
        runtime += perf_counter() - start
        block = frame.copy()
        for i in range(5):
            block[f"ROM_surface_zone_{i + 1}_C"] = state[:, i]
            block[f"surface_zone_{i + 1}_error_C"] = (
                state[:, i] - frame[f"surface_zone_{i + 1}_C"].to_numpy(float)
            )
            block[f"ROM_core_zone_{i + 1}_C"] = state[:, 5 + i]
        extremes = [
            model.extreme_temperatures(item, load)
            for item, load in zip(state, frame.current_A.to_numpy(float), strict=True)
        ]
        block[["ROM_Tmax_C", "ROM_Tmin_C", "ROM_DeltaT_C"]] = np.asarray(extremes)
        block["ROM_hotspot_zone_zero_based"] = np.argmax(state[:, :5], axis=1)
        rows.append(block)
    result = pd.concat(rows, ignore_index=True)
    scoring = result.time_s > 0.0
    error = result.loc[
        scoring, [f"surface_zone_{i}_error_C" for i in range(1, 6)]
    ].to_numpy(float)
    maximum_error = result.loc[scoring, "ROM_Tmax_C"].to_numpy(float) - result.loc[
        scoring, "Tmax_C"
    ].to_numpy(float)
    delta_error = result.loc[scoring, "ROM_DeltaT_C"].to_numpy(float) - result.loc[
        scoring, "DeltaT_C"
    ].to_numpy(float)
    predicted_mean = result.loc[
        scoring, [f"ROM_surface_zone_{i}_C" for i in range(1, 6)]
    ].to_numpy(float).mean(axis=1)
    target_mean = result.loc[
        scoring, [f"surface_zone_{i}_C" for i in range(1, 6)]
    ].to_numpy(float).mean(axis=1)
    terminal = []
    for _, block in result.loc[scoring].groupby("scenario_id"):
        count = max(int(np.ceil(len(block) / 4)), 1)
        terminal.append(
            float(
                np.max(
                    np.abs(
                        block[
                            [f"surface_zone_{i}_error_C" for i in range(1, 6)]
                        ].tail(count).mean(axis=0)
                    )
                )
            )
        )
    energy = 0.0
    for _, row in result.loc[scoring].iterrows():
        state = row[
            [
                *[f"ROM_surface_zone_{i}_C" for i in range(1, 6)],
                *[f"ROM_core_zone_{i}_C" for i in range(1, 6)],
            ]
        ].to_numpy(float)
        energy = max(
            energy,
            model.instantaneous_energy_audit(
                state, float(row.current_A), float(row.stack_voltage_V)
            )["energy_residual_percent"],
        )
    metrics = {
        "scenario_ids": sorted(result.scenario_id.unique().tolist()),
        "zone_overall_RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
        "each_zone_RMSE_C": [
            float(value) for value in np.sqrt(np.mean(np.square(error), axis=0))
        ],
        "Tmax_RMSE_C": float(np.sqrt(np.mean(np.square(maximum_error)))),
        "DeltaT_RMSE_C": float(np.sqrt(np.mean(np.square(delta_error)))),
        "mean_temperature_RMSE_C": float(
            np.sqrt(np.mean(np.square(predicted_mean - target_mean)))
        ),
        "hotspot_zone_accuracy_percent": float(
            np.mean(
                result.loc[scoring, "ROM_hotspot_zone_zero_based"].to_numpy(int)
                == result.loc[scoring, "hotspot_zone_zero_based"].to_numpy(int)
            )
            * 100.0
        ),
        "maximum_terminal_quarter_bias_C": float(max(terminal)),
        "maximum_energy_residual_percent": float(energy),
        "mean_runtime_per_step_s": float(runtime / max(scoring.sum(), 1)),
        "maximum_real_eigenvalue_per_s": float(
            np.max(np.real(np.linalg.eigvals(model.continuous_matrix)))
        ),
    }
    return metrics, result


def checks(metrics: dict[str, object], protocol: dict[str, object]) -> dict[str, bool]:
    threshold = protocol["acceptance_thresholds"]
    mapping = {
        "zone_overall_RMSE_C": "zone_overall_RMSE_C",
        "Tmax_RMSE_C": "Tmax_RMSE_C",
        "DeltaT_RMSE_C": "DeltaT_RMSE_C",
        "hotspot_zone_accuracy_percent": "hotspot_zone_accuracy_percent",
        "mean_temperature_RMSE_C": "mean_temperature_RMSE_C",
        "energy_residual_percent": "maximum_energy_residual_percent",
        "maximum_terminal_quarter_bias_C": "maximum_terminal_quarter_bias_C",
    }
    result = {}
    for name, limit in threshold.items():
        value = float(metrics[mapping[name]])
        result[name] = value >= float(limit) if name == "hotspot_zone_accuracy_percent" else value <= float(limit)
    return result


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    plant = FixedPWM60Plant(FROZEN_DYNAMIC, flow_rows=15, columns=5)

    # Calibration parent data are generated first.  Only the extrema output
    # reconstruction uses them; the ten-state dynamics are a direct projection.
    frames: dict[str, pd.DataFrame] = {}
    registry_rows = []
    for definition in definitions()[:2]:
        path = DATA / f"{definition['scenario_id']}.csv"
        frame = pd.read_csv(path) if path.exists() else generate(plant, definition)
        if not path.exists():
            frame.to_csv(path, index=False)
        frames[str(definition["scenario_id"])] = frame
        registry_rows.append((definition, path))

    if FROZEN.exists():
        parameters = load_parameters()
        action = "loaded existing projection freeze; no fitting"
    else:
        base = ProjectedTwoMassFiveZoneROM.project_from_v15(plant)
        maximum, minimum = fit_offsets(list(frames.values()))
        parameters = ProjectedParameters(
            base.parameters.capacity_J_K,
            base.parameters.thermal_operator_W_K,
            maximum,
            minimum,
        )
        write_json(
            FROZEN,
            {
                "schema_version": 1,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "parent_V15_sha256": PARENT_SHA,
                "protocol_sha256": sha256(PROTOCOL),
                "source_sha256": sha256(
                    ROOT / "models_v16/projected_two_mass_five_zone_rom.py"
                ),
                "dynamic_parameter_fit_count": 0,
                "extrema_reconstruction_fit_count": 4,
                "parameters": parameters.to_dict(),
                "no_retuning_after_freeze": True,
            },
        )
        action = "Galerkin projection plus calibration-only extrema offsets, then frozen"
    model = ProjectedTwoMassFiveZoneROM(parameters)
    model.assert_physical()
    frozen_hash = sha256(FROZEN)

    # Validation is generated and scored without changing the freeze.
    validation_definition = definitions()[2]
    validation_path = DATA / f"{validation_definition['scenario_id']}.csv"
    validation = (
        pd.read_csv(validation_path)
        if validation_path.exists()
        else generate(plant, validation_definition)
    )
    if not validation_path.exists():
        validation.to_csv(validation_path, index=False)
    frames[str(validation_definition["scenario_id"])] = validation
    registry_rows.append((validation_definition, validation_path))

    # The new seed-2601 scenario is executed once only after the model freeze.
    final_definition = definitions()[3]
    final_path = DATA / f"{final_definition['scenario_id']}.csv"
    first_final = not MARKER.exists()
    final = pd.read_csv(final_path) if final_path.exists() else generate(plant, final_definition)
    if not final_path.exists():
        final.to_csv(final_path, index=False)
    frames[str(final_definition["scenario_id"])] = final
    registry_rows.append((final_definition, final_path))

    calibration_metrics, calibration_rows = evaluate(
        model,
        [
            frames["two_mass_calibration_staircase"],
            frames["two_mass_calibration_triangular"],
        ],
    )
    validation_metrics, validation_rows = evaluate(model, [validation])
    final_metrics, final_rows = evaluate(model, [final])
    pd.concat((calibration_rows, validation_rows, final_rows), ignore_index=True).to_csv(
        OUT / "projected_two_mass_rom_predictions.csv", index=False
    )
    all_checks = {
        "calibration": checks(calibration_metrics, protocol),
        "validation": checks(validation_metrics, protocol),
        "synthetic_final_confirmation": checks(final_metrics, protocol),
    }
    method_pass = all(all_checks["validation"].values()) and all(
        all_checks["synthetic_final_confirmation"].values()
    )
    summary = {
        "action": action,
        "frozen_ROM_sha256": frozen_hash,
        "first_synthetic_final_confirmation_run": first_final,
        "parent_state_count": plant.state_count,
        "ROM_state_count": model.state_count,
        "dynamic_parameter_fit_count": 0,
        "calibration": calibration_metrics,
        "validation": validation_metrics,
        "synthetic_final_confirmation": final_metrics,
        "acceptance_checks": all_checks,
        "synthetic_method_gate_passed": bool(method_pass),
        "independent_experimental_validation_passed": False,
        "formal_ROM_claim_authorized": False,
        "reason": "The strict parent experimental gate is failed and the new final scenario is deterministic parent simulation, not fresh experiment.",
    }
    write_json(OUT / "projected_two_mass_rom_metrics.json", summary)
    write_json(
        OUT / "two_mass_dataset_registry.json",
        {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "parent_V15_sha256": PARENT_SHA,
            "scenarios": [
                {
                    **definition,
                    "file": str(path.relative_to(ROOT)).replace("\\", "/"),
                    "sha256": sha256(path),
                    "source_type": "frozen V15-equation deterministic simulation",
                    "fresh_experiment": False,
                }
                for definition, path in registry_rows
            ],
        },
    )
    if first_final:
        write_json(
            MARKER,
            {
                "executed_utc": datetime.now(timezone.utc).isoformat(),
                "frozen_ROM_sha256": frozen_hash,
                "synthetic_result": "pass" if method_pass else "fail",
                "no_retuning_after_marker": True,
            },
        )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
