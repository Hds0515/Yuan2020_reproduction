"""Execute the gated V16 fixed-duty ROM and observer research programme.

This script never changes V15 parameters.  The only experimental final
confirmation available in the repository has already been consumed, so V16
separates the strict experimental decision from a clearly labelled synthetic
method-development track derived from the frozen V15 reference.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.run_v13_fixed_pwm_confirmation import current_events, step_values  # noqa: E402
from models_v15 import Fixed60HighResolutionReference  # noqa: E402
from models_v16 import (  # noqa: E402
    FiveZoneParameters,
    PhysicsConstrainedFiveZoneROM,
    conservative_five_zone_average,
)
from models_v16.physics_constrained_five_zone_rom import (  # noqa: E402
    TOTAL_HEAT_CAPACITY_J_K,
)


OUT = ROOT / "outputs_v16"
DATA = OUT / "high_resolution_method_dataset"
DOCS = ROOT / "docs_v16"
FROZEN = OUT / "frozen_five_zone_rom.json"
FINAL_MARKER = OUT / "synthetic_final_confirmation_executed.json"
V15_SHA256 = "05566ac5903af63fdda10c1dd941a3bfbd6b7f607409007d341e320debf81ae7"
REGIONS = ("inlet", "intermediate", "outlet")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def voltage_proxy(current_A: np.ndarray) -> np.ndarray:
    return np.clip(33.5 - 0.34 * np.asarray(current_A, dtype=float), 19.5, 33.5)


def scenario_definitions() -> list[dict[str, object]]:
    return [
        {
            "scenario_id": "calibration_staircase",
            "purpose": "calibration",
            "duration_s": 800,
            "dt_s": 20,
            "profile": "5A/15A/30A/40A/10A deterministic staircase",
            "random_seed": None,
        },
        {
            "scenario_id": "calibration_triangular",
            "purpose": "calibration",
            "duration_s": 800,
            "dt_s": 20,
            "profile": "8A to 38A to 8A triangular load",
            "random_seed": None,
        },
        {
            "scenario_id": "validation_multisine",
            "purpose": "validation",
            "duration_s": 800,
            "dt_s": 20,
            "profile": "bounded deterministic two-frequency load",
            "random_seed": None,
        },
        {
            "scenario_id": "final_prbs_seed1601",
            "purpose": "synthetic_final_confirmation",
            "duration_s": 800,
            "dt_s": 20,
            "profile": "80 s block PRBS over [8,16,24,32,40] A",
            "random_seed": 1601,
        },
    ]


def scenario_inputs(definition: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dt = float(definition["dt_s"])
    duration = float(definition["duration_s"])
    times = np.arange(0.0, duration + dt, dt)
    name = str(definition["scenario_id"])
    if name == "calibration_staircase":
        current = np.select(
            [times < 100, times < 260, times < 440, times < 620],
            [5.0, 15.0, 30.0, 40.0],
            default=10.0,
        )
    elif name == "calibration_triangular":
        current = 8.0 + 30.0 * (1.0 - np.abs(times - 400.0) / 400.0)
    elif name == "validation_multisine":
        current = 23.0 + 10.0 * np.sin(2 * np.pi * times / 320.0) + 5.0 * np.sin(
            2 * np.pi * times / 140.0
        )
        current = np.clip(current, 5.0, 40.0)
    elif name == "final_prbs_seed1601":
        rng = np.random.default_rng(int(definition["random_seed"]))
        blocks = rng.choice(np.asarray((8.0, 16.0, 24.0, 32.0, 40.0)), size=10)
        current = blocks[np.minimum((times / 80.0).astype(int), len(blocks) - 1)]
    else:
        raise KeyError(name)
    return times, current.astype(float), voltage_proxy(current)


def generate_parent_scenario(
    plant: Fixed60HighResolutionReference, definition: dict[str, object]
) -> pd.DataFrame:
    times, current, voltage = scenario_inputs(definition)
    duty = np.full_like(times, 0.60)
    state = plant.simulate_anchors(
        times,
        current,
        voltage,
        ambient_C=23.0,
        initial_temperature_C=23.0,
        duties=duty,
    )
    surface = plant.surface_from_state(state)
    zones = conservative_five_zone_average(surface)
    maximum = surface.max(axis=(1, 2, 3))
    minimum = surface.min(axis=(1, 2, 3))
    flat_index = np.argmax(surface.reshape(len(times), -1), axis=1)
    unravelled = np.asarray(
        np.unravel_index(flat_index, (plant.cells, plant.flow_rows, plant.columns))
    ).T
    hot_row = unravelled[:, 1]
    hotspot_zone = np.minimum(
        (5.0 * (hot_row + 0.5) / plant.flow_rows).astype(int), 4
    )
    frame = pd.DataFrame(
        {
            "scenario_id": str(definition["scenario_id"]),
            "purpose": str(definition["purpose"]),
            "time_s": times,
            "current_A": current,
            "stack_voltage_V": voltage,
            "ambient_C": 23.0,
            "fan_duty": 0.60,
            **{f"T_zone_{index + 1}_C": zones[:, index] for index in range(5)},
            "Tmax_C": maximum,
            "Tmin_C": minimum,
            "DeltaT_C": maximum - minimum,
            "hotspot_zone_zero_based": hotspot_zone,
        }
    )
    return frame


def create_protocol() -> dict[str, object]:
    path = OUT / "preregistered_v16_protocol.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    protocol = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parent_V15_frozen_sha256": V15_SHA256,
        "scope": {
            "ambient_C": 23.0,
            "fan_duty": 0.60,
            "current_A": [2.0, 40.0],
            "measured_current_and_voltage_are_exogenous": True,
        },
        "experimental_data_status": {
            "Guo2024_2A_ramp": "already used for V13 calibration",
            "Guo2024_4A_ramp": "already consumed once as V13 final confirmation",
            "unused_fresh_independent_dynamic_experiment_available": False,
        },
        "method_dataset_scenarios": scenario_definitions(),
        "ROM_structure_frozen_before_validation": {
            "zones": "five equal-volume cathode-flow zones",
            "state_model": "stable positive RC chain",
            "total_heat_capacity_preserved": True,
            "heat_generation_preserved": True,
            "fit_data": "calibration scenarios only",
            "validation_retuning_allowed": False,
            "synthetic_final_retuning_allowed": False,
        },
        "ROM_acceptance_thresholds": {
            "zone_overall_RMSE_C": 0.75,
            "Tmax_RMSE_C": 0.75,
            "DeltaT_RMSE_C": 0.50,
            "hotspot_zone_accuracy_percent": 90.0,
            "mean_temperature_RMSE_C": 0.50,
            "energy_residual_percent": 0.5,
            "maximum_terminal_quarter_bias_C": 0.50,
        },
        "observer_preregistration": {
            "sensor_zone_indices_zero_based": [0, 4],
            "measurement_noise_seeds": list(range(10)),
            "nominal_sigma_C": 0.12,
            "pure_noise_sigma_C": 0.30,
            "model_bank": [
                "nominal",
                "adjacent_conductance_x0.75",
                "adjacent_conductance_x1.25",
                "ambient_loss_x0.85",
                "ambient_loss_x1.15",
            ],
            "gate_inherited_from_V6": {
                "probability_floor": 0.005,
                "new_probability_fraction_per_step": 0.35,
                "normalized_entropy_maximum": 0.85,
                "maximum_probability_minimum": 0.45,
                "sustained_excess_NIS_threshold": 1.50,
            },
            "acceptance": {
                "majority_improves": True,
                "median_Tmax_RMSE_improvement_percent_min": 10.0,
                "paired_95pct_CI_lower_above_zero": True,
                "nominal_and_pure_noise_loss_percent_max": 5.0,
                "hotspot_accuracy_must_not_decrease": True,
            },
        },
        "interpretation_rule": "Synthetic parent-model confirmation can test reduction and software, but cannot repair or replace the failed experimental high-fidelity gate.",
    }
    write_json(path, protocol)
    return protocol


def diagnose_v15_experimental_residuals() -> dict[str, object]:
    frame = pd.read_csv(ROOT / "outputs_v13/fixed_pwm_dynamic_validation.csv")
    frame = frame.loc[
        (frame.purpose == "final_confirmation") & frame.scored.astype(bool)
    ].copy()
    events, event_current = current_events("ramp_4A")
    frame["current_A"] = step_values(events, event_current, frame.time_s.to_numpy(float))
    frame["phase"] = np.where(frame.time_s < 1000.0, "increasing", "decreasing")
    frame["absolute_error_C"] = frame.error_C.abs()
    frame.to_csv(OUT / "v15_dynamic_residual_diagnosis.csv", index=False)

    region_rows: list[dict[str, object]] = []
    for region, block in frame.groupby("flow_region", sort=False):
        for phase, phase_block in block.groupby("phase", sort=False):
            error = phase_block.error_C.to_numpy(float)
            correlation = float(
                np.corrcoef(phase_block.current_A.to_numpy(float), error)[0, 1]
            )
            region_rows.append(
                {
                    "flow_region": region,
                    "phase": phase,
                    "mean_bias_C": float(np.mean(error)),
                    "RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
                    "maximum_absolute_error_C": float(np.max(np.abs(error))),
                    "error_current_correlation": correlation,
                }
            )
    pd.DataFrame(region_rows).to_csv(OUT / "v15_residual_by_region_phase.csv", index=False)
    worst = frame.loc[frame.absolute_error_C.idxmax()]
    target_order = []
    predicted_order = []
    order_failures: list[float] = []
    for time_s, block in frame.groupby("time_s", sort=True):
        ordered = block.set_index("flow_region").loc[list(REGIONS)]
        target_order.append(np.argsort(ordered.measured_C.to_numpy(float)))
        predicted_order.append(np.argsort(ordered.predicted_C.to_numpy(float)))
        if not np.array_equal(target_order[-1], predicted_order[-1]):
            order_failures.append(float(time_s))
    return {
        "source": "outputs_v13/fixed_pwm_dynamic_validation.csv",
        "data_already_consumed_as_final_confirmation": True,
        "worst_error": {
            "time_s": float(worst.time_s),
            "flow_region": str(worst.flow_region),
            "current_A": float(worst.current_A),
            "measured_C": float(worst.measured_C),
            "predicted_C": float(worst.predicted_C),
            "error_C": float(worst.error_C),
        },
        "regional_order_failure_times_s": order_failures,
        "interpretation": [
            "The largest error occurs at the high-load end of the increasing ramp in the outlet region.",
            "The intermediate region is biased low during much of the increasing ramp, while inlet/outlet errors change sign with load.",
            "This load- and position-dependent residual is consistent with missing heat/airflow spatial dynamics; it is not explained by the 0.0074 degC COMSOL-Python implementation difference.",
            "Sparse 100 s digitized anchors cannot uniquely identify sensor lag, fan lag, water state and local heat release.",
        ],
    }


def decode_dynamic(vector: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    logits = np.r_[0.0, vector[:4]]
    fractions = np.exp(logits - np.max(logits))
    capacity = TOTAL_HEAT_CAPACITY_J_K * fractions / fractions.sum()
    conductance = np.exp(vector[4:8])
    loss = np.exp(vector[8:13])
    return capacity, conductance, loss


def fit_extreme_offsets(calibration: list[pd.DataFrame]) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.concat(calibration, ignore_index=True)
    zones = frame[[f"T_zone_{index}_C" for index in range(1, 6)]].to_numpy(float)
    basis = np.column_stack((np.ones(len(frame)), frame.current_A.to_numpy(float) / 40.0))
    maximum_delta = np.maximum(frame.Tmax_C.to_numpy(float) - zones.max(axis=1), 0.0)
    minimum_delta = np.maximum(zones.min(axis=1) - frame.Tmin_C.to_numpy(float), 0.0)
    fit_max = least_squares(
        lambda values: basis @ values - maximum_delta,
        np.asarray((0.1, 0.5)),
        bounds=(0.0, np.inf),
    ).x
    fit_min = least_squares(
        lambda values: basis @ values - minimum_delta,
        np.asarray((0.1, 0.5)),
        bounds=(0.0, np.inf),
    ).x
    return fit_max, fit_min


def fit_rom(calibration: list[pd.DataFrame]) -> FiveZoneParameters:
    initial = np.r_[
        np.zeros(4),
        np.log(np.full(4, 10.0)),
        np.log(np.full(5, 10.0)),
    ]
    lower = np.r_[
        np.full(4, -2.5),
        np.log(np.full(4, 0.05)),
        np.log(np.full(5, 0.05)),
    ]
    upper = np.r_[
        np.full(4, 2.5),
        np.log(np.full(4, 100.0)),
        np.log(np.full(5, 100.0)),
    ]

    def residual(vector: np.ndarray) -> np.ndarray:
        capacity, conductance, loss = decode_dynamic(vector)
        parameters = FiveZoneParameters(
            capacity,
            conductance,
            loss,
            np.zeros(2),
            np.zeros(2),
        )
        rom = PhysicsConstrainedFiveZoneROM(parameters)
        blocks: list[np.ndarray] = []
        for frame in calibration:
            target = frame[[f"T_zone_{index}_C" for index in range(1, 6)]].to_numpy(float)
            prediction = rom.simulate(
                frame.time_s.to_numpy(float),
                frame.current_A.to_numpy(float),
                frame.stack_voltage_V.to_numpy(float),
                target[0],
            )
            blocks.append(((prediction[1:] - target[1:]) / 0.25).ravel())
        regularization = np.r_[
            0.03 * vector[:4],
            0.02 * (vector[4:8] - np.log(10.0)),
            0.02 * (vector[8:13] - np.log(10.0)),
        ]
        return np.concatenate((*blocks, regularization))

    result = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        max_nfev=160,
        xtol=1e-9,
        ftol=1e-9,
        gtol=1e-9,
        verbose=1,
    )
    if not result.success:
        raise RuntimeError(result.message)
    capacity, conductance, loss = decode_dynamic(result.x)
    maximum, minimum = fit_extreme_offsets(calibration)
    return FiveZoneParameters(capacity, conductance, loss, maximum, minimum)


def parameters_from_frozen() -> FiveZoneParameters:
    payload = json.loads(FROZEN.read_text(encoding="utf-8"))["parameters"]
    return FiveZoneParameters(
        np.asarray(payload["capacity_J_K"], dtype=float),
        np.asarray(payload["adjacent_conductance_W_K"], dtype=float),
        np.asarray(payload["ambient_loss_W_K"], dtype=float),
        np.asarray(payload["tmax_offset_C_intercept_and_current40"], dtype=float),
        np.asarray(payload["tmin_offset_C_intercept_and_current40"], dtype=float),
    )


def evaluate_rom(
    rom: PhysicsConstrainedFiveZoneROM, frames: list[pd.DataFrame]
) -> tuple[dict[str, object], pd.DataFrame]:
    rows: list[pd.DataFrame] = []
    elapsed = 0.0
    for frame in frames:
        target = frame[[f"T_zone_{index}_C" for index in range(1, 6)]].to_numpy(float)
        start = perf_counter()
        predicted = rom.simulate(
            frame.time_s.to_numpy(float),
            frame.current_A.to_numpy(float),
            frame.stack_voltage_V.to_numpy(float),
            target[0],
        )
        elapsed += perf_counter() - start
        maximum, minimum, delta = zip(
            *[
                rom.extreme_temperatures(state, current)
                for state, current in zip(
                    predicted, frame.current_A.to_numpy(float), strict=True
                )
            ],
            strict=True,
        )
        block = frame.copy()
        for index in range(5):
            block[f"ROM_T_zone_{index + 1}_C"] = predicted[:, index]
            block[f"zone_{index + 1}_error_C"] = predicted[:, index] - target[:, index]
        block["ROM_Tmax_C"] = maximum
        block["ROM_Tmin_C"] = minimum
        block["ROM_DeltaT_C"] = delta
        block["ROM_hotspot_zone_zero_based"] = np.argmax(predicted, axis=1)
        rows.append(block)
    result = pd.concat(rows, ignore_index=True)
    scoring = result.time_s > 0.0
    zone_error = result.loc[
        scoring, [f"zone_{index}_error_C" for index in range(1, 6)]
    ].to_numpy(float)
    tmax_error = result.loc[scoring, "ROM_Tmax_C"].to_numpy(float) - result.loc[
        scoring, "Tmax_C"
    ].to_numpy(float)
    delta_error = result.loc[scoring, "ROM_DeltaT_C"].to_numpy(float) - result.loc[
        scoring, "DeltaT_C"
    ].to_numpy(float)
    mean_error = result.loc[
        scoring, [f"ROM_T_zone_{index}_C" for index in range(1, 6)]
    ].to_numpy(float).mean(axis=1) - result.loc[
        scoring, [f"T_zone_{index}_C" for index in range(1, 6)]
    ].to_numpy(float).mean(axis=1)
    terminal_biases = []
    for _, frame in result.loc[scoring].groupby("scenario_id"):
        count = max(int(np.ceil(len(frame) / 4)), 1)
        errors = frame[[f"zone_{index}_error_C" for index in range(1, 6)]].tail(count)
        terminal_biases.append(float(np.max(np.abs(errors.mean(axis=0)))))
    energy_residual = 0.0
    for _, row in result.loc[scoring].iterrows():
        state = row[[f"ROM_T_zone_{index}_C" for index in range(1, 6)]].to_numpy(float)
        audit = rom.instantaneous_energy_audit(
            state, float(row.current_A), float(row.stack_voltage_V)
        )
        energy_residual = max(energy_residual, audit["energy_residual_percent"])
    metrics = {
        "scenario_ids": sorted(result.scenario_id.unique().tolist()),
        "zone_overall_RMSE_C": float(np.sqrt(np.mean(np.square(zone_error)))),
        "each_zone_RMSE_C": [
            float(value) for value in np.sqrt(np.mean(np.square(zone_error), axis=0))
        ],
        "Tmax_RMSE_C": float(np.sqrt(np.mean(np.square(tmax_error)))),
        "DeltaT_RMSE_C": float(np.sqrt(np.mean(np.square(delta_error)))),
        "mean_temperature_RMSE_C": float(np.sqrt(np.mean(np.square(mean_error)))),
        "hotspot_zone_accuracy_percent": float(
            np.mean(
                result.loc[scoring, "ROM_hotspot_zone_zero_based"].to_numpy(int)
                == result.loc[scoring, "hotspot_zone_zero_based"].to_numpy(int)
            )
            * 100.0
        ),
        "maximum_terminal_quarter_bias_C": float(max(terminal_biases)),
        "maximum_energy_residual_percent": float(energy_residual),
        "runtime_s": float(elapsed),
        "mean_runtime_per_step_s": float(elapsed / max(scoring.sum(), 1)),
        "maximum_real_stability_eigenvalue_per_s": float(
            np.max(np.real(rom.stability_eigenvalues()))
        ),
    }
    return metrics, result


def rom_checks(metrics: dict[str, object], protocol: dict[str, object]) -> dict[str, bool]:
    threshold = protocol["ROM_acceptance_thresholds"]
    metric_names = {
        "zone_overall_RMSE_C": "zone_overall_RMSE_C",
        "Tmax_RMSE_C": "Tmax_RMSE_C",
        "DeltaT_RMSE_C": "DeltaT_RMSE_C",
        "hotspot_zone_accuracy_percent": "hotspot_zone_accuracy_percent",
        "mean_temperature_RMSE_C": "mean_temperature_RMSE_C",
        "energy_residual_percent": "maximum_energy_residual_percent",
        "maximum_terminal_quarter_bias_C": "maximum_terminal_quarter_bias_C",
    }
    checks: dict[str, bool] = {}
    for threshold_name, value in threshold.items():
        metric = float(metrics[metric_names[threshold_name]])
        if threshold_name == "hotspot_zone_accuracy_percent":
            checks[threshold_name] = metric >= float(value)
        else:
            checks[threshold_name] = metric <= float(value)
    return checks


def observer_metrics(
    truth: np.ndarray,
    estimate: np.ndarray,
    true_tmax: np.ndarray,
    current: np.ndarray,
    rom: PhysicsConstrainedFiveZoneROM,
) -> dict[str, float]:
    score = np.arange(1, len(truth))
    error = estimate[score] - truth[score]
    estimated_tmax = np.asarray(
        [
            rom.extreme_temperatures(state, current[index])[0]
            for index, state in enumerate(estimate)
        ]
    )
    tmax_error = estimated_tmax[score] - true_tmax[score]
    return {
        "zone_RMSE_C": float(np.sqrt(np.mean(np.square(error)))),
        "Tmax_RMSE_C": float(np.sqrt(np.mean(np.square(tmax_error)))),
        "maximum_absolute_error_C": float(np.max(np.abs(error))),
        "q95_absolute_error_C": float(np.quantile(np.abs(error), 0.95)),
        "hotspot_zone_accuracy_percent": float(
            np.mean(np.argmax(estimate[score], axis=1) == np.argmax(truth[score], axis=1))
            * 100.0
        ),
    }


def build_model_bank(rom: PhysicsConstrainedFiveZoneROM) -> list[PhysicsConstrainedFiveZoneROM]:
    p = rom.parameters
    variants = [
        p,
        replace(p, adjacent_conductance_W_K=p.adjacent_conductance_W_K * 0.75),
        replace(p, adjacent_conductance_W_K=p.adjacent_conductance_W_K * 1.25),
        replace(p, ambient_loss_W_K=p.ambient_loss_W_K * 0.85),
        replace(p, ambient_loss_W_K=p.ambient_loss_W_K * 1.15),
    ]
    return [PhysicsConstrainedFiveZoneROM(parameters) for parameters in variants]


def run_gated_observer(
    rom: PhysicsConstrainedFiveZoneROM,
    truth: np.ndarray,
    current: np.ndarray,
    voltage: np.ndarray,
    times: np.ndarray,
    measurements: np.ndarray,
    sigma_C: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bank = build_model_bank(rom)
    count = len(bank)
    sensors = np.asarray((0, 4), dtype=int)
    observation = np.zeros((len(sensors), 5))
    observation[np.arange(len(sensors)), sensors] = 1.0
    states = np.repeat((truth[0] + np.linspace(0.25, -0.25, 5))[None, :], count, axis=0)
    covariance = np.repeat((np.eye(5) * 0.8)[None, :, :], count, axis=0)
    probability = np.r_[0.60, np.full(count - 1, 0.10)]
    process_noise = np.eye(5) * 0.004
    measurement_noise = np.eye(len(sensors)) * sigma_C**2
    single = np.empty_like(truth)
    gated = np.empty_like(truth)
    enabled = np.zeros(len(truth), dtype=bool)
    single[0] = states[0]
    gated[0] = states[0]
    sustained_evidence = 0.0
    for step, dt_s in enumerate(np.diff(times)):
        likelihood = np.empty(count)
        baseline_nis = 0.0
        for index, model in enumerate(bank):
            transition, forcing = model.exact_discrete_transition(
                current[step], voltage[step], float(dt_s)
            )
            states[index] = model.ambient_C + transition @ (
                states[index] - model.ambient_C
            ) + forcing
            covariance[index] = (
                transition @ covariance[index] @ transition.T + process_noise
            )
            innovation = measurements[step + 1] - observation @ states[index]
            innovation_covariance = (
                observation @ covariance[index] @ observation.T + measurement_noise
            )
            gain = covariance[index] @ observation.T @ np.linalg.inv(
                innovation_covariance
            )
            states[index] += gain @ innovation
            identity = np.eye(5)
            covariance[index] = (
                (identity - gain @ observation)
                @ covariance[index]
                @ (identity - gain @ observation).T
                + gain @ measurement_noise @ gain.T
            )
            nis = float(innovation @ np.linalg.solve(innovation_covariance, innovation))
            if index == 0:
                baseline_nis = nis
            _, logdet = np.linalg.slogdet(innovation_covariance)
            likelihood[index] = -0.5 * (
                nis + logdet + len(sensors) * np.log(2.0 * np.pi)
            )
        raw = np.log(np.maximum(probability, 0.005)) + likelihood
        raw -= raw.max()
        raw_probability = np.exp(raw)
        raw_probability /= raw_probability.sum()
        probability = 0.65 * probability + 0.35 * raw_probability
        probability = np.maximum(probability, 0.005)
        probability /= probability.sum()
        sustained_evidence = 0.90 * sustained_evidence + 0.10 * max(
            0.0, baseline_nis - len(sensors)
        )
        entropy = -float(np.sum(probability * np.log(np.maximum(probability, 1e-12))))
        normalized_entropy = entropy / np.log(count)
        gate = (
            sustained_evidence > 1.50
            and normalized_entropy < 0.85
            and probability.max() > 0.45
        )
        single[step + 1] = states[0]
        gated[step + 1] = probability @ states if gate else states[0]
        enabled[step + 1] = gate
    return single, gated, enabled


def observer_study(
    rom: PhysicsConstrainedFiveZoneROM, final_frame: pd.DataFrame
) -> tuple[dict[str, object], pd.DataFrame]:
    truth = final_frame[[f"T_zone_{index}_C" for index in range(1, 6)]].to_numpy(float)
    true_tmax = final_frame.Tmax_C.to_numpy(float)
    current = final_frame.current_A.to_numpy(float)
    voltage = final_frame.stack_voltage_V.to_numpy(float)
    times = final_frame.time_s.to_numpy(float)
    sensors = np.asarray((0, 4), dtype=int)
    rows: list[dict[str, object]] = []
    start = perf_counter()
    for seed in range(10):
        rng = np.random.default_rng(160100 + seed)
        measurements = truth[:, sensors] + rng.normal(0.0, 0.12, (len(truth), 2))
        single, gated, enabled = run_gated_observer(
            rom, truth, current, voltage, times, measurements, 0.12
        )
        single_metric = observer_metrics(truth, single, true_tmax, current, rom)
        gated_metric = observer_metrics(truth, gated, true_tmax, current, rom)
        improvement = 100.0 * (
            single_metric["Tmax_RMSE_C"] - gated_metric["Tmax_RMSE_C"]
        ) / max(single_metric["Tmax_RMSE_C"], 1e-12)
        rows.append(
            {
                "scenario": "V15_synthetic_final_parent",
                "seed": seed,
                **{f"single_{key}": value for key, value in single_metric.items()},
                **{f"gated_{key}": value for key, value in gated_metric.items()},
                "Tmax_RMSE_improvement_percent": improvement,
                "bank_enabled_fraction": float(np.mean(enabled)),
            }
        )

    nominal_truth = rom.simulate(times, current, voltage, truth[0])
    nominal_tmax = np.asarray(
        [rom.extreme_temperatures(state, load)[0] for state, load in zip(nominal_truth, current)]
    )
    for label, sigma, seed_offset in (
        ("nominal", 0.12, 170000),
        ("pure_measurement_noise", 0.30, 180000),
    ):
        for seed in range(10):
            rng = np.random.default_rng(seed_offset + seed)
            measurements = nominal_truth[:, sensors] + rng.normal(
                0.0, sigma, (len(nominal_truth), 2)
            )
            single, gated, enabled = run_gated_observer(
                rom, nominal_truth, current, voltage, times, measurements, sigma
            )
            single_metric = observer_metrics(
                nominal_truth, single, nominal_tmax, current, rom
            )
            gated_metric = observer_metrics(
                nominal_truth, gated, nominal_tmax, current, rom
            )
            loss = 100.0 * (
                gated_metric["Tmax_RMSE_C"] - single_metric["Tmax_RMSE_C"]
            ) / max(single_metric["Tmax_RMSE_C"], 1e-12)
            rows.append(
                {
                    "scenario": label,
                    "seed": seed,
                    **{f"single_{key}": value for key, value in single_metric.items()},
                    **{f"gated_{key}": value for key, value in gated_metric.items()},
                    "Tmax_RMSE_improvement_percent": -loss,
                    "bank_enabled_fraction": float(np.mean(enabled)),
                }
            )
    runtime = perf_counter() - start
    result = pd.DataFrame(rows)
    off_grid = result.loc[result.scenario == "V15_synthetic_final_parent"]
    improvements = off_grid.Tmax_RMSE_improvement_percent.to_numpy(float)
    bootstrap_rng = np.random.default_rng(160188)
    bootstrap = np.mean(
        improvements[
            bootstrap_rng.integers(0, len(improvements), (5000, len(improvements)))
        ],
        axis=1,
    )
    nominal_loss = -result.loc[
        result.scenario == "nominal", "Tmax_RMSE_improvement_percent"
    ].median()
    noise_loss = -result.loc[
        result.scenario == "pure_measurement_noise", "Tmax_RMSE_improvement_percent"
    ].median()
    single_hotspot = off_grid.single_hotspot_zone_accuracy_percent.mean()
    gated_hotspot = off_grid.gated_hotspot_zone_accuracy_percent.mean()
    checks = {
        "majority_improves": bool(np.mean(improvements > 0.0) > 0.5),
        "median_improvement_above_10pct": bool(np.median(improvements) > 10.0),
        "paired_CI_lower_above_zero": bool(np.quantile(bootstrap, 0.025) > 0.0),
        "nominal_and_pure_noise_loss_below_5pct": bool(max(nominal_loss, noise_loss) < 5.0),
        "hotspot_accuracy_not_decreased": bool(gated_hotspot >= single_hotspot),
    }
    summary = {
        "track": "synthetic method development; not an experimental confirmation",
        "fixed_sensor_zones_zero_based": [0, 4],
        "fixed_noise_seeds": list(range(10)),
        "median_Tmax_RMSE_improvement_percent": float(np.median(improvements)),
        "improved_seed_fraction": float(np.mean(improvements > 0.0)),
        "paired_mean_improvement_95pct_CI_percent": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "nominal_median_performance_loss_percent": float(nominal_loss),
        "pure_noise_median_performance_loss_percent": float(noise_loss),
        "single_hotspot_zone_accuracy_percent": float(single_hotspot),
        "gated_hotspot_zone_accuracy_percent": float(gated_hotspot),
        "runtime_s": float(runtime),
        "acceptance_checks": checks,
        "synthetic_numeric_gate_passed": bool(all(checks.values())),
        "core_innovation_authorized": False,
        "authorization_reason": "The parent experimental high-fidelity gate is not passed and no fresh observer confirmation experiment exists.",
    }
    return summary, result


def main() -> None:
    OUT.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    DOCS.mkdir(exist_ok=True)
    protocol = create_protocol()
    residual_diagnosis = diagnose_v15_experimental_residuals()
    write_json(OUT / "v15_dynamic_residual_diagnosis.json", residual_diagnosis)

    plant = Fixed60HighResolutionReference()
    plant.assert_in_scope(ambient_C=23.0, duty=0.60, current_A=40.0)
    definitions = scenario_definitions()
    frames: dict[str, pd.DataFrame] = {}
    registry_rows: list[dict[str, object]] = []
    for definition in definitions:
        path = DATA / f"{definition['scenario_id']}.csv"
        if path.exists():
            frame = pd.read_csv(path)
        else:
            frame = generate_parent_scenario(plant, definition)
            frame.to_csv(path, index=False)
        frames[str(definition["scenario_id"])] = frame
        registry_rows.append(
            {
                **definition,
                "data_file": str(path.relative_to(ROOT)).replace("\\", "/"),
                "data_file_sha256": sha256(path),
                "parent_model_sha256": V15_SHA256,
                "source_type": "deterministic simulation of frozen V15 parent",
                "eligible_as_fresh_experimental_confirmation": False,
            }
        )
    registry = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "fresh_experimental_confirmation_available": False,
        "scenarios": registry_rows,
    }
    write_json(OUT / "high_resolution_dataset_registry.json", registry)

    calibration = [
        frames["calibration_staircase"],
        frames["calibration_triangular"],
    ]
    if FROZEN.exists():
        parameters = parameters_from_frozen()
        fit_action = "loaded existing freeze; no optimization"
    else:
        parameters = fit_rom(calibration)
        payload = {
            "schema_version": 1,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "parent_V15_frozen_sha256": V15_SHA256,
            "protocol_sha256": sha256(OUT / "preregistered_v16_protocol.json"),
            "model_source_sha256": sha256(
                ROOT / "models_v16/physics_constrained_five_zone_rom.py"
            ),
            "calibration_scenarios": [
                "calibration_staircase",
                "calibration_triangular",
            ],
            "parameters": parameters.to_dict(),
            "no_retuning_after_freeze": True,
            "claim_scope": "V15 fixed-duty synthetic method-development ROM only",
        }
        write_json(FROZEN, payload)
        fit_action = "calibration-only optimization then frozen"
    frozen_hash = sha256(FROZEN)
    rom = PhysicsConstrainedFiveZoneROM(parameters)
    rom.assert_physical()

    calibration_metrics, calibration_output = evaluate_rom(rom, calibration)
    validation_metrics, validation_output = evaluate_rom(
        rom, [frames["validation_multisine"]]
    )
    first_final = not FINAL_MARKER.exists()
    final_metrics, final_output = evaluate_rom(rom, [frames["final_prbs_seed1601"]])
    combined = pd.concat(
        (calibration_output, validation_output, final_output), ignore_index=True
    )
    combined.to_csv(OUT / "five_zone_rom_predictions.csv", index=False)
    checks = {
        "calibration": rom_checks(calibration_metrics, protocol),
        "validation": rom_checks(validation_metrics, protocol),
        "synthetic_final_confirmation": rom_checks(final_metrics, protocol),
    }
    method_gate = all(checks["validation"].values()) and all(
        checks["synthetic_final_confirmation"].values()
    )
    rom_summary = {
        "fit_action": fit_action,
        "frozen_ROM_sha256": frozen_hash,
        "first_synthetic_final_confirmation_run": first_final,
        "parent_V15_experimental_gate_passed": False,
        "calibration": calibration_metrics,
        "validation": validation_metrics,
        "synthetic_final_confirmation": final_metrics,
        "acceptance_checks": checks,
        "synthetic_method_gate_passed": bool(method_gate),
        "independent_experimental_ROM_validation_passed": False,
        "formal_ROM_claim_authorized": False,
        "authorization_reason": "No unused independent dynamic experiment exists and the V15 strict experimental gate remains failed.",
    }
    write_json(OUT / "five_zone_rom_validation_metrics.json", rom_summary)
    if first_final:
        write_json(
            FINAL_MARKER,
            {
                "executed_utc": datetime.now(timezone.utc).isoformat(),
                "frozen_ROM_sha256": frozen_hash,
                "result": "pass" if method_gate else "fail",
                "data_type": "synthetic V15 parent scenario, not experiment",
                "no_retuning_after_marker": True,
            },
        )

    observer_summary, observer_rows = observer_study(
        rom, frames["final_prbs_seed1601"]
    )
    observer_rows.to_csv(OUT / "observer_method_study.csv", index=False)
    write_json(OUT / "observer_method_study_summary.json", observer_summary)

    strict_parent_gate = False
    control_authorized = bool(
        strict_parent_gate
        and rom_summary["independent_experimental_ROM_validation_passed"]
        and observer_summary["core_innovation_authorized"]
    )
    final_summary = {
        "branch_stage": "V16 physics-ROM-observer programme",
        "protected_parent": {
            "tag": "v15_fixed60_reference_confirmed",
            "frozen_reference_sha256": V15_SHA256,
            "V15_outputs_overwritten": False,
        },
        "step_1_scope_freeze": "complete",
        "step_2_preregistered_data_registry": {
            "complete": True,
            "fresh_independent_experiment_available": False,
        },
        "step_3_dynamic_discrepancy_diagnosis": {
            "complete": True,
            "primary_finding": residual_diagnosis["interpretation"][2],
            "worst_error": residual_diagnosis["worst_error"],
        },
        "step_4_high_resolution_confirmation": {
            "COMSOL_implementation_consistency_passed": True,
            "strict_experimental_gate_passed": False,
            "reason": "V15 maximum error 4.271 C exceeds 4.0 C and regional order 94.737% is below 95%; no fresh confirmation experiment is available.",
        },
        "step_5_five_zone_ROM": {
            "synthetic_method_study_completed": True,
            "synthetic_method_gate_passed": bool(method_gate),
            "independent_experimental_validation_passed": False,
            "formal_claim_authorized": False,
        },
        "step_6_sparse_observer": {
            "synthetic_method_study_completed": True,
            "synthetic_numeric_gate_passed": observer_summary[
                "synthetic_numeric_gate_passed"
            ],
            "preregistered_experimental_confirmation_passed": False,
            "core_innovation_authorized": False,
        },
        "step_7_control_decision": {
            "MPC_run": False,
            "control_research_authorized": control_authorized,
            "decision": "Do not run decisive MPC; variable-PWM high-fidelity plant is unavailable.",
        },
        "thesis_decision": {
            "recommended_core": [
                "physics-constrained energy-conserving fixed-duty five-zone ROM methodology",
                "structural-uncertainty-aware sparse thermal estimation, pending fresh confirmation",
            ],
            "claims_allowed_now": [
                "experimentally constrained fixed-60%-PWM high-resolution thermal reference",
                "COMSOL-Python implementation consistency",
                "synthetic V15-parent ROM/observer method results with explicit limitations",
            ],
            "claims_not_allowed": [
                "variable-PWM high-fidelity digital twin",
                "independently experimentally validated five-zone ROM",
                "experimentally confirmed gated MM-EKF core innovation",
                "MPC superiority",
            ],
            "next_external_action": "Acquire a newly reserved multi-point fixed-60% dynamic experiment; add fan speed/flow or pressure telemetry before any variable-PWM claim.",
        },
        "reference_temperature_C": 55.0,
        "reference_temperature_is_not_a_safety_limit": True,
        "draft_PR_only": True,
        "merge_to_main_allowed": False,
    }
    write_json(OUT / "final_summary.json", final_summary)
    print(json.dumps(final_summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
