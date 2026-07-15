"""Reproducible five-node consistency and nonlinear observer experiments."""

from __future__ import annotations

import argparse
from dataclasses import replace
import itertools
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from models.five_node_model import (  # noqa: E402
    FiveNodeParameters,
    interpolate_three_to_five,
    simulate_open_loop,
)
from models.three_node_model import load_frozen_parameters, simulate  # noqa: E402
from observers.ekf import ExtendedKalmanFilter  # noqa: E402
from observers.multiple_model_observer import MultipleModelObserver  # noqa: E402


def _variant(base: FiveNodeParameters, *, capacity: float = 1.0, cooling: float = 1.0,
             conduction: float = 1.0) -> FiveNodeParameters:
    return replace(
        base,
        total_heat_capacity_J_K=base.total_heat_capacity_J_K * capacity,
        cooling_coefficient_per_node_W_K=base.cooling_coefficient_per_node_W_K * cooling,
        conduction_edge_W_K=base.conduction_edge_W_K * conduction,
    )


def _new_filter(parameters: FiveNodeParameters, sensors: tuple[int, ...], initial: np.ndarray,
                measurement_sigma_C: float) -> ExtendedKalmanFilter:
    return ExtendedKalmanFilter(
        state_C=initial.copy(),
        covariance=np.eye(5) * 0.8,
        process_covariance=np.eye(5) * 0.0025,
        measurement_variance_C2=measurement_sigma_C**2,
        sensor_indices=sensors,
        parameters=parameters,
    )


def _run_filter(observer: ExtendedKalmanFilter, measurements: np.ndarray, current: np.ndarray,
                duty: np.ndarray, direction: np.ndarray, dt_s: float) -> np.ndarray:
    estimate = np.empty((len(current), 5))
    observer.update(measurements[0])
    estimate[0] = observer.state_C
    for k in range(len(current) - 1):
        observer.predict(current[k], duty[k], int(direction[k]), dt_s)
        observer.update(measurements[k + 1])
        estimate[k + 1] = observer.state_C
    return estimate


def _metrics(truth: np.ndarray, estimate: np.ndarray) -> dict[str, object]:
    error = estimate - truth
    node_rmse = np.sqrt(np.mean(error**2, axis=0))
    true_hotspot = np.max(truth, axis=1)
    estimate_hotspot = np.max(estimate, axis=1)
    hotspot_error = estimate_hotspot - true_hotspot
    position_accuracy = np.mean(np.argmax(truth, axis=1) == np.argmax(estimate, axis=1))
    return {
        **{f"node_{i + 1}_rmse_C": float(value) for i, value in enumerate(node_rmse)},
        "mean_node_rmse_C": float(np.sqrt(np.mean(error**2))),
        "hotspot_rmse_C": float(np.sqrt(np.mean(hotspot_error**2))),
        "hotspot_position_accuracy": float(position_accuracy),
        "maximum_absolute_error_C": float(np.max(np.abs(error))),
        "q95_absolute_error_C": float(np.quantile(np.abs(error), 0.95)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    root = args.root.resolve()
    output_five = root / "outputs" / "five_node"
    output_observer = root / "outputs" / "observer"
    output_five.mkdir(parents=True, exist_ok=True)
    output_observer.mkdir(parents=True, exist_ok=True)

    config = yaml.safe_load((root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    dt_s = float(config["project"]["dt_s"])
    heat_model = json.loads((root / "identification_v3" / "frozen_parameters_v3.json").read_text(encoding="utf-8"))["heat_model"]
    thermal = load_frozen_parameters(root / "identification_v3" / "frozen_parameters_v3.json")
    coolant = np.asarray(config["model"]["coolant_node_temperatures_forward_C"], dtype=float)
    five = FiveNodeParameters.from_three_node(thermal, coolant)

    fig15 = pd.read_csv(root / "digitization" / "regenerated" / "Fig15_digitized.csv")
    fig16 = pd.read_csv(root / "digitization" / "regenerated" / "Fig16_digitized.csv")
    current = np.interp(fig16["time_s"], fig15["time_s"], fig15["paper_simulation_current_A"])
    initial_three = fig16.loc[0, ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]].to_numpy(dtype=float)
    common = dict(
        dt_s=dt_s, bidirectional=True, initial_airflow_direction=-1, fan_initially_enabled=True,
        fan_enable_temperature_C=float(config["model"]["fan_enable_temperature_C"]),
        reference_temperature_C=float(config["model"]["reference_temperature_C"]),
        smc_deadband_C=float(config["model"]["smc_deadband_C"]),
        smc_update_period_s=float(config["model"]["smc_update_period_s"]),
        heat_model=heat_model, coolant_node_temperatures_forward_C=coolant,
    )
    three_result = simulate(thermal, current, initial_three, **common)
    truth, energy_residual = simulate_open_loop(
        five, current, three_result.fan_duty, three_result.airflow_direction,
        interpolate_three_to_five(initial_three), dt_s, heat_model,
    )

    consistency = pd.DataFrame({
        "time_s": fig16["time_s"],
        "three_node_mean_C": np.mean(three_result.temperature_C, axis=1),
        "five_node_mean_C": np.mean(truth, axis=1),
        "three_node_hotspot_C": np.max(three_result.temperature_C, axis=1),
        "five_node_hotspot_C": np.max(truth, axis=1),
        "five_node_energy_residual_W": energy_residual,
    })
    consistency.to_csv(output_five / "three_vs_five_node_consistency.csv", index=False)
    consistency_metrics = {
        "mean_temperature_rmse_C": float(np.sqrt(np.mean((consistency.three_node_mean_C - consistency.five_node_mean_C) ** 2))),
        "hotspot_rmse_C": float(np.sqrt(np.mean((consistency.three_node_hotspot_C - consistency.five_node_hotspot_C) ** 2))),
        "maximum_absolute_energy_residual_W": float(np.max(np.abs(energy_residual))),
        "five_node_total_capacity_J_K": five.total_heat_capacity_J_K,
        "area_preservation_ratio": float(5 * five.cooling_coefficient_per_node_W_K / (3 * thermal.K_cool_W_K)),
    }
    (output_five / "consistency_metrics.json").write_text(json.dumps(consistency_metrics, indent=2), encoding="utf-8")

    scenarios = {
        "nominal": (five, 0.12, 1.0),
        "measurement_noise": (five, 0.30, 1.0),
        "capacity_plus_12pct": (_variant(five, capacity=1.12), 0.12, 1.0),
        "cooling_minus_15pct": (_variant(five, cooling=0.85), 0.12, 1.0),
        "conduction_plus_20pct": (_variant(five, conduction=1.20), 0.12, 1.0),
        "load_preview_minus_10pct": (five, 0.12, 0.90),
    }
    rows: list[dict[str, object]] = []
    trajectory: dict[tuple[str, tuple[int, ...]], np.ndarray] = {}
    for scenario_index, (scenario, (model, sigma, current_scale)) in enumerate(scenarios.items()):
        rng = np.random.default_rng(20200715 + scenario_index)
        noise = rng.normal(0.0, sigma, size=truth.shape)
        for sensor_count in (1, 2, 3):
            for sensors in itertools.combinations(range(5), sensor_count):
                measurements = truth[:, sensors] + noise[:, sensors]
                observer = _new_filter(model, sensors, truth[0] + np.linspace(0.3, -0.3, 5), sigma)
                estimate = _run_filter(observer, measurements, current * current_scale,
                                       three_result.fan_duty, three_result.airflow_direction, dt_s)
                trajectory[(scenario, sensors)] = estimate
                rows.append({
                    "scenario": scenario,
                    "sensor_count": sensor_count,
                    "sensor_nodes_1_based": ",".join(str(i + 1) for i in sensors),
                    **_metrics(truth, estimate),
                })
    metrics = pd.DataFrame(rows).sort_values(["scenario", "sensor_count", "hotspot_rmse_C"])
    metrics.to_csv(output_observer / "sensor_layout_metrics.csv", index=False)
    summary = metrics.groupby(["sensor_count", "sensor_nodes_1_based"], as_index=False).agg(
        mean_hotspot_rmse_C=("hotspot_rmse_C", "mean"),
        worst_hotspot_rmse_C=("hotspot_rmse_C", "max"),
        mean_position_accuracy=("hotspot_position_accuracy", "mean"),
    ).sort_values(["sensor_count", "mean_hotspot_rmse_C"])
    summary.to_csv(output_observer / "sensor_layout_summary.csv", index=False)

    # A fair single-model versus model-bank comparison uses the best two-sensor
    # layout selected by aggregate error above, never by one test scenario.
    best_two_text = str(summary[summary.sensor_count == 2].iloc[0].sensor_nodes_1_based)
    best_two = tuple(int(value) - 1 for value in best_two_text.split(","))
    comparison_rows = []
    bank_models = [five, _variant(five, capacity=1.12), _variant(five, cooling=0.85), _variant(five, conduction=1.20)]
    for scenario_index, (scenario, (model, sigma, current_scale)) in enumerate(scenarios.items()):
        rng = np.random.default_rng(20200715 + scenario_index)
        measurements = truth[:, best_two] + rng.normal(0.0, sigma, size=(len(truth), len(best_two)))
        single = _new_filter(five, best_two, truth[0] + np.linspace(0.3, -0.3, 5), sigma)
        single_estimate = _run_filter(single, measurements, current * current_scale,
                                      three_result.fan_duty, three_result.airflow_direction, dt_s)
        filters = [_new_filter(item, best_two, truth[0] + np.linspace(0.3, -0.3, 5), sigma) for item in bank_models]
        bank = MultipleModelObserver(filters=filters, probabilities=np.full(len(filters), 1 / len(filters)))
        bank_estimate = np.empty_like(truth)
        bank.update(measurements[0]); bank_estimate[0] = bank.state_C
        for k in range(len(current) - 1):
            bank.predict(current[k] * current_scale, three_result.fan_duty[k], int(three_result.airflow_direction[k]), dt_s)
            bank.update(measurements[k + 1]); bank_estimate[k + 1] = bank.state_C
        for method, estimate in (("single_model_ekf", single_estimate), ("multiple_model_ekf", bank_estimate)):
            comparison_rows.append({"scenario": scenario, "method": method, **_metrics(truth, estimate)})
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(output_observer / "multiple_model_comparison.csv", index=False)

    nominal_estimate = trajectory[("nominal", best_two)]
    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    time = fig16["time_s"].to_numpy()
    axes[0].plot(time, np.max(truth, axis=1), label="true hotspot", lw=2)
    axes[0].plot(time, np.max(nominal_estimate, axis=1), "--", label="EKF hotspot")
    axes[0].set_ylabel("Temperature (°C)"); axes[0].legend(); axes[0].grid(alpha=0.25)
    axes[1].plot(time, np.argmax(truth, axis=1) + 1, label="true node")
    axes[1].plot(time, np.argmax(nominal_estimate, axis=1) + 1, "--", label="estimated node")
    axes[1].set(xlabel="Time (s)", ylabel="Hotspot node"); axes[1].legend(); axes[1].grid(alpha=0.25)
    figure.suptitle(f"Five-node EKF, sensors at nodes {best_two_text}")
    figure.tight_layout(); figure.savefig(output_observer / "ekf_hotspot_reconstruction.png", dpi=180); plt.close(figure)

    machine_summary = {
        "best_two_sensor_nodes_1_based": [i + 1 for i in best_two],
        "best_two_sensor_mean_hotspot_rmse_C": float(summary[(summary.sensor_count == 2) & (summary.sensor_nodes_1_based == best_two_text)].iloc[0].mean_hotspot_rmse_C),
        "multiple_model_mean_hotspot_rmse_C": float(comparison[comparison.method == "multiple_model_ekf"].hotspot_rmse_C.mean()),
        "single_model_mean_hotspot_rmse_C": float(comparison[comparison.method == "single_model_ekf"].hotspot_rmse_C.mean()),
        "multiple_model_improves_mean": bool(comparison[comparison.method == "multiple_model_ekf"].hotspot_rmse_C.mean() < comparison[comparison.method == "single_model_ekf"].hotspot_rmse_C.mean()),
    }
    (output_observer / "observer_summary.json").write_text(json.dumps(machine_summary, indent=2), encoding="utf-8")
    print(json.dumps({"consistency": consistency_metrics, "observer": machine_summary}, indent=2))


if __name__ == "__main__":
    main()
