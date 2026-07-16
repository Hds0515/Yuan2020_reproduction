"""Fair five-node sensor-layout and parameter-model-bank observer experiments."""

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


def _variant(
    base: FiveNodeParameters,
    *,
    capacity: float = 1.0,
    cooling: float = 1.0,
    conduction: float = 1.0,
) -> FiveNodeParameters:
    return replace(
        base,
        total_heat_capacity_J_K=base.total_heat_capacity_J_K * capacity,
        cooling_coefficient_per_node_W_K=base.cooling_coefficient_per_node_W_K * cooling,
        conduction_edge_W_K=base.conduction_edge_W_K * conduction,
    )


def _new_filter(
    parameters: FiveNodeParameters,
    sensors: tuple[int, ...],
    initial: np.ndarray,
    measurement_sigma_C: float,
    heat_model: str,
) -> ExtendedKalmanFilter:
    return ExtendedKalmanFilter(
        state_C=initial.copy(),
        covariance=np.eye(5) * 0.8,
        process_covariance=np.eye(5) * 0.0025,
        measurement_variance_C2=measurement_sigma_C**2,
        sensor_indices=sensors,
        parameters=parameters,
        heat_model=heat_model,
    )


def _run_filter(
    observer: ExtendedKalmanFilter,
    measurements: np.ndarray,
    current: np.ndarray,
    duty: np.ndarray,
    direction: np.ndarray,
    dt_s: float,
) -> np.ndarray:
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
    estimated_hotspot = np.max(estimate, axis=1)
    hotspot_error = estimated_hotspot - true_hotspot
    return {
        **{f"node_{i + 1}_rmse_C": float(value) for i, value in enumerate(node_rmse)},
        "mean_node_rmse_C": float(np.sqrt(np.mean(error**2))),
        "hotspot_rmse_C": float(np.sqrt(np.mean(hotspot_error**2))),
        "maximum_absolute_error_C": float(np.max(np.abs(error))),
        "q95_absolute_error_C": float(np.quantile(np.abs(error), 0.95)),
        "hotspot_position_accuracy": float(
            np.mean(np.argmax(truth, axis=1) == np.argmax(estimate, axis=1))
        ),
    }


def _paired_bootstrap_interval(
    improvements: np.ndarray, seed: int, replicate_count: int = 5000
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(improvements), size=(replicate_count, len(improvements)))
    means = np.mean(improvements[indices], axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output_root = (args.output_dir or root / "outputs_v4").resolve()
    output_five = output_root / "five_node"
    output_observer = output_root / "observer"
    output_five.mkdir(parents=True, exist_ok=True)
    output_observer.mkdir(parents=True, exist_ok=True)

    config = yaml.safe_load(
        (root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8")
    )
    seed = int(config["project"]["seed"])
    dt_s = float(config["project"]["dt_s"])
    frozen = json.loads(
        (root / "identification_v3" / "frozen_parameters_v3.json").read_text(encoding="utf-8")
    )
    heat_model = str(frozen["heat_model"])
    thermal = load_frozen_parameters(root / "identification_v3" / "frozen_parameters_v3.json")
    coolant = np.asarray(config["model"]["coolant_node_temperatures_forward_C"], dtype=float)
    five = FiveNodeParameters.from_three_node(thermal, coolant)

    digitized = output_root / "digitization" / "regenerated"
    if not digitized.exists():
        digitized = root / "digitization" / "regenerated"
    fig15 = pd.read_csv(digitized / "Fig15_digitized.csv")
    fig16 = pd.read_csv(digitized / "Fig16_digitized.csv")
    current = np.interp(fig16["time_s"], fig15["time_s"], fig15["paper_simulation_current_A"])
    initial_three = fig16.loc[
        0, ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]
    ].to_numpy(dtype=float)
    common = dict(
        dt_s=dt_s,
        bidirectional=True,
        initial_airflow_direction=-1,
        fan_initially_enabled=True,
        fan_enable_temperature_C=float(config["model"]["fan_enable_temperature_C"]),
        reference_temperature_C=float(config["model"]["reference_temperature_C"]),
        smc_deadband_C=float(config["model"]["smc_deadband_C"]),
        smc_update_period_s=float(config["model"]["smc_update_period_s"]),
        heat_model=heat_model,
        coolant_node_temperatures_forward_C=coolant,
    )
    three_result = simulate(thermal, current, initial_three, **common)

    combined = _variant(five, capacity=1.12, cooling=0.85, conduction=1.20)
    scenarios = {
        "nominal": (five, 0.12, 1.0),
        "capacity_plus_12pct": (_variant(five, capacity=1.12), 0.12, 1.0),
        "cooling_minus_15pct": (_variant(five, cooling=0.85), 0.12, 1.0),
        "conduction_plus_20pct": (_variant(five, conduction=1.20), 0.12, 1.0),
        "combined_mismatch": (combined, 0.12, 1.0),
        "load_input_minus_10pct": (five, 0.12, 0.90),
        "measurement_noise": (five, 0.30, 1.0),
    }
    truths: dict[str, np.ndarray] = {}
    measurements_full: dict[str, np.ndarray] = {}
    energy_residuals: dict[str, np.ndarray] = {}
    for scenario_index, (scenario, (plant, sigma, truth_input_scale)) in enumerate(scenarios.items()):
        truth, residual = simulate_open_loop(
            plant,
            current * truth_input_scale,
            three_result.fan_duty,
            three_result.airflow_direction,
            interpolate_three_to_five(initial_three),
            dt_s,
            heat_model,
        )
        rng = np.random.default_rng(seed + scenario_index)
        truths[scenario] = truth
        measurements_full[scenario] = truth + rng.normal(0.0, sigma, size=truth.shape)
        energy_residuals[scenario] = residual

    nominal_truth = truths["nominal"]
    consistency = pd.DataFrame(
        {
            "time_s": fig16["time_s"],
            "three_node_mean_C": np.mean(three_result.temperature_C, axis=1),
            "five_node_mean_C": np.mean(nominal_truth, axis=1),
            "three_node_hotspot_C": np.max(three_result.temperature_C, axis=1),
            "five_node_hotspot_C": np.max(nominal_truth, axis=1),
            "five_node_energy_residual_W": energy_residuals["nominal"],
        }
    )
    consistency.to_csv(output_five / "three_vs_five_node_consistency.csv", index=False)
    consistency_metrics = {
        "mean_temperature_rmse_C": float(
            np.sqrt(np.mean((consistency.three_node_mean_C - consistency.five_node_mean_C) ** 2))
        ),
        "hotspot_rmse_C": float(
            np.sqrt(np.mean((consistency.three_node_hotspot_C - consistency.five_node_hotspot_C) ** 2))
        ),
        "maximum_absolute_energy_residual_W": float(
            max(np.max(np.abs(value)) for value in energy_residuals.values())
        ),
        "five_node_total_capacity_J_K": five.total_heat_capacity_J_K,
        "area_preservation_ratio": float(
            5 * five.cooling_coefficient_per_node_W_K / (3 * thermal.K_cool_W_K)
        ),
    }
    (output_five / "consistency_metrics.json").write_text(
        json.dumps(consistency_metrics, indent=2) + "\n", encoding="utf-8"
    )

    layout_rows: list[dict[str, object]] = []
    trajectories: dict[tuple[str, tuple[int, ...]], np.ndarray] = {}
    for scenario, (_, sigma, _) in scenarios.items():
        truth = truths[scenario]
        for sensor_count in (1, 2, 3):
            for sensors in itertools.combinations(range(5), sensor_count):
                observer = _new_filter(
                    five,
                    sensors,
                    truth[0] + np.linspace(0.3, -0.3, 5),
                    sigma,
                    heat_model,
                )
                estimate = _run_filter(
                    observer,
                    measurements_full[scenario][:, sensors],
                    current,
                    three_result.fan_duty,
                    three_result.airflow_direction,
                    dt_s,
                )
                trajectories[(scenario, sensors)] = estimate
                layout_rows.append(
                    {
                        "scenario": scenario,
                        "sensor_count": sensor_count,
                        "sensor_nodes_1_based": ",".join(str(i + 1) for i in sensors),
                        **_metrics(truth, estimate),
                    }
                )
    layout_metrics = pd.DataFrame(layout_rows).sort_values(
        ["scenario", "sensor_count", "hotspot_rmse_C"]
    )
    layout_metrics.to_csv(output_observer / "sensor_layout_metrics.csv", index=False)
    layout_summary = (
        layout_metrics.groupby(["sensor_count", "sensor_nodes_1_based"], as_index=False)
        .agg(
            mean_hotspot_rmse_C=("hotspot_rmse_C", "mean"),
            worst_hotspot_rmse_C=("hotspot_rmse_C", "max"),
            mean_position_accuracy=("hotspot_position_accuracy", "mean"),
        )
        .sort_values(["sensor_count", "mean_hotspot_rmse_C"])
    )
    layout_summary.to_csv(output_observer / "sensor_layout_summary.csv", index=False)
    best_two_text = str(layout_summary[layout_summary.sensor_count == 2].iloc[0].sensor_nodes_1_based)
    best_two = tuple(int(value) - 1 for value in best_two_text.split(","))

    model_names = (
        "nominal",
        "capacity_plus_12pct",
        "cooling_minus_15pct",
        "conduction_plus_20pct",
        "combined_mismatch",
        "load_input_minus_10pct",
    )
    bank_models = [
        five,
        _variant(five, capacity=1.12),
        _variant(five, cooling=0.85),
        _variant(five, conduction=1.20),
        combined,
        five,
    ]
    input_scales = np.asarray([1.0, 1.0, 1.0, 1.0, 1.0, 0.90])
    method_rows: list[dict[str, object]] = []
    wide_rows: list[dict[str, object]] = []
    probability_frames: list[pd.DataFrame] = []
    bank_estimates: dict[str, np.ndarray] = {}
    for scenario, (_, sigma, _) in scenarios.items():
        truth = truths[scenario]
        measurements = measurements_full[scenario][:, best_two]
        initial = truth[0] + np.linspace(0.3, -0.3, 5)
        single = _new_filter(five, best_two, initial, sigma, heat_model)
        single_estimate = _run_filter(
            single,
            measurements,
            current,
            three_result.fan_duty,
            three_result.airflow_direction,
            dt_s,
        )
        filters = [
            _new_filter(model, best_two, initial, sigma, heat_model) for model in bank_models
        ]
        bank = MultipleModelObserver(
            filters=filters,
            probabilities=np.full(len(filters), 1.0 / len(filters)),
            model_names=model_names,
            input_scales=input_scales,
        )
        bank_estimate = np.empty_like(truth)
        probability = np.empty((len(truth), len(filters)))
        bank.update(measurements[0])
        bank_estimate[0] = bank.state_C
        probability[0] = bank.probabilities
        for k in range(len(current) - 1):
            bank.predict(
                current[k],
                three_result.fan_duty[k],
                int(three_result.airflow_direction[k]),
                dt_s,
            )
            bank.update(measurements[k + 1])
            bank_estimate[k + 1] = bank.state_C
            probability[k + 1] = bank.probabilities
        bank_estimates[scenario] = bank_estimate
        probability_frame = pd.DataFrame(probability, columns=[f"probability_{name}" for name in model_names])
        probability_frame.insert(0, "time_s", fig16["time_s"].to_numpy())
        probability_frame.insert(0, "scenario", scenario)
        probability_frames.append(probability_frame)

        single_metrics = _metrics(truth, single_estimate)
        bank_metrics = _metrics(truth, bank_estimate)
        for method, metrics in (
            ("single_model_ekf", single_metrics),
            ("parameter_model_bank_mmae", bank_metrics),
        ):
            method_rows.append({"scenario": scenario, "method": method, **metrics})
        improvement = 100.0 * (
            float(single_metrics["hotspot_rmse_C"]) - float(bank_metrics["hotspot_rmse_C"])
        ) / max(float(single_metrics["hotspot_rmse_C"]), 1e-12)
        wide_rows.append(
            {
                "scenario": scenario,
                **{f"single_{key}": value for key, value in single_metrics.items()},
                **{f"model_bank_{key}": value for key, value in bank_metrics.items()},
                "model_bank_hotspot_improvement_percent": improvement,
            }
        )

    method_comparison = pd.DataFrame(method_rows)
    method_comparison.to_csv(output_observer / "multiple_model_method_metrics.csv", index=False)
    comparison = pd.DataFrame(wide_rows)
    comparison.to_csv(output_observer / "multiple_model_comparison.csv", index=False)
    pd.concat(probability_frames, ignore_index=True).to_csv(
        output_observer / "model_probability_history.csv", index=False
    )

    mismatch_names = [
        "capacity_plus_12pct",
        "cooling_minus_15pct",
        "conduction_plus_20pct",
        "combined_mismatch",
        "load_input_minus_10pct",
    ]
    mismatch = comparison[comparison.scenario.isin(mismatch_names)]
    improvements = mismatch.model_bank_hotspot_improvement_percent.to_numpy(dtype=float)
    ci_lower, ci_upper = _paired_bootstrap_interval(improvements, seed + 991)
    improved_count = int(np.sum(improvements > 0.0))
    majority = improved_count > len(mismatch_names) / 2
    conclusion = (
        "The parameter-model-bank MM-EKF/MMAE improves hotspot RMSE in a majority of true mismatch scenarios."
        if majority
        else "当前多模型结构不构成明确创新优势 / the current multiple-model structure does not establish a clear innovation advantage."
    )
    summary = {
        "method_classification": "parameter-model-bank MM-EKF/MMAE; not an unknown-airflow-direction switching observer",
        "truth_generation_policy": "Each scenario uses its declared plant and load input to generate an independent truth trajectory.",
        "comparison_randomness_policy": "Within each scenario, single-model and model-bank observers use identical truth, measurements, and seed.",
        "declared_bank_models": [
            {"name": name, "input_scale": float(scale)}
            for name, scale in zip(model_names, input_scales, strict=True)
        ],
        "best_two_sensor_nodes_1_based": [i + 1 for i in best_two],
        "true_mismatch_scenario_count": len(mismatch_names),
        "model_bank_improved_scenario_count": improved_count,
        "model_bank_improves_majority_of_true_mismatches": majority,
        "mean_hotspot_improvement_percent": float(np.mean(improvements)),
        "worst_scenario_improvement_percent": float(np.min(improvements)),
        "paired_bootstrap_mean_improvement_95pct_CI_percent": [ci_lower, ci_upper],
        "paired_bootstrap_replicates": 5000,
        "conclusion": conclusion,
    }
    (output_observer / "observer_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    nominal_estimate = trajectories[("nominal", best_two)]
    figure, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    time = fig16["time_s"].to_numpy()
    axes[0].plot(time, np.max(nominal_truth, axis=1), label="true hotspot", lw=2)
    axes[0].plot(time, np.max(nominal_estimate, axis=1), "--", label="single EKF hotspot")
    axes[0].plot(time, np.max(bank_estimates["nominal"], axis=1), ":", label="model-bank hotspot")
    axes[0].set_ylabel("Temperature (°C)")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    axes[1].plot(time, np.argmax(nominal_truth, axis=1) + 1, label="true node")
    axes[1].plot(time, np.argmax(bank_estimates["nominal"], axis=1) + 1, "--", label="model-bank node")
    axes[1].set(xlabel="Time (s)", ylabel="Hotspot node")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    figure.suptitle(f"Five-node observer comparison, sensors at nodes {best_two_text}")
    figure.tight_layout()
    figure.savefig(output_observer / "ekf_hotspot_reconstruction.png", dpi=180)
    plt.close(figure)
    print(json.dumps({"consistency": consistency_metrics, "observer": summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
