"""V6 gated MM-EKF audit on the frozen V5 objects and noise seeds.

The V5 CSV is the test-set registry: object scales and object/seed identifiers
are read from it rather than sampled again.  The nominal model filter is always
maintained and the bank estimate is used only after sustained innovation and
probability-concentration evidence.  This is a post-hoc validation experiment,
not a new sensor-layout or model-bank design exercise.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from experiments.run_observer_decisive_validation import (  # noqa: E402
    metric,
    run_kf,
    simulate_truth,
    system_matrices,
)
from models.five_node_model import FiveNodeParameters, interpolate_three_to_five  # noqa: E402
from models.model_registry import load_primary_model, write_model_runtime_record  # noqa: E402
from models.three_node_model import load_frozen_parameters, simulate  # noqa: E402


def run_gated_bank(
    base, bank_scales, sensors, measurements, current, duty, direction,
    dt, heat_model, initial, sigma,
):
    """Run a rate-limited bank with innovation/entropy fallback to model zero."""
    count = len(bank_scales)
    states = np.repeat(initial[None, :], count, axis=0)
    covariance = np.repeat((np.eye(5) * 0.8)[None, :, :], count, axis=0)
    probability = np.r_[0.60, np.full(count - 1, 0.40 / (count - 1))]
    process_noise = np.eye(5) * 0.004
    measurement_noise = np.eye(len(sensors)) * sigma**2
    observation = np.zeros((len(sensors), 5))
    observation[np.arange(len(sensors)), sensors] = 1.0
    estimate = np.empty((len(current), 5))
    probability_history = np.empty((len(current), count))
    bank_enabled = np.zeros(len(current), dtype=bool)
    evidence_history = np.zeros(len(current))
    estimate[0] = initial
    probability_history[0] = probability
    sustained_evidence = 0.0
    for k in range(len(current) - 1):
        log_likelihood = np.empty(count)
        baseline_nis = 0.0
        for j, scales in enumerate(bank_scales):
            transition, forcing = system_matrices(
                base, scales, current[k], duty[k], int(direction[k]), dt, heat_model
            )
            states[j] = transition @ states[j] + forcing
            covariance[j] = transition @ covariance[j] @ transition.T + process_noise
            innovation = measurements[k + 1] - observation @ states[j]
            innovation_covariance = (
                observation @ covariance[j] @ observation.T + measurement_noise
            )
            gain = covariance[j] @ observation.T @ np.linalg.inv(innovation_covariance)
            states[j] += gain @ innovation
            covariance[j] = (np.eye(5) - gain @ observation) @ covariance[j]
            nis = float(innovation @ np.linalg.solve(innovation_covariance, innovation))
            if j == 0:
                baseline_nis = nis
            _, logdet = np.linalg.slogdet(innovation_covariance)
            log_likelihood[j] = -0.5 * (
                nis + logdet + len(sensors) * np.log(2.0 * np.pi)
            )

        raw_log_weight = np.log(np.maximum(probability, 0.005)) + log_likelihood
        raw_log_weight -= raw_log_weight.max()
        raw_probability = np.exp(raw_log_weight)
        raw_probability /= raw_probability.sum()
        # Limit a one-step probability collapse and preserve a meaningful floor.
        probability = 0.65 * probability + 0.35 * raw_probability
        probability = np.maximum(probability, 0.005)
        probability /= probability.sum()

        # For two sensors E[NIS]=2 under the nominal hypothesis.  Only persistent
        # excess innovation is evidence for activating the model mixture.
        excess_nis = max(0.0, baseline_nis - len(sensors))
        sustained_evidence = 0.90 * sustained_evidence + 0.10 * excess_nis
        entropy = -float(np.sum(probability * np.log(np.maximum(probability, 1e-12))))
        normalized_entropy = entropy / np.log(count)
        concentrated = normalized_entropy < 0.85 and probability.max() > 0.45
        enabled = sustained_evidence > 1.50 and concentrated
        estimate[k + 1] = probability @ states if enabled else states[0]
        probability_history[k + 1] = probability
        bank_enabled[k + 1] = enabled
        evidence_history[k + 1] = sustained_evidence
    return estimate, probability_history, bank_enabled, evidence_history


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output_root = (args.output_dir or root / "outputs_v6").resolve()
    output = output_root / "observer_gated"
    output.mkdir(parents=True, exist_ok=True)

    primary, primary_path, model_hash = load_primary_model(root)
    write_model_runtime_record(root, output, "observer_gated_v6")
    config = yaml.safe_load(
        (root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8")
    )
    settings = config["observer_v5"]
    seed = int(config["project"]["seed"])
    heat_model = str(primary["heat_model"])
    thermal = load_frozen_parameters(primary_path)
    coolant = np.asarray(config["model"]["coolant_node_temperatures_forward_C"], float)
    base = FiveNodeParameters.from_three_node(thermal, coolant)

    digitized = root / "outputs_v4" / "digitization" / "regenerated"
    fig15 = pd.read_csv(digitized / "Fig15_digitized.csv")
    fig16 = pd.read_csv(digitized / "Fig16_digitized.csv")
    current_1s = np.interp(
        fig16.time_s, fig15.time_s, fig15.paper_simulation_current_A
    )
    initial_three = fig16.loc[
        0, ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]
    ].to_numpy(float)
    three = simulate(
        thermal, current_1s, initial_three, dt_s=1, bidirectional=True,
        initial_airflow_direction=-1, fan_initially_enabled=True,
        fan_enable_temperature_C=53, reference_temperature_C=55,
        smc_deadband_C=0.5, smc_update_period_s=10,
        heat_model=heat_model, coolant_node_temperatures_forward_C=coolant,
    )
    stride = int(settings["evaluation_dt_s"])
    dt = float(stride)
    current = current_1s[::stride]
    duty = three.fan_duty[::stride]
    direction = three.airflow_direction[::stride]
    initial = interpolate_three_to_five(initial_three)
    sensors = (0, 4)  # frozen V5 design-only selection
    sigma = 0.12
    bank = np.array([
        [1, 1, 1, 1, 0], [.85, 1, 1, 1, 0], [1.15, 1, 1, 1, 0],
        [1, .75, 1, 1, 0], [1, 1, 1, .85, 0],
    ], float)

    registry = pd.read_csv(
        root / "outputs_v5" / "observer_decisive" / "continuous_off_grid_100x10.csv"
    )
    object_columns = [
        "C_th_scale", "K_cool_scale", "K_node_scale",
        "load_input_scale", "coolant_offset_C",
    ]
    objects = registry.groupby("object_index", sort=True)[object_columns].first()
    expected_pairs = registry[["object_index", "noise_seed"]].drop_duplicates()
    if len(objects) != 100 or len(expected_pairs) != 1000:
        raise RuntimeError("Frozen V5 registry must contain 100 objects x 10 seeds")

    rows: list[dict[str, float | int]] = []
    probability_rows: list[dict[str, float | int]] = []
    started = time.perf_counter()
    for object_index, object_row in objects.iterrows():
        scales = object_row.to_numpy(float)
        truth = simulate_truth(
            base, scales, current, duty, direction, dt, heat_model, initial
        )
        seeds = expected_pairs.loc[
            expected_pairs.object_index == object_index, "noise_seed"
        ].sort_values()
        for noise_seed in seeds:
            rng = np.random.default_rng(
                seed + 100000 + int(object_index) * 100 + int(noise_seed)
            )
            measurements = truth[:, sensors] + rng.normal(
                0, sigma, (len(truth), len(sensors))
            )
            observer_initial = initial + np.linspace(0.25, -0.25, 5)
            single = run_kf(
                base, bank[0], 1, sensors, measurements, current, duty,
                direction, dt, heat_model, observer_initial, sigma,
            )
            gated, probability, enabled, evidence = run_gated_bank(
                base, bank, sensors, measurements, current, duty, direction,
                dt, heat_model, observer_initial, sigma,
            )
            single_metrics = metric(truth, single)
            gated_metrics = metric(truth, gated)
            improvement = 100.0 * (
                single_metrics["hotspot_rmse_C"] - gated_metrics["hotspot_rmse_C"]
            ) / single_metrics["hotspot_rmse_C"]
            rows.append({
                "object_index": int(object_index), "noise_seed": int(noise_seed),
                **{name: float(value) for name, value in zip(object_columns, scales)},
                **{f"single_{key}": value for key, value in single_metrics.items()},
                **{f"gated_{key}": value for key, value in gated_metrics.items()},
                "hotspot_improvement_percent": improvement,
            })
            probability_rows.append({
                "object_index": int(object_index), "noise_seed": int(noise_seed),
                "final_max_probability": float(probability[-1].max()),
                "mean_normalized_entropy": float(np.mean(
                    -np.sum(probability * np.log(np.maximum(probability, 1e-12)), axis=1)
                    / np.log(len(bank))
                )),
                "bank_enabled_fraction": float(np.mean(enabled)),
                "maximum_sustained_innovation_evidence": float(evidence.max()),
            })
    runtime_s = time.perf_counter() - started
    off_grid = pd.DataFrame(rows)
    off_grid.to_csv(output / "continuous_off_grid_100x10_gated.csv", index=False)
    pd.DataFrame(probability_rows).to_csv(
        output / "gating_probability_audit.csv", index=False
    )

    nominal_rows: list[dict[str, float | int | str]] = []
    nominal_truth = simulate_truth(
        base, bank[0], current, duty, direction, dt, heat_model, initial
    )
    for noise_sigma, label, offset in [
        (0.12, "nominal", 0), (0.30, "pure_measurement_noise", 100)
    ]:
        for noise_seed in range(20):
            rng = np.random.default_rng(seed + 500000 + noise_seed + offset)
            measurements = nominal_truth[:, sensors] + rng.normal(
                0, noise_sigma, (len(nominal_truth), len(sensors))
            )
            observer_initial = initial + np.linspace(0.25, -0.25, 5)
            single = run_kf(
                base, bank[0], 1, sensors, measurements, current, duty,
                direction, dt, heat_model, observer_initial, noise_sigma,
            )
            gated, _, enabled, _ = run_gated_bank(
                base, bank, sensors, measurements, current, duty, direction,
                dt, heat_model, observer_initial, noise_sigma,
            )
            single_rmse = metric(nominal_truth, single)["hotspot_rmse_C"]
            gated_rmse = metric(nominal_truth, gated)["hotspot_rmse_C"]
            nominal_rows.append({
                "scenario": label, "seed": noise_seed,
                "single_hotspot_rmse_C": single_rmse,
                "gated_hotspot_rmse_C": gated_rmse,
                "gated_loss_percent": 100.0 * (gated_rmse - single_rmse) / single_rmse,
                "bank_enabled_fraction": float(np.mean(enabled)),
            })
    nominal = pd.DataFrame(nominal_rows)
    nominal.to_csv(output / "nominal_and_noise_loss_gated.csv", index=False)

    improvements = off_grid.hotspot_improvement_percent.to_numpy()
    bootstrap_rng = np.random.default_rng(seed + 88)
    bootstrap = np.mean(
        improvements[bootstrap_rng.integers(0, len(improvements), (5000, len(improvements)))],
        axis=1,
    )
    ci = [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))]
    summary = {
        "primary_model_sha256": model_hash,
        "test_registry": "outputs_v5/observer_decisive/continuous_off_grid_100x10.csv",
        "fixed_object_count": int(len(objects)),
        "fixed_noise_seed_count_per_object": 10,
        "paired_sample_count": int(len(off_grid)),
        "median_hotspot_rmse_improvement_percent": float(np.median(improvements)),
        "mean_hotspot_rmse_improvement_percent": float(np.mean(improvements)),
        "improved_sample_fraction": float(np.mean(improvements > 0)),
        "paired_mean_improvement_95pct_CI_percent": ci,
        "nominal_median_performance_loss_percent": float(
            nominal.loc[nominal.scenario == "nominal", "gated_loss_percent"].median()
        ),
        "pure_noise_median_performance_loss_percent": float(
            nominal.loc[
                nominal.scenario == "pure_measurement_noise", "gated_loss_percent"
            ].median()
        ),
        "runtime_s": runtime_s,
        "gate_configuration": {
            "initial_nominal_probability": 0.60,
            "probability_floor": 0.005,
            "new_probability_fraction_per_step": 0.35,
            "normalized_entropy_maximum": 0.85,
            "maximum_probability_minimum": 0.45,
            "sustained_excess_NIS_threshold": 1.50,
            "nominal_and_noise_fallback": "single model zero when evidence gate is closed",
        },
        "post_hoc_test_set_tuning_warning": True,
    }
    summary["core_innovation_acceptance"] = {
        "off_grid_majority_improves": summary["improved_sample_fraction"] > 0.5,
        "paired_CI_lower_above_zero": ci[0] > 0,
        "median_improvement_above_10pct": summary["median_hotspot_rmse_improvement_percent"] > 10,
        "nominal_and_noise_loss_below_5pct": max(
            summary["nominal_median_performance_loss_percent"],
            summary["pure_noise_median_performance_loss_percent"],
        ) < 5,
    }
    summary["MM_EKF_core_innovation_supported"] = all(
        summary["core_innovation_acceptance"].values()
    )
    (output / "observer_gated_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
