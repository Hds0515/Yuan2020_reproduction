"""Fair closed-loop comparison of PI-SMC and EKF-based hotspot MPC."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from controllers.hotspot_mpc import HotspotMPC  # noqa: E402
from controllers.pi_smc import PISMCController  # noqa: E402
from models.five_node_model import FiveNodeParameters, interpolate_three_to_five, step_five_node  # noqa: E402
from models.three_node_model import load_frozen_parameters  # noqa: E402
from observers.ekf import ExtendedKalmanFilter  # noqa: E402


def _variant(base: FiveNodeParameters, capacity: float = 1.0, cooling: float = 1.0,
             conduction: float = 1.0) -> FiveNodeParameters:
    return replace(base,
        total_heat_capacity_J_K=base.total_heat_capacity_J_K * capacity,
        cooling_coefficient_per_node_W_K=base.cooling_coefficient_per_node_W_K * cooling,
        conduction_edge_W_K=base.conduction_edge_W_K * conduction)


def _simulate(controller_name: str, true_model: FiveNodeParameters, controller_model: FiveNodeParameters,
              current_true: np.ndarray, current_preview: np.ndarray, initial: np.ndarray,
              sensors: tuple[int, ...], noise_sigma_C: float, seed: int, dt_s: float,
              heat_model: str, reference_C: float, actuator_limit: float) -> tuple[pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(seed)
    temperature = np.empty((len(current_true), 5)); estimate = np.empty_like(temperature)
    duty = np.zeros(len(current_true)); direction = np.full(len(current_true), -1, dtype=int)
    objective = np.full(len(current_true), np.nan)
    energy_residual = np.zeros(len(current_true))
    temperature[0] = initial
    observer = ExtendedKalmanFilter(
        state_C=initial + np.linspace(0.25, -0.25, 5), covariance=np.eye(5) * 0.8,
        process_covariance=np.eye(5) * 0.0036, measurement_variance_C2=noise_sigma_C**2,
        sensor_indices=sensors, parameters=controller_model, heat_model=heat_model,
    )
    observer.update(temperature[0, sensors] + rng.normal(0, noise_sigma_C, len(sensors)))
    estimate[0] = observer.state_C
    if controller_name in ("Original-PI-SMC", "FiveNodeMean-PI-SMC"):
        controller: object = PISMCController(
            controller_model.controller_parameters.Kp_1_K,
            controller_model.controller_parameters.Ki_1_Ks,
            reference_C=reference_C,
            feedback_mode="middle" if controller_name == "Original-PI-SMC" else "mean",
            maximum_duty=actuator_limit,
        )
    else:
        controller = HotspotMPC(controller_model, heat_model, reference_C=reference_C, maximum_duty=actuator_limit)
    for k in range(len(current_true) - 1):
        if controller_name in ("Original-PI-SMC", "FiveNodeMean-PI-SMC"):
            command, flow = controller.command(estimate[k], k * dt_s, dt_s)  # type: ignore[attr-defined]
        else:
            stride = max(1, int(round(controller.prediction_dt_s / dt_s)))  # type: ignore[attr-defined]
            preview = current_preview[k::stride][:controller.horizon_steps]  # type: ignore[attr-defined]
            command, flow, objective[k] = controller.command(estimate[k], preview)  # type: ignore[attr-defined]
        command = min(command, actuator_limit)
        result = step_five_node(temperature[k], current_true[k], command, flow, dt_s, true_model, heat_model)
        temperature[k + 1] = result.next_temperature_C
        energy_residual[k] = result.energy_residual_W
        duty[k] = command; direction[k] = flow
        observer.predict(current_preview[k], command, flow, dt_s)
        measurement = temperature[k + 1, sensors] + rng.normal(0, noise_sigma_C, len(sensors))
        observer.update(measurement); estimate[k + 1] = observer.state_C
    duty[-1] = duty[-2]; direction[-1] = direction[-2]; objective[-1] = objective[-2]; energy_residual[-1] = energy_residual[-2]
    hotspot = np.max(temperature, axis=1); gradient = np.max(temperature, axis=1) - np.min(temperature, axis=1)
    mean_temperature = np.mean(temperature, axis=1)
    dt_hours = dt_s / 3600.0
    switches = int(np.sum(direction[1:] != direction[:-1]))
    metrics = {
        "controller": controller_name,
        "maximum_hotspot_C": float(np.max(hotspot)),
        "maximum_hotspot_after_60s_C": float(np.max(hotspot[min(len(hotspot) - 1, int(round(60 / dt_s))):])),
        "hotspot_rmse_above_reference_C": float(np.sqrt(np.mean(np.maximum(hotspot - reference_C, 0.0) ** 2))),
        "time_above_reference_s": float(np.sum(hotspot > reference_C) * dt_s),
        "maximum_gradient_C": float(np.max(gradient)),
        "mean_gradient_C": float(np.mean(gradient)),
        "gradient_rmse_C": float(np.sqrt(np.mean(gradient**2))),
        "mean_temperature_tracking_rmse_C": float(np.sqrt(np.mean((mean_temperature - reference_C) ** 2))),
        "fan_energy_proxy_duty2_h": float(np.sum(duty**2) * dt_hours),
        "mean_absolute_duty_rate_per_s": float(np.mean(np.abs(np.diff(duty))) / dt_s),
        "airflow_direction_switches": switches,
        "mean_estimated_hotspot_error_C": float(np.mean(np.max(estimate, axis=1) - hotspot)),
        "hotspot_estimation_rmse_C": float(np.sqrt(np.mean((np.max(estimate, axis=1) - hotspot) ** 2))),
        "maximum_absolute_energy_residual_W": float(np.max(np.abs(energy_residual))),
        "actuator_saturation_fraction": float(np.mean(duty >= actuator_limit - 1e-12)),
        "thermal_constraint_infeasible_fraction": float(np.mean((hotspot > reference_C) & (duty >= actuator_limit - 1e-12))),
    }
    data = {"time_s": np.arange(len(current_true)) * dt_s, "current_A": current_true,
            "preview_current_A": current_preview, "fan_duty": duty, "direction": direction,
            "hotspot_C": hotspot, "gradient_C": gradient, "objective": objective}
    for i in range(5):
        data[f"true_T{i + 1}_C"] = temperature[:, i]
        data[f"estimated_T{i + 1}_C"] = estimate[:, i]
    return pd.DataFrame(data), metrics


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args(); root = args.root.resolve()
    output = root / "outputs" / "control"; output.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load((root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    frozen = json.loads((root / "identification_v3" / "frozen_parameters_v3.json").read_text(encoding="utf-8"))
    thermal = load_frozen_parameters(root / "identification_v3" / "frozen_parameters_v3.json")
    heat_model = frozen["heat_model"]; dt_s = float(config["project"]["dt_s"])
    five = FiveNodeParameters.from_three_node(thermal, np.asarray(config["model"]["coolant_node_temperatures_forward_C"]))
    fig15 = pd.read_csv(root / "digitization" / "regenerated" / "Fig15_digitized.csv")
    fig16 = pd.read_csv(root / "digitization" / "regenerated" / "Fig16_digitized.csv")
    base_current = np.interp(fig16.time_s, fig15.time_s, fig15.paper_simulation_current_A)
    initial_three = fig16.loc[0, ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]].to_numpy(dtype=float)
    initial = interpolate_three_to_five(initial_three)
    hot_ambient = replace(five, coolant_forward_C=five.coolant_forward_C + 5.0)
    scenarios = {
        "nominal": dict(true=five, current=base_current, preview=base_current, noise=0.12, limit=1.0),
        "critical_load_plus_25pct": dict(true=five, current=np.minimum(base_current * 1.25, 42.0), preview=np.minimum(base_current * 1.25, 42.0), noise=0.12, limit=1.0),
        "plant_mismatch": dict(true=_variant(five, capacity=0.88, cooling=0.82, conduction=1.20), current=base_current, preview=base_current, noise=0.12, limit=1.0),
        "measurement_noise": dict(true=five, current=base_current, preview=base_current, noise=0.30, limit=1.0),
        "load_preview_minus_10pct": dict(true=five, current=base_current, preview=base_current * 0.90, noise=0.12, limit=1.0),
        "ambient_plus_5C": dict(true=hot_ambient, current=base_current, preview=base_current, noise=0.12, limit=1.0),
        "actuator_limited": dict(true=_variant(five, cooling=0.65), current=np.minimum(base_current * 1.75, 42.0), preview=np.minimum(base_current * 1.75, 42.0), noise=0.12, limit=0.15),
    }
    rows = []; trajectories: dict[tuple[str, str], pd.DataFrame] = {}
    sensors = (1, 4)
    for index, (scenario, values) in enumerate(scenarios.items()):
        for controller_name in ("Original-PI-SMC", "FiveNodeMean-PI-SMC", "Hotspot-MPC"):
            frame, metrics = _simulate(
                controller_name, values["true"], five, values["current"], values["preview"], initial,
                sensors, values["noise"], 20200715 + index, dt_s, heat_model,
                float(config["model"]["reference_temperature_C"]), values["limit"],
            )
            frame.to_csv(output / f"{scenario}_{controller_name.lower().replace('-', '_')}.csv", index=False)
            trajectories[(scenario, controller_name)] = frame
            rows.append({"scenario": scenario, **metrics})
    metrics_frame = pd.DataFrame(rows)
    metrics_frame.to_csv(output / "controller_comparison_metrics.csv", index=False)
    aggregate = metrics_frame.groupby("controller", as_index=False).agg(
        mean_maximum_hotspot_C=("maximum_hotspot_C", "mean"),
        mean_maximum_hotspot_after_60s_C=("maximum_hotspot_after_60s_C", "mean"),
        worst_maximum_hotspot_C=("maximum_hotspot_C", "max"),
        mean_hotspot_violation_rmse_C=("hotspot_rmse_above_reference_C", "mean"),
        mean_maximum_gradient_C=("maximum_gradient_C", "mean"),
        mean_gradient_C=("mean_gradient_C", "mean"),
        mean_gradient_rmse_C=("gradient_rmse_C", "mean"),
        mean_temperature_tracking_rmse_C=("mean_temperature_tracking_rmse_C", "mean"),
        mean_fan_energy_proxy_duty2_h=("fan_energy_proxy_duty2_h", "mean"),
        total_direction_switches=("airflow_direction_switches", "sum"),
        mean_hotspot_estimation_rmse_C=("hotspot_estimation_rmse_C", "mean"),
    )
    aggregate.to_csv(output / "controller_aggregate_metrics.csv", index=False)

    figure, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for name, style in (("Original-PI-SMC", "-"), ("FiveNodeMean-PI-SMC", ":"), ("Hotspot-MPC", "--")):
        frame = trajectories[("critical_load_plus_25pct", name)]
        axes[0].plot(frame.time_s, frame.hotspot_C, style, label=name)
        axes[1].plot(frame.time_s, frame.gradient_C, style, label=name)
        axes[2].plot(frame.time_s, frame.fan_duty, style, label=name)
    axes[0].axhline(55.0, color="k", lw=1, alpha=0.5); axes[0].set_ylabel("Hotspot (°C)")
    axes[1].set_ylabel("Gradient (°C)"); axes[2].set(xlabel="Time (s)", ylabel="Fan duty")
    for axis in axes: axis.grid(alpha=0.25); axis.legend()
    figure.suptitle("Critical-load closed-loop comparison (EKF feedback)")
    figure.tight_layout(); figure.savefig(output / "critical_load_comparison.png", dpi=180); plt.close(figure)

    mpc = aggregate[aggregate.controller == "Hotspot-MPC"].iloc[0]
    pi = aggregate[aggregate.controller == "Original-PI-SMC"].iloc[0]
    summary = {
        "feedback_source": "two-sensor EKF estimated five-node state",
        "sensor_nodes_1_based": [2, 5],
        "scenario_count": len(scenarios),
        "actuator_limited_declared_infeasible": bool(
            metrics_frame[(metrics_frame.scenario == "actuator_limited") &
                          (metrics_frame.controller == "Hotspot-MPC")].iloc[0].thermal_constraint_infeasible_fraction > 0
        ),
        "maximum_hotspot_metric_initial_condition_limited": bool(mpc.mean_maximum_hotspot_C == pi.mean_maximum_hotspot_C),
        "mpc_lower_mean_hotspot_violation_rmse": bool(mpc.mean_hotspot_violation_rmse_C < pi.mean_hotspot_violation_rmse_C),
        "mpc_lower_mean_gradient": bool(mpc.mean_gradient_C < pi.mean_gradient_C),
        "mpc_mean_maximum_hotspot_C": float(mpc.mean_maximum_hotspot_C),
        "pi_mean_maximum_hotspot_C": float(pi.mean_maximum_hotspot_C),
        "mpc_mean_hotspot_violation_rmse_C": float(mpc.mean_hotspot_violation_rmse_C),
        "pi_mean_hotspot_violation_rmse_C": float(pi.mean_hotspot_violation_rmse_C),
        "mpc_mean_gradient_C": float(mpc.mean_gradient_C),
        "pi_mean_gradient_C": float(pi.mean_gradient_C),
        "mpc_mean_fan_energy_proxy_duty2_h": float(mpc.mean_fan_energy_proxy_duty2_h),
        "pi_mean_fan_energy_proxy_duty2_h": float(pi.mean_fan_energy_proxy_duty2_h),
    }
    (output / "control_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
