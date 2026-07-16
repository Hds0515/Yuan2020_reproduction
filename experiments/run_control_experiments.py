"""Closed-loop comparison with explicit reference and safety semantics."""

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

CONTROLLERS = (
    "AuthorMeasured-PI-SMC",
    "SparseSensor-PI-SMC",
    "FiveNodeMean-PI-SMC",
    "Hotspot-MPC",
)


def _variant(
    base: FiveNodeParameters,
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


def _safe_metrics(hotspot: np.ndarray, safe_C: float | None, dt_s: float) -> dict[str, object]:
    if safe_C is None:
        return {
            "safe_limit_configured": False,
            "safe_temperature_C": None,
            "time_above_safe_limit_s": None,
            "safe_limit_violation_rmse_C": None,
            "safety_conclusion_permitted": False,
        }
    violation = np.maximum(hotspot - safe_C, 0.0)
    return {
        "safe_limit_configured": True,
        "safe_temperature_C": safe_C,
        "time_above_safe_limit_s": float(np.sum(hotspot > safe_C) * dt_s),
        "safe_limit_violation_rmse_C": float(np.sqrt(np.mean(violation**2))),
        "safety_conclusion_permitted": True,
    }


def _simulate(
    controller_name: str,
    true_model: FiveNodeParameters,
    controller_model: FiveNodeParameters,
    current_true: np.ndarray,
    current_preview: np.ndarray,
    initial: np.ndarray,
    sensors: tuple[int, ...],
    noise_sigma_C: float,
    seed: int,
    dt_s: float,
    heat_model: str,
    reference_C: float,
    safe_C: float | None,
    actuator_limit: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rng = np.random.default_rng(seed)
    measurement_noise = rng.normal(0.0, noise_sigma_C, size=(len(current_true), len(sensors)))
    temperature = np.empty((len(current_true), 5))
    estimate = np.full_like(temperature, np.nan)
    duty = np.zeros(len(current_true))
    direction = np.full(len(current_true), -1, dtype=int)
    objective = np.full(len(current_true), np.nan)
    energy_residual = np.zeros(len(current_true))
    temperature[0] = initial

    uses_ekf = controller_name != "AuthorMeasured-PI-SMC"
    observer: ExtendedKalmanFilter | None = None
    if uses_ekf:
        observer = ExtendedKalmanFilter(
            state_C=initial + np.linspace(0.25, -0.25, 5),
            covariance=np.eye(5) * 0.8,
            process_covariance=np.eye(5) * 0.0036,
            measurement_variance_C2=noise_sigma_C**2,
            sensor_indices=sensors,
            parameters=controller_model,
            heat_model=heat_model,
        )
        observer.update(temperature[0, sensors] + measurement_noise[0])
        estimate[0] = observer.state_C

    if controller_name == "AuthorMeasured-PI-SMC":
        controller: object = PISMCController(
            controller_model.controller_parameters.Kp_1_K,
            controller_model.controller_parameters.Ki_1_Ks,
            reference_C=reference_C,
            feedback_mode="middle",
            maximum_duty=actuator_limit,
        )
    elif controller_name == "SparseSensor-PI-SMC":
        controller = PISMCController(
            controller_model.controller_parameters.Kp_1_K,
            controller_model.controller_parameters.Ki_1_Ks,
            reference_C=reference_C,
            feedback_mode="middle",
            maximum_duty=actuator_limit,
        )
    elif controller_name == "FiveNodeMean-PI-SMC":
        controller = PISMCController(
            controller_model.controller_parameters.Kp_1_K,
            controller_model.controller_parameters.Ki_1_Ks,
            reference_C=reference_C,
            feedback_mode="mean",
            maximum_duty=actuator_limit,
        )
    elif controller_name == "Hotspot-MPC":
        controller = HotspotMPC(
            controller_model,
            heat_model,
            reference_C=reference_C,
            maximum_duty=actuator_limit,
            control_dt_s=dt_s,
        )
    else:
        raise ValueError(f"Unknown controller: {controller_name}")

    for k in range(len(current_true) - 1):
        if controller_name == "AuthorMeasured-PI-SMC":
            feedback = temperature[k, [0, 2, 4]]
            command, flow = controller.command(feedback, k * dt_s, dt_s)  # type: ignore[attr-defined]
        elif controller_name in ("SparseSensor-PI-SMC", "FiveNodeMean-PI-SMC"):
            assert observer is not None
            command, flow = controller.command(observer.state_C, k * dt_s, dt_s)  # type: ignore[attr-defined]
        else:
            assert observer is not None
            stride = max(1, int(round(controller.prediction_dt_s / dt_s)))  # type: ignore[attr-defined]
            preview = current_preview[k::stride][: controller.horizon_steps]  # type: ignore[attr-defined]
            command, flow, objective[k] = controller.command(observer.state_C, preview)  # type: ignore[attr-defined]
        command = min(float(command), actuator_limit)
        result = step_five_node(
            temperature[k], current_true[k], command, int(flow), dt_s, true_model, heat_model
        )
        temperature[k + 1] = result.next_temperature_C
        energy_residual[k] = result.energy_residual_W
        duty[k] = command
        direction[k] = int(flow)
        if observer is not None:
            observer.predict(current_preview[k], command, int(flow), dt_s)
            observer.update(temperature[k + 1, sensors] + measurement_noise[k + 1])
            estimate[k + 1] = observer.state_C

    duty[-1] = duty[-2]
    direction[-1] = direction[-2]
    energy_residual[-1] = energy_residual[-2]
    if controller_name == "Hotspot-MPC":
        objective[-1] = objective[-2]
    hotspot = np.max(temperature, axis=1)
    gradient = np.max(temperature, axis=1) - np.min(temperature, axis=1)
    mean_temperature = np.mean(temperature, axis=1)
    reference_error = hotspot - reference_C
    dt_hours = dt_s / 3600.0
    switches = int(np.sum(direction[1:] != direction[:-1]))
    after_60 = min(len(hotspot) - 1, int(round(60.0 / dt_s)))
    metrics: dict[str, object] = {
        "controller": controller_name,
        "feedback_source": (
            "direct measured three-zone temperatures"
            if controller_name == "AuthorMeasured-PI-SMC"
            else "nodes 2 and 5 measurements with five-node EKF state"
        ),
        "reference_temperature_C": reference_C,
        "maximum_hotspot_C": float(np.max(hotspot)),
        "maximum_hotspot_after_60s_C": float(np.max(hotspot[after_60:])),
        "reference_tracking_error_mean_C": float(np.mean(reference_error)),
        "reference_tracking_error_rmse_C": float(np.sqrt(np.mean(reference_error**2))),
        "hotspot_violation_rmse_above_reference_C": float(
            np.sqrt(np.mean(np.maximum(reference_error, 0.0) ** 2))
        ),
        "time_above_reference_s": float(np.sum(hotspot > reference_C) * dt_s),
        "maximum_gradient_C": float(np.max(gradient)),
        "mean_gradient_C": float(np.mean(gradient)),
        "gradient_rmse_C": float(np.sqrt(np.mean(gradient**2))),
        "mean_temperature_tracking_rmse_C": float(
            np.sqrt(np.mean((mean_temperature - reference_C) ** 2))
        ),
        "fan_energy_proxy_duty2_h": float(np.sum(duty**2) * dt_hours),
        "mean_absolute_duty_rate_per_s": float(np.mean(np.abs(np.diff(duty))) / dt_s),
        "airflow_direction_switches": switches,
        "direction_switch_cost_proxy": float(switches * 2.0),
        "maximum_absolute_energy_residual_W": float(np.max(np.abs(energy_residual))),
        "actuator_saturation_fraction": float(np.mean(duty >= actuator_limit - 1e-12)),
        "thermal_constraint_infeasible_fraction": float(
            np.mean((hotspot > reference_C) & (duty >= actuator_limit - 1e-12))
        ),
        **_safe_metrics(hotspot, safe_C, dt_s),
    }
    if observer is None:
        metrics["mean_estimated_hotspot_error_C"] = None
        metrics["hotspot_estimation_rmse_C"] = None
    else:
        estimated_hotspot = np.max(estimate, axis=1)
        metrics["mean_estimated_hotspot_error_C"] = float(np.mean(estimated_hotspot - hotspot))
        metrics["hotspot_estimation_rmse_C"] = float(
            np.sqrt(np.mean((estimated_hotspot - hotspot) ** 2))
        )

    data: dict[str, object] = {
        "time_s": np.arange(len(current_true)) * dt_s,
        "current_A": current_true,
        "preview_current_A": current_preview,
        "fan_duty": duty,
        "direction": direction,
        "hotspot_C": hotspot,
        "mean_temperature_C": mean_temperature,
        "gradient_C": gradient,
        "reference_temperature_C": np.full(len(current_true), reference_C),
        "safe_temperature_C": np.full(len(current_true), np.nan if safe_C is None else safe_C),
        "objective": objective,
    }
    for i in range(5):
        data[f"true_T{i + 1}_C"] = temperature[:, i]
        data[f"estimated_T{i + 1}_C"] = estimate[:, i]
    return pd.DataFrame(data), metrics


def _percentage_improvement(baseline: float, candidate: float) -> float:
    return 100.0 * (baseline - candidate) / max(abs(baseline), 1e-12)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--safe-temperature-C", type=float)
    args = parser.parse_args()
    root = args.root.resolve()
    output_root = (args.output_dir or root / "outputs_v4").resolve()
    output = output_root / "control"
    output.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load(
        (root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8")
    )
    configured_safe = config.get("safety", {}).get("safe_temperature_C")
    safe_C = args.safe_temperature_C if args.safe_temperature_C is not None else configured_safe
    if safe_C is not None:
        safe_C = float(safe_C)
    reference_C = float(config["model"]["reference_temperature_C"])
    frozen = json.loads(
        (root / "identification_v3" / "frozen_parameters_v3.json").read_text(encoding="utf-8")
    )
    thermal = load_frozen_parameters(root / "identification_v3" / "frozen_parameters_v3.json")
    heat_model = str(frozen["heat_model"])
    dt_s = float(config["project"]["dt_s"])
    seed = int(config["project"]["seed"])
    five = FiveNodeParameters.from_three_node(
        thermal, np.asarray(config["model"]["coolant_node_temperatures_forward_C"])
    )
    digitized = output_root / "digitization" / "regenerated"
    if not digitized.exists():
        digitized = root / "digitization" / "regenerated"
    fig15 = pd.read_csv(digitized / "Fig15_digitized.csv")
    fig16 = pd.read_csv(digitized / "Fig16_digitized.csv")
    base_current = np.interp(fig16.time_s, fig15.time_s, fig15.paper_simulation_current_A)
    initial_three = fig16.loc[
        0, ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]
    ].to_numpy(dtype=float)
    initial = interpolate_three_to_five(initial_three)
    hot_ambient = replace(five, coolant_forward_C=five.coolant_forward_C + 5.0)
    scenarios = {
        "nominal": dict(true=five, current=base_current, preview=base_current, noise=0.12, limit=1.0, feasible=True),
        "critical_load_plus_25pct": dict(true=five, current=np.minimum(base_current * 1.25, 42.0), preview=np.minimum(base_current * 1.25, 42.0), noise=0.12, limit=1.0, feasible=True),
        "plant_mismatch": dict(true=_variant(five, capacity=0.88, cooling=0.82, conduction=1.20), current=base_current, preview=base_current, noise=0.12, limit=1.0, feasible=True),
        "measurement_noise": dict(true=five, current=base_current, preview=base_current, noise=0.30, limit=1.0, feasible=True),
        "load_preview_minus_10pct": dict(true=five, current=base_current, preview=base_current * 0.90, noise=0.12, limit=1.0, feasible=True),
        "ambient_plus_5C": dict(true=hot_ambient, current=base_current, preview=base_current, noise=0.12, limit=1.0, feasible=True),
        "actuator_limited": dict(true=_variant(five, cooling=0.65), current=np.minimum(base_current * 1.75, 42.0), preview=np.minimum(base_current * 1.75, 42.0), noise=0.12, limit=0.15, feasible=False),
    }
    rows: list[dict[str, object]] = []
    trajectories: dict[tuple[str, str], pd.DataFrame] = {}
    sensors = (1, 4)
    for index, (scenario, values) in enumerate(scenarios.items()):
        for controller_name in CONTROLLERS:
            frame, metrics = _simulate(
                controller_name,
                values["true"],
                five,
                values["current"],
                values["preview"],
                initial,
                sensors,
                values["noise"],
                seed + index,
                dt_s,
                heat_model,
                reference_C,
                safe_C,
                values["limit"],
            )
            slug = controller_name.lower().replace("-", "_")
            frame.to_csv(output / f"{scenario}_{slug}.csv", index=False)
            trajectories[(scenario, controller_name)] = frame
            rows.append({"scenario": scenario, "feasible_scenario": values["feasible"], **metrics})
    metrics_frame = pd.DataFrame(rows)
    metrics_frame.to_csv(output / "controller_comparison_metrics.csv", index=False)

    aggregate_spec = {
        "mean_maximum_hotspot_after_60s_C": ("maximum_hotspot_after_60s_C", "mean"),
        "mean_hotspot_violation_rmse_C": ("hotspot_violation_rmse_above_reference_C", "mean"),
        "mean_gradient_rmse_C": ("gradient_rmse_C", "mean"),
        "mean_temperature_tracking_rmse_C": ("mean_temperature_tracking_rmse_C", "mean"),
        "mean_fan_energy_proxy_duty2_h": ("fan_energy_proxy_duty2_h", "mean"),
        "total_direction_switches": ("airflow_direction_switches", "sum"),
        "mean_direction_switch_cost_proxy": ("direction_switch_cost_proxy", "mean"),
        "mean_absolute_duty_rate_per_s": ("mean_absolute_duty_rate_per_s", "mean"),
    }
    feasible_metrics = metrics_frame[metrics_frame.feasible_scenario]
    infeasible_metrics = metrics_frame[~metrics_frame.feasible_scenario]
    feasible_aggregate = feasible_metrics.groupby("controller", as_index=False).agg(**aggregate_spec)
    infeasible_aggregate = infeasible_metrics.groupby("controller", as_index=False).agg(**aggregate_spec)
    feasible_aggregate.to_csv(output / "controller_aggregate_feasible.csv", index=False)
    infeasible_aggregate.to_csv(output / "controller_aggregate_infeasible.csv", index=False)

    figure, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    styles = ("-", "-.", ":", "--")
    for name, style in zip(CONTROLLERS, styles, strict=True):
        frame = trajectories[("critical_load_plus_25pct", name)]
        axes[0].plot(frame.time_s, frame.hotspot_C, style, label=name)
        axes[1].plot(frame.time_s, frame.gradient_C, style, label=name)
        axes[2].plot(frame.time_s, frame.fan_duty, style, label=name)
    axes[0].axhline(reference_C, color="k", lw=1, alpha=0.5, label="reference")
    if safe_C is not None:
        axes[0].axhline(safe_C, color="tab:red", lw=1, alpha=0.6, label="configured safe limit")
    axes[0].set_ylabel("Hotspot (°C)")
    axes[1].set_ylabel("Gradient (°C)")
    axes[2].set(xlabel="Time (s)", ylabel="Fan duty")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.suptitle("Critical-load closed-loop comparison")
    figure.tight_layout()
    figure.savefig(output / "critical_load_comparison.png", dpi=180)
    plt.close(figure)

    mpc = feasible_aggregate[feasible_aggregate.controller == "Hotspot-MPC"].iloc[0]
    author = feasible_aggregate[
        feasible_aggregate.controller == "AuthorMeasured-PI-SMC"
    ].iloc[0]
    tradeoffs = {
        "hotspot_after_60s_improvement_percent": _percentage_improvement(
            float(author.mean_maximum_hotspot_after_60s_C),
            float(mpc.mean_maximum_hotspot_after_60s_C),
        ),
        "temperature_gradient_rmse_improvement_percent": _percentage_improvement(
            float(author.mean_gradient_rmse_C), float(mpc.mean_gradient_rmse_C)
        ),
        "fan_energy_proxy_change_percent": -_percentage_improvement(
            float(author.mean_fan_energy_proxy_duty2_h), float(mpc.mean_fan_energy_proxy_duty2_h)
        ),
        "direction_switch_count_change": int(mpc.total_direction_switches - author.total_direction_switches),
        "mean_temperature_tracking_rmse_change_percent": -_percentage_improvement(
            float(author.mean_temperature_tracking_rmse_C),
            float(mpc.mean_temperature_tracking_rmse_C),
        ),
    }
    limited_mpc = infeasible_metrics[
        (infeasible_metrics.scenario == "actuator_limited")
        & (infeasible_metrics.controller == "Hotspot-MPC")
    ].iloc[0]
    summary = {
        "reference_temperature_C": reference_C,
        "reference_semantics": "paper optimum/control reference only; not a safety upper limit",
        "safe_temperature_C": safe_C,
        "safety_conclusion": (
            "A configured external safe limit was evaluated."
            if safe_C is not None
            else "No externally supported safe limit was supplied; no safety conclusion is made."
        ),
        "controller_definitions": {
            "AuthorMeasured-PI-SMC": "direct true/measured temperatures at the three regional nodes; ideal paper baseline",
            "SparseSensor-PI-SMC": "nodes 2 and 5 plus EKF; middle-node PI-SMC feedback",
            "FiveNodeMean-PI-SMC": "nodes 2 and 5 plus EKF; five-node mean PI-SMC feedback",
            "Hotspot-MPC": "nodes 2 and 5 plus EKF; three-block hotspot MPC",
        },
        "control_blocking": {"horizon_steps": 12, "block_steps": [[1, 4], [5, 8], [9, 12]]},
        "direction_reversal": {
            "fixed_objective_penalty": 2.0,
            "minimum_dwell_s": 10.0,
            "reversal_dead_time_s": 2.0,
        },
        "feasible_scenario_count": int(feasible_metrics.scenario.nunique()),
        "infeasible_scenarios_excluded_from_normal_aggregate": ["actuator_limited"],
        "actuator_limited_declared_infeasible": bool(
            float(limited_mpc.thermal_constraint_infeasible_fraction) > 0.0
        ),
        "mpc_hotspot_improvement_still_holds": bool(
            tradeoffs["hotspot_after_60s_improvement_percent"] > 0.0
        ),
        "mpc_tradeoffs_relative_to_author_measured_pi_smc_feasible_scenarios": tradeoffs,
    }
    (output / "control_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
