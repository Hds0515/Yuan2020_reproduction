"""Exhaustive V6 critical-case resource search and frozen-controller audit."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import itertools
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from experiments.run_control_fair_comparison import simulate_case, variant  # noqa: E402
from models.five_node_model import FiveNodeParameters, interpolate_three_to_five  # noqa: E402
from models.model_registry import load_primary_model, write_model_runtime_record  # noqa: E402
from models.three_node_model import load_frozen_parameters  # noqa: E402


def _run_search_item(payload):
    config, true_model, control_model, current, preview, initial, sensors, sigma, seed, dt, heat = payload
    _, metrics = simulate_case(
        "Hotspot-MPC", true_model, control_model, current, preview, initial,
        sensors, sigma, seed, dt, heat, 55.0, 1.0, config,
    )
    return {**config, **metrics}


def _improvement(author: pd.Series, candidate: pd.Series) -> dict[str, float | int]:
    return {
        "maximum_hotspot_absolute_decrease_C": float(
            author.maximum_hotspot_C - candidate.maximum_hotspot_C
        ),
        "above_55C_RMSE_decrease_percent": float(
            100.0 * (author.hotspot_excess_rmse_C - candidate.hotspot_excess_rmse_C)
            / max(author.hotspot_excess_rmse_C, 1e-12)
        ),
        "gradient_RMSE_decrease_percent": float(
            100.0 * (author.gradient_rmse_C - candidate.gradient_rmse_C)
            / max(author.gradient_rmse_C, 1e-12)
        ),
        "fan_energy_change_percent": float(
            100.0 * (candidate.fan_energy_proxy_duty2_h / author.fan_energy_proxy_duty2_h - 1.0)
        ),
        "direction_switch_count_change": int(
            candidate.direction_switch_count - author.direction_switch_count
        ),
        "mean_temperature_tracking_RMSE_change_percent": float(
            100.0 * (candidate.mean_temperature_tracking_rmse_C
                     / max(author.mean_temperature_tracking_rmse_C, 1e-12) - 1.0)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--workers", type=int, default=min(6, os.cpu_count() or 1))
    args = parser.parse_args()
    root = args.root.resolve()
    output_root = (args.output_dir or root / "outputs_v6").resolve()
    output = output_root / "control_equal_resource"
    output.mkdir(parents=True, exist_ok=True)

    primary, primary_path, model_hash = load_primary_model(root)
    write_model_runtime_record(root, output, "control_equal_resource_v6")
    config = yaml.safe_load(
        (root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8")
    )
    thermal = load_frozen_parameters(primary_path)
    heat_model = str(primary["heat_model"])
    five = FiveNodeParameters.from_three_node(
        thermal, np.asarray(config["model"]["coolant_node_temperatures_forward_C"], float)
    )
    seed = int(config["project"]["seed"])
    dt = 10.0
    digitized = root / "outputs_v4" / "digitization" / "regenerated"
    fig15 = pd.read_csv(digitized / "Fig15_digitized.csv")
    fig16 = pd.read_csv(digitized / "Fig16_digitized.csv")
    base_current = np.interp(
        fig16.time_s, fig15.time_s, fig15.paper_simulation_current_A
    )[:: int(dt)]
    initial = interpolate_three_to_five(
        fig16.loc[
            0, ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]
        ].to_numpy(float)
    )
    sensors = (0, 4)
    critical_current = np.minimum(base_current * 1.25, 42.0)
    _, baseline_metrics = simulate_case(
        "AuthorMeasured-PI-SMC", five, five, critical_current, critical_current,
        initial, sensors, 0.12, seed, dt, heat_model, 55.0, 1.0,
    )

    combinations = list(itertools.product(
        [1.0, 2.0, 5.0, 10.0, 20.0, 50.0],
        [0.5, 1.0, 2.0, 4.0],
        [0.2, 0.5, 1.0, 2.0],
        [10.0, 30.0, 60.0],
        [1.0, 4.0, 10.0, 30.0],
    ))
    configs = [
        {
            "config_id": index,
            "fan_weight": fan,
            "hotspot_weight": hotspot,
            "gradient_weight": gradient,
            "dwell": dwell,
            "switch_weight": switch,
            "maximum_duty_step": 0.05,
            "duty_move_weight": 0.8,
        }
        for index, (fan, hotspot, gradient, dwell, switch) in enumerate(combinations)
    ]
    payloads = [
        (
            item, five, five, critical_current, critical_current, initial,
            sensors, 0.12, seed, dt, heat_model,
        )
        for item in configs
    ]
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        search_rows = list(executor.map(_run_search_item, payloads, chunksize=8))
    search = pd.DataFrame(search_rows).sort_values("config_id")
    search["energy_relative_difference_percent"] = 100.0 * (
        search.fan_energy_proxy_duty2_h / baseline_metrics["fan_energy_proxy_duty2_h"] - 1.0
    )
    search["within_5pct_equal_energy"] = (
        search.energy_relative_difference_percent.abs() <= 5.0
    )
    search.to_csv(output / "critical_full_factorial_1152.csv", index=False)

    equal = search[search.within_5pct_equal_energy]
    equal_found = not equal.empty
    if equal_found:
        frozen = equal.sort_values(
            ["hotspot_excess_rmse_C", "gradient_rmse_C", "direction_switch_count"]
        ).iloc[0]
        selection_basis = "best hotspot RMSE among critical points within ±5% fan energy"
    else:
        frozen = search.iloc[search.energy_relative_difference_percent.abs().argmin()]
        selection_basis = "nearest diagnostic point only; no equal-energy claim"
    same_dwell_pool = search[search.dwell == 10.0]
    equal_switch = same_dwell_pool.assign(
        switch_distance=(same_dwell_pool.direction_switch_count
                         - baseline_metrics["direction_switch_count"]).abs()
    ).sort_values(
        ["switch_distance", "hotspot_excess_rmse_C", "gradient_rmse_C"]
    ).iloc[0]
    selection = {
        "baseline_critical": baseline_metrics,
        "equal_energy_match_count": int(len(equal)),
        "equal_energy_match_found": bool(equal_found),
        "selection_basis": selection_basis,
        "frozen_resource_point": frozen.to_dict(),
        "same_10s_dwell_nearest_switch_point": equal_switch.to_dict(),
    }
    (output / "frozen_resource_point.json").write_text(
        json.dumps(selection, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    scenarios = {
        "nominal": (five, base_current, base_current, 0.12, 1.0, True),
        "critical": (five, critical_current, critical_current, 0.12, 1.0, True),
        "mismatch": (variant(five, .88, .82, 1.2), base_current, base_current, .12, 1.0, True),
        "noise": (five, base_current, base_current, .30, 1.0, True),
        "preview_minus10": (five, base_current, base_current * .9, .12, 1.0, True),
        "ambient_plus5": (variant(five, offset=5), base_current, base_current, .12, 1.0, True),
        "actuator_limited": (
            variant(five, cooling=.65), np.minimum(base_current * 1.75, 42),
            np.minimum(base_current * 1.75, 42), .12, .15, False,
        ),
    }
    frozen_config = {
        key: frozen[key] for key in [
            "fan_weight", "hotspot_weight", "gradient_weight", "dwell",
            "switch_weight", "maximum_duty_step", "duty_move_weight",
        ]
    }
    mean_config = {**frozen_config, "hotspot_weight": 0.0, "mean_weight": 8.0}
    result_rows = []
    for scenario_index, (scenario, values) in enumerate(scenarios.items()):
        true_model, current, preview, sigma, limit, feasible = values
        for controller, controller_config in [
            ("AuthorMeasured-PI-SMC", None),
            ("SparseSensor-PI-SMC", None),
            ("MeanTemperature-MPC", mean_config),
            ("Hotspot-MPC", frozen_config),
        ]:
            _, metrics = simulate_case(
                controller, true_model, five, current, preview, initial, sensors,
                sigma, seed + scenario_index, dt, heat_model, 55.0, limit,
                controller_config,
            )
            result_rows.append({"scenario": scenario, "feasible": feasible, **metrics})
    results = pd.DataFrame(result_rows)
    results.to_csv(output / "all_scenario_controller_metrics.csv", index=False)
    feasible = results[results.feasible]
    aggregate = feasible.groupby("controller", as_index=False).agg(
        maximum_hotspot_C=("maximum_hotspot_C", "max"),
        mean_above_55C_RMSE_C=("hotspot_excess_rmse_C", "mean"),
        mean_gradient_RMSE_C=("gradient_rmse_C", "mean"),
        mean_fan_energy_proxy_duty2_h=("fan_energy_proxy_duty2_h", "mean"),
        total_direction_switches=("direction_switch_count", "sum"),
        mean_temperature_tracking_RMSE_C=("mean_temperature_tracking_rmse_C", "mean"),
    )
    aggregate.to_csv(output / "feasible_scenario_aggregate.csv", index=False)
    author = aggregate[aggregate.controller == "AuthorMeasured-PI-SMC"].iloc[0]
    hotspot = aggregate[aggregate.controller == "Hotspot-MPC"].iloc[0]
    # Rename aggregate columns for the common comparison helper.
    author_common = pd.Series({
        "maximum_hotspot_C": author.maximum_hotspot_C,
        "hotspot_excess_rmse_C": author.mean_above_55C_RMSE_C,
        "gradient_rmse_C": author.mean_gradient_RMSE_C,
        "fan_energy_proxy_duty2_h": author.mean_fan_energy_proxy_duty2_h,
        "direction_switch_count": author.total_direction_switches,
        "mean_temperature_tracking_rmse_C": author.mean_temperature_tracking_RMSE_C,
    })
    hotspot_common = pd.Series({
        "maximum_hotspot_C": hotspot.maximum_hotspot_C,
        "hotspot_excess_rmse_C": hotspot.mean_above_55C_RMSE_C,
        "gradient_rmse_C": hotspot.mean_gradient_RMSE_C,
        "fan_energy_proxy_duty2_h": hotspot.mean_fan_energy_proxy_duty2_h,
        "direction_switch_count": hotspot.total_direction_switches,
        "mean_temperature_tracking_rmse_C": hotspot.mean_temperature_tracking_RMSE_C,
    })
    critical_author = pd.Series(baseline_metrics)
    critical_frozen = frozen
    critical_switch = equal_switch
    summary = {
        "primary_model_sha256": model_hash,
        "critical_search_configuration_count": int(len(search)),
        "critical_equal_energy_match_count": int(len(equal)),
        "equal_energy_match_found_within_5pct": bool(equal_found),
        "frozen_point_energy_relative_difference_percent": float(
            frozen.energy_relative_difference_percent
        ),
        "critical_equal_energy_comparison": _improvement(
            critical_author, critical_frozen
        ) if equal_found else None,
        "all_feasible_scenarios_frozen_point_comparison": _improvement(
            author_common, hotspot_common
        ) if equal_found else None,
        "same_10s_dwell_nearest_switch_comparison": {
            "MPC_switch_count": int(equal_switch.direction_switch_count),
            "author_switch_count": int(baseline_metrics["direction_switch_count"]),
            **_improvement(critical_author, critical_switch),
        },
        "reference_temperature_C": 55.0,
        "reference_is_safety_limit": False,
        "actuator_limited_excluded_from_feasible_aggregate": True,
        "Hotspot_MPC_equal_energy_advantage_holds": bool(
            equal_found
            and frozen.hotspot_excess_rmse_C < baseline_metrics["hotspot_excess_rmse_C"]
            and frozen.maximum_hotspot_C < baseline_metrics["maximum_hotspot_C"]
        ),
    }
    (output / "control_equal_resource_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
