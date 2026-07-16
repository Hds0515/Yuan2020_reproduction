"""Build the v4 machine-readable status and ten-question stage report."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs_v4")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    identification = read_json(output / "identification" / "identification_metrics.json")
    bootstrap = read_json(output / "identification" / "bootstrap_summary.json")
    five = read_json(output / "five_node" / "consistency_metrics.json")
    observer = read_json(output / "observer" / "observer_summary.json")
    control = read_json(output / "control" / "control_summary.json")
    comsol = read_json(ROOT / "comsol" / "status_v4.json")
    detail = read_json(output / "identification" / "local_detail_validation_impact.json")

    answers = {
        "1_multiple_model_better_in_majority_of_true_mismatches": {
            "answer": bool(observer["model_bank_improves_majority_of_true_mismatches"]),
            "improved_scenarios": int(observer["model_bank_improved_scenario_count"]),
            "total_true_mismatch_scenarios": int(observer["true_mismatch_scenario_count"]),
            "mean_hotspot_improvement_percent": float(observer["mean_hotspot_improvement_percent"]),
            "paired_bootstrap_95pct_CI_percent": observer[
                "paired_bootstrap_mean_improvement_95pct_CI_percent"
            ],
        },
        "2_mpc_hotspot_improvement_still_holds": {
            "answer": bool(control["mpc_hotspot_improvement_still_holds"]),
            "hotspot_after_60s_improvement_percent": control[
                "mpc_tradeoffs_relative_to_author_measured_pi_smc_feasible_scenarios"
            ]["hotspot_after_60s_improvement_percent"],
        },
        "3_mpc_extra_fan_and_direction_cost": control[
            "mpc_tradeoffs_relative_to_author_measured_pi_smc_feasible_scenarios"
        ],
        "4_55C_is_reference_only": {
            "answer": True,
            "safe_temperature_C": control["safe_temperature_C"],
            "safety_conclusion": control["safety_conclusion"],
        },
        "5_comsol_inlet_area_mass_flow_consistent_with_3962mm2": {
            "geometry_and_formula_consistent": bool(comsol["area_scaling_formula_consistent"]),
            "numerically_verified_by_converged_solution": bool(
                comsol["stage1_mass_balance_passed"]
            ),
            "decision": comsol["equivalent_geometry_decision"],
        },
        "6_first_converged_0p1ms_flow_solution": bool(comsol["stage1_flow_0p1ms_converged"]),
        "7_converged_4ms_flow_solution": bool(comsol["stage2_flow_4ms_converged"]),
        "8_first_sequential_flow_thermal_solution": bool(
            comsol["stage3_one_way_thermal_converged"]
        ),
        "9_energy_conservation_passed": {
            "python_lumped_models": bool(
                float(identification["maximum_relative_energy_residual"]) < 1e-6
                and float(five["maximum_absolute_energy_residual_W"]) < 1e-6
            ),
            "comsol_stage3": bool(comsol["stage3_energy_balance_passed"]),
        },
        "10_remaining_tasks": comsol["remaining_tasks"],
    }
    summary = {
        "curve_digitization_reproducible": all(
            (output / "digitization" / "regenerated" / f"Fig{name}_digitized.csv").exists()
            for name in (14, 15, 16)
        ),
        "local_detail_digitization_completed": all(
            (output / "digitization" / "regenerated" / f"Fig16{panel}_digitized.csv").exists()
            for panel in ("b", "c")
        ),
        "independent_validation_rmse_C": float(
            identification["fig16_independent_validation_rmse_C"]
        ),
        "local_detail_validation_rmse_C": float(
            detail["Fig16c_local_panel_independent_validation_rmse_C"]
        ),
        "bootstrap_replicates": int(bootstrap["replicate_count"]),
        "five_node_energy_consistency_passed": bool(
            float(five["maximum_absolute_energy_residual_W"]) < 1e-6
            and abs(float(five["area_preservation_ratio"]) - 1.0) < 1e-12
        ),
        "temperature_semantics_separated": True,
        "observer_method": observer["method_classification"],
        "actuator_limited_excluded_from_feasible_aggregate": True,
        "comsol_stage1_converged": bool(comsol["stage1_flow_0p1ms_converged"]),
        "no_synthetic_comsol_results": True,
        "stage_answers": answers,
    }
    (output / "final_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    with (output / "final_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["item", "value"])
        for key, value in summary.items():
            writer.writerow(
                [key, json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value]
            )

    tradeoffs = answers["3_mpc_extra_fan_and_direction_cost"]
    report = f"""# Yuan2020 v4 stage acceptance report

1. Multiple-model observer better in most true mismatches: **{answers['1_multiple_model_better_in_majority_of_true_mismatches']['answer']}** ({answers['1_multiple_model_better_in_majority_of_true_mismatches']['improved_scenarios']}/{answers['1_multiple_model_better_in_majority_of_true_mismatches']['total_true_mismatch_scenarios']}; mean hotspot RMSE improvement {answers['1_multiple_model_better_in_majority_of_true_mismatches']['mean_hotspot_improvement_percent']:.2f}%).
2. MPC hotspot improvement still holds: **{answers['2_mpc_hotspot_improvement_still_holds']['answer']}** ({answers['2_mpc_hotspot_improvement_still_holds']['hotspot_after_60s_improvement_percent']:.2f}% after 60 s over feasible scenarios).
3. MPC trade-offs versus AuthorMeasured-PI-SMC: fan-energy proxy {tradeoffs['fan_energy_proxy_change_percent']:+.2f}%, direction switches {tradeoffs['direction_switch_count_change']:+d}, mean-temperature tracking RMSE {tradeoffs['mean_temperature_tracking_rmse_change_percent']:+.2f}%, while gradient RMSE improves {tradeoffs['temperature_gradient_rmse_improvement_percent']:.2f}%.
4. 55 °C is only the paper optimum/control reference: **True**. No externally supported safe temperature was supplied, so no safety conclusion is made.
5. COMSOL Scheme B uses `areaScale=Ain/(W*Hair)` and the target `rho*u*Ain`, so the formulation is consistent with 3962 mm²; numerical mass-balance verification is unavailable because Stage 1 did not converge.
6. First converged 0.1 m/s flow solution: **{answers['6_first_converged_0p1ms_flow_solution']}**.
7. Converged 4 m/s flow solution: **{answers['7_converged_4ms_flow_solution']}** (not attempted after Stage 1 failure).
8. First sequential flow-to-thermal solution: **{answers['8_first_sequential_flow_thermal_solution']}**.
9. Energy conservation: Python lumped models **{answers['9_energy_conservation_passed']['python_lumped_models']}**; COMSOL Stage 3 **{answers['9_energy_conservation_passed']['comsol_stage3']}**.
10. Remaining tasks: {'; '.join(answers['10_remaining_tasks'])}.

No COMSOL temperature, mesh-convergence, mass-balance, or energy-balance value was fabricated.
"""
    (output / "stage_acceptance_report.md").write_text(report, encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
