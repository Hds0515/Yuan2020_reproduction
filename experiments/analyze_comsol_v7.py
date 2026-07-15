"""Create the frozen V7 COMSOL/Fig.7 validation table and gate decision."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
VELOCITIES = np.array([4, 5, 6, 8, 10, 12], dtype=float)
VALIDATION = {5.0, 6.0, 10.0, 12.0}


def load_comsol(name: str) -> pd.DataFrame:
    frame = pd.read_csv(
        ROOT / "comsol" / "v7_runtime" / f"physical_heat_metrics_{name}.csv",
        comment="%",
        header=None,
    )
    frame.columns = [
        "velocity_m_s", "T_inlet_region_C", "T_middle_region_C",
        "T_outlet_region_C", "Tmax_C", "Tmin_C", "DeltaT_C",
        "hotspot_location_normalized", "hotspot_y_m", "hotspot_z_m",
        "air_inlet_C", "air_outlet_C", "air_enthalpy_gain_W",
        "natural_convection_loss_W", "radiation_loss_W", "conduction_loss_W",
        "total_heat_input_W", "energy_residual_relative", "mdot_cell_kg_s",
        "mdot_stack_kg_s", "mdot_identity_relative_error",
    ]
    return frame


def main() -> None:
    target = pd.read_csv(
        ROOT / "outputs_v5" / "digitization" / "external" / "Fig7_region_summary.csv"
    )
    target = target.sort_values("velocity_m_s").reset_index(drop=True)
    for column in ["T_inlet_region_K", "T_middle_region_K", "T_outlet_region_K", "Tmax_K", "Tmin_K", "DeltaT_K"]:
        target[column.replace("_K", "_C")] = target[column] - (273.15 if column != "DeltaT_K" else 0.0)
    coarse, medium, fine = (load_comsol(name) for name in ["coarse", "medium", "fine"])
    rows = []
    for index, velocity in enumerate(VELOCITIES):
        predicted = fine.iloc[index]
        observed = target.iloc[index]
        region_errors = np.array([
            predicted.T_inlet_region_C - observed.T_inlet_region_C,
            predicted.T_middle_region_C - observed.T_middle_region_C,
            predicted.T_outlet_region_C - observed.T_outlet_region_C,
        ])
        medium_row = medium.iloc[index]
        mesh_region_difference = max(
            abs(medium_row.T_inlet_region_C - predicted.T_inlet_region_C),
            abs(medium_row.T_middle_region_C - predicted.T_middle_region_C),
            abs(medium_row.T_outlet_region_C - predicted.T_outlet_region_C),
        )
        mesh_tmax_relative = abs(medium_row.Tmax_C - predicted.Tmax_C) / abs(predicted.Tmax_C)
        role = "calibration" if velocity in {4.0, 8.0} else "independent_validation"
        rows.append({
            "velocity_m_s": velocity,
            "data_role": role,
            "T_inlet_region_C": predicted.T_inlet_region_C,
            "T_middle_region_C": predicted.T_middle_region_C,
            "T_outlet_region_C": predicted.T_outlet_region_C,
            "Tmax_C": predicted.Tmax_C,
            "Tmin_C": predicted.Tmin_C,
            "DeltaT_C": predicted.DeltaT_C,
            "hotspot_location_normalized": predicted.hotspot_location_normalized,
            "air_outlet_temperature_C": predicted.air_outlet_C,
            "air_enthalpy_gain_W": predicted.air_enthalpy_gain_W,
            "natural_convection_loss_W": predicted.natural_convection_loss_W,
            "radiation_loss_W": predicted.radiation_loss_W,
            "conduction_loss_W": predicted.conduction_loss_W,
            "total_heat_input_W": predicted.total_heat_input_W,
            "energy_residual_relative": predicted.energy_residual_relative,
            "mdot_cell_kg_s": predicted.mdot_cell_kg_s,
            "mdot_stack_kg_s": predicted.mdot_stack_kg_s,
            "mdot_identity_relative_error": predicted.mdot_identity_relative_error,
            "Fig7_T_inlet_region_C": observed.T_inlet_region_C,
            "Fig7_T_middle_region_C": observed.T_middle_region_C,
            "Fig7_T_outlet_region_C": observed.T_outlet_region_C,
            "Fig7_Tmax_C": observed.Tmax_C,
            "Fig7_Tmin_C": observed.Tmin_C,
            "Fig7_DeltaT_C": observed.DeltaT_C,
            "Fig7_hotspot_location_normalized": observed.hotspot_location_normalized,
            "region_RMSE_C": float(np.sqrt(np.mean(region_errors**2))),
            "hotspot_location_absolute_error": abs(
                predicted.hotspot_location_normalized - observed.hotspot_location_normalized
            ),
            "medium_fine_Tmax_relative_difference": mesh_tmax_relative,
            "medium_fine_max_region_difference_C": mesh_region_difference,
            "numerically_converged": bool(
                mesh_tmax_relative < 0.01
                and mesh_region_difference < 0.2
                and predicted.energy_residual_relative < 0.005
                and predicted.mdot_identity_relative_error < 1e-10
            ),
            "independent_RMSE_gate_passed": bool(
                role == "calibration" or np.sqrt(np.mean(region_errors**2)) <= 1.5
            ),
        })
    output = pd.DataFrame(rows)
    destination = ROOT / "outputs_v7"
    destination.mkdir(exist_ok=True)
    output.to_csv(destination / "comsol_fig7_validation.csv", index=False)
    correlation = float(spearmanr(
        output.hotspot_location_normalized,
        output.Fig7_hotspot_location_normalized,
    ).statistic)
    validation = output[output.data_role == "independent_validation"]
    summary = {
        "calibration_parameters": {
            "h_external_W_m2K": 30.0,
            "emissivity": 1.0,
            "both_at_upper_bound": True,
        },
        "calibration_region_RMSE_C": float(
            output.loc[output.data_role == "calibration", "region_RMSE_C"].mean()
        ),
        "independent_validation_region_RMSE_C": float(validation.region_RMSE_C.mean()),
        "independent_validation_worst_region_RMSE_C": float(validation.region_RMSE_C.max()),
        "all_independent_speeds_below_1p5C": bool(validation.independent_RMSE_gate_passed.all()),
        "hotspot_location_spearman_correlation": correlation,
        "hotspot_location_trend_consistent": bool(correlation >= 0.5),
        "temperature_scale_order_consistent_with_Fig7": bool(
            output.Tmin_C.min() >= 20.0 and output.Tmax_C.max() <= 70.0
        ),
        "all_numerically_converged": bool(output.numerically_converged.all()),
        "maximum_medium_fine_Tmax_relative_difference_percent": float(
            100 * output.medium_fine_Tmax_relative_difference.max()
        ),
        "maximum_medium_fine_region_difference_C": float(
            output.medium_fine_max_region_difference_C.max()
        ),
        "maximum_energy_residual_percent": float(100 * output.energy_residual_relative.max()),
        "physically_validated": bool(
            output.numerically_converged.all()
            and validation.independent_RMSE_gate_passed.all()
            and correlation >= 0.5
            and output.Tmin_C.min() >= 20.0
            and output.Tmax_C.max() <= 70.0
        ),
        "downstream_high_fidelity_dataset_authorized": False,
    }
    summary["downstream_high_fidelity_dataset_authorized"] = summary["physically_validated"]
    (destination / "comsol_fig7_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
