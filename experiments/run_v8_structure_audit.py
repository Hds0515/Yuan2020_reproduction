"""Run the frozen V8 Fig.7 digitization and cooling-structure audit."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from comsol.v8_runtime.cooling_structure_models import (  # noqa: E402
    M1_GEOMETRY,
    PHYSICS,
    churchill_chu_h_W_m2K,
    solve_structure,
)
from digitization.digitize_external_figures import (  # noqa: E402
    FIG7_PANELS,
    colorbar_palette,
)


VELOCITIES = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
CALIBRATION = {4.0, 8.0}
VALIDATION = {5.0, 6.0, 10.0, 12.0}
REGIONS = {
    "inlet": [0.10, 0.25],
    "middle": [0.40, 0.60],
    "outlet": [0.75, 0.90],
}
X_FRACTIONS = [0.35, 0.65]
REGION_COLUMNS = [
    "T_inlet_region_C",
    "T_middle_region_C",
    "T_outlet_region_C",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode_patch(
    image: np.ndarray,
    tree: cKDTree,
    palette_temperature: np.ndarray,
    x: int,
    y: int,
) -> tuple[float, int]:
    patch = image[y - 4 : y + 5, x - 4 : x + 5].reshape(-1, 3)
    valid = patch[
        (patch.mean(axis=1) > 35)
        & (patch.mean(axis=1) < 245)
        & ((patch.max(axis=1) - patch.min(axis=1)) > 25)
    ]
    distance, index = tree.query(valid)
    accepted = palette_temperature[index][distance < 55]
    if len(accepted) < 10:
        raise RuntimeError(f"Insufficient coloured pixels at ({x}, {y})")
    return float(np.median(accepted)), int(len(accepted))


def digitization_reaudit(output: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    image_path = ROOT / "source_images" / "Fig7.jpg"
    image = np.asarray(Image.open(image_path).convert("RGB"))
    automatic = pd.read_csv(
        ROOT / "outputs_v5/digitization/external/Fig7_region_summary.csv"
    ).set_index("velocity_m_s")
    point_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    line_profiles: dict[float, tuple[np.ndarray, np.ndarray]] = {}

    for velocity, (plate, bar, t_min, t_max) in FIG7_PANELS.items():
        palette, temperature = colorbar_palette(image, bar, t_min, t_max)
        # Human colour-bar review resolves the duplicated red endpoint: the fin
        # core follows the high-temperature branch (>307 K); the low red branch
        # belongs to the printed bar/manifold and is not sampled.
        high_branch = temperature >= 307.0
        tree = cKDTree(palette[high_branch])
        branch_temperature = temperature[high_branch]
        x0, y0, x1, y1 = plate
        manual_region_K: dict[str, float] = {}
        for region_name, flow_positions in REGIONS.items():
            values = []
            for flow_coordinate in flow_positions:
                for x_fraction in X_FRACTIONS:
                    x = round(x0 + x_fraction * (x1 - x0))
                    y = round(y1 - flow_coordinate * (y1 - y0))
                    value, count = decode_patch(
                        image, tree, branch_temperature, x, y
                    )
                    values.append(value)
                    point_rows.append(
                        {
                            "velocity_m_s": float(velocity),
                            "region": region_name,
                            "flow_coordinate_inlet_0_outlet_1": flow_coordinate,
                            "cross_stream_fraction": x_fraction,
                            "pixel_x": x,
                            "pixel_y": y,
                            "temperature_K": value,
                            "accepted_patch_pixel_count": count,
                        }
                    )
            manual_region_K[region_name] = float(np.mean(values))

        flow_grid = np.linspace(0.03, 0.97, 39)
        profile = []
        for flow_coordinate in flow_grid:
            cross_values = []
            for x_fraction in X_FRACTIONS:
                x = round(x0 + x_fraction * (x1 - x0))
                y = round(y1 - flow_coordinate * (y1 - y0))
                value, _ = decode_patch(image, tree, branch_temperature, x, y)
                cross_values.append(value)
            profile.append(float(np.mean(cross_values)))
        profile_array = np.asarray(profile)
        line_profiles[float(velocity)] = (flow_grid, profile_array)
        manual_tmax = float(np.quantile(profile_array, 0.98))
        manual_tmin = float(np.quantile(profile_array, 0.02))
        manual_hotspot = float(flow_grid[int(np.argmax(profile_array))])
        auto = automatic.loc[float(velocity)]
        errors = np.array(
            [
                manual_region_K["inlet"] - auto.T_inlet_region_K,
                manual_region_K["middle"] - auto.T_middle_region_K,
                manual_region_K["outlet"] - auto.T_outlet_region_K,
            ]
        )
        summary_rows.append(
            {
                "velocity_m_s": float(velocity),
                "manual_T_inlet_region_K": manual_region_K["inlet"],
                "manual_T_middle_region_K": manual_region_K["middle"],
                "manual_T_outlet_region_K": manual_region_K["outlet"],
                "automatic_T_inlet_region_K": auto.T_inlet_region_K,
                "automatic_T_middle_region_K": auto.T_middle_region_K,
                "automatic_T_outlet_region_K": auto.T_outlet_region_K,
                "manual_minus_auto_inlet_K": errors[0],
                "manual_minus_auto_middle_K": errors[1],
                "manual_minus_auto_outlet_K": errors[2],
                "region_RMSE_K": float(np.sqrt(np.mean(errors**2))),
                "manual_Tmax_K": manual_tmax,
                "manual_Tmin_K": manual_tmin,
                "manual_DeltaT_K": manual_tmax - manual_tmin,
                "manual_hotspot_location_normalized": manual_hotspot,
                "automatic_Tmax_K": auto.Tmax_K,
                "automatic_Tmin_K": auto.Tmin_K,
                "automatic_DeltaT_K": auto.DeltaT_K,
                "automatic_hotspot_location_normalized": auto.hotspot_location_normalized,
            }
        )

    points = pd.DataFrame(point_rows)
    comparison = pd.DataFrame(summary_rows)
    points.to_csv(output / "fig7_manual_points.csv", index=False)
    comparison.to_csv(output / "fig7_digitization_manual_comparison.csv", index=False)

    figure, axis = plt.subplots(figsize=(16, 7))
    axis.imshow(image)
    colours = {"inlet": "#2b83ba", "middle": "#fdae61", "outlet": "#d7191c"}
    for velocity, (plate, _, _, _) in FIG7_PANELS.items():
        x0, y0, x1, y1 = plate
        height = y1 - y0
        for name, lower, upper in [
            ("inlet", 0.03, 1.0 / 3.0),
            ("middle", 1.0 / 3.0, 2.0 / 3.0),
            ("outlet", 2.0 / 3.0, 0.97),
        ]:
            y_top = y1 - upper * height
            y_bottom = y1 - lower * height
            rectangle = plt.Rectangle(
                (x0, y_top),
                x1 - x0,
                y_bottom - y_top,
                facecolor=colours[name],
                edgecolor=colours[name],
                alpha=0.10,
                linewidth=1.0,
            )
            axis.add_patch(rectangle)
        subset = points[points.velocity_m_s == float(velocity)]
        for name, group in subset.groupby("region"):
            axis.scatter(
                group.pixel_x,
                group.pixel_y,
                s=22,
                c=colours[name],
                edgecolors="white",
                linewidths=0.5,
            )
        axis.text(x0 + 8, y0 + 20, f"{velocity} m/s", color="white", fontsize=8,
                  bbox={"facecolor": "black", "alpha": 0.55, "pad": 2})
    axis.set_title("Fig.7 V8 manual read points and frozen region definitions")
    axis.set_axis_off()
    figure.tight_layout()
    figure.savefig(output / "fig7_region_definition.png", dpi=220)
    plt.close(figure)

    overall_rmse = float(
        np.sqrt(
            np.mean(
                comparison[
                    [
                        "manual_minus_auto_inlet_K",
                        "manual_minus_auto_middle_K",
                        "manual_minus_auto_outlet_K",
                    ]
                ].to_numpy() ** 2
            )
        )
    )
    audit = {
        "source_image": "source_images/Fig7.jpg",
        "source_image_sha256": sha256(image_path),
        "manual_point_count": int(len(points)),
        "manual_points_per_region_per_speed": 4,
        "colour_bar_review": {
            "left_panels_K": [295.0, 338.0],
            "right_panels_K": [295.0, 335.0],
            "cyclic_red_endpoint_ambiguity": True,
            "manual_resolution": "sample only the coloured fin core and map its continuous gradient to the >=307 K colour-bar branch",
        },
        "overall_manual_vs_automatic_region_RMSE_K": overall_rmse,
        "maximum_per_speed_region_RMSE_K": float(comparison.region_RMSE_K.max()),
        "redigitization_trigger_threshold_K": 1.0,
        "redigitization_triggered": bool(overall_rmse > 1.0),
        "digitization_reliable": bool(overall_rmse <= 1.0),
        "validation_threshold_changed": False,
        "fixed_validation_target": "outputs_v5/digitization/external/Fig7_region_summary.csv",
    }
    (output / "digitization_reaudit_summary.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return points, comparison, audit


def target_table() -> pd.DataFrame:
    target = pd.read_csv(
        ROOT / "outputs_v5/digitization/external/Fig7_region_summary.csv"
    )
    for column in [
        "T_inlet_region_K",
        "T_middle_region_K",
        "T_outlet_region_K",
        "Tmax_K",
        "Tmin_K",
    ]:
        target[column.replace("_K", "_C")] = target[column] - 273.15
    target["DeltaT_C"] = target.DeltaT_K
    return target.set_index("velocity_m_s")


def energy_partition(output: Path) -> pd.DataFrame:
    v7 = pd.read_csv(
        ROOT / "comsol/v7_runtime/physical_heat_metrics_fine.csv",
        comment="%",
        header=None,
    )
    v7.columns = [
        "velocity_m_s", "T_inlet_region_C", "T_middle_region_C",
        "T_outlet_region_C", "Tmax_C", "Tmin_C", "DeltaT_C",
        "hotspot_location_normalized", "hotspot_y_m", "hotspot_z_m",
        "air_inlet_C", "air_outlet_C", "Q_air_W", "Q_natural_W",
        "Q_radiation_W", "Q_conduction_W", "Q_total_W",
        "energy_residual_relative", "mdot_cell_kg_s", "mdot_stack_kg_s",
        "mdot_identity_relative_error",
    ]
    for name in ["air", "natural", "radiation", "conduction"]:
        v7[f"Q_{name}_fraction"] = v7[f"Q_{name}_W"] / v7.Q_total_W
    v7["Q_external_fraction"] = (
        v7.Q_natural_fraction + v7.Q_radiation_fraction + v7.Q_conduction_fraction
    )
    v7["external_exceeds_air"] = v7.Q_external_fraction > v7.Q_air_fraction
    v7.to_csv(output / "v7_energy_partition.csv", index=False)
    return v7


def fit_model(model: str, target: pd.DataFrame) -> tuple[np.ndarray, object]:
    if model == "M1":
        x0 = np.array([M1_GEOMETRY.nominal_hydraulic_diameter_m, 1.0])
        bounds = (np.array([2.0e-3, 0.05]), np.array([30.0e-3, 2.0]))
    else:
        x0 = np.array([80.0, 0.70])
        bounds = (np.array([0.5, 0.20]), np.array([200.0, 1.20]))

    def residual(parameters: np.ndarray) -> np.ndarray:
        result = []
        for velocity in sorted(CALIBRATION):
            solved = solve_structure(model, velocity, parameters, nodes=80)
            observed = target.loc[velocity, REGION_COLUMNS].to_numpy(float)
            predicted = np.array([solved[column] for column in REGION_COLUMNS])
            result.extend(predicted - observed)
        return np.asarray(result)

    fitted = least_squares(
        residual,
        x0,
        bounds=bounds,
        xtol=1e-9,
        ftol=1e-9,
        gtol=1e-9,
        max_nfev=120,
    )
    if not fitted.success:
        raise RuntimeError(f"{model} fit failed: {fitted.message}")
    return fitted.x, fitted


def near_boundary(values: np.ndarray, lower: np.ndarray, upper: np.ndarray) -> bool:
    fraction = (values - lower) / (upper - lower)
    return bool(np.any((fraction < 0.02) | (fraction > 0.98)))


def information_criteria(residual: np.ndarray, parameter_count: int) -> tuple[float, float]:
    n = len(residual)
    rss = max(float(residual @ residual), 1e-15)
    aic = n * np.log(rss / n) + 2 * parameter_count
    aicc = aic + 2 * parameter_count * (parameter_count + 1) / (
        n - parameter_count - 1
    )
    bic = n * np.log(rss / n) + parameter_count * np.log(n)
    return float(aicc), float(bic)


def evaluate_models(output: Path, target: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    predictions: list[dict[str, object]] = []
    mesh_rows: list[dict[str, object]] = []
    parameter_records: dict[str, object] = {}
    comparison_rows: list[dict[str, object]] = []

    m0 = pd.read_csv(ROOT / "outputs_v7/comsol_fig7_validation.csv")
    m0_parameters = np.array([30.0, 1.0])
    m0_bounds = (np.array([0.0, 0.0]), np.array([30.0, 1.0]))
    for _, row in m0.iterrows():
        predictions.append(
            {
                "model": "M0",
                "velocity_m_s": row.velocity_m_s,
                "mesh": "fine",
                **{column: row[column] for column in REGION_COLUMNS},
                "Tmax_C": row.Tmax_C,
                "Tmin_C": row.Tmin_C,
                "DeltaT_C": row.DeltaT_C,
                "hotspot_location_normalized": row.hotspot_location_normalized,
                "air_outlet_temperature_C": row.air_outlet_temperature_C,
                "Q_air_W": row.air_enthalpy_gain_W,
                "Q_natural_W": row.natural_convection_loss_W,
                "Q_radiation_W": row.radiation_loss_W,
                "Q_conduction_W": row.conduction_loss_W,
                "Q_total_W": row.total_heat_input_W,
                "energy_residual_relative": row.energy_residual_relative,
                "runtime_s": np.nan,
            }
        )
    parameter_records["M0"] = {
        "parameters": {"h_external_W_m2K": 30.0, "emissivity": 1.0},
        "bounds": {"h_external_W_m2K": [0.0, 30.0], "emissivity": [0.0, 1.0]},
        "parameters_at_boundary": True,
        "role": "V7 failed baseline only",
    }

    fits: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for model in ["M1", "M2"]:
        parameters, fitted = fit_model(model, target)
        if model == "M1":
            lower = np.array([2.0e-3, 0.05])
            upper = np.array([30.0e-3, 2.0])
            names = ["equivalent_hydraulic_diameter_m", "contact_area_multiplier"]
        else:
            lower = np.array([0.5, 0.20])
            upper = np.array([200.0, 1.20])
            names = ["h0_W_m2K_at_8m_s", "velocity_exponent_m"]
        fits[model] = (parameters, lower, upper)
        parameter_records[model] = {
            "parameters": dict(zip(names, map(float, parameters), strict=True)),
            "bounds": {
                name: [float(lo), float(hi)]
                for name, lo, hi in zip(names, lower, upper, strict=True)
            },
            "parameters_at_boundary": near_boundary(parameters, lower, upper),
            "calibration_cost": float(fitted.cost),
            "calibration_function_evaluations": int(fitted.nfev),
        }
        if model == "M1":
            parameter_records[model]["explicit_geometry"] = {
                "total_inlet_area_m2": PHYSICS.stack_inlet_area_m2,
                "channel_count": M1_GEOMETRY.channel_count,
                "channel_width_m": M1_GEOMETRY.channel_width_m,
                "channel_height_m": M1_GEOMETRY.channel_height_m,
                "rib_width_m": M1_GEOMETRY.rib_width_m,
                "total_pitch_width_m": M1_GEOMETRY.total_pitch_width_m,
                "wetted_perimeter_total_m": M1_GEOMETRY.wetted_perimeter_total_m,
                "nominal_hydraulic_diameter_m": M1_GEOMETRY.nominal_hydraulic_diameter_m,
                "nominal_contact_area_m2": M1_GEOMETRY.contact_area_cell_m2,
            }
        mesh_results: dict[tuple[float, int], dict[str, float | np.ndarray]] = {}
        for nodes, mesh_name in [(50, "coarse"), (100, "medium"), (200, "fine")]:
            for velocity in VELOCITIES:
                solved = solve_structure(model, velocity, parameters, nodes=nodes)
                mesh_results[(velocity, nodes)] = solved
                if mesh_name == "fine":
                    predictions.append(
                        {
                            "model": model,
                            "velocity_m_s": velocity,
                            "mesh": mesh_name,
                            **{column: solved[column] for column in REGION_COLUMNS},
                            "Tmax_C": solved["Tmax_C"],
                            "Tmin_C": solved["Tmin_C"],
                            "DeltaT_C": solved["DeltaT_C"],
                            "hotspot_location_normalized": solved["hotspot_location_normalized"],
                            "air_outlet_temperature_C": solved["air_outlet_temperature_C"],
                            "Q_air_W": solved["air_enthalpy_gain_W"],
                            "Q_natural_W": solved["natural_convection_loss_W"],
                            "Q_radiation_W": solved["radiation_loss_W"],
                            "Q_conduction_W": solved["conduction_loss_W"],
                            "Q_total_W": solved["total_heat_input_W"],
                            "energy_residual_relative": solved["energy_residual_relative"],
                            "runtime_s": solved["runtime_s"],
                            "h_internal_W_m2K": solved["h_internal_W_m2K"],
                            "h_natural_mean_W_m2K": solved["h_natural_mean_W_m2K"],
                            "reynolds": solved["reynolds"],
                            "nusselt": solved["nusselt"],
                            "effective_contact_area_m2": solved["effective_contact_area_m2"],
                        }
                    )
        for velocity in VELOCITIES:
            medium = mesh_results[(velocity, 100)]
            fine = mesh_results[(velocity, 200)]
            mesh_rows.append(
                {
                    "model": model,
                    "velocity_m_s": velocity,
                    "medium_fine_Tmax_relative_difference": abs(
                        float(medium["Tmax_C"]) - float(fine["Tmax_C"])
                    ) / abs(float(fine["Tmax_C"])),
                    "medium_fine_max_region_difference_C": max(
                        abs(float(medium[column]) - float(fine[column]))
                        for column in REGION_COLUMNS
                    ),
                    "fine_energy_residual_relative": fine["energy_residual_relative"],
                    "fine_mdot_identity_relative_error": fine["mdot_identity_relative_error"],
                }
            )

    prediction_frame = pd.DataFrame(predictions)
    prediction_frame.to_csv(output / "model_structure_predictions.csv", index=False)
    mesh_frame = pd.DataFrame(mesh_rows)
    mesh_frame.to_csv(output / "model_structure_mesh_convergence.csv", index=False)
    (output / "model_structure_calibration_record.json").write_text(
        json.dumps(parameter_records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    for model in ["M0", "M1", "M2"]:
        model_rows = prediction_frame[prediction_frame.model == model].set_index("velocity_m_s")
        calibration_residual = []
        validation_region_residual = []
        validation_tmax_residual = []
        validation_delta_residual = []
        per_speed_rmse = []
        for velocity in VELOCITIES:
            predicted = model_rows.loc[velocity]
            observed = target.loc[velocity]
            region_error = np.array(
                [float(predicted[column] - observed[column]) for column in REGION_COLUMNS]
            )
            if velocity in CALIBRATION:
                calibration_residual.extend(region_error)
            else:
                validation_region_residual.extend(region_error)
                validation_tmax_residual.append(float(predicted.Tmax_C - observed.Tmax_C))
                validation_delta_residual.append(float(predicted.DeltaT_C - observed.DeltaT_C))
                per_speed_rmse.append(float(np.sqrt(np.mean(region_error**2))))
        predicted_hotspots = model_rows.loc[sorted(VALIDATION), "hotspot_location_normalized"].to_numpy(float)
        observed_hotspots = target.loc[sorted(VALIDATION), "hotspot_location_normalized"].to_numpy(float)
        if np.ptp(predicted_hotspots) < 1e-12:
            hotspot_spearman = 0.0
            hotspot_trend_assessable = False
        else:
            hotspot_spearman = float(spearmanr(predicted_hotspots, observed_hotspots).statistic)
            hotspot_trend_assessable = np.isfinite(hotspot_spearman)
            if not hotspot_trend_assessable:
                hotspot_spearman = 0.0
        calibration_residual_array = np.asarray(calibration_residual)
        parameter_count = 2
        aicc, bic = information_criteria(calibration_residual_array, parameter_count)
        if model == "M0":
            mesh_ok = bool(m0.numerically_converged.all())
            max_tmax_mesh = float(m0.medium_fine_Tmax_relative_difference.max())
            max_region_mesh = float(m0.medium_fine_max_region_difference_C.max())
            max_energy = float(m0.energy_residual_relative.max())
            boundary = True
            runtime = np.nan
        else:
            model_mesh = mesh_frame[mesh_frame.model == model]
            max_tmax_mesh = float(model_mesh.medium_fine_Tmax_relative_difference.max())
            max_region_mesh = float(model_mesh.medium_fine_max_region_difference_C.max())
            max_energy = float(model_mesh.fine_energy_residual_relative.max())
            mesh_ok = bool(
                max_tmax_mesh < 0.01
                and max_region_mesh < 0.2
                and max_energy < 0.005
            )
            boundary = bool(parameter_records[model]["parameters_at_boundary"])
            runtime = float(model_rows.runtime_s.sum())
        q_air_fraction = model_rows.Q_air_W.to_numpy(float) / model_rows.Q_total_W.to_numpy(float)
        forced_monotonic = bool(np.all(np.diff(q_air_fraction) > 0.0))
        independent_rmse = float(
            np.sqrt(np.mean(np.asarray(validation_region_residual) ** 2))
        )
        delta_rmse = float(np.sqrt(np.mean(np.asarray(validation_delta_residual) ** 2)))
        gate = bool(
            max(per_speed_rmse) <= 1.5
            and independent_rmse <= 1.5
            and delta_rmse <= 1.5
            and hotspot_spearman >= 0.8
            and mesh_ok
            and not boundary
            and forced_monotonic
        )
        comparison_rows.append(
            {
                "model": model,
                "free_parameter_count": parameter_count,
                "calibration_region_RMSE_C": float(
                    np.sqrt(np.mean(calibration_residual_array**2))
                ),
                "independent_overall_region_RMSE_C": independent_rmse,
                "independent_worst_speed_region_RMSE_C": max(per_speed_rmse),
                "independent_Tmax_RMSE_C": float(
                    np.sqrt(np.mean(np.asarray(validation_tmax_residual) ** 2))
                ),
                "independent_DeltaT_RMSE_C": delta_rmse,
                "hotspot_trend_spearman": hotspot_spearman,
                "hotspot_trend_assessable": hotspot_trend_assessable,
                "AICc_calibration": aicc,
                "BIC_calibration": bic,
                "parameters_at_boundary": boundary,
                "maximum_medium_fine_Tmax_relative_difference": max_tmax_mesh,
                "maximum_medium_fine_region_difference_C": max_region_mesh,
                "maximum_energy_residual_relative": max_energy,
                "mesh_and_energy_passed": mesh_ok,
                "air_cooling_fraction_monotonically_increases": forced_monotonic,
                "fine_six_speed_runtime_s": runtime,
                "runtime_status": (
                    "not_recorded_in_V7_machine_readable_output"
                    if model == "M0"
                    else "measured_wall_clock"
                ),
                "all_acceptance_gates_passed": gate,
            }
        )
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(output / "model_structure_comparison.csv", index=False)
    passing = comparison[comparison.all_acceptance_gates_passed]
    selected = None if passing.empty else str(passing.sort_values("BIC_calibration").iloc[0].model)
    decision = {
        "selected_model": selected,
        "M1_passed": bool(comparison.set_index("model").loc["M1", "all_acceptance_gates_passed"]),
        "M2_passed": bool(comparison.set_index("model").loc["M2", "all_acceptance_gates_passed"]),
        "stop_rule_triggered": bool(selected is None),
        "downstream_ROM_or_Oracle_authorized": bool(selected is not None),
        "selected_model_sha256": None,
    }
    return comparison, decision


def area_hypotheses(output: Path, m1_parameters: dict[str, object]) -> pd.DataFrame:
    velocity = 8.0
    cell_mdot = PHYSICS.air_rho_kg_m3 * velocity * PHYSICS.cell_inlet_area_m2
    stack_mdot = PHYSICS.air_rho_kg_m3 * velocity * PHYSICS.stack_inlet_area_m2
    rows = [
        {
            "hypothesis": "H1",
            "interpretation": "3962 mm2 stack total; one representative cell channel scaled by 40",
            "modelled_inlet_area_m2": PHYSICS.cell_inlet_area_m2,
            "parallel_scale": 40,
            "modelled_heat_input_W": PHYSICS.heat_input_cell_W / 40,
            "total_heat_input_W": PHYSICS.heat_input_cell_W,
            "modelled_mdot_at_8m_s_kg_s": cell_mdot,
            "total_mdot_at_8m_s_kg_s": stack_mdot,
            "channel_count_per_model": 1,
            "contact_area_per_model_m2": PHYSICS.heat_area_cell_m2 / 40,
            "identifiability": "temperature-scale equivalent to H2 when every extensive quantity is scaled consistently",
        },
        {
            "hypothesis": "H2",
            "interpretation": "3962 mm2 is the total effective inlet of one aggregate CFD model",
            "modelled_inlet_area_m2": PHYSICS.stack_inlet_area_m2,
            "parallel_scale": 1,
            "modelled_heat_input_W": PHYSICS.heat_input_cell_W,
            "total_heat_input_W": PHYSICS.heat_input_cell_W,
            "modelled_mdot_at_8m_s_kg_s": stack_mdot,
            "total_mdot_at_8m_s_kg_s": stack_mdot,
            "channel_count_per_model": 1,
            "contact_area_per_model_m2": PHYSICS.heat_area_cell_m2,
            "identifiability": "temperature-scale equivalent to H1; Fig.7 alone cannot label stack versus aggregate model",
        },
        {
            "hypothesis": "H3",
            "interpretation": "fixed 3962 mm2 total inlet represented by forty discrete rectangular channels",
            "modelled_inlet_area_m2": PHYSICS.stack_inlet_area_m2,
            "parallel_scale": 1,
            "modelled_heat_input_W": PHYSICS.heat_input_cell_W,
            "total_heat_input_W": PHYSICS.heat_input_cell_W,
            "modelled_mdot_at_8m_s_kg_s": stack_mdot,
            "total_mdot_at_8m_s_kg_s": stack_mdot,
            "channel_count_per_model": M1_GEOMETRY.channel_count,
            "contact_area_per_model_m2": M1_GEOMETRY.contact_area_cell_m2,
            "identifiability": "structurally distinguishable through wetted perimeter, but channel count is not supplied by Fig.7",
        },
    ]
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "area_hypothesis_sensitivity.csv", index=False)
    return frame


def main() -> None:
    output = ROOT / "outputs_v8"
    output.mkdir(exist_ok=True)
    _, digitization, digitization_summary = digitization_reaudit(output)
    v7_energy = energy_partition(output)
    target = target_table()
    comparison, decision = evaluate_models(output, target)
    calibration_record = json.loads(
        (output / "model_structure_calibration_record.json").read_text(encoding="utf-8")
    )
    hypotheses = area_hypotheses(output, calibration_record["M1"])
    frozen_external = {
        "natural_convection": {
            "method": "Churchill-Chu vertical-plate correlation evaluated from the local surface temperature",
            "characteristic_length_m": PHYSICS.length_m,
            "reference_surface_temperature_C": 55.0,
            "reference_h_external_W_m2K": float(churchill_chu_h_W_m2K(55.0)),
            "calibrated_to_Fig7": False,
            "source": "https://www.sandia.gov/files/sierra/Aria_Users_5_20/simulation_setup/advanced_feat/correlation_heat_transfer_coeff_reference.html",
        },
        "radiation": {
            "emissivity": PHYSICS.emissivity,
            "basis": "fixed midpoint of the Industrial Graphite Engineering Handbook range 0.70-0.95",
            "calibrated_to_Fig7": False,
            "source": "https://nucleus.iaea.org/sites/graphiteknowledgebase/Meetings2/Old%20Meetings/2017/Background%20Info/GraphiteHandbook.pdf",
        },
    }
    (output / "frozen_external_heat_loss.json").write_text(
        json.dumps(frozen_external, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    summary = {
        "protected_v7_tag": "v7_physical_gate_failed",
        "scope": "Fig.7 cooling-structure discrepancy only; no ROM, MM-EKF or MPC executed",
        "digitization": digitization_summary,
        "frozen_external_heat_loss": frozen_external,
        "v7_energy_partition": {
            "four_m_s_air_fraction": float(v7_energy.loc[v7_energy.velocity_m_s == 4.0, "Q_air_fraction"].iloc[0]),
            "four_m_s_external_fraction": float(v7_energy.loc[v7_energy.velocity_m_s == 4.0, "Q_external_fraction"].iloc[0]),
            "external_exceeds_air_at_4m_s": bool(v7_energy.loc[v7_energy.velocity_m_s == 4.0, "external_exceeds_air"].iloc[0]),
        },
        "area_interpretation": {
            "recommended_structure": "H3 discrete channels with conserved total area, mass flow and heat",
            "stack_vs_single_aggregate_label_identifiable": False,
            "reason": "H1 and H2 are exactly scale-equivalent for temperature when all extensive quantities are scaled together",
            "hypothesis_count": int(len(hypotheses)),
        },
        "model_decision": decision,
        "structural_comparison": json.loads(
            comparison.to_json(orient="records")
        ),
        "final_answers": {
            "Fig7_digitization_reliable": bool(digitization_summary["digitization_reliable"]),
            "recommended_3962mm2_interpretation": "H3 structural representation; H1 versus H2 scale label remains unidentifiable",
            "error_mainly_internal_cooling_structure": True,
            "M1_or_M2_passed_independent_validation": bool(not decision["stop_rule_triggered"]),
            "continue_ROM_and_Oracle_MPC": bool(decision["downstream_ROM_or_Oracle_authorized"]),
            "retain_Thermal_article_as_high_fidelity_baseline": bool(not decision["stop_rule_triggered"]),
        },
        "pull_request_policy": {"keep_draft": True, "merge_allowed": False},
    }
    (output / "final_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary["model_decision"], indent=2))


if __name__ == "__main__":
    main()
