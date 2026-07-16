"""Calibrate and independently validate the 30/50/80-node solid-air model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, brentq
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from models.distributed_thermal_air_model import DistributedParameters, dynamic_reversal, steady_state  # noqa: E402
from models.five_node_model import FiveNodeParameters, simulate_open_loop  # noqa: E402
from models.model_registry import load_primary_model, write_model_runtime_record  # noqa: E402
from models.three_node_model import ThermalParameters, heat_generation_W, load_frozen_parameters  # noqa: E402


def mass_flow(velocity: float, fan: pd.DataFrame) -> float:
    duty = np.clip(100.0 * velocity / 12.0, 0.0, 100.0)
    return float(np.interp(duty, fan.duty_ratio_percent, fan.corotation_mass_flow_kg_s))


def profile_target(profiles: pd.DataFrame, velocity: float) -> tuple[np.ndarray, np.ndarray]:
    frame = profiles[profiles.velocity_m_s == velocity]
    coordinate = frame.flow_coordinate_inlet_0_outlet_1.to_numpy()
    temperature_C = frame.temperature_K.to_numpy() - 273.15
    keep = (coordinate >= 0.03) & (coordinate <= 0.97)
    return coordinate[keep], temperature_C[keep]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output_root = (args.output_dir or root / "outputs_v5").resolve()
    output = output_root / "distributed_model"
    output.mkdir(parents=True, exist_ok=True)
    model_data, primary_path, model_hash = load_primary_model(root)
    write_model_runtime_record(root, output, "distributed_model")
    config = yaml.safe_load((root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    external = output_root / "digitization" / "external"
    if not external.exists():
        # V6 reuses the frozen V5 digitization instead of regenerating or
        # silently changing the external calibration/validation targets.
        external = root / "outputs_v5" / "digitization" / "external"
    fan = pd.read_csv(external / "Fig2_digitized.csv")
    profiles = pd.read_csv(external / "Fig7_spatial_profiles.csv")
    fig10 = pd.read_csv(external / "Fig10_digitized.csv")

    calibration_velocities = [4.0, 8.0]
    validation_velocities = [5.0, 6.0, 10.0, 12.0]

    def make_parameters(log_values: np.ndarray, nodes: int) -> DistributedParameters:
        ua, external_loss, axial, heat_fraction = np.exp(log_values)
        return DistributedParameters(nodes, ua, external_loss, axial, heat_fraction)

    def residual(log_values: np.ndarray) -> np.ndarray:
        values: list[float] = []
        parameters = make_parameters(log_values, 50)
        for velocity in calibration_velocities:
            result = steady_state(parameters, mass_flow(velocity, fan), 1)
            coordinate, target = profile_target(profiles, velocity)
            prediction = np.interp(coordinate, np.linspace(0.0, 1.0, 50), result.solid_temperature_C)
            values.extend((prediction - target) / 1.5)
        return np.asarray(values)

    fit = least_squares(
        residual,
        np.log([50.0, 0.2, 8.0, 0.25]),
        bounds=(np.log([0.01, 0.005, 0.01, 0.05]), np.log([200.0, 20.0, 100.0, 2.0])),
        max_nfev=160,
        x_scale="jac",
    )
    fitted = make_parameters(fit.x, 50)
    parameter_record = {
        "model": "independent_distributed_solid_air_v5",
        "primary_model_sha256_for_downstream_comparison_only": model_hash,
        "calibrated_with_velocities_m_s": calibration_velocities,
        "independent_validation_velocities_m_s": validation_velocities,
        "values": {
            "total_UA_W_K": fitted.total_UA_W_K,
            "total_external_loss_W_K": fitted.total_external_loss_W_K,
            "axial_edge_conductance_W_K": fitted.axial_edge_conductance_W_K,
            "effective_heat_fraction": fitted.effective_heat_fraction,
            "paper_heat_flux_W_m2": fitted.heat_flux_W_m2,
            "paper_active_area_m2": fitted.active_area_m2,
            "paper_mass_times_graphite_cp_J_K": fitted.total_heat_capacity_J_K,
        },
        "explicit_equivalent_assumptions": [
            "Fig. 7 velocity is mapped to Fig. 2 duty as duty=100*velocity/12 because the paper gives no direct map.",
            "The effective heat fraction represents heat-flux-area and unmodelled external-loss ambiguity.",
            "The model is one-dimensional along the cooling-flow direction and does not reproduce full 3D ribs.",
        ],
    }
    (output / "fitted_parameters.json").write_text(json.dumps(parameter_record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    rows = []
    profiles_out = []
    for velocity in [4.0, 5.0, 6.0, 8.0, 10.0, 12.0]:
        result = steady_state(fitted, mass_flow(velocity, fan), 1)
        coordinate, target = profile_target(profiles, velocity)
        prediction = np.interp(coordinate, np.linspace(0.0, 1.0, 50), result.solid_temperature_C)
        rmse = float(np.sqrt(np.mean((prediction - target) ** 2)))
        rows.append({
            "velocity_m_s": velocity,
            "data_role": "calibration" if velocity in calibration_velocities else "independent_validation",
            "profile_rmse_C": rmse,
            "predicted_Tmax_C": float(result.solid_temperature_C.max()),
            "predicted_Tout_air_C": float(result.air_node_outlet_C[-1]),
            "energy_residual_W": result.energy_residual_W,
            "mass_flow_kg_s": mass_flow(velocity, fan),
        })
        profiles_out.append(pd.DataFrame({
            "velocity_m_s": velocity,
            "flow_coordinate": np.linspace(0, 1, 50),
            "predicted_solid_temperature_C": result.solid_temperature_C,
            "predicted_air_outlet_of_node_C": result.air_node_outlet_C,
        }))
    validation = pd.DataFrame(rows)
    validation.to_csv(output / "six_velocity_validation.csv", index=False)
    pd.concat(profiles_out, ignore_index=True).to_csv(output / "six_velocity_profiles.csv", index=False)

    mesh_rows = []
    reference = None
    for nodes in [30, 50, 80]:
        parameters = DistributedParameters(
            nodes, fitted.total_UA_W_K, fitted.total_external_loss_W_K,
            fitted.axial_edge_conductance_W_K * (nodes - 1) / 49.0,
            fitted.effective_heat_fraction,
        )
        result = steady_state(parameters, mass_flow(8.0, fan), 1)
        if nodes == 80:
            reference = result
        mesh_rows.append({"nodes": nodes, "Tmax_C": float(result.solid_temperature_C.max()), "Tout_air_C": float(result.air_node_outlet_C[-1]), "energy_residual_W": result.energy_residual_W})
    mesh = pd.DataFrame(mesh_rows)
    t80 = float(mesh.loc[mesh.nodes == 80, "Tmax_C"].iloc[0])
    mesh["Tmax_difference_from_80_C"] = (mesh.Tmax_C - t80).abs()
    mesh.to_csv(output / "mesh_convergence_30_50_80.csv", index=False)

    # Fig. 10 is kept independent.  The paper does not report its exact fan speed;
    # evaluate the pre-registered six speeds and report the best/entire range without refitting.
    switch_time = 18.0
    observed_difference = fig10.surface_1_C.to_numpy() - fig10.surface_2_C.to_numpy()
    observed_cross_candidates = np.flatnonzero(np.sign(observed_difference[1:]) != np.sign(observed_difference[:-1])) + 1
    observed_cross = float(fig10.time_s.iloc[observed_cross_candidates[0]]) if len(observed_cross_candidates) else None
    reversal_rows = []
    for velocity in [4.0, 5.0, 6.0, 8.0, 10.0, 12.0]:
        steady = steady_state(fitted, mass_flow(velocity, fan), 1)
        time, history, outlet, energy = dynamic_reversal(fitted, mass_flow(velocity, fan), steady.solid_temperature_C, switch_time, 60.0, 0.1, 1)
        difference = history[:, -1] - history[:, 0]
        crossings = np.flatnonzero(np.sign(difference[1:]) != np.sign(difference[:-1])) + 1
        crossing = float(time[crossings[0]]) if len(crossings) else None
        predicted_1 = np.interp(fig10.time_s, time, history[:, -1])
        predicted_2 = np.interp(fig10.time_s, time, history[:, 0])
        rmse = float(np.sqrt(np.mean(np.r_[predicted_1 - fig10.surface_1_C, predicted_2 - fig10.surface_2_C] ** 2)))
        reversal_rows.append({"velocity_m_s": velocity, "predicted_cross_time_s": crossing, "observed_cross_time_s": observed_cross, "two_surface_rmse_C": rmse, "maximum_energy_residual_W": float(np.max(np.abs(energy)))})
    reversal = pd.DataFrame(reversal_rows)
    reversal.to_csv(output / "fig10_independent_reversal_validation.csv", index=False)

    # Frozen five-node comparison at the paper heat flux.  This is an external audit,
    # not a refit of the five-node parameters.
    thermal = load_frozen_parameters(primary_path)
    five = FiveNodeParameters.from_three_node(thermal, np.asarray(config["model"]["coolant_node_temperatures_forward_C"], dtype=float))
    target_heat = fitted.total_heat_W
    current_for_heat = brentq(lambda current: heat_generation_W(current, 55.0, thermal, str(model_data["heat_model"])) - target_heat, 0.01, 44.0)
    five_rows = []
    for velocity in [4.0, 5.0, 6.0, 8.0, 10.0, 12.0]:
        duty = velocity / 12.0
        samples = 5001
        temperature, residual = simulate_open_loop(
            five,
            np.full(samples, current_for_heat),
            np.full(samples, duty),
            np.ones(samples, dtype=int),
            np.full(5, 25.0),
            1.0,
            str(model_data["heat_model"]),
        )
        high = steady_state(fitted, mass_flow(velocity, fan), 1)
        high_aggregated = np.asarray([np.mean(chunk) for chunk in np.array_split(high.solid_temperature_C, 5)])
        frozen = temperature[-1]
        five_rows.append({
            "velocity_m_s": velocity,
            "five_region_rmse_C": float(np.sqrt(np.mean((frozen - high_aggregated) ** 2))),
            "Tmax_error_C": float(np.max(frozen) - np.max(high.solid_temperature_C)),
            "mean_temperature_error_C": float(np.mean(frozen) - np.mean(high.solid_temperature_C)),
            "five_node_energy_residual_W": float(np.max(np.abs(residual))),
        })
    five_comparison = pd.DataFrame(five_rows)
    five_comparison.to_csv(output / "frozen_five_node_external_comparison.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    for axis, velocity in zip(axes.flat, [4, 5, 6, 8, 10, 12], strict=True):
        target_x, target_y = profile_target(profiles, velocity)
        prediction = steady_state(fitted, mass_flow(velocity, fan), 1).solid_temperature_C
        axis.plot(target_x, target_y, label="digitized Fig. 7")
        axis.plot(np.linspace(0, 1, 50), prediction, "--", label="50-node model")
        axis.set_title(f"{velocity} m/s"); axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    figure = axes[0, 0].figure
    figure.suptitle("Independent distributed-model validation")
    figure.tight_layout(); figure.savefig(output / "six_velocity_overlay.png", dpi=190); plt.close(figure)

    summary = {
        "primary_model_sha256": model_hash,
        "calibration_profile_rmse_C": float(validation[validation.data_role == "calibration"].profile_rmse_C.mean()),
        "independent_validation_profile_rmse_C": float(validation[validation.data_role == "independent_validation"].profile_rmse_C.mean()),
        "maximum_absolute_energy_residual_W": float(validation.energy_residual_W.abs().max()),
        "mesh_50_vs_80_Tmax_difference_C": float(mesh.loc[mesh.nodes == 50, "Tmax_difference_from_80_C"].iloc[0]),
        "fig10_observed_cross_time_s": observed_cross,
        "fig10_best_unfitted_velocity_m_s": float(reversal.loc[reversal.two_surface_rmse_C.idxmin(), "velocity_m_s"]),
        "fig10_best_unfitted_rmse_C": float(reversal.two_surface_rmse_C.min()),
        "frozen_five_node_mean_region_rmse_C": float(five_comparison.five_region_rmse_C.mean()),
        "space_validation_passed_within_palette_uncertainty": bool(validation[validation.data_role == "independent_validation"].profile_rmse_C.mean() <= 1.5),
        "reversal_validation_passed": bool(reversal.two_surface_rmse_C.min() <= 0.75),
        "five_node_external_validation_passed": bool(five_comparison.five_region_rmse_C.mean() <= 1.0),
    }
    (output / "distributed_model_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
