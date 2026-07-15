"""Weighted multistart identification using calibration data only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from models.three_node_model import ThermalParameters, simulate  # noqa: E402

BASE_NAMES = ["C_th_J_K", "K_cool_W_K", "K_node_W_K", "Kp_1_K", "Ki_1_Ks"]
HEAT_NAMES = ["heat_a0_W", "heat_a1_W_A", "heat_a2_W_A2"]


def parameter_names(heat_model: str) -> list[str]:
    return BASE_NAMES + (HEAT_NAMES if heat_model == "B_quadratic_current" else [])


def vector_to_parameters(x: np.ndarray, heat_model: str) -> ThermalParameters:
    values = dict(zip(parameter_names(heat_model), map(float, x), strict=True))
    return ThermalParameters(**values)


def load_cases(root: Path) -> dict[str, object]:
    fig14 = pd.read_csv(root / "digitization" / "regenerated" / "Fig14_digitized.csv")
    fig15 = pd.read_csv(root / "digitization" / "regenerated" / "Fig15_digitized.csv")
    fig16 = pd.read_csv(root / "digitization" / "regenerated" / "Fig16_digitized.csv")
    fig14_target = fig14[
        ["paper_simulation_T1_C", "paper_simulation_T2_C", "paper_simulation_T3_C"]
    ].to_numpy()
    # State ordering follows the physical channel direction used by the V2 model.
    fig16_target = fig16[
        ["paper_simulation_T2_C", "paper_simulation_middle_C", "paper_simulation_T1_C"]
    ].to_numpy()
    return {
        "fig14_time": fig14["time_s"].to_numpy(),
        "fig14_current": np.full(len(fig14), 10.0),
        "fig14_target": fig14_target,
        "fig16_time": fig16["time_s"].to_numpy(),
        "fig16_current": np.interp(fig16["time_s"], fig15["time_s"], fig15["paper_simulation_current_A"]),
        "fig16_target": fig16_target,
    }


def simulate_both(
    x: np.ndarray,
    heat_model: str,
    cases: dict[str, object],
    config: dict[str, object],
) -> tuple[object, object]:
    parameters = vector_to_parameters(x, heat_model)
    model = config["model"]
    common = dict(
        dt_s=float(config["project"]["dt_s"]),
        fan_enable_temperature_C=float(model["fan_enable_temperature_C"]),
        reference_temperature_C=float(model["reference_temperature_C"]),
        smc_deadband_C=float(model["smc_deadband_C"]),
        smc_update_period_s=float(model["smc_update_period_s"]),
        heat_model=heat_model,
        coolant_node_temperatures_forward_C=np.asarray(model["coolant_node_temperatures_forward_C"]),
    )
    fig14 = simulate(
        parameters,
        cases["fig14_current"],
        cases["fig14_target"][0],
        bidirectional=False,
        initial_airflow_direction=1,
        fan_initially_enabled=False,
        **common,
    )
    fig16 = simulate(
        parameters,
        cases["fig16_current"],
        cases["fig16_target"][0],
        bidirectional=True,
        initial_airflow_direction=-1,
        fan_initially_enabled=True,
        **common,
    )
    return fig14, fig16


def residual_vector(
    x: np.ndarray,
    heat_model: str,
    cases: dict[str, object],
    config: dict[str, object],
    targets: tuple[np.ndarray, np.ndarray] | None = None,
) -> np.ndarray:
    fig14, fig16 = simulate_both(x, heat_model, cases, config)
    target14, target16 = targets or (cases["fig14_target"], cases["fig16_target"])
    mask14 = (cases["fig14_time"] >= 0.0) & (cases["fig14_time"] <= 1200.0)
    mask16 = (cases["fig16_time"] >= 0.0) & (cases["fig16_time"] <= 700.0)
    # Three-pixel line-width uncertainty is used as a measurement weight.
    sigma14_C = 3.0 * 20.0 / 1073.0
    sigma16_C = 3.0 * 20.0 / 1153.0
    return np.concatenate(
        [
            ((fig14.temperature_C[mask14] - target14[mask14]) / sigma14_C).ravel(),
            ((fig16.temperature_C[mask16] - target16[mask16]) / sigma16_C).ravel(),
        ]
    )


def multistart_fit(
    heat_model: str,
    cases: dict[str, object],
    config: dict[str, object],
    rng: np.random.Generator,
) -> tuple[object, pd.DataFrame]:
    names = parameter_names(heat_model)
    bounds_cfg = config["identification"]["parameter_bounds"]
    lower = np.asarray([bounds_cfg[n][0] for n in names], dtype=float)
    upper = np.asarray([bounds_cfg[n][1] for n in names], dtype=float)
    if heat_model == "A_thermoneutral_proxy":
        preferred = np.asarray([4819.765, 31.2555, 7.08419, 0.06918, 0.0014268])
    else:
        preferred = np.asarray([4819.765, 31.2555, 7.08419, 0.06918, 0.0014268, 20.0, 20.0, 0.25])
    starts = [np.clip(preferred, lower, upper)]
    count = int(config["identification"]["multistart_count"])
    for _ in range(count - 1):
        # Log sampling is used for positive parameters spanning several decades.
        starts.append(np.exp(rng.uniform(np.log(np.maximum(lower, 1e-8)), np.log(upper))))

    results = []
    rows = []
    for index, start in enumerate(starts):
        result = least_squares(
            residual_vector,
            start,
            args=(heat_model, cases, config),
            bounds=(lower, upper),
            x_scale="jac",
            max_nfev=180,
            ftol=1e-7,
            xtol=1e-7,
            gtol=1e-7,
        )
        results.append(result)
        rows.append(
            {
                "heat_model": heat_model,
                "start_index": index,
                "success": bool(result.success),
                "cost": float(result.cost),
                "optimality": float(result.optimality),
                "nfev": int(result.nfev),
                **{name: float(value) for name, value in zip(names, result.x, strict=True)},
            }
        )
        print(f"{heat_model} multistart {index + 1}/{count}: cost={result.cost:.6g}", flush=True)
    best = min(results, key=lambda item: item.cost)
    return best, pd.DataFrame(rows)


def rmse_metrics(result14: object, result16: object, cases: dict[str, object]) -> dict[str, float]:
    calibration16 = cases["fig16_time"] <= 700.0
    validation16 = cases["fig16_time"] >= 700.0
    return {
        "fig14_calibration_rmse_C": float(
            np.sqrt(np.mean((result14.temperature_C - cases["fig14_target"]) ** 2))
        ),
        "fig16_calibration_rmse_C": float(
            np.sqrt(np.mean((result16.temperature_C[calibration16] - cases["fig16_target"][calibration16]) ** 2))
        ),
        "fig16_independent_validation_rmse_C": float(
            np.sqrt(np.mean((result16.temperature_C[validation16] - cases["fig16_target"][validation16]) ** 2))
        ),
    }


def covariance_from_jacobian(best: object) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    jac = np.asarray(best.jac)
    dof = max(jac.shape[0] - jac.shape[1], 1)
    residual_variance = 2.0 * best.cost / dof
    covariance = np.linalg.pinv(jac.T @ jac, rcond=1e-12) * residual_variance
    standard = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denom = np.outer(standard, standard)
    correlation = np.divide(covariance, denom, out=np.zeros_like(covariance), where=denom > 0)
    singular = np.linalg.svd(jac, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
    return covariance, correlation, singular, condition


def block_resample(residual: np.ndarray, rng: np.random.Generator, block: int = 30) -> np.ndarray:
    n = len(residual)
    starts = rng.integers(0, max(n - block + 1, 1), size=int(np.ceil(n / block)))
    return np.concatenate([residual[s : s + block] for s in starts], axis=0)[:n]


def bootstrap(
    best: object,
    heat_model: str,
    cases: dict[str, object],
    config: dict[str, object],
    rng: np.random.Generator,
) -> pd.DataFrame:
    names = parameter_names(heat_model)
    bounds_cfg = config["identification"]["parameter_bounds"]
    lower = np.asarray([bounds_cfg[n][0] for n in names], dtype=float)
    upper = np.asarray([bounds_cfg[n][1] for n in names], dtype=float)
    fitted14, fitted16 = simulate_both(best.x, heat_model, cases, config)
    residual14 = cases["fig14_target"] - fitted14.temperature_C
    cal16 = cases["fig16_time"] <= 700.0
    residual16 = cases["fig16_target"][cal16] - fitted16.temperature_C[cal16]
    rows = []
    for replicate in range(int(config["identification"]["bootstrap_count"])):
        synthetic14 = fitted14.temperature_C + block_resample(residual14, rng)
        synthetic16 = cases["fig16_target"].copy()
        synthetic16[cal16] = fitted16.temperature_C[cal16] + block_resample(residual16, rng)
        fit = least_squares(
            residual_vector,
            best.x,
            args=(heat_model, cases, config, (synthetic14, synthetic16)),
            bounds=(lower, upper),
            x_scale="jac",
            max_nfev=100,
        )
        rows.append(
            {
                "replicate": replicate,
                "success": bool(fit.success),
                "cost": float(fit.cost),
                **{name: float(value) for name, value in zip(names, fit.x, strict=True)},
            }
        )
        print(f"bootstrap {replicate + 1}/{config['identification']['bootstrap_count']}: cost={fit.cost:.6g}", flush=True)
    return pd.DataFrame(rows)


def save_plots(
    root: Path,
    names: list[str],
    best: object,
    result14: object,
    result16: object,
    cases: dict[str, object],
    correlation: np.ndarray,
    singular: np.ndarray,
) -> None:
    out = root / "outputs" / "three_node"
    out.mkdir(parents=True, exist_ok=True)
    colors = ["tab:red", "tab:green", "tab:blue"]
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, color in enumerate(colors):
        ax.plot(cases["fig14_time"], cases["fig14_target"][:, i], color=color, alpha=0.45)
        ax.plot(cases["fig14_time"], result14.temperature_C[:, i], "--", color=color)
    ax.set(xlabel="Time (s)", ylabel="Temperature (°C)", title="Fig14 calibration: image trace (solid) vs v3 (dashed)")
    fig.tight_layout(); fig.savefig(out / "Fig14_overlay.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, color in enumerate(colors):
        ax.plot(cases["fig16_time"], cases["fig16_target"][:, i], color=color, alpha=0.45)
        ax.plot(cases["fig16_time"], result16.temperature_C[:, i], "--", color=color)
    ax.axvspan(700, 1200, color="0.9", label="independent validation")
    ax.set(xlabel="Time (s)", ylabel="Temperature (°C)", title="Fig16 calibration and untouched validation")
    ax.legend(); fig.tight_layout(); fig.savefig(out / "Fig16_overlay.png", dpi=180); plt.close(fig)

    mask = cases["fig16_time"] >= 700
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for i, color in enumerate(colors):
        ax.plot(cases["fig16_time"][mask], cases["fig16_target"][mask, i], color=color, alpha=0.45)
        ax.plot(cases["fig16_time"][mask], result16.temperature_C[mask, i], "--", color=color)
    ax.set(xlabel="Time (s)", ylabel="Temperature (°C)", title="Fig16 independent validation only (700-1200 s)")
    fig.tight_layout(); fig.savefig(out / "Fig16_validation_only.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=False)
    axes[0].plot(cases["fig14_time"], result14.temperature_C - cases["fig14_target"], lw=0.7)
    axes[0].set_title("Fig14 residuals"); axes[0].set_ylabel("°C")
    axes[1].plot(cases["fig16_time"], result16.temperature_C - cases["fig16_target"], lw=0.7)
    axes[1].axvline(700, color="k", ls="--"); axes[1].set_title("Fig16 residuals"); axes[1].set(xlabel="Time (s)", ylabel="°C")
    fig.tight_layout(); fig.savefig(out / "residuals.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.semilogy(cases["fig16_time"], np.maximum(np.abs(result16.total_energy_residual_W), 1e-18))
    ax.set(xlabel="Time (s)", ylabel="Absolute residual (W)", title="Discrete energy balance residual")
    fig.tight_layout(); fig.savefig(out / "energy_balance.png", dpi=180); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    image = axes[0].imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm")
    axes[0].set_xticks(range(len(names)), names, rotation=75, ha="right", fontsize=7)
    axes[0].set_yticks(range(len(names)), names, fontsize=7); axes[0].set_title("Local parameter correlation")
    fig.colorbar(image, ax=axes[0], fraction=0.046)
    axes[1].semilogy(range(1, len(singular) + 1), singular, "o-")
    axes[1].set(xlabel="Singular value index", ylabel="Weighted Jacobian singular value", title="Local sensitivity spectrum")
    fig.tight_layout(); fig.savefig(out / "parameter_identifiability.png", dpi=180); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    root = args.root.resolve()
    config = yaml.safe_load((root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    rng = np.random.default_rng(int(config["project"]["seed"]))
    cases = load_cases(root)
    out = root / "identification_v3"; out.mkdir(parents=True, exist_ok=True)

    candidate_rows = []
    best_by_model = {}
    multistarts = []
    for heat_model in ("A_thermoneutral_proxy", "B_quadratic_current"):
        best, starts = multistart_fit(heat_model, cases, config, rng)
        result14, result16 = simulate_both(best.x, heat_model, cases, config)
        metrics = rmse_metrics(result14, result16, cases)
        physical = True
        if heat_model == "B_quadratic_current":
            p = vector_to_parameters(best.x, heat_model)
            currents = np.linspace(0, 35, 100)
            heats = p.heat_a0_W + p.heat_a1_W_A * currents + p.heat_a2_W_A2 * currents**2
            physical = bool(np.all(heats >= -1e-9) and np.all(np.diff(heats) >= -1e-9))
        candidate_rows.append({"heat_model": heat_model, "physically_admissible": physical, **metrics})
        best_by_model[heat_model] = (best, result14, result16)
        multistarts.append(starts)
    comparison = pd.DataFrame(candidate_rows)
    admissible = comparison[comparison["physically_admissible"]]
    selected_model = str(admissible.sort_values("fig16_independent_validation_rmse_C").iloc[0]["heat_model"])
    best, result14, result16 = best_by_model[selected_model]
    names = parameter_names(selected_model)
    covariance, correlation, singular, condition = covariance_from_jacobian(best)
    boot = bootstrap(best, selected_model, cases, config, rng)
    standard = np.sqrt(np.maximum(np.diag(covariance), 0.0))

    units = {
        "C_th_J_K": "J/K", "K_cool_W_K": "W/K", "K_node_W_K": "W/K",
        "Kp_1_K": "1/K", "Ki_1_Ks": "1/(K s)", "heat_a0_W": "W",
        "heat_a1_W_A": "W/A", "heat_a2_W_A2": "W/A²",
    }
    identified = pd.DataFrame(
        {
            "parameter": names,
            "identified_value": best.x,
            "unit": [units[n] for n in names],
            "local_95pct_lower": best.x - 1.96 * standard,
            "local_95pct_upper": best.x + 1.96 * standard,
            "bootstrap_2p5": [boot[n].quantile(0.025) for n in names],
            "bootstrap_97p5": [boot[n].quantile(0.975) for n in names],
        }
    )
    identified.to_csv(out / "identified_parameters.csv", index=False)
    pd.DataFrame(covariance, index=names, columns=names).to_csv(out / "parameter_covariance.csv")
    pd.DataFrame(correlation, index=names, columns=names).to_csv(out / "parameter_correlation.csv")
    boot.to_csv(out / "profile_or_bootstrap_results.csv", index=False)
    comparison.to_csv(out / "heat_model_comparison.csv", index=False)
    pd.concat(multistarts, ignore_index=True).to_csv(out / "multistart_results.csv", index=False)

    split = {
        "calibration": [
            {"source": "digitization/regenerated/Fig14_digitized.csv", "time_s": [0, 1200]},
            {"source": "digitization/regenerated/Fig16_digitized.csv", "time_s": [0, 700]},
        ],
        "independent_validation": [
            {"source": "digitization/regenerated/Fig16_digitized.csv", "time_s": [700, 1200]}
        ],
        "validation_used_in_parameter_optimization": False,
        "validation_used_only_for_predeclared_heat_model_A_B_comparison": True,
    }
    (out / "calibration_validation_split.json").write_text(json.dumps(split, indent=2) + "\n", encoding="utf-8")
    metrics = rmse_metrics(result14, result16, cases)
    max_residual = float(max(np.max(np.abs(result14.total_energy_residual_W)), np.max(np.abs(result16.total_energy_residual_W))))
    max_heat = float(max(np.max(np.sum(result14.node_heat_generation_W, axis=1)), np.max(np.sum(result16.node_heat_generation_W, axis=1))))
    metrics.update(
        {
            "selected_heat_model": selected_model,
            "weighted_jacobian_condition_number": condition,
            "maximum_total_energy_residual_W": max_residual,
            "maximum_relative_energy_residual": max_residual / max(max_heat, 1e-12),
            "K_cool_times_Kp": float(best.x[names.index("K_cool_W_K")] * best.x[names.index("Kp_1_K")]),
            "K_cool_times_Ki": float(best.x[names.index("K_cool_W_K")] * best.x[names.index("Ki_1_Ks")]),
        }
    )
    (out / "identification_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    payload = {
        "model": "Yuan2020_three_node_control_oriented_v3",
        "heat_model": selected_model,
        "values": {name: float(value) for name, value in zip(names, best.x, strict=True)},
        "paper_reported_or_directly_derived": {
            "reference_temperature_C": 55.0,
            "fan_enable_temperature_C": 53.0,
            "smc_deadband_C": 0.5,
            "smc_update_period_s": 10.0,
            "coolant_node_temperatures_forward_C": [24.8, 30.05, 35.3],
            "number_of_cells": 40,
        },
        "assumptions": {
            "polarization_curve_is_proxy": True,
            "radiation_neglected": True,
            "uniform_node_heat_generation": True,
        },
        "metrics": metrics,
        "identifiability_note": "K_cool, Kp and Ki are structurally entangled when fan duty is unavailable; product parameters are more defensible than individual values.",
    }
    serialized = json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")
    payload["sha256_without_hash_field"] = hashlib.sha256(serialized).hexdigest()
    (out / "frozen_parameters_v3.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    save_plots(root, names, best, result14, result16, cases, correlation, singular)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
