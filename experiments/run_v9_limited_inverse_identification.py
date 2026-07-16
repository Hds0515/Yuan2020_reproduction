"""V9 limited equivalent-parameter inverse identification for Fig. 7.

The calibration split remains frozen at 4 and 8 m/s.  Speeds 5, 6, 10 and
12 m/s are never used to select the parameter vector.  The task estimates an
admissible *equivalent* parameter set; it does not recover the authors' CAD.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.spatial.distance import cdist
from scipy.stats import qmc, spearmanr


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from comsol.v9_inverse.limited_inverse_model import (  # noqa: E402
    DEFAULT_INITIAL,
    LOWER_BOUNDS,
    PARAMETER_NAMES,
    UPPER_BOUNDS,
    result_row,
    solve_inverse_model,
)


OUTPUT = ROOT / "outputs_v9"
DOCS = ROOT / "docs_v9"
CALIBRATION = (4.0, 8.0)
VALIDATION = (5.0, 6.0, 10.0, 12.0)
VELOCITIES = CALIBRATION + VALIDATION
REGION_COLUMNS = (
    "T_inlet_region_C",
    "T_middle_region_C",
    "T_outlet_region_C",
)
FIT_NODES = 32
FINAL_NODES = 200
RANDOM_SEED = 20260716


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def target_table() -> pd.DataFrame:
    source = ROOT / "outputs_v5/digitization/external/Fig7_region_summary.csv"
    target = pd.read_csv(source)
    for column in [
        "T_inlet_region_K",
        "T_middle_region_K",
        "T_outlet_region_K",
        "Tmax_K",
        "Tmin_K",
    ]:
        target[column.replace("_K", "_C")] = target[column] - 273.15
    target["DeltaT_C"] = target["DeltaT_K"]
    return target.set_index("velocity_m_s")


def calibration_residual(parameters: np.ndarray, target: pd.DataFrame) -> np.ndarray:
    residuals: list[float] = []
    for velocity in CALIBRATION:
        solved = solve_inverse_model(velocity, parameters, nodes=FIT_NODES)
        residuals.extend(
            getattr(solved, column) - float(target.loc[velocity, column])
            for column in REGION_COLUMNS
        )
    return np.asarray(residuals, dtype=float)


def multistart_fit(target: pd.DataFrame, starts: int = 16) -> tuple[np.ndarray, pd.DataFrame]:
    sampler = qmc.LatinHypercube(d=4, seed=RANDOM_SEED)
    unit = sampler.random(starts - 1)
    initial_points = qmc.scale(unit, LOWER_BOUNDS, UPPER_BOUNDS)
    initial_points = np.vstack([DEFAULT_INITIAL, initial_points])
    rows: list[dict[str, object]] = []
    solutions: list[tuple[float, np.ndarray]] = []
    for index, initial in enumerate(initial_points):
        fitted = least_squares(
            lambda values: calibration_residual(values, target),
            initial,
            bounds=(LOWER_BOUNDS, UPPER_BOUNDS),
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
            max_nfev=140,
        )
        rss = float(fitted.fun @ fitted.fun)
        solutions.append((rss, fitted.x.copy()))
        row: dict[str, object] = {
            "start_index": index,
            "success": bool(fitted.success),
            "rss_C2": rss,
            "calibration_RMSE_C": float(np.sqrt(np.mean(fitted.fun**2))),
            "nfev": int(fitted.nfev),
        }
        row.update({name: float(value) for name, value in zip(PARAMETER_NAMES, fitted.x, strict=True)})
        rows.append(row)
    frame = pd.DataFrame(rows).sort_values("rss_C2").reset_index(drop=True)
    best = min(solutions, key=lambda item: item[0])[1]
    return best, frame


def evaluate(parameters: np.ndarray, target: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float | bool]]:
    rows: list[dict[str, float | str]] = []
    for velocity in sorted(VELOCITIES):
        solved = solve_inverse_model(velocity, parameters, nodes=FINAL_NODES)
        observed = target.loc[velocity]
        row = result_row(solved)
        row["data_role"] = "calibration" if velocity in CALIBRATION else "independent_validation"
        for column in REGION_COLUMNS + (
            "Tmax_C",
            "Tmin_C",
            "DeltaT_C",
            "hotspot_location_normalized",
        ):
            row[f"Fig7_{column}"] = float(observed[column])
            row[f"error_{column}"] = float(row[column] - observed[column])
        region_errors = np.array([row[f"error_{column}"] for column in REGION_COLUMNS])
        row["region_RMSE_C"] = float(np.sqrt(np.mean(region_errors**2)))
        rows.append(row)
    frame = pd.DataFrame(rows).sort_values("velocity_m_s")
    validation = frame[frame.data_role == "independent_validation"]
    calibration = frame[frame.data_role == "calibration"]
    all_validation_region_errors = validation[
        [f"error_{column}" for column in REGION_COLUMNS]
    ].to_numpy(float)
    predicted_hotspots = validation.hotspot_location_normalized.to_numpy(float)
    observed_hotspots = validation.Fig7_hotspot_location_normalized.to_numpy(float)
    if np.ptp(predicted_hotspots) < 1e-12:
        hotspot_spearman = 0.0
        hotspot_assessable = False
    else:
        value = spearmanr(predicted_hotspots, observed_hotspots).statistic
        hotspot_spearman = float(value) if np.isfinite(value) else 0.0
        hotspot_assessable = bool(np.isfinite(value))
    metrics: dict[str, float | bool] = {
        "calibration_region_RMSE_C": float(
            np.sqrt(np.mean(calibration[[f"error_{c}" for c in REGION_COLUMNS]].to_numpy(float) ** 2))
        ),
        "independent_overall_region_RMSE_C": float(np.sqrt(np.mean(all_validation_region_errors**2))),
        "independent_worst_speed_region_RMSE_C": float(validation.region_RMSE_C.max()),
        "independent_Tmax_RMSE_C": float(np.sqrt(np.mean(validation.error_Tmax_C.to_numpy(float) ** 2))),
        "independent_DeltaT_RMSE_C": float(np.sqrt(np.mean(validation.error_DeltaT_C.to_numpy(float) ** 2))),
        "hotspot_trend_spearman": hotspot_spearman,
        "hotspot_trend_assessable": hotspot_assessable,
        "maximum_energy_residual_relative": float(frame.energy_residual_relative.max()),
        "six_speed_runtime_s": float(frame.runtime_s.sum()),
    }
    metrics["region_overall_gate_passed"] = bool(metrics["independent_overall_region_RMSE_C"] <= 1.5)
    metrics["region_each_speed_gate_passed"] = bool(metrics["independent_worst_speed_region_RMSE_C"] <= 1.5)
    metrics["DeltaT_gate_passed"] = bool(metrics["independent_DeltaT_RMSE_C"] <= 3.0)
    metrics["hotspot_gate_passed"] = bool(hotspot_assessable and hotspot_spearman >= 0.8)
    metrics["energy_gate_passed"] = bool(metrics["maximum_energy_residual_relative"] <= 0.005)
    return frame, metrics


def profile_metrics(parameters: np.ndarray) -> tuple[pd.DataFrame, dict[str, float]]:
    source = pd.read_csv(ROOT / "outputs_v5/digitization/external/Fig7_spatial_profiles.csv")
    rows: list[dict[str, float]] = []
    for velocity in sorted(VELOCITIES):
        observed = source[source.velocity_m_s == velocity].sort_values(
            "flow_coordinate_inlet_0_outlet_1"
        )
        solved = solve_inverse_model(velocity, parameters, nodes=FINAL_NODES)
        predicted_C = np.interp(
            observed.flow_coordinate_inlet_0_outlet_1.to_numpy(float),
            solved.flow_coordinate,
            solved.surface_C,
        )
        observed_C = observed.temperature_K.to_numpy(float) - 273.15
        for coordinate, obs, pred in zip(
            observed.flow_coordinate_inlet_0_outlet_1,
            observed_C,
            predicted_C,
            strict=True,
        ):
            rows.append(
                {
                    "velocity_m_s": velocity,
                    "flow_coordinate": float(coordinate),
                    "Fig7_temperature_C": float(obs),
                    "M3_temperature_C": float(pred),
                    "error_C": float(pred - obs),
                }
            )
    frame = pd.DataFrame(rows)
    calibration = frame[frame.velocity_m_s.isin(CALIBRATION)]
    validation = frame[frame.velocity_m_s.isin(VALIDATION)]
    return frame, {
        "calibration_full_profile_RMSE_C": float(np.sqrt(np.mean(calibration.error_C**2))),
        "independent_full_profile_RMSE_C": float(np.sqrt(np.mean(validation.error_C**2))),
    }


def jacobian_audit(parameters: np.ndarray, target: pd.DataFrame) -> tuple[np.ndarray, dict[str, object]]:
    base = calibration_residual(parameters, target)
    jacobian = np.empty((len(base), 4), dtype=float)
    for index in range(4):
        step = max(abs(parameters[index]) * 1e-4, (UPPER_BOUNDS[index] - LOWER_BOUNDS[index]) * 1e-6)
        plus = parameters.copy()
        minus = parameters.copy()
        plus[index] = min(plus[index] + step, UPPER_BOUNDS[index])
        minus[index] = max(minus[index] - step, LOWER_BOUNDS[index])
        jacobian[:, index] = (
            calibration_residual(plus, target) - calibration_residual(minus, target)
        ) / (plus[index] - minus[index])
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    condition = float(singular_values[0] / max(singular_values[-1], 1e-15))
    fisher = jacobian.T @ jacobian
    covariance_shape = np.linalg.pinv(fisher)
    standard = np.sqrt(np.maximum(np.diag(covariance_shape), 1e-30))
    correlation = covariance_shape / np.outer(standard, standard)
    report = {
        "singular_values": singular_values.tolist(),
        "jacobian_condition_number": condition,
        "rank": int(np.linalg.matrix_rank(jacobian)),
        "parameter_correlation_matrix": correlation.tolist(),
        "maximum_absolute_parameter_correlation": float(
            np.max(np.abs(correlation - np.eye(4)))
        ),
    }
    return jacobian, report


def bootstrap(parameters: np.ndarray, target: pd.DataFrame, replicates: int = 30) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED + 1)
    sigma_C = 0.662  # V8 manual-vs-automatic regional RMSE.
    rows: list[dict[str, float | int | bool]] = []
    for replicate in range(replicates):
        perturbed = target.copy()
        for velocity in CALIBRATION:
            perturbed.loc[velocity, list(REGION_COLUMNS)] = (
                target.loc[velocity, list(REGION_COLUMNS)].to_numpy(float)
                + rng.normal(0.0, sigma_C, size=3)
            )
        starts = [parameters]
        reflected = parameters.copy()
        reflected[2] *= -1.0
        reflected[3] *= -1.0
        starts.append(np.clip(reflected, LOWER_BOUNDS, UPPER_BOUNDS))
        candidates = []
        for initial in starts:
            fit = least_squares(
                lambda values: calibration_residual(values, perturbed),
                initial,
                bounds=(LOWER_BOUNDS, UPPER_BOUNDS),
                max_nfev=60,
                xtol=1e-9,
                ftol=1e-9,
                gtol=1e-9,
            )
            candidates.append(fit)
        fitted = min(candidates, key=lambda item: float(item.fun @ item.fun))
        row: dict[str, float | int | bool] = {
            "replicate": replicate,
            "success": bool(fitted.success),
            "rss_C2": float(fitted.fun @ fitted.fun),
        }
        row.update({name: float(value) for name, value in zip(PARAMETER_NAMES, fitted.x, strict=True)})
        rows.append(row)
    return pd.DataFrame(rows)


def profile_likelihood(parameters: np.ndarray, target: pd.DataFrame, grid_count: int = 7) -> pd.DataFrame:
    """One-at-a-time local profile scan with other parameters frozen.

    This is intentionally bounded and deterministic.  It diagnoses flatness but
    is not presented as a full re-optimised likelihood confidence interval.
    """
    rows: list[dict[str, float | str]] = []
    for fixed_index, name in enumerate(PARAMETER_NAMES):
        span = UPPER_BOUNDS[fixed_index] - LOWER_BOUNDS[fixed_index]
        local_low = max(LOWER_BOUNDS[fixed_index], parameters[fixed_index] - 0.08 * span)
        local_high = min(UPPER_BOUNDS[fixed_index], parameters[fixed_index] + 0.08 * span)
        for fixed_value in np.linspace(local_low, local_high, grid_count):
            candidate = parameters.copy()
            candidate[fixed_index] = fixed_value
            residual = calibration_residual(candidate, target)
            rows.append({
                "profile_parameter": name,
                "fixed_value": float(fixed_value),
                "rss_C2": float(residual @ residual),
                **{parameter_name: float(value) for parameter_name, value in zip(PARAMETER_NAMES, candidate, strict=True)},
            })
    return pd.DataFrame(rows)

def global_sensitivity(parameters: np.ndarray, target: pd.DataFrame, samples: int = 160) -> pd.DataFrame:
    # Local admissible design box around the calibrated point; the output is a
    # rank-correlation screen, not a Sobol variance decomposition.
    half_width = 0.15 * (UPPER_BOUNDS - LOWER_BOUNDS)
    lower = np.maximum(LOWER_BOUNDS, parameters - half_width)
    upper = np.minimum(UPPER_BOUNDS, parameters + half_width)
    unit = qmc.LatinHypercube(d=4, seed=RANDOM_SEED + 2).random(samples)
    design = qmc.scale(unit, lower, upper)
    outputs = []
    for sample in design:
        values = []
        for velocity in CALIBRATION:
            solved = solve_inverse_model(velocity, sample, nodes=40)
            values.extend(getattr(solved, column) for column in REGION_COLUMNS)
        outputs.append(values)
    output_array = np.asarray(outputs)
    rows = []
    output_names = [f"{int(v)}m_s_{region}" for v in CALIBRATION for region in ("inlet", "middle", "outlet")]
    for parameter_index, parameter_name in enumerate(PARAMETER_NAMES):
        for output_index, output_name in enumerate(output_names):
            rho = spearmanr(design[:, parameter_index], output_array[:, output_index]).statistic
            rows.append(
                {
                    "parameter": parameter_name,
                    "output": output_name,
                    "spearman_rank_correlation": float(rho),
                    "absolute_correlation": float(abs(rho)),
                }
            )
    return pd.DataFrame(rows)


def parameter_interval_summary(bootstrap_frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for index, name in enumerate(PARAMETER_NAMES):
        values = bootstrap_frame[name].to_numpy(float)
        lower, median, upper = np.quantile(values, [0.025, 0.5, 0.975])
        bound_span = UPPER_BOUNDS[index] - LOWER_BOUNDS[index]
        rows.append(
            {
                "parameter": name,
                "bootstrap_median": float(median),
                "bootstrap_2p5": float(lower),
                "bootstrap_97p5": float(upper),
                "CI_width_fraction_of_prior_bound": float((upper - lower) / bound_span),
            }
        )
    return pd.DataFrame(rows)


def write_docs(summary: dict[str, object], intervals: pd.DataFrame, multistart: pd.DataFrame) -> None:
    metrics = summary["validation_metrics"]
    params = summary["frozen_equivalent_parameters"]
    interval_lines = "\n".join(
        f"| {row.parameter} | {row.bootstrap_median:.5g} | [{row.bootstrap_2p5:.5g}, {row.bootstrap_97p5:.5g}] | {row.CI_width_fraction_of_prior_bound:.3f} |"
        for row in intervals.itertuples()
    )
    report = f"""# V9 受限等效参数逆向辨识报告

## 结论

本阶段完成了四参数等效逆向辨识，但没有恢复、也不声称恢复原论文真实 CAD。冻结模型只辨识等效内部换热强度、速度指数、下游冷却偏置和热源偏置。

- 独立区域总体 RMSE：**{metrics['independent_overall_region_RMSE_C']:.3f} °C**；
- 独立最差单风速区域 RMSE：**{metrics['independent_worst_speed_region_RMSE_C']:.3f} °C**；
- 独立温差 RMSE：**{metrics['independent_DeltaT_RMSE_C']:.3f} °C**；
- 独立 Tmax RMSE：**{metrics['independent_Tmax_RMSE_C']:.3f} °C**；
- 热点趋势 Spearman：**{metrics['hotspot_trend_spearman']:.3f}**，可评估={metrics['hotspot_trend_assessable']}；
- 最大能量残差：**{metrics['maximum_energy_residual_relative']:.3e}**。

区域平均温度相较 V8 M2 明显改善，但最差风速和热点位置门槛未通过。因此该模型可作为“与 Fig.7 区域温度相容的等效稳态模型”，不可作为空间高保真真值，ROM、MM-EKF 和 Oracle MPC 仍不得启动。

## 冻结等效参数

| 参数 | 数值 |
|---|---:|
"""
    for name in PARAMETER_NAMES:
        report += f"| {name} | {params[name]:.8g} |\n"
    report += f"""

其中 `h0` 综合了无法分离的水力直径和有效接触面积；`beta` 与 `eta` 分别代表等效冷却和热源空间偏置，不对应唯一真实歧管或热源 CAD。

## Bootstrap 区间

| 参数 | 中位数 | 95% 区间 | 区间宽度/先验范围 |
|---|---:|---:|---:|
{interval_lines}

## 可辨识性

- 雅可比秩：{summary['identifiability']['rank']} / 4；
- 雅可比条件数：{summary['identifiability']['jacobian_condition_number']:.3e}；
- 最大绝对参数相关系数：{summary['identifiability']['maximum_absolute_parameter_correlation']:.3f}；
- 近最优多起点解数量：{summary['multistart']['near_optimal_solution_count']}。

多起点结果存在不同符号的冷却偏置/热源偏置组合，说明 Fig.7 区域数据无法唯一分离“流量不均匀”和“热源不均匀”。因此输出的是可容许等效参数集合，而非真实几何参数。

## 门控判定

| 门槛 | 结果 |
|---|---|
| 独立总体区域 RMSE ≤ 1.5 °C | {metrics['region_overall_gate_passed']} |
| 每个独立风速区域 RMSE ≤ 1.5 °C | {metrics['region_each_speed_gate_passed']} |
| 独立 DeltaT RMSE ≤ 3 °C | {metrics['DeltaT_gate_passed']} |
| 热点 Spearman ≥ 0.8 | {metrics['hotspot_gate_passed']} |
| 能量误差 ≤ 0.5% | {metrics['energy_gate_passed']} |

最终结论：`high_fidelity_authorized={summary['decision']['high_fidelity_authorized']}`。
"""
    (DOCS / "limited_equivalent_inverse_identification.md").write_text(report, encoding="utf-8")

    methods = """# V9 方法与边界

1. 4、8 m/s 仍为唯一校准速度；5、6、10、12 m/s 不参与参数选择。
2. 总入口面积、总热输入、材料参数、Churchill–Chu 自然对流和石墨发射率均继承 V8 冻结值。
3. 只辨识四个等效参数，不辨识通道数、真实歧管 CAD、逐通道流量、接触热阻或堆叠层数。
4. 使用 16 个拉丁超立方多起点、30 次噪声 Bootstrap、有限差分雅可比、局部单参数剖面扫描和局部 LHS 秩相关敏感性分析。
5. 校准最优解按校准 RSS 冻结，不使用独立验证结果在多个局部解之间挑选。
6. 所有温度结论基于一维稳态等效模型；本环境未运行 COMSOL，后续若用于工程仿真，需在本地 COMSOL 中实现同一参数化并做一次独立复算。
"""
    (DOCS / "inverse_identification_methodology.md").write_text(methods, encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)
    DOCS.mkdir(exist_ok=True)
    target = target_table()
    best, multistart = multistart_fit(target)
    predictions, validation_metrics = evaluate(best, target)
    profiles, profile_summary = profile_metrics(best)
    jacobian, identifiability = jacobian_audit(best, target)
    bootstrap_frame = bootstrap(best, target)
    interval_frame = parameter_interval_summary(bootstrap_frame)
    profile_frame = profile_likelihood(best, target)
    sensitivity_frame = global_sensitivity(best, target)

    multistart.to_csv(OUTPUT / "multistart_solutions.csv", index=False)
    predictions.to_csv(OUTPUT / "inverse_model_predictions.csv", index=False)
    profiles.to_csv(OUTPUT / "full_profile_predictions.csv", index=False)
    bootstrap_frame.to_csv(OUTPUT / "bootstrap_parameters.csv", index=False)
    interval_frame.to_csv(OUTPUT / "parameter_intervals.csv", index=False)
    profile_frame.to_csv(OUTPUT / "profile_likelihood.csv", index=False)
    sensitivity_frame.to_csv(OUTPUT / "global_sensitivity_rank_correlations.csv", index=False)
    pd.DataFrame(jacobian, columns=PARAMETER_NAMES).to_csv(OUTPUT / "calibration_jacobian.csv", index=False)

    near_optimal = multistart[multistart.rss_C2 <= multistart.rss_C2.iloc[0] + 0.25]
    at_boundary = bool(
        np.any((best - LOWER_BOUNDS) / (UPPER_BOUNDS - LOWER_BOUNDS) < 0.02)
        or np.any((UPPER_BOUNDS - best) / (UPPER_BOUNDS - LOWER_BOUNDS) < 0.02)
    )
    all_gates = bool(
        validation_metrics["region_overall_gate_passed"]
        and validation_metrics["region_each_speed_gate_passed"]
        and validation_metrics["DeltaT_gate_passed"]
        and validation_metrics["hotspot_gate_passed"]
        and validation_metrics["energy_gate_passed"]
        and not at_boundary
    )
    parameter_mapping = {
        "h0_W_m2K_at_8m_s": "lumped effective internal heat transfer; hydraulic diameter and contact area are not separately identifiable",
        "velocity_exponent_m": "velocity sensitivity of effective forced convection",
        "downstream_cooling_bias_beta": "equivalent downstream cooling imbalance; not a unique manifold geometry",
        "heat_source_skew_eta": "equivalent linear heat-source skew; not a measured electrochemical heat map",
    }
    frozen_parameters = {name: float(value) for name, value in zip(PARAMETER_NAMES, best, strict=True)}
    frozen_path = OUTPUT / "frozen_equivalent_parameters.json"
    frozen_payload = {
        "model": "M3_limited_equivalent_inverse_model",
        "parameters": frozen_parameters,
        "bounds": {
            name: [float(lower), float(upper)]
            for name, lower, upper in zip(PARAMETER_NAMES, LOWER_BOUNDS, UPPER_BOUNDS, strict=True)
        },
        "calibration_velocities_m_s": list(CALIBRATION),
        "independent_validation_velocities_m_s": list(VALIDATION),
        "parameter_interpretation": parameter_mapping,
        "not_original_CAD": True,
    }
    frozen_path.write_text(json.dumps(frozen_payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    summary: dict[str, object] = {
        "scope": "limited equivalent inverse identification; no claim of original geometry recovery",
        "source_target": {
            "path": "outputs_v5/digitization/external/Fig7_region_summary.csv",
            "sha256": sha256(ROOT / "outputs_v5/digitization/external/Fig7_region_summary.csv"),
            "digitization_manual_audit_RMSE_C": 0.662,
        },
        "frozen_equivalent_parameters": frozen_parameters,
        "frozen_parameter_file_sha256": sha256(frozen_path),
        "validation_metrics": {**validation_metrics, **profile_summary},
        "identifiability": identifiability,
        "multistart": {
            "start_count": int(len(multistart)),
            "best_calibration_RSS_C2": float(multistart.rss_C2.iloc[0]),
            "near_optimal_solution_count": int(len(near_optimal)),
            "near_optimal_threshold_delta_RSS_C2": 0.25,
        },
        "bootstrap": {
            "replicates": int(len(bootstrap_frame)),
            "noise_sigma_C": 0.662,
            "interval_file": "outputs_v9/parameter_intervals.csv",
        },
        "decision": {
            "equivalent_regional_temperature_model_improved": bool(
                validation_metrics["independent_overall_region_RMSE_C"] < 2.845862878019487
            ),
            "actual_geometry_recovered": False,
            "high_fidelity_authorized": all_gates,
            "ROM_MM_EKF_Oracle_MPC_authorized": all_gates,
            "use_allowed": "steady regional-temperature sensitivity and bounded equivalent simulation only",
            "use_prohibited": "original CAD claim, hotspot truth, high-fidelity ROM training or decisive control validation",
        },
    }
    (OUTPUT / "final_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_docs(summary, interval_frame, multistart)
    print(json.dumps(summary["decision"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
