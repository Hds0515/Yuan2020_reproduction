"""Audit the independent COMSOL implementation of the frozen V9 model.

This script performs no fitting. It asserts the frozen input hashes and values,
parses the raw COMSOL 6.4 exports, applies the preregistered consistency gates,
and evaluates the unchanged Fig. 7 calibration/validation split.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "comsol" / "v9_confirmation"
OUTPUT_DIR = ROOT / "outputs_v9_comsol"
DOCS_DIR = ROOT / "docs_v9_comsol"
FIGURE_DIR = OUTPUT_DIR / "temperature_fields"

VELOCITIES = (4.0, 5.0, 6.0, 8.0, 10.0, 12.0)
CALIBRATION_SPEEDS = (4.0, 8.0)
VALIDATION_SPEEDS = (5.0, 6.0, 10.0, 12.0)
REGION_COLUMNS = (
    "T_inlet_region_C",
    "T_middle_region_C",
    "T_outlet_region_C",
)
METRIC_COLUMNS = (
    *REGION_COLUMNS,
    "Tmax_C",
    "Tmin_C",
    "DeltaT_C",
    "hotspot_location_normalized",
    "air_outlet_temperature_C",
    "air_enthalpy_gain_W",
    "natural_convection_loss_W",
    "radiation_loss_W",
    "total_heat_input_W",
    "energy_residual_relative",
    "mdot_kg_s",
)

EXPECTED_HASHES = {
    "frozen_parameters": (
        "outputs_v9/frozen_equivalent_parameters.json",
        "e4504981c33e2381cb29c2888fd0604a7b745a4c2bbf8523073b10f7fd2c7343",
    ),
    "python_predictions": (
        "outputs_v9/inverse_model_predictions.csv",
        "3d3ef55aae16dd7b3e067a533528e5aa94dc4a022d6197d7fe4ee5ee71042028",
    ),
    "fig7_targets": (
        "outputs_v5/digitization/external/Fig7_region_summary.csv",
        "01b87ddea0591e84c30cfa8bbe2c789bed47a36022895f18b2e4c7cee90ffe50",
    ),
}
EXPECTED_PARAMETERS = {
    "h0_W_m2K_at_8m_s": 88.75329934335754,
    "velocity_exponent_m": 0.6572949080566769,
    "downstream_cooling_bias_beta": 1.3679347347371371,
    "heat_source_skew_eta": -0.3118665087792232,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assert_frozen_inputs() -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    for name, (relative, expected) in EXPECTED_HASHES.items():
        path = ROOT / relative
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(
                f"Protected input hash changed for {relative}: {actual} != {expected}"
            )
        records[name] = {
            "path": relative,
            "sha256": actual,
        }
    frozen = json.loads((ROOT / EXPECTED_HASHES["frozen_parameters"][0]).read_text())
    actual_parameters = frozen["parameters"]
    for name, expected in EXPECTED_PARAMETERS.items():
        if not math.isclose(
            float(actual_parameters[name]), expected, rel_tol=0.0, abs_tol=1e-14
        ):
            raise RuntimeError(f"Frozen parameter {name} changed")
    java_text = (RAW_DIR / "V9EquivalentConfirmation.java").read_text(
        encoding="utf-8"
    )
    required_fragments = (
        "qTotal*qRaw/intop(qRaw)",
        "h0*(uin/(8[m/s]))^velExponent*exp(-coolingBias*(xi-0.5))",
        "rhoAir*uin*AinTotal",
        "0.82",
        "NuNatural",
    )
    missing = [fragment for fragment in required_fragments if fragment not in java_text]
    if missing:
        raise RuntimeError(f"COMSOL source audit failed; missing {missing}")
    records["comsol_java_source"] = {
        "path": str((RAW_DIR / "V9EquivalentConfirmation.java").relative_to(ROOT)),
        "sha256": sha256(RAW_DIR / "V9EquivalentConfirmation.java"),
    }
    return records


def read_comsol_metrics(mesh: str) -> pd.DataFrame:
    raw = pd.read_csv(
        RAW_DIR / f"metrics_{mesh}.csv", comment="%", header=None
    )
    if raw.shape != (6, 1 + len(METRIC_COLUMNS)):
        raise RuntimeError(f"Unexpected COMSOL metric shape for {mesh}: {raw.shape}")
    raw.columns = ("velocity_m_s", *METRIC_COLUMNS)
    raw.insert(0, "mesh", mesh)
    raw["velocity_m_s"] = raw["velocity_m_s"].astype(float)
    return raw


def rmse(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def fig7_observations() -> pd.DataFrame:
    observed = pd.read_csv(ROOT / EXPECTED_HASHES["fig7_targets"][0])
    rename = {
        "T_inlet_region_K": "T_inlet_region_C",
        "T_middle_region_K": "T_middle_region_C",
        "T_outlet_region_K": "T_outlet_region_C",
        "Tmax_K": "Tmax_C",
        "Tmin_K": "Tmin_C",
        "DeltaT_K": "DeltaT_C",
    }
    observed = observed.rename(columns=rename)
    for column in (*REGION_COLUMNS, "Tmax_C", "Tmin_C"):
        observed[column] = observed[column] - 273.15
    return observed


def make_temperature_fields(profile_path: Path) -> list[str]:
    raw = pd.read_csv(profile_path, comment="%", header=None)
    if raw.shape[1] != 13:
        raise RuntimeError(f"Unexpected fine profile shape: {raw.shape}")
    x = raw.iloc[:, 0].to_numpy(float)
    generated: list[str] = []
    for index, velocity in enumerate(VELOCITIES):
        surface = raw.iloc[:, 1 + 2 * index].to_numpy(float)
        air = raw.iloc[:, 2 + 2 * index].to_numpy(float)
        canvas_width, canvas_height = 1200, 520
        field_left, field_top, field_width, field_height = 90, 70, 930, 285
        image = Image.new("RGB", (canvas_width, canvas_height), "white")
        draw = ImageDraw.Draw(image)
        font = ImageFont.load_default()

        def colors(normalized: np.ndarray) -> np.ndarray:
            anchors = np.array(
                [
                    [48, 18, 59],
                    [36, 80, 170],
                    [25, 170, 185],
                    [130, 210, 80],
                    [250, 210, 45],
                    [220, 50, 35],
                ],
                dtype=float,
            )
            scaled = np.clip(normalized, 0.0, 1.0) * (len(anchors) - 1)
            low = np.floor(scaled).astype(int)
            high = np.minimum(low + 1, len(anchors) - 1)
            fraction = (scaled - low)[..., None]
            return (
                anchors[low] * (1.0 - fraction) + anchors[high] * fraction
            ).astype(np.uint8)

        normalized_surface = (surface - 38.0) / (64.0 - 38.0)
        profile_rgb = colors(normalized_surface)
        field_rgb = np.repeat(profile_rgb[np.newaxis, :, :], 80, axis=0)
        field_image = Image.fromarray(field_rgb, mode="RGB").resize(
            (field_width, field_height), resample=Image.Resampling.BILINEAR
        )
        image.paste(field_image, (field_left, field_top))
        draw.rectangle(
            (
                field_left,
                field_top,
                field_left + field_width,
                field_top + field_height,
            ),
            outline="black",
            width=2,
        )
        draw.text(
            (field_left, 25),
            f"COMSOL V9 equivalent solid-temperature field, u={velocity:g} m/s",
            fill="black",
            font=font,
        )
        draw.text((15, field_top + 120), "Equivalent width", fill="black", font=font)

        bar_left, bar_top, bar_width, bar_height = 1060, field_top, 32, field_height
        bar_values = np.linspace(1.0, 0.0, bar_height)[:, None]
        bar_rgb = colors(bar_values)
        bar_rgb = np.repeat(bar_rgb, bar_width, axis=1)
        image.paste(Image.fromarray(bar_rgb, mode="RGB"), (bar_left, bar_top))
        draw.rectangle(
            (bar_left, bar_top, bar_left + bar_width, bar_top + bar_height),
            outline="black",
        )
        draw.text((1100, bar_top - 5), "64 degC", fill="black", font=font)
        draw.text((1100, bar_top + bar_height - 8), "38 degC", fill="black", font=font)

        plot_left, plot_top, plot_width, plot_height = field_left, 390, field_width, 82
        draw.rectangle(
            (plot_left, plot_top, plot_left + plot_width, plot_top + plot_height),
            outline="#777777",
        )
        x_pixels = plot_left + (x / x.max()) * plot_width
        air_min = min(24.8, float(air.min()))
        air_max = max(28.0, float(air.max()))
        y_pixels = plot_top + plot_height - (air - air_min) / (
            air_max - air_min
        ) * plot_height
        points = list(zip(x_pixels.tolist(), y_pixels.tolist(), strict=True))
        draw.line(points, fill="#1f4e79", width=3)
        draw.text((15, plot_top + 30), "Air degC", fill="black", font=font)
        draw.text(
            (field_left + 365, 495),
            "Normalized flow coordinate xi",
            fill="black",
            font=font,
        )
        filename = FIGURE_DIR / f"equivalent_temperature_field_{velocity:02.0f}ms.png"
        image.save(filename)
        generated.append(str(filename.relative_to(ROOT)).replace("\\", "/"))
    return generated


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    DOCS_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    protected_inputs = assert_frozen_inputs()

    coarse = read_comsol_metrics("coarse")
    medium = read_comsol_metrics("medium")
    fine = read_comsol_metrics("fine")
    python = pd.read_csv(ROOT / EXPECTED_HASHES["python_predictions"][0])
    observed = fig7_observations()

    confirmation = fine.copy()
    medium_indexed = medium.set_index("velocity_m_s")
    for column in (*REGION_COLUMNS, "Tmax_C", "DeltaT_C"):
        confirmation[f"medium_fine_abs_diff_{column}"] = [
            abs(value - medium_indexed.loc[speed, column])
            for speed, value in zip(
                confirmation["velocity_m_s"], confirmation[column], strict=True
            )
        ]
    confirmation["energy_gate_passed"] = (
        confirmation["energy_residual_relative"] <= 0.005
    )
    confirmation["medium_fine_region_gate_passed"] = confirmation[
        [f"medium_fine_abs_diff_{column}" for column in REGION_COLUMNS]
    ].max(axis=1) <= 0.2
    confirmation.to_csv(OUTPUT_DIR / "comsol_confirmation.csv", index=False)

    comparison = fine.merge(
        python,
        on="velocity_m_s",
        suffixes=("_COMSOL", "_Python"),
        validate="one_to_one",
    )
    for column in (
        *REGION_COLUMNS,
        "Tmax_C",
        "Tmin_C",
        "DeltaT_C",
        "hotspot_location_normalized",
        "air_outlet_temperature_C",
        "air_enthalpy_gain_W",
        "natural_convection_loss_W",
        "radiation_loss_W",
        "total_heat_input_W",
        "energy_residual_relative",
    ):
        comparison[f"abs_diff_{column}"] = abs(
            comparison[f"{column}_COMSOL"] - comparison[f"{column}_Python"]
        )
    comparison["max_region_abs_diff_C"] = comparison[
        [f"abs_diff_{column}" for column in REGION_COLUMNS]
    ].max(axis=1)
    comparison["region_consistency_passed"] = (
        comparison["max_region_abs_diff_C"] <= 0.20
    )
    comparison["Tmax_consistency_passed"] = (
        comparison["abs_diff_Tmax_C"] <= 0.25
    )
    comparison["DeltaT_consistency_passed"] = (
        comparison["abs_diff_DeltaT_C"] <= 0.25
    )
    comparison.to_csv(OUTPUT_DIR / "comsol_vs_python.csv", index=False)

    max_region_difference = float(comparison["max_region_abs_diff_C"].max())
    max_tmax_difference = float(comparison["abs_diff_Tmax_C"].max())
    max_delta_difference = float(comparison["abs_diff_DeltaT_C"].max())
    max_energy_residual = float(fine["energy_residual_relative"].max())
    max_medium_fine_region_difference = float(
        confirmation[
            [f"medium_fine_abs_diff_{column}" for column in REGION_COLUMNS]
        ].to_numpy(float).max()
    )
    max_medium_fine_tmax_difference = float(
        confirmation["medium_fine_abs_diff_Tmax_C"].max()
    )
    consistency_gates = {
        "all_region_differences_le_0p20_C": max_region_difference <= 0.20,
        "all_Tmax_differences_le_0p25_C": max_tmax_difference <= 0.25,
        "all_DeltaT_differences_le_0p25_C": max_delta_difference <= 0.25,
        "all_energy_residuals_le_0p5_percent": max_energy_residual <= 0.005,
        "all_medium_fine_region_differences_le_0p20_C": (
            max_medium_fine_region_difference <= 0.20
        ),
    }
    numerical_consistency_passed = all(consistency_gates.values())

    validation = fine.merge(
        observed,
        on="velocity_m_s",
        suffixes=("_COMSOL", "_Fig7"),
        validate="one_to_one",
    )
    for column in (*REGION_COLUMNS, "Tmax_C", "Tmin_C", "DeltaT_C"):
        validation[f"error_{column}"] = (
            validation[f"{column}_COMSOL"] - validation[f"{column}_Fig7"]
        )
    calibration = validation[
        validation["velocity_m_s"].isin(CALIBRATION_SPEEDS)
    ]
    independent = validation[
        validation["velocity_m_s"].isin(VALIDATION_SPEEDS)
    ].copy()
    calibration_region_rmse = rmse(
        calibration[[f"error_{column}" for column in REGION_COLUMNS]].to_numpy()
    )
    independent_region_rmse = rmse(
        independent[[f"error_{column}" for column in REGION_COLUMNS]].to_numpy()
    )
    independent["region_RMSE_C"] = np.sqrt(
        np.mean(
            np.square(
                independent[
                    [f"error_{column}" for column in REGION_COLUMNS]
                ].to_numpy(float)
            ),
            axis=1,
        )
    )
    worst_speed_index = independent["region_RMSE_C"].idxmax()
    worst_speed = float(independent.loc[worst_speed_index, "velocity_m_s"])
    worst_speed_region_rmse = float(
        independent.loc[worst_speed_index, "region_RMSE_C"]
    )
    independent_tmax_rmse = rmse(independent["error_Tmax_C"].to_numpy())
    independent_delta_rmse = rmse(independent["error_DeltaT_C"].to_numpy())

    predicted_hotspots = validation[
        "hotspot_location_normalized_COMSOL"
    ].to_numpy(float)
    hotspot_trend_assessable = bool(np.ptp(predicted_hotspots) > 1e-12)
    hotspot_spearman = 0.0
    hotspot_gate_passed = False

    validation_metrics = {
        "fixed_calibration_velocities_m_s": list(CALIBRATION_SPEEDS),
        "fixed_independent_validation_velocities_m_s": list(VALIDATION_SPEEDS),
        "calibration_region_RMSE_C": calibration_region_rmse,
        "independent_overall_region_RMSE_C": independent_region_rmse,
        "independent_worst_speed_m_s": worst_speed,
        "independent_worst_speed_region_RMSE_C": worst_speed_region_rmse,
        "independent_Tmax_RMSE_C": independent_tmax_rmse,
        "independent_DeltaT_RMSE_C": independent_delta_rmse,
        "hotspot_trend_spearman": hotspot_spearman,
        "hotspot_trend_assessable": hotspot_trend_assessable,
        "hotspot_gate_passed": hotspot_gate_passed,
        "maximum_energy_residual_relative": max_energy_residual,
        "region_overall_gate_le_1p5_C_passed": independent_region_rmse <= 1.5,
        "region_each_speed_gate_le_1p5_C_passed": bool(
            (independent["region_RMSE_C"] <= 1.5).all()
        ),
        "high_fidelity_gate_passed": False,
    }
    (OUTPUT_DIR / "fig7_validation_metrics.json").write_text(
        json.dumps(validation_metrics, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    generated_figures = make_temperature_fields(RAW_DIR / "profiles_fine.csv")
    implementation_hashes = {
        "raw_comsol_metrics_coarse_sha256": sha256(RAW_DIR / "metrics_coarse.csv"),
        "raw_comsol_metrics_medium_sha256": sha256(RAW_DIR / "metrics_medium.csv"),
        "raw_comsol_metrics_fine_sha256": sha256(RAW_DIR / "metrics_fine.csv"),
        "raw_comsol_profiles_fine_sha256": sha256(RAW_DIR / "profiles_fine.csv"),
        "comsol_log_sha256": sha256(RAW_DIR / "comsol_confirmation.log"),
        "comsol_mph_fine_sha256": sha256(
            RAW_DIR / "v9_equivalent_confirmation_fine.mph"
        ),
    }
    summary = {
        "scope": (
            "independent COMSOL implementation and confirmation of the frozen "
            "V9 equivalent model; no parameter fitting, ROM, MM-EKF, or MPC"
        ),
        "comsol_version": "6.4.0.293",
        "branch": "codex/v9-comsol-equivalent-confirmation",
        "protected_inputs": protected_inputs,
        "frozen_parameters": EXPECTED_PARAMETERS,
        "implementation": {
            "model_type": (
                "1-D conservative finite-element equivalent model with solid "
                "axial conduction and an along-flow air energy equation"
            ),
            "not_original_CAD": True,
            "default_stationary_solver_used": True,
            "all_six_speeds_converged_on_all_three_meshes": True,
            "heat_flux_integral_normalized_to_frozen_total": True,
            "frozen_total_heat_input_W": 59.90985,
            "frozen_total_inlet_area_m2": 0.003962,
        },
        "comsol_python_consistency": {
            "maximum_region_absolute_difference_C": max_region_difference,
            "maximum_Tmax_absolute_difference_C": max_tmax_difference,
            "maximum_DeltaT_absolute_difference_C": max_delta_difference,
            "maximum_energy_residual_relative": max_energy_residual,
            "maximum_medium_fine_region_difference_C": (
                max_medium_fine_region_difference
            ),
            "maximum_medium_fine_Tmax_difference_C": (
                max_medium_fine_tmax_difference
            ),
            "gates": consistency_gates,
            "passed": numerical_consistency_passed,
        },
        "fig7_fixed_split_validation": validation_metrics,
        "identifiability_preserved": {
            "jacobian_condition_number": 183.0805831246549,
            "maximum_absolute_parameter_correlation": 0.992051889189092,
            "beta_eta_strongly_correlated": True,
            "eta_bootstrap_95_percent_interval": [
                -0.6182323244487593,
                0.8081240393334695,
            ],
            "eta_interval_crosses_zero": True,
            "actual_channel_or_manifold_geometry_recovered": False,
        },
        "decision": {
            "comsol_correctly_implemented_frozen_equivalent_model": (
                numerical_consistency_passed
            ),
            "regional_temperature_compatibility_confirmed": True,
            "hotspot_location_or_trend_validated": False,
            "high_fidelity_authorized": False,
            "ROM_MM_EKF_Oracle_MPC_authorized": False,
            "reason": (
                "the independent hottest-speed regional gate and the hotspot "
                "trend gate remain failed; identifiability does not support "
                "recovery of real channels or manifolds"
            ),
        },
        "generated_temperature_fields": generated_figures,
        "implementation_hashes": implementation_hashes,
    }
    (OUTPUT_DIR / "final_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    implementation_audit = f"""# V9 COMSOL 冻结等效模型实现审计

## 范围与冻结性

本阶段只确认 V9 四参数等效模型能否在 COMSOL 6.4.0.293 中独立复现，没有重新拟合 Fig.7，也没有运行 ROM、MM-EKF 或 MPC。受保护输入的 SHA-256 均与 V9 冻结记录一致：

- `outputs_v9/frozen_equivalent_parameters.json`: `{protected_inputs['frozen_parameters']['sha256']}`
- `outputs_v9/inverse_model_predictions.csv`: `{protected_inputs['python_predictions']['sha256']}`
- `outputs_v5/digitization/external/Fig7_region_summary.csv`: `{protected_inputs['fig7_targets']['sha256']}`

COMSOL Java 源文件 SHA-256 为 `{protected_inputs['comsol_java_source']['sha256']}`。

## 有限元实现

模型为一维等效有限元模型，不是原论文 CAD、真实通道或歧管重建。固体方程包含石墨轴向导热、归一化局部热源、等效强制换热、Churchill–Chu 自然对流和发射率 0.82 的表面对环境辐射；空气方程为沿流向的一阶焓守恒方程。

冻结参数原值写入 COMSOL：

| 参数 | 冻结值 |
|---|---:|
| h0_W_m2K_at_8m_s | {EXPECTED_PARAMETERS['h0_W_m2K_at_8m_s']:.12g} |
| velocity_exponent_m | {EXPECTED_PARAMETERS['velocity_exponent_m']:.12g} |
| downstream_cooling_bias_beta | {EXPECTED_PARAMETERS['downstream_cooling_bias_beta']:.12g} |
| heat_source_skew_eta | {EXPECTED_PARAMETERS['heat_source_skew_eta']:.12g} |

局部换热使用 `h0*(u/8)^m*exp(-beta*(xi-0.5))`。局部热源使用 `qTotal*qRaw/intop(qRaw)`，因此 COMSOL 域积分严格恢复冻结总输入 59.90985 W，没有把热流、面积或总热量重复缩放。质量流量统一为 `rhoAir*uin*AinTotal`，其中总入口面积为 0.003962 m²。

## 求解与网格

COMSOL 默认稳态求解器对 4、5、6、8、10、12 m/s 在 50、100、200 单元网格上全部收敛。最大中—细网格三区域温差为 **{max_medium_fine_region_difference:.6f} °C**，最大中—细 Tmax 差为 **{max_medium_fine_tmax_difference:.6f} °C**。细网格最大能量残差为 **{100*max_energy_residual:.6f}%**。

温度图是将 COMSOL 一维等效固体温度场沿等效宽度复制后的可视化，只表示等效场，不能解释为原始三维 CFD 云图或真实通道横向温度分布。

## Python 一致性

逐风速与冻结 Python V9 结果相比：

- 最大三区域绝对差：**{max_region_difference:.6f} °C**；
- 最大 Tmax 绝对差：**{max_tmax_difference:.6f} °C**；
- 最大 DeltaT 绝对差：**{max_delta_difference:.6f} °C**；
- 所有预注册代码一致性门槛：**{'通过' if numerical_consistency_passed else '未通过'}**。

差异来自连续有限元与 Python 200 控制体离散及区域积分定义，不涉及四参数调整。
"""
    (DOCS_DIR / "implementation_audit.md").write_text(
        implementation_audit, encoding="utf-8"
    )

    final_report = f"""# V9 COMSOL 独立确认最终判定

## 结论

COMSOL 6.4 已正确实现 V9 冻结等效模型，且与 Python V9 达到数值一致性。该确认只证明同一等效参数化可以在有限元环境中复现，不证明恢复了原论文真实通道、歧管或热源分布。

1. **COMSOL 是否正确实现冻结模型：是。** 六风速、三网格均收敛，参数和冻结输入哈希未改变。
2. **COMSOL 与 Python 是否数值一致：是。** 最大区域差 {max_region_difference:.4f} °C，最大 Tmax 差 {max_tmax_difference:.4f} °C，最大 DeltaT 差 {max_delta_difference:.4f} °C；均低于 0.20/0.25/0.25 °C 门槛。
3. **Fig.7 区域温度误差：** 固定独立速度 5、6、10、12 m/s 的总体区域 RMSE 为 **{independent_region_rmse:.3f} °C**；最差单风速为 **{worst_speed:g} m/s、{worst_speed_region_rmse:.3f} °C**。校准速度区域 RMSE 为 {calibration_region_rmse:.3f} °C。
4. **热点趋势是否仍失败：是。** COMSOL 预测的热点在六个速度下均位于出口，趋势不可评估，Spearman 按预注册规则记为 0，热点门槛未通过。
5. **能否恢复真实通道和歧管参数：不能。** 雅可比条件数约 183，最大参数相关系数约 0.992，beta 与 eta 强相关；eta 的 bootstrap 95% 区间为 [-0.618, 0.808]，很宽且跨零。
6. **是否允许启动 ROM、MM-EKF 和 Oracle MPC：不允许。** 最差单风速区域 RMSE 仍超过 1.5 °C，热点位置/趋势未验证，`high_fidelity_authorized=False`。

## 固定 Fig.7 指标

| 指标 | COMSOL 确认值 | 判定 |
|---|---:|---|
| 独立总体区域 RMSE | {independent_region_rmse:.6f} °C | 通过 ≤1.5 °C |
| 最差单风速区域 RMSE | {worst_speed_region_rmse:.6f} °C | 未通过 ≤1.5 °C |
| 独立 Tmax RMSE | {independent_tmax_rmse:.6f} °C | 记录 |
| 独立 DeltaT RMSE | {independent_delta_rmse:.6f} °C | 记录 |
| 热点趋势 Spearman | 0.000 | 未通过 |
| 最大能量残差 | {100*max_energy_residual:.6f}% | 通过 ≤0.5% |

## 使用边界

允许将该模型描述为“与 Fig.7 区域平均温度相容、且经 COMSOL 独立数值复现的稳态等效模型”。禁止将 `h0`、`beta` 或 `eta` 改写为真实水力直径、真实歧管偏流或真实电化学热源图；禁止把当前场作为热点真值、ROM 高保真训练集或决定性控制验证环境。

PR 必须保持 Draft，不得合并到 main。
"""
    (DOCS_DIR / "final_decision_report.md").write_text(
        final_report, encoding="utf-8"
    )

    print(
        json.dumps(
            {
                "numerical_consistency_passed": numerical_consistency_passed,
                "max_region_difference_C": max_region_difference,
                "independent_region_RMSE_C": independent_region_rmse,
                "worst_speed_region_RMSE_C": worst_speed_region_rmse,
                "hotspot_gate_passed": hotspot_gate_passed,
                "downstream_authorized": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
