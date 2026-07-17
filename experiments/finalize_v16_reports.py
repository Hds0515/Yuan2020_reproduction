"""Build the final V16 integrity audit and human-readable decision reports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs_v16"
DOCS = ROOT / "docs_v16"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def write_text(path: Path, text: str) -> None:
    path.write_text(text.strip() + "\n", encoding="utf-8")


def main() -> None:
    DOCS.mkdir(exist_ok=True)
    v15 = read_json(ROOT / "outputs_v15/final_summary.json")
    residual = read_json(OUT / "v15_dynamic_residual_diagnosis.json")
    single = read_json(OUT / "five_zone_rom_validation_metrics.json")
    projected = read_json(OUT / "projected_two_mass_rom_metrics.json")
    observer = read_json(OUT / "observer_method_study_summary.json")
    single_freeze = read_json(OUT / "frozen_five_zone_rom.json")
    projected_freeze = read_json(OUT / "frozen_projected_two_mass_rom.json")
    integrity = {
        "V15_frozen_reference_matches_declared_sha256": sha256(
            ROOT / "outputs_v15/frozen_high_resolution_reference.json"
        )
        == v15["frozen_reference_sha256"],
        "single_mass_current_protocol_matches_historical_freeze": sha256(
            OUT / "preregistered_v16_protocol.json"
        )
        == single_freeze["protocol_sha256"],
        "single_mass_protocol_note": "The protocol timestamp was inadvertently regenerated after the negative freeze; the mismatch is disclosed and this failed baseline is not used for a positive claim.",
        "two_mass_protocol_matches_freeze": sha256(
            OUT / "preregistered_two_mass_rom_protocol.json"
        )
        == projected_freeze["protocol_sha256"],
        "two_mass_model_source_matches_freeze": sha256(
            ROOT / "models_v16/projected_two_mass_five_zone_rom.py"
        )
        == projected_freeze["source_sha256"],
    }
    write_json(OUT / "integrity_audit.json", integrity)

    final = {
        "branch_stage": "V16 seven-step physics-ROM-observer programme",
        "protected_parent": {
            "tag": "v15_fixed60_reference_confirmed",
            "frozen_reference_sha256": v15["frozen_reference_sha256"],
            "V15_outputs_overwritten": False,
        },
        "step_1_scope_freeze": {
            "complete": True,
            "valid_scope": "23 C ambient, fixed 60% PWM, 2-40 A, measured current and stack voltage inputs",
        },
        "step_2_preregistered_data_registry": {
            "complete": True,
            "fresh_independent_experiment_available": False,
            "experimental_registry_template": "outputs_v16/fresh_experiment_registry_TEMPLATE.csv",
            "synthetic_method_datasets_are_labelled_as_simulation": True,
        },
        "step_3_dynamic_discrepancy_diagnosis": {
            "complete": True,
            "worst_error": residual["worst_error"],
            "regional_order_failure_times_s": residual[
                "regional_order_failure_times_s"
            ],
            "conclusion": "Load- and position-dependent structural discrepancy dominates implementation error; available anchors cannot uniquely separate fan lag, sensor lag, water state and local heat release.",
        },
        "step_4_high_resolution_confirmation": {
            "COMSOL_implementation_consistency_passed": v15["COMSOL"][
                "implementation_consistency_passed"
            ],
            "strict_experimental_gate_passed": False,
            "maximum_dynamic_error_C": v15["experimental_evidence"][
                "fixed60_regional_dynamic"
            ]["maximum_absolute_error_C"],
            "regional_order_accuracy_percent": v15["experimental_evidence"][
                "fixed60_regional_dynamic"
            ]["region_order_accuracy_percent"],
            "reason": "4.271 C exceeds the preregistered 4.0 C limit and 94.737% is below the 95% ordering limit; no unused confirmation experiment exists.",
        },
        "step_5_five_zone_ROM": {
            "single_mass_ablation_completed": True,
            "single_mass_synthetic_gate_passed": single[
                "synthetic_method_gate_passed"
            ],
            "two_mass_projection_completed": True,
            "two_mass_validation_mean_temperature_RMSE_C": projected[
                "validation"
            ]["mean_temperature_RMSE_C"],
            "two_mass_validation_zone_RMSE_C": projected["validation"][
                "zone_overall_RMSE_C"
            ],
            "two_mass_validation_Tmax_RMSE_C": projected["validation"][
                "Tmax_RMSE_C"
            ],
            "two_mass_validation_DeltaT_RMSE_C": projected["validation"][
                "DeltaT_RMSE_C"
            ],
            "two_mass_synthetic_final_gate_passed": projected[
                "synthetic_method_gate_passed"
            ],
            "independent_experimental_validation_passed": False,
            "formal_claim_authorized": False,
            "conclusion": "The two-mass structure preserves mean temperature, stability, hotspot zone and energy, but five-zone extrema reconstruction is not accurate enough.",
        },
        "step_6_sparse_observer": {
            "synthetic_method_study_completed": True,
            "synthetic_numeric_gate_passed": observer[
                "synthetic_numeric_gate_passed"
            ],
            "median_Tmax_RMSE_improvement_percent": observer[
                "median_Tmax_RMSE_improvement_percent"
            ],
            "preregistered_experimental_confirmation_passed": False,
            "core_innovation_authorized": False,
            "reason": "The observer used the failed single-mass method model and synthetic parent truth; it is an algorithm check only.",
        },
        "step_7_control_decision": {
            "MPC_run": False,
            "control_research_authorized": False,
            "decision": "Do not run decisive MPC. The variable-PWM plant and independently validated ROM are unavailable.",
        },
        "research_decision": {
            "what_is_available": "An experimentally constrained fixed-60% high-resolution reference and a reproducible ROM structure study.",
            "what_is_not_available": "A general high-fidelity plant or an independently validated five-zone hotspot ROM.",
            "next_required_action": "Collect a newly reserved multi-point fixed-60% dynamic confirmation experiment; measure fan speed and preferably airflow or pressure. Variable-PWM work requires separate PWM40/60/80 telemetry.",
            "thesis_writing": "Methods, reproduction, discrepancy and negative-result chapters may start; final innovation claims must wait for fresh confirmation data.",
        },
        "integrity": integrity,
        "reference_temperature_C": 55.0,
        "reference_temperature_is_not_a_safety_limit": True,
        "draft_PR_only": True,
        "merge_to_main_allowed": False,
    }
    write_json(OUT / "final_summary.json", final)

    scope = f"""
# V16研究范围与数据注册

V15标签 `v15_fixed60_reference_confirmed` 保持不变，冻结文件SHA-256为
`{v15['frozen_reference_sha256']}`。V16没有覆盖任何V15文件。

当前可辩护范围仅为23°C环境、固定60% PWM、2–40 A，并使用实测电流和
电堆电压作为外生输入。Guo 2024的2 A阶梯已用于辨识，4 A阶梯已经执行过
一次最终确认；仓库中不存在尚未使用的新动态实验，因此不能再次把4 A数据
称为独立确认。

`high_resolution_dataset_registry.json`和`two_mass_dataset_registry.json`中的
数据均明确标记为冻结V15父模型模拟。它们只用于ROM和观测器软件方法验证。
下一次真实实验使用`fresh_experiment_registry_TEMPLATE.csv`，必须在采集前冻结
calibration、validation和final-confirmation用途，并在采集后填写SHA-256。
"""
    write_text(DOCS / "scope_and_data_registry.md", scope)

    worst = residual["worst_error"]
    discrepancy = f"""
# V15动态失配诊断

最大误差出现在{worst['time_s']:.0f} s、40 A升载末端的出口区域：实验
{worst['measured_C']:.3f}°C，模型{worst['predicted_C']:.3f}°C，误差
{worst['error_C']:.3f}°C。唯一的区域排序失败点为
{', '.join(f'{v:.0f} s' for v in residual['regional_order_failure_times_s'])}。

残差随负载、流向位置和升降载阶段改变符号。该结构性特征远大于COMSOL—
Python最大0.00739°C实现差，不能归因于数值求解器。现有100 s间隔数字化锚点
也不足以唯一分离传感器热惯性、风扇/流量惯性、膜水状态和局部产热分配。
因此本阶段没有利用旧确认集重新调参。
"""
    write_text(DOCS / "dynamic_discrepancy_diagnosis.md", discrepancy)

    high_resolution = f"""
# 高分辨率参考确认

COMSOL 6.4与Python冻结方程的一致性已经通过，但这只证明实现一致。当前独立
空间RMSE为{v15['experimental_evidence']['independent_spatial_RMSE_C']:.3f}°C，
固定60% PWM区域动态RMSE为
{v15['experimental_evidence']['fixed60_regional_dynamic']['regional_RMSE_C']:.3f}°C。

严格门槛仍失败：最大误差4.271°C高于4.0°C，区域排序94.737%低于95%。
所以模型名称保持为“实验约束的固定占空比高分辨率热参考”，不能改写成变PWM
高保真数字孪生，也不能用模拟场景替代新的实验确认。
"""
    write_text(DOCS / "high_resolution_confirmation.md", high_resolution)

    sv = single["validation"]
    pv = projected["validation"]
    pf = projected["synthetic_final_confirmation"]
    rom_doc = f"""
# 五区域ROM结构验证

## 单热容五状态消融

冻结单热容RC链在验证集的区域RMSE为{sv['zone_overall_RMSE_C']:.3f}°C，
Tmax RMSE为{sv['Tmax_RMSE_C']:.3f}°C，DeltaT RMSE为
{sv['DeltaT_RMSE_C']:.3f}°C。辨识容量集中到一个区域，部分导热参数触及上界，
说明把表面—核心双时间尺度压缩成五个温度状态并不可辨识。

## 五区双热容投影

第二种结构采用五个流向区域，每区一个表面状态和一个核心状态，共10状态；
动态矩阵由6000状态父模型容量加权Galerkin投影得到，不辨识动态参数。验证集
平均温度RMSE为{pv['mean_temperature_RMSE_C']:.3f}°C、热点区准确率
{pv['hotspot_zone_accuracy_percent']:.1f}%，但区域RMSE为
{pv['zone_overall_RMSE_C']:.3f}°C、Tmax RMSE为{pv['Tmax_RMSE_C']:.3f}°C、
DeltaT RMSE为{pv['DeltaT_RMSE_C']:.3f}°C。新seed 2601模拟确认的区域RMSE为
{pf['zone_overall_RMSE_C']:.3f}°C，仍未通过。

结论：双热容结构适合平均温度快速预测，但当前五区划分不能可信重构局部极值
和温差。两种模型均严格守恒、稳定且满足在线计算速度；这不等于精度门槛通过。
不得继续通过增加自由参数追逐模拟确认集。
"""
    write_text(DOCS / "five_zone_ROM_validation.md", rom_doc)

    observer_doc = f"""
# 稀疏观测器方法验证

固定首末两个区域传感器、10个噪声seed和V6继承门控参数后，模拟方法测试中
门控MM-EKF相对单EKF的Tmax RMSE中位改善为
{observer['median_Tmax_RMSE_improvement_percent']:.2f}%，95%区间为
[{observer['paired_mean_improvement_95pct_CI_percent'][0]:.2f}%,
{observer['paired_mean_improvement_95pct_CI_percent'][1]:.2f}%]，标称和纯噪声
中位损失均为0%。数值门槛通过。

但该测试依赖冻结父模型模拟，并使用未通过ROM门槛的单热容方法模型。因此它
只证明门控代码在规定模拟压力测试中工作，不能确定为论文核心创新。需要新的
预注册实验确认集，并在通过的ROM上重新确认。
"""
    write_text(DOCS / "observer_method_study.md", observer_doc)

    control_doc = """
# 控制研究停止决定

本阶段没有运行MPC。原因不是计算资源不足，而是三个前置条件不成立：V15严格
实验门失败、变PWM高保真输入模型不存在、五区域ROM未通过局部极值验证。
在这些条件下调MPC只能得到依赖错误模型的结论。

控制研究保持降级。只有新实验确认高分辨率模型、独立确认ROM，并取得PWM40、
PWM60、PWM80下风扇转速/流量或压差遥测后，才允许重新预注册公平控制实验。
55°C始终只是参考温度，不是安全上限。
"""
    write_text(DOCS / "control_decision.md", control_doc)

    final_report = f"""
# V16最终决策报告

七步程序已经执行到证据允许的边界：V15被保护，数据用途和SHA-256已注册，
动态残差已定位，两种五区ROM结构及稀疏观测器方法均已冻结运行，MPC按停止规则
未运行。

1. 高分辨率参考：COMSOL—Python一致，但严格实验高保真门失败。
2. 新独立实验：不存在；模拟数据没有被冒充为实验。
3. 单热容五状态ROM：失败。
4. 五区双热容投影ROM：平均温度有效，但区域极值和DeltaT失败。
5. 门控MM-EKF：模拟数值门通过，但没有获得实验核心创新授权。
6. MPC：未授权、未运行。
7. 论文：可以开始方法、复现、结构失配和负面结果章节；最终创新结论必须等待
   新的多点动态确认实验。

当前能交付的是“固定60% PWM实验约束高分辨率参考”，不是通用高保真模型。
下一项不可由软件替代的工作是采集新的多点温度、风扇转速以及流量或压差数据。
Draft PR不得合并到main。
"""
    write_text(DOCS / "final_decision_report.md", final_report)
    write_text(OUT / "final_decision_report.md", final_report)

    manifest: dict[str, str] = {}
    for root in (OUT, DOCS):
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            if path.name == "sha256_manifest.json":
                continue
            manifest[str(path.relative_to(ROOT)).replace("\\", "/")] = sha256(path)
    write_json(OUT / "sha256_manifest.json", manifest)


if __name__ == "__main__":
    main()
