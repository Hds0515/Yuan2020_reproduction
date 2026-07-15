# Yuan 2020 PEM 燃料电池热管理复现 v3

这是对 V2 工程的审计式重建。受保护的论文 PDF、三张原始 JPG、原始 raw CSV 和 `identification/frozen_parameters_v2.json` 未被修改；V2 的逐文件 SHA-256 副本保存在本仓库基线提交/tag `baseline_v2` 中。

当前结论：图片数字化、三节点辨识/独立验证、五节点模型、EKF、多模型观测器和热点 MPC 均已实际运行；COMSOL 6.4 Java 已实际编译，但稳态求解未收敛，因此网格无关性、CFD 能量守恒和 COMSOL—三节点交叉验证均明确标记为未完成。

## 一键复现

建议先在工程根目录建立环境并安装锁定依赖：

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\reproduce_python_only.bat
```

Python-only 入口会从 `source_images/Fig14.jpg`、`Fig15.jpg`、`Fig16.jpg` 开始，依次重新数字化、辨识、验证、运行五节点/观测器/控制实验和测试，并把本次完整产物收集到全新的 `outputs_v3`。它不把 V2 cleaned CSV、PNG 或结果 JSON 当作模型输入；V2 CSV 仅在图片提取完成后用于审计性比较。

若需要同时重试 COMSOL：

```powershell
.\reproduce_all.bat
```

COMSOL 求解可能耗时较长；当前失败证据和安全中止点见 `comsol/status.json` 与 `comsol/logs/solve_twelfth.log`。

## 核心结果

- Fig14 校准 RMSE：0.2491 °C。
- Fig16 0–700 s 校准 RMSE：0.2506 °C。
- Fig16 700–1200 s 独立验证 RMSE：0.2960 °C。
- 三节点最大相对能量残差：3.57e-14。
- 三/五节点平均温度 RMSE：0.00449 °C；五节点最大绝对能量残差：1.18e-11 W。
- 两传感器推荐位置：节点 2、5；跨六类失配场景的热点 RMSE 为 0.0681 °C。
- 多模型 EKF 平均热点 RMSE 0.0467 °C，单模型 EKF 为 0.0512 °C。
- 七场景中热点 MPC 使用两传感器 EKF 状态；与论文式中点 PI–SMC 相比降低热点越界 RMSE 和平均温差，但风机能耗代理更高。不可行场景单独标记，不混入常规可行性结论。

## 导航

- 总报告：`docs/final_reproduction_report_CN.md`
- 工程审计：`audit/project_audit.md`
- 数字化：`docs/digitization_audit.md`
- 模型假设：`docs/model_assumptions.md`
- 可辨识性：`docs/parameter_identifiability.md`
- COMSOL：`docs/comsol_cross_validation.md`
- 五节点与观测器：`docs/five_node_observer_results.md`
- MPC：`docs/hotspot_mpc_results.md`
- 机器可读状态：`outputs/final_summary.json`

本工程复现的是独立空气冷却流道方向热管理，不能称为开放阴极反应空气—冷却空气耦合模型，也不支持氧传输、水淹或膜含水量结论。
