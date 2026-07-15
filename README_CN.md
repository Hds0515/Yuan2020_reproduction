# Yuan2020 热管理复现：v4 方法学修正版

本分支修复了 v3 的三个关键问题：多模型观察器场景未真正进入 truth、MPC 在整个预测域只用单一控制量，以及 55 °C 被混称为安全上限。受保护的原始 JPG、V2 raw/cleaned CSV、`frozen_parameters_v2.json` 和 `frozen_parameters_v3.json` 不会被改写。

## 复现

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_all.py --clean-output
```

如需重试 COMSOL Stage 1：

```powershell
.\.venv\Scripts\python.exe run_all.py --clean-output --with-comsol
```

新结果只写入 `outputs_v4/`、`comsol/generated_v4/`、`comsol/results_v4/` 和 `comsol/logs_v4/`。旧结果不会复制进 fresh bundle；V2 状态由 `baseline_v2` 标签和 `baseline_v2_outputs/README.md` 只读引用。

## 当前结论

- 55 °C 仅为论文最优工作温度/控制参考，不是安全上限。`T_safe` 默认为空；未提供外部可信值时不作安全结论。
- 参数模型库属于 MM-EKF/MMAE，不是未知气流方向的切换观察器。它在 5/5 个真实失配场景降低热点 RMSE，平均改善 30.48%，配对 bootstrap 95% 区间为 9.27%–51.69%。
- Hotspot-MPC 使用 12 步、三控制块（1–4、5–8、9–12），包含方向反转代价、10 s 最小驻留时间和 2 s 反转无效风量时间。
- 在可行场景中，MPC 相对 AuthorMeasured-PI-SMC 的 60 s 后热点峰值改善 1.44%，温差 RMSE 改善 46.78%；代价为风机能耗代理增加 30.30%、方向切换增加 316 次、平均温度跟踪 RMSE 增加 121.37%。`actuator_limited` 不计入正常均值。
- Fig.16(b)/(c) 已直接从局部放大图重新提取；参数与模型选择不使用 700–1200 s 验证区反向调节。
- COMSOL 采用代表性并联通道方案 B，`areaScale=Ain/(W*Hair)` 实际进入总质量流量。0.1 m/s Stage 1 仍未收敛，因此 0.5–4 m/s 延续、传热、网格无关性和 CFD 能量守恒均未执行，也未生成虚假结果。

机器可读汇总见 `outputs_v4/final_summary.json`，人工验收结论见 `outputs_v4/stage_acceptance_report.md`，COMSOL 证据见 `comsol/status_v4.json` 和 `comsol/logs_v4/`。
