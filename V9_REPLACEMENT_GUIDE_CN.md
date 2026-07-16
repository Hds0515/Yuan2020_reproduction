# V9 文件覆盖与运行说明

## 覆盖位置

将补丁 ZIP 解压到 `Yuan2020_reproduction_v3` 工程根目录，并保留目录结构。新增文件不会覆盖 V7/V8 结果。

新增内容：

- `comsol/v9_inverse/`：四参数等效稳态模型；
- `experiments/run_v9_limited_inverse_identification.py`：多起点辨识、独立验证、Bootstrap 和可辨识性审计；
- `run_v9_inverse.py`：运行入口；
- `outputs_v9/`：冻结参数、预测、图和机器汇总；
- `docs_v9/`：中文方法与结果报告；
- `tests/test_v9_inverse_identification.py`：V9测试；
- `audit_v9/`：测试记录和SHA-256清单；
- `CODEX_V9_COMSOL_CONFIRMATION_PROMPT_CN.md`：交给本地Codex继续做COMSOL一致性确认。

## 运行

在工程根目录执行：

```powershell
python run_v9_inverse.py
python -m pytest -q
```

当前已验证结果：

```text
24 passed
```

## 科学结论

V9将独立区域总体RMSE由V8 M2的2.846°C降低到1.157°C，温差RMSE降低到1.648°C；但最差单风速RMSE为1.622°C，热点始终位于出口附近，热点趋势没有通过。

因此：

- 可以用于区域平均温度的等效稳态仿真和敏感性分析；
- 不能声称恢复了原论文真实通道、歧管、热源或接触结构；
- 不能作为高保真真值训练ROM、验证MM-EKF或得出Oracle MPC决定性结论。
