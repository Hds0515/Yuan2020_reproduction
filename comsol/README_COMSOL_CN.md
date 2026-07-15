# COMSOL 6.4 最小复杂度收敛任务

当前模型选择方案 B：代表性并联通道。完整决定和尺度证明见 `docs/equivalent_geometry_decision.md`。总入口面积为 `Ain=3962 mm²`，代表性通道面积为 `W*Hair`，所有总量使用 `areaScale=Ain/(W*Hair)` 缩放；Stage 1 的目标总质量流量为 `rho*uin*Ain`。

运行：

```powershell
.\run_comsol_windows.bat
```

脚本只执行 Stage 1：关闭传热、入口 0.1 m/s、Laminar Flow、映射/扫掠网格（厚度方向 4 层）、Stokes 初始化、全耦合 PARDISO，再延续到惯性层流。编译物、结果和日志分别写入 `generated_v4/`、`results_v4/` 和 `logs_v4/`。

当前 Stage 1 未收敛。压力出口、显式压力点试验、截面细化和扫掠网格均已尝试；PARDISO 的压力—速度线性系统仍出现病态/相对残差失败并达到 Newton 迭代上限。按预定停止规则，Stage 2 和 Stage 3 未启动。

可靠状态以 `status_v4.json` 为准。没有收敛的 MPH、质量守恒验收和能量账本时，不得声称速度场、压力场、温度场、网格收敛或 CFD 能量守恒已经完成。
