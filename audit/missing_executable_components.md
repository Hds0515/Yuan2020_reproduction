# V2 缺失或不可执行组件

## 必须重建

- 原图颜色分割、像素轨迹提取、坐标转换、清洗、覆盖图和不确定度报告。
- 仅使用规定校准集的参数辨识入口，包括热源模型 A/B、加权残差、多起点、边界、局部协方差、bootstrap/profile 和可辨识性分析。
- 显式 `dt`、显式控制配置、抗饱和以及逐项能量账本的三节点模型。
- 可清空输出目录并从原图重新生成结果的一键入口和测试。
- 可编译、可求解、包含 `maxop1`、完整能量积分、三种网格、稳态速度扫描和瞬态换向研究的 COMSOL 6.4 模型。
- 真正调用三节点冻结模型的 COMSOL 交叉验证。
- 能量一致的五节点动力学模型。
- EKF、多模型观测器、传感器布局搜索和不确定性实验。
- 使用观测状态的热点约束 MPC，以及与 PI+SMC 的公平多场景比较。

## V2 不能重建的已有文件

- `outputs/01_fig15_digitized_current.png` 至 `outputs/05_parameter_correlation.png`
- `outputs/three_node_fig14_calibration.csv`
- `outputs/three_node_fig16_validation.csv`
- `identification/identified_parameters.csv` 及相关统计文件
- `outputs/five_node_EKF_timeseries.csv`
- `outputs/five_node_MPC_timeseries.csv`
- `outputs/five_node_observer_mpc_metrics.json`

这些文件可作为基线对照，但不得作为 v3 生成流程的输入。
