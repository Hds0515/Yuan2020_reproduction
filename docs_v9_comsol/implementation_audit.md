# V9 COMSOL 冻结等效模型实现审计

## 范围与冻结性

本阶段只确认 V9 四参数等效模型能否在 COMSOL 6.4.0.293 中独立复现，没有重新拟合 Fig.7，也没有运行 ROM、MM-EKF 或 MPC。受保护输入的 SHA-256 均与 V9 冻结记录一致：

- `outputs_v9/frozen_equivalent_parameters.json`: `e4504981c33e2381cb29c2888fd0604a7b745a4c2bbf8523073b10f7fd2c7343`
- `outputs_v9/inverse_model_predictions.csv`: `3d3ef55aae16dd7b3e067a533528e5aa94dc4a022d6197d7fe4ee5ee71042028`
- `outputs_v5/digitization/external/Fig7_region_summary.csv`: `01b87ddea0591e84c30cfa8bbe2c789bed47a36022895f18b2e4c7cee90ffe50`

COMSOL Java 源文件 SHA-256 为 `185e8a6610dc4d342a0abd878efb9b742aca78c27b4188f0486c07614ed7caa7`。

## 有限元实现

模型为一维等效有限元模型，不是原论文 CAD、真实通道或歧管重建。固体方程包含石墨轴向导热、归一化局部热源、等效强制换热、Churchill–Chu 自然对流和发射率 0.82 的表面对环境辐射；空气方程为沿流向的一阶焓守恒方程。

冻结参数原值写入 COMSOL：

| 参数 | 冻结值 |
|---|---:|
| h0_W_m2K_at_8m_s | 88.7532993434 |
| velocity_exponent_m | 0.657294908057 |
| downstream_cooling_bias_beta | 1.36793473474 |
| heat_source_skew_eta | -0.311866508779 |

局部换热使用 `h0*(u/8)^m*exp(-beta*(xi-0.5))`。局部热源使用 `qTotal*qRaw/intop(qRaw)`，因此 COMSOL 域积分严格恢复冻结总输入 59.90985 W，没有把热流、面积或总热量重复缩放。质量流量统一为 `rhoAir*uin*AinTotal`，其中总入口面积为 0.003962 m²。

## 求解与网格

COMSOL 默认稳态求解器对 4、5、6、8、10、12 m/s 在 50、100、200 单元网格上全部收敛。最大中—细网格三区域温差为 **0.005177 °C**，最大中—细 Tmax 差为 **0.000004 °C**。细网格最大能量残差为 **0.028160%**。

温度图是将 COMSOL 一维等效固体温度场沿等效宽度复制后的可视化，只表示等效场，不能解释为原始三维 CFD 云图或真实通道横向温度分布。

## Python 一致性

逐风速与冻结 Python V9 结果相比：

- 最大三区域绝对差：**0.039817 °C**；
- 最大 Tmax 绝对差：**0.018636 °C**；
- 最大 DeltaT 绝对差：**0.014915 °C**；
- 所有预注册代码一致性门槛：**通过**。

差异来自连续有限元与 Python 200 控制体离散及区域积分定义，不涉及四参数调整。
