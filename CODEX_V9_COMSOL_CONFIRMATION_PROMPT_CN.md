继续处理 `Yuan2020_reproduction_v3` 工程。

当前 V9 已在 Python 等效模型中完成受限四参数逆向辨识。禁止重新调整参数以迎合验证集；本阶段只做 COMSOL 独立实现与一致性确认。

## 一、建立分支与保护

从当前 V9 工程建立：

`codex/v9-comsol-equivalent-confirmation`

不得覆盖：

- `outputs_v7/`、`outputs_v8/`、`outputs_v9/`
- `docs_v7/`、`docs_v8/`、`docs_v9/`
- 历史 MPH、原始 Fig.7 数字化数据

新结果写入：

- `comsol/v9_confirmation/`
- `outputs_v9_comsol/`
- `docs_v9_comsol/`

PR 保持 Draft，不合并。

## 二、冻结参数，不得重新拟合

读取：

`outputs_v9/frozen_equivalent_parameters.json`

冻结值为：

- `h0_W_m2K_at_8m_s = 88.753299`
- `velocity_exponent_m = 0.65729491`
- `downstream_cooling_bias_beta = 1.3679347`
- `heat_source_skew_eta = -0.31186651`

参数含义：

1. `h0`：8 m/s 下等效内部换热强度，综合未解析的水力直径和有效接触面积；
2. `m`：速度指数；
3. `beta`：沿流向等效冷却偏置；
4. `eta`：沿流向等效热源偏置。

这些参数不是原论文真实 CAD 参数，不得改写为真实通道尺寸、真实歧管或真实热流分布。

## 三、在 COMSOL 中实现同一等效模型

使用 V8 冻结的：

- 总入口面积；
- 总热输入；
- Table 2 材料参数；
- 环境温度；
- Churchill–Chu 自然对流；
- 石墨发射率 0.82。

定义归一化流向坐标：

`xi = x/L`

定义等效内部换热系数：

`h_forced(x,u) = h0*(u/8[m/s])^m*exp(-beta*(xi-0.5))`

定义局部热流：

`q_raw(x) = 1 + eta*(2*xi-1)`

必须通过积分归一化，使总热输入严格等于 V8 冻结总热输入：

`q_flux(x) = q_total*q_raw(x)/integral(q_raw dA)`

禁止因单位、面积或方向定义造成重复缩放。

该模型可以采用“等效边界换热 + 沿程空气能量方程”，不要求伪装成显式原始通道 CFD。

## 四、六风速独立复算

运行：

4、5、6、8、10、12 m/s。

输出：

- 入口、中部、出口区域平均温度；
- Tmax、Tmin、DeltaT；
- 热点位置；
- 空气出口温度；
- 空气焓升；
- 自然对流损失；
- 辐射损失；
- 总热输入；
- 能量残差；
- 温度云图。

## 五、先验证代码一致性，不重新拟合 Fig.7

将 COMSOL 结果与：

`outputs_v9/inverse_model_predictions.csv`

进行逐风速比较。

一致性门槛：

- 三区域温度差均 ≤ 0.20 °C；
- Tmax 差 ≤ 0.25 °C；
- DeltaT 差 ≤ 0.25 °C；
- 能量误差 ≤ 0.5%；
- 中细网格区域温度差 ≤ 0.2 °C。

若不一致，只允许修复：

- 单位；
- 坐标方向；
- 热流归一化；
- 面积尺度；
- 边界选择；
- 后处理区域定义。

不得调整四个冻结参数。

## 六、再按冻结验证集评价 Fig.7

仍使用：

- 校准速度：4、8 m/s；
- 独立速度：5、6、10、12 m/s。

不得改变数据划分。

预期 Python V9 指标为：

- 独立总体区域 RMSE 约 1.157 °C；
- 最差单风速区域 RMSE 约 1.622 °C；
- Tmax RMSE 约 1.637 °C；
- DeltaT RMSE 约 1.648 °C；
- 热点位置始终靠近出口，热点趋势门槛失败。

COMSOL 确认的目标不是制造更优结果，而是确认同一等效参数化在有限元环境中可复现。

## 七、可辨识性结论必须保留

V9 当前结果：

- 雅可比条件数约 183；
- 最大参数相关系数约 0.992；
- `beta` 和 `eta` 强相关；
- Bootstrap 中 `eta` 置信区间很宽并跨零。

因此最终报告必须明确：

- 区域平均温度可以较好拟合；
- 流量偏置和热源偏置不能唯一分离；
- 真实几何未恢复；
- 热点位置未验证；
- 不能用于训练高保真 ROM 或做决定性 MPC 结论。

## 八、最终文件

生成：

- `outputs_v9_comsol/comsol_confirmation.csv`
- `outputs_v9_comsol/comsol_vs_python.csv`
- `outputs_v9_comsol/fig7_validation_metrics.json`
- `outputs_v9_comsol/final_summary.json`
- `docs_v9_comsol/implementation_audit.md`
- `docs_v9_comsol/final_decision_report.md`

最终回答：

1. COMSOL 是否正确实现了 V9 冻结等效模型；
2. COMSOL 与 Python 是否达到数值一致性；
3. Fig.7 区域温度误差是多少；
4. 热点趋势是否仍失败；
5. 能否恢复真实通道和歧管参数；
6. 是否允许启动 ROM、MM-EKF 和 Oracle MPC。

只有全部高保真门槛通过才允许进入后续控制研究；预期当前热点门槛仍不通过，因此不得自动授权。
