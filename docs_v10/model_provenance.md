# V10 文献约束参考模型：来源与参数审计

## 1. 模型身份

V10 是替代性研究基线，不是 Yuan 2020 未公开 CAD 的重建，也不改变 V7–V9 的负面结论。模型以 Shahsavari 等人的公开空冷 PEMFC 模型和实验表格为主要来源，建立一个守恒的代表性阴极通道。论文原文及本次存档副本的 SHA-256 为
`3E5CC71356AEAC3F726DCCB212D74842B59E7EB0DCF2BD21E20012539C47E768`。

主要来源：Shahsavari, Desouza, Bahrami, Kjeang, *Thermal analysis of air-cooled PEM fuel cells*, International Journal of Hydrogen Energy 37 (2012) 18261–18271, [DOI: 10.1016/j.ijhydene.2012.09.075](https://doi.org/10.1016/j.ijhydene.2012.09.075)，[作者公开 PDF](https://www.sfu.ca/~mbahrami/pdf/2012/S.%20Shahsavari%2C%20A.%20Desouza%2C%20M.%20Bahrami%2C%20E.%20Kjeang%20-%20Thermal%20analysis%20of%20air-cooled%20PEM%20fuel%20cells.pdf)。

## 2. 参数分类

### 论文直接给定

- 80 个阴极通道；电池宽 280 mm；流向长度 60 mm；
- 阴极通道高 2.5 mm、梯形长边 2.5 mm、侧壁角 80°；
- 双极板厚 5.5 mm，阳极/阴极 GDL 各 0.2 mm，CCM 厚 0.05 mm；
- 面内导热率：双极板 60 W/(m·K)，GDL 10 W/(m·K)，CCM 1.5 W/(m·K)；
- Table 1 的六组入口速度、总电池热量、实验 Tmax 与论文模拟 Tmax；入口温度 21°C；
- 论文说明完整模型包含入口/出口 plenum，并使用实验局部电流密度与电压作为输入；这些空间输入未在 Table 1 中公开。

### 文献补充闭合关系

- 完全发展矩形管、均匀热流的 Shah–London Nusselt 相关式；声明的等效矩形宽高比为 1，得到 `Nu=3.610224`。相关式形式可核查于[矩形微通道换热研究](https://www.sciencedirect.com/science/article/pii/S2214157X21000198)。
- 21°C 附近干空气常物性：`rho=1.204 kg/m³`、`Cp=1006 J/(kg·K)`、`k=0.02514 W/(m·K)`、`mu=1.825e-5 Pa·s`。
- 方形管近似的 Darcy `f Re=56.91` 仅用于压降诊断，不参与温度验证。

### 由上述量确定的几何量

- 通道节距 3.5 mm；
- 梯形短边 1.618365 mm；
- 单通道入口面积 `5.147956e-6 m²`；
- 湿周 `9.195498e-3 m`；
- 水力直径 `2.239338e-3 m`；
- 代表通道轴向 `kA=1.1692625e-3 W·m/K`。

### 校准参数

无。六组实验温度未用于拟合 `Nu`、几何、材料参数、热量或空气物性。

### 未知或未公开参数

- 入口/出口 plenum 的全部 CAD 和逐通道流量分配；
- Table 1 六工况的 96 点局部电流密度输入；
- 实验夹具、端板、接触热阻和外部漏热的完整边界；
- 六工况的 24 点温度原始矩阵；
- Yuan 2020 Fig.7 的真实内部通道和歧管几何。

热接触电阻不能随意忽略为“已知真实值”。Sadeghifar 等人的研究表明，BPP–GDL 接触热阻会随压紧力、PTFE/MPL 和表面平整度显著变化：[Journal of Power Sources 273 (2015) 96–104](https://www.sciencedirect.com/science/article/pii/S0378775314014700)。本阶段没有用该不确定项拟合验证集。

## 3. 冻结与可追溯性

全部常数、派生量、验证划分和源文件哈希记录在 `outputs_v10/frozen_parameter_registry.json`；其 SHA-256 为
`bbb5818382a7d25f33eebe70f3f1f29dc0cda5ca2892583e089922dd6fdef9d0`。Python 和 COMSOL 读取同一物理定义，任何后续修改必须建立新版本，不得覆盖 V10。
