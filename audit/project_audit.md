# V2 工程完整性审计

审计对象：`Yuan2020_full_reproduction_v2/Yuan2020_full_reproduction_v2`。审计遵循先检查、后修改原则；V2 原目录未被修改。其 51 个文件已逐文件 SHA-256 一致地复制到独立 v3 仓库，基线提交为 `98af306`，标签为 `baseline_v2`。

## 结论

V2 是“冻结数据和已有结果包”，不是从原图到控制实验均可执行的研究工程。现有三节点函数可运行，并能从已有 cleaned CSV 复算 Fig16 的 700-1200 s RMSE 0.306303 °C；但数字化、辨识、五节点、EKF 和 MPC 均缺少真实可执行入口，因此已有 PNG/CSV/JSON 不能从当前源码完整重建。

## 审计项目

1. `requirements.txt` 不足：版本只有宽松下限，未包含 `PyYAML`、`pytest`，也未锁定解释器；本机 Python 3.13.5 环境缺少 Matplotlib。需增加可复现环境记录与依赖锁定。
2. 所有 15 个 CSV、5 个 JSON、7 个 PNG 和 13 个 JPG 均可读取；原始三图尺寸分别为 Fig14 1634x1221、Fig15 1634x1225、Fig16 3374x2962。
3. `frozen_parameters_v2.json` 的内部哈希有效：删除哈希字段后，以 UTF-8、`ensure_ascii=False`、缩进 2、无末尾换行序列化，SHA-256 为 `8f669482afc5818feba0019673b07063f713e3ebe2b6b7ea31608f1533d84262`。
4. `models/three_node_model.py` 的函数可独立导入和调用，但没有 CLI；时间步被隐含为 1 s，控制参数大量硬编码，能量项未输出。
5. `digitization/digitize_figures.py` 只打印已有冻结数据位置，是占位脚本。
6. 工程没有参数辨识程序，只有辨识结果文件；不能复算多起点优化、协方差、相关性或置信区间。
7. `models/five_node_observer_mpc.py` 只读取并打印已有 JSON，是占位脚本；不存在可执行五节点、EKF 或 MPC 实现。
8. V2 COMSOL Java 文件在已安装 COMSOL 6.4.0.293 上实际编译失败：`main` 未处理 `model.save` 抛出的 `IOException`。
9. Java 文件在结果表达式中调用 `maxop1(T)`，但没有创建 `maxop1`。
10. `run_comsol_windows.bat` 指向不存在的默认 C 盘安装路径，并把 `.java` 直接交给 `comsolbatch -inputfile`；COMSOL 6.4 要求先用 `comsolcompile` 生成 `.class`，再由 `comsolbatch` 执行 `.class`。脚本还包含阻塞式 `pause`。
11. Java 几何选择依赖自动生成的 `geom1_air_dom`、`geom1_solid_dom`，尚未经过成功运行验证；热传递接口、耦合、出口无黏性应力、外表面散热和能量积分均不完整。
12. 论文表 3 的 `A_inlet=3962 mm²` 是控制模型冷却入口总截面积，而 `A_node=1420 mm²` 是节点控制体截面积；V2 又令 `L*W=247 cm²` 并构造单个均质空气块，这属于等效几何假设，不能称为原 CAD。
13. 论文表 2 的空气黏度写作 `1.849x10^-5 m²/s`，单位对应运动黏度；V2 将数值直接写入 `dynamicviscosity`，量纲错误。等效模型应明确采用约 `1.849e-5 Pa*s` 的动态黏度，或用运动黏度乘密度换算，并记录选择。
14. `validate_comsol_reduction.py` 没有调用三节点模型，只计算 COMSOL 三区温差统计。
15. 当前输出不可完整重建：数字化、辨识、EKF、MPC 的生成代码缺失，COMSOL 未编译和求解。

## 论文边界核对

已对原 PDF 的第 4-10 页进行文本提取和页面渲染核查。论文明确使用独立空气冷却流道、稳态 Laminar Flow 与 Heat Transfer 3D 模型，以及沿冷却流道离散的三节点控制模型。不能将其扩展解释为开放阴极反应空气-冷却空气耦合模型，也不能从该模型推出氧传输、水淹或膜含水量结论。

## 审计布尔状态

详见 `audit/project_audit.json`。只有三节点数值函数被判定为可执行；COMSOL 编译、求解及所有端到端复现状态均为 false。
