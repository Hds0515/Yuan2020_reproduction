# COMSOL 6.4 等效三维模型

模型名称固定为 `Yuan2020_equivalent_3D`。论文未公开完整 CAD 和流道尺寸，所以它不是作者原始几何的精确复刻。

`Yuan2020Equivalent3D.java` 已在 `E:\COMSOL\COMSOL64\Multiphysics` 的 COMSOL 6.4 环境中实际编译成功，并完成几何、材料、物理场、网格、区域算子和参数化研究的 API 初始化。最后一次 Stokes 初始化稳态求解在 12% 处长期振荡，人工安全中止；没有得到收敛 MPH 或可信的场量导出。

运行 `run_comsol_windows.bat` 可重新编译并求解。仅当 `comsol_zone_temperatures.csv` 存在且包含六个有限工况时，`validate_comsol_reduction.py` 才执行真实三节点稳态比较；否则它返回状态码 2，绝不生成虚构验证值。

当前可靠状态以 `status.json` 为准。`logs/solve_twelfth.log` 是最后一次求解证据。当前不得声称：求解完成、网格收敛、CFD 能量守恒通过、温度云图已生成或三节点交叉验证通过。
