# Yuan 2020 PEMFC 三节点模型复现 v2

本工程完成了：

1. Fig. 14、15、16 曲线数字化；
2. `C_th`、`K_cool`、`K_node`、`Kp`、`Ki` 的正则化辨识；
3. 论文曲线与复现曲线叠加图；
4. 校准集与独立验证集分离；
5. 参数冻结和哈希记录；
6. COMSOL 6.4 等效三维建模脚本；
7. COMSOL 三区域温度验证管线；
8. 五节点、双传感器 EKF 和热点约束 MPC 原型。

## 首先查看

- `docs/completion_report_CN.md`
- `identification/frozen_parameters_v2.json`
- `outputs/02_fig14_calibration_overlay.png`
- `outputs/03_fig16_calibration_validation_overlay.png`
- `comsol/README_COMSOL_CN.md`

## 已验证的复现精度

- Fig. 14 总体 RMSE：0.258 °C
- Fig. 16 独立验证 RMSE：0.306 °C

## 重要状态

三节点控制导向模型已经达到“曲线级、分离验证的论文功能复现”。

COMSOL 脚本已经生成，但由于当前环境没有 COMSOL，尚未实际求解 `.mph`。
