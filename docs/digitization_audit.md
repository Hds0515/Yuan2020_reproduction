# 曲线数字化重建审计

## 数据来源与隔离

v3 的曲线由 `source_images/Fig14.jpg`、`Fig15.jpg`、`Fig16.jpg` 直接生成。提取程序不读取 V2 cleaned CSV；只有在三个 regenerated CSV 已写出后，`validate_digitization.py` 才把它们与 V2 CSV 作对照。原始 JPG、V2 raw CSV 与 V2 cleaned CSV 均未覆盖。

## 算法

1. 使用 `axis_calibration.json` 的数据范围和图框作为初始标定。
2. 按每条曲线的 RGB 原型计算颜色距离，并用饱和度阈值排除黑色坐标轴、文字和白色背景。
3. 使用 `retained_y_pixels` 排除图例。Fig16 只保留主图 (a)，不使用放大图 (b)(c)。
4. 保存每个通过阈值的原始像素及其颜色距离。
5. 对每列以颜色距离倒数加权中位数取得线芯，仅剔除孤立且与两侧轨迹不连续的点。
6. 只对缺失列插值，再转换到数据坐标并重采样为 1 s；没有从 V2 CSV 平移、缩放或滤波。

## Fig16 标定修正

V2 的 Fig16 像素矩形 `[608, 10, 1562, 711]` 对应约 2048 像素宽的缩放副本，而受保护原图宽 3374 像素。直接应用会切到错误区域。v3 在原图上重新检测黑色边框，主图 (a) 的矩形为 `[1006, 18, 2572, 1171]`，并在配置中保留了这一修正说明。

## 像素分辨率不确定度

以下为纯量化的半像素下限，不包含 JPEG 压缩、抗锯齿、线宽和曲线重叠；实际曲线读数建议至少按约 3 像素线宽评估。

| 图 | x 分辨率 | 半像素 x 不确定度 | y 分辨率 | 半像素 y 不确定度 |
|---|---:|---:|---:|---:|
| Fig14 | 0.82079 s/pixel | ±0.41040 s | 0.018639 °C/pixel | ±0.009320 °C |
| Fig15 | 0.82305 s/pixel | ±0.41152 s | 0.032588 A/pixel | ±0.016294 A |
| Fig16 | 0.76628 s/pixel | ±0.38314 s | 0.017346 °C/pixel | ±0.008673 °C |

按 3 像素线宽计，温度轨迹的读数尺度约为 Fig14 ±0.056 °C、Fig16 ±0.052 °C，电流轨迹约为 Fig15 ±0.098 A。曲线重叠处的结构性误差可能更大，因此参数辨识同时报告 bootstrap 区间与相关性，不能把半像素值当作完整置信区间。

## 可审计输出

- `digitization/regenerated/Fig14_digitized.csv`
- `digitization/regenerated/Fig15_digitized.csv`
- `digitization/regenerated/Fig16_digitized.csv`
- `digitization/regenerated/intermediate/*_pixel_points.csv`：所有颜色匹配像素。
- `digitization/regenerated/intermediate/*_preclean.csv`：每列原始线芯、匹配数、剔除标记。
- `digitization/regenerated/intermediate/*_data_coordinates.csv`：清洗后但尚未 1 s 重采样的数据坐标。
- `outputs/digitization/*_pixel_overlay.png`：回贴原图检查图。
- `outputs/digitization/digitization_metrics.json`：有效列比例、匹配像素数和阈值。
- `outputs/digitization/v2_comparison.json`：只在提取完成后生成的 V2 对照。

## 当前限制

Fig16 的多条平滑仿真线在约 55 °C 附近重叠，颜色核心会被上层线遮挡，因而紫色中点仿真曲线的直接有效列比例约 25%；缺失段由相邻同色像素连续插值。覆盖图显示轨迹与原图一致，但该区间的结构性不确定度高于不重叠曲线，后续结论不得忽略这一点。
