# V8 / V7 能量分配诊断

## V7 六风速分配

| 风速 (m/s) | Q_air/Q_total | Q_natural/Q_total | Q_radiation/Q_total | Q_conduction/Q_total | 外部散热合计 |
|---:|---:|---:|---:|---:|---:|
| 4 | 30.01% | 56.44% | 13.56% | 0% | 69.99% |
| 5 | 34.10% | 53.24% | 12.66% | 0% | 65.90% |
| 6 | 37.35% | 50.69% | 11.96% | 0% | 62.65% |
| 8 | 42.23% | 46.85% | 10.92% | 0% | 57.77% |
| 10 | 45.74% | 44.05% | 10.18% | 0% | 54.23% |
| 12 | 48.46% | 41.90% | 9.62% | 0% | 51.52% |

4 m/s 时外部散热是通道空气散热的2.33倍；即使到12 m/s，外部散热仍略高于强制空气冷却。这不是由能量不守恒造成的：V7 最大能量残差只有0.02794%。根因是 `h_external=30 W/(m²·K)` 和 `emissivity=1.0` 均被4、8 m/s校准推到上界，外部边界实际承担了连续1 mm空气层没有解释的冷却能力。

## V8 冻结外部边界

- 自然对流不再拟合，采用 Churchill–Chu 竖直平板相关式，特征长度0.249369 m并按局部表面温度计算；55°C参考表面对应 `h_external=4.859 W/(m²·K)`。[Sandia Aria相关式说明](https://www.sandia.gov/files/sierra/Aria_Users_5_20/simulation_setup/advanced_feat/correlation_heat_transfer_coeff_reference.html)
- 石墨总发射率固定为0.82，是《Industrial Graphite Engineering Handbook》给出的0.70–0.95范围中点，不用5、6、10、12 m/s调整。[IAEA Graphite Handbook](https://nucleus.iaea.org/sites/graphiteknowledgebase/Meetings2/Old%20Meetings/2017/Background%20Info/GraphiteHandbook.pdf)

V8 结构模型中，外部散热只由上述相关式和固定发射率决定。结构模型仍无法同时复现温度水平与8–11°C沿程温差，因此不能再把误差转移给外部散热参数。
