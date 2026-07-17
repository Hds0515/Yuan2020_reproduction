# V16研究范围与数据注册

V15标签 `v15_fixed60_reference_confirmed` 保持不变，冻结文件SHA-256为
`05566ac5903af63fdda10c1dd941a3bfbd6b7f607409007d341e320debf81ae7`。V16没有覆盖任何V15文件。

当前可辩护范围仅为23°C环境、固定60% PWM、2–40 A，并使用实测电流和
电堆电压作为外生输入。Guo 2024的2 A阶梯已用于辨识，4 A阶梯已经执行过
一次最终确认；仓库中不存在尚未使用的新动态实验，因此不能再次把4 A数据
称为独立确认。

`high_resolution_dataset_registry.json`和`two_mass_dataset_registry.json`中的
数据均明确标记为冻结V15父模型模拟。它们只用于ROM和观测器软件方法验证。
下一次真实实验使用`fresh_experiment_registry_TEMPLATE.csv`，必须在采集前冻结
calibration、validation和final-confirmation用途，并在采集后填写SHA-256。
