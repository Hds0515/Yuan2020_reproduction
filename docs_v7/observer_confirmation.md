# V7 门控 MM-EKF 预注册确认

## 状态：预注册框架冻结，确认实验未执行

V6 的门控 MM-EKF 在固定 V5 压力测试上通过数值门槛，但属于看到 V5 失败模式后的后验开发。V7 要求使用新的高分辨率确认场景；由于 COMSOL 物理验证失败，这些场景不存在，因此不能声称通过 V7 预注册确认。

`outputs_v7/observer_confirmation_registry.json` 冻结了拟沿用的 V6 传感器位置、模型库和门控阈值，并明确记录空的高分辨率场景集合和 `blocked_by_physical_validation` 状态。没有优化传感器、调整模型库、调整门控阈值或重新抽取 seed。

V7 结论：门控 MM-EKF仍是候选创新，但尚未成为经独立预注册确认的核心创新。
