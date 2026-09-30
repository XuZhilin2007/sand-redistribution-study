# Erratum：M1A / M1A.1 的 D_max 数值口径（2026-09-18）

> 性质：报告口径勘误（reporting erratum）。**canonical 数据本身完全一致，无需修改任何结果文件。**
>
> 审计范围：`experiments/m1a_adaptive_boundary/results/m1a_baseline_trajectory.csv`、`experiments/m1a_1_long_horizon/results/long_horizon_trajectory_K100.csv`、M1A 与 M1A.1 的两份完成报告及已提交文档。

## 1. 冲突描述

- M1A 完成报告（对话记录）称：T=100 运行的 "D_max 终值 ≈ 0.0497"。
- M1A.1 完成报告与检查点表称：t=100 的 D_max = 3.0919×10⁻²。
- M1A.1 同时声明其 K=100 轨迹与 M1A canonical 轨迹前 101 行逐位一致。

三者表面冲突。

## 2. 审计结果

| 事实 | 数值 |
|---|---|
| M1A canonical CSV（101 行，t=0..100）row t=99 的 D_max | 4.966640989722404×10⁻² |
| M1A canonical CSV row t=100（最终状态）的 D_max | 3.091936013113139×10⁻² |
| M1A.1 K=100 CSV 相同行的 D_max | 与上表一致至 12 位有效数字（该 CSV 存储格式为 12 位科学计数；两文件相对差最大 4.4×10⁻¹³，纯属存储精度差） |
| U_L2 / a_t / U_density 逐行对照 | 同样一致（a_t 精确相等；U_L2 值相等，仅字符串格式不同） |

**结论：canonical data 一致。** 差异只是：

1. **时间索引口径**：0.0497 是 **t=99**（倒数第二个状态）的 D_max；t=100（最终状态）的 D_max 是 **3.0919×10⁻²**。M1A 完成报告把 t=99 的值误写为"终值"。
2. **存储格式**：M1A CSV 用全精度 repr，M1A.1 CSV 用 12 位科学计数——数值相同。

## 3. 勘误与影响

1. 正确数值（后续所有文档统一采用）：**D_max(t=100, K=100, canonical) = 3.0919×10⁻²**；D_max(t=99) = 4.9666×10⁻²。
2. M1A 完成报告中"D_max 终值 ≈ 0.0497"应读作"D_max(t=99) ≈ 0.0497"。
3. 已提交文档中"波动衰减至 ≈0.05"的表述（M1A_MODEL.md、M1A_REPORT.md、RESEARCH_LOG.md）描述的是末段量级（t=99 处 0.0497），作为阶段量级大致成立，但**最终状态**是 0.031——以本勘误为准。
4. `docs/M1A_LONG_HORIZON.md` 第 19 行原句"T=100 时 D_max 仍为 0.0497"已就地加勘误标记修正为 t=100 → 3.09×10⁻²（t=99 才是 0.0497）。
5. 按治理，冻结的历史 M1A report（`experiments/m1a_adaptive_boundary/M1A_REPORT.md`）与其结果文件**不做修改**；本勘误文档与 RESEARCH_LOG 记录为准。
6. 该勘误不影响 M1A/M1A.1 的任何结论（不停止的 T=100 观察、T=5000 停止触发等均与精确数值无关）。

## 4. 后续口径

- 引用 T=100 运行的最终 D_max 时，一律写 3.09×10⁻²（t=100）；
- 引用末段波动水平时，注明"t≈99 附近 D_max≈0.05"；
- 行索引约定（重申）：canonical M1A/M1A.1 CSV 的 row t 对应**状态 t**（t 次更新后），其 D_max/a_t 为该状态上的决策诊断。
