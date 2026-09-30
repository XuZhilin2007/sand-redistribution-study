# M1B.1 实验报告：Uniformity Metric & Observation-Scale Audit

> 日期：2026-09-18
>
> 分析文档：[docs/M1B_METRIC_SCALE_AUDIT.md](../../docs/M1B_METRIC_SCALE_AUDIT.md)；前置：[M1B_TOLERANCE.md](../../docs/M1B_TOLERANCE.md)；代码：`src/sand_m0/adaptive.py`（未修改）；配置：[config.json](config.json)（运行前冻结）
>
> 性质：纯诊断轮。成功标准不是证明 M1B 好或坏，而是确定此前"高 regret"究竟代表真实丢失大量均匀化收益，还是 residual 分母与 fine-grid 观察尺度的表象。

## 1. 开始前的两项治理审计

1. **测试计数勘误**：M1B 完成报告"共 154 项"为算术错误（其加数和为 174）；实测当前运行时断言 172（37/25/54/18/18/20，个别套件因轮内测试修订与早先记录差 ±1–3）。纯计数勘误，不影响数据与结论，已记录于 RESEARCH_LOG。
2. **Deadband 文献连接**：Crossref API 独立核实引用（Luceño 2003，*Handbook of Statistics* Vol. 22, pp. 695–727, DOI 10.1016/S0169-7161(03)22021-6）；书目元数据独立验证，结构描述来自章节标题 + GPT 阅读与本轮未读全文的边界已注明。记录为 Literature-supported modeling connection（不支持 human perception / epsilon 数值 / novelty）。已同步至 M1B_TOLERANCE.md §7 与 RESEARCH_LOG。

## 2. 方法

- 复用 canonical dynamics 重新生成 exact reference 状态史（未改模型）；回归锚点：停止轮 vs 已提交 M1B 结果 20/20 一致；K=100 轨迹 vs 已提交 M1A.1 一致（t=0..100）。
- 固定观察尺度 B ∈ {5,10,25,50}（均整除所有 K，精确聚合无插值）：`U_B = B·Σ_b(P_b − 1/B)²`；B=K 时精确退化为 fine-grid `U_density`（测试验证）。fine-grid metric 保留不动。
- 三读数：absolute_gap、residual_regret（小分母警告）、**improvement_loss / improvement_retention = 1 − loss**（相对"初始→该尺度最优"的总可获改善）。
- 测试：`tests/test_m1b_1_metric_scale.py`（17 项）通过；六套既有测试（37+25+54+18+18+20 项）保持通过。

## 3. 结果一：fixed-B 最优时刻不收敛——漂移属于动力学

| 观察尺度 | K=50 | K=100 | K=200 | K=400 |
|---|---:|---:|---:|---:|
| B=5（20% 区） | 190 | 602 | 1111 | 730 |
| B=10（10% 区） | 292 | 561 | 1363 | 2627 |
| B=25（4% 区） | 292 | 653 | 1349 | 2627 |
| B=50（2% 区） | 292 | 653 | 1349 | 2627 |
| fine（B=K） | 292 | 653 | 1349 | 2627 |

B=10/25/50 的 t_best 与 fine-grid **完全同步 ∝K**；B=5 例外（平坦、非单调：20% 区的 U_B 在 t≈200 后改善幅度极小，argmin 由微小平移决定）。**"U 最优随 K 漂移"不是细网格伪影——固定粗观察尺度下动力学同样持续改善到 ∝K 的时刻。** U_best(K,B) 本身随 K 下降（更细的模拟在同等粗尺度上做到更均匀，如 B=10：1.17→0.56×10⁻³）。

## 4. 结果二：retention 全面修正 residual-regret 读数（K=100；K=200/400 同型）

| epsilon | B=5 | B=10 | B=25 | B=50 | fine |
|---:|---:|---:|---:|---:|---:|
| 0（exact） | 0.9911 | 0.9913 | 0.9907 | 0.9908 | 0.9901 |
| 0.005 | 0.9948 | 0.9951 | 0.9935 | 0.9907 | 0.9857 |
| 0.01 | 0.9935 | 0.9882 | 0.9845 | 0.9777 | 0.9740 |
| 0.02 | 0.9951 | 0.9885 | 0.9710 | 0.9628 | 0.9605 |
| 0.04 | 0.9766 | 0.9514 | 0.9200 | 0.9142 | 0.9128 |

同一停止态的两种读数对照（K=100、eps=0.04、B=5）：residual_regret = **103.2**，improvement_retention = **0.9766**。residual 口径除以趋零的 U_best，把"扔掉 2.3% 可获改善"放大成"差 103 倍"。

## 5. 对 M1B split verdict 的重新判断（Case A）

> **正式修正**：M1B 原结论"uniformity near-optimality 没有 robust region"过强。在 improvement-retention 语义下存在宽而稳定的 near-optimal region：epsilon ≤ 0.02 时全部 (K, B) 的 retention ≥ 0.96；epsilon = 0.005 时 ≥ 0.986；即便 epsilon = 0.04 仍达 0.91–0.98。原高 regret 数值主要来自 residual 口径的小分母。

保留的 nuance（Case B 成分，如实记录）：retention 在细网格略低于粗尺度（eps=0.04：fine 0.913 vs B=5 0.977）——细尺度结构被牺牲得略多；epsilon=0.005 的 retention 在部分 (K,B) 上甚至高于 exact 停止（exact 停在 U 最优之后，轻微过修正）。

## 6. Scale dependence of epsilon

停止时刻与 K、B 无关（controller 未改，锚点验证）；停止质量随 B 温和变化且方向一致：coarse ≥ fine，差幅从 eps ≤ 0.02 的 ≤2.5 个百分点到 eps=0.04 的 ~6 个百分点。回答任务问题：M1B 的停止**不是**"在所有尺度都明显太早"——eps ≤ 0.01 时连细网格都保留 ≥97.4%；但也不能声称"肉眼大尺度已足够均匀"（B 未与视觉校准）。

## 7. Multiscale 汇总（diagnostic-only）

`multiscale_retention_summary.csv` 给出按 (K, epsilon) 的 B∈{5,10,25,50} retention mean/min/max。例：K=100、eps=0.005 → mean 0.9910 [0.9857, 0.9951]；eps=0.04 → mean 0.9550 [0.9128, 0.9766]。仅作趋势参考，不升级为 primary metric，不用于优化 epsilon。

## 8. Evidence Classification

- **Model-derived Mathematical Result**：B=K 时 U_B 与 fine-grid U_density 精确同一；精确聚合保持总质量（测试验证）。
- **Verified Simulation Result**：16 个 (K,B) 最优表、95 个 (K,eps,B) 三读数、两个回归锚点、retention 的跨 K/跨 B 矩阵。
- **Model-dependent Observation**：fixed-B（B≥10）t_best ∝K（漂移属动力学）；B=5 平坦例外；retention near-optimal 区域（eps ≤ 0.02 → ≥0.96 全尺度）；coarse-fine retention 差的方向。
- **Working Hypothesis**：人式粗观察尺度（B≈5–10 量级）下，有限容差停止与"接近最优"几乎不可区分——可能是真实操作"无明显过剩即停"在功能上成立的原因；待感知校准。
- **Literature-supported terminology/modeling connection**：dead-band adjustment（Luceño 2003，Crossref 核实）；one-sided KS functional（承 M1A.1）。
- 无 Real-world calibrated result；无 novelty 声明。

## 9. What This Does NOT Prove

- 不证明任何 B 对应人类视觉分辨率/JND；
- 不定义"停止足够好"的阈值（90/95/99% 线仅为读图辅助）；
- 不证明 retention 关系跨 q/alpha/其他控制器泛化（单 q、单 alpha、M1A 边界策略）；
- 不证明 fixed-B t_best ∝K 的解析机制；
- 不实现 U-derivative stopping、combined stopping、optimized epsilon、hysteresis、smoothing、noisy perception、coarse boundary selection——全部等待 Owner + GPT review M1B.1。

## 10. Reproducibility

- **Python** 3.14.6（仓库 `.venv`）；**依赖** numpy 2.5.3、matplotlib 3.11.2。
- **命令**：`.venv/Scripts/python.exe experiments/m1b_1_metric_scale/run_experiment.py`
- **测试**：`tests/test_m1b_1_metric_scale.py`（17 项）+ 既有六套（37+25+54+18+18+20）全部通过。
- **随机性**：deterministic，无种子（np.random 仅用于聚合恒等测试的合成状态）。
- **结果文件**（`experiments/m1b_1_metric_scale/results/`，SHA-256 见 metadata.json）：`fixed_scale_optimum.csv`、`stop_quality_by_scale.csv`（95 行）、`fixed_B_t_best_vs_K.csv`、`multiscale_retention_summary.csv`、四张 PNG、`metadata.json`。
- **可复现实证**：干净树上重跑，数据输出位级一致（metadata 哈希比对）。
- **冻结声明**：B 集、阈值、读数定义运行前写入 config.json；未修改任何 controller 或既有结果文件。

## 11. Recommended Next Research Question

诊断显示真正的取舍不是"停止早 vs 晚"，而是**观察尺度**：细尺度均匀化收益在停止后仍在缓慢积累（∝K），粗尺度收益早已饱和。待 Owner + GPT 决定：

1. **感知尺度的正式化（谨慎）**：若接受"操作者按粗尺度评价"，M1B 的 epsilon ≤ 0.02 已在功能上近优——下一步是把"粗尺度评价"作为显式模型成分（evaluation-scale 与 controller 解耦的正式化），而非继续调 epsilon；
2. **Phase 2B**（Beta family / r=U²）：检验 M0 判据、M1B 容差与本次尺度结论在分布维度是否闭合。
