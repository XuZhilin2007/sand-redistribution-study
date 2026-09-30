# M1B.1 诊断文档：均匀度指标与观察尺度审计

> 状态：纯诊断（2026-09-18）。M1A/M1B dynamics、adaptive boundary policy、epsilon stopping rule、throw distribution、alpha、controller **一概未改**；观察尺度 B 只是诊断仪器，从不进入 policy。
>
> 实验：[M1B_1_REPORT.md](../experiments/m1b_1_metric_scale/M1B_1_REPORT.md)；前置：[M1B_TOLERANCE.md](M1B_TOLERANCE.md)、[M1A_LONG_HORIZON.md](M1A_LONG_HORIZON.md)
>
> 核心问题：M1B 的"停止太早 / uniformity regret 很大"究竟是真实的控制器—目标错位，还是 residual-regret 的小分母问题 + fine-grid 观察尺度错配？

## 0. 两个治理事项（先于结果）

1. **测试计数勘误**：M1B 完成报告写"共 154 项"，其自身加数（20+37+25+53+18+21）和为 **174**——算术错误。实测当前运行时断言数为 **172**（37/25/54/18/18/20；M0.2 与 M1A.1 两套因轮内测试修订与早先报告相差 ±1–3，均为各轮 commit 内记录的编辑）。纯计数勘误，不影响任何实验数据或结论。
2. **Deadband 文献连接升级**：dead-band adjustment 文献引用已经 Crossref API 独立核实（2026-09-18）：Alberto Luceño, "Ch. 19. Dead-band adjustment schemes for on-line feedback quality control", *Handbook of Statistics* Vol. 22, Elsevier, 2003, pp. 695–727, DOI `10.1016/S0169-7161(03)22021-6`。**provenance**：书目元数据经 Crossref 独立验证；dead-band 结构描述（周期测量偏差、仅在超出 action limits 时调整）来自章节标题与 GPT 提供的阅读，本轮未读全文。记录为 **Literature-supported modeling connection**；不支持 human perception、不支持任何 epsilon 数值、不构成 novelty。（已同步至 [M1B_TOLERANCE.md](M1B_TOLERANCE.md) §7。）

## 1. 方法

- 复用 canonical dynamics 重新生成 exact M1A reference 状态史（K ∈ {50,100,200,400}，T=5000，tol=10⁻¹²）；回归锚点：停止轮与已提交 M1B 结果 20/20 一致，K=100 轨迹与已提交 M1A.1 轨迹一致（t=0..100）。
- 观察尺度：把 K 个内部 bin 精确聚合为 B 个 macro-bin（B ∈ {5,10,25,50}，均整除所有 K，无插值），`U_B = B·Σ_b(P_b − 1/B)²`。B=K 时该式**精确退化为**既有 fine-grid `U_density`（已测试）。B=5/10/25/50 分别对应约 20%/10%/4%/2% 面积的观察区。
- 停止态质量三读数：absolute_gap = U_stop − U_best；residual_regret = (U_stop − U_best)/U_best；**improvement_loss = (U_stop − U_best)/(U_initial − U_best)**；retention = 1 − loss。不预设任何"near-optimal"阈值（90/95/99% 线仅为读图辅助）。

## 2. 核心发现一：fixed-B 的 U 最优时刻并不收敛——漂移属于动力学

| 观察尺度 | t_best(K=50) | t_best(K=100) | t_best(K=200) | t_best(K=400) | 跨 K |
|---|---:|---:|---:|---:|---|
| B=5（20% 区） | 190 | 602 | 1111 | 730 | 非单调（平坦） |
| B=10（10% 区） | 292 | 561 | 1363 | 2627 | ≈9.0× |
| B=25（4% 区） | 292 | 653 | 1349 | 2627 | ≈9.0× |
| B=50（2% 区） | 292 | 653 | 1349 | 2627 | ≈9.0× |
| fine（B=K） | 292 | 653 | 1349 | 2627 | ≈9.0× |

**fixed-B（B≥10）的 t_best 与 fine-grid 完全同步 ∝K**——"U 最优随 K 漂移"不是细网格伪影，而是动力学性质：粗到 10% 宽的观察区也在持续改善到 ∝K 的时刻。B=5（20% 区）例外且非单调：最粗尺度的 U_B 在 t≈200 后已近乎平坦（argmin 由微小波动决定，190/602/1111/730）。

**证据等级**：Verified Simulation Result。

## 3. 核心发现二：residual regret 是小分母伪影；retention 全面修正 M1B 结论

K=100 的 improvement retention（每个观察尺度；K=200/400 形状相同）：

| epsilon | B=5 | B=10 | B=25 | B=50 | fine (B=K) |
|---:|---:|---:|---:|---:|---:|
| 0（exact） | 0.9911 | 0.9913 | 0.9907 | 0.9908 | 0.9901 |
| 0.005 | 0.9948 | 0.9951 | 0.9935 | 0.9907 | 0.9857 |
| 0.01 | 0.9935 | 0.9882 | 0.9845 | 0.9777 | 0.9740 |
| 0.02 | 0.9951 | 0.9885 | 0.9710 | 0.9628 | 0.9605 |
| 0.04 | 0.9766 | 0.9514 | 0.9200 | 0.9142 | 0.9128 |

K=400 同型（exact 0.989–0.993；eps=0.005 ≥0.9868；eps=0.01 ≥0.9729；eps=0.02 ≥0.9520；eps=0.04 ≥0.9148）。

**同一个停止态，两种读数天差地别**：K=100、epsilon=0.04、B=5 的 residual_regret = **103.2**（U_best 仅 7.2×10⁻⁵），而 improvement retention = **0.9766**——即"扔掉了可获得改善的 2.3%"却被 residual 口径读成"差 103 倍"。M1B 报告中 0.82–7.9 的 regret 数值全部属于此类小分母放大。

**对 M1B 结论的正式修正（Case A）**：

> M1B 的停止状态在 improvement-retention 语义下**存在宽而稳定的 near-optimal region**：epsilon ≤ 0.02 时所有观察尺度、所有 K 的 retention ≥ 0.96；epsilon = 0.005 时 ≥ 0.986；即便 epsilon = 0.04 仍有 0.91–0.98。原先"uniformity near-optimality 没有 robust region"的说法过强，其高 regret 数值主要来自 residual 口径的小分母（U_best → 0）。

**证据等级**：Verified Simulation Result（retention 数据）+ Model-dependent Observation（near-optimal 区域的刻画）。

## 4. 核心发现三：epsilon 的尺度依赖温和且方向一致

同一 epsilon 的停止时刻不变（controller 未改——回归锚点验证），但停止质量随观察尺度**温和变化、方向一致**：coarse 尺度的 retention 略高于 fine（eps=0.04：B=5 为 0.977 vs fine 0.913，差 6 个百分点；eps ≤ 0.02 时各尺度差 ≤2.5 个百分点）。解读边界：M1B 的停止**不是**"在所有尺度都明显太早"（eps ≤ 0.01 时连细网格都保留 ≥97%）；也**不能**解读为"肉眼大尺度已经足够均匀"——B 只是诊断尺度，未与视觉校准。

## 5. Multiscale 汇总（diagnostic-only）

`U_multiscale` 类的跨 B 均值 retention 已作为诊断列输出（`multiscale_retention_summary.csv`：按 (K, epsilon) 给出 B∈{5,10,25,50} 上 retention 的 mean/min/max）。不升级为 primary metric，不用于优化 epsilon。

## 6. 证据等级汇总

- **模型推导数学结果**：B=K 时 U_B 与 fine-grid U_density 精确同一（同一公式）；精确聚合保持总质量。
- **已验证模拟结果**：16 个 (K,B) 的 reference 最优表；95 个 (K,eps,B) 的三读数；两个回归锚点（M1B 停止轮 20/20、K=100 轨迹 vs M1A.1）。
- **模型依赖观察**：fixed-B t_best ∝K（B≥10）；B=5 平坦例外；retention 的 near-optimal 区域（epsilon ≤ 0.02 / ≥0.96）；coarse-fine retention 差的方向与幅度。
- **工作猜想**：人式粗观察（类似 B≈5–10 的尺度）下，有限容差停止与"接近最优"几乎不可区分——这可能是真实操作中"无明显过剩即停"在功能上成立的原因；待未来感知校准。
- **Literature-supported terminology/modeling connection**：one-sided KS functional（承 M1A.1）；dead-band adjustment（Luceño 2003，Crossref 核实）。

## 7. What This Does NOT Prove

- 不证明任何 B 对应人类视觉分辨率或 JND；
- 不证明 retention ≥ 95% 的停止"足够好"——本轮不定义成功阈值；
- 不证明 coarse-fine retention 差在其它 q/alpha 下保持（单 q、单 alpha）；
- 不证明 fixed-B t_best ∝K 的机制（未做解析推导）；
- 不修改 controller、不引入 smoothing/hysteresis/noisy perception——全部等待 Owner + GPT review M1B.1。
