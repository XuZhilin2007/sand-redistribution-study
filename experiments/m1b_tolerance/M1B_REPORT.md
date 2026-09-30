# M1B 实验报告：感知容差 / Deadband 停止

> 日期：2026-09-18
>
> 模型定义：[docs/M1B_TOLERANCE.md](../../docs/M1B_TOLERANCE.md)；前置：[M1A_MODEL.md](../../docs/M1A_MODEL.md)、[M1A_LONG_HORIZON.md](../../docs/M1A_LONG_HORIZON.md)；代码：`src/sand_m0/adaptive.py`（未修改——M1B 只是同一 runner 的不同 stop_tolerance）；配置：[config.json](config.json)（运行前冻结）
>
> 前置事项：本轮开始前完成了一项数据审计（D_max 数值口径冲突），结论为 canonical 数据一致、系完成报告把 t=99 值（0.0497）误写为终值（t=100 实为 3.09×10⁻²）。勘误见 [ERRATUM_2026-09-18_M1A_DMAX.md](../../docs/ERRATUM_2026-09-18_M1A_DMAX.md)，冻结文件未修改。

## 1. 研究问题

> 把 exact-zero stopping 改为有限 tolerance 后，能否形成更稳定、分辨率不敏感、同时更接近原始现实操作动机的 feedback stopping rule？

三个子问题：(1) 有限容差是否消除 exact-stop 的 K scaling？(2) 停止得太早还是太晚？(3) 是否存在宽而稳定的 near-optimal tolerance region？

冻结设计：epsilon ∈ {0（=10⁻¹²）, 0.005, 0.01, 0.02, 0.04} × K ∈ {50, 100, 200, 400}，T_max=5000，alpha=0.25，q_near。5×4=20 个 deterministic case，结构实验而非参数优化。

## 2. 数学结构（先于运行推导）

- **First-passage 恒等（模型推导数学结果）**：M1B(epsilon) 与 exact reference 在停止前逐步相同，故 `t_eps = min{t : D*_max(t) ≤ epsilon}` 是 reference D_max 序列的首达时间；停止态 = reference 在 t_eps 的状态。**运行前即可由 reference 序列预测全部 epsilon 的停止轮**——20/20 个 case 的 canonical 运行与预测一致。
- **epsilon 单调性**：`epsilon₁ ≤ epsilon₂ ⟹ t_{eps₂} ≤ t_{eps₁}`（集合包含，平凡但严格）；200 点对数网格 × 4 K 数值验证全部通过。
- **分辨率极限**：epsilon ≥ 0.005 的首达时间 K 稳定（数值）；未证明定理（Working Hypothesis，见 [M1B_TOLERANCE.md](../../docs/M1B_TOLERANCE.md) §2.4）。

## 3. 结果一：有限容差消除 K scaling（问题 1：是）

| epsilon | K=50 | K=100 | K=200 | K=400 | 跨 K 比值 |
|---:|---:|---:|---:|---:|---:|
| 0（=10⁻¹²） | 324 | 699 | 1445 | 2876 | **8.9×** |
| 0.005 | 253 | 340 | 338 | 321 | 1.34× |
| 0.01 | 165 | 159 | 173 | 169 | **1.09×** |
| 0.02 | 76 | 80 | 76 | 82 | 1.08× |
| 0.04 | 34 | 35 | 34 | 34 | **1.03×** |

**有限容差彻底移除 exact-stop 的 K 依赖**：epsilon ≥ 0.005 时停止轮跨 K 变化 ≤1.34 倍（epsilon ≥ 0.01 时 ≤1.09 倍），而 exact 标准 8.9 倍。过渡相当陡：0.5% 的容差就足以坍缩大部分 K 敏感性。epsilon=0.005 在 K=50 处有轻微残余依赖（253 vs ~330）——容差与 1/K 离散尺度可比时 K 依赖部分回归，与理论预期一致。

## 4. 结果二：停止太早（问题 2）

Reference 的 U 最优时刻 t_U_best = 292 / 653 / 1349 / 2627（K=50/100/200/400，自身 ∝K）；U_density_best = 3.71 / 3.82 / 3.69 / 3.56 ×10⁻³。

| epsilon | delta_t（K=100） | delta_t（K=400） | regret 范围（K=50..400） |
|---:|---:|---:|---|
| 0 | +46 | +249 | 0.85 – 1.43 |
| 0.005 | −313 | −2306 | **0.82 – 1.65** |
| 0.01 | −494 | −2458 | 2.00 – 2.51 |
| 0.02 | −573 | −2545 | 3.41 – 4.45 |
| 0.04 | −618 | −2593 | 7.35 – 7.90 |

- **exact 规则轻微过修正**（delta_t = +46@K=100，regret 0.85–1.43，即停止态 ≈2×U 最优——与 M1A.1 的 1.85 倍一致）；
- **所有有限 epsilon 都停得过早**（delta_t < 0，且差距随 K 增大：t_U_best∝K 而 t_stop(eps) 与 K 无关）；
- **regret 随 epsilon 单调上升**：没有任何冻结 epsilon 同时做到"早停"与"接近 U 最优"。

公平的语境：regret 是模型内相对量。M1B 在 K=100、epsilon=0.005 的停止态 raw U ≈ 8.5×10⁻⁵，仍是 M0 最佳（1.68×10⁻³）的 1/20。

## 5. 结果三：robust tolerance region（问题 3）

**分 verdict 回答**：

- **分辨率稳健性上：存在 robust region**。epsilon ≥ 0.005 的停止轮跨 K 稳定（≤1.34×；≥0.01 时 ≤1.09×），且停止轮数大幅缩短（699→159–340@K=100）。
- **均匀度近优性上：不存在**。区域内所有 epsilon 的 regret ≥ 0.82（最好情形 epsilon=0.005），且 regret 随 epsilon 单调上升；唯一的低 regret 点（exact，低 K）恰好不属于有限容差的稳健区。

**结构性原因（模型依赖观察）**：D_max（sup-of-CDF）衰减快、U（L2）改善慢且持续到 ∝K 的时刻——两个指标的时间尺度错位，使"容差停止"与"U 近优停止"在本模型中渐近不相容。**这不是缺陷，而是 M1B 的核心发现**：perceptual-tolerance 停止天然是"足够好就停"，不是"最优才停"。

**结论措辞（按治理）**：M1B has a robust tolerance region for resolution-robust stopping (epsilon ≥ 0.005, 尤其 ≥0.01)，但 uniformity regret 在该区域内不趋于零；不声称任何 epsilon 是"最优"或"视觉阈值"。

## 6. Controls

q_uniform 与 q_far（K=100、epsilon=0.01）：均在 t=0 立即停止（D_max(0) = 1.11×10⁻¹⁶ / 0.00）——M1B 未破坏"无近侧过剩即无动作"的正确行为。

## 7. Deadband / perceptual 连接（文献定位）

- **Deadband connection**（Literature-supported modeling connection，标准控制/质量工程概念）：偏差处于动作限内不执行调整——与 M1B 的 `D_max ≤ epsilon → 停止` 结构相同。未做控制论文献检索原文核实，仅记录概念联系。
- **Perceptual motivation**：真实观察是"没有明显值得修正的区域就停止"（Real-world Observation）；`D_max ≤ epsilon` 是其最小代理（Model Assumption）。epsilon 未由人类实验测量：**不得称 JND、不得称视觉阈值、不得声称某 epsilon 是真实最优**；psychophysical calibration 留给未来。

## 8. Evidence Classification

- **Real-world Observation**：Owner 实际操作中"没有明显值得处理的区域就停止"（承自 checkpoint）。
- **Model Assumption**：`D_max ≤ epsilon` 停止代理；冻结 epsilon 集；epsilon=0 用 10⁻¹² 实现。
- **Model-derived Mathematical Result**：first-passage 恒等（M1B 停止轮 = reference 序列首达，停止态 = reference 状态）；epsilon 单调性（含 200 点网格 × 4 K 验证）。
- **Verified Simulation Result**：20 个 (epsilon, K) canonical case 的全部数值（停止轮、D_max@stop、U@stop、regret、near-half mass、停止前边界）；controls 立即停止。
- **Model-dependent Observation**：epsilon ≥ 0.005 消除 K scaling；有限 epsilon 一律早停且 regret 随 epsilon 单调上升；稳健区存在但不含低 regret 点。
- **Working Hypothesis**：固定正 epsilon 的 t_eps 在 K→∞ 收敛到有限值；D_max 与 U 的时间尺度错位是停止-最优权衡的结构根源。
- **Literature-supported modeling connection**：deadband/action-limit 概念（未做原文核实检索，仅概念联系）。
- **无 Real-world calibrated result**。

## 9. Unexpected Results

1. **K-scaling 坍缩得非常陡**：0.5% 容差（相对 D_max 末段水平 ~0.02–0.05 并不小）就把 8.9× 压到 1.34×——原本预期需要 epsilon 与离散误差同量级的细致分析，实际一档就够。
2. **exact 规则的 regret（0.85–1.43）与最小有限容差（0.82–1.65）相当**：容差并没有在 regret 上吃亏太多，真正的分歧在更高 epsilon；而 eps=0.005 在 K=50 上甚至比 exact 更接近 U 最优（因为 exact 在 K=50 也轻微过修正）。
3. **t_U_best 本身 ∝K**：U 的最优点也随分辨率线性后移——"何时最均匀"与"何时停止"都带 K 足迹，但只有后者被容差消除了。
4. 无其他异常；controls 与全部单调性/恒等测试一次通过。

## 10. What This Does NOT Prove

- 不证明 epsilon 的任何取值对应人类感知（无校准）；
- 不证明 t_eps(epsilon) 的 K→∞ 极限存在（Working Hypothesis）；
- 不证明 regret-epsilon 关系在 q、alpha、K 之外泛化（单 q、单 alpha）；
- 不证明"robust region"在其它分布下存在（q_near only）；
- 不引入 hysteresis/noisy perception/visual smoothing/2D/替代 q（等待 Owner + GPT review M1B）；
- epsilon=0.005–0.04 的任何数值无现实操作含义。

## 11. Reproducibility

- **Python** 3.14.6（仓库 `.venv`）；**依赖** numpy 2.5.3、matplotlib 3.11.2。
- **命令**：`.venv/Scripts/python.exe experiments/m1b_tolerance/run_experiment.py`
- **测试**：`tests/test_m1b_tolerance.py`（20 项）+ 既有五套（37+25+53+18+21 项）全部通过。
- **随机性**：全部 deterministic，无种子。
- **结果文件**（`experiments/m1b_tolerance/results/`，SHA-256 见 metadata.json）：`epsilon_k_summary.csv`（20 case 全量）、`first_passage_monotonicity.csv`、`controls.csv`、三张 PNG、`metadata.json`。
- **冻结声明**：epsilon 集、K 集、T_max、阈值实现方式均运行前写入 config.json；无事后调参、无事后追加 epsilon。

## 12. Recommended Next Research Question

按本轮证据，最有信息量的下一问（待 Owner + GPT 决定）：

1. **停止判据与 U 最优的错位能否用第二个量修复？**——D_max 快、U 慢的时间尺度错位是 regret 的根源；可研究以"U 的变化率"或混合判据作为停止信号（仍属停止规则研究，不引入视觉机制）；
2. 或回到 **Phase 2B**（Beta family / r=U² 上检验 M0 判据与 M1B 容差行为），把 M0–M1B 的结论链在分布维度闭合。
