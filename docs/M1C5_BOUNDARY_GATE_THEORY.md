# M1C.5 理论文档：Boundary-Gated Rebound Theory——churn 切换门的解析化

> 状态：理论 + existing-trajectory verification 轮（2026-09-19）。M1C.4 的直接后续；**无新模型、无新参数、无新实验**——只对既有 canonical Variant A 确定性轨迹（K∈{50,100,200,400}）做推导核验。
>
> 实验：[M1C5_REPORT.md](../experiments/m1c5_boundary_gate/M1C5_REPORT.md)；前置：[M1C4_CHURN_REGIME_AUDIT.md](M1C4_CHURN_REGIME_AUDIT.md)（a-门控 observation）、[M1C3_REDISTRIBUTION_LAW_ORDER.md](M1C3_REDISTRIBUTION_LAW_ORDER.md)（精确反弹判据）、[M1C2_TARGET_MATCHED_KERNEL_THEORY.md](M1C2_TARGET_MATCHED_KERNEL_THEORY.md)
>
> North Star：把 M1C.4 的经验门（a<0.55 必收缩、a≥0.8 必反弹、apparent gate ≈0.58）提升为**由 exact rebound criterion 推出的解析可能边界**。

## 0. 裁决（Outcome B1：analytic gate 成功）

**边界门控被解析证明为必要条件，且闭式门槛定量落在观察门上**：

1. **Region-1 化简（精确）**：在扫掠前缀上，反弹残差严格化为
   `M·Δ(j) − γ(j) = (1−α)·D(j) − D_max + α·h(j/K)`，其中 `h(x) = x[F_a(2−x)−1]`、`F_a = F(a) = a + D_max`。全部 5344 个 active 轮核验，最大残差 **4.0×10⁻¹⁶**。
2. **State-reduced 必要条件（定理）**：由 `(1−α)D(j) ≤ (1−α)D_max`，
   `M·Δ(j) − γ(j) ≤ α·(h(j/K) − D_max)`，故 **rebound ⟹ max_{j≤j\*} h(j/K) > D_max**——只依赖 (a, D_max)，不依赖完整 D-profile。
3. **解析最大化**：h 是开口向下抛物线。`F_a ≤ 1/2 ⟹ h ≡ 0`（精确，rebound 不可能）；`F_a > 1/2` 时内部最大在 `x\* = 1 − 1/(2F_a)`，值 `(2F_a−1)²/(4F_a)`。x\* ≤ a 的条件是 `(a+D_max)(1−a) ≤ 1/2`，对 canonical A 全部状态成立（D_max ≤ 0.25 < √2−1；实测 0 轮违反）。
4. **闭式门槛**：`(2F_a−1)²/(4F_a) > D_max` 的解给出
   **a_crit(D) = [(1−D) + √(D(D+2))]/2**（大根；小根落在 F_a≤1/2 惰性分支）。
   **rebound 在 a ≤ a_crit(D_max) 时数学上不可能。**
5. **经验验证**：全部 1958 个 rebound 轮满足 a > a_crit（**FN = 0**，定理要求）；a_crit 的 churn-window 分布 p50 = **0.593–0.604**（跨 K 几乎不变）、p5–p95 = 0.539–0.644——**解析门槛自然落在 M1C.4 观察门 ≈0.58 区域**。
6. **必要 ≠ 充分（混合区解释）**：gate-open 但仍收缩的轮每 K 有 21–177 个；其反弹驱动前缀的 profile 因子 P(x_r) = D(x_r)/D_max 分布显著区别于 rebound 轮——**边界位置开启可能门，D-profile 形状决定是否真的反弹**。

## 1. Canonical specialization（离散精确关系）

Variant A：T(j) = j/K；G_near 离散 CDF；Δ(j) = G(j) − j/K。a = j\*/K 是 D 的 argmax ⟹ 精确关系：

```text
D(a) = D_max        （argmax 定义）
F(a) = a + D_max    （F = D + T 在 j* 处）
M = alpha·F(a) = alpha·(a + D_max)
```

（离散形式；continuum 记号仅为表达。）

## 2. Region-1 化简（引理 1，精确）

Region 1（j ≤ j\*）的 target-matched 更新（M1C.2）：`D′_TM(j) = (1−α)D(j) − α·T(j)·(1−F_a)`。代入 γ(j) = D_max − D′_TM(j) 与 M·Δ(j) = α·F_a·(j/K)(1−j/K)：

```text
M·Δ(j) − γ(j) = (1−α)·D(j) − D_max + α·h(j/K)
h(x) = x[F_a(2−x) − 1]        （Region 1，精确）
```

（Prompt 候选公式确认正确。）suffix（j > j\*）由 M1C.3 永不反弹，故必要条件只需 Region 1。核验：5344 轮 × 全前缀，残差 ≤ 4.0×10⁻¹⁶。

## 3. State-reduced 必要条件（定理，引理 2）

D(j) ≤ D_max ⟹ (1−α)D(j) − D_max ≤ −α·D_max，故

```text
M·Δ(j) − γ(j) ≤ α·(h(j/K) − D_max)    （j ≤ j*，精确上界）
```

**rebound（∃j 残差 > 0）⟹ max_{j≤j\*} h(j/K) > D_max。** 这是一个 **state-reduced necessary condition**：只依赖 (a, D_max) 两个标量。上界在 D(j) = D_max（即 j = j\*）时最紧，但 j\* 处残差恒负（M1C.3：等价于 G(j\*)>1 不可能）——bound 不因这一点失效（仍是有效上界）。

## 4. 解析最大化（引理 3）

h(x) = (2F_a−1)x − F_a·x²，h(0) = 0，开口向下：

| 情形 | max_{[0,a]} h | 结论 |
|---|---|---|
| F_a ≤ 1/2 | 0（h 单调降，h(0)=0） | 0 > D_max 假 ⟹ **rebound 不可能（精确）** |
| F_a > 1/2，x\* ≤ a | h(x\*) = (2F_a−1)²/(4F_a)，x\* = 1−1/(2F_a) | 必要条件 (2F_a−1)²/(4F_a) > D_max |
| F_a > 1/2，x\* > a | h(a) = a[D(2−a) − (1−a)²] | 见下 |

**x\* ≤ a 条件**：x\* ≤ a ⟺ (a+D_max)(1−a) ≤ 1/2。该式对 a 的最大值 ((1+D)/2)² 在 **D ≤ √2−1 ≈ 0.414** 时 ≤ 1/2。canonical A 全程 D_max ≤ 0.25（初值即全轨迹最大，实测 4 K 全部 0 轮违反）⟹ **x\* ≤ a 恒成立，内部最大公式适用于每一轮**。

（x\* > a 分支的 h(a) = a[D(2−a)−(1−a)²]：因 a(2−a)−1 = −(1−a)² ≤ 0，该值 ≤ a·D(2−a) < D 当 a < 1 且 D < D/(1)……直接验证：h(a) > D 要求 a(2−a) > 1 − D(2−a)/a，而 2a−a² ≤ 1 恒成立（(1−a)² ≥ 0），故 h(a) ≤ a·D(2−a)；对 a<1、D≤0.25 有 a(2−a) ≤ 1，故 h(a) ≤ D(2−a)，h(a)>D 需要 2−a>1 即 a<1 且 a·D(2−a)>D——只在 a→1 极限时接近。该分支在 canonical 轨迹上不出现（x\*≤a 恒成立），记录完整性而已。）

## 5. 闭式门槛（定理，主结果）

`(2F_a−1)²/(4F_a) > D_max`，F_a = a + D_max。解等式 `(2(a+D)−1)² = 4D(a+D)`：

```text
4a² + 4(D−1)a + (1−4D) = 0
a = [(1−D) ± √(D(D+2))] / 2
```

**大根即门槛**：`a_crit(D) = [(1−D) + √(D(D+2))]/2`。小根 a = [(1−D)−√(D(D+2))]/2 落在 F_a < 1/2 的惰性分支（该分支 h≡0，rebound 无论 a 如何都不可能），无 gate 意义。φ(F) = (2F−1)²/(4F) 在 F > 1/2 严格递增 ⟹ 必要条件等价于 **a > a_crit(D_max)**。

**逻辑方向（重要）**：

```text
a ≤ a_crit(D_max)  ⟹  rebound 不可能          （已证明）
a > a_crit(D_max)  ⟹  rebound 可能，未必发生    （necessary only）
```

**数值样例**：a_crit(0.005) = 0.548，a_crit(0.02) = 0.590，a_crit(0.025) = 0.600，a_crit(0.05) = 0.635，a_crit(0.25) = 0.750。D 越小门越低——D_max 收得越紧，越小的前缀也能产生反弹。

## 6. 离散 vs 连续诚实记录

正式陈述是**离散**的：rebound ⟹ max_{j≤j\*} h(j/K) > D_max（格点上取值）。闭式 a_crit 是 **continuum-style envelope**，与离散格点 max 的差别有两种 regime：

- 一般（x\* 离格网近）：O(F_a/(4K²))——实测 K=100 worst 2.4×10⁻⁵、K=200 6.2×10⁻⁶；
- x\* 落在最前格点之前（F_a→1/2⁺，x\* < 1/K，如 K=50 t=314：F_a=0.5016、x\*=0.0032）：格点 max 已在 h 下降段，gap 达 O(1/K)——实测 K=50 worst 1.4×10⁻⁴。

两种 regime 的 gap 都远小于 D_max 量级（0.005–0.25），**对本轮全部结论（FN=0、门槛位置）无实际影响**；离散精确条件以 max_j h(j/K) 为准（`state_reduced_gate` 同时返回两者）。

## 7. Existing trajectory verification（Verified Numerical Analysis）

| 检查 | 结果 |
|---|---|
| Region-1 化简残差 | ≤ 4.0×10⁻¹⁶（5344 轮全前缀） |
| State-reduced bound 违反 | 0（逐轮逐前缀） |
| F_a ≤ 1/2 分支 h ≤ 0 | 精确成立 |
| x\* ≤ a | 0 轮违反（4 K 全部轨迹） |
| **必要门槛假阴性** | **FN = 0**（1958 rebound 轮全部 a > a_crit） |
| gate 曲线下收缩轮（被理论正确排除） | churn window 内 170 / 379 / 809 / 1624（各 K 全部收缩轮减去 gate-open 收缩） |

### 跨 K a_crit 分位数（churn window）

| K | p5 | p25 | p50 | p75 | p95 | rebound frac | FN | gate-open 收缩 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 0.5528 | 0.5823 | 0.6043 | 0.6222 | 0.6438 | 0.3654 | 0 | 21 |
| 100 | 0.5451 | 0.5784 | 0.5992 | 0.6203 | 0.6387 | 0.3624 | 0 | 52 |
| 200 | 0.5418 | 0.5733 | 0.5952 | 0.6165 | 0.6366 | 0.3692 | 0 | 88 |
| 400 | 0.5391 | 0.5692 | 0.5926 | 0.6136 | 0.6349 | 0.3687 | 0 | 177 |

**与 0.58 的关系**：a_crit p50 ≈ 0.59–0.60、整个 p5–p95 带 0.54–0.64——**解析门槛 tracks the observed gate**（M1C.4：a<0.55 零反弹、[0.55,0.6) 近零、0.8+ 全反弹）。措辞：analytic necessary boundary tracks the observed gate；**不是** "theorem proves rebound threshold equals 0.58"。

## 8. 必要 vs 充分与混合区（§10/§11 的回答）

- **gate-open contraction 轮**（B_t > 0 却收缩）：每 K 21/52/88/177 个（占 gate-open 轮的 16%/15%/13%/13%）。它们证明 necessary ≠ sufficient。
- **Profile 因子**（§11）：反弹驱动前缀 x_r（残差 argmax，恒在扫掠前缀内）处 P(x_r) = D(x_r)/D_max。正常 churn 带（D_max > 0.005）内：rebound 轮 P 分布集中于 [−1, 1] 高段（h 峰位 prefix 的 D 相对高度决定余量被穿越），gate-open contraction 轮显著散布于负 P 区（D(x_r) 低于均匀份额 ⟹ (1−α)D 项贡献为负 ⟹ 即使 gate 开着也补不回 margin）。`fig_profile_factor.png`。
- **§10 same-a 反例**：中间区 [0.55, 0.8) 的混合行为由 (a, D_max) 二维_gate（a 相同、D_max 不同 ⟹ B 不同）+ profile 形状共同解释——单变量 a 不可能充分。

## 9. Evidence Classification

- **Model-derived Mathematical Result**：Region-1 化简恒等式；state-reduced 上界；h 的解析最大化（含 F_a≤1/2 精确不可能分支与 x\*≤a 条件）；闭式门槛 a_crit(D) 及"rebound ⟹ a > a_crit(D_max)"必要条件定理。
- **Verified Numerical Analysis**：5344 轮全部恒等式/bound/门槛核验（FN=0、残差 ≤4×10⁻¹⁶）；a_crit 分位数表；P(x_r) 分布；确定性复现（输出位级一致）。
- **Model-dependent Observation**：a_crit 分布跨 K 稳定且集中于 0.54–0.64（track 观察门）；gate-open contraction 的占比与 P 分布；D 越小门越低的结构。
- **Working Hypothesis**：反复的边界漂移穿越该门产生长期 churn switching regime（本轮不证明 boundary-return dynamics）。

## 10. What this adds to the churn story

```text
边界位置 a（相对 a_crit(D_max)）→ rebound 是否"可能"        [定理]
完整 D-profile（x_r 处 D 高度）→ 可能性是否兑现为反弹        [精确残差判据 + P 因子观察]
收缩后边界为何重新漂回门上                                   [未解决——下一阶段]
```

## 11. Limitations

1. 必要条件非充分（by design；充分条件未寻找——Prompt §12 允许停止在 necessary + exact criterion）；
2. x\* ≤ a 恒成立依赖 D_max ≤ √2−1（canonical 满足；一般状态需逐轮检查）；
3. closed-form 是 continuum-style envelope（离散两种 rounding regime 已记录）；离散精确条件用格点 max；
4. 单一 q/alpha/一维/deterministic/4 K；a_crit 分布的跨 K 稳定性是 4 点观察；
5. P(x_r) 是简单诊断（deep-valley 轮比值爆炸，已过滤 D_max > 0.005 并记录）；
6. 不证明 boundary-return dynamics、不证明长期 churn。

## 12. Next Research Gate（交回 GPT + Owner）

1. **boundary-return dynamics**：为什么收缩后边界重新漂回门上（C-C-R 节奏的"回到门上"半边）；
2. **充分条件/半减小界**：gate-open 轮中预测哪些 rebound（P 因子的解析化）；
3. **长期化**：门控 + 切换结构是否足以支撑 A 不收敛定理；
4. 依据 Prompt §21，推导 + 验证完成后停止。

## 13. Repository 资产

- `src/sand_m0/kernel_theory.py` 扩展：`rebound_residual` / `state_reduced_gate` / `a_crit_gate`（pure functions，含定理表述）；
- `experiments/m1c5_boundary_gate/`：config / run_analysis.py / per_round_gate.csv（5344 行）/ gate_summary.csv / 3 图 / M1C5_REPORT.md；
- `tests/test_m1c5_boundary_gate.py`（10 项）。
