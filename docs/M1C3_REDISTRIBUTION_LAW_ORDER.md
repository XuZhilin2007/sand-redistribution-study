# M1C.3 理论文档：Redistribution-Law Order & Instability Criterion——CDF 序定理与精确单步反弹判据

> 状态：数学推广轮（2026-09-19）。M1C.2 的直接后续；**theory + existing-trajectory verification，无新模拟、无新变体、无参数扫描**。全部推导沿用 M1C.2 的 canonical 离散记号。
>
> 前置：[M1C2_TARGET_MATCHED_KERNEL_THEORY.md](M1C2_TARGET_MATCHED_KERNEL_THEORY.md)（分解恒等式与 target-matched 定理）、[M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md](M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md)
>
> North Star：把"near-biased q 会重新注入误差"升级为精确数学陈述——**哪些 redistribution laws 保证每一步稳定改善；不满足该序关系的 laws，在什么精确条件下 biased reinjection 压过 controller 本来的纠偏，造成 D_max 反弹。**

## 1. Motivation

M1C.2 已证明（Model-derived Mathematical Result）：`G = T` 时每个 active step 严格收缩 D_max，且实际更新相对 target-matched counterfactual 的**唯一**额外项是失配再注入 `M·Δ(j)`，`Δ(j) := G(j) − T(j)`。对 canonical q_near vs 均匀目标 `Δ(x) = x(1−x) ≥ 0`（内部恒正）。本轮回答三个问题：

1. **序定理**：`G ≤ T`（逐前缀）是否足以保证每个 active step 的 D_max 严格下降？
2. **精确反弹判据**：`G ≰ T` 时，单轮 D_max 上升的 exact necessary-and-sufficient 条件是什么？
3. **验证**：该判据能否逐轮精确预测现有 canonical A 轨迹中约 36% 的回升轮？

## 2. Canonical decomposition（承 M1C.2，未改）

记号与 M1C.2 完全一致：离散 K bins、Σp=1、T(j)=j/K（一般化时 0≤T(j)≤1、T(K)=1）、D(j)=F(j)−T(j)、D_max=max_j D(j) ≥ 0、j\* 为 argmax（最小下标）、active ⟺ D_max > tol（canonical 1e-12）、M = alpha·F(j\*)、G(j)=Σ_{i≤j}g_i：

```text
D′(j) = D′_TM(j) + M·Δ(j)        （两个 region 同式）
Region 1 (j ≤ j\*): D′_TM(j) = (1−alpha)D(j) − alpha·T(j)·(1 − F(j\*))
Region 2 (j >  j\*): D′_TM(j) = D(j) − M(1 − T(j))
```

## 3. CDF-order theorem（本轮主定理）

> **定理（cumulatively target-dominated redistribution）**
>
> 设 G、T 为离散 CDF（非降、取值 [0,1]、G(K)=T(K)=1），T(j) < 1 对一切 j < K；alpha ∈ (0,1)；当前步 active（D_max > 0，canonical 中由 tol ≥ 0 自动保证）。若
>
> `G(j) ≤ T(j)` 对每个前缀 j 成立
>
> 则 **D_max(t+1) < D_max(t)（strict）**。
>
> **术语**：满足 `G ≤ T` 的 law 称为 **cumulatively target-dominated redistribution**（累计目标受控再分配：任何前缀收回的份额都不超过目标分配给它的份额）。Uniform target 下的 G=T（M1C.2 情形）与 G_far(x)=x²≤x 都是该类的成员——后者仅作 CDF 序的示例说明，**本轮无任何 far-biased 模拟**。

**证明**：active 步 M = alpha·F(j\*) > 0（F(j\*) ≥ j\*/K > 0）。由 G ≤ T 得 M·Δ(j) ≤ 0，故 D′(j) ≤ D′_TM(j) 对一切 j。再由 M1C.2 定理的逐点形式：Region 1 的 D′_TM(j) ≤ (1−alpha)D_max < D_max；Region 2 的 D′_TM(j) ≤ D_max − M(1−T(j)) < D_max 对 j < K（用 T(j) < 1），j=K 处 D′_TM(K) = 0 < D_max。三项合并：每个前缀 D′(j) < D_max ⟹ max_j D′(j) < D_max。∎

**最小假设清单**（避免未检查假设）：不需要 T 严格递增、不需要 bins 严格为正、不需要光滑性、无 continuum limit；strictness 只在退化目标（某 j < K 使 T(j)=1）处失效——canonical 均匀目标以余量满足 T(j)=j/K<1。

**Scope note（对 M1C.2 的澄清，非更正）**：M1C.2 定理的 strict 部分同样隐含需要 T(j) < 1（j < K）；其 canonical 实例（均匀目标）不受影响。

## 4. 三个结构性推论（均为精确结论）

1. **Suffix 不可能性**：对**任意** fixed law（无论 Δ 符号），Region 2 永不产生反弹：`D′(j) = D(j) − M(1−G(j)) ≤ D(j) ≤ D_max`（因 G(j) ≤ 1）。**D_max 反弹只可能发生在被扫前缀内部（j ≤ j\*）**。推论：反弹轮的新 argmax 必落在扫掠前缀内——在 A 全部实际反弹轮上验证成立（§7）。
2. **前缀 margin 下界**：γ(j) ≥ alpha·D_max 对一切 j ≤ j\*（由 D(j) ≤ D_max 与 T(j) ≥ 0）。于是**任何反弹的必要条件**：∃ j ≤ j\* 使 M·Δ(j) > alpha·D_max；再由 Δ(j) ≤ 1−T(j) ≤ 1−1/K，均匀目标下必要条件化为 **M > alpha·D_max·K/(K−1)**。
3. **旧 argmax 处不可能反弹**：在 j = j\* 处，`M·Δ(j\*) > γ(j\*)` 等价于 `G(j\*) > 1`——不可能。反弹永远来自非 argmax 的更短前缀。

## 5. Exact one-step rebound criterion（本轮第二个主结果）

定义 pointwise contraction margin `γ(j) := D_max − D′_TM(j)`（active 步上逐点为正）。由分解式：

> **单步反弹判据（exact iff）**
>
> `D_max(t+1) > D_max(t)  ⟺  ∃ j：M·Δ(j) > γ(j)`
>
> 即：本轮 biased redistribution 在某个前缀重新塞回的累计过剩，**超过**该前缀在 target-matched 更新下本来能获得的纠偏幅度。

**判据的精确性**：这是恒等变形（max 泛函的定义），无附加条件。配合 §4.1，满足判据的 j 必在扫掠前缀内。注意这是**单步**判据——不得称为 long-run instability（见 §11 边界）。

## 6. Global sufficient bound 与幅度界

记 `Δ₊ := max_j max(Δ(j),0)`、全局 margin `Γ := D_max − max_j D′_TM(j) = min_j γ(j) > 0`（active 步）：

- **充分不增条件**：`M·Δ₊ ≤ Γ ⟹ D_max′ ≤ D_max`（strict `<` 时严格下降）。
- **幅度上界**：`D_max′ − D_max ≤ M·Δ₊ − Γ`。

这两个界给出"再注入强度 vs 纠偏余量"的最简竞争关系，但**非**精确判据（精确判据是逐点的 γ）。

## 7. Canonical q_near specialization

`Δ(x) = x(1−x)`，网格峰值 Δ₊ = 1/4（j = K/2；在全部 4 个 canonical 网格上精确成立）。故每轮最大可能正再注入为 **M/4**——它是判据的"可能性上界"：`M/4 ≤ Γ` 的轮必然不反弹；`M/4 > Γ` 只说明存在反弹的**可能**，实际是否反弹由逐点 γ(j) 决定。

## 8. Existing A trajectory verification（本轮核心验证）

**无新模拟**——直接在 canonical Variant A 既有轨迹（确定性重算，位级与已提交结果一致）上逐轮验证：

| K | active 轮 | 实际反弹轮 | 判据预测轮 | FP | FN |
|---:|---:|---:|---:|---:|---:|
| 50 | 324 | 117（36.1%） | 117 | 0 | 0 |
| 100 | 699 | 251（35.9%） | 251 | 0 | 0 |
| 200 | 1445 | 531（36.8%） | 531 | 0 | 0 |
| 400 | 2876 | 1059（36.8%） | 1059 | 0 | 0 |
| **合计** | **5344** | **1958** | **1958** | **0** | **0** |

数值残差：判据 profile（D′_TM + MΔ）与实际下一状态 D_max 之差 ≤ **2.0×10⁻¹⁵**；与 one_step_excess 恒等式之差 ≤ 2.8×10⁻¹⁷。判据逐轮 iff 匹配（浮点 tolerance 1e−12，未触发边界情形）。

同时验证的精确预测：全部反弹轮的新 argmax 都在扫掠前缀内（§4.1）；suffix 前缀在全部 5344 轮从不升破旧 D_max；global bound 与幅度界逐轮成立；γ(j) ≥ alpha·D_max 抽样成立。

## 9. Diagnostic decomposition（解释性，非参数研究）

全部 5344 轮按实际反弹/下降分类，在**实际新 argmax** 处取值：

| 类别 | n | M 中位 | Δ@new 中位 | γ@new 中位 | R = MΔ/γ 中位 | R 范围（p5–p95） |
|---|---:|---:|---:|---:|---:|---:|
| 反弹轮 | 1958 | 0.2123 | 0.2436 | 0.0266 | **1.70** | 1.05 – 5.87（全部 > 1） |
| 下降轮 | 3386 | 0.1087 | 0.1688 | 0.0319 | **0.55** | 0.16 – 0.91（全部 < 1） |

解读：反弹轮的特征组合是**更大的动作**（M ≈ 0.21 vs 0.11）、**失配峰值位置的新 argmax**（Δ ≈ 0.244 ≈ 峰值 0.25）与**更小的局部纠偏余量**（γ ≈ 0.027）。R = MΔ/γ 在两类轮上的分离是判据的逐点重述（在新 argmax 处 R>1 ⟺ 该点破 D_max），与理论完全一致；真正的预测力在于判据对全部 1958/3386 轮的零误分。反弹占比跨 K 稳定在 0.36–0.37。

## 10. Counterexamples（converse 不成立的精确刻画）

- **Counterexample A（G > T somewhere 仍可逐步收缩）**：构造近端微倾斜 law（G>T 但 Δ₊ < 10⁻³），在 A 的全部抽样状态上手工执行 canonical 步：判据判"无反弹"、实际 D_max 逐轮严格下降。**正失配的存在只是该位置产生反弹的必要条件（对该 j 而言），不是充分条件**——target-matched 收缩可以完全吸收小幅失配。
- **Counterexample B（G ≤ T 反弹）**：按定理不可能。tiny-K 确定性扫描（K∈{2,3,4}，2955 个 manual active steps，多个 G≤T laws 含 uniform/far/混合/随机受控 law）：无一例 D_max 上升。属于数学核验（deterministic check），不构成研究轨迹声明。
- **未证明的方向**："若某 law 使 D_max 在一切状态的一切 active 步都不反弹，则 G ≤ T"——**未建立**（本轮不声称 iff）。已证明的 iff 只在单步判据层面（§5）。

## 11. 与 long-run churn 的边界（重要）

**one-step rebound criterion ≠ long-run churn theorem.** 本轮判据精确回答"什么条件下某一轮 D_max 反弹"，它不推出：无限多次反弹、不变带、周期轨道、不收敛性、exact stopping ∝K。A 的长期 churn（M1B.2 实测的平稳波动带）作为"收缩项 + 反弹机制反复触发"的联合长期行为仍是 Working Hypothesis。

## 12. Evidence Classification

- **Model-derived Mathematical Result**：CDF-order 定理（G ≤ T ⟹ active 步严格收缩）；精确单步反弹判据（iff）；suffix 不可能性；旧 argmax 不可能性；前缀 margin 下界 γ ≥ alpha·D_max 与必要条件 M > alpha·D_max·K/(K−1)；global sufficient bound；幅度上界。
- **Verified Numerical Check**：A 轨迹 5344 轮判据验证（FP=FN=0，残差 ≤ 2.0×10⁻¹⁵）；反弹位置预测；tiny-K 扫描（2955 步）；global/magnitude bound 逐轮核验。
- **Model-dependent Observation**：反弹占比 ≈ 0.36–0.37 跨 K 稳定；反弹/下降轮的诊断分解（M、Δ、γ、R 的中位差异）。
- **Working Hypothesis**：反复触发的单步反弹机制 + target-matched 收缩的联合长期行为 = 观测 churn 带（未证明）。

## 13. Limitations / What This Does NOT Establish

1. 长期结果一个都没有证明（§11）；
2. G ≤ T 的"普适必要性"未证明（只有单步 iff 与反例 A）；
3. 判据依赖 selection 的 j\* 与 M——它是给定 controller 下的判据，不含 selection 几何的变化分析；
4. 判据与界都是离散精确陈述；continuum limit 未建立；
5. 反弹占比 ≈ 0.36 的不变性未理论化（为什么是这个量级完全未知）；
6. tiny-K 扫描与 A 轨迹核验是数值检查，不是证明（证明部分独立成立于 §3/§5 的推导）。

## 14. Next Research Gate（交回 GPT + Owner）

1. **长期化的最小一步**：能否把单步判据 + 再注入的统计结构升级为"反弹无穷次/带不变性"的定理（A 不收敛证明）；
2. **G ≤ T 类的刻画**：普适必要性（是否 ∃ 状态使任何 Δ(j₀)>0 的 law 在某步反弹）；
3. **反弹占比 0.36 的来源**：判据在 A 轨迹上的触发率理论；
4. 依据 Prompt §21，本轮在三项交付（序定理、精确判据、A 轨迹验证）完成后停止。

## 15. Repository 资产

- `src/sand_m0/kernel_theory.py` 扩展：`cumulative_mismatch` / `contraction_margin` / `rebound_criterion` / `cumulative_excess_`（pure functions，含定理表述 docstring）；
- `tests/test_m1c3_order_criterion.py`（11 项断言；含 5344 轮 iff 验证与 tiny-K 扫描）；
- 本文档；无新实验目录、无新模拟数据。
