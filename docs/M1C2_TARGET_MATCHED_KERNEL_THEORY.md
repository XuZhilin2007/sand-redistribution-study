# M1C.2 理论文档：Target-Matched Kernel Theory——为什么 target-matched fixed redistribution 没有 churn

> 状态：数学机制审计轮（2026-09-19）。M1C.1 的直接后续；**本轮只做理论推导 + 恒等式数值核验，无新模拟、无新变体、无参数扫描**。全部推导从 repository 当前 canonical 离散更新规则（`src/sand_m0/adaptive.py::run_m1a_deterministic`）出发。
>
> 前置：[M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md](M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md)（U 无 churn 的实证）、[M1B_MICRO_CORRECTION.md](M1B_MICRO_CORRECTION.md)（A 的 churn 实证）、[M1A_LONG_HORIZON.md](M1A_LONG_HORIZON.md)（D_max 与 one-sided KS functional 的恒等）
>
> 核心问题（North Star）：
>
> > **为什么 target-matched fixed redistribution（U）没有 churn，而 near-biased fixed redistribution（A）进入 churn？**
>
> 一句话回答（本轮确立）：**一轮更新的 D-profile 变化可以严格分解为一个"在任何状态下都严格收缩 D_max 的 target-matched 分量"加上唯一的失配项 `M·D_g(j)`（被移除质量 × redistribution law 自己对目标的累计过剩）；U 的失配项恒为零，A 的失配项恒为正、峰值 M/4，与整个 churn 波动带同量级。**

## 1. Question

M1C.1 确立：selection rule 不变时，fixed + deficit-unaware 但 target-matched 的 respray law（U）无 churn（t=6 exact stop，全部 canonical K），而 near-biased law（A）产生 persistent churn（∝K）。本轮把该 simulation-level observation 提升为 mathematical mechanism：

- 主问题（Q1/Q2）：`G = T` 时 `D_max` 是否必然 non-increasing？何时 strict？
- 结构问题（Q3）：为什么 `D_max = 0` 不等价于 bin 级均匀（U 停止态 U_density ≈ 5.4×10⁻³）？
- 分离问题（Q4）：为什么 deficit-fill oracle B 能同时清零 D_max 与 bin 级偏差？
- 失配问题（Q5/Q6）：`G ≠ T` 时一步更新多出了什么项？near-biased 与 uniform 的 export 几何差多少？

## 2. Exact canonical assumptions（全部来自当前代码）

- **空间**：K 个等宽 bin（等面积坐标），质量向量 p ≥ 0，Σp = 1。
- **目标**：均匀 T(j) = j/K（selection rule 与 U 指标编码的目标）。
- **累计分布与过剩**：F(j) = Σ_{i≤j} p_i；D(j) = F(j) − T(j)；D(K) = 0（质量守恒）⟹ **D_max := max_j D(j) ≥ 0 恒成立**。
- **Selection**：j\* = argmax_j D(j)，并列取最小下标（`np.argmax` 约定）；推导不依赖 tie-breaking。
- **Active 条件**：D_max > stop_tolerance（canonical 1e-12）；active ⟹ D_max > 0。
- **Removal**：前缀 1..j\* 每 bin 乘 (1−alpha)；被移除质量 M = alpha·F(j\*)；canonical alpha = 0.25。
- **Respray**：`p_next = p_minus + M·g`，g 为固定 bin 概率（sum 1），G(j) = Σ_{i≤j} g_i。
- **记号警告**：continuum 记号（F(x)、G(x)、T(x)）只是表达捷径；**canonical 模型是离散的，continuum limit 未建立**。D_max 与连续 one-sided KS functional 的离散/连续恒等已在 M1A.1 确立，但本轮一切证明都在离散设定。

## 3. One-step update derivation（引理 1：region identity）

一轮 active sweep 后：

- j ≤ j\*：F′(j) = (1−alpha)F(j) + M·G(j)；
- j > j\*：F′(j) = F(j) − M + M·G(j)（前缀失去的 M 从一切延伸进扫掠区的累计质量中扣除）。

减去 T(j)、代入 F = D + T、用 M = alpha·F(j\*)：

> **引理 1（一步 region identity，精确）**
>
> - Region 1（j ≤ j\*）：`D′(j) = (1−alpha)·D(j) + alpha·[F(j\*)·G(j) − T(j)]`
> - Region 2（j > j\*）：`D′(j) = D(j) − M·(1 − G(j))`

**引理 2（mismatch 分解，精确）**。把 G 与 T 分离（两侧 Region 的失配项**相同**，因为 alpha·F(j\*) = M）：

> - Region 1：`D′(j) = (1−alpha)·D(j) − alpha·T(j)·(1 − F(j\*)) + M·[G(j) − T(j)]`
> - Region 2：`D′(j) = D(j) − M·(1 − T(j)) + M·[G(j) − T(j)]`
>
> **失配项 = M·D_g(j)**，其中 `D_g(j) := G(j) − T(j)` 是 redistribution law **自己对目标的累计过剩**。即：每轮被移除的质量按 law 自身的失配 profile 把累计过剩**再注入**回每个前缀。

**引理 3（边界 export 恒等式，精确）**：`F′(j\*) − F(j\*) = −M·(1 − G(j\*))`——M1B.2 net-export 恒等式对任意固定 respray law 的推广。

数值核验：引理 1 在 A、U 全部 **5368 个 active 轮**（4 个 K）上与实际下一状态 D 逐位对照，最大残差 **5.4×10⁻¹⁵**；引理 2 残差 ≤ 2.5×10⁻¹⁵；引理 3 残差 ≤ 3.4×10⁻¹⁶（`tests/test_m1c2_kernel_theory.py`）。

## 4. Target-matched case 与单调性定理（本轮主结果）

设 G = T（canonical U：g_i = 1/K），alpha ∈ (0,1)，active 轮（D_max > 0）。引理 2 中失配项恒为零：

- Region 1：`D′(j) = (1−alpha)·D(j) − alpha·T(j)·(1 − F(j\*))`。两项符号：第一项 ≤ (1−alpha)·D_max（因 D(j) ≤ D_max 且 D_max ≥ 0）；第二项 ≤ 0（因 T(j) ≥ 0、F(j\*) ≤ 1）。故 **D′(j) ≤ (1−alpha)·D_max < D_max**。
- Region 2：`D′(j) = D(j) − M·(1 − T(j)) ≤ D(j) ≤ D_max`（M > 0、T(j) ≤ 1），且对每个 j < K **严格**小于 D(j)（per-prefix 单调下降）。
- D′(K) = F′(K) − T(K) = 1 − 1 = 0。

> **定理（target-matched 单调性；Model-derived Mathematical Result）**
>
> 在 canonical 离散 M1 模型中，若 respray law 与目标精确匹配（G(j) = T(j) 对一切 j），则每个 active 轮满足 **D_max(t+1) < D_max(t)**（strict），且 D_max(t) ≥ 0 始终。因此 D_max 严格递减，直至触发停止条件（D_max ≤ tol）后冻结。
>
> **假设与适用域**：离散模型；0 ≤ T(j) ≤ 1、T(K) = 1（一般目标即可，均匀目标只是实例）；G = T 精确成立；alpha ∈ (0,1)；M = alpha·F(j\*) > 0（active 轮自动满足）；与 tie-breaking、tol 具体数值无关。
>
> **注意事项（不 overclaim）**：
> 1. **无一致收缩率**：Region 2 的单轮下降量 M(1−T(j)) 在 j\* 接近 K 时可任意小——定理不给出几何收敛率；
> 2. **不证明有限时间停止**：严格递减 + 下有界只保证收敛到某个 L ≥ 0；"t=6 到达 tol"在全部 4 个 canonical K 上是实测（Model-dependent Observation），非定理；
> 3. **不覆盖 U_density**：U 轨迹上 U_density 严格下降是另一个泛函的实测行为，本定理不涉及；
> 4. **不证明收敛到均匀目标**（见 §6——D_max = 0 的状态一般不是均匀态）。

**推论（反事实）**：证明只用了当前状态，未用历史——因此 target-matched 一步映射对**任意** D_max > 0 的状态严格收缩。在 A 的全部 699 个（K=100）状态上验证：把 A 的状态代入 target-matched 公式，D_max 每轮严格下降；而 A 的**实际** D_max 在 36% 的轮次上升。**selection/状态不是 A churn 的来源——失配项是**。

## 5. 失配机制：A 的更新 = 收缩项 + 正的再注入（Q5/Q6）

对 A（g = q_near，T 均匀）：`D_g(j) = G_near(j) − j/K = (j/K)(1 − j/K) ∈ [0, 1/4]`，峰值恰在 j = K/2。

> **A 的每轮更新 = target-matched 严格收缩项 + M·(j/K)(1−j/K) ≥ 0 的系统性正再注入。**

量级（K=100 实测）：每轮峰值再注入 M/4 的中位数 **0.0376**（范围 0.0039–0.0623）；在下一轮 argmax 处实际"被看到"的再注入中位数 **0.0274**（p95 0.0590）。而 A 的整个 D_max 波动带为 p5–p95 = **0.0046–0.0596**——**每轮再注入与整个 churn 带同量级**。U 的失配项恒为零（实测 max |G−T| ≤ 1.7×10⁻¹⁷）。

**Export 几何（Q6）**：净 export = M(1 − G(a))。

| a | A：1 − G_near(a) = (1−a)² | U：1 − T(a) = 1 − a | 比值（A/U）= (1−a) |
|---:|---:|---:|---:|
| 0.50 | 0.25 | 0.50 | 0.50 |
| 0.90 | 0.01 | 0.10 | 0.10 |
| 0.94（U 的消灭扫掠） | 0.0036 | 0.06 | **0.06（约 17× 弱）** |

即：near-biased kernel 下大前缀扫掠的净 export 比 target-matched 弱一个 (1−a) 因子；a → 1 时几乎完全不 export（removed 质量以 2× 近端密度灌回）。这从几何上解释了 M1B.2 的观察（A 的消灭必须用小前缀 F(a) ≤ 0.5）与 M1C.1 的观察（U 的消灭恰是大前缀 a = 0.94）。

**与 churn 的关系（等级：机制解释，非证明）**：A 的失配项恒为正（near-bias 对均匀目标只"多灌"不"少灌"），每轮把与波动带同量级的累计过剩注回系统，恰好抵消 target-matched 分量的收缩——观测到的平稳 churn regime 与"收缩–再注入平衡"一致。但本框架**没有证明 A 不收敛**（无 A 的 D_max 下界定理）；它精确指出了破坏定理假设的那一项。

## 6. D_max = 0 ≠ bin 级均匀（Q3）

`D_max = 0` ⟺ F(j) ≤ j/K 对一切 j ⟺ one-sided KS functional ≤ 0——只约束**运行积分**，不约束单 bin。反例（构造）：`p = (0, 2/K, 0, 2/K, …)`：每个偶前缀恰为均匀份额、奇前缀低于份额 ⟹ D_max = 0 精确，而 U_density = 1（≈3× q_near 初值）。canonical U 停止态是温和版本：D_max ≤ 0（±tol）而 U_density = 5.42×10⁻³（bin 偏差 −0.12u…+0.23u，40 bin 高于均匀）。**controller 的判据泛函（累计过剩）与评价泛函（bin 级 L2）不同，前者的零点集远大于后者。**

## 7. 为什么 B 能同时清零两者（Q4）

B 的停止构型是**饱和构型**：某轮扫掠后无任何 bin 高于目标（E_minus = 0）。此时由 M1C 恒等式 D_def = M 精确成立，deficit-比例回填把每个 deficit bin **恰好**填到 u_i；"所有 bin ≤ u 且总和 = 1"强制成均匀态 ⟹ U = 0 且 D_max = 0（机器精度）。U 的停止集 {D_max ≤ 0} 不约束 bin。**两种 controller 的停止判据是不同泛函的零点：B 的判据（逐 bin deficit）与点态均匀一致，U 的判据（累计过剩）不一致**——这是 10⁻³¹ 与 5.4×10⁻³ 之间 ~28 个数量级差距的结构来源，而非 U"收敛慢"。

## 8. Numerical verification（汇总）

全部由 `tests/test_m1c2_kernel_theory.py`（14 项断言）在 canonical 轨迹上核验：

| 恒等式/定理 | 核验范围 | 最大残差 |
|---|---|---|
| 引理 1 一步 region identity | A、U 全部 5368 active 轮 × 4 K | 5.4×10⁻¹⁵ |
| 引理 2 mismatch 分解（含 M·D_g 精确形） | A 全部 active 轮 × 4 K | 2.5×10⁻¹⁵ |
| 引理 3 边界 export 恒等式 | A、U 全部 active 轮 × 4 K | 3.4×10⁻¹⁶ |
| 定理 strict 单调性（U 实轨迹） | U 全部 active 轮 × 4 K | 严格 <（逐轮） |
| Region 1 界 D′ ≤ (1−alpha)D_max、Region 2 per-prefix 单调 | U 全部 active 轮 × 4 K | ≤ 1e−12（j=K 处 D(K)=0±ulp 的实现态浮点残差） |
| 反事实（定理对任意状态成立） | A 全部 699 轮（K=100）+ 4 K | 严格 <（逐轮） |
| D_max=0 非均匀构造态 | (0, 2/K, 0, 2/K, …) | D_max=0 精确、U_density=1 精确 |
| U 停止态分离 | canonical U t=6 状态 | D_max ≤ 1e−12、U_density ∈ (1e−3, 1e−2) |
| export 几何 (1−a) 因子 | a ∈ {0.5, 0.9, 0.94} | 精确 |

## 9. Evidence Classification

- **Model-derived Mathematical Result**：引理 1/2/3（一步 region identity、mismatch 分解 M·D_g、边界 export 恒等式）与定理（G = T ⟹ D_max 严格单调下降直至停止）——从 canonical 更新规则严格推导，全部经数值核验。定理对一般目标 T（0 ≤ T ≤ 1、T(K)=1）成立，均匀目标为实例。
- **Verified Numerical Check**：§8 的全部数值核验（在既有 canonical 轨迹上，非新实验）。
- **Model-dependent Observation**：U 在 4 个 K 上 t=6 有限时间停止；U_density 严格下降；A 的再注入量级（0.027–0.038/轮）与 churn 带同量级；A 消灭几何为小前缀、U 为大前缀。
- **Working Hypothesis**：A 的 churn 是"严格收缩项 + 恒正失配再注入"的平衡态——解释了机制但**不是** A 不收敛的证明。

## 10. Implications

1. M1C.1 的实证结论获得机制基础：churn 的关键承载结构是 respray law 的失配项 M·D_g；target-matched ⟹该项消失 ⟹收缩定理成立；
2. "deficit-unawareness 本身导致 churn"被理论排除：U 完全 deficit-unaware 却满足收缩定理；
3. churn 的理论化路径清晰化：要么证 A 的下界/不变集定理（失配项阻止 D_max → 0），要么证明 A 的 D_max 收敛到正极限——两者都超出本轮；
4. U 的有限时间停止（t=6）与 B 的饱和构型成为下一步分析对象（§12）。

## 11. Limitations / What This Does NOT Establish

1. 无收敛率、无有限时间停止定理（t=6 是观察）；D_max 的收敛极限 L 是否为 0 未证明（仅 U 实测到达 tol）；
2. A 的非收敛（churn 持久性）未证明——只有失配项与波动带同量级的定量事实；
3. U_density 单调性未理论化；停止态锯齿结构（U_density ≈ 5.4×10⁻³）未做到达/稳定性分析；
4. 离散设定；continuum limit 未建立；
5. 只考虑均匀目标下的 target-matched（一般 T 的定理是形式推广，未对非均匀 T 数值实验——本轮禁止新模拟）；
6. B 侧分析沿用 M1C 已有恒等式，本轮未新增 B 的理论。

## 12. Next Research Gate（交回 GPT + Owner）

1. **U 的停止时间**：t=6 是否有界（K 无关）？能否从收缩定理 + 初值结构推导停止轮数的显式界？
2. **A 的持久性**：失配项 M·D_g ≥ 0 能否升级为 D_max 的下界/不变集定理（即证明 A 的 churn 是必然而非观察）？
3. **一般 target-matched law**：非均匀目标 + 匹配 law 的行为（定理已覆盖，模拟留待决策）；
4. **三种停止态的泛函几何**：{D_max = 0} 的 U_density 上确界、B 的饱和构型到达率；
5. 依据 Prompt §13，本轮在理论 + 最小核验完成后停止。

## 13. Repository 资产

- 数学 helper：`src/sand_m0/kernel_theory.py`（`one_step_excess` / `target_matched_excess` / `mismatch_reinjection` / `prefix_cdf`，含定理表述的 docstring）；
- 测试：`tests/test_m1c2_kernel_theory.py`（14 项）；
- 文档：本文件；无新实验目录、无新模拟数据。
