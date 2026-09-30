# Sand Redistribution Study — Stage 2 Closure
## Feedback Mismatch Mechanism

> 日期：2026-09-20
>
> 性质：**Stage 2 closure checkpoint**。目标是让完全不了解当前聊天的新 GPT/Zcode 窗口仅凭 repository reality 就能继续：Stage 1 与 Stage 2 研究了什么、哪些是定理、哪些只是 canonical observation、哪些解释已被推翻、为什么下一阶段必须转向 generality。
>
> 当前 reality（写作时核验）：commit `361acd4`，branch `main`，working tree clean，测试 **288 项运行时断言**全部通过（`python tests/count_tests.py`）。
>
> 前置 checkpoint：Stage 1 闭合记录（2026-09-19，M0 → M1B.2；保留在私人研究档案）。本文取代其"下一步"部分，不改写其历史内容。

---

## 1. Research origin（现实来源）

一切源于 Owner 的真实操作（Real-world Observation）：

```text
观察当前沙面
→ 扫回明显过密 / 过近区域
→ 重新泼出
→ 再观察
→ 适时停止
```

**核心原始问题**：

> 在带反馈的随机再分配过程中，为什么有限时间内可能先改善、随后反复纠偏，并最终在某种停止条件下结束？

Stage 2 的全部理论都是对这一现象的最小机制追问，不淹没现实来源；也尚未与现实数据接触（real throw distribution 未测量）。

## 2. Research North Star（原则，未变）

1. **minimal mechanism first**：每轮只加一个机制、只回答一个问题；
2. theory 与 reality-inspired validation 严格分离；
3. falsification 正式保留（见 §10）；
4. **repository reality > chat history**；
5. 不为 novelty 强行复杂化；
6. **不把 canonical 特例结果伪装成普适理论**（Stage 2 最重要的纪律，见 §5/§6 的 scope 划分）。

## 3. 当前模型与记号（Stage 2 通用对象）

一维等面积坐标 x∈[0,1]，K 个等宽 bin，质量 p≥0、Σp=1。**目标分布 T**（canonical：均匀，T(j)=j/K）；**redistribution law G**（fixed、外生，canonical q_near：G(x)=2x−x²）；**cumulative excess** D(j)=F(j)−T(j)、D_max=max_j D(j)（= one-sided KS functional，M1A.1 恒等）；**selection**：a=j\*/K 为 D 的 argmax；**removal**：扫掠前缀每 bin 移除 α 比例（canonical α=0.25），M=α·F(a)；**redistribution**：移除质量按 fixed G 重泼。全部结论限于 **deterministic expected-mass path、finite K、exact/finite-deadband stopping**；MC 路径、continuum limit、随机性均在 scope 之外（M1C 系列 scope 声明见各 config）。

## 4. Canonical research ladder（紧凑版）

### M0 — fixed selection + fixed redistribution（Stage 1，理论闭合）

- finite optimum / overshoot（q_near 最佳 t≈4.46 → t=4/5 平台）；
- 单状态变量精确解 z=(1−α/4)^t；alpha 只是时间尺度参数；
- **一般 fixed-q 内部最优判据** `Q_f·S_n > Q·S_f`（Model-derived Mathematical Result，幂族 6/6 预测命中）；
- M0 支线（fixed-q/fixed-zone）理论基本闭合。

### M1A / M1B — state-dependent cumulative-excess feedback（Stage 1）

observe → decide → act：a_t = argmax D_t(j)；D_max ≤ tol 停止。

- exact stopping 存在但 stopping time ∝ K（324/699/1445/2876）；
- finite deadband（eps ∈ {0.005..0.04}）使停止时间跨 K 稳定（≤1.34×）且 improvement retention 91–99.9%；
- late stage **不是** micro-correction（动作恒定 ~15%/轮，M1B.2 证伪）；
- canonical 系统进入 substantial-redistribution **churn** regime（D_max 在 0.003–0.08 带内震荡数百至数千轮）；
- exact stop 来自低频 **annihilable stopping configuration**（一次扫掠清零全部正累计过剩；F(a)≤0.5 必要条件，4/4）。

## 5. Stage 2 core story — M1C 系列（问题 → 结果 → 证据等级）

### M1C — Correction-kernel ablation（oracle 对照）

**问题**：churn 是否来自 selection 还是 redistribution？

- A（canonical selection + fixed near-biased q）：churn（699 轮 @K=100）；
- B（同 selection + deficit-aware oracle）：全部 canonical K **t=3 达机器精度均匀**并 exact stop；总搬动质量 0.32 vs A 的 103。

**结论**：cumulative-prefix selection 本身并不必然导致 churn；redistribution mechanism 是关键变量。
**等级**：Verified Simulation Result / Model-dependent Observation（baseline reproduction 112/112）。

### M1C.1 — Fixed-uniform kernel control（拆开两个因素）

**问题**：M1C 同时改变了 adaptive+target-aware+deficit-aware；churn 究竟来自 near-bias 还是 deficit-unawareness？

- U（同 selection、同 q_near 初始状态、**fixed uniform** redistribution——仍 fixed/state-unaware/deficit-unaware）：全部 canonical K **t=6 exact stop**，D_max/U_density 严格单调下降，无 churn；
- U 停止态 D_max=0 但 **U_density≈5.4×10⁻³ 非 binwise 均匀**（锯齿构型）——controller 目标（累计过剩）清零 ≠ 逐 bin 均匀。

**结论**：deficit-unawareness 本身不是 canonical churn 的充分解释（F5，见 §10）；**near-bias / target mismatch** 成为主要嫌疑。
**等级**：Verified Simulation Result / Model-dependent Observation（baseline integrity 132/132）。

### M1C.2–M1C.6 — 从对照到解析理论（见 §6/§7/§8/§9）

- M1C.2：一步分解恒等式 + target-matched 严格单调定理；
- M1C.3：G≤T 序定理 + 精确单步 rebound 判据（5344 轮 FP=FN=0）；
- M1C.4：churn regime 审计（跨 K 稳定带、0.36=P(R>1)、C-C-R 节奏、边界门控 observation）；
- M1C.5：解析反弹门 a_crit(D)（1958 rebound 轮 FN=0）；
- M1C.6：boundary return + peak competition（反馈回路闭合；证伪简单两峰交换）。

## 6. General theory established（与 canonical observations 严格分开）

以下为 **Model-derived Mathematical Result**：从 canonical 离散更新规则严格推导、对声明的假设类成立（非 real-world、非 continuum、非全模型族）。

### 6.1 一步分解恒等式（Stage 2 最核心的一般结构）

```text
D_new(j) = D_TM_new(j) + M·[G(j) − T(j)]
```

- T = target CDF，G = fixed redistribution law CDF，Δ(j)=G(j)−T(j)；
- D_TM_new = 同状态、同 selection、同 removal、**G=T** 时的 counterfactual；
- **mismatch term = M·Δ**：每轮被移除质量按 law 自身对目标的累计过剩 profile 重新注入（Region 1/2 同式）；
- 推论（M1C.5 一般化 net-export 恒等式）：prefix change = −M(1−G(a))，对任意 fixed law 成立。

（M1C.2；核验残差 ≤5.4×10⁻¹⁵，5368 轮 × 4 K。）

### 6.2 Target-matched monotonicity theorem

**假设**：离散模型；G(j)=T(j) ∀j；α∈(0,1)；active 步（D_max>0）；T(K)=1、0≤T≤1、T(j)<1（j<K，strictness 所需）。

```text
G = T 且 active step ⟹ D_max(t+1) < D_max(t)   （strict）
```

**边界（不过度声称）**：只针对 one-sided cumulative discrepancy；D_max=0 不等价于 binwise 等于 target（允许 U_density 高达 1 的锯齿态，构造 (0,2/K,0,2/K,…) 已验证）；不自动证明 finite stopping time（U 的 t=6 是观察）。无一致收缩率。

（M1C.2；对一般 T 成立，均匀目标为实例。）

### 6.3 CDF-order theorem（sufficient stability class）

**假设**：G、T 为离散 CDF（值域 [0,1]、G(K)=T(K)=1）、T(j)<1（j<K）、α∈(0,1)、active 步。

```text
G(j) ≤ T(j) ∀j（cumulatively target-dominated redistribution）
⟹ 每个 active step 严格收缩 D_max
```

**明确**：这是 **sufficient stability class，不是 iff characterization**（G>T somewhere 的 law 仍可能逐步收缩——反例已构造；普适必要性未证明）。

（M1C.3。）

### 6.4 精确单步 rebound criterion

```text
D_max(t+1) > D_max(t)  ⟺  ∃j: M·Δ(j) > γ(j)
γ(j) = D_max(t) − D_TM_new(j)   （pointwise contraction margin，active 步逐点>0）
```

直觉：**mismatch reinjection 超过该位置 target-matched update 原本能获得的纠偏幅度**。附带结构定理：suffix（j>j\*）永不反弹（任意 law）；旧 argmax 处不可能反弹；前缀 margin 下界 γ≥α·D_max ⟹ 反弹必要条件 M > α·D_max·K/(K−1)（均匀目标）。

**canonical A 轨迹验证**：5344 active 轮、1958 rebound 轮、**FP=0、FN=0**（残差 ≤2.0×10⁻¹⁵）。

**明确**：one-step rebound theorem ≠ long-run churn theorem。

（M1C.3。）

## 7. Canonical q_near specialization（canonical specialization, not general theorem）

以下只对 canonical T(x)=x、G_near(x)=2x−x² 成立，正式划分为 canonical 特例：

```text
Δ(x) = x(1−x)，max Δ = 1/4（j=K/2）
正 mismatch reinjection 每轮至多 M/4（实测中位 0.027–0.038/轮，与整个 churn 波动带 0.005–0.060 同量级）
export geometry：1−G_near(a) = (1−a)²  vs  target-matched 1−T(a) = 1−a
```

含义：near-biased redistribution 对 large-prefix sweep 的净 export correction 弱一个 (1−a) 因子（a=0.94 时约 17× 弱）——几何上解释 A 的小前缀消灭构型（F(a)≤0.5，M1B.2）与 U 的大前缀消灭（a=0.94，M1C.1）。

## 8. Analytic rebound gate（M1C.5 定理，canonical specialization）

```text
rebound ⟹ a > a_crit(D_max)
a_crit(D) = [(1−D) + √(D(D+2))] / 2
```

推导链：Region-1 化简 `MΔ−γ = (1−α)D − D_max + α·h`（h(x)=x[F_a(2−x)−1]，残差 ≤4×10⁻¹⁶）→ state-reduced 上界 `≤ α(max_j h − D_max)` → h 最大化（F_a≤1/2 ⟹ 精确不可能分支；F_a>1/2 ⟹ 内部 max；x\*≤a 条件 (a+D)(1−a)≤1/2 在 canonical 全程成立，D_max≤0.25<√2−1）→ 闭式门槛。

**必须明确**：

- **必要条件，非充分**（gate-open contraction 轮每 K 21–177 个；P 因子 P(x_r)=D(x_r)/D_max 分布与 rebound 轮分离——profile 形状决定兑现）；
- continuum-style closed-form envelope 与离散格点条件差两种 rounding regime（一般 O(F_a/4K²)、x\*<1/K 相位 regime O(1/K)，实测 ≤1.4×10⁻⁴）；离散精确条件以格点 max_j h(j/K) 为准；
- canonical 轨迹全部 1958 rebound 轮 **FN=0**；
- a_crit churn-window 分位数 p50 = **0.593–0.604**（跨 K 稳定）、p5–p95 = 0.539–0.644——**解析门槛 tracks M1C.4 观察门 ≈0.58**。

一句话故事：**boundary 决定 rebound 是否进入可能区，full D-profile 决定可能性是否真正兑现。**

## 9. Canonical churn dynamics（全部为 Model-dependent Observation / Verified Numerical Analysis）

不得升级为 general theorem。

### 9.1 Cross-K stable operating band（M1C.4）

K={50,100,200,400} 的 churn window（onset=首个 D_max<0.05 的轮，=23，4 K 相同；M1B.2 divider）：

- **R_t = max_j M·max(Δ,0)/γ 分布 collapse**：p5/p50/p95 = 0.40–0.43 / 0.689–0.708 / 3.89–4.49；
- normalized boundary、moved mass band、|Δa|、switch freq（0.936–0.944）全部近似稳定；
- **rebound fraction ≈0.362–0.369 = P(R_t>1)**——0.36 不是独立 universal constant，而是稳定 R 分布穿阈值的 empirical frequency。

### 9.2 Boundary-gated switching（M1C.6）

```text
gate-open large boundary → rebound → large boundary switch
→ contraction（gap-closing share = 1.000）→ gate return → next rebound
```

- post-rebound gate exit：P(B_post<0) = 0.93–0.97；
- **gate never-return = 0 in 1958 episodes**；typical gate return = 1–2 contraction rounds（median tau_return=2）；
- next rebound median ≈3 rounds（max 5–7）；
- boundary map **两分支**（近端 a∈[0.05,0.65]→远端 a∈[0.6,1.0] 与反向；几乎无点在对角线）；
- near/far boundary large hopping（|Δa| 中位 0.39–0.40）；
- B 逐轮交替震荡跨 K collapse（快弹跳 + a_crit 慢降的两时间尺度）。

**措辞边界**：只说"当前 canonical trajectories 中全部 observed episodes return"；**不说**"已证明永远 return"。

### 9.3 Terminal behavior

terminal episode 与普通 episode 结构同构（len 2–5 vs median 3；gate 在 4 K 中 2 个曾重开 1 轮未兑现）；exact stopping = 普通 churn episode 恰好进入 annihilable geometry（衔接 Stage 1 M1B.2），无独立 approach-to-stop 动力学。

## 10. Peak/profile mechanism——被证伪的简化图像（M1C.6）

**原候选**：canonical switching 主要是两个固定峰交换主导（two-peak exchange）。

**数据否定其作为主机制**：既有 secondary peak 接管（P1）占比随 K **7.8% → 0.85%** 衰减；reinjection 重塑/新生峰（P2_new）**31% → 86%** 递增；第三峰接管（P_multi）56% → 12%。

**当前解释**：mismatch reinjection M·Δ(x) 对**整个** cumulative-excess profile 重新塑形；boundary switching 不是简单两个固定峰轮换。（已列入 §10 falsified list F7。）

**保留的精确内容**：逐对峰高更新恒等式 `ΔH = ΔH_TM + M[Δ(j1)−Δ(j2)]`（M1C.2 推论，残差 ≤2.1×10⁻¹⁵）；contraction 轮 gap-closing share = 1.000（active 峰受 (1−α)+export 双重压制 −0.036/轮，competing 峰仅 −0.021/轮）。

## 11. Stopping mechanism（整合 Stage 1 + M1C.4/6）

finite-K canonical A **已 finite-time exact stop（324/699/1445/2876）——一定不要描述成"长期不收敛"**。

当前事实：

- churn 是 **long pre-stopping transient**；
- terminal episode 无独立动力学阶段（与普通 episode 同构：len 2–5 vs median 3；gate 在 4 K 中 2 个曾重开 1 轮未兑现 rebound）；
- stopping 发生于普通 churn episode 中 trajectory 进入 annihilable geometry（D_max 谷 + 小前缀消灭构型，M1B.2）；
- finite deadband 不需要等待这种 exact geometry（在普通波动谷停手，retention 91–99.9%）。

严格避免（除非未来有独立证据）：infinite churn、non-convergence、invariant distribution、attractor、chaos、ergodicity。

## 12. Explicitly falsified / superseded interpretations（禁止复活）

| # | 已废弃解释 | 推翻证据 | 现行解释 |
|---|---|---|---|
| F1 | "M1A 不会停止"（T=100 窗口结论） | M1A.1：T=5000 下 exact stop 于 324/699/1445/2876 | 会停止；等待时间 ∝K；stop = 稀有 annihilable configuration |
| F2 | "large boundary jumps 主要是 grid plateau artifact" | M1A.1：649 大跳全为真实 peak switching（96% near-tie 为相邻格点） | 大跳是真实 profile 切换；M1C.4/6 进一步给出两分支 boundary map 与门控关联 |
| F3 | "finite epsilon near-optimality 很差" | M1B.1：small-denominator amplification；retention 91–99.9% | retention 语义下存在宽 near-optimal region |
| F4 | "late-stage moved mass shrinks as O(1/K)（micro-correction regime）" | M1B.2：moved/L1/net_export 跨 K 恒定（CV≤3.5%） | 平稳 churn 波动 regime + 稀有 annihilable stopping configuration |
| F5 | "deficit-unawareness 本身是 churn 主要原因" | M1C.1：U 仍 fixed+deficit-unaware 却无 churn（t=6） | near-bias / target mismatch 是关键驱动；deficit-awareness 买 bin 级精度（5.4×10⁻³ → 10⁻³¹） |
| F6 | "churn 主要由 cumulative-prefix selection 本身造成" | M1C：同 selection 换 kernel 后 churn 消失（B t=3）；M1C.1：U 无 churn | redistribution law 与 target 的 mismatch 是关键变量 |
| F7 | "canonical switching 主要是两个 persistent peaks 轮换主导" | M1C.6：P1 7.8%→0.85% 随 K 衰减，P2_new 31%→86% | mismatch reinjection M·Δ(x) 重塑整个 profile；switching 非两峰轮换 |
| F8 | "A 是不收敛系统 / 应证明 infinitely many rebounds" | 方向性错误：finite-K canonical A 已 exact stop（324/699/1445/2876） | churn 是 long pre-stopping transient；研究对象是 pre-stopping 动力学，不是 non-convergence |

另有 Stage 1 口径勘误三条（D_max 数值、测试计数 ×2）与 Stage 2 实现勘误一条（M1C.6 的 0/1-based 索引错位，pair-identity 残差 5.6e-3 → 2.1e-15，已留痕 RESEARCH_LOG）。

## 13. Evidence hierarchy snapshot

| Result | Scope | Evidence class |
|---|---|---|
| M0 fixed-q interior optimum condition（`Q_f·S_n > Q·S_f`） | general fixed-q M0 | Model-derived Mathematical Result |
| 一步分解 `D_new = D_TM_new + M·Δ` | general M1 update class | Model-derived Mathematical Result |
| target-matched strict D_max decrease | general M1 update class（G=T） | Model-derived Mathematical Result |
| G≤T order theorem（cumulatively target-dominated） | general M1 update class | Model-derived Mathematical Result |
| exact rebound iff criterion（M·Δ > γ） | general fixed-G/T update | Model-derived Mathematical Result |
| 逐对峰高更新恒等式 `ΔH = ΔH_TM + M·Δdiff` | general M1 update class | Model-derived Mathematical Result |
| a_crit gate `a > [(1−D)+√(D(D+2))]/2` | canonical q_near specialization | Model-derived Mathematical Result |
| rebound fraction ≈0.36 = P(R>1) | canonical A, 4 K | Model-dependent Observation |
| gate return never failed in 1958 episodes | canonical A | Verified Numerical Analysis / Model-dependent Observation |
| two-branch boundary map；C-C-R rhythm | canonical A | Model-dependent Observation |
| t_stop ≈ K；stopping = annihilable configuration | canonical finite-K observations | Model-dependent Observation |
| deadband retention 91–99.9% | canonical M1B family | Verified Simulation Result |
| B/U t=3/t=6 机器精度收敛 | canonical A vs B/U 对照 | Verified Simulation Result |

**scope 纪律**：general M1 update class 的定理不自动适用于其他 selection 规则、其他目标族、continuum 或随机版本；canonical specialization 不外推。

## 14. Literature positioning（Stage 2 前已核验的定位；本轮不新增检索）

Neighboring fields（Stage 1 literature positioning 冻结的线索，均未做 novelty 声明）：discrepancy reduction / point movement / deletion；online thinning；rank-driven selection-resampling processes；nonlinear Markov / selection–mutation particle systems；density-feedback swarm control；deadband / threshold control（已核实连接：Luceño 2003 dead-band adjustment schemes, Handbook of Statistics 22, Crossref-verified——仅支持结构与 M1B 的相似性）。

**Honest positioning**：大框架已有；state-dependent excess selection + stochastic redistribution 与 density-feedback control 工作高度近邻。**本项目不能 claim "invented a new process family"，novelty 未建立。**

当前项目更准确的潜在价值：

> **当 feedback selection 已知哪里存在 cumulative excess、但 correction 被约束为一个 fixed exogenous redistribution law 时，target mismatch（Δ = G − T）如何决定 finite-time rebound、switching 与 churn-like pre-stopping transients**——即 mismatch-structured 的 stability/instability 理论（M1C.2/3/5 的定理链）+ 其 canonical 动力学实现。

## 15. What remains unresolved

### Generality

- current long-churn dynamics 是否只属于 q_near；
- positive cumulative mismatch 的幅度/形状到什么程度产生 sustained pre-stopping churn；
- sign-changing mismatch（far / target-dominated 之外）的 dynamics；
- target-dominated laws 的 dynamics 定理已有、轨迹行为未系统观察。

### Theory

- long pre-stopping regime 无 general theorem（只有单步判据 + 门控）；
- gate-return 在 canonical trajectories 稳定（never=0），但无 general return theorem；tau_return≈2 无解析来源；
- rebound fraction ≈0.36 无解析来源；
- t_stop ~ K 无解析证明；
- continuous-space limit 未解决；P(x_r) profile factor 未解析化。

### Reality connection

- real throw q 未测量；human perception / epsilon 未校准（当前 eps 无现实数值含义）；
- 2D model 未建立；local patches / noisy perception 未建；finite-particle feedback 未系统研究。

### Novelty

> 未建立（见 §14）。

## 16. Stage 2 closure statement

> **Stage 2 — Feedback Mismatch Mechanism is closed for the current canonical model.**

含义：

- 当前 q_near canonical dynamics 不再继续无限细抠；
- 不继续证明 tau_return≈2、不继续拟合 boundary map、不继续构造越来越复杂的 canonical reduced model；
- 除非未来 generality study 暴露必须回头处理的问题。

这不是项目结束，而是：**canonical mechanism phase 已够完整**——机理链（mismatch → gate → rebound → reshaping → return）已经从对照实验走到解析定理再到逐轮验证，继续在同一个特例上打磨的边际收益低于转向 generality。

## 17. Next stage（冻结方向）

# Stage 3 / M2 — Redistribution-Mismatch Generality Gate

核心研究问题：

> **Are the rebound / switching / churn mechanisms specific to canonical q_near, or do they arise systematically across a broader class of fixed redistribution laws with positive cumulative mismatch relative to the target?**

中文：canonical q_near 的 churn 是一个特殊函数的偶然现象，还是 fixed redistribution law 相对 target 存在 positive cumulative mismatch 时的一类更普遍动力学？

现有理论（M1C.2/3 的分解、序定理、单步判据）直接提供 Stage 3 的预测工具——每个新 law 的行为可在运行前由 Δ=G−T 的 sign/shape/amplitude 部分预测。

## 18. Stage 3 principles（提前冻结）

**Not a blind parameter sweep**：不做 20 个 q、大规模网格、arbitrary powers、tuning until desired result。

**Theory-guided kernel selection**：只选少量 structural controls，例如：

1. target-matched `G=T`（已有：M1C.1 U）；
2. target-dominated `G<T`（定理已保证稳定；轨迹行为未观察）；
3. weak positive mismatch；
4. canonical q_near（已有）；
5. stronger / differently shaped positive mismatch；
6. optional sign-changing mismatch。

**具体函数在 Stage 3 开启时再冻结，本 closure task 不决定。**

## 19. Stage 3 evaluation questions（每个 law 的检查单）

- sign / shape / amplitude of Δ = G − T；
- whether rebound is possible（a_crit gate 的适用性）；
- rebound frequency（P(R>1)）；
- whether sustained churn-like transient appears；
- whether boundary gate exists；
- whether current one-step theory predicts behavior（分解/判据的逐轮验证）；
- exact stopping behavior；
- dependence on K；
- whether the canonical q_near story survives qualitatively。

**closure task 不运行这些实验。**

## 20. Repository state（Stage 2 closure 时）

```text
docs/        RESEARCH_CHECKPOINT_2026-09-19（Stage 1，冻结）、本文件（Stage 2 closure）、
             M1C_CORRECTION_KERNEL_ABLATION / M1C1_FIXED_UNIFORM_KERNEL_ABLATION /
             M1C2_TARGET_MATCHED_KERNEL_THEORY / M1C3_REDISTRIBUTION_LAW_ORDER /
             M1C4_CHURN_REGIME_AUDIT / M1C5_BOUNDARY_GATE_THEORY / M1C6_BOUNDARY_RETURN、
             M0/M1A/M1B 系列、MODEL_ASSUMPTIONS（§11/§12 kernel 登记）、
             RESEARCH_LOG（append-only）、NEXT_STEPS（顶部指向 Stage 3）
src/         sand_m0：model / simulate / mechanism / adaptive / diagnostics / plotting /
             kernel_theory（M1C.2–6 数学 helpers，全部 pure functions）
tests/       15 套，288 项运行时断言（count_tests.py 为唯一口径）
experiments/ m0_baseline … m1c_kernel_ablation / m1c1_uniform_kernel_ablation /
             m1c2…m1c6（M1C.2–M1C.6 为 theory/audit 轮：m1c2 与 m1c5 为 verification
             脚本目录，m1c4/m1c6 为 audit 目录）
archive/     2026-09-13 原始 checkpoint（只读）
```

## 21. New-window bootstrap read order

1. 根 [README](../README.md)（当前阶段、证据标签、禁止误读）；
2. 本文件（Stage 2 closure：研究逻辑、定理/observation 分层、falsified 清单、Stage 3 冻结）；
3. 当时的私人 `NEXT_STEPS.md`（Stage 3 / M2 优先；未收入公开快照）；
4. [MODEL_ASSUMPTIONS.md](MODEL_ASSUMPTIONS.md)（登记的假设）；
5. 需要细节时再回读本快照中的 M1C 系列文档；追加式 `RESEARCH_LOG.md` 留在私人研究档案；
6. Stage 1 checkpoint 仅作历史冻结节点回读。

原则：不需要读取任何旧聊天；被证伪解释以 §12 为准；canonical observation 不得当 theorem 引用。

## 22. Stage 2 交付清单（资产索引）

| 轮 | 交付 | 类型 |
|---|---|---|
| M1C | kernel ablation A/B；baseline reproduction 112/112；deficit-fill kernel | simulation |
| M1C.1 | A/U/B 三方对照；integrity 132/132；U t=6 无 churn | simulation |
| M1C.2 | 分解恒等式 + 单调定理；kernel_theory.py；残差 ≤5.4e-15 | theory + verification |
| M1C.3 | 序定理 + rebound iff 判据；5344 轮 FP=FN=0 | theory + verification |
| M1C.4 | churn regime 审计（R 带、a-门控、C-C-R、terminal 温和） | dynamics audit |
| M1C.5 | a_crit 闭式门；1958 rebound 轮 FN=0 | theory + verification |
| M1C.6 | episode 审计（never-return=0、gap-closing=1.0、P1/P2 分类、两分支 map） | dynamics audit |
