# Sand Redistribution Study — Stage 3 Closure
## Redistribution-Mismatch Generality（M2）

> 日期：2026-09-20
>
> 性质：**Stage 3 closure checkpoint**。目标是让完全不了解旧聊天的新 GPT / researcher / external mathematical reviewer 仅凭 repository reality 就能准确恢复：Stage 3 为什么启动、Gate 1 / Gate 2 实际发现了什么、哪些是 theorem、哪些只是 numerical / model-dependent observation、哪些 generality hypothesis 被否定、canonical q_near 到底特殊在哪里、Stage 2 theory 泛化到了什么程度、以及为什么现在应进入 External Mathematical Review 而不是继续自行加 kernel。
>
> 当前 reality（写作时核验）：commit `77b6f7f`，branch `main`，working tree clean，测试 **330 项运行时断言**全部通过（`python tests/count_tests.py`）。
>
> 前置 checkpoint：[RESEARCH_CHECKPOINT_2026-09-20.md](RESEARCH_CHECKPOINT_2026-09-20.md)（Stage 2 closure，历史冻结节点，不覆盖）。本文不修改任何科学代码、frozen kernels 或 raw experiment outputs（本轮 documentation-only）。

---

## 1. 正式路线与执行记录

```text
Stage 2 Closure (a544eee)
↓
Stage 3 / M2 Gate 1        feat d418ae2 / docs 7cf0fa3   （PASS）
↓
Gate 1 Interpretation Audit feat 795c559 / docs 7629dc4
↓
Stage 3 / M2 Gate 2        feat af4ece5 / docs 77b6f7f   （PASS）
↓
Stage 3 Closure            ← 本轮（documentation-only）
↓
External Mathematical Review Gate（下一步，未启动）
```

三个执行轮均为 Owner + GPT 冻结 prompt 的严格实现：无 post-hoc tuning、无中途增删 kernel、无 horizon 延长。关键资产：

| 轮 | 实验/分析 | 文档 |
|---|---|---|
| Gate 1 | `experiments/m2_gate1_generality/`（K0/K-/KW/KC/KS；5510 active 轮） | [M2_GATE1_REPORT.md](../experiments/m2_gate1_generality/M2_GATE1_REPORT.md)、[M2_GATE1_GENERALITY_GATE.md](M2_GATE1_GENERALITY_GATE.md) |
| Audit | `run_interpretation_audit.py`（5483 轮再生核验零失配） | [M2_GATE1_INTERPRETATION_AUDIT.md](../experiments/m2_gate1_generality/M2_GATE1_INTERPRETATION_AUDIT.md) |
| Gate 2 | `experiments/m2_gate2_sign_placement/`（XA/XO；7524 active 轮） | [M2_GATE2_REPORT.md](../experiments/m2_gate2_sign_placement/M2_GATE2_REPORT.md)、[M2_GATE2_SIGN_PLACEMENT.md](M2_GATE2_SIGN_PLACEMENT.md) |

## 2. Stage 3 core question（原始冻结问题）

> **canonical q_near 的 rebound / switching / churn 是否只是一个特殊 redistribution function 的行为，还是 fixed redistribution law 相对于 target 的 cumulative mismatch 所产生的一类更一般动力学？**

Stage 3 采用 **theory-guided minimal controls**（sign / amplitude / shape / sign-placement 四个结构维度、每维度一至两个 kernel），不是 parameter sweep。

## 3. Preserved Stage 2 general theory（未被 Stage 3 推翻；scope 不变）

以下仍是 **Model-derived Mathematical Result**（general M1 update class；非 real-world、非 continuum、非全模型族）：

1. **一步分解恒等式**：`D_new(j) = D_TM_new(j) + M·[G(j) − T(j)]`，`Δ(j) := G(j) − T(j)`（M1C.2）。
2. **Target-matched 严格收缩**：`G = T` ⟹ active step 严格下降 D_max（M1C.2）。
3. **CDF-order sufficient stability class**：`G ≤ T` ⟹ active step 严格下降 D_max（M1C.3）——**sufficient，不是 iff**。
4. **精确单步 rebound 判据**：`rebound ⟺ ∃j: M·Δ(j) > γ(j)`，`γ(j) = D_max − D_TM_new(j)`（M1C.3）。
5. canonical `a_crit` 闭式门（M1C.5）**仍只属于 canonical q_near specialization**，未在 Stage 3 泛化或使用。

Stage 3 在 7 个新/对照 kernels（K0、K-、KW、KS、XA、XO + canonical KC）上逐轮 numerical verification：分解残差 max ≤ 7.6×10⁻¹⁵（Gate 1）/ 5.3×10⁻¹⁵（Gate 2），non-boundary 判据 FP=FN=0——**实现与理论持续一致，无一例非边界违例**。

## 4. Gate 1 synthesis（sign / amplitude / shape）

| kernel | Δ(x) | 结果（全部 K∈{50,100,200,400}，除 K- 仅 K=100） | 结论 |
|---|---|---|---|
| K0 | 0 | t=6 stop、零反弹 | target-matched reference（复现 M1C.1 U；integrity 132/132） |
| K- | −x(1−x) | t=3 stop、零反弹、严格单调 | 负控制与 G≤T 序定理完全一致（5/5 验收） |
| KW | (1/4)x(1−x) | **t=8 stop、零反弹**（max R_t 0.68–0.71 不穿 1） | **W0**：positive mismatch existence alone is not sufficient to generate rebound along these trajectories；amplitude materially affects dynamics（不可写 general necessity） |
| KC（anchor） | x(1−x) | t=324/699/1445/2876、N_R=1958、f_R≈0.36（复现 Stage 2） | canonical 长瞬态仅为 anchor |
| KS | flat-top（同 max=1/4、同 ∫=1/6、同对称/端点） | **t=26–29、5–10 次 rebound、短促复发** | **S1**：one-step 判据幸存于形状改变，长复发动力学 shape-sensitive；**peak mismatch 与 integrated mismatch 不足以决定 canonical-like long churn**；不得说 "flat top itself causes the difference"（正则性混淆） |

Gate 1 理论验证：5510 active 轮，判据 TP=1983/TN=3527/**FP=FN=0**，numerical-boundary 0 例。

## 5. Gate 1 Interpretation Audit synthesis（witness geometry）

- **KC canonical witnesses**（1958 rebound 轮）：稳定 interior 带，median x≈0.455–0.48、p05–p95≈[0.20, 0.62]；Middle 75–80%；**Far (x>2/3) = 6/1958 = 0.306%**（corrected count；per-K 0.19–0.85%）；median |x−a|≈0.38–0.39（远在扫掠边界左侧）。
- **KS witnesses**：钉在 flat-top 的两个斜率转折边缘 **x≈1/3（19/25）与 x≈2/3（6/25）**，全部距 edge ≤0.013——不在幅度平台内部。
- **KW criterion-relevant risk**：prefix-restricted argmax 靠 **Near** 侧（全部轮 Near share 59%；top-25%-R_t 轮 100%）；raw argmax 24/32 轮落在结构性不可能的 suffix（E(K)=−D_max 平凡值）——由此确立 prefix-restricted 补充口径。
- **Multi-step finding**：rebound 后 risk/witness 位置大幅重建（KC median |Δx|=0.300；KS=0.335，两 edge 间交替）；**无 persistent hot spot**——one-step rebound generation 与 multi-step state reconstruction 是不同问题。这一发现直接动机化 Gate 2。
- 数据来源纪律：committed per-round CSV 无 state profile，audit 从冻结 config 确定性再生轨迹并对全部 5483 轮逐字符串核对 committed 字段（零失配）后才派生；无新 dynamics。

## 6. Gate 2 design logic（sign placement）

核心 pair：**XA（Witness-Aligned）与 XO（Witness-Opposed）**，满足逐点符号反射 `Δ_O(x) = −Δ_A(x)` ⟹ |Δ| 处处相同（max/min ±1/6 @ x=1/3、5/6）：

| | XA | XO |
|---|---|---|
| Δ 符号 | >0 on (0, 2/3)（覆盖历史 witness 带）；<0 on (2/3, 1) | <0 on (0, 2/3)；>0 on (2/3, 1) |
| G | 3/2x \| x/2+1/3 \| 2/3 \| 2x−1 | x/2 \| 3/2x−1/3 \| 2x−2/3 \| 1 |

同断点 {1/3, 2/3, 5/6}、同 regularity、同 piecewise complexity、同逐点 |Δ|、同初态、同控制器、同停止规则。**唯一操纵变量 = sign placement。** 这是 Stage 3 中最干净的 sign-placement comparison（XA-vs-KC 不是 sign effect identification）。

## 7. Gate 2 results

### XA（aligned）

全部 4 K：**t=6 exact stop、N_R=0**（max R_t 0.52–0.54）；risk 迁移到 x≈1/3（XA 自身正失配峰）。

> placing substantial positive mismatch over the historical canonical witness-support region is not sufficient to produce rebound or recover canonical long churn.

### XO（opposed）

全部 4 K：**`not stopped within pre-registered Gate-2 horizon`（finite-horizon censoring；Stage 3 首次）**。N_R=41/581/56/2380（K=50/100/200/400），f_R=0.082/0.581/0.028/0.595。呈现：

- **narrow D_max plateau**（≈0.1344 / ≈0.13778，近 K 无关，K=400 p5–p95=[0.1378,0.1380]）；
- **boundary locking**：整个 horizon 只取 ~4 个离散值（0.5、0.64、≈0.83、≈0.835），几乎无 large boundary switching（sw_freq≈0 vs canonical 0.92）；
- risk 锁定在 x≈5/6（XO 自身正失配峰）；actual rebound witnesses 全部在正支持区、Δ 严格 positive；
- **大量 numerical-boundary / exact-tie 轮**（3492/7524，全部在 XO；见 §9）；
- **强分辨率敏感性**（K=200 Short vs K=100/400 Long-span；f_R 差 20×）。

**措辞纪律（§11 核心）**：不得写 "XO reproduces canonical churn"。正确表述：

> **XO produces a long finite-horizon recurrent / marginal plateau regime under the current finite-K model.** Stage 3 uncovered a qualitatively different recurrent regime, not a replication of the canonical mechanism.

对照：canonical q_near churn = 跨 K 稳定宽运行带、substantial redistribution、大幅 near/far 边界跳转、gate exit/return、最终 annihilable stopping（finite stopping times 已在多个 K 观察）；XO = 窄平台、边界锁定、边际/exact-tie 反弹、几乎无大切换、强分辨率敏感、4/4 censoring。XO 的"最终是否停止"在本轮 scope 内**未知**——禁止 non-convergence / infinite churn / infinitely many rebounds 词汇。

### Gate 2 theory verification

7524 active 轮：分解残差 max 5.3×10⁻¹⁵；TP=2231、TN=4059；FP=407、FN=827 **全部 floating-boundary**（|worst| ≤ 1e-12；**non-boundary FP=FN=0**）；全部 3058 个 actual rebound witness 处 Δ 严格 positive。

## 8. Gate 2 main falsification（X4）

原 hypothesis：positive mismatch aligned with canonical historical witness-support region should be more conducive to rebound / recurrence。实际：

```text
XA (aligned):    zero rebound, t=6 stop, all K
XO (opposed):    long finite-horizon recurrent plateau, all K censored
```

> **static canonical witness alignment failed as a transferable predictor in the controlled XA/XO pair.**
> （不得写：historical witness localization has no predictive value in general。）

进一步 observation：**risk geometry rapidly relocates to each redistribution law's own positive mismatch support**（XA→1/3、XO→5/6）。当前更合理的 Working Hypothesis（非 theorem）：

> **witness/risk geometry is endogenous to the law-induced trajectory.**

## 9. Resolution / discretization concern（unresolved，External Review 高优先级）

XO 在 K=50/100/200/400 的 rebound frequency 与 recurrence signatures 高度不同，且 law 使用 rational breakpoints（1/3、2/3、5/6）。因此：

> **whether the XO plateau/tie-riding regime reflects a robust dynamical structure or a finite-grid arithmetic/resonance effect remains unresolved.**

本轮不解决；列入 External Review 优先问题（§13 Q6–Q8）。

## 10. Exact-tie / numerical-boundary concern

XO 大量轮次 |M·Δ − γ| ≪ 1e-12，并存在 real-arithmetic marginal ties（K=50 的 380 个 FP 轮中 377 个 actual `d_D_max` **逐位 = 0**）。必须区分三层：

1. numerical implementation 在冻结 tolerance 下**保持一致**（non-boundary FP=FN=0；分解残差 ≤1e-12；再生逐位核验）；
2. one-step theorem **没有发现任何非边界违例**；
3. 但 XO 的 qualitative long-horizon interpretation 受大量 marginal states 影响。

正确表述：**"numerically reproducible, but mathematically delicate."** 不得称当前结果"数值不可靠"。专业数学分析（exact / finite-state dynamics）具有高边际价值。

## 11. Stage 3 overall synthesis（five-layer story）

**Layer 1 — General one-step theory survived。** 跨 target-matched / target-dominated / weak positive / alternative shape / sign-changing XA / sign-changing XO 全部 kernels：one-step decomposition 与 exact rebound criterion 与实现/数值检验完全一致（约 1.3 万 active 轮、non-boundary 零误分）。

> the mechanism explaining why a particular step rebounds is genuinely broader than canonical q_near.

**Layer 2 — Positive mismatch alone does not organize long churn。** KW（正、零反弹）、KS（正、复发但短促）、XA（正、零反弹）：`Δ > 0` 既不保证 rebound 也不保证 long recurrent dynamics。

**Layer 3 — Simple scalar mismatch strength is insufficient。** KS 与 KC 匹配 max Δ 与 ∫Δ 但 long dynamics 极不相同：

> peak mismatch and integrated mismatch are not sufficient statistics for long-run transient organization.

**Layer 4 — Static witness geometry is not transferable。** XA/XO 直接 falsify 简单 historical alignment hypothesis：

> risk geometry must be treated as trajectory-dependent / endogenous.

**Layer 5 — Multiple recurrent regimes exist。** 目前至少两种定性不同的 recurrent regime：canonical q_near churn（hopping、宽 D_max 带、最终 annihilable stopping）与 XO marginal plateau / tie-riding regime（窄平台、边界锁定、censored）。研究焦点由此重构（**research framing / working hypothesis，非 theorem**）：

> **how redistribution law and evolving feedback state jointly organize multi-step recurrent transient regimes.**

## 12. What Stage 3 did NOT establish

- positive mismatch 的 general churn theorem；
- long recurrence 的 necessary/sufficient condition；
- amplitude threshold / phase transition；
- universal witness band / universal sign-placement rule；
- infinite recurrence / XO non-convergence（XO 最终行为未知）；
- continuum behavior；
- resolution-independent XO regime；
- real-world validity；
- novelty（见 §14）。

## 13. Consolidated External Review questions（供 Review Pack 压缩）

**Correctness / theory status**

- **Q1** Stage 2 的一步分解、G≤T 序定理、exact rebound criterion 是否属于已有标准结果/corollary？最自然的已有框架归属（discrepancy theory、selection-mutation/replicator、disturbance-rejection/ISS、deadband control 等）？
- **Q2** 这些定理是否存在 hidden assumptions 或更自然的定理表述（如对一般 target T、一般 removal 几何）？

**Long dynamics**

- **Q3** 什么 law/state 量能预测 long recurrence（而非仅 one-step rebound）？Gate 1/2 证据：R_t 分布不区分 KS 与 KC（p50 0.75 vs 0.70）但瞬态长度差 25×；risk 锁定各自 Δ 峰却不决定复发与否。
- **Q4** 是否存在 Lyapunov / operator / low-dimensional / finite-state representation？（witness 每轮大幅重建、无 hot spot；γ_t(j) 演化结构。）

**XO plateau / resolution**

- **Q5** XO 的 exact-tie plateau（real-arithmetic E≡0、边界锁定 ~4 值、D_max 平台 ≈0.134/0.138 近 K 无关）是否存在解析描述？
- **Q6** rational breakpoints（1/3、2/3、5/6）与 K 的组合是否产生 arithmetic/grid resonance？plateau 水平的近 K 无关性是结构还是巧合？
- **Q7** XO 是 robust dynamical regime 还是 finite-resolution artifact / hybrid？与 canonical churn 的本质区别是什么？
- **Q8** piecewise-flat G（密度含零段）在什么条件下产生精确边际 tie 与平台锁定？

**Positioning**

- **Q9** 最自然的领域归属与关键遗漏文献（Stage 1/2 定位线索：discrepancy minimization、online thinning、rank-driven selection-resampling、density-feedback swarm control、deadband control）？
- **Q10** 当前结果（Layer 1–5）是否值得进一步数学严格化？若值得，最小严格化路径是什么？

（合并自 Gate 1 ERC-1–4、audit ERC-5、Gate 2 ERC-6–7，避免重复编号。）

## 14. Novelty position

> **novelty not established**（不变）。

Stage 3 的 strong falsification / 新 regime 不自动构成 novelty claim。更准确的 potential contribution（potential positioning，非 verdict）：

> a carefully governed model study showing that one-step mismatch theory generalizes broadly while long recurrent transients depend on richer law–state interaction and can organize into qualitatively distinct finite-K regimes.

## 15. Reality connection remains open

real throw distribution 未测量；perceptual epsilon 未校准；2D 未建立；human selection noise 未建模；finite-particle stochasticity 未系统研究。Stage 3 是数学 generality study，**不改变也不推进**以上任何一项。

## 16. Evidence hierarchy snapshot（Stage 3 增量）

| Result | Scope | Evidence class |
|---|---|---|
| 一步分解 / G=T 收缩 / G≤T sufficient class / exact rebound criterion | general M1 update class | **Model-derived Mathematical Result**（Stage 2 建立并经 Stage 3 全 kernel 逐轮验证；a_crit 仍属 canonical specialization） |
| K-：t=3、零反弹、严格单调 | K- @ K=100 | Verified Simulation Result（+定理后盾） |
| KW：零反弹、t=8（4 K）；KS：5–10 rebound、t=26–29（4 K） | Gate 1 frozen configs | Verified Simulation Result |
| XA：零反弹、t=6（4 K）；XO：4 K 10K censoring + plateau；XA/XO contrast | Gate 2 frozen configs | Verified Simulation Result |
| Gate 1 witness localization（KC band [0.20,0.62]、KS edges、KW Near risk；Far=6/1958） | audited kernels × 4 K | Verified Numerical Analysis |
| support migration（risk 锁定自身正峰）；decomposition residuals；non-boundary FP/FN=0；XO tie-riding diagnostics | Stage 3 | Verified Numerical Analysis |
| canonical churn regime；KS edge-pinned witnesses；XO marginal plateau；resolution sensitivity | canonical/相应 kernels | Model-dependent Observation |
| multi-step recurrence 需要 law × evolving-state geometry；risk geometry endogenous；one-step 标量分布不决定 transient lifetime | Stage 3 | Working Hypothesis |

## 17. Explicitly falsified / superseded（承接 Stage 2 F1–F8）

| # | 已废弃解释 | 推翻证据 | 现行解释 |
|---|---|---|---|
| F9 | positive cumulative mismatch alone broadly produces canonical-like churn | KW（零反弹）+ KS（短促复发）+ XA（零反弹） | positive mismatch 可产生无 rebound、短促复发或其他 dynamics，取决于 amplitude/profile/state interaction |
| F10 | peak mismatch 与 integrated mismatch 足以刻画 long transient behavior | KS vs KC（同峰同积分、瞬态长度差 ~25×） | 标量失配强度不是 long-run 组织的充分统计量；profile/regularity 与 law–state 交互参与 |
| F11 | canonical historical rebound-witness locations 构成对新 laws 可迁移的 intervention map | XA/XO controlled pair（aligned 零反弹 vs opposed 长复发） | witness/risk geometry 是 trajectory-dependent 的、在新 laws 下重组（endogenous，Working Hypothesis） |
| （guardrail，非 falsified claim） | 任何 long recurrent sign-changing behavior 都应解读为 canonical-like churn——旧文档从未如此声称；列为今后解读护栏 | XO 的形态学差异（平台/锁定/边际 rebound/censoring vs hopping/宽带/annihilable stopping） | 多个定性不同的 recurrent regimes 可能并存；不得混称 "churn" |

另有本轮 A2 bookkeeping 更正（audit 文档 KC Far share prose 笔误 31/1958 → 6/1958）与 Gate 1 报告措辞收紧（否定对象限定为 "positive mismatch alone broadly produces canonical-like churn is not supported"，不升级为否定 model-class 现象）——均已留痕。

## 18. Stage 3 closure statement

> **Stage 3 / M2 — Redistribution-Mismatch Generality: CLOSED.**

含义：theory-guided generality gate 已完成（positive / negative / weak / alternative-shape / sign-changing controls 全覆盖）；one-step theory generality 与 long-transient non-generality 已清晰区分；不再自行追加 kernel；§9–§10 的 resolution 与 exact-tie 问题、§13 的理论问题**转入 External Mathematical Review Gate**。Stage 2 checkpoint 与全部 raw results 保持冻结。

## 19. Next stage

> **External Mathematical Review Gate**（下一任务：prepare Expert Review Pack；不联系 reviewer、不启动 review 本身）。

Stage 4：**NOT FROZEN**。候选路线（由 external review 结果决定）：further general mathematical theory；collaborator；literature / contribution reframing；undergraduate research/modeling report；reality-inspired validation / 2D / perception。

## 20. New-window bootstrap read order

1. 本文件（Stage 3 closure）；
2. [RESEARCH_CHECKPOINT_2026-09-20.md](RESEARCH_CHECKPOINT_2026-09-20.md)（Stage 2 closure，theory 定义）；
3. [M2_GATE1_GENERALITY_GATE.md](M2_GATE1_GENERALITY_GATE.md) + [M2_GATE2_SIGN_PLACEMENT.md](M2_GATE2_SIGN_PLACEMENT.md)（Stage 3 两轮）；
4. [M2_GATE1_INTERPRETATION_AUDIT.md](../experiments/m2_gate1_generality/M2_GATE1_INTERPRETATION_AUDIT.md)（witness geometry）；
5. `MODEL_ASSUMPTIONS.md`（§11–§15 kernel 登记）、`NEXT_STEPS.md`（当前优先）、`RESEARCH_LOG.md`（append-only 历史）。
