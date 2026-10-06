# S3 — Canonical XO Boundary-Reduction Theorem：正式登记文档

> 登记 / proof deposit 日期：2026-09-30。性质：已完成、已审计证明的仓库入库与证据分层 reconciliation；无新研究、无新定理范围。
>
> **登记状态**：`REGISTERED — PROOF DEPOSITED / INDEPENDENT AUDIT PASS`。
>
> **S3：CLOSED — PASS**。证明本体：[S3_XO_BOUNDARY_THEOREM_PROOF.md](S3_XO_BOUNDARY_THEOREM_PROOF.md)。独立审计：[S3_XO_BOUNDARY_THEOREM_AUDIT.md](S3_XO_BOUNDARY_THEOREM_AUDIT.md)，GPT-6.1 Sol，**PASS — no substantive gap found**。
>
> **状态迁移 provenance**：`04dab39` 是私人来源仓库中的登记提交，当时仅有定理陈述，证明尚未入库。后来先存入 Astra 的完整 P1–P4 证明及独立审计 provenance，再完成登记升级；不把历史实验改写为证明。私人提交编号与 RESEARCH_LOG 仅为来源记录，不是本公开独立历史中的提交或证明依赖；公开证明 §2 已完整重录既有不变量证明。引用本版时使用公开标签与文件路径，见[阅读说明](READING_NOTES.md)。 / The commit ID belongs to the private source history. Its log records provenance, not a dependency of the deposited public proof; cite this snapshot using a public tag or commit and file path.

## 1. 冻结适用范围

全部限定同时适用：finite-dimensional deterministic M1 mass-feedback model；canonical initialization（CDF 2x − x² 的精确 bin 差分）；XO CDF（x/2 ｜ 3x/2 − 1/3 ｜ 2x − 2/3 ｜ 1，断点 1/3、2/3、5/6）；alpha = 1/4；smallest-index argmax；exact arithmetic；integer K ≥ 6。

不扩展至其他 laws、其他初态、其他 alpha、2D、float/mp 实现或一般 sand redistribution 系统。完整模型定义与更新式见证明 §1。

## 2. 已核对符号

| 符号 | 定义 | 证明位置 |
|---|---|---|
| h | 1/K，均匀 bin 份额 | §1 |
| r | K mod 6 | §1 |
| m | ⌈5K/6⌉ | §1 |
| ℓ | m − 1，边界左指标 | §1 |
| b | m/K | §1 |
| c | 2b − b²，守恒质量 F_{t,m} | §§1–2 |
| a | b(1−b) = c−b，守恒超额 D_{t,m} | §§1–2 |
| q | q_K = (6−r)/(3K)，边界 bin 重撒份额 | §1 的闭式推导 |
| x_t | p_{t,m} | §§1、6 |
| t* | min{t : s_t ∈ {ℓ,m}} | §§1、3–5 |

## 3. 主定理 T1–T5 与证明对应

**T1 — Global invariant [Existing Math]**：对所有 t，s_t ≤ m、F_{t,m}=c、D_{t,m}=a，且 d_t ≥ a ≥ 10/121 > 10⁻¹²，故 tolerance stopping rule 不触发。原证明已完整重录于[公开证明 §2](S3_XO_BOUNDARY_THEOREM_PROOF.md#2-existing-invariant-dependency--t1-existing-math)，不依赖未公开的旧评审草稿。F_{t,m}=c 是既有不变量的 CDF 形式重述，不另算新发现。

**T2 — Boundary entry [New Mathematical Result]**：t* 存在且有限，t* ≤ ⌊81K/20⌋+1。完整证明见 proof §3（P1）。该宽松上界不宣称 S2 观察到的 t*∈{1,2} 对所有 K 成立。

**T3 — Boundary persistence [New Mathematical Result]**：对每个 t≥t*，argmax_j D_{t,j} ⊆ {ℓ,m}。证明见 proof §4（K≥7 的端点不等式）及 §5（K=6 canonical 两步桥接）。结论从第一次边界 selection 起成立，不只从较强不变区域的到达时间起成立。

**T4 — Scalar boundary closure [New Mathematical Result]**：对 t≥t*，

- x_{t+1} = (3/4)x_t + cq/4，当 x_t>h；
- x_{t+1} = x_t + (c−x_t)q/4，当 x_t≤h；
- s_t=m 当且仅当 x_t>h，否则 s_t=ℓ；d_t=a+(h−x_t)⁺。

证明见 proof §6（P3）。闭包决定边界坐标、selection、discrepancy、平局、移除质量与 discrepancy increments；**不重构完整质量向量**。精确 threshold equality x_t=h 归 ℓ 支，但 d_{t+1}=d_t=a，是 exact tie；“B 支严格收缩”需要 x_t<h。

**T5 — Residue-class dichotomy [New Mathematical Result]**：

- K mod 6 ∈ {0,1,2}：eventual permanent m-selection，即存在 T，使所有 t≥T 均有 s_t=m；
- K mod 6 ∈ {3,4,5}：s_t=ℓ 与 s_t=m 均出现无限多次。

证明见 proof §7（P4）：标量判据为 cq≥h 与 cq<h；canonical XO 中等号不出现。跨向下支包括落在 x=h。**不推出周期性、周期轨道、渐近 dwell 常数或 full-state asymptotics。**

## 4. 独立审计与证据分层

独立审计方 **GPT-6.1 Sol**；结论 **PASS — no substantive gap found**。来源为 Owner 在 2026-09-30 deposit 请求中的明确报告。审计范围包括 dependency assumptions、P1–P4、K=6 exceptional bridge、indexing 与 threshold equality，沿用先前登记的范围。完整审计原始对话未另行提供；仓库 [audit provenance record](S3_XO_BOUNDARY_THEOREM_AUDIT.md) 如实保存已提供的身份、verdict、范围、记录日期与来源，不虚构逐条 reviewer 意见。

| 层级 | 内容 | 仓库证据 |
|---|---|---|
| **[Existing Math]** | 原 XO support-endpoint invariant，T1 | 公开证明 §2 完整重录原证明 |
| **[New Mathematical Result]** | Canonical XO boundary-reduction theorem，T2–T5 | 已入库 proof §§3–7 |
| **[Independent Audit]** | GPT-6.1 Sol：PASS — no substantive gap found | 已入库 audit provenance record（Owner 提供结论） |
| **[VCR]** | S1 exact replay；S2 51-case census 与 C1–C7 | 原 committed outputs 不变；发现与验证支撑，不代替证明 |
| **[Obs]/[Hyp] / open** | 未覆盖的 law、初态、alpha、2D、周期性、full-state asymptotics、dwell 结构 | 本轮不研究、不升级 |

## 5. Deposit completion 与验证记录

1. 证明本体已入库：含 P1–P4、h/ℓ 定义、K=6 桥接、threshold equality、entry bound。
2. 独立审计 provenance 已入库：reviewer、记录日期、verdict、范围与来源明确；不冒充完整审计 transcript。
3. T1–T5 与证明章节逐项对应，h=1/K、ℓ=m−1 已正式核对。
4. Scope 保持 canonical XO / canonical initialization / alpha=1/4 / smallest-index / exact arithmetic / K≥6。
5. Deposit round 对已有 51 K CSV 做机械一致性核对，并运行既有完整测试入口；2026-09-30 的记录为 20 suites / 382 项运行时断言通过。原始实验与冻结 configs 不变。公开版另行验证，不把本条当作新测试。

**REGISTERED — PROOF DEPOSITED / INDEPENDENT AUDIT PASS；S3 CLOSED — PASS。**

## 6. 定理与实验标签的分界

定理使用 eventual permanent m-selection / both occur infinitely often。LOCK/CYCLE 仍仅表示冻结分类器在 H=10K 内的实验分类。此证明**没有额外证明**该有限 cutoff 分类器对每个未测 K 都与无限时间二分一致。

S1/S2 历史报告保留其发现当时的 [VCR] 范围；其 open proof targets 由本登记与证明覆盖的部分取代。无 full-state 一维化、无周期性结论、无新的“共振”或 dwell 极限主张。

## 7. Closure

发现与验证链：既有 invariant → S1 exact precision correction → S2 scalar boundary evidence。证明链：S3 P1 entry → P2 persistence（含 K=6）→ P3 closure → P4 residue dichotomy。独立审计 PASS 后，证明现已完成仓库 deposit。

Canonical XO 的 P1–P4 理论缺口闭合。下一步是 synthesis / documentation；任何其他研究分支由 Owner 另行决定。本 deposit 不自动启动 S4/S8、2D、噪声、dwell 渐近或新的 proof round。
