# 沙子再分配研究：公开版研究报告 / Sand Redistribution Study: Public Research Report

> 发布快照说明 / Release snapshot note (2026-10-01): This is a research record for a deterministic one-dimensional model, not a physical sand simulator or a peer-reviewed paper. The Chinese and English accounts below describe the same evidence. The [formal XO proof](S3_XO_BOUNDARY_THEOREM_PROOF.md) is the theorem-level source; dated experimental reports retain their original, time-limited wording.

## 中文

### 问题与模型

项目起于一个现实观察：人看到沙面局部堆积后，把一段沙扫回、重新撒出，再判断是否继续。仓库研究的是由此抽象出的**一维、有限维、确定性质量反馈模型**，不是沙粒碰撞、人的动作或真实沙面的数值仿真。

把总质量 1 分在 \(K\) 个格子。第 \(j\) 个前缀的质量为 \(F_{t,j}\)，相对均匀目标的**单侧累计超额**为

\[
D_{t,j}=F_{t,j}-j/K,\qquad d_t=\max_j D_{t,j}.
\]

每轮选择达到最大超额的最小指标 \(s_t\)，若 \(d_t\le10^{-12}\) 则按模型的**容差规则**停机；否则从前 \(s_t\) 格各取走四分之一质量，按固定喷撒分布重新投入所有格子。\(d_t=0\) 只表示没有前缀超过其目标份额，**不等于整个质量向量均匀**。模型、更新式及假设见 [MODEL_ASSUMPTIONS](MODEL_ASSUMPTIONS.md) 和 [XO 证明 §1](S3_XO_BOUNDARY_THEOREM_PROOF.md)。

### 已得到什么

**一般单步数学（[Math]）**：在模型所定义的更新类内，新旧分布的差别可精确拆成“目标匹配更新”加上“本轮移走的质量 × 喷撒律相对目标的累计失配”。目标匹配时，活跃一步的最大超额严格下降；喷撒律累计分布处处不高于目标时，也有严格下降的充分条件。精确的单步反弹判据需要当前完整状态，**不是多轮寿命定理**。正式推导见 [M1C.2](M1C2_TARGET_MATCHED_KERNEL_THEORY.md) 与 [M1C.3](M1C3_REDISTRIBUTION_LAW_ORDER.md)。

**对照实验（限定范围的计算证据）**：Stage 2/3 比较不同的固定喷撒律，发现“有正失配”“失配峰值相同”“失配积分相同”等单一特征不足以决定是否出现 canonical 轨迹那样的长时程反复修正。它们是指定格数、初态、规则与运行时限内的结果，不构成所有喷撒律的分类定理。见 [Stage 2 checkpoint](RESEARCH_CHECKPOINT_2026-09-20.md) 与 [Stage 3 checkpoint](RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md)。

**canonical XO 定理（[Existing Math] T1；[New Mathematical Result] T2–T5）**：固定 canonical 初态、XO 喷撒律、移除比例 \(1/4\)、最小指标平局规则、精确算术和整数 \(K\ge6\)。设

\[
m=\lceil5K/6\rceil,\quad \ell=m-1,\quad h=1/K,\quad x_t=p_{t,m}.
\]

原不变量 T1 说明 \(F_{t,m}\) 守恒、\(s_t\le m\)，且 \(d_t\) 有严格正下界，因此这条精确轨迹不会触发当前容差停机规则。新证明 T2–T3 说明轨迹有限步后进入并持续留在只有 \(\ell,m\) 能取得最大超额的区域。此后 \(x_t\) 满足一个**精确的双分支仿射递推**，由它可恢复边界选择、\(d_t\)、平局和本轮移动质量；它**不能**恢复整个质量向量。特别地，\(x_t=h\) 时选较小的 \(\ell\)，而 \(d_t\) 在该步保持相等。

定理 T5 进一步给出无限时间的选择二分：\(K\bmod6\in\{0,1,2\}\) 时，最终永远选择 \(m\)；\(K\bmod6\in\{3,4,5\}\) 时，\(\ell\) 和 \(m\) 都会被无限次选择。**这不是周期性结论**。见[英文正式证明](S3_XO_BOUNDARY_THEOREM_PROOF.md)、[中文完整译本](S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md)、[定理登记](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)及简短的 [XO 定理导读](XO_THEOREM_GUIDE_ZH.md)。

**精确计算（[VCR]）**：S1 以有理数重算指定的四个 \(K\)，发现一些浮点“反弹”其实发生在精确持平轮。S2 在 \(K=6,\ldots,53\) 及 \(100,200,400\) 共 51 个配置上，逐步验证边界递推，并在冻结的 \(H=10K\) 时限内得到 51/51 的余数类分类。这些计算帮助发现并检验机制，不能替代上述一般 \(K\) 证明。实验标签 **LOCK/CYCLE 只属于该有限时限分类器**；定理说的是无限时间选择行为，未额外证明未测试的每个 \(K\) 在 \(H=10K\) 时已经取得对应标签。见 [S1](S4_W1_S1_XO_PRECISION.md) 和 [S2](S4_W1_S2_BOUNDARY_MICROSCOPE.md)。

### 为什么有意义，什么仍未知

意义在于：一个由完整质量分布驱动的选择过程，在这套特定规则及轨迹上，最终有一组重要**边界观测量**可由单个数精确决定；格子数与 XO 支集终点的对齐，又把选择行为分成六个余数类中的两组。它提供了可核查的结构，而不只是图上的平台形状。

该结论没有覆盖其他初态、其他移除比例、其他喷撒律、二维、随机扰动、真实沙坑，也没有给出完整状态的渐近描述或周期轨道定理。一般文献中的新颖性尚未系统确立。项目目前完成了 canonical XO 的 S3 证明与仓库内审计记录；后续研究分支没有自动启动。[审计来源说明](S3_XO_BOUNDARY_THEOREM_AUDIT.md)如实记录：独立审计的 PASS 结论由项目所有者提供，完整审计对话未收入仓库。

## English

### Question and model

The project began with a physical observation: a person sees an uneven sand surface, sweeps back an over-dense region, redistributes that material, and decides whether to continue. The repository studies a **one-dimensional, finite-dimensional, deterministic mass-feedback abstraction** of that process. It does not simulate grain collisions, body motion, or a real sand surface.

Unit mass is placed in \(K\) bins. For prefix \(j\), let \(F_{t,j}\) be cumulative mass and let \(D_{t,j}=F_{t,j}-j/K\) be its **one-sided cumulative excess** over a uniform target; \(d_t=\max_jD_{t,j}\). Each round selects the smallest maximizer \(s_t\). The model stops when \(d_t\le10^{-12}\); otherwise it removes one quarter of the mass from each selected-prefix bin and redistributes that mass according to a fixed law. This is a **tolerance stopping rule**. In particular, \(d_t=0\) does not imply a uniform mass vector. See [model assumptions](MODEL_ASSUMPTIONS.md) and [XO proof §1](S3_XO_BOUNDARY_THEOREM_PROOF.md).

### Established results and evidence levels

**General one-step mathematics ([Math]).** Within the defined update class, the updated discrepancy is exactly a target-matched update plus removed mass times the redistribution law's cumulative mismatch from the target. Target-matched active steps strictly decrease maximal excess; cumulative laws below the target satisfy a sufficient strict-descent condition. The exact one-step rebound criterion depends on the full current state and is not a lifetime theorem. See [M1C.2](M1C2_TARGET_MATCHED_KERNEL_THEORY.md) and [M1C.3](M1C3_REDISTRIBUTION_LAW_ORDER.md).

**Controlled computational comparisons.** Stage 2/3 controls show, within their specified grids and horizons, that the existence, peak, or integral of positive mismatch alone does not determine a long correction transient like the canonical trajectory. These are scoped computational results, not a classification of every redistribution law. See the [Stage 2](RESEARCH_CHECKPOINT_2026-09-20.md) and [Stage 3](RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md) checkpoints.

**Canonical XO theorem ([Existing Math] T1; [New Mathematical Result] T2–T5).** The theorem fixes canonical initialization, the XO law, removal fraction \(1/4\), smallest-index tie handling, exact arithmetic, and integer \(K\ge6\). Put \(m=\lceil5K/6\rceil\), \(\ell=m-1\), \(h=1/K\), and \(x_t=p_{t,m}\). The earlier invariant T1 conserves \(F_{t,m}\), bounds selection by \(m\), and keeps discrepancy strictly above the current stopping tolerance. The new proof T2–T3 gives finite entry into, and persistence of, a regime where only \(\ell\) and \(m\) can maximize discrepancy. In that regime, an exact two-branch affine recurrence for \(x_t\) determines boundary selection, discrepancy, ties, and moved mass. It does **not** reconstruct the full mass vector. At \(x_t=h\), the smaller index \(\ell\) is selected and discrepancy remains tied.

T5 establishes an infinite-time selection dichotomy: for \(K\bmod6\in\{0,1,2\}\), selection eventually remains at \(m\); for residues \(\{3,4,5\}\), both \(\ell\) and \(m\) are selected infinitely often. It does **not** assert periodicity. See the [formal English proof](S3_XO_BOUNDARY_THEOREM_PROOF.md), its [full Chinese translation](S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md), the [registration and exact scope](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md), and the shorter [Chinese theorem guide](XO_THEOREM_GUIDE_ZH.md).

**Exact computation ([VCR]).** S1 replayed four selected \(K\) values with rational arithmetic and separated genuine rebounds from floating-point signals on exact-tie rounds. S2 verified the boundary recurrence step by step for 51 values, \(K=6,\ldots,53\) plus \(100,200,400\), and found a 51/51 residue classification under the frozen horizon \(H=10K\). This evidence supported discovery and validation, not the universal proof. The experimental labels **LOCK/CYCLE** refer only to that finite-horizon classifier. The theorem's infinite-time dichotomy does not additionally establish the classifier's label at \(H=10K\) for every untested \(K\). See [S1](S4_W1_S1_XO_PRECISION.md) and [S2](S4_W1_S2_BOUNDARY_MICROSCOPE.md).

### Value and limits

The result is a precise structural reduction of **boundary observables** along one canonical trajectory. After entry, one bin mass controls the relevant boundary choices, while the alignment of the XO support endpoint with the \(K\)-bin grid yields the residue split. This is stronger than a visual plateau pattern and narrower than a reduction of the full system.

Other initial states, removal fractions, laws, 2D models, stochastic effects, real sand, full-state asymptotics, and periodic trajectories are outside the theorem. Novelty relative to the mathematical literature has not been established. Canonical XO S3 is closed in the repository; other research branches have not started automatically. The [audit provenance](S3_XO_BOUNDARY_THEOREM_AUDIT.md) states that the independent PASS verdict was supplied by the Owner and that the full audit transcript is not deposited.
