# 公开研究报告 / Public research report

> v0.2.0，2026-10-07。当前补充成果见[登记](XO_EXTENSION_REGISTRATION.md)；原 S3、S1/S2 报告和冻结材料保留原范围。

## 中文

### 问题与模型

研究起点是对不均匀沙面“扫回一部分再撒开”的观察。仓库采用一维有限维确定性质量反馈抽象，不模拟真实沙粒碰撞或沙坑。K 个格子内总质量为 1，前 j 格的累计质量为 F_j，单侧累计超额 D_j=F_j−j/K，每轮取最小指标的最大超额。更新前 d≤10⁻¹² 则停止；否则从所选前缀每格移除 α 份额，再依固定律撒回。原实验默认 α=1/4；A 在原 canonical XO 前提下推广 α。d=0 不意味着完整质量向量均匀。模型和原假设见 [MODEL_ASSUMPTIONS](MODEL_ASSUMPTIONS.md)。

### 结果与证据

**一般单步数学：** 更新可分解为目标匹配更新，加上“移走质量 × 再分配律累计失配”。目标匹配时活跃步最大超额严格下降；累计律不高于目标也给出充分下降条件。完整状态的一步反弹判据不是多轮寿命定理。见 [M1C.2](M1C2_TARGET_MATCHED_KERNEL_THEORY.md)、[M1C.3](M1C3_REDISTRIBUTION_LAW_ORDER.md)。

**历史对照计算：** Stage 2/3 的指定规则、初态、格数与时限比较说明，仅有失配峰值或积分等特征不足以决定 canonical 轨迹式的长瞬态。这是有限计算，不是所有再分配律的分类。

**原 canonical XO S3：** canonical 初态、XO、α=1/4、最小指标平局、精确算术、K≥6。T1 给出前 m=ceil(5K/6) 格质量守恒及正偏差下界；T2/T3 给出有限进入并持续留在 {m−1,m} 边界区；T4 给出 x_t=p_{t,m} 的精确双分支闭包，决定边界选择、偏差和移动质量，不能重建完整质量向量。T5 的余数类二分为 0/1/2 最终永久选择 m，3/4/5 两边界被无限次选择，没有周期性声明。原[英文证明](S3_XO_BOUNDARY_THEOREM_PROOF.md)、[完整中文译本](S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md)及[登记](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)保持原范围。

**补充 A：** 同样的 canonical XO、平局和精确算术前提，K≥6、任意 0<α<1。A 证明上述长期边界结构，首次进入有保守界 floor(81K/(80α))+1；K=6 的首次进入持久性单独处理。见 [A 正式证明](XO_ALPHA_EXTENSION_PROOF.md)。

**补充 T6：** 仅 α=1/4、H=10K，窗口含状态 t=0,…,H−1 的 H 项选择，使用原 S2 冻结分类函数。对所有 K≥6，有限窗 LOCK/CYCLE 与长期余数二分一致，没有 OTHER。T6 依赖原 S3，自身驻留与计数论证不需要 A。见 [T6 正式证明](XO_T6_CLASSIFIER_PROOF.md)。

**A/T6 范围分开：** K=9、α=1/8、H=90 时，分类器因 R_m=10≥K 优先输出 LOCK，但长期两边界仍无限切换；K=10、α=1/1000、H=100 时可未进入并输出 OTHER。A 不能扩展成任意 α 的 10K 分类保证。

**计算核验：** 原 S1/S2 的 51 个 K 记录保持不变。补充脚本核对 A 的 255 组精确矩阵（大 K 实际只跑 300 步），T6 的 295 个实例和 15 项条件；独立整数完整状态算法另跑 29 条 α=1/4 轨迹，含 K=200、400 的完整 10K 窗口。有限浮点标量代理与 exact 完整状态分别登记，均不能替代全称证明。实际整合结果见[验收](XO_EXTENSION_VALIDATION.md)。

**辅助附录：** F 的单调地板/条件守恒、首次选择公式、限定有效初态与强制规则探针、径向聚合/半共轭和文献线索。S12 的 200000 步统计仅为标量 float 代理；完整状态 exact 33 步核对另列。C1 只对径向选择、固定非负归一核、正整数 R/T、0<α<1 及径向停机成立；角向偏差按各环被扫次数衰减。详见[辅助索引](XO_SUPPLEMENTARY_RESULTS.md)。

### 意义与仍未知的范围

成果说明，特定质量反馈轨迹的一组边界观测量可由单个质量坐标精确决定，支集终点与网格对齐给出长期余数二分，原 1/4 冻结窗口也有一般 K 的可靠性证明。长期结构的操作力度推广和有限窗分类保证是两项不同结论。

周期轨道、完整状态渐近、其他初态/再分配律、一般噪声、一般二维和现实模型仍未由这些结果解决。书目线索没有完成内容级查新，学术新颖性尚未确立。原 S3 的[审计来源](S3_XO_BOUNDARY_THEOREM_AUDIT.md)与本次[独立补充审查](XO_EXTENSION_REVIEW.md)分开；两者均不是期刊同行评审声明。

## English

### Question and model

The motivating observation was a person correcting an uneven sand surface by sweeping some material back and redistributing it. The repository studies a one-dimensional, finite-dimensional, deterministic mass-feedback abstraction, not grain collisions or a real sandpit. Unit mass occupies K bins. Prefix discrepancy is D_j=F_j−j/K; each round selects the smallest maximizer. The process stops before updating when maximal discrepancy is at most 10⁻¹²; otherwise it removes an α fraction from the selected prefix and redistributes it according to a fixed law. Historical experiments use α=1/4. A extends the removal fraction under canonical XO assumptions. Zero one-sided discrepancy does not imply a uniform full mass vector.

### Results and evidence

**General one-step mathematics.** The update separates into a target-matched update and removed mass times cumulative redistribution mismatch. Target matching gives strict descent on active steps; cumulative laws below the target give another sufficient condition. A full-state one-step rebound condition is not a lifetime theorem. See M1C.2 and M1C.3 above.

**Historical controlled computations.** Stage 2/3 comparisons show that mismatch peaks or integrals alone do not determine a canonical-style long transient within the tested laws, initial states, grids, and horizons. They do not classify every redistribution law.

**Original canonical XO S3.** With canonical initialization, XO, α=1/4, smallest-index ties, exact arithmetic, and K≥6, T1 gives prefix conservation and a positive discrepancy floor. T2/T3 give finite entry into and persistence of the {m−1,m} boundary regime, where m=ceil(5K/6). T4 gives an exact two-branch recurrence for x_t=p_{t,m}, determining boundary choices, discrepancy, and moved mass, but not the full mass vector. T5 gives eventual permanent selection of m for residues 0/1/2 and infinitely many selections of both boundaries for residues 3/4/5. It does not assert periodicity. The original proof, translation, registration, and audit retain their scope.

**Supplement A.** Under the same canonical XO and exact-tie assumptions, A extends the long-term structure to every 0<α<1 and K≥6. The conservative entry bound is floor(81K/(80α))+1; the K=6 first-entry case is treated separately. See the deposited A proof and registration.

**Supplement T6.** At α=1/4 and H=10K only, the original frozen S2 classifier observes H pre-update choices at t=0,…,H−1. For every K≥6, LOCK/CYCLE agrees with the long-term residue dichotomy, with no OTHER. The proof depends on original S3 and its own dwell/counting bounds, not on A.

**Separate scopes.** At K=9, α=1/8, H=90, a tail of ten m choices makes the classifier return LOCK although both boundaries occur infinitely often in the long term. At K=10, α=1/1000, H=100, no entry occurs and the label is OTHER. A does not supply T6's finite-window guarantee at arbitrary α.

**Computational validation.** The original 51-case S1/S2 record is preserved. Supplement checks include A's 255 exact cases, whose large-K windows are capped at 300, and T6's 295 cases plus 15 exact conditions. A separate common-denominator integer algorithm advances 29 full mass-vector trajectories at α=1/4, including complete 10K windows for K=200 and 400. Finite scalar float proxies and exact full-state checks are separately labelled. They support, rather than replace, the analytical universal statements. The validation page records integrated results.

**Scoped appendices.** The release also includes a flat-tail monotone floor/conditional-conservation proposition, a first-update selection formula, valid finite initial-state probes, a restricted forced-rule experiment, radial aggregation/semiconjugacy, and bibliographic leads. S12's 200000-step statistics advance a scalar float proxy; the exact 33-state full-vector check is separate. C1 assumes radial-only selection/stopping, a fixed normalized nonnegative kernel, positive integer dimensions, and 0<α<1; angular decay is counted by sweeps of each ring.

### Value and limits

The work gives a precise reduction of boundary observables along a specified mass-feedback trajectory, a residue-class long-term dichotomy, and a reliable finite-window classifier at the original α=1/4. The arbitrary-α long-term extension and the finite-window theorem are distinct.

Periodicity, full-state asymptotics, other initial states/laws, general stochastic or 2D theory, and real-world validation remain open. Literature leads do not establish novelty. Original S3 audit provenance and the separate A/T6 review are recorded individually; neither is a journal peer-review claim.
