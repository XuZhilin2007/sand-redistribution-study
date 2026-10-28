# 阅读说明、术语与勘误 / Reading notes, terminology, and errata

> 更新 / Updated: 2026-10-27：追加下列 E1–E3 勘误入口与引用清单。原 2026-10-07、v0.2.0 阅读说明保留。 / Added the E1–E3 errata index below; the earlier v0.2.0 reading notes are retained.

## 2026-10-27 定向勘误 / Targeted corrections

当前数学更正以 [2026-10-27 评审勘误 E1–E3](review_2026-10-27/ERRATA.md) 为准。以下明确其适用位置；日期明确的历史文稿、实验脚本、原始输出和发布快照保留原样，不把历史原句作为未经限定的当前结论引用。勘误中的反例与推导已经给出，本次只补足维护入口和定向复验，见[修正验收记录](corrections_2026-10-27/VERIFICATION.md)。

The [review errata E1–E3](review_2026-10-27/ERRATA.md) give the corrected mathematics. The index below identifies the historical passages to read with those corrections. Dated documents, experiment sources and outputs, released snapshots and the completed review remain unchanged. The [targeted verification record](corrections_2026-10-27/VERIFICATION.md) records this maintenance pass.

| 勘误 / Erratum | 正确范围 / Correct statement | 原文与重复引用 / Historical locations |
|---|---|---|
| E1：M0 整数轮最优性 | `c_n > c_f` 等价于连续路径顶点 `z* < 1`；存在严格优于初态的正整数轮，当且仅当 `z* < (1+ρ)/2`，其中 `ρ=1−α(1−Q)`、`0<Q<1`、`0<α≤1`。等号时 t=0、1 并列，α=0 单独为静止轨道。 | [M0_GENERAL_Q](M0_GENERAL_Q.md) 引言、§1.2、§1.4、§A.5–A.6；[M0.2 报告](../experiments/m0_2_distribution_family/M0_2_REPORT.md) §2.1、§2.2、§7、§8；[09-19 checkpoint](RESEARCH_CHECKPOINT_2026-09-19.md) 的一般判据及幂族条目；[09-20 checkpoint](RESEARCH_CHECKPOINT_2026-09-20.md) 的 M0 结果表；[研究日志](RESEARCH_LOG.md) 的 Phase 2A 判据。上述宽泛表述均按 E1 限定。 |
| E2：iid 直方图损失期望 | 对历史未乘 K 的损失 `V(p)=Σ(p_i−1/K)²`，`E V(p̂)=V(p)+(1−Σp_i²)/N=(1−1/N)V(p)+(1−1/K)/N`。固定均匀噪声项只在均匀处精确成立，在一般状态是上界近似。 | [M0_MECHANISM](M0_MECHANISM.md) §2.7；[M0.1 实验脚本](../experiments/m0_1_mechanism/run_experiment.py) 的 MC 均值比较及其精确期望注释。历史统计值保留，不能据该近似比较声称已验证错误的精确等式。 |
| E3：B 持久性证明的斜率上界 | 由 `w≤C<1` 得 `w/2−1<−1/2<0`；由 `w≥29/42` 得 `3w/2−1≥1/28>0`、`2w−1≥8/21>0`。这给出所需的先降后升，与 `C>Θ(K)` 的端点比较一起恢复原论证。 | [V03_B_ROBUSTNESS_PROOF](V03_B_ROBUSTNESS_PROOF.md) §2(3)。错误数值上界不使用；定理的假设和结论保持原范围。 |

E1 同样限定 [M0.2 历史预测脚本](../experiments/m0_2_distribution_family/run_experiment.py) 中的 `predicted_interior_optimum`：该字段实际检测连续顶点 `z*<1`（代码含浮点容差），不能当作任意 α 下整数轮改善的充分条件；其另行在记录窗口内逐轮比较二次损失所得的 `predicted_best_round_discrete` 不受此标签问题影响。[既有测试](../tests/test_m0_2_family.py) 中 `Q_f*S_n>Q*S_f ⇔ z*<1` 的检查仍正确，六个冻结案例的数值观察也不因此改变；它们不承担一般整数轮定理。

E2 不修改 [noise_floor_l2](../src/sand_m0/model.py) 的均匀样本公式 `(K−1)/(KN)`，该函数已明确只讨论 iid 均匀样本。以上三项均不进入 v0.4.0 主定理的证明链；定向复验不扩展模型、参数范围或成果声明。

E3 的通用上界写作 `w≤C<1`：原定理允许 `x=p_m=0`，此时 `w=C−x=C`。既有勘误的 `w<C` 简写在这个边界应放宽为 `w≤C`；所需 `w<1` 及斜率结论不变，定向检查包含该零质量边界。

## 当前与历史材料 / Current and historical material

canonical XO 原 S3 的[正式证明](S3_XO_BOUNDARY_THEOREM_PROOF.md)和[登记](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)保持原范围，S3 为 **CLOSED — PASS**。v0.2.0 的 A/T6 分别见[补充登记](XO_EXTENSION_REGISTRATION.md)与[独立审查记录](XO_EXTENSION_REVIEW.md)。[S1](S4_W1_S1_XO_PRECISION.md)、[S2](S4_W1_S2_BOUNDARY_MICROSCOPE.md)及更早 checkpoint 保留各自形成时的实验范围、术语与开放问题；旧规划不构成当前待执行任务。

The original S3 [proof](S3_XO_BOUNDARY_THEOREM_PROOF.md) and [registration](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md) retain their scope; S3 remains **CLOSED — PASS**. The v0.2.0 A/T6 supplements have separate [registration](XO_EXTENSION_REGISTRATION.md) and [review provenance](XO_EXTENSION_REVIEW.md). S1, S2, and dated checkpoints retain their original experimental scope and historical questions; their old plans are not active instructions.

原 S3 固定 canonical 初态、XO、α=1/4、最小指标平局、精确算术及 K≥6。补充 A 将长期结论推广到 0<α<1；补充 T6 在原 α=1/4、H=10K 下证明所有 K≥6 的冻结分类可靠性，只依赖原 S3。A 不把这个有限窗保证推广到其他 α。边界闭包不恢复完整质量向量，所有定理均不宣称周期性。

Original S3 fixes canonical XO, α=1/4, smallest-index ties, exact arithmetic, and K≥6. A extends the long-term result to 0<α<1; T6 proves universal frozen-classifier reliability only at α=1/4 and H=10K, depending on S3 rather than A. Scalar closure does not reconstruct the full mass vector, and no periodicity is claimed.


## 新补充的证据边界 / Supplement evidence boundaries

[验收入口](XO_EXTENSION_VALIDATION.md)区分精确完整状态、进入后的精确标量闭包和有限浮点标量代理。S12 的 200000 步统计只推进标量 x，d 由闭包推算；完整质量向量的 33 步 exact 核对另列。C1 为径向聚合/半共轭，参数及角向被扫次数公式见[正式附录](XO_RADIAL_AGGREGATION.md)。

The validation record separates exact full-state evolution, exact scalar closure after entry, and finite scalar float proxies. S12 does not measure a complete mass vector over 200000 steps; its exact 33-state full-vector check is separate. C1 is radial aggregation/semiconjugacy within its explicitly stated parameter and stopping scope.

## 证据名称 / Evidence names

| 名称 / Name | 阅读方式 / Meaning |
|---|---|
| Existing Math；New Mathematical Result | 分别指既有证明结果与新入库证明结果；范围由各自陈述和证明限定。 / Previously proved results and newly deposited proved results, within their stated assumptions. |
| Independent Audit | 已记录的独立审计结论与来源。原 S3 的 PASS 由 Owner 提供，完整审计对话未发布；A/T6 的工具协助独立技术审查另行登记；它不是期刊同行评审声明。 / The original S3 verdict was supplied by the Owner; A/T6 have separate tool-assisted technical-review provenance. It is not a journal peer-review claim. |
| VCR — Verified Computational Result | 指定代码、输入、算术方式及运行时限内的计算核验。旧文使用 Verified Simulation Result、Verified Numerical Analysis 等名称；这些计算证据没有因此变成一般定理。 / Computational evidence for specified code, inputs, arithmetic, and horizon. Older names such as Verified Simulation Result and Verified Numerical Analysis describe computational evidence, not universal proofs. |
| Working Hypothesis；Model-dependent Observation | 分别为待检验猜想与限定模型、数据范围内的观察。 / A hypothesis awaiting evaluation, or an observation within a specified model and dataset. |

## 易混淆的术语 / Terms that need care

| 历史词或标签 / Historical term or label | 准确含义 / Precise meaning |
|---|---|
| exact stop | 历史实验名称；实际为容差停机 $`d_t\le10^{-12}`$，不能据此宣称精确等于零。 / A historical experiment name for tolerance stopping, not a proof that discrepancy becomes exactly zero. |
| LOCK / CYCLE | S2 冻结有限窗标签。T6 仅在 α=1/4、H=10K 保证全 K 一致；LOCK 的尾长条件优先于 CYCLE。其他 α 下标签可能不符合长期二分。 / Finite-window labels; T6 guarantees all-K reliability only at α=1/4, H=10K. LOCK has precedence; other α values may disagree with the long-term dichotomy. |
| “周期-2”旧表述 / Old period-2 wording | 已由 2026-09-30 reconciliation 取代。正确描述为两个边界之间的切换，允许同支连续驻留；S3 无周期性结论。 / Superseded by two-boundary switching with same-branch dwell; S3 proves no periodicity. |
| exact tie（轮变化指标 / round-change metric） | 通常指 $`d_{t+1}=d_t`$。它不等于选取指标发生平局。 / Equality of successive discrepancy values, not necessarily a tie between selection candidates. |
| selection tie（边界区 / boundary regime） | $`x_t=h`$ 时两个边界超额相等，按规则选较小指标。 / At threshold equality the two boundary discrepancies tie and the smaller index is selected. |
| finite-K | 有限维；质量状态仍取连续实数值。 / Finite-dimensional, with continuously valued mass coordinates. |

冻结配置和原始输出中的历史名称保留。S1 配置字段 `exact_tolerence_fraction` 的拼写也是冻结来源的一部分；其值为 1/10^12。该字段与结果 metadata 保持原样，以便核对原始配置和输出。

Historical names in frozen configurations and raw outputs are preserved. The S1 key `exact_tolerence_fraction` is a legacy spelling for the exact tolerance 1/10^12; the frozen configuration and metadata retain it for reproducibility.

## 记号按文档局部定义 / Notation is local to each document

| 记号 / Symbol | 定义与出处 / Definition and source |
|---|---|
| $`D_{t,j}`$ | 当前状态的累计超额；不是喷撒律的失配。 / Cumulative excess of the current state, not mismatch of the redistribution law. |
| $`D_g(j)`$；$`\Delta(j)`$ | M1C2 与后续单步文稿的同一失配量 $`G(j)-T(j)`$，以各文目标 T 为准。 / Two names for the same law-minus-target cumulative mismatch, using the document's target T. |
| h（XO 证明 / XO proof） | 网格份额 1/K。 / Uniform bin share 1/K. |
| h(x)（M1C5） | 残差辅助函数 $`x[F_a(2-x)-1]`$；不是 XO 的网格份额。 / Auxiliary residual function, distinct from the XO grid share. |
| a（M1C5） | 被选前缀的位置 j*/K。 / Position of the selected prefix. |
| a（XO 证明 / XO proof） | 守恒的边界超额 b(1−b)。 / Conserved boundary discrepancy. |
| q（XO 边界证明 / XO boundary proof） | 第 m 格的喷撒质量份额；不是完整分布向量。 / Redistribution share of bin m, not the entire distribution vector. |

## M1C5 文稿更正 / M1C5 text corrections

更正位置为 [M1C5](M1C5_BOUNDARY_GATE_THEORY.md) §0 第 3 条、§5 以及 §4 端点说明。旧句“h 恒为零”应为：F_a≤1/2 时 h≤0，连续区间最大值 h(0)=0，正格点上 h<0。例如 F_a=1/2、x=1/2 时，h(x)=−1/8。

The corrected locations are M1C5 §0 item 3, §5, and the endpoint discussion in §4. The former assertion that h vanishes identically was incorrect. For F_a≤1/2, h is nonpositive, its continuous maximum is h(0)=0, and its values at positive grid points are strictly negative. At F_a=x=1/2, h(x)=−1/8.

§4 的端点分支现在用直接展开的恒等式，D=D_max≥0、0≤a≤1：

$$
h(a)-D=-(a+D)(1-a)^2\le0.
$$

This identity replaces the inaccurate intermediate argument in the endpoint branch. Both corrections preserve the existing necessary-gate conclusion and its scope. They do not alter any experiment or reopen canonical XO S3.

## 来源提交与公开版本 / Source commits and public versions

早期文稿或实验 metadata 中的某些提交 hash 属于私人研究仓库，是来源标识，不一定能在公开仓库解析。公开快照有独立 Git 历史；请用公开版本标签、公开提交和本仓库文件路径引用本版内容。私人日志、评审草稿和聊天记录不是公开证明的必要依赖。

Some commit hashes in dated documents or experiment metadata identify the private source history and cannot be resolved in this repository. The public snapshot has an independent history. Cite a public version tag or commit and the deposited files; private logs, review drafts, and chat records are not dependencies of the public proof.

版本与引用信息见[公开 Releases](https://github.com/XuZhilin2007/sand-redistribution-study/releases)。Apache LICENSE 附录中的版权占位符是官方应用示例；项目版权与文件许可范围另行说明，不能把示例当作项目漏填的声明。

See the [public releases](https://github.com/XuZhilin2007/sand-redistribution-study/releases) for versioned citation. Copyright placeholders in the Apache LICENSE appendix are official application examples; project attribution and licensing are specified separately.
