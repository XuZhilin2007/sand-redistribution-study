# 沙子再分配研究（Sand Redistribution Study）

## v0.4.0：统一参数区间上的全 K 两轮最优停止（2026-10-24，本地成果冻结）

在 canonical XO、精确算术和最小指标平局规则下，全部整数 K≥6、
全部实数 α∈[1999/8000,2001/8000]、全部整数 t≥0 均满足
U(p_t)≥U(p_2)。本版将 α=1/4 的结果扩展到与 K 无关的双侧参数区间。
完整证明由前缀、D₂ 进入、定理 A 驻留、W0、统一尾部及有限补集组成。

当前成果入口：[v0.4.0 登记](docs/V040_REGISTRATION.md)、
[验收与复验](docs/V040_VALIDATION.md)、[固定文件清单](docs/V040_MANIFEST.json)。
本地研究成果已收敛，公开发布仍以已有 v0.3.0 Release 为准。
下方各版本与阶段记录保留形成时的范围。

## v0.3.0：全 K 两轮最优停止与周期结构（2026-10-21）

在 canonical XO、α=1/4、精确算术和最小指标平局规则下，全部整数 K≥6、
全部整数 t≥0 满足 U(p_t)≥U(p_2)，其中 U=KΣ(p_i−1/K)²。
本版同时收录 T7 的最终精确选择字与渐近周期状态、B/B-II 初态结果、
一般 α 分类窗口与词典、平台和端点证书。完整范围见
[v0.3.0 成果登记](docs/V030_REGISTRATION.md)，实际验收见
[发布验收](docs/V030_VALIDATION.md)，文献定位见
[继承说明](docs/V030_RELATED_WORK.md)。

这是模型内的解析与计算辅助证明，不是实际沙坑的最优操作次数。
公开发布沿用[独立快照仓库](https://github.com/XuZhilin2007/sand-redistribution-study)，
私有过程材料与完整研究历史继续保留在本仓库。
正式公开版本已发布：[v0.3.0 Release](https://github.com/XuZhilin2007/sand-redistribution-study/releases/tag/v0.3.0)。
私有完整历史见[归档 Release](https://github.com/XuZhilin2007/sand-redistribution-study-private/releases/tag/v0.3.0-history)
的 Git bundle 附件，保留全部原提交 ID。因历史证书超过 GitHub 单文件限制，
私有远端 main 保留原基线；本机 main 为最新收尾树。恢复完整 v0.3.0 应下载 bundle，
不要把归档索引的自动源码包当作完整研究树。后续研究从最新本机 main 另建分支。

下方 v0.2.0 及 Stage 1–4 的“当前”“尚未证明”“未启动”等描述按形成时的
历史范围读取；本版新增结论以以上登记为准。

## v0.2.0 历史补充成果（2026-10-07）

A 与 T6 已分别正式登记。A 覆盖 canonical 初态、XO、最小指标平局、精确算术、K≥6 与 0<α<1 的长期边界结论；T6 只保证 α=1/4、H=10K 的原冻结分类器。两者范围分开。
当前入口为[补充成果登记](docs/XO_EXTENSION_REGISTRATION.md)、[独立审查记录](docs/XO_EXTENSION_REVIEW.md)、[候选验收](docs/XO_EXTENSION_VALIDATION.md)及[辅助附录](docs/XO_SUPPLEMENTARY_RESULTS.md)。原 S3、S1/S2 冻结材料和以下历史路线按原范围保留。
实验室现行稿与旧稿保存在 `archive/extension_lab_2026-10-07/`；完整独立审查记录保存在 `release_review_2026-10-07/`，均不进入公开仓库。公开仓库保持独立提交历史。版本发布对应关系另见 `release_records/v0.2.0/`（完成发布时登记）。

```text
python -B tests/count_tests.py
python -B experiments/xo_extensions/validate.py --scope full --out-dir ../xo-v020-run
```

有限验证、标量代理与数学证明的证据等级按各自登记读取；新版本不宣称周期性、完整状态渐近、一般噪声/二维理论或文献新颖性。下方原阶段记录中的“未启动”不否定本段已登记的 A/T6 及限定附录。


> 当前阶段：Stage 4 — Canonical XO S3 CLOSED — PASS；证明已入库，独立审计 PASS（2026-09-30）
>
> 文档维护日期：2026-10-06；基线整理日期：2026-09-13
>
> 仓库性质：私有研究主仓库；保留完整历史与内部材料；公开发布使用独立快照仓库

## 项目是什么

本项目源于一次真实的沙坑操作：人在粗糙水泥地附近用两根棍子把沙向外泼出，观察落沙分布，再选择性地把明显需要处理的沙扫回身前并重新泼出。

一次普通泼沙通常呈扇形，近处较密、远处较稀。操作者并不会机械地无限重复，而是根据当前沙面作出反馈：近处太密时回收近处；中远处出现明显密集块时，可能从该处一路向身体方向扫回；当近处已可接受且中远处没有明显密集块时停止。

本项目研究的是这种带观察、选择性回收和停止判断的随机再分配过程。

目前正式研究采用由该问题启发的一维确定性质量反馈模型；真实泼沙中的随机性尚未经过实测校准。

## 原始问题

本项目的核心问题不是“无限迭代后是否收敛到均匀状态”，而是：

> 在“观察当前分布 → 回收明显需要处理的沙 → 重新泼出”的有限反馈操作中，是否会出现一个最佳或近最佳的均匀状态？什么样的反馈策略与停止规则能够到达并识别它？

这里的“均匀”必须与明确的空间尺度和评价指标一起定义，不能直接等同于肉眼观感，也不能默认存在唯一的现实最优状态。

## 当前研究阶段

本项目源于一次真实沙坑操作（观察沙面 → 扫回过密区域 → 重新泼出 → 再观察 → 适时停止），研究这种带反馈的随机再分配过程中有限时间最佳状态的出现与停止条件。

**Stage 1、Stage 2、Stage 3 均已闭合，Stage 4 的 canonical XO S3 已 CLOSED — PASS。** 当前 XO 数学结论以[入库证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)和[定理登记](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)为准；各阶段 checkpoint 保留为历史节点。

- **Stage 1（M0 → M1B.2，2026-09-17 至 09-19）**：canonical baseline → 解析机制 → 一般 fixed-q 理论 → 自适应反馈（M1A）→ 长时程诊断（M1A.1）→ 容差停止（M1B）→ 指标审计（M1B.1）→ 机制证伪（M1B.2）。核心现象：state-dependent feedback 下进入 substantial-redistribution churn regime，数值容差停机（d ≤ 10⁻¹²，历史名称 exact stop）来自低频 annihilable configuration（tested-range stopping time ∝K），finite deadband 恢复 tested-range K 稳健性。
- **Stage 2（M1C → M1C.6，2026-09-19 至 09-20，Feedback Mismatch Mechanism）**：对照实验（oracle B 与 fixed-uniform U 均无 churn → churn 关键驱动是 near-bias/target mismatch 而非 selection 或 deficit-unawareness）→ 一般理论（一步分解 `D_new = D_TM_new + M·Δ`、target-matched 与 G≤T 单调性定理、精确单步 rebound 判据——canonical A 轨迹 5344 轮 FP=FN=0）→ canonical 解析门 `a_crit(D)`（1958 反弹轮 FN=0，tracks 观察门 ≈0.58）→ churn 动力学闭合（gate exit/return cycle：never-return=0/1958 episodes、gap-closing=1.0、switch 主机制为 mismatch 重塑 profile）。

**Stage 3 / M2 — Redistribution-Mismatch Generality 已闭合**（2026-09-20，[closure checkpoint](docs/RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md)；[M2 Gate 1](docs/M2_GATE1_GENERALITY_GATE.md) / [M2 Gate 2](docs/M2_GATE2_SIGN_PLACEMENT.md)）。五层结论：general one-step theory 在 7 个 kernels 上计算核验；positive mismatch alone 不组织 tested long churn（F9）；peak/integrated 失配强度不是充分统计量（F10）；静态 canonical witness geometry 不可迁移（F11——XA 零反弹 t=6 vs XO 全 K 10K censoring 的窄平台复发瞬态，witness geometry is endogenous）；存在多个定性不同的 recurrent regimes。External Review Pack v0.4 完成并已获两位真人 reviewer 第一轮反馈（2026-09-29 登记）。

**Stage 4 / canonical XO**：Wave 1 S1（[XO 高精度重放](docs/S4_W1_S1_XO_PRECISION.md)）与 S2（[边界切换显微镜](docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md)）均已完成。S1 修正 XO rebound 的浮点口径；S2 给出 51/51 tested-range mod-6 分类和逐步 Fraction 边界坐标闭包。LOCK/CYCLE 保持冻结有限 horizon 分类标签。**S3 CLOSED — PASS**；[canonical XO 证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)已入库，[定理登记](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)为 **REGISTERED — PROOF DEPOSITED / INDEPENDENT AUDIT PASS**；[独立审计 provenance](docs/S3_XO_BOUNDARY_THEOREM_AUDIT.md)如实记录 Owner 提供的 GPT-6.1 Sol PASS 结论。T2–T5 证明有限边界进入、持久性、边界标量闭包与无限时间 residue dichotomy；T1 保持 [Existing Math]，S1/S2 保持 [VCR]。**当前优先维护公开快照；S8、S9 及其他研究分支暂不启动。**

所有结论只对各自声明的一维 toy model 成立；canonical XO 已有入库并获独立审计 PASS 的数学证明，不构成一般沙子再分配定理或现实物理结论；novelty 未建立。

## 信息等级

后续文档、实验记录和结论必须使用下列等级。标签描述的是证据性质，不描述观点听起来是否合理。

| 中文标签 | 英文标签 | 使用条件 |
|---|---|---|
| **真实观察** | Real-world Observation | Owner 在真实沙坑或水泥地操作中实际观察到或实际执行过的事实。 |
| **模型假设** | Model Assumption | 为使模拟可运行而人为加入的规则、参数、几何或测量定义。 |
| **模型推导数学结果** | Model-derived Mathematical Result | 从版本化 canonical model 的明确假设和更新规则中严格推导，并通过独立数值实现核验的数学结论。只对声明的模型成立；"数学上精确"不等于真实世界正确；不自动构成文献支持结果、novelty，也不自动推广到其他模型族。（2026-09-18 由 Owner + GPT 批准新增） |
| **已验证计算结果** | Verified Computational Result（VCR） | 已在明确代码版本、完整输入、算术方式、运行时限和评价方法下实际复现的结果；只对声明范围成立。旧文的 Verified Simulation Result / Verified Numerical Analysis 保留其历史名称，不因此升级为一般定理。 |
| **工作猜想** | Working Hypothesis | 由观察或多组探索现象提出，但尚未证明、尚未充分验证的研究猜想。 |
| **文献支持结果** | Literature-supported Result | 经正式文献检索并由可追溯来源支持的结论。当前仓库尚未建立此类结果。 |

另设历史状态标签 **探索性模拟发现（Exploratory Simulation Finding）**：用于描述早期聊天探索中反复出现的现象。该标签不自动满足“已验证计算结果”的复现要求。2026-09-13 初始 checkpoint 中的探索性发现属于这一历史状态；后续已复现计算和已入库证明按各自证据等级登记。

详细治理规则见 [研究 checkpoint](docs/RESEARCH_CHECKPOINT_2026-09-13.md) 和 [模型假设台账](docs/MODEL_ASSUMPTIONS.md)。当前术语、局部记号与勘误见[中英双语阅读说明](docs/READING_NOTES.md)。

## 当前探索性认识

以下内容全部是**探索性模拟发现**，不是理论定理，也不是现实物理结论：

1. 多种简化模型的初始分布都表现为近密远疏和明显不均匀。
2. 选择性回收并重新泼出后，模拟中的空间均匀度通常能够明显改善。
3. 多种简化模型会经历有限的最佳或近最佳阶段。
4. 如果无视状态继续操作，常出现近处被掏空、稀疏沟或远处堆积等过度修正（over-correction）。
5. “没有明显值得修正的区域就停止”可能是核心机制之一。
6. 人的视觉容忍度可能具有功能意义；对任意微小偏差持续修正反而容易过度操作。

这些描述只概括历史探索的共同现象。具体轮数、阈值、扫线和参数值不能跨模型比较，更不能外推为真实世界最优值。

## 最关键的未知量

真实的单次泼沙分布（throw distribution）尚未测量。

早期模型用过 `r = U^2` 来制造“近处更容易落、远处较少”的径向分布。它是**高影响临时假设（high-impact provisional assumption）**，没有真实实验数据支持，并会直接影响：

- 初始近密远疏的强度；
- 每轮质量向外转移的速度；
- 最佳或近最佳阶段出现的时机；
- 模拟中的最终空间结构。

因此，后续应先比较一族候选泼沙分布，再决定是否值得采集真实照片做粗略校准。详见 [后续步骤](docs/NEXT_STEPS.md)。

## 仓库结构

```text
sand-redistribution-study/
├── README.md
├── .gitattributes
├── .gitignore
├── docs/
│   ├── RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md
│   ├── RESEARCH_CHECKPOINT_2026-09-20.md
│   ├── RESEARCH_CHECKPOINT_2026-09-19.md
│   ├── RESEARCH_CHECKPOINT_2026-09-13.md
│   ├── RESEARCH_SCOPE.md
│   ├── MODEL_ASSUMPTIONS.md
│   ├── M0_MODEL.md
│   ├── M0_MECHANISM.md
│   ├── M0_GENERAL_Q.md
│   ├── M1A_MODEL.md
│   ├── M1A_LONG_HORIZON.md
│   ├── M1B_TOLERANCE.md
│   ├── M1B_METRIC_SCALE_AUDIT.md
│   ├── M1B_MICRO_CORRECTION.md
│   ├── M1C_CORRECTION_KERNEL_ABLATION.md
│   ├── M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md
│   ├── M1C2_TARGET_MATCHED_KERNEL_THEORY.md
│   ├── M1C3_REDISTRIBUTION_LAW_ORDER.md
│   ├── M1C4_CHURN_REGIME_AUDIT.md
│   ├── M1C5_BOUNDARY_GATE_THEORY.md
│   ├── M1C6_BOUNDARY_RETURN.md
│   ├── M2_GATE1_GENERALITY_GATE.md
│   ├── M2_GATE2_SIGN_PLACEMENT.md
│   ├── STAGE4_RESEARCH_MAP.md
│   ├── S4_W1_S1_XO_PRECISION.md
│   ├── S4_W1_S2_BOUNDARY_MICROSCOPE.md
│   ├── RESEARCH_CHECKPOINT_2026-09-29_STAGE4_POST_S2.md
│   ├── READING_NOTES.md
│   ├── POST_CLOSURE_MATHEMATICAL_REVIEW_NOTES.md
│   ├── ASTRA_EXTERNAL_REVIEW_READING_MAP.md
│   ├── ERRATUM_2026-09-18_M1A_DMAX.md
│   ├── S3_XO_BOUNDARY_THEOREM_PROOF.md
│   ├── S3_XO_BOUNDARY_THEOREM_AUDIT.md
│   ├── S3_XO_BOUNDARY_THEOREM_REGISTRATION.md
│   ├── NEXT_STEPS.md
│   └── RESEARCH_LOG.md
├── src/
│   └── sand_m0/            # canonical simulation source（M0/M1A 类模型；2026-09-17 起）
├── tests/
│   ├── test_m0.py
│   ├── test_m0_mechanism.py
│   ├── test_m0_2_family.py
│   ├── test_m1a_adaptive.py
│   ├── test_m1a_1_diagnostics.py
│   ├── test_m1b_tolerance.py
│   ├── test_m1b_1_metric_scale.py
│   ├── test_m1b_2_micro_scaling.py
│   ├── test_m1c_kernel_ablation.py
│   ├── test_m1c1_uniform_kernel.py
│   ├── test_m1c2_kernel_theory.py
│   ├── test_m1c3_order_criterion.py
│   ├── test_m1c4_churn_dynamics.py
│   ├── test_m1c5_boundary_gate.py
│   ├── test_m1c6_boundary_return.py
│   ├── test_m2_gate1_generality.py
│   ├── test_m2_gate1_audit.py
│   ├── test_m2_gate2_sign_placement.py
│   ├── test_s4_s1_precision.py
│   └── test_s4_s2_microscope.py
├── experiments/
│   ├── m0_baseline/        # canonical M0 正式实验（config / runner / results / 报告）
│   ├── m0_1_mechanism/     # M0.1 机制分析：alpha 时间尺度实验（2026-09-18）
│   ├── m0_2_distribution_family/  # M0.2 Phase 2A：分布族与运行前预测（2026-09-18）
│   ├── m1a_adaptive_boundary/     # M1A：自适应扫掠边界，最小反馈模型（2026-09-18）
│   ├── m1a_1_long_horizon/        # M1A.1：长时程与分辨率诊断（2026-09-18）
│   ├── m1b_tolerance/             # M1B：感知容差 / deadband 停止（2026-09-18）
│   ├── m1b_1_metric_scale/        # M1B.1：指标与观察尺度审计（2026-09-18）
│   ├── m1b_2_micro_scaling/       # M1B.2：late-stage 微修正标度诊断（2026-09-19）
│   ├── m1c_kernel_ablation/       # M1C：correction-kernel ablation A vs B（2026-09-19）
│   ├── m1c1_uniform_kernel_ablation/  # M1C.1：fixed-uniform kernel ablation A/U/B（2026-09-19）
│   ├── m1c4_churn_dynamics/       # M1C.4：churn-regime 动力学审计（2026-09-19）
│   ├── m1c5_boundary_gate/        # M1C.5：边界门控反弹理论验证（2026-09-19）
│   ├── m1c6_boundary_return/      # M1C.6：boundary return 与峰竞争审计（2026-09-19）
│   ├── m2_gate1_generality/       # M2 Gate 1：mismatch sign/amplitude/shape 泛性对照（2026-09-20）
│   ├── m2_gate2_sign_placement/   # M2 Gate 2：sign placement 符号翻转对照 XA/XO（2026-09-20）
│   ├── s4_w1_s1_xo_precision/     # S4 Wave 1 S1：XO 高精度重放（exact/mp/float 四层，2026-09-29）
│   ├── s4_w1_s2_boundary_microscope/  # S4 Wave 1 S2：XO 边界切换显微镜（mod-6 分类 51/51 + 边界标量闭包精确成立，2026-09-29）
│   └── raw_chat_runs/
│       ├── README.md
│       └── [12 组历史探索目录]
└── archive/
    ├── README.md
    └── checkpoint_2026-09-13/
        └── sand_redistribution_checkpoint_2026-09-13/
            └── [原始 checkpoint 内容]
```

`src/` 与 `experiments/m0_baseline/` 自 2026-09-17 起存放首个 canonical model M0 及其正式实验；详见 [M0 模型定义](docs/M0_MODEL.md) 与 [M0 实验报告](experiments/m0_baseline/M0_REPORT.md)。`experiments/raw_chat_runs/` 仍为只读历史证据，两者严格分离。

## 文档导航

- [阅读说明、术语与勘误（中英双语）](docs/READING_NOTES.md)：当前结论与历史材料、VCR 名称、容差停机、边界切换、局部记号和 M1C5 文稿更正；公开版本见[独立快照仓库](https://github.com/XuZhilin2007/sand-redistribution-study)。
- [Stage 4 Research Map v1.0](docs/STAGE4_RESEARCH_MAP.md)：**Stage 4 历史规划（PROPOSED；S3 已 CLOSED — PASS）**——开放问题矩阵、S1–S9 branch + parked 清单、三波次预算（≤12 轮进 Final Synthesis）、Owner+GPT 决策门 G-S4-A…F（2026-09-29）。
- [S4 Wave 1 S1 报告](docs/S4_W1_S1_XO_PRECISION.md)：**XO 高精度重放（Wave 1 首轮，已完成）**——exact rational tier + mpmath ladder 四层对照；XO d/selection/stopping 为真实算术性质（理论解禁）；exact N_R = 0/398/0/1598 与锁定平台/两边界切换分类（H 内）；float rebound 口径修正（2026-09-29）。
- [S4 Wave 1 S2 报告](docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md)：**XO 边界切换显微镜（Wave 1 第二轮，已完成）**——冻结 config 全量执行：mod-6 分类预测 51/51 全中（LOCK ⟺ K mod 6 ∈ {0,1,2}，tested K=6..53,100,200,400，H 内冻结分类器）；边界坐标标量闭包——p_m 递推与 full exact 轨道的边界坐标逐步 Fraction 恒等（C1–C7 全过、零 scientific events）；三类 residue 切换指纹（2026-09-29）。
- [Post-S2 研究状态简报](docs/RESEARCH_CHECKPOINT_2026-09-29_STAGE4_POST_S2.md)：2026-09-29 历史 checkpoint；已加 S3 proof-deposit 状态注，原 S1/S2 [VCR] 范围保持。
- [S3 定理登记](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)：**REGISTERED — PROOF DEPOSITED / INDEPENDENT AUDIT PASS**；T1–T5 与证据分层。
- [S3 canonical XO 证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)：P1–P4 完整证明，含 K=6、threshold equality 与 t* ≤ ⌊81K/20⌋+1；仅边界闭包，无周期性结论。
- [S3 独立审计记录](docs/S3_XO_BOUNDARY_THEOREM_AUDIT.md)：GPT-6.1 Sol，PASS — no substantive gap found；Owner 提供的审计 provenance。
- [阶段 checkpoint 2026-09-20 Stage 3](docs/RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md)：**Stage 3 closure（历史冻结节点）**——五层结论、F9–F11、resolution/exact-tie concerns、当时的 External Review 问题清单 Q1–Q10；canonical XO 的后续证明状态以 S3 入库证明为准。
- [阶段 checkpoint 2026-09-20](docs/RESEARCH_CHECKPOINT_2026-09-20.md)：Stage 2 closure（历史冻结节点）——研究逻辑、定理/observation 分层、falsified 清单。
- [阶段 checkpoint 2026-09-19](docs/RESEARCH_CHECKPOINT_2026-09-19.md)：Stage 1 closure（历史冻结节点）。
- [研究 checkpoint](docs/RESEARCH_CHECKPOINT_2026-09-13.md)：问题来源、历史探索、当前认识、限制和交接状态。
- [研究范围](docs/RESEARCH_SCOPE.md)：现在研究什么、暂时不研究什么。
- [模型假设台账](docs/MODEL_ASSUMPTIONS.md)：集中记录临时假设、依据、可信度和校准需求（含 M0 新增假设）。
- [M0 模型定义](docs/M0_MODEL.md)：首个 canonical model 的正式定义与证据等级边界。
- [M0 机制解释](docs/M0_MECHANISM.md)：为什么第 4–5 轮最佳、alpha 为何只是时间尺度参数（直觉 + 数学附录）。
- [一般 q 机制分析](docs/M0_GENERAL_Q.md)：任意 throw distribution 下的结构与 interior optimum 判据（Phase 2A 理论）。
- [M0 实验报告](experiments/m0_baseline/M0_REPORT.md)：EXP-M0-A/B/C 的可复现结果（2026-09-17）。
- [M0.1 机制实验报告](experiments/m0_1_mechanism/M0_1_REPORT.md)：correction-strength 时间尺度实验（2026-09-18）。
- [M0.2 分布族报告](experiments/m0_2_distribution_family/M0_2_REPORT.md)：冻结幂族、运行前预测与对照（2026-09-18）。
- [M1A 模型定义](docs/M1A_MODEL.md)：自适应扫掠边界——首个最小状态依赖反馈模型（模型假设）。
- [M1A 实验报告](experiments/m1a_adaptive_boundary/M1A_REPORT.md)：M0 vs M1A、边界轨迹、自然停止与 alpha 结构检查（2026-09-18）。
- [M1A.1 长时程诊断](docs/M1A_LONG_HORIZON.md)：数值容差停机在 T=5000 内触发、tested-range 停止时间 ∝K、边界歧义细化与 one-sided KS 术语注记（2026-09-18）。
- [M1B 模型定义](docs/M1B_TOLERANCE.md)：感知容差 / deadband 停止——只改停止规则的最小扩展。
- [M1B.1 指标尺度审计](docs/M1B_METRIC_SCALE_AUDIT.md)：高 regret 是小分母伪影；retention 语义下存在宽 near-optimal region（2026-09-18）。
- [M1B.2 微修正标度诊断](docs/M1B_MICRO_CORRECTION.md)：micro-correction 假说被证伪；数值容差停机对应稀有单轮消灭事件，tested-range 尾巴 ∝K 在等待时间（2026-09-19）。
- [M1C correction-kernel ablation](docs/M1C_CORRECTION_KERNEL_ABLATION.md)：selection 不变、只换 redistribution kernel 的单机制对照——churn 完全消失，B 三轮精确均匀收敛（2026-09-19）。
- [M1C 实验报告](experiments/m1c_kernel_ablation/M1C_REPORT.md)：A vs B 可复现结果与 baseline reproduction 门禁（2026-09-19）。
- [M1C.1 fixed-uniform kernel ablation](docs/M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md)：中间 control Variant U（fixed uniform、仍 deficit-unaware）同样无 churn——churn 关键驱动定位到 near-bias/target mismatch（2026-09-19）。
- [M1C.1 实验报告](experiments/m1c1_uniform_kernel_ablation/M1C1_REPORT.md)：A/U/B 可复现结果与 baseline integrity 门禁（2026-09-19）。
- [M1C.2 target-matched kernel theory](docs/M1C2_TARGET_MATCHED_KERNEL_THEORY.md)：一步更新严格分解为收缩的 target-matched 分量 + 失配项 M·D_g；G=T ⟹ D_max 严格单调下降定理（2026-09-19）。
- [M1C.3 redistribution-law order](docs/M1C3_REDISTRIBUTION_LAW_ORDER.md)：G≤T 序定理 + 精确单步反弹判据 M·Δ>γ——A 轨迹 5344 轮 iff 验证零误分（2026-09-19）。
- [M1C.4 churn-regime audit](docs/M1C4_CHURN_REGIME_AUDIT.md)：R_t 跨 K 稳定带（0.36=P(R>1)）、边界位置门控反弹、C-C-R 切换节奏、终端温和低谷（2026-09-19）。
- [M1C.5 boundary-gate theory](docs/M1C5_BOUNDARY_GATE_THEORY.md)：反弹必要门槛解析闭式 a_crit(D)=[(1−D)+√(D(D+2))]/2——1958 反弹轮 FN=0，门槛 tracks 观察门 ≈0.58（2026-09-19）。
- [M1C.6 boundary return](docs/M1C6_BOUNDARY_RETURN.md)：churn 反馈回路闭合——gate exit/return cycle（never-return=0）、峰竞争 gap 收缩（share=1.0）、switch 主机制为 reinjection 重塑（2026-09-19）。
- [M2 Gate 1 generality gate](docs/M2_GATE1_GENERALITY_GATE.md)：Stage 3 首实验——sign/amplitude/shape 三对照（K-/KW/KS）；单步判据 5510 轮零误分泛化，churn-like transient 未在 positive 对照出现；contract W0/S1（2026-09-20）。
- [M2 Gate 2 sign placement](docs/M2_GATE2_SIGN_PLACEMENT.md)：逐点符号反射对 XA/XO——XA 零反弹 t=6、XO 全 K 10K censoring 窄平台复发；witness geometry endogenous（X4 强证伪）（2026-09-20）。
- [Post-closure mathematical review notes](docs/POST_CLOSURE_MATHEMATICAL_REVIEW_NOTES.md)：Stage 3 closure 后的数学卫生登记——stopping=tolerance 术语、finite-K≠finite-state、M1C5 h≤0 勘误、XO near-marginal 软化、invariant-anchor analytic clue、Review Pack 措辞规则（2026-09-21）。
- [Expert Review Pack Draft v0.4](docs/external_review/EXPERT_REVIEW_PACK_DRAFT_V0.4.md)：保留于私有仓库的历史评审材料；**2026-09-29 已登记两位真人 reviewer 第一轮反馈，External Review Gate 已关闭**，不再等待旧稿所列的 final verification。v0.1–v0.3、配套 CLAIM_LEDGER / SOURCE_INDEX 和导出件作为历史来源保留；当前 canonical XO 结论以 S3 入库证明与登记为准。
- [后续步骤](docs/NEXT_STEPS.md)：从 canonical baseline 到可选真实实验的分阶段路线。
- [研究日志](docs/RESEARCH_LOG.md)：只追加的决策与进展记录。
- [历史实验说明](experiments/raw_chat_runs/README.md)：raw runs 的来源、边界与使用规则。
- [归档说明](archive/README.md)：原始 checkpoint 的完整性记录和作废文件说明。

## 禁止误读

不得把任一历史模拟中的以下数值解释为真实世界最优参数：

- “最佳第 4 轮”“最佳第 15 轮”或任何具体轮数；
- 某个固定扫线半径，例如 `R = 0.4`；
- 某个视觉或停止阈值，例如 `1.8`；
- 某个回收比例或 recovery efficiency，例如 `90%`。

70%–100% 回收效率在早期敏感性探索中对核心定性现象影响较小、更多改变速度；这仍只是特定简化模型下的探索性结果，不能据此宣称真实回收效率存在 90% 最优点。

当前也不得声称：已证明收敛、已证明普适的有限最优时机、已确认学术新颖性，或已有 publication-quality data。

## 历史 prompt 状态

旧版 `CODEX_BOOTSTRAP_PROMPT.md` 只保存在不可改写的历史 checkpoint 中，已明确作废。它要求立即实现 baseline simulator 的旧指令不再有效；当前研究推进应以本 README、`docs/` 下的治理文档以及 Owner 的最新明确指令为准。

## 给下一次研究窗口的最短交接路径

按以下顺序阅读即可恢复上下文：

1. 本 README；
2. [READING_NOTES.md](docs/READING_NOTES.md)（当前术语与勘误）；
3. [S3_XO_BOUNDARY_THEOREM_REGISTRATION.md](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)与[入库证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)（当前 canonical XO 数学结论）；
4. [Stage 3 checkpoint](docs/RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md)与[Stage 2 checkpoint](docs/RESEARCH_CHECKPOINT_2026-09-20.md)（历史证据与单步理论）；
5. [MODEL_ASSUMPTIONS.md](docs/MODEL_ASSUMPTIONS.md)；
6. [RESEARCH_SCOPE.md](docs/RESEARCH_SCOPE.md)与[NEXT_STEPS.md](docs/NEXT_STEPS.md)；
7. 最后在确有需要时查看 `experiments/raw_chat_runs/`，不要先从历史图片中的“最佳值”反推研究结论。
