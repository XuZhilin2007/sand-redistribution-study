# Stage 4 Post-S2 Research State Brief（研究状态简报）

> **当前状态注（2026-09-30）**：本文件主体为 Post-S2 历史快照。Canonical XO **S3 CLOSED — PASS**；[证明已入库](S3_XO_BOUNDARY_THEOREM_PROOF.md)，[独立审计 PASS](S3_XO_BOUNDARY_THEOREM_AUDIT.md)。主体中“约化未证明/无限时间行为开放”的当时状态，在 P1–P4 的精确 scope 内已由该证明取代；原 S1/S2 [VCR] 分类不变。

> 日期：2026-09-29。性质：**consolidation 文档**——S2 收口后的研究状态汇总，供 Owner 评审、未来数学协作者与外部 reviewer 阅读定位。不含新实验、新 kernel、新证明主张。
>
> 证据分级（全文严格使用）：**[Math]** = Mathematical result（既有推导/证明）；**[VCR]** = Verified computational result（ stated scope 内的精确计算结果）；**[Obs]** = Observation（数值/实验模式）；**[Hyp]** = Hypothesis（候选解释或未来定理目标）。四类不互相升级。
>
> 写作前核验：本简报的全部关键数字在写作当日对照 committed CSV/JSON 重新验证（脚本核验，见 §8）；发现并修正一处 prose 过度概括（S2 报告的 LOCK 尾部公式，见 §8c）。HEAD 于写作时为 `102c86a`。

---

## 0. 模型一分钟（自包含所需的最小定义）

K 个等宽 bin 上的质量向量 p（总和 1）。每轮：

1. 算累计超额 D_j =（前 j 个 bin 的质量）− j/K，j = 1..K；
2. 选择 s_t = 最小 argmax D_j（平局取小指标），d_t = D_{s_t}；
3. 若 d_t ≤ 10⁻¹² 停机；否则扫前缀 1..s_t：每 bin 失去 1/4 自身质量，全部移除质量按固定喷洒律 q 重撒到**所有** K 个 bin。

**XO 喷洒律**（Stage 3 Gate 2 冻结）：CDF G_O(x) = x/2（x ≤ 1/3）｜3x/2 − 1/3（≤ 2/3）｜2x − 2/3（≤ 5/6）｜1（x ≥ 5/6）——负失配覆盖历史 witness 带、正失配只在远端；断点 {1/3, 2/3, 5/6}。初态 = canonical q_near（CDF 2x − x²）。记 m = ⌈5K/6⌉、b = m/K、p_m = bin m 的质量、F_m = 前 m 个 bin 的总质量。

---

## 1. Executive Summary

**S2 改变了什么**：S1 证明 XO 的"长瞬态/反复反弹"现象学需要用精确算术重读（一半 tested 配置在精确算术下、H=10K 内**零反弹**，float 的 rebound 多数是舍入噪声）；S2 进一步把 XO 的**边界观测量**（selection 与 discrepancy 的演化）在 tested 范围内**精确闭包于一个标量有理递推**——单个变量 p_m（bin m 质量）围绕阈值 1/K 的两支仿射更新——并在 51 个 K 的 exact census 上验证了 **regime 分类由 K mod 6 决定**（51/51 应验，tested 范围、冻结有限 horizon 分类器）。

**当前的数学对象**：映射 φ_K：p_m ↦ p′_m（p_m > 1/K 时 p′ = ¾p + ¼F_m q_K；否则 p′ = p + ¼(F_m − p_m)q_K，其中 F_m = 2b − b²、q_K = (6 − (K mod 6))/(3K)）。LOCK/CYCLE 是**冻结的有限 horizon（H=10K）实验分类标签**，不是定理级渐近态；φ_K 的无限时间行为是开放问题。φ_K 是显式、可证伪、每步可 Fraction 验证的标量 switched affine 映射族——S3 的定理化目标（bounded proof attempt，见 §6）。

**范围限定**：以上为 tested-range（K = 6..53, 100, 200, 400；canonical 初态；XO kernel；H = 10K）内的 [VCR]。mod-6 规则**不是**一般定理；约化**不是**已被证明的一般命题。

## 2. Research Evolution（不夸大的版本）

1. **现实问题**（Stage 1）：沙坑反复"喷沙—扫沙"何时收敛、何时不收敛 → 抽象为 1D deterministic mass-feedback 模型（M0/M1A，finite-dimensional，非物理仿真）。
2. **Stage 2–3**：one-step 数学（分解恒等式、rebound 判据、boundary gate）+ 7-kernel 对照；发现两个反事实 law（XA/XO）行为截然不同；**XO 出现全 K censoring 的窄平台复发瞬态**（当时为 float 观测 [Obs]）；falsification 清单 F9–F11 排除了三类简单解释；XO selection-bound invariant 获得归纳证明（Pack v0.3/v0.4，[Math]）。reviewer 共同指向 switching/invariant 结构。
3. **S1（precision correction）**：exact rational 重放（可行——CDF 在有理网格上严格有理）分离了真动力学与浮点效应：d/selection/stopping 是真算术性质；**K=50/200 在 H 内零反弹（分类 LOCK）**（float 的 41/56 个"rebound"全是 tie 轮噪声），K=100/400 为两边界切换（分类 CYCLE）[VCR]。
4. **S2（boundary reduction discovery）**：设计评审的纸面推导（当时 [Hyp]）——F_m 守恒 + argmax 限制 ⟹ 边界动力学 = 一维 p_m 映射，且二分判据坍缩为 K mod 6——被冻结 census 在 51 个 K 上**全部精确证实** [VCR]；C1–C7 零失败、零 scientific event。

每一步的证据等级变化都被显式登记；没有把 [VCR] 写成 [Math]，也没有把 [Hyp] 写成结论。

## 3. Current Established Results

### 3a. Existing mathematical results（[Math]，Pack v0.3/v0.4 §6.2）

对 canonical 初态 + XO law、K ≥ 6（m = ⌈5K/6⌉、b = m/K）：

1. s_t ≤ m 对一切 t（selection 永不越过 m）；
2. D_{t,m} = b(1−b) 对一切 t（不变量坐标精确守恒）；
3. d_t ≥ b(1−b) ≥ 10/121 > 10⁻¹² 一致成立 ⟹ **tolerance stopping rule 在 exact real-arithmetic 模型中不可触发**。

边界：不证明 state 收敛性、不证明周期性、floating rebound 计数不在其内（S1 已另行裁定）。

### 3b. S1 verified computational results（[VCR]，tested K = 50/100/200/400，H = 10K）

1. **精确 rebound 计数**：N_R^exact = 0 / 398 / 0 / 1598（float：41 / 581 / 56 / 2380）。K=50/200 的 float rebound 现象整体是 exact-tie 轮上的舍入噪声伪造；K=100/400 的 recurrent 现象真实存在但 float 计数偏高 46%–49%。
2. **float boundary 旗标 = exact tie 的逐轮 1:1 检测器**（7500/7500 轮）——exact |worst| 严格双峰（0 vs ≥5.9×10⁻⁵），1e-12 诊断落在 ≥7 个数量级空隙中。
3. **锁定平台**（K=50/200）：selection 3 轮下降后锁定 s = m 至 horizon，d_t = b(1−b) **精确常值**（K=50: 84/625 = 0.1344；K=200: 0.137775），状态仍缓慢演化（非不动点）。
4. **两边界切换**（K=100/400，分类 CYCLE，H 内）：selection 限制于 {m−1, m}，**非严格交替**（同支驻留 dwell 1–3），exact 窄带 |e| ≤ 7.9×10⁻⁴ / 1.9×10⁻⁴。
5. float64 的 d_t 全程追踪 exact 真值至 ≤ 5.3×10⁻¹⁵；KC sanity 精确复现（t_stop = 324、N_R = 117 @ K=50）。

### 3c. S2 verified computational results（[VCR]，tested K = 6..53 ∪ {100, 200, 400}，H = 10K）

1. **边界坐标标量闭包精确成立**（C1–C7 全过，51/51）：post-entry 边际动力学（selection/discrepancy）与标量递推 φ_K 逐步 Fraction 恒等（full 状态向量继续多维演化，本结果限于边界观测量）；F_m = 2b − b² 逐状态守恒（C1）；argmax 限制零违例（C2）；分支方程、selection 阈值一致性、d 重构全部精确（C3–C5）；q_K 闭式 (6 − (K mod 6))/(3K) 与定义值全等（C7）。
2. **mod-6 分类 51/51**（冻结分类器，H=10K 内）：LOCK ⟺ K mod 6 ∈ {0,1,2}（25 K）；CYCLE ⟺ K mod 6 ∈ {3,4,5}（26 K）；零 OTHER、零 mismatch、零 scientific event。
3. **三类 residue 切换指纹**（由已验证递推精确生成，H 内；selection 均限制于 {m−1, m}，非严格交替）：mod-3 慢切换（mean rounds per rebound ~10.2–12.6，A-dwell 长达 11）；mod-4 紧凑切换（~2.5，maxA=2/maxB=3）；mod-5 B-重（~3.3，B-dwell 长达 8）。
4. LOCK K 的尾部游程：R_m = H − t* − 𝟙[入口轮为 m−1]（15 个 K 经一轮 m−1 进入、10 个 K 直接以 m 进入）。
5. descent 观察：t* = 2（50/51；K=6 为 1）；s_1/K **非常数**（4/7, 5/8, 2/3, 7/10, …, 16/25）——S1 时代"0.64K"是其四个 K 的精确事实但不外推（范围修正，已登记）。

**Scope 声明**：以上全部为 tested-range 计算证据。51 个 K ≠ 一般 K；canonical 初态 + XO kernel 以外未测；H = 10K 截断语义不变（censoring）。

**[2026-09-30 proof deposit]** **S3 CLOSED — PASS**；[证明本体](S3_XO_BOUNDARY_THEOREM_PROOF.md)已入库，[登记](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)为 **REGISTERED — PROOF DEPOSITED / INDEPENDENT AUDIT PASS**。T1 保持 [Existing Math]；T2–T5 为 [New Mathematical Result]；[独立审计](S3_XO_BOUNDARY_THEOREM_AUDIT.md) GPT-6.1 Sol：PASS — no substantive gap found。§3a–3c 保留当时的数学/计算证据分层；一般 canonical XO 的 P1–P4 现由新证明覆盖，S1/S2 的有限 H 分类并未自动升级。

## 4. Current XO Mechanism Picture（可理解版，非形式证明）

**为什么不变量坐标守恒**：当 selection s ≤ m 时，被扫走的只有前缀质量，而 XO 喷洒律在 x = m/K 处 CDF 已达 1——**全部重撒质量都落回前缀 1..m**。所以 F_m（前 m bin 总质量）一轮不少地保恒：F_m ≡ 2b − b²。这就是 [Math] 定理 D_{t,m} = b(1−b) 的机制化读法（同一事实的 CDF 形式）。

**边界变量与开关**：d（最大超额）在尾部冻结后只剩两个竞争者：D_m（= anchor，恒定）与 D_{m−1} = anchor + (1/K − p_m)。于是**整个边界博弈只看一个数**：bin m 的质量 p_m 与均匀份额 1/K 的比较——p_m > 1/K 时选 m（d = anchor，平台），p_m ≤ 1/K 时选 m−1（d 抬高 (1/K − p_m)，出现"反弹"）。

**两个 regime 的来源**：s = m 轮里 bin m 被扫掉 ¼ 又按 q_K 回填；s = m−1 轮里 bin m 不被扫、只按漏率 q_K 进账。回填速率与均匀份额之比 F_m·q_K / (1/K) 决定 tested 轨迹的走向：**≥ 1 时 p_m 持续高于阈值 → 分类为 LOCK（H 内）**；**< 1 时 p_m 衰减下穿、再被漏回 → 分类为 CYCLE（H 内）**。而 F_m·q_K ≥ 1/K 经闭式化简恰好等价于 **K mod 6 ∈ {0,1,2}**——residue 类依赖的来源是喷洒律**支集截断 5/6 与网格 1/K 的对齐**：截断使边界 bin 的接收质量 q_K = (6 − (K mod 6))/(3K) 随余数改变，从而产生 tested 范围内的 residue 类分类差（"共振"一词仅作类比，非已建立的机制名称）。

**状态并非不动点**：锁定期间 (s, d) 冻结，但前缀内部质量继续按 (1−α) 型收缩演化——"平台"指观测量，不是状态。

## 5. Post-S2 当时的开放问题（历史；P1–P4 已由 S3 证明覆盖）

**Proof-level（S3 的定理化清单）**：
- 内部边际引理：post-entry 是否恒有 max_{j ≤ m−2} D_j < max(D_{m−1}, D_m)（C2 在 51 K 上零违例的解释）——判据定理化的关键前置；
- 判据定理：对一般 K ≥ 6，LOCK ⟺ F_m q_K ≥ 1/K ⟺ K mod 6 ∈ {0,1,2}（目前是 51 K 的 [VCR] + [Hyp]）；
- 约化命题的形式化（给定边际引理，C6 恒等的证明是短代数）；
- descent 级联的规律（s_1/K 非常数——descent 是否有自己的精确结构）。

**Dynamics-level**：
- φ_K 的 dwell 结构与 mean-rounds-per-rebound 趋势（tested 范围：mod-3 12.57→10.16、mod-4 ~2.5、mod-5 3.60→3.26；趋势观察，无极限值主张）的解析解释；
- 锁定态状态演化是否收敛到前缀仿射映射的显式不动点（|Δp| ~ (1−α)^t 衰减样 [Obs]，无渐近断言）；
- H = 10K 之外的行为（censoring 纪律：不可外推；"mod-3 K 末端短 m-尾"是轨迹结束于 A-dwell 中段，非周期衰亡证据）。

**Modeling-level**：
- 机制对非 canonical law 的可迁移性（mod-6 的"6"来自 XO 断点分母；一般 law 的断点分母/漏率结构——S11 方向）；
- α 依赖（不变量对 α 无关是现成 corollary 候选——S5 登记项）；
- 初态依赖（suffix-dominance 区域的吸引性——S6）；
- 2D 回归（S7）：1D 机制是否在 radial-collapse/真 2D selection 下存活、变形或消失。

## 6. Recommended Next Research Stage（评估，不启动）

**[2026-09-30 状态]** **S3 CLOSED — PASS，证明已入库，独立审计 PASS**（见[登记](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)）。下一步为 synthesis / documentation，其他分支不自动启动。以下为 deposit 前的定位建议，仅作历史保留。

**S3 方向（推荐主线下一步）——重新定位为 Bounded XO Boundary-Reduction Proof Attempt**：(1) 证明 canonical 轨迹到达两边界 regime（t\* 存在且 post-entry argmax 限于 {m−1, m}）；(2) 证明内部坐标此后不能超越边界对（内部边际引理）；(3) 在此条件下形式化标量递推与 locking-vs-repeated-crossing 判据（含精确 tie 处理：p_m = 1/K 归 m−1 支）。**不列入 S3**：dwell 渐近、full symbolic-dynamics 分类、新的大规模计算 sweep、2D 实现——S3 是 bounded proof attempt，非开放理论期。收益：把 51 K 的 [VCR] 升级为一般 [Math]（若成功），完成"物理问题 → 精确机制"的闭环；风险：边际引理可能需要真实工作（若失败，失败模式本身 informative——定位哪一段归纳断裂）；成本：中等（归纳模板已有，数据完备）。Map 允许 S3a 提前（"视 wave 1 结果提前"）——S2 的结果正是该条款预期的那种触发。

**S4 / S8（Wave 1 剩余 cheap probe）**：S4 诊断面板须以 exact 口径执行（S1 裁定）；S8 文献核验无风险。二者可与 S3 并行或先后，均为 Wave 1 收口（G-S4-C）所需。

**2D 方向（S7-C1）**：意义 = 1D 机制现已完全显式，radial-collapse 同构候选可以对照一个"已知答案"的 1D 基线来陈述——这是 2D 回归的最佳时机窗口；但成本高于 S3，且其同构 candidate 尚未核验。建议排在 S3 之后（与 Map 的 wave-2 顺序一致）。

**Stop / synthesis 方向**：已可成书的完整弧线 = 物理问题 → 1D 模型 → one-step 数学 → 反事实 controls → XO 精确机制（S1+S2）。现在收口的代价：mod-6 未证明（只有 51 K 证据）、2D 未回答、S4/S8 未做——三条都是廉价可补的缺口。综合判断：**不建议现在终稿**；建议 S3 先行（定理化是最高信息增益、直接延续），S4/S8 随后收口 Wave 1，G-S4-C 再定 S7-C1。最终决定权在 Owner + GPT。

## 7. Owner-Level Explanation（非数学专业版）

**发现了什么**：你们沙坑模型里那个"永远扫不完"的怪现象，现在有了一个精确的、经过验证的解释（在测试过的范围内）。第一，电脑没有骗人——平台高度是真的；但它确实"编造"了大部分反弹：用精确分数算术重算后，一半的测试配置在模拟范围内**一次真正的反弹都没有**——扫沙位置走到一个格子后直到模拟终点都没再离开，表面高度纹丝不动（底下沙子还在缓慢流动）。第二，哪种行为会发生，由一个简单得惊人的规则决定：**看格子数 K 除以 6 的余数**——余 0、1、2 的所有测试格子都一路锁定到模拟终点（分类 LOCK），余 3、4、5 则在两个格子之间来回切换直到终点（分类 CYCLE）——这是有限模拟内的行为，"永远"尚未证明。第三，**边界上发生了什么**可以被压缩成**一个数字的更新规则**：某个特定格子的沙量围绕"平均值"上下穿越，就产生切换——注意这只是边界上的博弈，系统其余部分仍在许多维里继续演化。

**为什么重要**：(1) 它证明这些现象是真实算术的性质，不是计算机误差——这解除了此前的最大解读风险；(2) K mod 6 这种规则意味着现象有精确的结构原因（喷洒律 5/6 支集截断与网格对齐决定边界格子的回填量），不是杂乱复杂性；(3) 边界机制被化简成一行公式，意味着后面的数学工作（证明、推广、加噪声）都有了明确的靶子。

**还有什么不确定**：规则在测试过的 51 种格子（最大 400）上全部应验，但**还不是对所有格子大小的数学证明**；只测了一种喷洒律、一种初始沙堆；模型是一维的，而真实沙坑是二维的——2D 里这套机制是否还在，完全未知（这正是 Map 里 2D 回归分支要回答的）。

## 8. Validation appendix（本轮 consolidation 做了什么）

1. **脚本核验**（对照 committed artifacts，非对话记忆）：S1——exact ties 497/199/1997/799、exact N_R 0/398/0/1598 vs float 41/581/56/2380、lock e≡0（K=50/200）、boundary==tie 7500/7500、KC 复现、d 追踪 ≤5.3e-15；S2——census 51 行（25 LOCK + 26 CYCLE、51/51 match、每 K 满足 mod-6 规则）、C1–C7/C6 全 51/51、t* 分布、三类指纹（maxA/maxB）、零 events、冻结 config 与 `e706a0d` 逐字节一致。全部 PASS 除下述一项。
2. **发现并修正一处 prose 过度概括**：S2 报告曾写"LOCK K 尾部 R_m = H − t* − 1"——精确核验显示正确公式为 R_m = H − t* − 𝟙[入口轮为 m−1]（15 K 经 m−1 进入、10 K 直接进入；引用例 K=50/200/6 均属前者故原句在其上成立）。已按"措辞最小更正"纪律就地修正 canonical 报告一句；原始 CSV 与其余结论不受影响。
3. 测试口径：`tests/count_tests.py` = **382 项运行时断言全部通过**（本轮无代码改动）。

## 附：关键文档与提交索引

| 内容 | 位置 | 提交 |
|---|---|---|
| Stage 4 规划（内部历史文件，未收入公开快照） | 私人研究档案 | a41fd49 |
| S1 canonical 报告 | docs/S4_W1_S1_XO_PRECISION.md | f2e23db / 1514a3c |
| S2 冻结 config | experiments/s4_w1_s2_boundary_microscope/config.json | e706a0d |
| S2 canonical 报告 | docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md | c5e09e4 / 102c86a |
| 既有 XO 不变量 T1（原证明已重录） | docs/S3_XO_BOUNDARY_THEOREM_PROOF.md §2 | 私人原记录 a872169 |
| 措辞纪律（内部历史文件，未收入公开快照） | 私人研究档案 | 15f62eb |
