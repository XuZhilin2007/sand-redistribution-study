# S4 Wave 1 — S2：XO 边界切换显微镜执行报告（canonical）

> **状态更新（2026-09-30）**：Canonical XO **S3 CLOSED — PASS**；[P1–P4 证明已入库](S3_XO_BOUNDARY_THEOREM_PROOF.md)，[独立审计 PASS](S3_XO_BOUNDARY_THEOREM_AUDIT.md)，[登记完成](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)。下文保留本文件形成时的实验/规划范围；其中 canonical XO 的 entry、persistence、scalar closure 与 residue dichotomy 证明目标已由 S3 覆盖，不再是当前待执行任务。S1/S2 仍为 [VCR]，LOCK/CYCLE 仍为有限 H 标签；未覆盖的推广问题不自动启动。
> **Tie 说明**：S3 一般陈述中，x_t=1/K 归 m−1 支且 discrepancy 为 exact tie；B 支严格收缩要求 x_t<1/K。原 S2 census 无此 boundary-coordinate equality，冻结检查与原产物保留。

> 日期：2026-09-29。性质：Stage 4 Wave 1 S2 执行轮完成报告（冻结 config `e706a0d` 的首次执行；config 零改动）。
>
> 冻结 config：`experiments/s4_w1_s2_boundary_microscope/config.json`（sha256 `a94508f4…`，冻结于任何 S2 计算之前）。本报告遵守其 evidence_discipline 与措辞规则。
>
> 证据分级：**[Theorem]** = 既有 Model-derived Mathematical Result（本轮只核验、不扩张）；**[VCR]** = Verified Computational Result（本轮 exact 运行的计算结果，限 tested K）；**[机制状态]** = 候选机制的存活/证伪判定；**[假设更新]** = mod-6 假设的判定。全文不使用 "proved"（除指称既有 [Theorem] 的正式地位）。

---

## 1. 执行摘要（What ran）

- **对象**：冻结 config 的完整执行——Step 0（committed-CSV 门）→ Step 1（exact K census）→ Step 2（约化核验 C1–C7）。
- **算术**：exact rational（`fractions.Fraction`，S1 引擎及其六层核验继承）；无任何 float/mp 参与证据。
- **K 覆盖**：census K = 6..53（48 值）+ anchors {100, 200, 400}（扩展仪表化重跑），共 51 个 K，H = 10K。
- **运行**：1275.8 s（K=400 主导）；无 hard halt；**scientific events = 0**；`results/metadata.json` 记录 git SHA 与全部输出 sha256。
- **冻结件合规**：K 集、分类器、horizon、算术、核验判据零改动；分类器/检查语义与冻结文本逐一对应。

## 2. Step-0 结果（committed-CSV 门）

全部通过（详见 `results/step0_s1_csv_checks.csv`）：

| 检查 | 结果 |
|---|---|
| S0a 锚点分类复现 | LOCK/CYCLE/LOCK/CYCLE @ K=50/100/200/400 ✓ |
| S0b F_m = 4·moved（s=m 轮） | 4 锚点全部 ✓（坏轮 0） |
| S0c lock 签名 e ≡ 0（K=50/200） | ✓ |
| S0d 分解恒等式（A = ties+rebounds；contractions = B + t*） | 4 锚点全部精确成立（497=497+0、597=199+398、1997=1997+0、2397=799+1598；403=401+2 等） |
| S0e 分支 B 方程（CSV 恢复，容差 1e-50） | B→B 对 4/4、B→A 上穿 1998/1998；最大残差 3.3×10⁻⁶³（= 60 位存储截断尺度，≪ 容差） |

## 3. Census 表（Step 1）

**主结果：mod-6 预测 51/51 全中，零 OTHER，零 mismatch**（逐 K 明细见 `results/per_k_census.csv`）：

| K mod 6 | 预测 | 观测 | K 值 |
|---|---|---|---|
| 0, 1, 2 | LOCK | **LOCK（25/25）** | 6,7,8,12,…,50（含 200） |
| 3, 4, 5 | CYCLE | **CYCLE（26/26）** | 9,10,11,15,…,53（含 100,400） |

结构事实（全部为 [VCR]，限 tested K）：

- **t\***：首次进入 {m−1, m} 在 50/51 个 K 上为 t\*=2（K=6 为 t\*=1）；descent 行程逐 K 记录在案。
- **LOCK K（冻结分类器标签，H=10K 内）**：入口后 selection 不再离开 m，尾部游程延伸至 horizon 末端——R_m = H − t\* − 𝟙[入口轮为 m−1]（15 个 K 经一轮 m−1 进入：K=50: 497、K=200: 1997、K=6: 58；10 个 K 直接以 m 进入：R_m = H − t\*。精确公式由 post-S2 consolidation 核验所定；原句"H − t\* − 1"在所引三例上成立但对全部 25 个 LOCK K 过度概括），H 内零 rebound、d = b(1−b) 精确常值（tie-tail = R_m）。H 之外的行为未测（censoring 纪律）。
- **CYCLE K（冻结分类器标签，H=10K 内）的三类切换指纹**（全部由已验证的标量边界递推精确生成，见 §4 C6；selection 限制于 {m−1, m}，**非严格交替**——存在同支连续驻留，如 K=100 的序列含 BAA 型段）：
  - **mod-3（K=9,15,…,51）**：慢切换——mean rounds per rebound ~10.2–12.6（A-dwell 长达 9–11 轮的连续精确 tie，B-dwell ≤2）；该均值在 tested K 内由 12.57（K=9）降至 10.16（K=51）——趋势观察，不作极限值主张；
  - **mod-4（K=10,16,…,400）**：紧凑切换——mean rounds per rebound ~2.5（max branch dwell：A=2、B=3），即 S1 在 K=100/400 看到的指纹；
  - **mod-5（K=11,17,…,53）**：B-重切换——mean rounds per rebound ~3.3（maxA=1, maxB=7–8，长 B-dwell）。
- **descent 修正记录**：s_1/K 在 census 中**非常数**（4/7, 5/8, 2/3, 7/10, …, 16/25）——S1 在其四个 tested K 上观察到的 s_1 = 0.64K 是该四点的精确事实但不外推；descent 级联的真正规律成为新问题（§6e）。
- mod-3 K 的 horizon 末端出现 R_m = 3–7 的短 m-尾（仍判 CYCLE，R_m < K）——由映射语言：轨迹恰结束于一次 A-dwell 中段；非"周期衰亡"证据，不作渐近解读。

## 4. 约化核验报告（Step 2，C1–C7）

**全部通过，51/51，零 scientific events**（逐轮明细 `results/reduction_K{K}.csv`）：

| 检查 | 内容 | 结果 |
|---|---|---|
| C1 | F_{t,m} = 2b − b² 于全部 H+1 个状态（Fraction 相等） | 51/51 ✓（违例 0） |
| C2 | post-entry argmax ⊆ {m−1, m}（内部坐标严格低于边界对） | 51/51 ✓（违例 0；p_m = 1/K 精确 tie 轮另行计数） |
| C3 | 分支方程（A: p′=¾p+¼F_m q_K；B: p′=p+¼(F_m−p)q_K，实测前缀和）逐轮 Fraction 相等 | 51/51 ✓（0 mismatch） |
| C4 | s=m ⟺ p_m > 1/K（含 tie 归属） | 51/51 ✓ |
| C5 | d_t = b(1−b) + (1/K − p_m)⁺ | 51/51 ✓ |
| C6 | 标量边界递推（常值 F_m 的 p_m 映射）与 full 轨道的边界坐标（p_m 及逐轮 tie/rebound/收缩类别）逐步 Fraction 恒等 | **applicable 51/51，exact 恒等 51/51** |
| C7 | q_K = (6 − (K mod 6))/(3K) vs 定义值；descent 事实 | 51/51 ✓ |

无首次失败位置（无失败）；无 hard halt（H1–H4 未触发）。

## 5. 解读更新（Interpretation update）

**Existing theorem**：[Theorem]（s_t ≤ m、D_{t,m} = b(1−b)、tolerance stop 不触发）地位不变；C1/C2 构成其在 51 个 K 上的逐轮计算核验实例（非新证明）。

**Verified computational result（[VCR]，限 tested K = 6..53, 100, 200, 400）**：
1. **边界坐标标量闭包在全部 tested K 上严格成立**：post-entry 边际动力学（selection 与 discrepancy 的演化）与常值 F_m 的标量 p_m 递推逐步 Fraction 恒等（C6）——XO 的全部**边界观测量**（LOCK/CYCLE 分类、tie/rebound/收缩分类、dwell 结构、窄带幅度）由该显式有理递推精确生成；full 状态向量继续多维演化（非不动点），本结果不主张整个系统一维化；
2. **mod-6 分类 51/51**（冻结分类器，H=10K 内）：LOCK ⟺ K mod 6 ∈ {0,1,2}（25 K）、CYCLE ⟺ K mod 6 ∈ {3,4,5}（26 K），含 K=6..12 的设计期预测逐一应验（6,7,8 → LOCK；9,10,11 → CYCLE；12 → LOCK）；
3. 三类 residue 切换指纹（§3，H 内）；
4. 机制代数的全部辅助恒等式（分解恒等式、q_K 闭式、F_m 守恒、selection 阈值）精确成立。

**Candidate mechanism status**：**存活**——候选机制在冻结判据下经受住全部 51 个 K 的精确核验，未出现任何 scientific event。"XO 边界切换机制是否真实到值得更深层理论"的裁定：**是**（在 tested 范围内；约化仍非定理）。

**Hypothesis update**：mod-6 规则作为 **tested-range 规则性存活**（51/51）。按冻结纪律它**不是**"被发现的定律"——tested 范围是 K ≤ 400 的 51 个点，一般性属于 S3 的证明目标。

**新 unlocked 问题（喂给 S3，非本轮任务）**：
- (a) **内部边际引理**：C2 在 51/51 个 K 上零违例——为什么内部坐标在 post-entry 永远追不上边界对？这是把 mod-6 规则定理化的关键引理；
- (b) **判据定理化**：LOCK ⟺ F_m·q_K ≥ 1/K ⟺ K mod 6 ∈ {0,1,2}（K ≥ 6）能否在 (a) 之上完成归纳证明；
- (c) **residue 类 dwell 结构**：三类 mean-rounds-per-rebound 趋势（tested 范围：mod-3 12.57→10.16、mod-4 ~2.5、mod-5 3.60→3.26）与分支驻留分布的解析解释（无极限值主张）；
- (d) **descent 级联规律**：s_1/K 非常数——descent 的真正结构（S1 时代 0.64K 表述的范围修正已登记）；
- (e) **S12 良定化**：LOCK/CYCLE 分类对应标量递推在阈值 1/K 处的行为（H 内）——噪声注入 probe 现在有精确的实验对象。

## 6. Repository state

- **新增**：`experiments/s4_w1_s2_boundary_microscope/`：`run_s2.py`（执行器：Step-0 门 + census + C1–C7 + H2 锚点逐字符串核验 + 冻结分类器）、`results/`（step0_s1_csv_checks.csv、per_k_census.csv、reduction_K{6..53,100,200,400}.csv ×51、metadata.json）；`tests/test_s4_s2_microscope.py`（T1–T4）。
- **未动**：冻结 config（sha256 校验一致）、S1 全部产物、历史文档。
- **测试**：`tests/count_tests.py` = **382 项运行时断言全部通过**（349 + 本轮 33）。
- **提交**：feat + docs 两段式（见 commit SHA）；**push = no**。

## 附：边界与措辞合规

censoring 语义不变（H=10K 截断；mod-3 K 的短 m-尾不作"周期衰亡"解读）；无新 kernel、无 K 扩展、无 S3 证明工作、无 2D、无物理实验（scope reminder 全部遵守）；"proved" 仅指 [Theorem] 的既有地位；mod-6 全文为 hypothesis/规则性表述。
