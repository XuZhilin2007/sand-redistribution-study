# S4 Wave 1 — S1：XO 高精度重放验证报告（canonical）

> **状态更新（2026-09-30）**：Canonical XO **S3 CLOSED — PASS**；[P1–P4 证明已入库](S3_XO_BOUNDARY_THEOREM_PROOF.md)，[独立审计 PASS](S3_XO_BOUNDARY_THEOREM_AUDIT.md)，[登记完成](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)。下文保留本文件形成时的实验/规划范围；其中 canonical XO 的 entry、persistence、scalar closure 与 residue dichotomy 证明目标已由 S3 覆盖，不再是当前待执行任务。S1/S2 仍为 [VCR]，LOCK/CYCLE 仍为有限 H 标签；未覆盖的推广问题不自动启动。

> 日期：2026-09-29。性质：Stage 4 Wave 1 首轮（cheap probe S1）的验证任务完成报告。
>
> 任务来源：Stage 4 Wave 1 冻结 prompt（XO High-Precision Replay）；当时的私人 Stage 4 Research Map S1 为规划出处，未收入公开快照。
> 预注册：本轮全部 config（層级、K 集、可行性门、比较量与判类规则）在**任何计算运行之前**冻结于 `experiments/s4_w1_s1_xo_precision/config.json`（`pre_registration` 字段），运行后未做任何调整。
>
> 证据分级（本文所有结论逐条标注）：
> **[Theorem]** = Model-derived Mathematical Result（既有登记，本轮不升级、只做计算核验）；
> **[VNA-X]** = Verified Numerical Analysis（exact rational tier，exact Fraction 逐轮计算）；
> **[VNA-HP]** = Verified Numerical Analysis（mpmath 50/100 dps ladder）；
> **[RCO]** = Robust Computational Observation（跨全部算术层一致的数值结构）；
> **[MDO]** = Model-dependent Observation（依赖模型的解读性观察）；
> **[WH]** = Working Hypothesis（未验证解释）。
> 本轮**不产生任何新 theorem**；exact tier 的逐轮核验是计算实例，不是证明。

---

## 1. 测试了什么（What was tested）

- **对象**：XO kernel（Witness-Opposed，Gate 2 冻结）在 canonical 初态、α=1/4、exact tolerance 1e-12、H(K)=10K horizon 下的完整动力学；K ∈ {50, 100, 200, 400}（与 Gate 2 完全相同的冻结矩阵）。
- **问题**（按冻结 prompt）：高精度算术下 XO 的异常数值行为是否保持——discrepancy 轨迹 d_t、不变量坐标 D_{t,m}、selection 序列、stopping 行为、rebound 相关诊断；哪些观察 robust / precision-sensitive / unknown；switching 结构是否稳定。
- **Sanity 对照**：KC（canonical kernel）同参数重放，应复现 committed Gate 1 的 t_stop = 324/699/1445/2876 与逐轮分类。
- **明确不做**（冻结边界，已遵守）：无新 redistribution law、无参数 sweep、无 2D、无物理实验、无一般渐近理论、无 symbolic dynamics 形式化。

## 2. 精度方法（Precision methodology）

四层算术阶梯（tier ladder），全部实现同一冻结语义（`run_m1a_deterministic` 的数学规则：min-argmax selection、α 前缀 removal、M=αF_s 全域 respray、d ≤ 1e-12 停机）：

| Tier | 算术 | 说明 |
|---|---|---|
| F (float64) | numpy float64 | 现行 committed 实现，**经 committed runner 代码路径再生**并与 committed CSV 逐字符串核验后作为基线（不另写重复 CSV） |
| P50 | mpmath dps 50 + 15 guard digits | 同一数学语义；bin probs 不做 renormalization |
| P100 | mpmath dps 100 + 15 guard digits | 同上 |
| X (exact) | `fractions.Fraction` | **exact real-arithmetic 轨迹**：所有 bin probs 在执行网格上严格有理（near CDF = 2x−x²、XO CDF 分段线性有理系数，在有理边 j/K 上取严格有理值；telescoping sum 逐例断言 == 1）；tolerance = 1/10¹²；tie 判定 worst == 0 严格可决 |

**与 Research Map 的偏差（已登记 RESEARCH_LOG）**：Map S1 曾判定 exact Fraction 路线不可行（"near-CDF 有根式"）。仓库 reality：根式只存在于 ppf（Monte-Carlo 采样路径），deterministic expected-mass 模型从不使用；CDF 在有理网格点上严格有理，故 X tier 可行并被加入。Map 的 mpmath ladder 原样保留。

**六层核验协议**（全部通过后才采信结果）：

1. **Stage V 字符串核验** [VNA-HP]：committed per-round CSV 的 sha256 与 metadata 记录一致（Gate 1 + Gate 2）；经 committed runner 模块再生 XO 全部 7500 active rounds + KC 全部 5344 rounds，**逐行逐字段字符串零失配**。此后 F tier = committed record。
2. **Exact tier 内嵌逐轮自检** [VNA-X]：每轮以 Fraction 相等验证 (a) 单步分解恒等式推论 **max_j E_t(j) = d_{t+1} − d_t**（M1C.2 恒等式的直接推论）；(b) **D_{t,m} = b(1−b) 逐轮严格成立**；(c) s_t ≤ m；(d) 未停机轮 d_t > 1/10¹²。四个 K 全部通过（K=50: 500 轮、K=100: 1000、K=200: 2000、K=400: 4000）。
3. **Repo helper 交叉验证**：exact 诊断 worst/γ 与 `kernel_theory.rebound_residual` / `target_matched_excess`（float）在配对状态上一致至 ≤ 1e-15。
4. **mp ↔ exact 互检** [VNA-HP]：K=50 前 40 轮 selection 序列逐轮一致；d_t 一致至 < 1e-38（P50）/ < 1e-55（P100）。
5. **独立重实现核验**：一个**完全独立编写**（未读 replay 引擎）的 Fraction 重实现（仅依据冻结 spec）对 K=50 全程 500 轮复算——selection 序列（25, 32, 41, 42, 42, …）、0 rebound / 3 contraction / 497 tie、d_3 起 = 84/625 严格常值、D_{t,42} = 84/625 于全部 501 个状态成立——**全部吻合，零差异**。
6. **小 K 手算校验**：K=6 exact 重放 50 轮，初始 D profile = (5,8,9,8,5,0)/36、s_0=3、d_0=1/4、不变量 D_{t,5} = 5/36 逐轮成立（[Theorem] 的 K=6 非平凡实例核验）。

运行成本：exact tier 每 K 0.67 / 7.7 / 21.8 / 941.2 s（可行性门 300/900/2700 s 全部通过，四个 K 均取得 exact tier）；全程 1152 s。

## 3. 数值结果（Numerical results）

### 3.1 主结果总表（XO，全部 tiers 对照）

| K | m | exact 平台后结构 | exact ties | exact N_R | float N_R（其中 tie 噪声伪造） | mp100 N_R | exact 平台带 \|e\| (t≥10) |
|---|---|---|---|---|---|---|---|
| 50 | 42 | **selection 锁定 s=42 至 horizon 结束（497 连续轮）** | 497/500 | **0** | 41（41 全部伪造） | 191 | **0（d = 0.1344 严格常值）** |
| 100 | 84 | **selection 两边界切换 {83, 84}**（H=10K 内；非严格交替——同支驻留 dwell 1–3） | 199/1000 | **398** | 581（183 伪造） | 479 | ≤ 7.88×10⁻⁴ |
| 200 | 167 | **锁定 s=167（1997 连续轮）** | 1997/2000 | **0** | 56（56 全部伪造） | 376 | **0（d = 0.137775 严格常值）** |
| 400 | 334 | **两边界切换 {333, 334}**（H=10K 内；非严格交替） | 799/4000 | **1598** | 2380（782 伪造） | 1719 | ≤ 1.94×10⁻⁴ |

exact N_R = 严格 d_{t+1} > d_t 轮数（exact rational 比较）；伪造 = 该 tier 的 rebound 旗标落在 exact-tie 轮上（见 3.4）。

### 3.2 Selection 行程（itinerary）——跨全部算术层逐轮一致 [RCO]

四个 K、四种算术（F/P50/P100/X）的 selection 序列 s_t **逐轮完全相同，零分歧**（共 7500 轮 × 3 组 tier 对）：

- K=50：25 → 32 → 41 → 42（此后锁定，497 轮不动）；
- K=100：50 → 64 → {83↔84} 交替至 horizon；
- K=200：100 → 128 → 166 → 167（锁定）；
- K=400：200 → 256 → {333↔334} 交替至 horizon。

共性结构：前 3 轮单向下降 s_0=K/2、s_1=0.64K、s_2≈0.83K，随后进入 m=⌈5K/6⌉（[Theorem] 的 s_t ≤ m 逐轮核验 [VNA-X]）。**distinct s 恰为 4 个**（每个 K、每个 tier）。LOCK/CYCLE 二分（冻结分类器语义，H=10K 内）：K=50、200 分类为 LOCK；K=100、400 分类为 CYCLE（两边界切换 {m−1, m}，非严格交替）。注意 K=50 与 K=100 共享 b=21/25、K=200 与 K=400 共享 b=167/200——**b 单独不决定二分**（K 本身进入判据）[RCO；机制 = WH]。

### 3.3 d_t 轨迹与平台

- **float64 全程追踪 exact 真值至机器精度**：max_t |d_t^F − d_t^X| = 7.5e-16 / 1.9e-15 / 2.0e-15 / 5.2e-15（K=50/100/200/400，跨 500–4000 轮）[VNA-X + VNA-HP]。
- **锁定 K（50/200）**：exact 轨迹 3 轮收缩后 d_t 严格常值且**恰好等于不变量坐标**：K=50: 1/4 → 401/2500 → 8579/62500 → **84/625 = 0.1344 = b(1−b) 精确**（e_t ≡ 0 自 t=3）；K=200: → **5511/40000 = 0.137775 精确**。即平台不是"贴近 anchor"而是**d 与不变量坐标严格重合**（e_t = 0，exact Fraction 相等）[VNA-X]。
- **循环 K（100/400）**：d 在窄带内真实振荡，exact |e_t| ≤ 7.88×10⁻⁴（K=100）/ 1.94×10⁻⁴（K=400）[VNA-X]。
- **锁定期间状态并非不动点**：K=50 锁定后 max|Δp| 逐轮仍 ~3.0e-3 → 2.3e-3 → 1.7e-3（比值 ≈ 1−α），状态继续演化而 (s_t, d_t) 冻结 [VNA-X；渐近含义不做断言]。
- **Stopping**：全部 tier、全部 K 均在 H=10K censoring（无停机）；exact tier 逐轮 d_t > 1/10¹² 核验通过——与 [Theorem]（tolerance stop 在 exact real-arithmetic 模型中不可触发）一致，本轮为该结果的**计算核验实例**（非新证明）。

### 3.4 Rebound 诊断的分解（本轮核心新事实）

按 exact 轮类（TIE: d_{t+1}=d_t 精确；REB: 严格升；CON: 严格降）分解各 tier 的 rebound 旗标 [VNA-X]：

| K | exact (TIE/REB/CON) | F 旗标落点 | mp100 落点 | mp50 落点 |
|---|---|---|---|---|
| 50 | 497 / 0 / 3 | REB 0/0，TIE 41，CON 0 | TIE 191 | TIE 186 |
| 100 | 199 / 398 / 403 | **REB 398/398**，TIE 183，CON 0 | REB 398/398，TIE 81 | REB 398/398，TIE 106 |
| 200 | 1997 / 0 / 3 | REB 0/0，TIE 56，CON 0 | TIE 376 | TIE 215 |
| 400 | 799 / 1598 / 1603 | **REB 1598/1598**，TIE 782，CON 0 | REB 1598/1598，TIE 121 | REB 1598/1598，TIE 65 |

三条精确规律（每个 K、每个有限精度 tier）：

1. **真实 rebound 零漏报**：exact REB 轮全部被该 tier 旗标（这些轮的 |worst| ≥ 5.9×10⁻⁵，远高于噪声）；
2. **严格收缩轮零误报**（CON 轮无一被旗标）；
3. **全部伪造集中于 exact-tie 轮**，且伪造率随算术实现剧烈变化：float 在 tie 轮上的伪造率 K=50/100/200/400 = 8.2% / 91.9% / 2.8% / 97.9%，mp100 = 38.4% / 40.7% / 18.8% / 15.1%，mp50 = 37.4% / 53.3% / 10.8% / 8.1%——**tie 轮上的 strict 有限精度比较是四舍五入噪声的方向游戏，任何有限精度都无权威性**（F 旗标集 ⊋ exact 旗标集在全部 K 成立：float 只多报、不漏报、不误报收缩轮）。

### 3.5 exact tie 的可判定性与 1e-12 诊断 [VNA-X]

exact |worst| 严格双峰：TIE 轮 = **0（精确）**；非 TIE 轮最小 |worst| = 2.9×10⁻³ / 2.4×10⁻⁴ / 7.4×10⁻⁴ / 5.9×10⁻⁵（K=50/100/200/400）。1e-12 诊断容差落在 ≥ 7 个数量级的空隙中——因此 **float numerical-boundary 旗标与 exact tie 逐轮完全一一对应（4×100%）**。committed Gate 2 的 boundary 统计（497/500、199/1000、1997/2000、799/4000，46.6% 总体）是对 exact tie 结构的忠实检测。

### 3.6 KC sanity [VNA-X + VNA-HP]

- t_stop：exact K=50 = **324**（= float）；mp100 = 324 / 699 / 1445 / 2876（= float，全部 K）。
- N_R @ K=50：exact = 117 = float（KC 无 numerical-boundary 轮，其计数本就非边界敏感；现在 exact 确认）。
- canonical 参考行为在精度阶梯下完全稳定。

## 4. 与既往观察的对照（Comparison with previous observations）

| 既往观察（出处） | 本轮裁定 |
|---|---|
| XO 窄平台 D_max ≈ 0.134–0.138（Gate 2 当时的浮点观察） | **确认并锐化**：锁定 K 的平台 = anchor **精确值**（e=0）；循环 K 窄带 exact ≤ 8×10⁻⁴。"平台 ≈ anchor" 的贴近关系在锁定 K 上是**严格相等**（[VNA-X]；[Theorem] 本体不变） |
| selection 集中于少量离散边界值（~4 个） | **确认**：每 K 恰 4 个 distinct s，且行程跨全部算术层逐轮一致（升级为 [RCO]） |
| XO 数值边界轮 46.6%、强 K 双峰 99.4/19.9/99.9/20.0%（Gate 2） | **确认其真实性并给出 exact 含义**：= exact tie 份额（逐轮 1:1 对应）；双峰 = 锁定/循环二分的 tie 份额差 |
| XO rebound 计数 N_R = 41/581/56/2380、f_R = 0.082/0.581/0.028/0.595（Gate 2 float） | **重大修正**：exact N_R = **0 / 398 / 0 / 1598**。K=50/200 的 float rebound 现象整体是**浮点伪造**（exact 无任何 rebound）；K=100/400 的 recurrent 现象真实存在但 float 计数被 tie 噪声抬高 +46% / +49%；任何有限精度计数在 tie 轮上都由噪声决定 |
| Gate 2 / Stage 3 对 XO 的逐 K 定性标签：K=50 "Long-span recurrent / churn-like transient"、K=200 "Short recurrent transient" | **就 exact 实在而言被取代**：在 H=10K 内 K=50/200 的 exact 轨迹是**零 rebound 的锁定平台**（selection 锁定于 m、d = anchor 精确常值、状态仍缓慢演化）。历史文档不改写（append-only 纪律）；替代解读以本轮为准并已登记 RESEARCH_LOG。注意 §7 措辞规则同时禁止反向越界：这不证明 state 收敛、不证明周期性、不加任何 horizon 之外断言 |
| "XO 骑在 real-arithmetic 精确边际上，numerically reproducible but mathematically delicate"（Gate 2 措辞） | **确认并升级为可判定**：骑的就是 exact tie（worst = 0 严格成立，[VNA-X]）；"delicate" 部分现在有了精确刻画（tie 轮的有限精度旗标 = 噪声） |
| "d_t 贴近 anchor 的窄平台（宽度 ~1e-4）"（早期浮点观察） | **修正**：宽度不是统一的 ~1e-4——锁定 K 上精确为 0，循环 K 上 ≤ 8×10⁻⁴；observed 层的表述应由本轮 exact 值细化 |
| KC t_stop / N_R（Gate 1） | **确认**：精度阶梯下完全复现 |
| 不变量 [Theorem]（[公开证明 §2](S3_XO_BOUNDARY_THEOREM_PROOF.md) 重录原证明） | **未动**；获得四 K 全程逐轮 Fraction 级计算核验 + K=6 手算实例（不构成新证明） |

## 5. 解读更新（Interpretation update）

**Confirmed（跨全部算术层成立）**：

1. Selection itinerary 与其全部结构（4 值集中、前 3 轮下降模板、s ≤ m、LOCK/CYCLE 二分（H 内）的**形态**）[RCO]；
2. d_t 数值轨迹（float64 = exact 至 ~1e-15；平台值、窄带宽度、4 值/76 值结构）[RCO]；
3. Stopping/censoring 行为（全 tier 全 K 无停机；[Theorem] 计算核验）[RCO]；
4. boundary 旗标 = exact tie 检测器（逐轮 1:1）[RCO]；
5. KC 全部参考量 [RCO]。

**Revised（相对既往数值观察）**：

1. **XO rebound 现象学（最重要）**：exact N_R = 0/398/0/1598。K=50/200：无 rebound 的锁定平台，float 计数与定性标签是浮点伪造；K=100/400：真实 recurrent 但精确计数只能来自 exact tier。**下游一切使用 XO float N_R/f_R/S_R/定性标签的地方必须切换到 exact 口径或显式标注**。历史 raw CSV 不动。
2. 平台带宽度按 K 二分（0 vs ≤8×10⁻⁴），取代统一 "~1e-4" 表述。
3. 有限精度下的 rebound 旗标集本身（含 mp 层互相不一致）——precision-sensitive by construction。

**Unresolved（诚实保留）**：

1. 锁定/循环二分的**机制**（为何 K=50/200 锁定、K=100/400 循环；b 共享证明 b 不单独决定）——S3 的头号目标 [WH：与 Q6 "rational breakpoints × K 的 arithmetic/grid resonance" 直接相关，但机制未验证]；
2. 锁定态状态漂移（max|Δp| ~ (1−α)^t 衰减样）是否趋于 0 / 状态是否收敛——asymptotic 开放（censoring 纪律：不做 horizon 外断言）；
3. 两边界切换（CYCLE 分类）在 H=10K 之外的持续性与最终命运；
4. 二分结构的一般 K 刻画（仅 4 个 tested K）；
5. e_t 在循环 K 上是否永不逃逸窄带（观测限于 horizon）。

## 6. 建议下一步（Recommended next action）

**主裁定：继续 XO 理论（Map 的三分支结果中的"解禁"分支，附带一项重大解读修正）。**

- d_t / selection / stopping 全部 robust → XO 结构是**真实算术性质**，非浮点效应 → **XO 理论投资（S3a/S3b）解禁**，与 Map 的依赖设计一致；
- 但 rebound 口径必须修正：**exact tier 结果成为 XO rebound/分类的 canonical 参照**（float 记录降级为"边界检测器 + 伪造率已刻画的历史数据"）；
- **S3 的头号具体目标由本轮给出**：解释锁定/循环二分（K=50/200 vs 100/400）——一个精确、可证伪、数据完备的问题（exact 轨迹与全部轮类已入库）；
- **S4 面板执行前必须吸收本轮**：XO 行使用 exact 口径 N_R/轮类，或显式 caveat；
- **S12（tie-break 噪声 probe）问题变得良定**：向锁定平台注入 selection 噪声是否将其踢入循环——本轮提供了 exact 基线；
- Wave 1 其余 probe（S2、S8）按 Map 预算继续，不受本轮影响。

## 7. 措辞与边界合规

- 未使用 "finite-state"（模型始终为 finite-dimensional piecewise-affine switched 系统）；未写 "proved plateau level"（平台=anchor 相等仅锁定 K、仅 exact tier 计算核验，标注 [VNA-X]）；XO 不写 non-convergent/churn（K=50/200 就 exact 实在改称 locked plateau，限 H 内）；censoring 语义未动（H=10K 预注册，无外推）；"exact stop" 表述规则未触及（本轮无停机）。
- 本轮无任何对历史冻结文件/committed CSV 的改写；全部修正以 append 方式登记（RESEARCH_LOG 2026-09-29 S1 条目 + 本文档 §4）。
- 冻结 config 未调；无 sweep；无新 kernel。

## 附：产物与复现

- 目录：`experiments/s4_w1_s1_xo_precision/`（config.json 预注册、replay.py 引擎、run_s1.py 编排、results/ 全部 CSV + metadata.json）；
- 复现：`.venv/Scripts/python.exe experiments/s4_w1_s1_xo_precision/run_s1.py`（Stage V → 阶梯 → 比较，~19 min，输出位级可复现）；
- 测试：`tests/test_s4_s1_precision.py` 19 项断言（`tests/count_tests.py` 总计 **349** 项全部通过）；
- 运行时 metadata 记录的 git SHA 为运行时 HEAD（Stage 3 既有两段式惯例：feat 提交承载代码+结果，docs 提交承载 canonical 文档与日志）。
