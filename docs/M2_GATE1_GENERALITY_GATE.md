# M2 Gate 1 文档：Redistribution-Mismatch Generality Gate——sign / amplitude / shape 三对照

> 状态：**Stage 3 / M2 Gate 1 execution 轮**（2026-09-20）。Stage 2 closure（[RESEARCH_CHECKPOINT_2026-09-20.md](RESEARCH_CHECKPOINT_2026-09-20.md)）后的第一个 Stage 3 实验；Owner + GPT 冻结 prompt 执行。
>
> 性质：新增 2 个解析 kernel（KW/KS）+ 复用 generic fixed-G 支持；canonical simulator 零修改；全部结果 deterministic 一次运行可复现（`python experiments/m2_gate1_generality/run_experiment.py`）。
>
> 实验报告（数值细节）：[experiments/m2_gate1_generality/M2_GATE1_REPORT.md](../experiments/m2_gate1_generality/M2_GATE1_REPORT.md)。

## 1. Research question（冻结）

Stage 2 已确立一般理论（一步分解、G≤T 序定理、精确单步 rebound 判据）与 canonical q_near 的 long pre-stopping churn regime。Stage 3 / M2 的核心问题：

> **Are the rebound / switching / churn mechanisms specific to canonical q_near, or do they arise systematically across a broader class of fixed redistribution laws with positive cumulative mismatch relative to the target?**

Gate 1 只切三个变量：

```text
A. mismatch sign        （K-：Δ ≤ 0 镜像）
B. mismatch amplitude   （KW：同形 1/4 幅度）
C. mismatch spatial shape（KS：同 peak / 同 integral / 同对称性的 flat-top）
```

对 dynamics 的影响：rebound、repeated rebound、switching、long-span churn-like transient、exact stopping / finite-horizon censoring。**不研究**：continuum、Monte Carlo、finite particles、2D、perception、alpha sweep、deadband、new stopping rules。

Stage 2 的 canonical 数字（rebound fraction ≈0.36、boundary gate、a_crit、C-C-R、two-branch map、gate return、stopping geometry）**一律不得当成 general truth**；本轮不重新证明任何定理，只在新 kernels 上逐轮 numerical verification。

## 2. Frozen kernel set（运行前冻结）

| kernel | G(x) | Δ(x) | 关键结构性质 | 角色 |
|---|---|---|---|---|
| K0 | x | 0 | G=T | zero-mismatch integrity reference（复用 M1C.1 Variant U） |
| K- | x² | −x(1−x) | Δ≤0，canonical 的精确符号镜像 | sign control + 定理后盾 negative control |
| KW | 5/4·x − 1/4·x² | x(1−x)/4 | 同 canonical 形状，λ=1/4 幅度（冻结，不做 family sweep） | amplitude ablation |
| KC | 2x−x² | x(1−x) | canonical | calibration anchor（最低限度 integrity reproduction） |
| KS | 分段（见下） | 分段 flat-top | Δ≥0、max=1/4、∫=1/6、对称、端点为 0——与 canonical 的 peak/integral/symmetry/endpoints 全对齐 | spatial-shape control |

```text
Delta_S(x) = 3/4·x        , 0 <= x <= 1/3
             1/4          , 1/3 <= x <= 2/3
             3/4·(1-x)    , 2/3 <= x <= 1
```

**Scope boundary（冻结）**：KS 同时改变 mismatch 正则性（分段线性 vs 光滑抛物线）。KC vs KS 的任何差异只能表述为对 mismatch spatial profile / regularity 的敏感性。

## 3. Exact configuration（全部继承 canonical，未动）

- target `T(x)=x`；`alpha=0.25`；cumulative-excess selection（argmax，tie 取小下标）；prefix removal；deterministic expected-mass path；exact stopping `D_max ≤ 1e-12`；**canonical initial state = q_near bin probabilities，同一 K 下五个 kernel 完全共享**（`run_m1a_deterministic(dist=q_near, redistribution_dist=kernel law)` 的现有 M1C.1 分角色机制，simulator 未改一行）。
- **Discretization**：`q_j = G(j/K) − G((j−1)/K)`，解析 CDF 在 bin 边界上精确求值（`ThrowDistribution.bin_probs`）；不做 midpoint density sampling、不做近似积分、不归一化（仅 canonical fp 残差校正）、**不把 KS 的 1/3、2/3 断点取整到网格**（K=100 时 1/3 落在 bin 33/34 之间，直接求值）。
- **Horizon**：`H(K)=10K` active rounds 或 exact stop 取先；到界未停记为 `not stopped within pre-registered Gate-1 horizon`（finite-horizon censored observation），禁止自动延长或解读为 non-convergence。
- `tau_num = 1e-12` 仅用于 floating validation（mass normalization、分解残差、CDF 检查、numerical-boundary 分类），**不是新的 stopping tolerance**。
- per-round 诊断（每个 active 轮）：t、D_max_before、j_t/a_t、moved mass、D_max_after、rebound boolean、boundary jump（|Δa|，switch 阈值 0.2 = M1C.4 惯例）、exact-stop 状态、分解残差、判据结果、criterion margin / worst point、R_t（M1C.4 语义 `max_{j:γ>0} M·max(Δ,0)/γ`，helper 未重定义）。

## 4. Execution matrix（冻结顺序，全部执行，无 selection bias）

```text
A2 static validation   K-/KW/KS 在各自全部执行 K 上（105/105 required 通过；173/173 含参照行）
A3 K0 integrity        K = 50,100,200,400   vs committed M1C.1 U rows（16 字段逐字符串比较）
A4 KC integrity        K = 50,100,200,400   vs committed M1C.1 A rows + M1B.2 annihilation t_stop
A5 K- acceptance       K = 100              （零反弹/严格下降/分解/判据/质量守恒 五项验收）
B1 KW @ K=100          B2 KS @ K=100
C  KW, KS @ K = 50,200,400（无论 K=100 结果如何固定执行）
```

**Baseline integrity：132/132 全部通过**（t_stop 324/699/1445/2876 与 D_max/U_density/moved 等字段逐字符串一致）。**K- 验收 5/5 通过**（t=3 exact stop、零反弹、最小降幅 −0.130、残差 3.1e-16）。因此无任何"暂停科学解释"情形。

## 5. Theory verification（5510 active 轮）

| 项 | 结果 |
|---|---|
| 分解恒等式 `D_new = D_TM_new + M·Δ`（对实际下一状态逐轮验证） | max residual **7.55×10⁻¹⁵** ≤ τ_num |
| 精确 rebound 判据（iff） | **TP=1983 / TN=3527 / FP=0 / FN=0** |
| numerical-boundary cases（|worst| ≤ 1e-12） | **0** |
| R_t（M1C.4 语义）> 1 ⟺ 判据反弹 | 逐轮等价，0 例外 |
| 质量守恒 | max 1.64×10⁻¹⁴ ≤ τ_num；min bin 7.4×10⁻⁶ > 0 |

TP=1983 = KC 的 1958（Stage 2 全部 canonical 反弹轮，复现）+ KS 的 25。**单步理论在新 kernels 上零误分泛化。**

## 6. Results（证据标签按 README 信息等级）

### 6.1 K0（G=T）— integrity reference

全部 K：t=6 exact stop、无反弹、moved ≈1.12——与 M1C.1 Variant U 一致（[Model-derived Mathematical Result]（M1C.2 定理 G=T ⟹ 严格收缩）+ [Verified Simulation Result]）。

### 6.2 K-（G=x²≤T）— negative control

K=100：t=3 exact stop（比 U 快：远偏 law 的 net export 1−a² > 1−a）、零反弹、moved 0.579。序定理预测与轨迹完全一致（[Model-derived Mathematical Result] + [Verified Simulation Result]）。**K- 通过是"实现/离散化/数值回归"的排除证明，不是科学发现。**

### 6.3 KW（λ=1/4）— W0 分支

全部 K：t=8 exact stop、**N_R=0**、R_t 单调爬升但 max 0.68–0.71 始终 <1、boundary 活跃（switch freq 0.571）但不兑现、D_max 严格单调降到 ~2×10⁻¹⁶、moved ≈1.42（K0 为 1.12：mismatch 使停止推迟 2 轮）。跨 K 几乎不动（t_stop 恒 8）。

**Frozen contract W0（[Verified Simulation Result] / [Model-dependent Observation]）**：

> positive mismatch existence alone is not sufficient to overcome the target-matched contraction margin along these trajectories.

**不得**写：positive mismatch never causes rebound（单点幅度、单初态）。

### 6.4 KS（flat-top）— S1 分支

| K | t_stop | N_R | f_R | S_R/K | rebound 轮 | max R_t |
|---|---|---|---|---|---|---|
| 50 | 26 | 5 | 0.192 | 0.220 | 13,15,21,22,24 | 2.56 |
| 100 | 26 | 5 | 0.192 | 0.110 | 13,15,21,22,24 | 2.79 |
| 200 | 26 | 5 | 0.192 | 0.055 | 13,15,21,22,24 | 2.60 |
| 400 | 29 | 10 | 0.345 | 0.048 | 8,12,13,17,18,20,22,23,25,27 | 5.20 |

K=100 签名：前 12 轮收缩（D_max 0.25→0.064，R_t 0.55→0.96 逼近 1）→ t=13/15 两次反弹 → 收缩 → t=21 大反弹（+3.7×10⁻²，R=2.79）→ t=22/24 两次边际反弹（R≈1.00–1.02）→ t=26 exact stop。rebound inter-arrival 中位 2 轮；longest contraction run 13；switch freq 0.80；**全部 rebound 轮 a ≥ 0.66**（R_t p50 ≈0.75，接近 canonical 的 0.69）。

**Frozen contract S1（[Verified Simulation Result] / [Model-dependent Observation]）**：

> one-step rebound mechanism survives the shape change, while long recurrent dynamics remains shape-sensitive.

**不得**写：flat top causes the difference（正则性混淆）；mismatch shape does not matter（单一 alternative shape）。

### 6.5 KC — canonical reproduction

t_stop 324/699/1445/2876、N_R 117/251/531/1059（合计 1958）、f_R 0.361–0.368、S_R/K 6.3–7.2、switch freq 0.92–0.94、moved 49–415——与 Stage 2 全部 committed 数字一致（[Verified Simulation Result]，仅作 anchor，未做新的 canonical 深挖）。

## 7. Cross-K / resolution sensitivity

- **类别跨 K 稳定**（KW=No rebound、KS=Short recurrent transient、KC=Long-span churn-like、K0/K-=No rebound，全部 4 K 一致）；**无任何 censoring**（最长轨迹 2876 < 10K）；**理论验证跨 K 稳定**（FP=FN=0，残差 ≤7.6×10⁻¹⁵）。
- **定量细节分辨率敏感**：KS 的 rebound 轮次集合在 K≤200 完全相同（{13,15,21,22,24}——近乎 resolution-invariant 的早期瞬态），K=400 增至 10 次；S_R/K 随 K 下降主要是有界 S_R 对增长 K 的标度效应。KW 的 t_stop=8 与 R_t 分位跨 K 不变。
- 第一轮冻结问题的回答：**两个新现象（KW 的无反弹快速停止、KS 的短促复发后快速停止）都 survives resolution changes；没有出现 resolution-specific 现象。**

## 8. Negative findings / falsification（完整保留）

1. **Churn-like transient 在两个 positive 对照上均未出现**——本轮结果**不支持 `positive mismatch alone broadly produces canonical-like churn` 这一强 generality claim**（同 peak、同 integral、同对称性的 flat-top mismatch（KS）只产生 19–34% 的短促反弹后快速停止；1/4 幅度（KW）连一次反弹都没有）。当前证据表明，canonical long-span churn 依赖于比正 mismatch 的存在、峰值和积分强度更细的 mismatch 结构；amplitude 与 spatial profile/regularity 均会显著改变观察到的 dynamics，但维持 long transient 的具体多步机制仍未解决。（措辞经 2026-09-20 interpretation audit 收紧；不得升级为"mismatch-driven churn 不是 model-class phenomenon"。）
2. **Rebound 本身也不是 positive mismatch 的必然伴随**（KW 零反弹）：mismatch 存在性 ≠ rebound；兑现还要求 reinjection 强度与逐点 γ 的对齐（S0 contract 提示的 Δ/γ pointwise alignment）。
3. canonical 的 `a_crit` 闭式门**未**在 KW/KS 上使用（canonical specialization）；仅 descriptive 记录：KS 全部 25 个 rebound 轮 a≥0.66，定性上与 canonical 的"大边界才反弹"观察一致（[Model-dependent Observation]，非 gate 验证）。
4. 无 finite-horizon censoring 观测（全部轨迹远早于 10K exact stop）——"10K 可能不够"这一顾虑本轮无数据支持。

## 9. Interpretation boundaries（禁止越界）

- 全部结果只对：canonical initial state、α=0.25、K∈{50,100,200,400}、H=10K、deterministic expected-mass path 成立。
- W0/S1 的措辞边界见 §6；分类 label（Short recurrent transient 等）是 frozen reporting convention（S_R/K≥1 分界是记录约定，非 claimed transition）。
- 禁止词汇（Stage 2 纪律延续）：non-convergent、infinite churn、attractor、chaos、ergodic、invariant distribution——本轮所有轨迹都 exact stop，研究对象始终是 finite-time transient。
- KS 的形状/正则性混淆未解除；λ 只有 1/4 一个非 canonical 点。

## 10. Unresolved questions（记录，不在本轮展开）

1. **什么泛函预测"长瞬态 law vs 快停止 law"？** KW（max R_t≈0.7）不穿阈、KS（R_t p50≈0.75、p95≈1.5）短促穿越、canonical（R_t 分布与 KS 接近！p50≈0.70）却持续数百轮——R_t 分布本身不区分 KS 与 KC，需要 Δ/γ 的 pointwise alignment 或 profile 演化层面的量。
2. **λ 临界问题**：λ=1/4 无反弹、λ=1（canonical）长 churn；1/4 与 1 之间是否存在 λ*（以及它是否依赖初态/K）——**possible follow-up question，仅记录，本轮禁止 sweep**。
3. KS 的 K=400 定量跳变（N_R 5→10）是网格采样 flat-top 边缘的离散效应还是趋势？
4. a≥0.66 的 KS 反弹边界观察能否得到类似 a_crit 的解析门槛（需要把 M1C.5 的 h-最大化推广到一般 G——**new theorem，属后续 gate 决策**）。

## 11. External Review Candidate Questions（本轮新增，供 Stage 3 closure 后评审）

| # | 问题 | 来源 |
|---|---|---|
| ERC-1 | 对固定 Δ 幅度/峰值的 law，"长瞬态 vs 快速 exact stop"是否存在 Δ(j) 与 γ(j) 逐点对齐的解析刻画（而非数值判据逐轮分类）？ | §8.2 + frozen S0 hint |
| ERC-2 | KS 的 shape isolation 数学上是否成立（flat-top 与分段线性正则性混淆）？是否需要 smooth-approximation 对照才能声称"shape matters"？ | §2 scope boundary |
| ERC-3 | 本轮 KS/KC 的 R_t 分布几乎相同（p50 0.75 vs 0.70）但 transient 长度差 25 倍——单步判据的"穿越频率"显然不是瞬态长度的充分统计量，长期动力学是否需要多步/几何论证？ | §10.1 |
| ERC-4 | 已知框架对应：一步分解 + margin/reinjection 竞争是否与 control theory 的 disturbance-rejection / ISS（input-to-state stability）或 discrepancy theory 已知结果同构？ | 文献定位待 Stage 3 closure 后统一做 |
| ERC-5 | rebound witness localization 与 evolving contraction margin γ_t(j) 是否存在低维或 operator-level 表示，可预测 multi-step recurrence（witness 每轮大幅迁移、无 persistent hot spot、S_R 与 stopping 的联合结构），而不仅是逐轮 rebound？ | 2026-09-20 interpretation audit（ERC-1 的轨迹级互补问题） |

## 12. Explicit next gate（冻结交接，本轮不做任何实现）

**Gate 1 已按 frozen stop condition 完成（PASS）**：K0/KC integrity ✓、K- @K=100 ✓、KW/KS @ {50,100,200,400} ✓、theory verification ✓、summary metrics ✓、minimal plots ✓、tests（306 项运行时断言，官方口径）✓、canonical report ✓、research-log ✓、next-step checkpoint ✓。

下一步（Owner + GPT 决策，Zcode 不得自行启动）：

1. **GPT/Owner interpretation of Gate 1**（本报告 + docs 交由其裁决 W0/S1 的采纳与 negative finding 的定位）；
2. **Gate 2 — Sign-changing mismatch** 的 kernel 冻结（Stage 3 原则中的 optional sign-changing；本轮未实现、未运行、只允许作为 follow-up 记录）；
3. 或按 §10 unresolved questions 决定是否插入 amplitude/shape 中间点（需要新的冻结 prompt）。

## 13. 资产索引

| 交付 | 路径 |
|---|---|
| 冻结配置 | `experiments/m2_gate1_generality/config.json` |
| Runner（一次运行复现全部） | `experiments/m2_gate1_generality/run_experiment.py` |
| 实验报告 | `experiments/m2_gate1_generality/M2_GATE1_REPORT.md` |
| 逐轮诊断（5510 行） | `experiments/m2_gate1_generality/results/per_round_diagnostics.csv` |
| 汇总表（Table 2） | `experiments/m2_gate1_generality/results/rebound_summary.csv` |
| 静态验证 / baseline integrity | 同目录 `static_kernel_validation.csv` / `baseline_integrity.csv` |
| 图 A/B/C + 对比图 | 同目录 `fig_gate1_*.png` |
| 新增 kernel 构造器 | `src/sand_m0/model.py`（`weak_positive_distribution`、`flat_top_distribution`） |
| 测试 | `tests/test_m2_gate1_generality.py`（18 项；全库 306 项官方口径通过） |
