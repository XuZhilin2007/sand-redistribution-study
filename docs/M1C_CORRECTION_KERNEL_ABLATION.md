# M1C 诊断文档：Correction-Kernel Ablation——固定 selection、只换 redistribution kernel 的单机制对照

> 状态：机制对照实验轮（2026-09-19）。研究问题由 GPT + Owner 在 Focused Literature Positioning 之后冻结。
>
> 实验：[M1C_REPORT.md](../experiments/m1c_kernel_ablation/M1C_REPORT.md)；前置：[M1B_MICRO_CORRECTION.md](M1B_MICRO_CORRECTION.md)（churn 机制）、[M1A_LONG_HORIZON.md](M1A_LONG_HORIZON.md)、[M1B_TOLERANCE.md](M1B_TOLERANCE.md)、[M1B_METRIC_SCALE_AUDIT.md](M1B_METRIC_SCALE_AUDIT.md)
>
> 正式研究问题（冻结）：
>
> > **Does churn persist when the M1 selection rule is held fixed but the fixed biased redistribution law is replaced by a target-aware corrective kernel?**
>
> 待检验 Working Hypothesis（文献定位后冻结）：M1 后期 persistent churn 的来源，是 controller 虽能识别"哪里需要处理"，但被扫回的质量只能经过一个固定、外生、带偏的 redistribution law `q` 重新进入系统，无法针对当前 deficit 做纠正。

## 0. 裁决

**Churn 完全消失（不是减轻）。** 在 canonical M1 configuration（q_near、alpha=0.25、exact stopping、K∈{50,100,200,400}、T_max=5000）下：

- **Variant A**（canonical baseline，redistribution = fixed biased `q_near`）：完整复现 M1B.2 的全部 churn 行为（baseline reproduction 112/112 逐位通过）——exact stop 在 t=324/699/1445/2876（∝K），churn 每轮搬动 ~15% 总质量，D_max 在 0.003–0.08 带内震荡数百至数千轮。
- **Variant B**（selection rule 逐位不变，仅把 `q` 换成 target-aware deficit-filling kernel）：**3 轮内达到机器精度均匀并触发 exact stopping**（全部 4 个 K 都是 t=3；停止态 U_density ≈ 10⁻³¹–10⁻³³，max|p−u| ≈ 10⁻¹⁷），D_max 与 U_density 逐轮严格下降，总搬动质量 0.32（A 为 49–415，相差 150–1300 倍），无任何震荡期。

**结论等级**：本实验**支持**（supported by this experiment）当前 Working Hypothesis——在当前 canonical setup 中，"redistribution 与当前 deficit 脱钩"是 M1 persistent churn 的关键承载结构之一；把它替换成 deficit-aware kernel 后 churn 不复存在。**不**声称 fixed biased `q` 是 churn 的唯一原因、不声称该结果推广到其他 q/alpha/模型族、不声称 deficit-filling 是现实可行 controller（它是 oracle 式 ablation control，见 §12）。

## 1. Research Question 与动机

M1B.2 确立了当前机制解释：M1A 的 exact controller 在快速宏观清理（~前 100 轮）后进入**平稳 churn 波动 regime**（每轮搬动 ~15% 总质量、D_max 在 K 稳健带内不规则震荡、永不收敛），exact stopping 靠低频进入"annihilable stopping configuration"（等待轮数 ∝K）。

Focused Literature Positioning 之后，当前 Working Hypothesis 把 churn 的嫌疑集中在 redistribution law 上：cumulative-excess selection 能定位"哪里多余"，但被移除的质量按固定带偏的 `q` 重泼，与当前 deficit 完全脱钩——q_near 在近端的回填密度是均匀的 2 倍，会系统性再造它刚移走的近端过剩。

M1C 只检验一个最小问题：**保持 selection rule 逐位不变，只替换 redistribution kernel，churn 是否仍然存在？** 这是 falsification-oriented 的单变量机制对照，不是参数扫描。

## 2. Variant 定义

### Variant A — canonical M1 baseline（未动）

`run_m1a_deterministic(..., redistribution="throw")`（默认值），即当前 M1A/M1B/M1B.2 canonical 实现：

1. `D_t(j) =` 前缀累计质量 − `j/K`；`j_t = argmax_j D_t(j)`；
2. 若 `D_max > tol`：前缀 `1..j_t` 每bin移除 `alpha` 比例，`M_removed = alpha·F_t(j_t)`；
3. 移除质量按 fixed `q_near` 重泼：`p_{t+1} = p_minus + M_removed·q`；
4. `D_max ≤ tol` 停止，状态冻结。

### Variant B — target-aware corrective redistribution（唯一改变点）

`run_m1a_deterministic(..., redistribution="deficit_fill")`。第 1、2、4 步逐位不变，仅第 3 步替换为（prompt 冻结定义；目标分布 `u_i = 1/K` 即 selection rule 与 U 指标已编码的均匀目标）：

```text
deficit_i = max(u_i − p_minus_i, 0)
D_def     = Σ_i deficit_i
add_i     = M_removed · deficit_i / D_def      （D_def > 0 时）
p_next_i  = p_minus_i + add_i
```

这是一个 **deficit-filling oracle kernel**：它精确看到每个 bin 的当前 deficit，并把被移除质量按 deficit 比例放回。它不是现实沙坑模型，不是任何文献中的原始 controller，只是 mechanism ablation control：给 controller 一个能看到 deficit 并针对性分配的理想化纠正机制。

实现为 `src/sand_m0/adaptive.py::deficit_fill_redistribution` + `run_m1a_deterministic` 的新关键字参数（默认 `"throw"` 保持 canonical 行为位级不变；MC 路径未改，本 ablation 全部在 deterministic exact reference path 上——Stage 1 churn 证据都在该路径）。

## 3. Invariants（全部成立并经测试覆盖）

| Invariant | 状态 |
|---|---|
| Mass conservation：`sum(p_next) = sum(p_current)` | 成立（实测 |Σ−1| ≤ 4.4×10⁻¹⁶，全部 K 全部轮） |
| Non-negativity | 成立（fill 只加不减；最小 bin 质量 1.7×10⁻³ 量级，K=400） |
| Same selection：与 A 共用 canonical decision logic | 成立（同状态同 decision row；B 自身轨迹的记录边界与 argmax 规则逐轮一致） |
| No local overfill：fill 不把任何 deficit bin 填到 target 之上 | 成立（由 §5 恒等式 `add_i ≤ deficit_i` 保证；测试验证 `p_next ≤ max(p_minus, u)`） |
| Degenerate fallback：`D_def = 0` 或 `M_removed > D_def` 必须显式失败、不得静默回退到 q | 成立（kernel 内显式 raise；全部运行中从未触发——它们在数学上不可能，见 §5） |

## 4. Canonical configuration（继承 Stage 1，未新增）

| 项 | 值 |
|---|---|
| 模型 / dynamics | M1A deterministic exact reference（未改） |
| q / alpha / tol / T_max | `q_near`；alpha=0.25；exact tol=1e-12；T_max=5000 |
| K | 50 / 100 / 200 / 400（主 canonical K=100 做完整轨迹对照；其余 3 个为 Stage 1 已有 diagnostic resolutions，用于排除单 K 偶然） |
| 有限 deadband overlay | epsilon ∈ {0.005, 0.01, 0.02, 0.04}（M1B 冻结集） |
| 指标 | U_density = K·U_L2；D_max；moved_mass；L1_state_step；delta_U_density；improvement retention（M1B.1/M1B.2 原指标，未新造指标族） |
| 冻结声明 | config.json 运行前写入；无 post-hoc 调整（§15 遵守） |

## 5. 数学结果（Model-derived Mathematical Result）

**恒等式（deficit 总量 = 移除后正过剩 + 移除质量）**：记 `d_i = max(u_i − p_minus_i, 0)`、`e_i = max(p_minus_i − u_i, 0)`，则 `u_i − p_minus_i = d_i − e_i`。逐项求和并注意 `Σp_minus = 1 − M_removed`：

```text
D_def − E_minus = Σ(u_i − p_minus_i) = 1 − (1 − M_removed) = M_removed
⟹  D_def = E_minus + M_removed ≥ M_removed > 0   （ whenever M_removed > 0 ）
```

推论（全部由更新规则严格推出并数值核验）：

1. **Degenerate case 不可能**：`D_def ≥ M_removed > 0`（A/B 移除规则下 `M_removed = alpha·prefix_mass > 0`），故 Variant B 不存在"deficit 不够放"或"无处可放"的情形；实现中的两个 raise 只是 bug 检测器，从未触发。
2. **Fill 永不过填**：`add_i = M_removed·deficit_i/D_def ≤ deficit_i`，故 deficit bin 填后 `≤ u_i`；高于 target 的 bin 不接收质量。
3. **饱和情形 = 一步精确均匀**：若某轮扫掠后 `E_minus = 0`（无任何 bin 高于均匀），则 `D_def = M_removed` 精确成立，fill 把每个 deficit bin **恰好**填到 `u_i` → 状态精确均匀 → `D_max(t+1) = 0` → 下一轮 exact stop。这给出 Variant B 自己的 stopping configuration 概念，与 M1B.2 的 annihilable configuration 对应，但在 B 中它是被 kernel **快速构造**出来的，不是等待出来的。

## 6. Baseline reproduction（Variant A，Stage 1 门禁）

对 4 个 canonical K 重跑 canonical dynamics，与已提交 M1B.2/M1A.1 结果做 112 项逐位对照（字符串级，含 6dp/5dp 舍入字段）：exact stop 轮、pre-stop 边界与 moved、late-stage（D_max<0.005 带）moved/net_export/L1 中位数、D_max/U_density 波动带分位数、dip 频次与首达、16 个 (K,eps) first-passage 与 retention——**112/112 全部通过**（`results/baseline_reproduction.csv`）。

即：Variant A 成功复现 M1B.2 Stage 1 的核心 churn behavior（early rapid improvement → churn regime → substantial late-stage redistribution → D_max/U 震荡 → delayed exact stopping ∝K → finite-eps 在普通波动谷停手）。基线未漂移，A-vs-B 比较有效。

## 7. 结果

### 7.1 Variant A（canonical churn，对照）

| K | t_stop | t_U_best | U_density@stop | retention@stop | moved 中位数 | moved 总量 | D_max 带内上升轮占比 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 324 | 292 | 9.00×10⁻³ | 0.9840 | 0.1544 | 49.02 | 0.362 |
| 100 | 699 | 653 | 7.07×10⁻³ | 0.9901 | 0.1502 | 103.30 | 0.360 |
| 200 | 1445 | 1349 | 7.88×10⁻³ | 0.9873 | 0.1480 | 210.59 | 0.368 |
| 400 | 2876 | 2627 | 7.14×10⁻³ | 0.9892 | 0.1466 | 415.27 | 0.368 |

A 的 churn 特征全部在场：moved_mass 全程 ~0.15/轮（late band 0.231–0.238）；D_max 带内 36% 的轮次上升、U_density 43–46% 的轮次上升（震荡而非收敛）；best state 比 stop 早 32–249 轮；stop 态 U_density（7–9×10⁻³）比 best（3.6–3.8×10⁻³）差近 2 倍。

### 7.2 Variant B（deficit-fill kernel）

**全部 4 个 K：t_stop = 3，停止态机器精度均匀。**

| K | t_stop | t_U_best | U_density@stop | retention | moved 总量 | D_max 上升轮占比 | U 上升轮占比 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 3 | 3 | 3.6×10⁻³³ | 1.000000 | 0.3251 | 0 | 0 |
| 100 | 3 | 3 | 1.8×10⁻³¹ | 1.000000 | 0.3201 | 0 | 0 |
| 200 | 3 | 3 | 3.3×10⁻³³ | 1.000000 | 0.3213 | 0 | 0 |
| 400 | 3 | 3 | 1.6×10⁻³³ | 1.000000 | 0.3207 | 0 | 0 |

K=100 逐轮轨迹（4 个 K 几乎相同，只差格点尺度）：

| t | a_t | D_max | moved | U_density | D_def | E_minus | fill ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.50 | 0.250000 | 0.187500 | 3.33×10⁻¹ | 0.270825 | 0.083325 | 0.6923 |
| 1 | 0.33 | 0.083325 | 0.103331 | 4.39×10⁻² | 0.110275 | 0.006944 | 0.9370 |
| 2 | 0.11 | 0.006944 | 0.029236 | 6.59×10⁻⁴ | 0.029236 | 0（精确） | 1.0000 |
| 3 | — | 8.9×10⁻¹⁶ ≤ tol → **stop** | 0 | 1.84×10⁻³¹ | — | — | — |

（停止态 max|p_i − u_i| = 1.6×10⁻¹⁷；near zone 质量 = 0.5000 精确，A 为 0.477–0.493。）

**机制（数据确立的收敛路径）**：

- fill ratio 逐轮上升 0.692 → 0.937 → 1.000：被移除质量相对当前 deficit 的比例逐轮增大，因为 q 的"再造过剩"效应被消除后，正过剩单调集中于越来越短的近端前缀；
- 高于均匀的 bin 数按 **⌈K/3⌉（state 1）→ ⌈K/9⌉（state 2）→ 0（state 3）** 收缩（17/6、33/11、67/22、133/44 实测）；
- 最后一扫的条件：state 2 时最高于均匀 bin 的质量 0.0112 已低于 `u/(1−alpha) = 4u/3 = 0.01333`（state 1 时 0.0149 高于该阈值），故 t=2 的一扫把全部剩余正过剩压到均匀以下 → `E_minus = 0` → §5 推论 3 的饱和 → 精确均匀 → exact stop；
- **与 M1B.2 的连接**：A 中 exact stop 需要等待低频 annihilable configuration（∝K 轮）；B 中 deficit-aware kernel 主动、快速（3 轮）构造出该构型的饱和版。等待尾巴整个消失。

### 7.3 Epsilon overlay（有限 deadband stopping）

| variant | eps=0.005 | eps=0.01 | eps=0.02 | eps=0.04 |
|---|---|---|---|---|
| A（4 K） | t=253/340/338/321，retention 0.982–0.991 | t=165/159/173/169 | t=76/80/76/82 | t=34/35/34/34 |
| B（4 K） | **t=3**（=exact stop，retention 1.0） | **t=2**（U=6.59×10⁻⁴，retention 0.9980） | t=2 | t=2 |

A 的 first-passage 与已提交 M1B 值逐位一致；B 的停止轮跨 K 完全稳定（2–3 轮），停止态 U_density 6.6×10⁻⁴ 已低于 A 的整个活跃期波动带下沿（p5 ≈ 5.6×10⁻³）一个量级。M1B 的"在普通波动谷停手"机制在 B 中退化为平凡情形：谷在 2 轮内就深到 machine-precision。

## 8. A vs B 判读（对照 Prompt §13 的检查单）

| 检查项 | Variant A | Variant B |
|---|---|---|
| early rapid improvement | ✓（U_density 0.333 → ~10⁻²，前 ~100 轮） | ✓（0.333 → 6.6×10⁻⁴，2 轮，更陡） |
| substantial late-stage redistribution | ✓（moved ~0.15/轮直至 stop） | **✗**（3 轮后无任何动作） |
| oscillatory / churn regime | ✓（D_max 带 0.003–0.08，36% 轮上升） | **✗**（D_max、U_density 逐轮严格下降，0 轮上升） |
| delayed exact stopping | ✓（t=324–2876，∝K） | **✗**（t=3，跨 K 恒定） |
| uniformity 更接近单调改善 | — | ✓（严格单调下降，无一次回升） |
| D_max 更接近单调下降 | — | ✓（严格单调下降） |
| moved mass 衰减 | ✗（恒定 ~0.15） | ✓（0.188 → 0.103 → 0.029 → 0） |
| best state vs stopping state | 错开 32–249 轮 | **重合**（best = stop = 3） |
| exact / finite-deadband stopping 结构变化 | — | ✓（∝K 等待尾巴消失；eps stopping 提前至 2–3 轮且跨 K 恒定） |

注意（诚实记录）：B 的"单调"是本实验 4 条轨迹上的实测事实（Verified Simulation Result），不是推导出的定理——本文档不声称 deficit-fill kernel 在一切状态下必然单调（见 §11 未证明项）。

## 9. Interpretation（对照 Prompt §14 的结果树）

结果属于 **Outcome 1（B 基本消除 churn）的最强形式**：A 明显 churn、B 转为严格单调收敛并在 3 轮内 exact stop、late-stage moved mass 归零、best time 与 stopping 重合、stopping 的 ∝K 尾巴消失。

允许的表述：

> 在当前 canonical M1 setup（q_near、alpha=0.25、cumulative-excess selection、exact/finite-deadband stopping、deterministic 路径、K∈{50,100,200,400}）中，仅把 redistribution kernel 从 fixed biased `q` 换成 target-aware deficit-filling 后，persistent churn 完全消失。**该实验支持当前 Working Hypothesis：M1 churn 与"redistribution 与当前 deficit 脱钩"这一结构有重要关系。**

不允许的表述（Overclaim 禁令）：

- ✗ "已证明 fixed biased q 是 churn 的唯一原因"——本实验只确立了一个方向的机制证据：移除该结构即移除 churn（在此 setup 内）；未检验其他 churn 源是否可独立产生 churn；
- ✗ "deficit-filling 是可行/最优 controller"——B 是 oracle（精确知道 u 与每个 bin 的 deficit），现实操作者没有该信息；
- ✗ 任何关于其他 q / alpha / 2D / 随机环境的推广。

## 10. Evidence Classification

- **Model-derived Mathematical Result**：§5 恒等式 `D_def = E_minus + M_removed` 及三条推论（degenerate 不可能、fill 永不过填、E_minus=0 ⟹ 一步精确均匀）——由更新规则严格推出并数值核验（identity 残差 ≤ 1.8×10⁻¹⁶）。
- **Verified Simulation Result**：baseline reproduction 112/112；A/B 全部轨迹量（t_stop、t_best、U_density@stop、retention、moved 统计、带分位数、eps overlay）；B 的 3 轮收敛路径与 fill ratio 序列；invariants 全部实测界（|Σ−1|≤4.4×10⁻¹⁶ 等）。
- **Model-dependent Observation**：churn 在 B 下完全消失；B 的 exact stop 跨 K 恒定（t=3，4 个 K）；⌈K/3⌉→⌈K/9⌉→0 的收缩模式；fill ratio 单调升至 1；above-u bin 质量 vs `u/(1−alpha)` 阈值的穿越顺序。
- **Working Hypothesis（本轮后状态更新）**：M1 persistent churn 与 deficit-unaware redistribution 有重要关系——由本轮**升级为"supported by this experiment"**；其推广性（其他 q/alpha、部分感知、噪声）仍为开放问题。
- 无新文献声明；无 Real-world calibrated result。

## 11. Limitations / What This Does NOT Establish

1. **Oracle 信息**：B 精确知道 `u_i = 1/K` 与每个 bin 的 deficit。结果只说明"给 controller 完整 deficit 信息时 churn 消失"，不说明部分感知/噪声感知下消失；
2. 单一 q（q_near）、单一 alpha（0.25）、一维、deterministic path（MC 未实现/未测）、4 个 K；
3. t=3 跨 K 恒定是 4 点观察，未推导；"⌈K/3⌉→⌈K/9⌉"收缩与 `u/(1−alpha)` 穿越的轮数未证明（未尝试一般 K 的归纳）；
4. B 的单调性是实测，未证明；不排除在其他参数点出现非单调；
5. 不证明 churn 的唯一原因（见 §9）；不涉及 exact stop 等待 ∝K 的解析证明（A 侧机制未变）；
6. 未做任何敏感性/鲁棒性扫描（Prompt §7 禁止）——所有"是否稳健"问题留待 Owner + GPT 决策。

## 12. 与文献定位语境的关系

本实验回应的是文献定位后冻结的 Working Hypothesis（rank-driven selection / online thinning / selection–mutation particle systems / density-feedback control / deadband control 的定位背景见 RESEARCH_CHECKPOINT_2026-09-19 §6/§9）。本轮**不新增**文献声明、不扩张 novelty claim；M1C 在该语境中的角色是：把"redistribution law 是否 load-bearing"从猜想变成一个受控对照的数据点，供后续定位与理论化使用。

## 13. Repository 资产

- 实现：`src/sand_m0/adaptive.py`（新增 `deficit_fill_redistribution`；`run_m1a_deterministic` 新增 `redistribution` 参数，默认不变）；
- 测试：`tests/test_m1c_kernel_ablation.py`（20 项；kernel 单元、invariants、shared selection、reproducibility、committed anchors）；
- 实验：`experiments/m1c_kernel_ablation/`（config.json / run_experiment.py / results / M1C_REPORT.md；一条命令复现）；
- 结果文件：`baseline_reproduction.csv`（112 项）、`per_round_diagnostics.csv`（A+B 全部活动轮）、`variant_summary.csv`、`epsilon_overlay.csv`、两张图、metadata.json（SHA-256）。

## 14. Decision Gate / Unresolved Questions（交回 GPT + Owner）

1. **推广性 gate**：deficit-aware kernel 的结论是否跨 q（q_far/q_uniform/幂族）、跨 alpha、跨 K 更细网格成立？（本轮禁止扫描，未做。）
2. **部分感知版本**：把 oracle deficit 换成有噪/粗化/有限分辨率的 deficit 估计，churn 在多大感知误差下重新出现？——这直接连接 perception model 方向与 deadband 文献；
3. **理论化**：t=3 收敛与 ⌈K/3⌉→⌈K/9⌉ 收缩能否对一般 K 证明？B 的 stopping time 是否有闭式？
4. **A 侧既有开放问题不变**：exact stop 等待 ∝K 的解析化、连续极限、churn regime 的长期动力学；
5. 依据 Prompt §19，本轮在完成 A vs B 后停止，不进入 Variant C/D、不引入 mixture/interpolation、不做 sweep。
