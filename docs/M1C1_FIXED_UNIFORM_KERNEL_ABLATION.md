# M1C.1 诊断文档：Fixed-Uniform Kernel Ablation——把 near-bias 与 deficit-unawareness 拆开

> 状态：机制对照实验轮（2026-09-19）。研究问题由 GPT + Owner 冻结，是 M1C 的直接后续 refinement。
>
> 实验：[M1C1_REPORT.md](../experiments/m1c1_uniform_kernel_ablation/M1C1_REPORT.md)；前置：[M1C_CORRECTION_KERNEL_ABLATION.md](M1C_CORRECTION_KERNEL_ABLATION.md)、[M1B_MICRO_CORRECTION.md](M1B_MICRO_CORRECTION.md)
>
> 正式研究问题（冻结）：
>
> > **If the M1 selection rule remains unchanged and redistribution remains fixed and deficit-unaware, does removing the near-bias eliminate or substantially weaken churn?**
>
> 中文：如果完全保持 M1 的 selection rule，也仍然使用固定、完全不看当前 deficit 的 redistribution，只把 near-biased `q` 换成 Uniform `q`，churn 会不会消失或明显减弱？

## 0. 裁决

**Churn 消失（Outcome U1），且 M1C 的机制解释被精化。** Variant U（fixed uniform respray，仍然 fixed + state-unaware + deficit-unaware）在全部 4 个 canonical K 下 **t=6 exact stop**：D_max 与 U_density 逐轮严格单调下降（0 次回升），无 churn 波动期，无 ∝K 等待尾巴，total moved mass 1.12（A 为 49–415）。

同时出现两个必须如实报告的结构性事实：

1. **U 的停止态不是 bin 级均匀**：U_density = 5.25–5.55×10⁻³（跨 K 稳健），但 D_max = 0（精确）——controller 自己的目标（cumulative excess 清零）精确达成，而 bin 级 U 指标仍停在 A 的 churn 波动带下沿量级。只有 B（deficit-aware oracle）达到 bin 级机器均匀。
2. **U 的动作不衰减**：6 轮全部是全尺寸动作（moved 0.097–0.236，中位数 0.193 > A 的 0.150）；终止不是"动作缩到零"，而是第 6 轮的一次性大前缀消灭（a=0.94，F(a)=0.94 > 0.5——与 A 的小前缀消灭几何 F(a)≤0.5 **相反**）。

## 1. Research Question 与 M1C 为什么不足以隔离 bias

M1C 把 redistribution 从 "fixed near-biased q" 一次性换成 "adaptive + target-aware + deficit-aware oracle"，同时改变了多个性质，因此无法区分 churn 的来源是：

1. redistribution law 本身的 near-bias / target mismatch；
2. law 是 fixed、不看当前 deficit；
3. 两者共同作用。

M1C.1 增加唯一缺失的中间 control **Variant U**：仍然 fixed、state-unaware、deficit-unaware、每轮同一个 law——但这个 fixed law 恰好等于均匀目标。

## 2. Variant 定义

| Variant | Selection | 初始状态 | Redistribution law | 性质 |
|---|---|---|---|---|
| **A** | canonical cumulative-excess argmax（未动） | q_near bins（未动） | fixed `q_near`（未动） | fixed / target-mismatched / deficit-unaware |
| **U**（新） | **与 A 完全相同** | **与 A 完全相同（q_near）** | fixed `q_uniform[i] = 1/K`：`add_i = M_removed/K` | fixed / target-matched / deficit-unaware |
| **B** | 与 A/U 完全相同 | 与 A/U 相同 | M1C deficit-fill oracle | adaptive / target-matched / deficit-aware |

- **A vs U 唯一机制差异**：respray 概率从 `q_near`（近端密度 2× 均匀）换成 `q_uniform`（逐 bin 1/K）。selection、初始状态、alpha、stopping、指标全部逐位不变。
- **U vs B 唯一机制差异**：respray 是否读取当前 deficit。
- U 走**与 A 完全相同的 canonical fixed-throw code path**（`redistribution="throw"`），只通过新可选参数 `redistribution_dist=DISTRIBUTIONS["uniform"]` 替换 respray 概率（`dist` 仍提供初始状态）。默认 `None` 时行为与 M1C 完全位级一致。无新 simulator。

## 3. Controlled-variable logic

```text
A vs U：fixed + deficit-unaware + near-biased  vs  fixed + deficit-unaware + uniform
        → 隔离 near-bias / target mismatch
U vs B：fixed uniform + deficit-unaware  vs  adaptive deficit-aware oracle
        → 隔离 state / deficit awareness
```

注意：该因子表是机制解释的工具，**不**构成"三个性质可独立操纵"的证明。

## 4. Canonical configuration（完全继承 Stage 1 / M1C）

K∈{50,100,200,400}（主 K=100）、alpha=0.25、exact tol=1e-12、T_max=5000、初始状态 q_near、同一指标与 churn 诊断、eps overlay {0.005,0.01,0.02,0.04}、deterministic exact reference path。未新增 K，未做 sweep，config.json 运行前冻结，无 post-hoc 调整。

## 5. Baseline integrity（Stage 1 门禁）

A、B 在冻结配置下重跑，与已提交 M1C `variant_summary.csv` 做 16 字段 × 8 行（4 K × 2 变体）字符串级对照，另加 A 的 4 个 t_stop 对照 M1B.2 annihilation CSV：**132/132 全部通过**。历史结果无 drift，A/U/B 比较有效。

## 6. 数学结果（Model-derived Mathematical Result）

1. **Uniform respray 的 D-profile 恒等式**：`p_next = p_minus + M_removed/K`（逐 bin 常数），故
   `D_{t+1}(j) = D_minus(j) + (j/K)·M_removed`（对一切 j；j=K 处 −M_removed + M_removed = 0 精确回零）。uniform respray 对 D-profile 加一个**随 j 线性上升的 tilt**——这解释了 U 的 argmax 向远端前缀漂移（a_t: 0.50→0.59→0.69）与大前缀消灭的自然可达性。
2. **一般 fixed respray law 的 net-export 恒等式（M1B.2 恒等式的推广）**：对任意固定 respray law `g`，
   `prefix_mass_change = −moved·(1 − F_g(a_t))`。A（g=q_near）与 U（g=q_uniform，F_g(a)=a）逐轮数值核验，残差 ≤ 1.0×10⁻¹⁴。
3. **U 的更新不变量**：每轮加入质量 = `M_removed/K`（精确常数向量，|added − M/K| ≤ 2.5×10⁻¹⁷），总和 = M_removed（≤ 5.6×10⁻¹⁷）；质量守恒 |Σ−1| ≤ 2.2×10⁻¹⁶；无非负违规。
4. **Selection invariance**：同状态 → 同 decision（A/U t=0 行逐位相同；U 自身轨迹的边界与 argmax 规则逐轮一致）。

## 7. 主结果

### 7.1 三变体总表（4 个 canonical K）

| metric | A | U | B |
|---|---|---|---|
| t_exact_stop | 324 / 699 / 1445 / 2876（∝K） | **6 / 6 / 6 / 6** | 3 / 3 / 3 / 3 |
| t_U_best | 292 / 653 / 1349 / 2627 | 6（= stop） | 3（= stop） |
| U_density@stop | 8.99 / 7.07 / 7.88 / 7.14 ×10⁻³ | **5.55 / 5.42 / 5.25 / 5.41 ×10⁻³** | 10⁻³¹ – 10⁻³³（机器均匀） |
| U_density@best | 3.71 / 3.82 / 3.69 / 3.56 ×10⁻³ | = stop | = stop |
| retention@stop | 0.984–0.990 | 1.000000 | 1.000000 |
| total moved mass | 49.0 / 103.3 / 210.6 / 415.3 | **1.12（跨 K 几乎相同）** | 0.32 |
| moved 中位数/轮 | 0.147–0.154 | 0.193 | 0.104 |
| late band（D_max<0.005）moved | 0.231–0.238 | 0.2358（仅 1 轮 = 末轮消灭） | n/a（带为空） |
| D_max 上升轮占比 | 0.36 | **0** | 0 |
| U 上升轮占比 | 0.43–0.46 | **0** | 0 |
| near zone @stop | 0.477–0.493 | 0.4744–0.4750 | 0.5000（精确） |

### 7.2 U 的 6 轮轨迹（K=100；4 个 K 结构相同）

| t | a_t | D_max | moved | U_density |
|---:|---:|---:|---:|---:|
| 0 | 0.50 | 0.250000 | 0.187500 | 3.333×10⁻¹ |
| 1 | 0.59 | 0.165025 | 0.188756 | 1.627×10⁻¹ |
| 2 | 0.69 | 0.097261 | 0.196815 | 7.000×10⁻² |
| 3 | 0.34 | 0.049016 | 0.097254 | 2.498×10⁻² |
| 4 | 0.84 | 0.027148 | 0.216787 | 1.428×10⁻² |
| 5 | 0.94 | 0.003173 | 0.235793 | 5.968×10⁻³ |
| 6 | **stop**：D_max = −2.2×10⁻¹⁶ ≤ tol；U_density = 5.42×10⁻³ | | | |

停止态结构：全部前缀累计质量 ≤ 均匀份额（D_max 精确 ≤ 0），但 bin 级偏差达 −0.123u…+0.233u（40 个 bin 高于均匀、60 个低于）——一种"波动但从不累计超出"的构型。**controller 的目标（D_max）精确清零，bin 级均匀度（U）仍有限**。这是 A（stop 态 U_density 7.1×10⁻³，同样 D_max=0）与 U 共有的 stopping 结构；B 的饱和填充是唯一同时清零两者的机制。

### 7.3 Epsilon overlay（K=100）

| variant | eps=0.005 | 0.01 | 0.02 | 0.04 |
|---|---|---|---|---|
| A | t=340（ret 0.986） | t=159 | t=80 | t=35 |
| U | **t=5**（U=5.97×10⁻³，ret 0.9983） | t=5 | t=5 | **t=4**（U=1.43×10⁻²，ret 0.973） |
| B | t=3（ret 1.0） | t=2 | t=2 | t=2 |

对 A，有限 eps 的意义是"提前退出 churn 窗口"（340 vs 699）。对 U，**deadband 不再具有该意义**——没有 churn 窗口可退，eps stop 只比 exact stop 早 1–2 轮、停在几乎同样质量的状态。U 的 finite-eps 语义退化为"快速收敛过程上的提前 1–2 步截止"。

## 8. Q1–Q8 逐问回答（U）

| 问题 | 答案 | 依据 |
|---|---|---|
| Q1 early improvement → sustained churn → delayed stopping？ | **只有第一段**。改善 0.333→2.5×10⁻²（3 轮）；无 churn 期；t=6 exact stop | §7.2 |
| Q2 D_max 长期震荡？ | **否**——严格单调下降，0 次上升，6 轮后终止 | share_D_max_increases = 0 |
| Q3 uniformity 反复改善/恶化？ | **否**——U_density 逐轮严格下降，0 次回升 | share_delta_U_positive = 0 |
| Q4 late-stage moved mass substantial 还是衰减？ | **不适用/不衰减**：无 late stage；6 轮全部全尺寸（0.097–0.236），末轮消灭动作 0.236 是全程最大——终止靠消灭，不靠动作收缩 | §7.1–7.2 |
| Q5 best state 明显早于 stop？ | **否**——best = stop = 6 | §7.1 |
| Q6 t_stop 随 K 增长？ | **否**——6/6/6/6（K 无关；Model-dependent Observation，4 点） | §7.1 |
| Q7 finite-eps 仍有"提前退出 churn"意义？ | **否**——eps stop 只早 1–2 轮，语义退化为快速收敛上的提前截止 | §7.3 |
| Q8 U 定性上更接近 A 还是 B？ | **宏观动力学像 B**（快速单调收敛、无 churn、stop 跨 K 恒定）；**stopping 机制像 A**（一次性 cumulative-excess 消灭、停止态 D_max=0 而 U 有限），但消灭几何相反（大前缀 F(a)=0.94 vs A 的小前缀 F(a)≤0.5） | §7.2、§9 |

## 9. Mechanism interpretation（A 为什么 churn、U 为什么不 churn）

- **A 的 churn 来源**：q_near 在近端的回填密度是均匀的 2 倍，每轮系统性再造它刚移走的近端过剩——D_max 无法持续下降，进入平稳波动；且 q_near 的回填几何使消灭事件必须用小前缀扫掠（M1B.2 必要条件 F(a)≤0.5），该构型在波动中低频出现 → 等待 ∝K。
- **U 为什么不 churn**：uniform respray 对每个前缀恰好加入其均匀份额（§6 恒等式 1 的 tilt），不再再造超出累计结构所能吸收的过剩——D_max 单调下降到零。近端再创造的消失同时使**大前缀消灭**（F(a)=0.94）成为自然构型，6 轮内达成，K 无关。
- **三变体机制分解**：
  - **A vs U**：churn 的有无由 near-bias / target mismatch 决定 → **near-bias 是 canonical M1 churn 的关键驱动**（本轮 setup 内）；
  - **U vs B**：deficit-awareness 不改变"是否 churn"（U 已无 churn），它改变的是**收敛的精度与停止态的结构**——U 停在 D_max=0 的"累计均匀"锯齿态（U_density ≈ 5.4×10⁻³），B 达到 bin 级机器均匀（10⁻³¹）。deficit 信息买到的是 bin 级精度，不是 churn 的消除。

## 10. Evidence Classification

- **Model-derived Mathematical Result**：§6 恒等式 1–3（D-profile tilt、一般 fixed-law net-export 恒等式、U 常数加入向量）——由更新规则严格推出并数值核验。
- **Verified Simulation Result**：integrity 132/132；A/U/B 全部轨迹量（t_stop、t_best、U@stop/U@best、retention、moved 统计、带分位数、eps overlay、U 不变量实测界）。
- **Model-dependent Observation**：U 无 churn 且 t=6 跨 K 恒定（4 点）；U 停止态的 D_max=0 / U≈5.4×10⁻³ 分离结构；大前缀消灭几何；deadband 语义变化；argmax 向远端漂移。
- **Working Hypothesis（本轮后更新）**：见 §12——near-bias/target mismatch 是 churn 的关键驱动（supported by this experiment）；"deficit-unawareness 本身是 churn 主因"的强解读被**削弱**。

## 11. Limitations / What This Does NOT Establish

1. 只测了一个 uniform control（恰好等于目标的 fixed law）；未测任意其他 target-matched fixed law（如远偏、其他幂族 p≠1）；
2. 单一 q_near 初始状态、单一 alpha=0.25、一维、deterministic path（MC 未实现 U）、4 个 K；
3. t=6 与 K 无关是 4 点观察，未证明；U 的单调性是实测（6 轮），未证明一般性；
4. U 停止态的锯齿结构未做稳定性/再现性分析（它是否是 attractor 结构未知）；
5. 不证明 near-bias 是 churn 的**唯一**来源（只证明：去掉它，churn 在此 setup 消失）；
6. 未做任何敏感性/鲁棒性扫描（Prompt §11 禁止）。

## 12. Updated Hypothesis Status 与 Next Decision Gate

**对核心问题（§0）的回答**：在 selection rule 不变、redistribution 仍 fixed 且 deficit-unaware 的条件下，**去掉 near-bias 本身就消除了 churn**（不是减弱——是消失）。

**假设状态更新**：

- M1C 后的解释："churn 与 fixed, deficit-unaware redistribution 有重要关系"——其**宽泛解读（deficit-unawareness 本身）被本轮削弱**：U fixed + deficit-unaware 却无 churn。
- **被支持的新解释**：churn 的关键承载结构是 **fixed redistribution law 与目标之间的 near-bias / target mismatch**；deficit-awareness 的作用是 bin 级精度（U：5.4×10⁻³ vs B：10⁻³¹），不是 churn 的消除。
- 以上均为 Model-dependent Observation / Working Hypothesis 级（single uniform control、deterministic path、4 K）。

**Decision gate（交回 GPT + Owner）**：

1. **target-matched fixed law 的一般性**：U 的无 churn 是否对其他 target-matched fixed law（p≠1 幂族、远偏）成立？（新对照轮，待冻结）
2. **A 侧旧问题保留**：exact-stop 等待 ∝K 解析化、连续极限、churn regime 长期动力学——A 机制未变；
3. **停止态结构**：U 的 D_max=0 / U 有限锯齿态与 A 的 annihilable 构型、B 的饱和态三者关系值得理论化（controller 目标与 U 指标的分离）；
4. 依据 Prompt §17，本轮在 A/U/B canonical comparison 完成后停止，不进入 interpolation/adaptive-q/noise/2D。

## 13. Repository 资产

- 实现：`src/sand_m0/adaptive.py`（`run_m1a_deterministic` 新增可选 `redistribution_dist`；默认 None 位级不变；无新 simulator）；
- 测试：`tests/test_m1c1_uniform_kernel.py`（19 项）；
- 实验：`experiments/m1c1_uniform_kernel_ablation/`（config.json / run_experiment.py / results / M1C1_REPORT.md）；
- 结果：`baseline_integrity.csv`（132 项）、`per_round_diagnostics.csv`（A/U/B 逐轮）、`variant_summary.csv`、`epsilon_overlay.csv`、两张图、metadata.json。
