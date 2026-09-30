# M1C.1 实验报告：Fixed-Uniform Kernel Ablation（A vs U vs B）

> 日期：2026-09-19
>
> 分析文档：[docs/M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md](../../docs/M1C1_FIXED_UNIFORM_KERNEL_ABLATION.md)（研究问题、controlled-variable logic、数学结果、Q1–Q8 判读、假设更新）；前置：[M1C_REPORT.md](../m1c_kernel_ablation/M1C_REPORT.md)
>
> 代码：`src/sand_m0/adaptive.py`（`run_m1a_deterministic` 新增可选 `redistribution_dist`；默认 `None` 与 M1C 行为位级一致；U 走与 A 完全相同的 fixed-throw code path）；配置：[config.json](config.json)（运行前冻结）
>
> 性质：单变量机制对照。Variant U 是唯一新机制：fixed、state-unaware、deficit-unaware，但 respray law 恰为均匀目标（`q_uniform[i] = 1/K`）。

## 1. 裁决（headlines）

| 问题 | 答案 |
|---|---|
| A/B 是否与 M1C 提交结果一致（无 drift）？ | **是**：132/132 字符串级对照通过（16 字段 × 4 K × 2 变体 + A 的 4 个 t_stop vs M1B.2） |
| 去掉 near-bias（保持 fixed + deficit-unaware）后 churn 是否消失？ | **是——完全消失**。U 在全部 4 个 K 下 t=6 exact stop（A：324–2876 ∝K）；D_max 与 U_density 逐轮严格下降、0 次回升；无 ∝K 尾巴 |
| U 的停止态是什么？ | D_max = 0（精确）但 **U_density ≈ 5.4×10⁻³ 非 bin 均匀**——controller 目标精确达成而 U 指标有限；与 A 的 stopping 结构同类、与 B 的机器均匀不同 |

## 2. Canonical configuration（继承，未新增）

与 M1C 完全相同：K∈{50,100,200,400}（主 K=100）、alpha=0.25、exact tol=1e-12、T_max=5000、初始状态 q_near、同一指标/诊断/eps overlay。U 通过 `redistribution_dist=DISTRIBUTIONS["uniform"]` 启用；`add_i = M_removed/K` 对所有 bin 相同，不读取任何状态信息。

## 3. 主结果

### 3.1 三变体总表（`results/variant_summary.csv`）

| metric | A | U | B |
|---|---|---|---|
| t_exact_stop | 324 / 699 / 1445 / 2876 | **6 / 6 / 6 / 6** | 3 / 3 / 3 / 3 |
| t_U_best | 292 / 653 / 1349 / 2627 | 6（= stop） | 3（= stop） |
| U_density@stop | 7.1–9.0×10⁻³ | **5.25–5.55×10⁻³** | 10⁻³¹–10⁻³³ |
| U_density@best | 3.6–3.8×10⁻³ | = stop | = stop |
| retention@stop | 0.984–0.990 | 1.000000 | 1.000000 |
| total moved | 49.0–415.3 | **1.12** | 0.32 |
| moved 中位数 | 0.147–0.154 | 0.193 | 0.104 |
| D_max 上升轮占比 | 0.36 | **0** | 0 |
| U 上升轮占比 | 0.43–0.46 | **0** | 0 |

### 3.2 U 的逐轮轨迹（K=100）

| t | a_t | D_max | moved | U_density | 累计 moved |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.50 | 0.250000 | 0.187500 | 3.333×10⁻¹ | 0.188 |
| 1 | 0.59 | 0.165025 | 0.188756 | 1.627×10⁻¹ | 0.376 |
| 2 | 0.69 | 0.097261 | 0.196815 | 7.000×10⁻² | 0.573 |
| 3 | 0.34 | 0.049016 | 0.097254 | 2.498×10⁻² | 0.670 |
| 4 | 0.84 | 0.027148 | 0.216787 | 1.428×10⁻² | 0.887 |
| 5 | 0.94 | 0.003173 | 0.235793 | 5.968×10⁻³ | 1.123 |
| 6 | stop：D_max = −2.2×10⁻¹⁶；U_density = 5.42×10⁻³ | | | | |

- 全部 6 轮是全尺寸动作；末轮是**大前缀消灭**（a=0.94，F(a)=0.94）——与 A 的小前缀消灭（F(a)≤0.5，M1B.2）几何相反。
- 停止态：全部前缀 ≤ 均匀份额（D_max 精确 ≤ 0），bin 偏差 −0.123u…+0.233u（40 高 / 60 低）——"波动但从不累计超出"的构型。
- argmax 向远端漂移（0.50→0.59→0.69）与 uniform respray 的 D-profile 线性 tilt 恒等式一致（分析文档 §6）。

### 3.3 Epsilon overlay（`results/epsilon_overlay.csv`）

- A、B 与已提交 M1C 值逐位一致。
- U：eps=0.005/0.01/0.02 → t=5（U=5.97×10⁻³，ret 0.9983）；eps=0.04 → t=4（U=1.43×10⁻²，ret 0.973）。**deadband 不再是"提前退出 churn 窗口"**（无窗口可退），只是快速收敛上的提前 1–2 轮截止。

### 3.4 U 不变量（`results/metadata.json` + 测试）

质量守恒 |Σ−1| ≤ 2.2×10⁻¹⁶；无非负违规（min bin 4.8×10⁻⁴，K=400）；逐轮加入向量 = M_removed/K 精确常数（|added − M/K| ≤ 2.5×10⁻¹⁷）；removed 全额回填（≤ 5.6×10⁻¹⁷）；q_uniform 每bin = 1/K（≤ 1e-15）且 sum=1；selection invariance 与确定性复现逐位通过。

## 4. Churn classification（U）

**Churn absent**（Prompt §9 Outcome U1，非 U2/U3）：早期改善（0.333→2.5×10⁻²，3 轮）后无任何震荡期，D_max/U_density 严格单调至 exact stop，无 ∝K 尾巴，best=stop。两个如实报告的结构性 nuance：(a) 停止态非 bin 均匀（U_density 5.4×10⁻³，处于 A 的 churn 带下沿量级）；(b) 动作不衰减——终止靠第 6 轮全尺寸消灭，不靠动作收缩。

## 5. Evidence Classification

见分析文档 §10。本报告数值为 **Verified Simulation Result**（冻结 config、一条命令复现、SHA-256 入 metadata）；D-profile tilt、一般 fixed-law net-export 恒等式（A/U 逐轮核验 ≤1.0×10⁻¹⁴）与 U 常数加入向量为 **Model-derived Mathematical Result**；t=6 跨 K 恒定、锯齿停止态、大前缀消灭几何、deadband 语义变化为 **Model-dependent Observation**；near-bias 为 churn 关键驱动、deficit-awareness 买 bin 级精度为 **Working Hypothesis（supported by this experiment）**。

## 6. Reproducibility

- Python 3.14.6（`.venv`）；numpy 2.5.3、matplotlib 3.11.2。
- 命令：`.venv/Scripts/python.exe experiments/m1c1_uniform_kernel_ablation/run_experiment.py`
- 测试：`tests/test_m1c1_uniform_kernel.py`（19 项）+ 既有九套全部通过；`python tests/count_tests.py` = 233。
- 随机性：deterministic，无种子。
- 结果文件：`baseline_integrity.csv`、`per_round_diagnostics.csv`、`variant_summary.csv`、`epsilon_overlay.csv`、两张 PNG、`metadata.json`。
- 冻结声明：config.json 运行前写入；模型仅新增可选 `redistribution_dist` 参数（默认不变）；无 post-hoc 调整（frozen kernel 集 {q_near, q_uniform, deficit-fill}）。

## 7. What This Does NOT Prove

见分析文档 §11（单 uniform control、单 q/alpha/1D/deterministic、t=6 为 4 点观察、单调性未证明、锯齿态结构未知、不证明唯一原因、无敏感性）。依据 Prompt §17，A/U/B canonical comparison 完成后停止；后续由 GPT + Owner 决策。
