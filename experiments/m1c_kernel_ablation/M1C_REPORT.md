# M1C 实验报告：Correction-Kernel Ablation（A vs B）

> 日期：2026-09-19
>
> 分析文档：[docs/M1C_CORRECTION_KERNEL_ABLATION.md](../../docs/M1C_CORRECTION_KERNEL_ABLATION.md)（研究问题、Variant 定义、invariants、数学结果、判读与证据分级）；前置：[M1B_MICRO_CORRECTION.md](../../docs/M1B_MICRO_CORRECTION.md)
>
> 代码：`src/sand_m0/adaptive.py`（新增 `deficit_fill_redistribution` kernel 与 `redistribution` 参数；默认值保持 canonical 行为位级不变）；配置：[config.json](config.json)（运行前冻结）
>
> 性质：单变量机制对照实验。Variant A = canonical M1（fixed biased `q_near` redistribution，未动）；Variant B = selection rule 逐位不变，仅把 redistribution 换成 target-aware deficit filling（`u_i = 1/K`）。

## 1. 裁决（ headlines ）

| 问题 | 答案 |
|---|---|
| Variant A 是否复现 M1B.2 Stage 1 canonical churn behavior？ | **是**：112/112 逐位对照通过（见 §3） |
| Variant B 是否仍出现 churn？ | **否**：churn 完全消失。B 在全部 4 个 canonical K 下 **3 轮**达到机器精度均匀并 exact stop（A：324–2876 轮） |
| 研究问题答案 | **Churn does NOT persist.** 替换 redistribution kernel 后 persistent churn 不复存在——支持（supported by this experiment）当前 Working Hypothesis |

## 2. Canonical configuration（与 Stage 1 完全一致）

`q_near`、alpha=0.25、exact tol=1e-12、T_max=5000、K∈{50,100,200,400}（主 K=100）、deterministic exact reference path、epsilon overlay {0.005, 0.01, 0.02, 0.04}。Variant B 目标 `u_i = 1/K`。无新增参数、无 sweep、无 post-hoc 调整。

## 3. Baseline reproduction（Stage 1 门禁）

对每个 K 重跑 canonical dynamics，与已提交 M1B.2/M1A.1 结果做 112 项对照（`results/baseline_reproduction.csv`）：exact stop 轮与 pre-stop 结构、late-stage（D_max<0.005）moved/net_export/L1 中位数、D_max/U_density 波动带分位数、dip 频次与首达、16 个 (K,eps) first-passage t 与 retention——**全部通过，舍入字段逐位一致**。若失败则按 config.baseline_policy 中止（未触发）。

## 4. 主结果

### 4.1 轨迹量总表（`results/variant_summary.csv`）

| variant | K | t_stop | t_U_best | U_density@stop | retention | moved 中位数 | moved 总量 | D_max 上升轮占比 | U 上升轮占比 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 50 | 324 | 292 | 9.00×10⁻³ | 0.9840 | 0.1544 | 49.02 | 0.362 | 0.432 |
| A | 100 | 699 | 653 | 7.07×10⁻³ | 0.9901 | 0.1502 | 103.30 | 0.360 | 0.431 |
| A | 200 | 1445 | 1349 | 7.88×10⁻³ | 0.9873 | 0.1480 | 210.59 | 0.368 | 0.451 |
| A | 400 | 2876 | 2627 | 7.14×10⁻³ | 0.9892 | 0.1466 | 415.27 | 0.368 | 0.463 |
| B | 50 | **3** | **3** | **3.6×10⁻³³** | **1.000000** | 0.1058 | **0.3251** | **0** | **0** |
| B | 100 | **3** | **3** | **1.8×10⁻³¹** | **1.000000** | 0.1033 | **0.3201** | **0** | **0** |
| B | 200 | **3** | **3** | **3.3×10⁻³³** | **1.000000** | 0.1046 | **0.3213** | **0** | **0** |
| B | 400 | **3** | **3** | **1.6×10⁻³³** | **1.000000** | 0.1040 | **0.3207** | **0** | **0** |

A 侧保留全部 churn 指纹（moved ~0.15/轮、D_max 带 0.003–0.08 震荡、best 早于 stop 32–249 轮、∝K 尾巴）；B 侧 D_max 与 U_density 逐轮严格下降、无任何回升轮、best = stop。

### 4.2 Variant B 的 3 轮收敛路径（K=100；4 个 K 只差格点尺度）

| t | a_t | D_max | moved | U_density | D_def | E_minus | fill ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.50 | 0.250000 | 0.187500 | 3.33×10⁻¹ | 0.270825 | 0.083325 | 0.6923 |
| 1 | 0.33 | 0.083325 | 0.103331 | 4.39×10⁻² | 0.110275 | 0.006944 | 0.9370 |
| 2 | 0.11 | 0.006944 | 0.029236 | 6.59×10⁻⁴ | 0.029236 | 0（精确） | 1.0000 |
| 3 | stop：D_max = 8.9×10⁻¹⁶ ≤ 1e-12；U_density = 1.84×10⁻³¹；max\|p−u\| = 1.6×10⁻¹⁷ | | | | | | |

机制：高于均匀的 bin 数按 ⌈K/3⌉ → ⌈K/9⌉ → 0 收缩；state 2 时最高于均匀 bin 的质量（0.0112 < 4u/3 = 0.01333）已被最后一扫压到均匀以下 → `E_minus = 0` → 恒等式 `D_def = E_minus + M_removed` 退化为饱和填充 → 精确均匀。推导见分析文档 §5。

### 4.3 Epsilon overlay（`results/epsilon_overlay.csv`）

- A：与已提交 M1B 值逐位一致（eps=0.005：t=253/340/338/321；eps=0.04：t≈34）。
- B：eps=0.005 → t=3（retention 1.0）；eps=0.01/0.02/0.04 → t=2（U_density=6.59×10⁻⁴，retention 0.9980）。跨 K 恒定；停止态均匀度已低于 A 整个活跃期波动带下沿（p5≈5.6×10⁻³）一个量级。

### 4.4 Invariants（B，全部 K 全部活动轮）

| 检查 | 实测最坏值 |
|---|---|
| 质量守恒 \|Σp−1\| | ≤ 4.4×10⁻¹⁶ |
| 非负性 | 最小 bin 质量 1.7×10⁻³（K=400） |
| fill 引起的过填 | 无（`add_i ≤ deficit_i` 由恒等式保证；测试验证 `p_next ≤ max(p_minus, u)`） |
| deficit 恒等式残差 \|D_def − E_minus − M_removed\| | ≤ 1.8×10⁻¹⁶ |
| fill ratio = M_removed/D_def | ≤ 1.000000（=1 恰在 t=2 饱和轮） |
| degenerate guards（D_def=0 / M>D_def raise） | 从未触发（数学上不可能，见分析文档 §5） |

## 5. Unexpected / Notable

1. **churn 不是被抑制而是被消灭**：B 没有更慢的收敛或残余微震荡——3 轮精确均匀，机器精度验证。
2. **exact stopping 跨 K 恒定（t=3）**：M1A/M1B 的 ∝K 等待尾巴在 B 下整个消失。B 主动构造 A 中需要等待 ∝K 轮的 annihilable 构型（的饱和版）。
3. **fill ratio 单调升至 1**（0.692 → 0.937 → 1.000）：q 的"再造过剩"效应被消除后，正过剩单调集中于缩短的近端前缀。
4. A 的 exact-stop 态 U_density（7–9×10⁻³）与 B 的（10⁻³¹）相差 ~28 个数量级——A 的"annihilation"清零的是 cumulative excess，不是逐 bin 偏差；B 的饱和填充把两者同时清零。

## 6. Evidence Classification

见分析文档 §10。本报告的全部数值为 **Verified Simulation Result**（固定代码版本、冻结 config、一条命令复现、输出 SHA-256 入 metadata）；§5 恒等式及推论为 **Model-derived Mathematical Result**；跨 K 恒定 t=3、收缩模式等为 **Model-dependent Observation**；机制解释的推广性为 **Working Hypothesis**。

## 7. Reproducibility

- Python 3.14.6（仓库 `.venv`）；numpy 2.5.3、matplotlib 3.11.2。
- 命令：`.venv/Scripts/python.exe experiments/m1c_kernel_ablation/run_experiment.py`
- 测试：`tests/test_m1c_kernel_ablation.py`（20 项）+ 既有八套（37+25+54+18+18+20+11+11）全部通过；`python tests/count_tests.py` = 214。
- 随机性：deterministic，无种子。
- 结果文件（`results/`，SHA-256 见 metadata.json）：`baseline_reproduction.csv`、`per_round_diagnostics.csv`（A+B 全部活动轮逐轮诊断）、`variant_summary.csv`、`epsilon_overlay.csv`、两张 PNG、`metadata.json`。
- 冻结声明：config.json（含 variant 定义、invariants、K/T/alpha/tol/eps、baseline 门禁策略）运行前写入；模型除新增 kernel 参数外未改；无 post-hoc 调整。

## 8. What This Does NOT Prove

见分析文档 §11（oracle 信息、单 q/单 alpha/1D/deterministic、t=3 为 4 点观察、单调性未证明、不证明唯一原因、未做敏感性）。依据 Prompt 冻结的停止条件，本实验在完成 A vs B 后停止；后续方向（部分感知 kernel、跨 q/alpha gate、理论化）由 GPT + Owner 决策。
