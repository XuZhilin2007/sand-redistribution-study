# M0.1 实验报告：correction strength 是 M0 的时间尺度参数吗

> 日期：2026-09-18
>
> 模型定义：[docs/M0_MODEL.md](../../docs/M0_MODEL.md)；机制推导：[docs/M0_MECHANISM.md](../../docs/M0_MECHANISM.md)；代码：`src/sand_m0/`（新增 `mechanism.py`）；配置：[config.json](config.json)
>
> 证据等级声明：本文数值结论均可由一条命令精确复现，属 **已验证模拟结果**（仅对 M0 成立）；由模型定义直接推导、并经独立数值核验的解析结论见第 5 节的等级讨论。本文不产生真实沙坑结论，alpha 的一切数值都不是现实回收效率。

## 1. 实验问题

只回答两个问题：

1. 为什么 EXP-M0-A 在第 4–5 轮达到最佳均匀状态？
2. 改变 `alpha`，主要改变"最佳状态本身"，还是"到达最佳状态需要多少轮"？

本实验**不是参数优化**：`alpha ∈ {0.10, 0.25, 0.50, 0.75, 1.00}` 是运行前固化的声明集合，不从中挑选"最好 alpha"。其余一切（`q=q_near`、`a=0.5`、`K=100`、`T=30`、指标、更新规则）与 M0 canonical baseline 完全一致。

## 2. 方法

- **Deterministic**：六个 alpha（含边界 0.0）各跑完整轨迹（`run_deterministic`）。
- **解析机制**：`src/sand_m0/mechanism.py` 从 M0 定义精确导出闭式解（近区 `z·q_i`、远区 `(4−3z)·q_i`，`z=(1−alpha/4)^t`）、near-mass recurrence、精确二次函数 `U(z)` 与最优点 `z*`。全部为离散 bin 上的精确代数，无数值近似。
- **MC spot check**：alpha ∈ {0.25, 0.75, 1.00}，各 5 种子（20260917–21，N=10⁵），确认粒子版本无系统偏离；另做一次性 40 种子收敛探针（不入库，记录于机制文档 §2.7）。
- **测试**：`tests/test_m0_mechanism.py`（25 项断言）在实验运行前全部通过；原 M0 测试（37 项）保持通过。

复现命令：

```bash
.venv/Scripts/python.exe experiments/m0_1_mechanism/run_experiment.py
```

## 3. 结果

### 3.1 alpha 对比（deterministic）

| alpha | 每轮因子 ρ=1−alpha/4 | 理论 t* | 实际最佳轮 | 最佳 U | gap = U(最佳) − U(z*) | 最佳时近区质量 | max\|U_sim − U(z) 二次式\| |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0（边界） | 1.000000 | — | 0（状态不动） | 3.333000e-03 | 1.667e-03 | 0.7500 | 8.2e-18 |
| 0.10 | 0.975000 | 11.3643 | 11 | 1.667350e-03 | 1.29e-06 | 0.5677 | 9.1e-18 |
| 0.25 | 0.937500 | 4.4581 | 4 | 1.679567e-03 | 1.35e-05 | 0.5794 | 1.2e-17 |
| 0.50 | 0.875000 | 2.1547 | 2 | 1.672596e-03 | 6.53e-06 | 0.5742 | 2.8e-17 |
| 0.75 | 0.812500 | 1.3857 | 1 | 1.770316e-03 | **1.04e-04** | 0.6094 | 3.1e-17 |
| 1.00（边界） | 0.750000 | 1.0001 | 1 | 1.666063e-03 | **2.11e-11** | 0.5625 | 2.1e-17 |

读法：

1. **理论 t* 与实际最佳轮逐项吻合**（11.36→11、4.46→4、2.15→2、1.39→1、1.00→1）。
2. **六个 alpha 的 U(z) 是同一条曲线**：每个 alpha 的模拟 U 与同一二次式的偏差都在浮点级（≤8.2×10⁻¹⁸）；理论最佳轮由 `t* = ln(z*)/ln ρ` 给出，与模拟 argmin 一致。
3. **最佳状态本身不随 alpha 改变**：连续意义的最优状态由 `z*` 唯一决定（近区质量 0.75·z* ≈ 0.5625），六个 alpha 的"实际能停下的最好状态"仅因轮次离散化而不同——gap 从 alpha=1.00 的 2×10⁻¹¹（第一拍精确命中 z=0.75）到 alpha=0.75 的 1.04×10⁻⁴（z 网格跳过 z* 最远）。
4. **alpha=0 边界**：状态永不移动（U(t)=U(0) 恒定，测试验证），不存在时间推进，也无从谈最佳轮。
5. **alpha=1 边界**：每轮近区全部重泼，z 每轮 ×0.75，第 1 轮即达连续最优；过程仍然健康（重泼的 75% 落回近区，系统不会一步清空近区）。

### 3.2 时间重合性（time-rescaling）检查

把横轴从轮次 t 换成状态变量 `z = ρ^t`（等价地，近区质量 = 0.75·z）：

- 六个 alpha 的 (z, U) 点**精确落在同一条二次曲线上**（合并残差 ≤ 10⁻¹⁷，见 `fig_rescaling_u_vs_z.png` 右幅）；
- 理论曲线与连续极限（z*=3/4 精确）几乎重合，离散修正 −2.8×10⁻⁵。

**结论：在 M0 该结构下，这是明确的 time-rescaling property——alpha 只重参数化时间，不改变状态路径。** 该性质是 M0 特定结构（固定 q 重采样 + 固定近区 + 线性守恒系统）的推论，不得推广到其他模型或真实过程。

### 3.3 MC spot check

| alpha | 单种子 max\|U−U(z)\| | 5 种子均值 max\|mean−(U(z)+底线)\| |
|---:|---:|---:|
| 0.25 | 1.42e-04 | 6.50e-05 |
| 0.75 | 1.32e-04 | 6.00e-05 |
| 1.00 | 1.36e-04 | 3.65e-05 |

单种子偏差是多批次抽样噪声的尾部（理论 sd ≈ 2–4×10⁻⁵）；三个 alpha 的均值偏差方向相关是因为共用同一批种子。40 种子收敛探针（alpha=0.25）：各轮 `mean(U_MC − U(z))` 全部落在理论底线 9.87×10⁻⁶ 的 ±1σ 内。**粒子版本无系统偏离。**

## 4. 对两个研究问题的最终回答

1. **为什么第 4–5 轮最佳**：M0-A 的 deterministic 状态被单一变量 z 完全描述；U 是 z 的精确二次函数，最低点 z*=3/4（近区质量 56.25%）；baseline alpha=0.25 下 z=(15/16)^t，z=3/4 对应 t≈4.46；整数轮 4、5 从两侧夹住它且 U(4)、U(5) 仅差 0.25%，形成浅平台，第 4 轮因 z₄ 略近 z* 而胜出。完整推导见 [M0_MECHANISM.md](../../docs/M0_MECHANISM.md)。
2. **alpha 改变什么**：alpha 只通过 ρ=1−alpha/4 决定 z 的每轮衰减速度——**主要改变"多少轮到达"，不改变"最佳状态本身"**。精确表述：状态路径 {m(z)} 与连续最优 (z*, U(z*)) 对一切 alpha 相同；alpha 影响的只有 (a) 到达时间 t*，(b) 因轮次离散化而产生的可达最优 gap（≤10⁻⁴ 量级）。

## 5. 证据等级分类

**A. M0 mathematical structure（由定义直接推出的精确性质）**：闭式解 m(z)、near-mass recurrence、U(z) 精确二次式、z* 公式、连续 z*=3/4、time-rescaling。这些是 M0 定义的代数推论，本应由"数学推导"等级承载。

**B. Verified Simulation Result（本实验的数值验证）**：上述结构在六个 alpha 下与逐步模拟器一致至浮点级（≤8.2×10⁻¹⁸）；理论 t* 与实际最佳轮逐项吻合；离散 gap 数值表；MC spot check 无系统偏离。全部可一条命令复现（commit、config、seeds 见 metadata.json）。

**C. Working Hypothesis（待其他设定检验）**：
1. time-rescaling 可能对"固定 throw distribution + 固定回收几何"类线性过程普遍成立，而对状态依赖的反馈策略（Phase 3 的 adaptive/density-based policy）不一定成立——那是它最可能失效的地方。
2. z* 由形状集中度比 S_n:S_f 决定（q_near 下 7:1 → 3/4）；不同 throw distribution 会给出不同 z* 与不同最佳轮——Phase 2 的分布族比较可直接检验该预测。

**治理建议（单独提出，待 Owner + GPT 决定，本文不自行改动治理体系）**：本轮首次出现"由 canonical model 定义精确推导并经独立数值核验"的结果。它们目前只能挂在 Verified Simulation Result 之下，但二者性质不同：后者是"模拟观察到的现象"，前者是"模型内在的数学事实"。建议正式增设标签：

> **Model-derived Mathematical Result（模型推导数学结果）**：由已版本化的 canonical model 定义出发、以精确推导（非近似模拟）得到，并已对照同一版本的模拟器做过独立数值核验的结论。它仍然只对该模型成立，不因"数学上精确"而获得任何现实有效性。

在 Owner 批准前，本文将上述结果保守标注为 Verified Simulation Result，不启用新标签。

## 6. 本文不证明的事情

- 不证明真实沙坑中"回收效率只改变速度"——M0 的一维等比例回收与现实动作存在已知 model mismatch（见假设台账 §5）。
- 不证明所有 redistribution model 都有 time-rescaling 性质——该性质依赖本轮明确的线性结构。
- 不证明 finite-time optimum 普适存在——M0-A 的最优点存在性是特定 q、特定 a、特定指标的代数事实。
- 不声明学术新颖性；无文献支持结果。
- alpha 的任何数值（含 alpha=1 的"一步最优"）都不具有现实操作含义。

## 7. 复现信息

- **Python** 3.14.6；**依赖** numpy 2.5.3、matplotlib 3.11.2（仓库专用 `.venv`）。
- **命令**：`.venv/Scripts/python.exe experiments/m0_1_mechanism/run_experiment.py`
- **测试**：`.venv/Scripts/python.exe tests/test_m0_mechanism.py`（25 项）与 `tests/test_m0.py`（37 项）。
- **种子**：MC spot check 用种子 20260917–20260921（M0 种子表前 5 个，显式记录于 config 与 metadata）；deterministic 部分无随机性。
- **结果文件**（`experiments/m0_1_mechanism/results/`，SHA-256 见 metadata.json）：
  - `alpha_comparison.csv` — 每 alpha 的 ρ、t*、实际最佳轮、最佳 U、gap、理论偏差
  - `deterministic_trajectories.csv` — 6×31 行逐轮 (z, near_mass, U, TV, 理论值)
  - `theory_u_of_z.csv` — 共同理论曲线（离散精确 + 连续极限）
  - `mc_spotcheck.csv` — 3 alpha × 5 seeds × 31 轮逐点记录
  - `fig_alpha_uniformity_vs_round.png`、`fig_rescaling_u_vs_z.png`、`fig_mc_spotcheck.png`
  - `metadata.json` — 参数、理论摘要、git 状态、输出哈希
- **对既有 M0 输出的影响**：本轮为机制分析在 `simulate.py` 中新增 `run_deterministic_states`（重构 `run_deterministic` 为其包装）。重构后重跑 M0 官方实验，16/16 个数据输出**位级一致**（以 SHA-256 验证），canonical 结果未受影响。

## 8. 局限与下一步

1. 单一分布（q_near）下的机制；z* 对分布族的依赖（S_n:S_f 比值如何随 q 变化、最佳轮如何移动）是 Phase 2 的直接可检验预测。
2. 线性、时齐、固定几何是 time-rescaling 的前提；引入状态依赖回收（Phase 3）后该性质预计失效，失效方式本身就是有价值的观察。
3. 本轮未触碰 stopping rule；M0-A 的 z-路径为将来评估"在线停止能否停在 z* 附近"提供解析基准。
