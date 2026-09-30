# M1A.1 实验报告：Long-horizon & Resolution Diagnostics

> 日期：2026-09-18
>
> 理论/分析文档：[docs/M1A_LONG_HORIZON.md](../../docs/M1A_LONG_HORIZON.md)；模型定义：[docs/M1A_MODEL.md](../../docs/M1A_MODEL.md)；代码：`src/sand_m0/adaptive.py`（未修改）+ `src/sand_m0/diagnostics.py`（新增，仅分析）；配置：[config.json](config.json)（运行前冻结）
>
> 证据等级：见第 8 节。核心修正：**M1A 的 exact stopping criterion 在全部四个分辨率下都被触发**——T=100 时"不停止"是预算表象。

## 1. 目的与约束

回答唯一问题：M1A 在 T=100 下的"持续改善、边界振荡、不停止、最佳贴预算末端"是长期反馈动力学，还是有限预算/分辨率/argmax 离散的表象？

约束：M1A feedback rule、stopping rule（tol=10⁻¹²）、throw distribution、更新规则**一概未改**；无感知阈值、平滑、迟滞、粗网格、替代 q、噪声修正、2D。K 的 bin 概率由连续 CDF 精确积分生成（`bin_probs`），从不跨 K 复制。

## 2. 配置

| 项 | 值 |
|---|---|
| q / alpha | `q_near = 2(1−x)` / 0.25 |
| T_max | 5000（每 K） |
| K | 50, 100, 200, 400（基线 100） |
| stop_tolerance | 1e-12（冻结） |
| 诊断阈值 | A_1% = 0.99·D_max，A_5% = 0.95·D_max，大跳阈值 |Δa|>0.2（仅诊断，不改 policy） |
| recurrence | lag 1..50，仅活动段，L1 |
| 复现命令 | `.venv/Scripts/python.exe experiments/m1a_1_long_horizon/run_experiment.py` |

测试：`tests/test_m1a_1_diagnostics.py`（21 项）通过；四套既有测试（37+25+53+18 项）保持通过。

## 3. 主结果

### 3.1 精确停止在所有 K 下触发（对 T=100 结论的正式修正）

| K | 停止轮 | 比值 | 停止态 U_density | 停止态 D_max | 近半区终值 |
|---:|---:|---:|---:|---:|---:|
| 50 | 324 | — | 9.00×10⁻³ | 6.7×10⁻¹⁶ | 0.4790 |
| 100 | 699 | ×2.16 | 7.07×10⁻³ | 1.6×10⁻¹⁴ | 0.4872 |
| 200 | 1445 | ×2.07 | 7.88×10⁻³ | −3.3×10⁻¹⁶ | 0.4929 |
| 400 | 2876 | ×1.99 | 7.14×10⁻³ | −1.0×10⁻¹⁵ | 0.4772 |

停止后状态逐位冻结（最后 100/500 轮窗口 std = 0）。停止态满足 `F(x) ≤ x ∀x`（每前缀不过剩），近半区质量略低于理想 0.5。

**停止时间近似 ∝ K**（比值趋近 2）。**停止态 U_density 跨 K 稳定**（0.0071–0.0090）。这是 Model-dependent Observation——尤其 K→∞ 的外推（连续极限是否有限时间停止）未定。

### 3.2 U 的长期行为（K=100）

| t | U | U_density | D_max | a_t | 近半区 |
|---:|---:|---:|---:|---:|---:|
| 0 | 3.333×10⁻³ | 0.3333 | 0.2500 | 0.50 | 0.7500 |
| 100 | 1.704×10⁻⁴ | 0.01704 | 3.092×10⁻² | 0.73 | 0.5153 |
| 500 | 8.595×10⁻⁵ | 0.008595 | 2.297×10⁻² | 0.19 | 0.5112 |
| 653（最佳） | **3.821×10⁻⁵** | **0.003821** | 1.119×10⁻² | 0.58 | 0.5058 |
| 699（停止）→5000 | 7.067×10⁻⁵ 冻结 | 0.007067 | 1.64×10⁻¹⁴ | 冻结 | 0.4872 |

- 总体单调改善（t=100 时已完成 95% 的降幅），t=653 达 U 最优，随后**回升 85%** 至停止；
- **模型自身判据（D_max→0）与 U 指标轻微分歧**：停止态比 U 最优态差 1.85 倍，且近半区低于 0.5（判据允许远端过满）；
- 无长期波动、无周期、无持续恶化；
- **跨 K 的 U_density 最优值稳定**：3.71 / 3.82 / 3.69 / 3.56 ×10⁻³（K=50/100/200/400）。

### 3.3 Recurrence（活动段，lag 1..50，L1）

单步典型 L1 = 7.63×10⁻²；全部 50 个 lag 的最小值 = 2.74×10⁻²（lag 1）——**无短周期近重复状态**（最近的"重访"仍差典型单步的 36%）。不规则切换，无 limit cycle 证据；不声称 chaos。

### 3.4 边界歧义重新刻画（699 个活动决策）

| 诊断 | 数值 | 占比 |
|---|---:|---:|
| top-2 值差 ≤ 1%·D_max | 601 | 86% |
| top-2 为相邻格点 | 671 | 96% |
| 远处第二峰（sep>0.1） | 8 | 1.1% |
| A_1% 多分量 | 75 | 10.7% |
| farthest A_1% 点中位数 | 0.02 | — |

**近平局是局部格点级现象**（runner-up 几乎总是相邻 bin）；"多个远离区域同样合理"基本不存在。MC 边界高频翻转的解释由此收窄为相邻格点间的噪声选择——未来感知尺度研究只需把相邻并列合并（更粗的决策分辨率），而非远处消歧。

### 3.5 边界大跳：全部是真实峰切换

活动段 649/699（93%）的决策 |Δa|>0.2；**649 个大跳全部 classified 为 peak switch**（新旧边界位置并非都在双方 A_5% 集合内），0 个 flat-top drift。D(x) 峰每轮被扫掠真实重排——大跳是反馈动力学的组成而非平顶漂移。

## 4. 精确停止的措辞与可达性

- 正确措辞（按治理）：**"在当前 deterministic M1A、该初始条件和 T≤5000 内，exact stopping criterion 于 t=324/699/1445/2876（K=50/100/200/400）触发。"**
- 不得声称："M1A 永远会停止"（无 K→∞ 证明）或"收敛"（停止≠U 最优）。
- 可达性分析尝试：停止要求 `m_0 ≤ 1/K`；扫掠回填使 `m_0' = 0.75·m_0 + alpha·M_prefix·q_0`，q_near 的 `q_0·K≈1.99`——扫后仍 ≤ 1/K 需要该轮 `M_prefix ≲ 0.5025`，而边界大部分时间 >0.5。这解释了停止缓慢，但不构成不可达证明（完整证明未获得，如实止步）。详见 [M1A_LONG_HORIZON.md](../../docs/M1A_LONG_HORIZON.md) §8。

## 5. One-sided CDF discrepancy 术语注记

`D_max(t) = max_j [F_t(j/K) − j/K]` 与针对 Uniform(0,1) 的 **one-sided Kolmogorov–Smirnov functional** `sup_x [F(x) − G(x)]`（G(x)=x）是同一数学对象；阶梯 CDF 下离散 max **精确等于**连续 sup（已证明：sup 在 bin 右端点取到；数值验证于测试）。这是 Literature-supported terminology 级别的数学联系：**不表示 M1A 在做 KS 检验**，无假设检验成分，不构成 novelty claim。

## 6. 对 M1A 解读的修正

完整修正表见 [M1A_LONG_HORIZON.md](../../docs/M1A_LONG_HORIZON.md) §9。要点：不停止→停止（t=699）；持续振荡→停止前瞬态+冻结；最佳贴末端→U 最优 t=653、停止态 1.85×最优；"高度模糊"→局部格点级近平局；大跳→真实峰切换。M1A_REPORT.md 保留为 T=100 历史记录，本报告为其正式修正与延伸。

## 7. 文件

`experiments/m1a_1_long_horizon/results/`（SHA-256 见 metadata.json）：

- `long_horizon_trajectory_K100.csv` — 5001 行逐轮：U、U_density、D_max、连续 sup 校验列、a_t、near、total mass、active、margin、全部歧义诊断
- `checkpoints_and_windows.csv` — 检查点 + 末 100/500 窗口统计
- `k_sensitivity_summary.csv` — 四 K 汇总
- `recurrence_l1_K100.csv` — lag 1..50
- `boundary_ambiguity_summary.csv`、`largest_boundary_jumps.csv`
- `fig_long_horizon_K100.png`、`fig_k_sensitivity.png`、`fig_d_profiles.png`
- `metadata.json`（含 headline、git 状态、输出哈希）

## 8. Evidence Classification

- **Model-derived Mathematical Result**：`max_j D(j) = sup_x [F_t(x) − x]`（阶梯 CDF 精确性）——证明 + 数值核验；停止态刻画 `D_max=0 ⟺ F(x)≤x ∀x`（承自 M1A）。
- **Verified Simulation Result**：四 K 的停止轮与冻结态、U/D_max 检查点与窗口统计、recurrence 表、歧义统计、大跳分类、K=100 与 M1A canonical 前 101 行逐位一致。
- **Model-dependent Observation**：T_stop ∝ K 标度及其 K→∞ 含义未定；停止态 ≠ U 最优态（1.85×）；近平局的局部性；U_density 最优值跨 K 稳定。
- **Literature-supported terminology**：one-sided KS functional 术语联系（标准定义；非检验、非 novelty）。
- **Working Hypothesis**：连续极限停止时间可能发散；停止判据与 U 指标的分歧是"感知容差"（M1B）的模型内证据。
- 无新的 Real-world Observation；adaptive rule 保持 Model Assumption。

## 9. What This Does NOT Prove

见 [M1A_LONG_HORIZON.md](../../docs/M1A_LONG_HORIZON.md) §10：不证明连续极限有限时间停止、不证明可达性条件、不证明周期性缺失的一般性（仅 lag≤50 L1）、不做 chaos/收敛率声明、不改变 M1A 的 Model Assumption 地位。

## 10. Recommended Next

**M1B — Perceptual Tolerance / Deadband Feedback**（本轮诊断的直接后继，待 Owner + GPT 批准）：诊断显示 (a) exact stop 可达但缓慢（T∝K），(b) 停止态比 U 最优态差 1.85 倍（判据与指标分歧），(c) 边界近平局是格点级平坦。三者共同指向"以感知容差替代精确零过剩"的正式建模。其次：T_stop ∝ K 的解析解释。
