# M1C.4 诊断文档：Churn-Regime Scaling & Switching Audit——既有 canonical A 轨迹的确定性动力学审计

> 状态：dynamics audit 轮（2026-09-19）。M1C.3 的直接后续；**无新模型、无新参数、无随机性**——只对既有 canonical Variant A 确定性轨迹（K∈{50,100,200,400}，q_near、alpha=0.25、exact stopping）做结构分析。
>
> 实验：[M1C4_REPORT.md](../experiments/m1c4_churn_dynamics/M1C4_REPORT.md)；前置：[M1C3_REDISTRIBUTION_LAW_ORDER.md](M1C3_REDISTRIBUTION_LAW_ORDER.md)（单步反弹判据）、[M1C2_TARGET_MATCHED_KERNEL_THEORY.md](M1C2_TARGET_MATCHED_KERNEL_THEORY.md)、[M1B_MICRO_CORRECTION.md](M1B_MICRO_CORRECTION.md)
>
> North Star：为什么"收缩占上风"与"mismatch reinjection 占上风"会在很长时间里反复切换，形成跨 resolution 的 churn regime？——**simple structure if it exists; honest absence of simple structure if it does not.**

## 0. 裁决

**核心发现（全部为确定性结构，非统计推断）**：

1. **反弹由边界位置门控（boundary-position gating）**：R_t 与扫掠边界 a_t 强单调相关——churn window 内 **a < 0.55 的轮 100% 收缩、a ≥ 0.8 的轮 100% 反弹**（K=100 与 K=400 分箱几乎一致；过渡带 [0.55, 0.8)）。这是 M1C.3 判据在 (a, R) 平面上的几何：**rebound ≡ controller 走进大前缀区域**。
2. **跨 K 稳定的运行带**：R_t 的 churn-window 分位数跨 K 几乎不变（p50 0.689–0.708、p95 3.89–4.49）；rebound fraction 0.362–0.369 就是 P(R_t > 1)——**0.36 不是新常数，而是稳定 R 分布穿阈值的结果**。
3. **简单 switching rhythm**：R 几乎从不连续（全部 5344 轮仅 9 次 R→R；R 游程中位 1、最大 2）；收缩常成对（C 游程中位 2）；RC ≈ CR（逐轮交替）；**反弹后下一轮 96–99% 发生大边界切换**。主导循环：`C → C → R → (switch) → C → C → R …`——定性上即 §11 的反馈回路。
4. **无精确周期、无低维精确回归**：无明确周期峰（lag-2 agreement 0.59–0.64 只是短 R 游程的代数结果）；标准化 4 维摘要（D_max, a, M, R）最近回归距离 0.002–0.018（非零、非接近零）。
5. **终端温和**：最后 5–20 轮无 R 的系统性剧变；D_max 中位数降至 bulk 的一半以下（0.005–0.017 vs 0.021–0.027）——**停止是普通波动探底后穿过 exact-zero 几何，而非 approach-to-stop 的独立动力签名**（与 M1B.2 的"普通波动谷 + 稀有 annihilable configuration"一致，以非随机语言表述）。

对应 Prompt §20 的结果分类：**Outcome C1（跨 K 结构稳定）+ C2（无低周期 simple rhythm）+ C4 的弱版本（终端只有温和低谷签名）**。

## 1. Churn window 定义（复用 Stage 1，非新发明）

Operational window（analysis convention，写入 config 冻结）：**onset = 首个 D_max < 0.05 的 active 轮**（M1B.2 文档已用其作为宏观清理期/churn 期的近似划分并注明"无明确突变点"），**end = t_stop − 1**（停止态本身是 stop event）。实测 onset = **23 轮（4 个 K 完全相同）**；window 长度 301/676/1422/2853。早期段（t < 23）完全被排除在 churn 统计之外。

## 2. R_t（rebound pressure ratio）定义与等价性

```text
R_t = max_{j: γ_t(j) > 0} [ M_t · max(Δ(j), 0) / γ_t(j) ]
γ_t(j) = D_max(t) − D′_TM,t(j)    （pointwise contraction margin，active 步逐点为正，M1C.2）
Δ(j) = G_near(j) − T(j) = (j/K)(1 − j/K)
```

**R_t > 1 ⟺ D_max 反弹（M1C.3 精确判据）**——ratio 形式与 M1C.3 excess 形式 `∃j: MΔ(j) > γ(j)` 等价（MΔ > γ > 0 自动要求 Δ > 0）。边界处理：γ(j) 在 active 步逐点 > 0 是 M1C.2 定理（精确算术），故不存在 γ=0 除法边界；实现仅作浮点防护。**等价性在全部 5344 个 active 轮上逐轮核验：0 mismatches**（tests/test_m1c4_churn_dynamics.py）。

## 3. 主结果 1：跨 K 稳定的运行带（Question A/B）

Churn window 分位数（`results/churn_window_summary.csv`）：

| K | window | rebound frac | R p5 | R p50 | R p95 | a p50 | M p50 | \|Δa\| p50 | switch freq |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 301 | 0.3654 | 0.432 | 0.708 | 3.890 | 0.560 | 0.1517 | 0.380 | 0.940 |
| 100 | 676 | 0.3624 | 0.421 | 0.702 | 4.391 | 0.555 | 0.1492 | 0.400 | 0.936 |
| 200 | 1422 | 0.3692 | 0.403 | 0.701 | 4.487 | 0.547 | 0.1470 | 0.390 | 0.944 |
| 400 | 2853 | 0.3687 | 0.408 | 0.689 | 4.432 | 0.545 | 0.1460 | 0.393 | 0.944 |

- **R 分布跨 K collapse**（p50 差 < 3%，p5/p95 差 < 15%）；a、M、|Δa| 中位数同样稳定（与 M1B.2 的 moved/D_max 带稳定性一致）。
- **0.36 的来源**：rebound fraction = P(R_t > 1)（判据恒等式的直接结果），而 R 分布跨 K 稳定 ⟹ 0.36–0.37 不是新的独立常数，是稳定运行带穿过阈值的频率。**仍是 Model-dependent Observation**（4 个 K、一个 q/alpha），不得当作 universal constant。

## 4. 主结果 2：边界位置门控（本轮最重要新发现）

Churn window 内按边界分箱的 P(rebound | a)：

| a 区间 | K=100 | K=400 |
|---|---:|---:|
| < 0.50（5 箱合并） | 0.000（297 轮） | 0.000（1275 轮） |
| [0.50, 0.55) | 0.000 | 0.000 |
| [0.55, 0.60) | 0.050 | 0.006 |
| [0.60, 0.70) | 0.219 | 0.371 |
| [0.70, 0.80) | 0.892 | 0.930 |
| ≥ 0.80 | 1.000 | 1.000 |

**确定性硬边界**：a < 0.55 ⟹ 必收缩；a ≥ 0.8 ⟹ 必反弹；最小反弹边界 0.57–0.59、最大收缩边界 0.76（重叠过渡带 [0.55, 0.76]）。(a, R) 操作平面（`fig_phase_a_vs_R_K100.png`）显示 R_t 随 a 单调上升、曲线自身很薄（a 几乎决定 R）。

**机制直觉（Working Hypothesis 级）**：Δ(j) = (j/K)(1−j/K) 的峰值在 j = K/2。当 controller 扫掠**大前缀**（a ≳ 0.6）时，(i) 被移除质量 M ≈ alpha·F(a) 大，(ii) 失配最强的中段 bin 落在扫掠区内且该处 margin 变薄——reinjection M·Δ 压过 γ ⟹ 反弹。小前缀扫掠时 M 小、强 Δ 区落入 suffix（那里 margin 有 M(1−T(j)) 保护项，M1C.3 结构推论）⟹ 必收缩。精确的 a\* ≈ 0.58 分界公式未推导（Next gate）。

**这就是 switching 的驱动**（Prompt §11 候选回路的定量确认）：大前缀边界 → 强 reinjection → 反弹 → **argmax 跳走**（P(switch | rebound) = 0.96–0.99）→ 新边界大概率进入小前缀收缩区 → 连续收缩对（C 游程中位 2）→ D_max 慢降、margin 消耗 → 边界再次漂回大前缀 → 下一次反弹。

## 5. 主结果 3：switching structure（Question C）

| K | R→R | R→C | C→R | C→C | R 游程 med/max | C 游程 med/max |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 3 | 107 | 106 | 84 | 1 / 2 | 2 / 4 |
| 100 | 2 | 243 | 242 | 188 | 1 / 2 | 2 / 5 |
| 200 | 2 | 523 | 522 | 374 | 1 / 2 | 2 / 5 |
| 400 | 2 | 1050 | 1049 | 751 | 1 / 2 | 2 / 6 |

（empirical transition frequencies，**非 Markov probabilities**——过程确定性。）结构：**R 几乎从不连续**（4 K 合计仅 9 次 R→R）；R→C 与 C→R 几乎相等（逐轮交替的主干）；收缩倾向成对（C 游程中位 2）。主导符号模式：`C C R (switch) C C R …`。大边界切换频率 0.94（churn window 内绝大多数轮都是真实 peak switching，与 M1A.1 一致）。

## 6. 周期性 / 回归审计（无简单低维结构的诚实记录）

- **无精确周期**：符号 lag agreement 最大 0.59–0.64 出现在 lag 2 附近——这是"R 游程 ≤ 2、C 游程 ~2"的直接代数结果，不是周期峰；R_t/D_max 自相关无孤立显著峰（`lag_agreement.csv`）。
- **无低维精确回归**：标准化 4 维摘要 (D_max, a, M, R) 的 stride 采样最近回归距离 min 0.0022–0.018、median 0.089–0.225（stride 采样使这些距离是全密度最近距离的上界）。**结论：no simple reduced switching law / low-period cycle was identified**——churn 是 bounded deterministic switching regime，不是低周期循环（Prompt §15 遵守；不使用 attractor/chaos/ergodicity 语言）。

## 7. 终端行为（Distance-to-stop + cross-K alignment）

| K | last-5 D_max med | last-20 D_max med | bulk D_max med | last-5 R med | last-5 M med | last-5 a range |
|---:|---:|---:|---:|---:|---:|---|
| 50 | 0.0051 | 0.0208 | 0.0275 | 0.524 | 0.117 | [0.12, 0.82] |
| 100 | 0.0055 | 0.0146 | 0.0246 | 0.617 | 0.139 | [0.07, 0.86] |
| 200 | 0.0093 | 0.0164 | 0.0224 | 0.698 | 0.150 | [0.195, 0.85] |
| 400 | 0.0102 | 0.0159 | 0.0210 | 0.544 | 0.100 | [0.035, 0.7375] |

- **温和签名**：终端 5–20 轮 D_max 中位系统性低于 bulk（波动探底），R 中位不高于 bulk、rebound 频率无显著变化（0.2–0.4，小样本噪声），M 略低，a 范围收窄并触及小区间（K=400 的 0.035 正是 annihilator 边界）。
- **判读**：stop 之前**没有**独立的 deterministic approach-to-stop 动力学（无 R 骤降、无 M 收缩到零的过程）——exact stop 是普通 churn 波动**偶然**（非随机意义：确定性轨迹低频地）把状态送入 annihilable configuration。与 M1B.2 的机制结论一致并加以动力学补充：等待发生在波动带内，而非某个渐近出口。
- Cross-K alignment（τ = t_stop − t，最后 40 轮，`fig_terminal_alignment.png`）：四个 K 的 D_max 都在 10⁻³–10⁻² 低谷段震荡、R 围绕 1 以下、M 保持全尺寸——终端 pattern 跨 K 相似但**温和**（无尖锐共同序列）。

## 8. Evidence Classification

- **Model-derived Mathematical Result**：仅 M1C.3 已证明的 `R_t > 1 ⟺ rebound`（本轮逐轮重新核验，0 mismatches）及 R_t 定义的良定性（γ > 0）。本轮**无新定理**，不制造。
- **Verified Numerical Analysis**：5344 轮 R_t 等价性核验；churn-window 分位数表；分箱 P(rebound | a)；run/transition 表；lag agreement 与自相关；回归审计；终端统计与 τ 对齐。全部确定性复现（两次运行输出位级一致）。
- **Model-dependent Observation**：R/a/M/|Δa| 的跨 K 稳定带；0.36–0.37 = P(R>1)；a-门控结构（a<0.55 必收缩、a≥0.8 必反弹）；`C C R` 主导模式与 RR 罕见性；终端温和低谷。
- **Working Hypothesis**：churn = "边界漂移进入大前缀区 → mismatch reinjection 胜过 margin → 反弹 → 边界切换 → 收缩对 → margin 再消耗"的确定性切换循环的长期行为；a\* ≈ 0.58 分界的解析式。均未证明。

## 9. What churn now means（当前最严格版本）

**已解释（theorem/精确层级）**：单轮反弹的充要条件（M1C.3）；rebound fraction 恒等于 P(R_t > 1)（判据恒等式）；本轮新增的确定性事实——反弹被边界位置门控、R 游程 ≤ 2、反弹后必切换、终端无独立签名。

**已解释（observation 层级）**：0.36 的跨 K 稳定性（R 分布稳定）；switching 的主导节奏（C-C-R）；stop 与普通波动谷的关系。

**未解释/未证明**：为什么 R 分布稳定在该形状（a\* 分界的解析来源）；为什么边界会反复漂回大前缀区（margin 消耗动力学无定量定理）；长期不收敛性（churn 永不终止到带外）——**single-step rebound theorem ≠ long-lived churn dynamics**，后者仍是无证明的 descriptive regime + Working Hypothesis。

## 10. Limitations

1. 单一 q/alpha/一维/deterministic path/4 个 K；onset=0.05 是沿用 M1B.2 的近似 divider（无突变点声明）；
2. R_t 压缩了全部逐前缀信息为单一 max——门控结构的精确定量（为什么是 0.55/0.8）未推导；
3. 回归审计是 stride 采样上界、4 维摘要；lag 分析 ≤ 100；
4. 终端窗口小样本（5/10/20）统计噪声大，只作定性比较；
5. 无 chaos/ergodicity/attractor 分析，也不声称。

## 11. Next Research Gate（交回 GPT + Owner）

1. **a\* 分界的解析式**：从 γ(j) 与 M·Δ(j) 的显式公式推导反弹/收缩的边界阈值（本轮只有 observation）；
2. **margin 消耗动力学**：为什么边界反复漂回大前缀区（CC 对与 C→R 转换的机制）；
3. **长期化**：本轮结构是否足以支撑"A 不收敛/带不变性"定理；
4. 依据 Prompt §22，四项交付（window audit、R 跨 K、switching、terminal alignment）完成后停止。

## 12. Repository 资产

- 分析：`experiments/m1c4_churn_dynamics/`（config / run_analysis.py / results 12 文件 / M1C4_REPORT.md）；一条命令复现；
- 绘图：`plotting.py` 新增 4 函数（R threshold、cross-K quantiles、(a,R) 平面、terminal alignment）；
- 测试：`tests/test_m1c4_churn_dynamics.py`（9 项）；无新模拟机制、模型代码零修改（plotting 除外）。
