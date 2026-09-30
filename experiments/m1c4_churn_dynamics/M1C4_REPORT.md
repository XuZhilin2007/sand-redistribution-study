# M1C.4 实验报告：Churn-Regime Scaling & Switching Audit

> 日期：2026-09-19；分析文档：[docs/M1C4_CHURN_REGIME_AUDIT.md](../../docs/M1C4_CHURN_REGIME_AUDIT.md)；前置：[M1C3_REDISTRIBUTION_LAW_ORDER.md](../../docs/M1C3_REDISTRIBUTION_LAW_ORDER.md)
>
> 性质：对既有 canonical Variant A 确定性轨迹的 dynamics audit。无新模型/参数/随机性；模型代码零修改。

## 1. 裁决

| 问题 | 答案 |
|---|---|
| R_t > 1 ⟺ rebound 是否逐轮成立？ | 是——全部 5344 active 轮 0 mismatches |
| churn diagnostics 是否跨 K 稳定？ | 是——R 分位数 p50 0.689–0.708、rebound fraction 0.362–0.369 |
| 0.36 是什么？ | = P(R_t > 1)，稳定 R 分布穿阈值的频率（非新常数） |
| 是否有简单 rhythm？ | 有——`C C R (switch)` 主导：R 游程 ≤ 2（R→R 全程仅 9 次）、C 游程中位 2、反弹后 96–99% 大切换 |
| 是否有精确周期/低维回归？ | **无**（诚实记录） |
| 终端签名？ | 温和——D_max 低谷、无 R 剧变；stop = 普通波动探底穿过 exact-zero 几何 |

**本轮最重要的新发现**：**反弹由边界位置门控**——churn window 内 a<0.55 必收缩（297+1275 轮零反弹）、a≥0.8 必反弹（100%），R_t 随 a 近似单调（`fig_phase_a_vs_R_K100.png`）。

## 2. Churn window 与 R_t 定义

- Window（analysis convention，复用 M1B.2 divider）：onset = 首个 D_max < 0.05 的 active 轮（= 23，4 K 相同），end = t_stop − 1。
- `R_t = max_{j:γ>0} M·max(Δ(j),0)/γ(j)`，γ = D_max − D′_TM；等价性精确（M1C.3），γ>0 逐点成立故无除零边界。

## 3. 数值表

Churn-window 分位数与 switching（见分析文档 §3/§5 的两张表）；终端统计（§7 表）。关键量：R p50 0.689–0.708、p95 3.89–4.49；a p50 0.545–0.560；M p50 0.146–0.152；|Δa| p50 0.38–0.40；switch freq 0.936–0.944；P(switch|rebound prev) 0.964–0.987；RR/RC/CR/CC（K=400）= 2/1050/1049/751。

## 4. Evidence Classification

Model-derived Mathematical Result：R_t ⟺ rebound（既有）。Verified Numerical Analysis：全部逐轮核验与汇总表（确定性，重跑位级一致）。Model-dependent Observation：跨 K 稳定带、a-门控、`C C R` 节奏、终端温和低谷。Working Hypothesis：确定性切换循环的长期行为、a\*≈0.58 分界解析式。**one-step rebound theorem ≠ long-lived churn dynamics**（后者的不收敛性未证明）。

## 5. Reproducibility

- 命令：`.venv/Scripts/python.exe experiments/m1c4_churn_dynamics/run_analysis.py`（1.4 s，deterministic）。
- 测试：`tests/test_m1c4_churn_dynamics.py`（9 项）+ 既有 11 套全过；`python tests/count_tests.py` = 267。
- 结果文件（SHA-256 入 metadata.json）：`per_round_churn_diagnostics.csv`（5344 行）、`churn_window_summary.csv`、`boundary_switch_summary.csv`、`lag_agreement.csv`、`recurrence_summary.csv`、`terminal_summary.csv`、`terminal_alignment.csv`、4 张 PNG、metadata。
- 冻结声明：config.json（window 定义、R_t 定义、switch 阈值 0.2= M1A.1 既有、终端窗口、回归审计参数）运行前写入；无 post-hoc 调整。

## 6. What This Does NOT Prove

见分析文档 §10：不证明长期不收敛/带不变性；不解释 a\* 的解析来源；回归审计是采样上界；不使用 chaos/ergodicity/attractor 语言。依据 Prompt §22 停止条件，四项交付完成后停止，交回 GPT + Owner。
