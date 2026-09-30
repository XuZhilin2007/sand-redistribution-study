# M1C.5 实验报告：Boundary-Gated Rebound Theory

> 日期：2026-09-19；分析文档：[docs/M1C5_BOUNDARY_GATE_THEORY.md](../../docs/M1C5_BOUNDARY_GATE_THEORY.md)；前置：[M1C4_CHURN_REGIME_AUDIT.md](../../docs/M1C4_CHURN_REGIME_AUDIT.md)
>
> 性质：theory + existing-trajectory verification。无新模型/参数/实验；canonical Variant A 确定性轨迹（K∈{50,100,200,400}）。

## 1. 裁决（Outcome B1）

**M1C.4 的经验边界门被解析推导为必要条件，闭式门槛定量命中观察门。**

- **定理链（全部精确、逐轮核验）**：Region-1 化简 `M·Δ − γ = (1−α)D − D_max + α·h`（残差 ≤ 4.0×10⁻¹⁶）→ state-reduced 上界 `≤ α(max_j h − D_max)`（0 违反）→ h 最大化（F_a≤1/2 ⟹ 精确不可能；F_a>1/2 ⟹ 内部 max (2F_a−1)²/(4F_a)，x\*≤a 条件在 canonical 全程成立）→ **必要门槛 `a > a_crit(D_max)`，a_crit(D) = [(1−D)+√(D(D+2))]/2**。
- **经验验证**：全部 **1958 个 rebound 轮满足 a > a_crit（FN = 0）**；a_crit churn-window p50 = **0.593–0.604**（跨 K 稳定）、p5–p95 = 0.539–0.644——解析门槛 tracks the observed ≈0.58 gate。
- **必要 ≠ 充分**：gate-open 但收缩的轮 21/52/88/177 个；其 P(x_r) = D(x_r)/D_max 分布与 rebound 轮显著分离——**边界开启可能门，profile 形状决定是否反弹**。

## 2. 逻辑方向（防误读）

```text
a ≤ a_crit(D_max)  ⟹  rebound 不可能          【已证明】
a > a_crit(D_max)  ⟹  rebound 可能，未必发生    【necessary only】
```

不声称 "a > a_crit ⟹ rebound"；不声称 0.58 是定理值（a_crit 依赖 D_max，其分布跨轮变化）。

## 3. 数值汇总

见分析文档 §7 两张表（a_crit 分位数 + 检查汇总）与三张图（`fig_gate_scatter.png`：a–D 平面 + a_crit 曲线，rebound 全部在曲线上方；`fig_gate_margin_B.png`：B≤0 纯收缩、B>0 混合；`fig_profile_factor.png`：P 因子分离）。

## 4. Evidence Classification

- **Model-derived Mathematical Result**：引理 1（Region-1 化简）、引理 2（state-reduced 上界与必要条件）、引理 3（解析最大化，含精确惰性分支）、定理（闭式门槛）。
- **Verified Numerical Analysis**：5344 轮全链核验；FN=0；a_crit 分位数表；确定性复现（SHA-256）。
- **Model-dependent Observation**：a_crit 带跨 K 稳定于观察门区域；gate-open contraction 占比 13–16%；P 因子分离。
- **Working Hypothesis**：边界反复穿越解析门产生长期 churn switching（未证明）。

## 5. Reproducibility

- 命令：`.venv/Scripts/python.exe experiments/m1c5_boundary_gate/run_analysis.py`（1.3 s，deterministic；重跑输出位级一致）。
- 测试：`tests/test_m1c5_boundary_gate.py`（10 项）+ 既有 12 套全过；`python tests/count_tests.py` = **277**。
- 结果文件：`per_round_gate.csv`（5344 行）、`gate_summary.csv`、3 张 PNG、metadata.json（SHA-256）。
- 冻结声明：config.json（理论链、门槛闭式、churn window、P 因子定义与 D_max>0.005 过滤）运行前写入；无 post-hoc 调整。

## 6. What This Does NOT Prove

见分析文档 §11：必要条件非充分；x\*≤a 条件依赖 D_max ≤ √2−1（canonical 满足）；closed-form 为 continuum-style envelope（两种离散 rounding regime 已记录，离散以格点 max 为准）；不证明 boundary-return dynamics 与长期 churn。依据 Prompt §21 停止。
