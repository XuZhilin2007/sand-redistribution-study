# M1C.6 实验报告：Boundary Return & Peak-Competition Dynamics

> 日期：2026-09-19；分析文档：[docs/M1C6_BOUNDARY_RETURN.md](../../docs/M1C6_BOUNDARY_RETURN.md)；前置：[M1C5_BOUNDARY_GATE_THEORY.md](../../docs/M1C5_BOUNDARY_GATE_THEORY.md)
>
> 性质：对既有 canonical Variant A 确定性轨迹的 dynamics audit。无新模型/参数/随机性；模型代码零修改。

## 1. 裁决（反馈回路闭合）

```text
① B>0（门上）→ rebound 可能【M1C.5 定理】→ profile 兑现 → rebound
② rebound 大前缀扫掠 → mismatch reinjection 重塑 profile → argmax 弹回近/中部
   → B 退出 gate（P(B_post<0)=0.93–0.97）
③ contraction：active 峰 (1−α) 压制 + export；competing 峰弱 export + mismatch 补偿
   → peak gap 每轮收缩（gap-closing share = 1.000，全部 4 K）
④ D_max 下降 → a_crit 下降（慢变量）
⑤ boundary 弹跳（|Δa|≈0.4）+ 门缓降 → 1–2 轮内 B 回正（tau_return median 2）
⑥ gate 重开 → 下一次 rebound（tau_reb median 3）
```

**C-C-R 符号节奏 = gate exit → return cycle。** 两个新结构事实：(a) **gate 从不失效**（1958 个 episode 中 never-return = 0）；(b) **switch 主机制随 K 从峰交换（P1 7.8%→0.85%）转向 reinjection 重塑新峰（P2_new 31%→86%）**。

## 2. 关键数字

| K | episodes | P(B_post<0) | return 2 轮 | never | tau_reb med/max | P1 | P_multi | P2_new | gap-closing |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 116 | 0.931 | 0.509 | **0** | 3/5 | 0.078 | 0.560 | 0.310 | **1.000** |
| 100 | 250 | 0.968 | 0.600 | **0** | 3/6 | 0.036 | 0.488 | 0.440 | **1.000** |
| 200 | 530 | 0.962 | 0.545 | **0** | 3/6 | 0.009 | 0.270 | 0.704 | **1.000** |
| 400 | 1058 | 0.962 | 0.546 | **0** | 3/7 | 0.009 | 0.120 | 0.855 | **1.000** |

B 对齐（K=100 中位）：+0.07（s=−2）→ −0.29（s=−1）→ **+0.26（s=0）** → −0.17 → **+0.05（s=2，59.8% 重开）** → −0.10 → −0.03。Boundary map：近端/远端两分支，几乎无点在对角线。

## 3. Evidence Classification

- **Model-derived Mathematical Result**：逐对峰高更新恒等式 `ΔH = ΔH_TM + M[Δ(j1)−Δ(j2)]`（M1C.2 推论，残差 ≤ 2.1×10⁻¹⁵）。
- **Verified Numerical Analysis**：1958 episodes 全表；B 对齐；return/switch 统计；gap-closing 逐轮；terminal 比较；位级复现。
- **Model-dependent Observation**：never-return=0；2 轮 return 主导；gap-closing 100%；P2_new 随 K 增；两分支 map；B 交替 collapse。
- **Working Hypothesis**：repeated gate exit/return through mismatch-rescaled peak competition = 长 churn regime 的生成机制（长期行为无定理）。

## 4. Reproducibility

- 命令：`.venv/Scripts/python.exe experiments/m1c6_boundary_return/run_analysis.py`（2.0 s，deterministic）。
- 测试：`tests/test_m1c6_boundary_return.py`（11 项）+ 既有 13 套全过；`python tests/count_tests.py` = **288**。
- 结果文件（SHA-256 入 metadata）：`episode_table.csv`（1958 行）、`boundary_map.csv`、`gate_return_summary.csv`、`gate_alignment.csv`、`terminal_comparison.csv`、4 张 PNG、metadata.json。
- 冻结声明：config.json（episode/peak/plateau/sep/switch 分类定义）运行前写入；无 post-hoc 调整。

## 5. What This Does NOT Prove

见分析文档 §10：不证明长期不收敛/带不变性；boundary-return 的解析化（tau_return≈2、B 交替）未定理化；terminal 比较为每 K 单样本。依据 Prompt §24 停止。
