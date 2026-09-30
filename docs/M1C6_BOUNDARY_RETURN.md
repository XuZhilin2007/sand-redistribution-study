# M1C.6 诊断文档：Boundary Return & Peak-Competition Dynamics——闭合确定性 churn 反馈回路

> 状态：dynamics audit 轮（2026-09-19）。M1C.5 的直接后续；**无新模型、无新参数、无随机性**——只分析既有 canonical Variant A 确定性轨迹（K∈{50,100,200,400}）。
>
> 实验：[M1C6_REPORT.md](../experiments/m1c6_boundary_return/M1C6_REPORT.md)；前置：[M1C5_BOUNDARY_GATE_THEORY.md](M1C5_BOUNDARY_GATE_THEORY.md)、[M1C4_CHURN_REGIME_AUDIT.md](M1C4_CHURN_REGIME_AUDIT.md)、[M1C3_REDISTRIBUTION_LAW_ORDER.md](M1C3_REDISTRIBUTION_LAW_ORDER.md)
>
> North Star：解释反馈回路缺失的一半——**rebound + boundary switch 之后，是什么让 sweep boundary 在 1–2 个 contraction 轮内重新回到解析反弹门之上？**

## 0. 裁决（Outcome R1 + R4 弱版：反馈回路闭合）

**回路闭合。** 每一半都有精确方程 + 轨迹证据支撑：

```text
① 边界在门上（B>0）→ rebound 可能（M1C.5 定理）→ profile 兑现 → rebound
② rebound 扫掠（大前缀）→ mismatch reinjection 重塑 profile
   → argmax 大幅弹回近/中部（|Δa|≈0.4）→ B 强制退出 gate（P(B_post<0)=0.93–0.97）
③ contraction 轮：active 峰被 (1−α) 压制 + export，competing 峰只受弱 export
   且被 mismatch 项补偿 → peak gap 每轮收缩（gap-closing share = 1.000）
④ 同时 D_max 下降 → a_crit(D) 下降（慢变量）
⑤ boundary 弹跳 + 门缓降 → 1–2 轮内 B 回正（tau_return median 2；never = 0.0000）
⑥ gate 重开 → profile 兑现 → 下一次 rebound（tau_reb median 3）
```

**C-C-R 符号节奏就是 gate exit → return cycle**。新发现两个结构性事实：

1. **B 逐轮交替震荡**：对齐后 B 的中位路径为 +0.07（s=−2）→ **−0.29**（s=−1）→ +0.26（s=0，rebound）→ −0.17（s=1）→ +0.05（s=2，59.8% 重开）→ −0.10 → −0.03——boundary 不是渐进漂移回门，而是**每轮大幅弹跳**（boundary map 近端/远端两分支，几乎无点在对角线上），B 随弹跳交替变号，a_crit 随 D_max 缓降是慢变量。
2. **switch 的主机制随 K 变化**：既有 secondary peak 接管（P1）占 7.8%（K=50）→ **0.85%（K=400）**；reinjection 重塑出的新位置（P2_new）占 31% → **86%**；第三峰接管（P_multi）56% → 12%。即**大 K 下 switch 不是两峰交换，而是 mismatch reinjection 对整个 profile 的重新塑形**。

## 1. Episode 定义与总体统计

Episode = rebound transition t_r → t_r+1 开启，至下一 rebound 前（terminal episode 至 stop）。扫掠边界 a_{t_r} = state t_r 的 argmax。全部 4 K 合计 **1958 个 rebound episodes**（K=50/100/200/400：116/250/530/1058 个非 terminal + 1 terminal）。

## 2. Gate exit → return（Question B）

| K | episodes | P(B_post<0) | return 1 轮 | 2 轮 | ≥3 轮 | **never** | tau_reb median/max |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | 116 | 0.9310 | 0.069 | 0.509 | 0.422 | **0.0000** | 3 / 5 |
| 100 | 250 | 0.9680 | 0.032 | 0.600 | 0.368 | **0.0000** | 3 / 6 |
| 200 | 530 | 0.9623 | 0.038 | 0.545 | 0.417 | **0.0000** | 3 / 6 |
| 400 | 1058 | 0.9622 | 0.038 | 0.546 | 0.416 | **0.0000** | 3 / 7 |

- rebound 后下一轮 B 中位 −0.17（93–97% 在门下）——**gate exit 是 rebound 的必然后果**（D_max 被抬高 ⟹ a_crit 上升 + argmax 弹离大前缀）；
- **gate 从不失效**：1958 个 episode 中 0 个 never-return（gate 在 stop 前总会重开）；
- return 集中在 2 轮（51–60%），3+ 轮 37–42%；
- tau_rebound median 3（p25=2）：gate 重开（s=2）后下一次 rebound 通常立即或隔 1 轮发生。

（empirical deterministic episode frequencies——非 stochastic model。）

## 3. Peak competition（§7/§9/§10）

**精确逐对峰高更新恒等式（T1，Model-derived Mathematical Result）**：对任意固定 index pair (j1, j2)，由 M1C.2 分解直接得到

```text
D_new(j1) − D_new(j2) = [D_TM_new(j1) − D_TM_new(j2)] + M[Δ(j1) − Δ(j2)]
```

top-2 meaningful peaks 上逐轮核验，最大残差 **2.1×10⁻¹⁵**（4 K）。

**Contraction 轮的相对运动**（rebound 前后各 contraction 轮，gap-closing share）：

| K | d(active peak) 中位 | d(competing peak) 中位 | gap 收缩轮占比 |
|---:|---:|---:|---:|
| 50 | −0.0361 | −0.0181 | **1.0000** |
| 100 | −0.0355 | −0.0214 | **1.0000** |
| 200 | −0.0356 | −0.0224 | **1.0000** |
| 400 | −0.0354 | −0.0228 | **1.0000** |

**机制（由 exact equations 支持）**：active peak j1 = 扫掠边界自身，受双重压制——(1−α) 缩放 + net export（Region 1 公式）；competing peak j2 在 suffix，只降 M(1−T(j2)) 且被 mismatch 项 M·Δ(j2) 部分补偿（secondary peak 集中在 Δ 峰区 x≈0.5–0.6，全 churn 轮 secondary 位置 x 中位 0.57–0.58）。净效果：active 每轮多降 ~0.013–0.018，gap 单调收缩，1–2 轮内换位/弹回。

## 4. Rebound-induced switch 的真实来源（§8）

| K | P1（既有 secondary 接管） | P_multi（第三峰接管） | P2_new（reinjection 重塑新峰） | no_switch |
|---:|---:|---:|---:|---:|
| 50 | 0.0776 | 0.5603 | 0.3103 | 0.0517 |
| 100 | 0.0360 | 0.4880 | 0.4400 | 0.0360 |
| 200 | 0.0094 | 0.2698 | 0.7038 | 0.0170 |
| 400 | 0.0085 | 0.1200 | **0.8554** | 0.0161 |

- **P1（简单两峰交换）是少数且随 K 衰减**；P2_new 随 K 强增。rebound 大前缀扫掠的 mismatch 注入 M·Δ（Δ 峰在 x=1/2）系统性地把 profile 重新塑形，新 argmax 多出现在**扫掠后 profile 的新位置**——Prompt §8 的 candidate mechanism 以 **P2/mixed 为主**，纯 P1 两峰交换图像在大 K 不成立。
- no_switch 1.6–5.2% 与 M1C.4 的 P(switch|rebound)=0.96–0.99 一致。
- operational definitions（peak、sep=max(2,K//50)、plateau→首索引、P1/P_multi/P2 分类）冻结于 config。

## 5. Boundary map（§12）

a_t → a_{t+1} 散点呈**两分支结构**（`fig_boundary_map_K100.png`）：近端分支（a_t∈[0.05,0.65] → a_{t+1}∈[0.6,1.0]，几乎全部 contraction）与远端分支（a_t∈[0.55,1.0] → a_{t+1}∈[0.05,0.68]，含全部 rebound 轮）；**几乎无点在对角线上**——boundary 每轮在近/远两区之间大幅弹跳（|Δa| 中位 0.39–0.40），从不停留。rebound 轮从门上（a>a_crit）发起、落点中位 a≈0.46。

## 6. 可预测性（§13，简单 diagnostics）

B 对齐的跨 K collapse（图 `fig_gate_episode_alignment.png`：四个 K 的中位路径与四分位带几乎重合）+ 两分支 boundary map 表明低维结构存在；但 B 的逐轮交替由弹跳相位主导（a 跳 |Δa|≈0.4 ≫ a_crit 慢漂 ~0.04/轮），**单靠 (a_t, D_max) 不决定 a_{t+1} 的分支归属**（两分支在同 a 值区间重叠）。未训练任何 predictor；未做高阶拟合。

## 7. Terminal episode（§15）

| K | terminal len | tau_return | episode 内 B>0 轮 | ordinary len median | ordinary tau_return median |
|---:|---:|---:|---:|---:|---:|
| 50 | 2 | 2 | 0 | 3 | 2 |
| 100 | 5 | 3 | 1 | 3 | 2 |
| 200 | 2 | 2 | 0 | 3 | 2 |
| 400 | 4 | 2 | 1 | 3 | 2 |

terminal episode 与普通 episode 结构相同（len、tau_return 同量级）；K=50/200 的 terminal 中 gate 未及重开即停止，K=100/400 中 gate 重开过 1 轮但未兑现 rebound——**exact stopping = 一次普通 episode 恰好进入 annihilable geometry（D_max 波动谷 + 小前缀消灭构型，M1B.2），而非 gate 机制的失效**（Outcome R4 弱版本）。

## 8. Evidence Classification

- **Model-derived Mathematical Result**：逐对峰高更新恒等式（M1C.2 分解的直接推论；残差 ≤ 2.1×10⁻¹⁵）。无新不等式定理（T2/T3 未强求）。
- **Verified Numerical Analysis**：episode 分割与 1958 个 episode 表；B 对齐统计；tau_return/tau_rebound 分布；gap-closing 逐轮核验；switch 机制分类；boundary map；terminal 比较。确定性复现（输出位级一致）。
- **Model-dependent Observation**：gate never-return = 0；return 集中于 2 轮；gap-closing 100%；P2_new 随 K 增、P1 随 K 减；boundary map 两分支；B 交替震荡跨 K collapse。
- **Working Hypothesis**：**repeated gate exit/return through mismatch-rescaled peak competition generates the observed long churn regime**——反馈回路的两半（M1C.5 门控 + 本轮 return）均有方程与数据支撑，长期行为（永不提前停止、带不变性）仍无定理。

## 9. Reduced-model viability（§23 G）

**Qualified yes**：两分支 boundary map + gap 单调收缩 + a_crit 慢漂给出清晰的 reduced 描述骨架；但"两峰"实为 **active peak vs mismatch reinjection 在 Δ 峰区重塑的竞争结构**（P2/P_multi 主导，纯两峰交换仅 0.8–7.8%），reduced model 必须包含 M·Δ(x) 的 reinjection 几何，而非简单两峰动力学。是否值得建——交 GPT + Owner。

## 10. Limitations

1. 单一 q/alpha/一维/deterministic/4 K；episode 统计为确定性频率（非概率模型）；
2. peak/sep/plateau/机制分类为冻结的 operational conventions（config 已记录）；3+ 峰情形归入 P_multi；
3. 可预测性只做了 binwise/duplicate 级 diagnostics，无系统建模；
4. B 的交替震荡与弹跳相位的解析关系未推导；tau_return=2 的主导性未定理化；
5. terminal 比较为 4 个样本（每 K 1 个 terminal episode）；
6. 不声称 limit cycle/attractor/chaos/ergodicity/infinite rebounds（Prompt §17 遵守；有限 K 已 finite-time stop）。

## 11. Next Research Gate（交回 GPT + Owner）

1. **two-peak + reinjection reduced model**（qualified yes 的具体化）；
2. tau_return≈2 与 B 交替的解析化（快弹跳 + 慢门漂的两时间尺度近似）；
3. 长期化：门控（M1C.5）+ return（本轮）合并不收敛定理的可行性；
4. 依据 Prompt §24，四项交付（gate-return、peak-competition、boundary map、terminal comparison）完成后停止。

## 12. Repository 资产

- `experiments/m1c6_boundary_return/`：config / run_analysis.py / episode_table.csv（1958 行）/ boundary_map.csv / gate_return_summary.csv / gate_alignment.csv / terminal_comparison.csv / 4 图 / M1C6_REPORT.md；
- `tests/test_m1c6_boundary_return.py`（11 项）；
- plotting.py 新增 4 函数；模型代码零修改。
