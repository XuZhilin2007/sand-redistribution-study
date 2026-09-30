# M1B.2 实验报告：Late-stage Micro-correction Scaling（假说证伪与机制修正）

> 日期：2026-09-19
>
> 分析文档：[docs/M1B_MICRO_CORRECTION.md](../../docs/M1B_MICRO_CORRECTION.md)；前置：[M1A_LONG_HORIZON.md](../../docs/M1A_LONG_HORIZON.md)、[M1B_TOLERANCE.md](../../docs/M1B_TOLERANCE.md)、[M1B_METRIC_SCALE_AUDIT.md](../../docs/M1B_METRIC_SCALE_AUDIT.md)；代码：`src/sand_m0/adaptive.py`（未修改）；配置：[config.json](config.json)（运行前冻结）
>
> 性质：纯诊断/理论解释轮。核心问题：exact stopping time ∝K 与 finite-epsilon stopping 的 K 稳健性，是否因为 exact feedback 在后期进入"每轮动作 O(1/K)"的 micro-correction regime？

## 1. 裁决

**假说被证伪；机制被替换。**

| 待检验预期（假说） | 实测 |
|---|---|
| 每轮 moved_mass 后期收缩到 O(1/K) | **否**：late-stage（D_max<0.005 段）中位数 0.231–0.238，跨 K CV=1.1%，β=+0.007；全局中位数 0.147–0.154（CV≈2%） |
| L1 状态步长收缩到 O(1/K) | **否**：late-stage 0.107–0.117（CV=3.5%，β=+0.031）；全局 0.075（CV≈0.5%） |
| 净输出收缩到 O(1/K) | **部分复杂化**：全局 0.0278（CV=0.4%）K 无关；深谷段 0.0007–0.0015 且 β=−0.216（噪声大、无干净 1/K 标度） |
| exact stop 前动作越来越小 | **停止轮本身是小前缀动作**（moved 0.009–0.060，a=0.035–0.24），但通往它的全部日常动作都是全尺寸 churn |
| finite epsilon 在 micro 区之前停止 | **问题本身不成立**：动作从不缩到 micro；M1B 停止轮 moved=0.224–0.249（moved·K=11–100），与日常 churn 同规模 |

**替代机制（由数据确立）**：宏观快速清理（前 ~100 轮）→ **平稳 churn 波动 regime**（动作恒定、D_max 在 0.003–0.08、U_density 在 p5≈0.0056 的带内不规则震荡）→ **exact stop = 稀有的单轮消灭事件**（一次扫掠同时清零全部正过剩）→ 等待该事件的轮数 ∝K。M1B 的有限 epsilon 把"等待稀有事件"换成"在普通波动谷停手"。

## 2. 逐轮 action 诊断（新数据）

每个活动轮记录：sweep_mass、moved_mass、net_export_mass（=moved×(1−F_q(a))，与 canonical 更新精确一致的恒等式）、L1_state_step、delta_U_density（带符号）、delta_U_per_moved。文件：`action_diagnostics_per_round.csv`（K=50: 324 行 / K=100: 699 / K=200: 1445 / K=400: 2876）。

恒等式验证（测试）：moved = alpha·sweep_mass ✓；前缀质量变化 = −net_export ✓（从 states 直接复核）；L1 ≥ net_export ✓。

## 3. 相位与 D_max 分箱

- **相对相位分箱**（0–25/25–50/50–75/75–90/90–99/99–100%）：各相位中位数 moved 0.14–0.16、net_export ≈0.027、L1 ≈0.075——无明显趋势（`phase_summary.csv`）。
- **D_max 分箱**（<0.005 / 0.005–0.01 / 0.01–0.02 / 0.02–0.05 / ≥0.05）：moved 与 L1 的中位数在各带几乎不变；net_export 随带变化大（高带 ~0.03，深谷带 ~0.001，由 (1−a)² 因子驱动）（`D_max_band_summary.csv`）。

## 4. 跨 K scaling

| 量（late-stage 中位数） | K=50 | K=100 | K=200 | K=400 | CV(X) | CV(K·X) | CV(K²·X) | β |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| moved_mass | 0.2357 | 0.2312 | 0.2374 | 0.2377 | 0.011 | 0.72 | 1.19 | +0.007 |
| net_export_mass | 0.000849 | 0.001478 | 0.000718 | 0.000655 | 0.354 | 0.521 | 1.06 | −0.216 |
| L1_state_step | 0.1106 | 0.1070 | 0.1170 | 0.1153 | 0.035 | 0.73 | 1.19 | +0.031 |

**X 本身 collapse（CV 最小），K·X 与 K²·X 都不 collapse**——三个量都不是 O(1/K)。4 点拟合仅作诊断（β 不确定度大，尤其 net_export）。

## 5. 尾巴长度与归一化时间

| K | t_U_best/K | t_stop/K | tail/K（≤0.02 后） | tail/K（≤0.01 后） | tail/K（≤0.005 后） |
|---:|---:|---:|---:|---:|---:|
| 50 | 5.84 | 6.48 | 4.96 | 3.18 | 1.42 |
| 100 | 6.53 | 6.99 | 6.19 | 5.40 | 3.59 |
| 200 | 6.75 | 7.23 | 6.85 | 6.36 | 5.54 |
| 400 | 6.57 | 7.19 | 6.99 | 6.77 | 6.39 |

归一化时间趋于常数（t_U_best/K≈6.5–6.7；t_stop/K≈6.5–7.2）；**经过 D_max=0.005 后的尾巴本身 ≈ O(K)**（1.42→6.39 ×K）——micro-scale 存在于等待时间。

## 6. Epsilon overlay

16 个 (K, epsilon) 首达点：moved_mass = 0.224–0.249（全部 gross churn，moved·K = 11–100）；retention 0.912–0.991。**M1B 停止轮就是普通的全尺寸 churn 轮，恰好落在 D_max 波动谷里**。D_max ≤ eps 的低谷频次（每 1000 活动轮）：eps=0.005 → 27.8/58.7/74.0/91.1（K=50/100/200/400，温和递增）；eps=0.01 → 111/143/183/220；首次低谷时刻跨 K 稳健（0.005：253/340/338/321）。

## 7. Annihilation 事件（exact stop 的真实形态）

| K | t_stop | pre-stop a | F(a) | F(a)≤0.5 | moved | 正过剩前缀 pre→post |
|---:|---:|---:|---:|---|---:|---|
| 50 | 324 | 0.24 | 0.422 | ✓ | 0.0604 | 8 → 0 |
| 100 | 699 | 0.07 | 0.135 | ✓ | 0.0182 | 25 → 0 |
| 200 | 1445 | 0.195 | 0.352 | ✓ | 0.0511 | 123 → 0 |
| 400 | 2876 | 0.035 | 0.069 | ✓ | 0.0089 | 28 → 0 |

消灭扫掠是小前缀动作；充要条件 (i)/(ii) 与必要条件 F(a)≤0.5 的推导见 [M1B_MICRO_CORRECTION.md](../../docs/M1B_MICRO_CORRECTION.md) §5。

## 8. Evidence Classification

- **Model-derived Mathematical Result**：net_export 恒等式（前缀质量变化 = −moved×(1−F_q(a))）；边界外单调性；消灭事件充要条件 (i)/(ii) 与 F(a)≤0.5 必要条件推论。
- **Verified Simulation Result**：逐轮 action 诊断表（4 K）；churn 中位数跨 K CV≤2%；D_max/U_density 波动带分位数；尾巴长度表；epsilon overlay 16 点；消灭事件 4/4 结构；回归锚点（M1B 20/20、M1A.1 t=0..100）。
- **Model-dependent Observation**：平稳 churn regime 与两期划分；波动带 K 稳健性；tail ∝K；消灭扫掠为小前缀动作；深谷段 net_export 变小且无干净标度。
- **Working Hypothesis**：消灭等待 ∝K 的概率机制（需同时满足的格点约束随 K 增长）；连续极限解读。
- 无新文献声明；无 Real-world calibrated result。

## 9. Unexpected / Negative Results

1. **micro-correction 假说被证伪**——本轮预设其主要检验对象，结果毛动作规模全程恒定（这是本轮最有价值的结果）。
2. **delta_U_density 晚期正负兼有且量级不小**（K=100：每轮 −0.11…+0.04）——U 的"最优点"是波动带内的有利采样，不是趋势终点；强化 M1B.1 的 retention 解释。
3. **消灭扫掠是小前缀动作**（moved 0.009–0.060）而通往它的动作是全尺寸 churn——"最后一击"与"日常动作"规模完全不同。
4. net_export 深谷段无干净标度（β=−0.216）——(1−a)² 因子随 argmax 摆动主导，如实报告为噪声/无标度。
5. 4/4 消灭事件满足 F(a)≤0.5——必要条件得到验证（而非证伪项）。

## 10. What This Does NOT Prove

- 不证明消灭等待 ∝K 的解析规律（4 点 + 假说）；
- 不证明连续极限 exact stopping 不可达；
- 不证明平稳 churn regime 的遍历性/混合性/混沌性（无此类分析）；
- 不证明波动带 K 稳健性在其他 q/alpha/控制器下保持；
- 不引入任何 controller 修改（perceptual-B、hysteresis、noisy perception、U-derivative stopping 等全部留待 Owner + GPT review）。

## 11. Reproducibility

- **Python** 3.14.6（仓库 `.venv`）；**依赖** numpy 2.5.3、matplotlib 3.11.2。
- **命令**：`.venv/Scripts/python.exe experiments/m1b_2_micro_scaling/run_experiment.py`
- **测试**：`tests/test_m1b_2_micro_scaling.py`（13 项）+ 既有七套（37+25+54+18+18+20+17）全部通过。
- **随机性**：deterministic，无种子。
- **结果文件**（`experiments/m1b_2_micro_scaling/results/`，SHA-256 见 metadata.json）：`action_diagnostics_per_round.csv`（5344 行）、`phase_summary.csv`、`D_max_band_summary.csv`、`cross_K_scaling.csv`、`band_and_dip_statistics.csv`、`tail_lengths.csv`、`epsilon_overlay.csv`、`annihilation_events.csv`、四张 PNG、`metadata.json`。
- **可复现实证**：干净树上重跑，数据输出位级一致。
- **冻结声明**：相位分箱、D_max 带、tail 阈值、eps overlay 集合运行前写入 config.json；模型与 controller 未改。

## 12. Recommended Next Research Question

机制链已闭合到"稀有消灭事件"，最有价值的下一问（待 Owner + GPT 决定）：

1. **消灭事件的形状统计**：characterize "可消灭状态"（全部正过剩限于短近前缀）在波动中的出现频率随 K 的增长规律——把 t∝K 从观察推进为可计算的到达时间模型；
2. 或 **Phase 2B**（Beta family / r=U²）闭合分布维度（M0 判据、M1B 容差、消灭机制三个结论的分布稳健性）。
