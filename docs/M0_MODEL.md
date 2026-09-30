# Canonical Model M0 定义（Model Definition）

> 状态：canonical baseline，已实现并运行（见 [实验报告](../experiments/m0_baseline/M0_REPORT.md)）
>
> 建立日期：2026-09-17
>
> 代码：`src/sand_m0/`；实验：`experiments/m0_baseline/`
>
> 证据等级声明：本文档定义的是一个 toy model。所有由 M0 产生的结果只对该模型成立，不得外推为真实沙坑规律或真实最优参数。

## 1. M0 回答的唯一问题

> 一个极简的"近处偏多 → 选择性回收近区 → 按原有偏置重新分配"的过程，本身是否可以产生有限时间的最佳均匀状态？

即检查是否出现轨迹：

`初始不均匀 → 均匀度改善 → 某一有限轮达到最好 → 继续操作后再次恶化`

M0 不研究真实人体操作，不寻找真实参数，不加停止规则（先看完整轨迹，才能判断"先改善再恶化"是否真的存在）。

## 2. 空间定义

一维归一化坐标 `x ∈ [0, 1]`，作为**等面积坐标**使用：

- `x = 0` 表示最靠近操作者，`x = 1` 表示最远端；
- 相同长度的 x 区间代表相同地面面积；
- 因此"真正均匀铺沙"对应 x 上的均匀分布 `Uniform(0,1)`。

M0 不引入真实二维半径、扇角或颗粒动力学（这些属于历史探索与未来模型的范围）。

## 3. 三种 throw distribution

每个实验的初始分布与重泼分布是同一个 `q(x)`。三个分布均为明确、可精确采样的一维 toy distribution（有解析 CDF 与反函数，见 `src/sand_m0/model.py`）：

| 名称 | 密度 | 性质 | 角色 |
|---|---|---|---|
| `q_near` | `2(1-x)` | 近处概率高，越远越低 | 主实验对象（EXP-M0-A） |
| `q_uniform` | `1` | 已完全均匀 | 负对照（EXP-M0-B） |
| `q_far` | `2x` | 远处比近处密 | 反例/边界对照（EXP-M0-C） |

三者不代表真实泼沙数据；`q_near` 与历史上的 `r = U^2`（径向坐标）没有函数关系，只是同类"近密远疏"toy 偏置在本模型坐标下的最简表达。

## 4. 每轮更新规则

固定参数（写入 `experiments/m0_baseline/config.json`，不得事后修改后仍称 baseline）：

- 近区边界 `a = 0.50`（near zone：`x < 0.50`，恰为前 50 个 bin）
- 修正强度 `alpha = 0.25`
- bin 数 `K = 100`
- 轮数 `T = 30`（`t = 0 … 30`），无停止规则

每一轮 `t → t+1`：

1. 找到当前所有 `x < a` 的沙（near zone）；
2. 每一份 near-zone 沙，以概率 `alpha = 0.25` 被独立选中；
3. 被选中的沙从当前位置移除；
4. 每一份被移除的沙，独立地从本实验的 throw distribution `q(x)` 重新抽样落位；
5. 未被选中的沙保持不动；
6. 得到下一轮状态。

对 `alpha` 的解释边界：`alpha` 只表示"每轮 near-zone 中约 25% 的质量被重新分配"，**不是**真实水泥地的 recovery efficiency，更不是现实扫沙回收率。

## 5. 均匀度指标

主指标（primary metric）— histogram squared-L2 distance：

把 `[0,1]` 均分为 `K = 100` 个等宽（=等面积）bin，设 `p_i(t)` 为第 `i` 个 bin 在第 `t` 轮的质量比例，理想均匀值为 `1/K`：

```text
U(t) = Σ_{i=1}^{K} ( p_i(t) − 1/K )²
```

- 数值越小越均匀，完全均匀时 `U = 0`；
- 等价的 MSE 读法为 `U/K`（每 bin 平均平方偏差）；
- 有限样本噪声底线：N 个粒子的完全均匀直方图期望值为 `(K−1)/(K·N)`；N = 100000 时约 `9.9×10⁻⁶`。低于该量级的差异无法与抽样噪声区分。

辅助 sanity-check 指标（secondary，仅用于交叉核对，不用于下结论）— total variation distance：

```text
TV(t) = ½ · Σ_{i=1}^{K} | p_i(t) − 1/K |
```

## 6. 两条独立验证路径

M0 要求每条轨迹同时用两种方式得到：

### A. Deterministic expected-mass version

不模拟单个粒子，直接跟踪 100 个 bin 的理论平均质量。每轮：

- near-zone 各 bin 质量 × `(1 − alpha)`；
- 被移走的总质量 `alpha · M_near(t)` 按 `q` 的 bin 概率精确重新分配；
- 总质量恒为 1。

它没有 Monte Carlo 噪声，代表粒子数无限多时的平均行为。

### B. Monte Carlo particle version

真实随机粒子模拟：`N = 100000`，20 个独立种子（`20260917 … 20260936`，全部记录于 config 与 metadata），`t = 0 … 30`，每个种子独立运行。输出每轮 `U` 的 mean 与 std（跨种子）、每个种子的最好轮次分布、以及与 deterministic 曲线的逐轮差值。

由更新规则可直接推出的解析事实（已写入 `tests/test_m0.py` 并验证）：near-zone 质量 `M_near` 每轮按因子 `(1 − alpha + alpha·Q_near)` 几何衰减，其中 `Q_near = ∫₀^a q(x)dx`。三个实验分别为 0.9375（near）、0.875（uniform）、0.8125（far）。因此 M0 中 near zone 最终必然被掏空——这是一个可以从规则直接读出的结构性结局。

## 7. 三个正式实验

| 实验 | q | 检验问题 |
|---|---|---|
| EXP-M0-A | `q_near` | 是否出现有限 interior minimum（`t_best > 0`），且之后继续操作使 U 再次上升？ |
| EXP-M0-B | `q_uniform` | `t = 0` 是否已最佳？继续选择性回收是否反而破坏均匀性？（负对照） |
| EXP-M0-C | `q_far` | 偏置方向与回收动作不匹配时，A 的改善模式是否消失或显著改变？（反例/边界） |

参数固定为第 4 节的值；无论结果如何都不得为了"漂亮曲线"事后调参并冒充 baseline。若需改变参数，属于下一轮研究，由 GPT + Owner 决定并在研究日志记录。

## 8. 证据等级适用范围

- M0 的 deterministic 轨迹与 MC 汇总在固定代码、参数、种子下可精确复现，因此相关结论可标记 **已验证模拟结果（Verified Simulation Result）**，但**只对 M0 模型成立**；
- "近偏置 + 匹配回收动作才出现 interior minimum"这类跨实验比较结论属于 **模型依赖观察（Model-dependent Observation）**；
- 任何向真实沙坑、真实操作者或文献理论的推广都是 **工作猜想（Working Hypothesis）** 或更低等级，M0 本身不产生现实证据与文献证据。

## 9. 与历史材料的关系

- `experiments/raw_chat_runs/` 与 `archive/` 中的历史探索结果不作为本实验数据，未被修改；
- M0 的最近区比例回收规则与历史 `sand_model_v0` 的"最近固定比例回收"同属 toy 回收策略，但 M0 是独立实现，参数（`a`、`alpha`、`K`、`T`、`N`、seeds）与本模型坐标均与任何历史 run 不同，数值不可与历史 CSV 横向比较；
- M0 属于当时私人研究计划的 Phase 1（canonical baseline）的实现，并顺带覆盖了 Phase 2 中"分布方向敏感性"的最小版本（near/uniform/far 三个极端），完整的分布族比较仍留给后续阶段。
