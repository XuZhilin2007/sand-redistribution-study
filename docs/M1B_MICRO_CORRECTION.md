# M1B.2 诊断文档：Late-stage Micro-correction Scaling——假说检验与机制修正

> 状态：纯诊断/理论解释轮（2026-09-19）。M1A dynamics、adaptive boundary policy、epsilon stopping rule、q、alpha、metric、controller **全部未改**。
>
> 实验：[M1B_2_REPORT.md](../experiments/m1b_2_micro_scaling/M1B_2_REPORT.md)；前置：[M1A_LONG_HORIZON.md](M1A_LONG_HORIZON.md)、[M1B_TOLERANCE.md](M1B_TOLERANCE.md)、[M1B_METRIC_SCALE_AUDIT.md](M1B_METRIC_SCALE_AUDIT.md)
>
> 待检验假说（Working Hypothesis，本轮**主动证伪**）：
>
> > exact feedback 在后期进入 spatial-discretization-scale micro-correction regime：每轮实际修正质量越来越小（O(1/K)），因此 exact stopping time 随 K 增长；positive epsilon 在该 regime 之前截断过程。

## 0. 裁决

**假说被证伪。** 逐轮 action 诊断显示：从 t=0 到 exact stop，controller 的**毛动作规模完全不缩水**——每个 K 下每轮搬动约 15% 的总沙量（moved_mass 中位数 0.147–0.154，跨 K CV≈2%），状态 L1 步长 ≈0.075（CV≈0.5%），净输出 ≈0.028（CV≈0.4%）。不存在"动作收缩到 O(1/K)"的 micro-correction regime。

**替代机制（由数据确立）**：exact controller 在快速宏观清理（~前 100 轮）后进入**平稳 churn 波动 regime**——D_max 在 K 稳健的带区（约 0.003–0.08）内不规则震荡，U_density 在其自身带区（p5≈0.0056，全 K）内震荡；**exact stopping 是"一次扫掠同时清零全部正过剩"的稀有消灭事件（annihilation event）**，等待时间随 K 增长（≈∝K）。M1B 的有限 epsilon 之所以 K 稳健且 retention 高，是因为它**不再等待稀有消灭事件，而是在普通波动谷（D_max ≤ epsilon 的低谷）停手**。

以下分节给出证据。所有数值来自 canonical dynamics 的 exact reference 轨迹（K=100 前 101 行与已提交 M1A.1 轨迹一致、停止轮与已提交 M1B 结果 20/20 一致——回归锚点测试通过）。

## 1. 证伪证据：毛动作规模不随 K 缩水

**全局 churn 中位数（全部活动轮）**：

| 量 | K=50 | K=100 | K=200 | K=400 | 跨 K CV |
|---|---:|---:|---:|---:|---:|
| moved_mass | 0.1544 | 0.1502 | 0.1480 | 0.1466 | ≈2% |
| net_export_mass | 0.0278 | 0.0278 | 0.0277 | 0.0279 | ≈0.4% |
| L1_state_step | 0.0757 | 0.0751 | 0.0751 | 0.0747 | ≈0.5% |

**late-stage（D_max < 0.005 深谷段）中位数**：

| 量 | K=50 | K=100 | K=200 | K=400 | log-log 斜率 β |
|---|---:|---:|---:|---:|---:|
| moved_mass | 0.2357 | 0.2312 | 0.2374 | 0.2377 | +0.007 |
| L1_state_step | 0.1106 | 0.1070 | 0.1170 | 0.1153 | +0.031 |
| net_export_mass | 0.000849 | 0.001478 | 0.000718 | 0.000655 | −0.216 |

- moved_mass 与 L1：β ≈ 0（K 无关），CV_X ≪ CV_{K·X}——**X 本身 collapse，K·X 不 collapse**，与 1/K 假说（β=−1）相反。
- net_export 全局 K 无关；深谷段内变小（~0.001）且噪声大、β=−0.216 远离 −1——无干净 1/K 标度。结构性原因见 §5：净输出 = alpha·F(a)·(1−F_q(a))，深谷轮的 argmax 前缀在大小间摆动（a=0.06–0.99），(1−a)² 因子使其在 0.0002–0.045 间大幅摆动。
- **注意**：4 个 K 点的拟合不是定理；β 与 CV 仅作诊断。

## 2. 替代机制：平稳 churn + K 稳健波动带

活动段（exact stop 之前）的 D_max 与 U_density 分位数：

| K | D_max p5/p50/p95 | U_density p5/p50/p95 |
|---:|---|---|
| 50 | 0.0066 / 0.0293 / 0.0783 | 0.00595 / 0.01315 / 0.0658 |
| 100 | 0.0046 / 0.0251 / 0.0596 | 0.00554 / 0.01083 / 0.0386 |
| 200 | 0.0038 / 0.0226 / 0.0531 | 0.00564 / 0.00942 / 0.0234 |
| 400 | 0.0033 / 0.0212 / 0.0510 | 0.00556 / 0.00864 / 0.0182 |

- **D_max 波动带跨 K 稳健**（p50 0.021–0.029，p5 0.0033–0.0066）——不是逐渐收紧的收敛带，而是长期驻留的震荡带；
- **U_density 带的下沿 p5 ≈ 0.0056 跨 K 几乎不变**——这直接解释 M1B.1 的发现（停止态 retention 高且 K 稳健）：波动带的下限就是"停止能拿到"的水平；
- 每轮 delta_U_density 在晚期**正负兼有**（K=100：−0.11 … +0.04），U 的"最优点"是波动带内的有利采样，不是趋势终点。

**regime 描述（只描述、不定义 sharp transition）**：宏观清理期（t≲100，delta_U 大且为正）→ 平稳 churn 波动期（其后全部轮次：动作恒定、D_max/U 震荡）。操作上可用"phase 相对位置"或"D_max 首次 <0.05"近似划分，但两者间无明确突变点。

## 3. exact stopping 的真实机制：稀有消灭事件

exact stop 轮的前一轮（pre-stop）诊断（4/4 一致）：

| K | t_stop | pre-stop a_t | F(a) | F(a)≤0.5 | moved | 正过剩前缀 pre→post |
|---:|---:|---:|---:|---|---:|---|
| 50 | 324 | 0.24 | 0.422 | ✓ | 0.0604 | 8 → 0 |
| 100 | 699 | 0.07 | 0.135 | ✓ | 0.0182 | 25 → 0 |
| 200 | 1445 | 0.195 | 0.352 | ✓ | 0.0511 | 123 → 0 |
| 400 | 2876 | 0.035 | 0.069 | ✓ | 0.0089 | 28 → 0 |

- **一次扫掠清零全部正过剩**（post 全为 0）；pre-stop 的 D_max 只有 6.5×10⁻⁴–9.4×10⁻³；
- 消灭扫掠本身是**小前缀动作**（moved 0.009–0.060，远小于日常 churn 的 ~0.15）；
- 所有事件满足必要条件 **F(a) ≤ 0.5**（推导见 §5；直观：q_near 在 x≈0 的密度是均匀的 2 倍，扫大前缀会把近端重新填到均匀之上，故清零必须用落在近半区的小前缀扫掠）。

**为什么等待时间 ∝K（Working Hypothesis，未证明）**：消灭事件要求波动把状态送入一个特定形状——全部正过剩集中在满足 §5 包络条件的短近前缀内。D_max 的**带底**（~0.003）跨 K 稳健，但"带内波动恰好给出可消灭形状"的每轮概率随 K 下降（需同时满足的约束点数随 K 增长）→ 等待轮数 ∝K。

## 4. 尾巴长度：微尺度体现在等待时间，不在动作规模

| K | t_U_best/K | t_stop/K | tail/K（D_max≤0.02 后） | tail/K（≤0.01 后） | tail/K（≤0.005 后） |
|---:|---:|---:|---:|---:|---:|
| 50 | 5.84 | 6.48 | 4.96 | 3.18 | 1.42 |
| 100 | 6.53 | 6.99 | 6.19 | 5.40 | 3.59 |
| 200 | 6.75 | 7.23 | 6.85 | 6.36 | 5.54 |
| 400 | 6.57 | 7.19 | 6.99 | 6.77 | 6.39 |

- `t_U_best/K ≈ 6.5–6.7`、`t_stop/K ≈ 6.5–7.2`——归一化后趋于稳定常数（T∝K 的规整化表述）；
- **经过 D_max=0.005 之后到 exact stop 的尾巴本身 ≈ O(K)**（1.42×K → 6.39×K，随 K 增长趋稳）——micro-scale 确实存在于**等待时间**中；
- 对照：毛动作规模全程不缩（§1）——尾巴长的不是动作小，而是**等得到那个特殊形状**。

## 5. 数学分析（可严格推导的部分）

**恒等式 1（净输出，精确）**：一次扫掠后前缀质量变化

```text
M_prefix(t+1) − M_prefix(t) = −alpha·M_prefix(t)·(1 − F_q(a_t)) = −net_export
```

即 `net_export = alpha·F(a_t)·(1 − F_q(a_t))`，其中 argmax 前缀满足 `F(a_t) = a_t + D_max(t)`。对 q_near：`1 − F_q(a) = (1−a)²`，故 `net_export = alpha·(a + D_max)·(1−a)²`——a→1 或 a→0、或 D_max→0 都使其变小；这解释其 0.0002–0.045 的大幅摆动（已数值验证恒等式）。

**恒等式 2（边界外单调，承 M1A.1）**：`D_{t+1}(x) = D_t(x) − alpha·F(a_t)·(1 − F_q(x)) ≤ D_t(x) ∀ x > a_t`。

**消灭事件的充要条件（精确）**：一次边界为 a 的扫掠实现清零（∀x: D'(x) ≤ 0）当且仅当

```text
(i)   ∀ x > a:  D(x) ≤ alpha·F(a)·(1 − F_q(x))          （边界外：过剩须被本次输出覆盖）
(ii)  ∀ x ≤ a:  (1−alpha)·D(x) + alpha·(F(a)·F_q(x) − x) ≤ 0   （边界内：回填后不得转正）
```

**推论（必要条件）**：q_near 在 x→0⁺ 处 `F_q(x) ≈ 2x`；若近端存在正过剩（小 x 处 D(x) ≥ 0），则 (ii) 在小 x 要求 `x·(2F(a) − 1) ≲ −(1−alpha)/alpha·D(x) ≤ 0`，即 **F(a) ≤ 0.5**——扫掠边界必须落在近半区。4 个观测事件全部满足（F(a) = 0.069–0.422）。

**未证明**：消灭事件等待时间 ∝K 的解析证明（需要 D-profile 形状的波动统计，超出本轮）；连续极限下 exact stopping 是否有限时间可达。

## 6. Epsilon overlay：有限容差停在普通波动谷

16 个 (K, epsilon) 首达点处的动作诊断：moved_mass = 0.224–0.249（**全部是日常 churn 规模**，moved·K = 11–100 ≫ 1），retention 0.912–0.991。即：**M1B 的停止轮本身也是一次全尺寸 churn 动作**——它停 in a dip（D_max ≤ epsilon 的波动低谷），不是在动作缩到微观时停。

D_max ≤ eps 的低谷在活动段内多次出现（eps=0.005：9–262 次；eps=0.04：226–2416 次），首次出现时刻跨 K 稳健（0.005：253/340/338/321）——这正是 M1B 停止轮 K 稳健的机制。

## 7. 与 deadband 及 M1B.1 retention 的关系

- **修正后的 deadband 解释（Model-dependent Observation / Working Hypothesis）**：M1B 的有限 epsilon 不是"截断 micro-correction 尾巴"（该尾巴不存在于动作规模），而是"**不再等待稀有消灭事件，接受波动带下沿**"。停止态落在 U_density 带下沿（p5≈0.0056）附近，故 retention 高（M1B.1）。
- **不声称**：deadband universally regularizes all adaptive redistribution systems——本轮只确立 M1A/M1B 这一个系统内的机制。

## 8. Evidence Classification

- **模型推导数学结果**：net_export 恒等式；边界外单调性；消灭事件充要条件 (i)/(ii) 及 F(a) ≤ 0.5 必要条件推论——全部由更新规则严格推出并数值核验。
- **已验证模拟结果**：四 K 的逐轮 action 诊断表、churn 中位数的 K 无关性（CV≤2%）、D_max/U 波动带的跨 K 分位数、消灭事件 4/4 结构、尾巴长度表、epsilon overlay 16 点。
- **模型依赖观察**：平稳 churn regime 的存在与两期划分；波动带 K 稳健性；tail ∝K；消灭扫掠为小前缀动作。
- **工作猜想**：消灭事件等待 ∝K 的概率机制（约束点数随 K 增长）；连续极限解读。
- 无新文献声明；无 Real-world calibrated result。

## 9. What This Does NOT Prove

- 不证明消灭事件等待 ∝K 的解析规律（4 点观察 + 机制假说）；
- 不证明连续极限下 exact stopping 不可达（K→∞ 外推未定）；
- 不证明平稳 churn regime 的遍历性/混合性（无此类分析）；
- 不证明波动带 K 稳健性在其他 q/alpha 下保持；
- 不引入任何 controller 修改（perceptual B、hysteresis、noisy perception 等全部留待 Owner + GPT review）。
