# M1B 模型定义：感知容差 / Deadband 停止

> 状态：M1B canonical 定义（2026-09-18）。相对 M1A **只改一个机制**：停止判据从 `D_max ≤ 10⁻¹²` 改为 `D_max ≤ epsilon`（有限容忍度）。边界规则、更新规则、q、alpha、指标、初始条件全部不变。
>
> 实验：[M1B_REPORT.md](../experiments/m1b_tolerance/M1B_REPORT.md)；前置：[M1A_MODEL.md](M1A_MODEL.md)、[M1A_LONG_HORIZON.md](M1A_LONG_HORIZON.md)
>
> 证据等级声明：`D_max ≤ epsilon` 停止是 **模型假设（Model Assumption）**——"没有明显值得处理的区域就停止"（真实观察）的最小数学代理。epsilon **不是**人眼阈值（JND），未做任何心理物理校准。

---

# 第 1 层：直觉解释（不需要高数）

## 1.1 为什么"精确零"不现实

M1A 的停止条件是"任何一段都不比均匀**多出哪怕一丝**"（D_max ≤ 10⁻¹²）。M1A.1 发现这可以实现，但要 324–2876 轮（K=50–400），而且**网格越细、耗时越长**（约 ∝K）——就像要求地面绝对一粒沙不差：分辨率越高，越难宣布"完全没差"。现实中操作者的标准从来不是零，而是"**没有明显值得处理的多余沙**"。M1B 把这个"明显"写成一个小正数 epsilon。

## 1.2 epsilon 是什么

epsilon = 0.01 的含义：**从脚边往前任何一段，累计沙量比完全均匀时多出的部分，不超过总沙量的 1%**。就这一句话。它作用在确定性 D_max 信号上，不是视觉模型，也没有噪声；"人眼能否分辨 1%"是完全另一个（未做）的问题。

## 1.3 一个重要发现：容差解决了网格敏感，但停止会偏早

冻结实验（epsilon ∈ {0, 0.5%, 1%, 2%, 4%} × K ∈ {50,100,200,400}）的结果：

- **epsilon ≥ 0.5% 时，停止轮数几乎与 K 无关**（K=50→400 的停止轮变化 ≤1.34 倍；epsilon≥1% 时 ≤1.09 倍），而精确零标准的停止轮在同样范围内差 8.9 倍。**网格敏感被有限容差消除了。**
- **但所有有限 epsilon 都在均匀度最优点之前停**：精确零标准在 U 最优之后 ~46 轮才停（轻微过修正），而 epsilon=0.5%/1% 在 U 最优之前 313/494 轮（K=100）就停了。原因很结构化：**U 的改善一直持续到 ∝K 的时刻，而容差停止的时刻与 K 无关**——两者只在个别点相遇。
- 于是"停止态离 U 最优的差距"（regret）随 epsilon 单调变大：epsilon=0.5% 时 0.8–1.65 倍，1% 时 2.0–2.5 倍，4% 时 7.3–7.9 倍。

直觉总结：**D_max（最重前缀的过剩）衰减得快，U（整体均匀度）改善得慢**——两个指标的节奏不同。容差截止的是 D_max 的尾声，但 U 还在慢速爬坡。想要"又快又网格无关"要付出"停止时离最均匀还差一截"的代价；这就是 M1B 的核心取舍，不是实现缺陷。

## 1.4 与控制论 deadband 的关系

控制/质量控制里早有同样想法：偏差在动作限内就不调整（deadband / action limit），避免在阈值附近频繁抖动。M1B 的"D_max ≤ epsilon 就不再扫"与它是同一个结构。这是建模传统的联系，不是校准。

**文献核实更新（2026-09-18，M1B.1）**：引用已通过 Crossref API 独立核实——Alberto Luceño, "Ch. 19. Dead-band adjustment schemes for on-line feedback quality control", *Handbook of Statistics* Vol. 22, Elsevier, 2003, pp. 695–727, DOI `10.1016/S0169-7161(03)22021-6`。provenance：书目元数据经 Crossref 验证；dead-band 结构描述（周期测量偏差、仅在超出 action limits 时调整）来自章节标题与 GPT 阅读建议，本轮未读全文。记录为 **Literature-supported modeling connection**：不支持 human perception、不支持任何 epsilon 数值、不构成 novelty。

---

# 第 2 层：Mathematical Definition

## 2.1 记号（承 M1A）

状态 `m(t) ∈ R^K`，累计过剩 `D_t(j) = Σ_{i≤j} m_i(t) − j/K`，`D_max(t) = max_j D_t(j)`。由阶梯 CDF 的右连续性，`D_max(t) = sup_{x∈[0,1]} [F_t(x) − x]` **精确成立**（one-sided CDF discrepancy / one-sided KS functional，证明与核验见 [M1A_LONG_HORIZON.md](M1A_LONG_HORIZON.md)）。

## 2.2 epsilon-deadband 停止

```text
M1B: 于状态 t，若 D_max(t) ≤ epsilon：停止（状态冻结，实验结束）；
     否则按 M1A 规则扫掠（边界 = argmax D_t(j)，前缀比例 alpha 移除，q 重泼）。
```

epsilon = 0 的 case 用 canonical 数值容差 10⁻¹² 实现，即精确 M1A 基线。

## 2.3 First-passage 解释与单调性（Model-derived Mathematical Result）

**命题 1（first-passage 恒等）**：记 reference 轨迹为 epsilon=0 的 M1A 轨迹（在停止前不被截断的完整演化），`D*_max(t)` 为其过剩序列。则 M1B(epsilon) 的停止时刻

```text
t_eps = min{ t ≥ 0 : D*_max(t) ≤ epsilon }
```

且 M1B(epsilon) 在 t ≤ t_eps 的每一步都与 reference 逐位相同（归纳：停止前更新规则相同）。因此停止态 = reference 在 t_eps 的状态——**不需要单独模拟，reference 序列预测一切 epsilon 的停止**（20/20 个 case 数值验证）。

**命题 2（epsilon 单调性）**：`epsilon₁ ≤ epsilon₂ ⟹ t_{eps₂} ≤ t_{eps₁}`。证明：`{t : D*≤eps₁} ⊆ {t : D*≤eps₂}`，首达时间对扩张的集合不增 ∎（200 点网格 × 4 K 数值验证）。

推论：`D_max(t) ≤ epsilon` 的停止就是 reference 信号对水平 epsilon 的**首达时间（first-passage time）**——该定义是良构的，且全部 M1B 量（停止态、regret）都是 reference 轨迹的读数。

## 2.4 分辨率极限（Working Hypothesis）

M1A.1 观察到 t_stop(exact) ∝ K。M1B 数值显示：对固定 epsilon ≥ 0.005，`t_eps` 在 K=50–400 间变化 ≤1.34 倍（epsilon ≥ 0.01 时 ≤1.09 倍）——**数值迹象表明 t_eps(epsilon) 对固定正 epsilon 收敛到 K 无关的有限值**。机制：D_max(t) 作为 CDF 差距的连续泛函，其离散近似误差为 O(1/K)；当 epsilon ≫ O(1/K) 时首达发生在离散/连续信号一致的阶段。**未证明定理**；且 epsilon 与 1/K 可比时（如 epsilon=0.005, K=400：0.005 vs ~0.0025）K 依赖部分回归（观测到 K=50 的 253 vs K≥100 的 ~330）。

## 2.5 冻结的 epsilon 集与禁止事项

`epsilon ∈ {0（=10⁻¹²）, 0.005, 0.01, 0.02, 0.04}` 运行前冻结；是简单的小到大尺度，非校准。禁止事后追加"恰好停在 U 最优"的 epsilon 作为正式结果；epsilon 不解释为 JND；hysteresis（双阈值/重启记忆）不在本轮（当前只有 continue→stop 单向转移，无重启机制）。
