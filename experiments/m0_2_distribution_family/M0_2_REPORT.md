# M0.2 Phase 2A 实验报告：Throw Distribution Family

> 日期：2026-09-18
>
> 理论文档：[docs/M0_GENERAL_Q.md](../../docs/M0_GENERAL_Q.md)；模型：[docs/M0_MODEL.md](../../docs/M0_MODEL.md)；代码：`src/sand_m0/`；配置：[config.json](config.json)（运行前冻结）
>
> 证据等级：解析结论为 **模型推导数学结果（Model-derived Mathematical Result）**，数值对照为 **已验证模拟结果**，跨分布的现象比较为 **模型依赖观察**。均只对 fixed-q / fixed-zone 的 M0 类模型成立。

## 1. 唯一核心问题

> M0/M0.1 中发现的 finite-time optimum 与 time-rescaling，是 `q(x)=2(1-x)` 的特殊巧合，还是一类不同 throw distributions 中都会出现的结构？

**回答：是一整类模型的结构性质，不是 q_near 的巧合。** 依据：§2 的一般推导（先于模拟完成并独立核验）+ §3 的六个分布对照（预测在运行前写出，运行后 6/6 命中）。

## 2. 理论先行（运行模拟之前完成）

### 2.1 一般 q 的结构（Model-derived Mathematical Result）

对**任意**固定 q（初始分布 = 重泼分布 = q）、固定近区、deterministic 期望质量动力学：

1. **形状锁定**：近区 bin = `z·q_i`，远区 bin = `λ·q_i`，`λ=(1−Qz)/Q_f`——两区永远保持 q 的形状，只改幅度；
2. **一维轨迹**：整个系统由单一变量 `z = ρ^t`（`ρ = 1−alpha(1−Q)`）描述；
3. **精确二次函数**：`U(z) = S_n z² + S_f((1−Qz)/Q_f)² − 1/K`；
4. **time-rescaling**：路径、U(z)、z* 与 alpha 无关，alpha 只决定 z 的前进速度（成立条件：固定 q、固定区、线性期望动力学、`0<Q<1`、`alpha∈(0,1]`）；
5. **interior optimum 一般条件**：存在当且仅当 **`Q_f·S_n > Q·S_f`（即每单位质量集中度 `c_n = S_n/Q > c_f = S_f/Q_f`）**；`c_n ≤ c_f` 时 t=0 最佳、单调恶化。

以上全部由 M0 更新规则严格推导，并独立核验：一个**不属于任何分布族**的非对称分段线性 q 与模拟器一致至 1.4×10⁻¹⁷；六个族内成员一致至 ≤1.1×10⁻¹⁷。完整推导见 [M0_GENERAL_Q.md](../../docs/M0_GENERAL_Q.md) 数学附录。

### 2.2 冻结的分布族（先于任何模拟写入 config.json）

**幂族** `q_p(x) = (p+1)(1−x)^p`，p > −1：p 是"朝向操作者的幂律倾斜"，p=0 均匀、p=1 恰为 M0 baseline q_near、p>0 近偏、p<0 远偏；CDF 与反函数解析，采样与 bin 概率稳定。闭式（连续、a=0.5）：`c_n/c_f = (2^(2p+1)−1)/(2^(p+1)−1)`，`z*_cont = 2(1−u)/(3−4u)`，`u=2^{−(p+1)}`，故 **interior ⟺ p > 0**。历史 `r=U²` 按计划保持为未来单独 sensitivity case，未被强行纳入。

六个 case（含均匀对照、远偏反例、近边界 mild case）在运行前写入预测表 `predictions.csv`：

| p | 角色 | Q | c_n/c_f | z*_disc | 理论 t*（cont） | 预测最佳轮 | 预测最佳 U |
|---:|---|---:|---:|---:|---:|---:|---:|
| −0.25 | 远偏反例 | 0.4054 | 0.617 | 1.29451 | —（无 interior） | 0 | 1.1232e-03 |
| 0 | 均匀边界 | 0.5000 | 1.000 | 1.00000 | —（z*=1 边界） | 0 | ≈0 |
| 0.25 | 近边界 mild | 0.5796 | 1.327 | 0.87924 | 1.16 | 1 | 2.1924e-04 |
| 0.5 | 温和近偏 | 0.6464 | 1.641 | 0.81526 | 2.21 | 2 | 6.4528e-04 |
| 1.0 | baseline | 0.7500 | 2.334 | 0.74997 | 4.46 | 4 | 1.6796e-03 |
| 2.0 | 强近偏 | 0.8750 | 4.429 | 0.69995 | 11.23 | 11 | 3.9496e-03 |

预测要点：**最佳轮次随近偏强度单调后移**（1→2→4→11），因为近偏越强、重泼落回近区的份额越大（Q 越大 → ρ=1−alpha·Q_f 越接近 1），z 前进越慢；同时**最佳状态本身也越差**（U(z*) 随 p 上升）。

## 3. 模拟对照（canonical simulator，baseline alpha=0.25）

| p | interior 预测/观察 | 最佳轮 预测/观察 | 最佳 U 预测/观察 | z* 预测 / 观察最佳轮处 z | U(0)→U(30) | max\|U_sim−U(z)\| |
|---:|---|---|---|---|---|---:|
| −0.25 | 无/无 | 0/0 | 1.1232e-03 / 1.1232e-03 | 1.29451 / 1.00000 | 1.12e-3 → 1.20e-2 单调恶化 | 3.5e-18 |
| 0 | 无/无 | 0/0 | ≈0 / ≈0 | 1.00000 / 1.00000 | 0 → 9.64e-3 单调恶化 | 9.5e-18 |
| 0.25 | 有/有 | 1/1 | 2.1924e-04 / 2.1924e-04 | 0.87924 / 0.89489 | 4.16e-4 → 9.98e-3 | 3.7e-18 |
| 0.5 | 有/有 | 2/2 | 6.4528e-04 / 6.4528e-04 | 0.81526 / 0.83104 | 1.25e-3 → 1.08e-2 | 6.9e-18 |
| 1.0 | 有/有 | 4/4 | 1.6796e-03 / 1.6796e-03 | 0.74997 / 0.77248 | 3.33e-3 → 1.14e-2 | 1.0e-17 |
| 2.0 | 有/有 | 11/11 | 3.9496e-03 / 3.9496e-03 | 0.69995 / 0.70523 | 8.00e-3 → 8.39e-3 | 5.2e-18 |

- **最佳轮次 6/6 精确命中**；interior 存在性 6/6 正确；最佳 U 6/6 与理论一致（打印精度内）。
- "观察最佳轮处 z" 与 z* 的差（如 p=1: 0.77248 vs 0.74997）是**轮次离散化**：观察值是离 z* 最近的整数轮位置，不是 z* 本身——这与 M0.1 的浅平台解释一致，不是预测失败。
- 六条 U 轨迹与各自精确二次函数的最大偏差 ≤1.0×10⁻¹⁷：**整个族都在一维路径上精确演化**。
- 轨迹形态：三个 interior case 都是"先改善→有限轮最佳→恶化"；p=2 的终态（8.39e-3）只比初始（8.00e-3）差 4.9%——强近偏时终态几乎回到起点（理论条件 `S_n/S_f < 1/Q_f²−1` 预测该行为，31 < 63 ✓）。

## 4. Alpha time-rescaling spot checks（非参数优化）

p ∈ {0.5, 2.0, −0.25} × alpha ∈ {0.10, 0.50, 1.00}（`alpha_rescaling.csv`）：

| p | z*（三档 alpha 下恒定） | 理论 t* → 观察最佳轮 |
|---:|---:|---|
| 0.5 | 0.815260（恒定） | 5.67→**6**；1.05→**1**；0.47→**1** |
| 2.0 | 0.699953（恒定） | 28.36→**28**；5.53→**6**；2.67→**3** |
| −0.25 | 1.294513（>1，恒定） | 无 interior：三个 alpha 下 t=0 恒最佳 |

- **z* 完全不随 alpha 改变**；改变 alpha 移动的是最佳轮次，且理论 t* 9/9 命中观察（p=2 在 alpha=0.10 下理论预测第 28 轮、模拟正好第 28 轮）。
- 全部 9 个 (p, alpha) 组合的 (z, U) 点落在各自 q 的精确二次曲线上（残差 ≤1.1×10⁻¹⁷）。
- 反例侧：p=−0.25 的 c_n<c_f 由 q 唯一决定，**任何 alpha 都不能制造 interior optimum**。

## 5. Monte Carlo spot checks

p ∈ {0.5, 2.0, −0.25}，alpha=0.25，各 5 种子（20260917–21，N=10⁵）。逐轮 MC 均值对"U(z)+噪声底线"的最大偏差：1.80×10⁻⁵（p=0.5）、4.06×10⁻⁵（p=2.0）、1.80×10⁻⁵（p=−0.25）——与 5 种子均值的抽样散布量级一致（单种子 sd ≈ 3×10⁻⁵），无随 p 或 t 的系统漂移。**有限粒子版本围绕 deterministic 预测波动。**

## 6. Counterexamples / 边界 case（主动保留）

1. **p=−0.25（远偏反例）**：`c_n/c_f = 0.617 < 1` → 无 interior optimum，31 轮单调恶化，且**任何 alpha 都无法改变**（§4）。这是"finite-time optimum 不是普适现象"的族内实证。
2. **p=0（均匀边界）**：`c_n = c_f` 精确成立，z* = 1 恰好落在路径起点——理论上的临界情形，观察为 t=0 最佳、单调恶化（`U(t) = 0.01(1−0.875^t)²` 精确成立，偏差 6.9×10⁻¹⁸）。
3. **p=0.25（近边界 mild）**：interior 存在但很浅（最佳 U 距初始仅降 47%后迅速恶化，margin to runner-up 8.1e-5）——靠近边界的 interior case，表明"p 刚过 0"时现象存在但微弱。
4. **连续闭式的适用边界**：连续积分 ∫q² 在 p ≤ −0.5 发散，故连续闭式要求 p > −0.5；离散实现只需 p > −1。本轮远偏 case 取 p=−0.25 保证两者都可计算；有限-bin 修正最大处也在该 case（c 比值 +1.6%），已如实记录（[M0_GENERAL_Q.md](../../docs/M0_GENERAL_Q.md) §A.6）。

## 7. 对"q_near 是否特殊"的明确回答

`q_near`（p=1）**不是**一个特殊的漂亮例外，而是一般理论中的一个普通点：

- 形状锁定、单变量轨迹、精确二次 U(z)、time-rescaling 对**任意固定 q** 成立（在 M0 的规则假设下）——M0.1 的结果从"单个分布的性质"升级为 **fixed-q / fixed-zone M0 class 的结构性质**；
- interior optimum 的存在性由不等式 `c_n > c_f` 完全决定，q_near（比值 2.33）只是舒适地落在 interior regime 内部的一个成员；远偏与均匀分布落在另一侧，边界在 `c_n = c_f`；
- q_near 的具体数值（z*=3/4、t*≈4.46）依赖它的 7:1 集中度比——族内其他成员给出不同 z* 与最佳轮，但**都由同一公式预测**。

## 8. 证据分类

**模型推导数学结果（Model-derived Mathematical Result）**：§2.1 的五条一般结构（形状锁定、一维轨迹、精确二次式、time-rescaling 及其条件、interior 条件 `Q_f·S_n > Q·S_f`）与幂族闭式——严格推导 + 独立数值核验（含族外分布），只对 fixed-q / fixed-zone M0 类模型成立。

**已验证模拟结果**：六个分布的 deterministic 轨迹与二次式一致至 ≤1.1×10⁻¹⁷；最佳轮 6/6、interior 存在性 6/6、alpha spot check 9/9 命中理论预测；MC spot check 无系统偏离。全部由一条命令复现（commit/config/seeds 见 metadata.json）。

**模型依赖观察**：最佳轮次随近偏强度单调后移（1→2→4→11）且最佳状态随近偏变差；interior case 的浅平台；强近偏终态几乎回到起点。这些依赖本冻结族、`a=0.5`、K=100 与 histogram 指标。

**工作猜想**：`c_n > c_f` 判据可能在任何"固定选择几何 + 固定重泼形状"的选择-重采样过程中成立；真实沙坑的单次泼沙若落在 interior regime，将呈现同类有限最优——两条都待 Phase 3（状态依赖策略预计破坏 time-rescaling）与未来现实粗校准检验。

**Preliminary literature positioning**（仅记录，未核实原文，不构成 Literature-supported Result）：rank-driven Markov processes（代表工作侧重 long-time / large-N threshold 与 limiting distribution）、repeated balls-into-bins（反复移除-重分配机制，侧重 equilibrium / mixing / self-stabilization）、balanced allocations（侧重 load gap / maximum load）——均与本项目机制相邻但研究问题不同。本轮未做原文核实，故不升级为文献支持结果。

## 9. 本文不证明的事情

- 不证明所有 redistribution model 都有 interior optimum 或 time-rescaling（远偏反例与边界 case 恰恰画出了"没有"的区域）；
- 不证明真实沙坑的单次泼沙分布属于 interior regime——真实 throw distribution 尚未测量；
- 不证明 phase diagram（6 个 case 不足以划分行为区域，本轮仅做 ordered comparison / condition map）；
- 不证明 novelty；无 Literature-supported Result；`r=U²` 未纳入本轮。

## 10. 复现信息

- **Python** 3.14.6；**依赖** numpy 2.5.3、matplotlib 3.11.2（仓库专用 `.venv`）。
- **命令**：`.venv/Scripts/python.exe experiments/m0_2_distribution_family/run_experiment.py`（runner 内部先写 predictions.csv 再运行模拟，次序即认识论次序）。
- **测试**：`tests/test_m0_2_family.py`（53 项断言；另 M0 37 项、M0.1 25 项保持通过）。
- **种子**：MC spot check 用 20260917–20260921；deterministic 无随机性。
- **结果文件**（`experiments/m0_2_distribution_family/results/`，SHA-256 见 metadata.json）：`predictions.csv`（运行前理论表）、`comparison.csv`（预测 vs 观察）、`deterministic_trajectories.csv`（6×31 行）、`alpha_rescaling.csv`、`mc_spotcheck.csv`、四张 PNG、`metadata.json`。
- **冻结声明**：family 定义、case 列表、alpha 集合、seeds 全部在 config.json 中于任何 M0.2 模拟运行之前写定；本轮无事后调参。

## 11. 下一步建议

1. **Phase 2B（若 Owner 同意）**：把判据推广到非幂族分布（Beta family、历史 `r=U²` 作为独立 case），检验 `c_n > c_f` 条件在新族上的预测力；
2. **Phase 3**：引入状态依赖回收策略，观察 time-rescaling 如何被破坏（预测：路径不再 alpha-无关）；
3. 在任何 phase-diagram 工作之前，先积累足够多的确认 case。
