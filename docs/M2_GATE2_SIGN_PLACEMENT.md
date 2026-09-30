# M2 Gate 2 文档：Sign-Placement Generality Gate——witness-aligned vs witness-opposed 符号翻转对照

> 状态：**Stage 3 / M2 Gate 2 execution 轮**（2026-09-20）。Gate 1（PASS）与 Interpretation Audit 之后的第二个 Stage 3 实验；Owner + GPT 冻结 prompt 执行。
>
> 性质：新增 2 个逐点符号反射的 sign-changing kernel（XA/XO）；canonical simulator 零修改；deterministic 一次运行可复现（`python experiments/m2_gate2_sign_placement/run_experiment.py`）。
>
> 实验报告（数值细节）：[experiments/m2_gate2_sign_placement/M2_GATE2_REPORT.md](../experiments/m2_gate2_sign_placement/M2_GATE2_REPORT.md)。

## 1. Research question（冻结）

Gate 1 已确立：positive mismatch existence alone 不足以产生 rebound（KW W0）；peak/integrated mismatch 不足以解释 canonical churn（KS S1）；audit 观察到 canonical rebound witness 集中在 interior 带 [0.20, 0.62]。Gate 2 只问：

> **保持同一个 sign-changing mismatch profile 的绝对值不变（Δ_O = −Δ_A 逐点），把 positive mismatch 放在历史 witness-support region（XA）与其外侧（XO），如何改变 rebound / recurrence / churn-like dynamics？**

主要因果对照 **XA ↔ XO**；KC 仅为 contextual reference（committed Gate 1 results，未重跑）。

## 2. Frozen kernel pair

| | XA（Witness-Aligned） | XO（Witness-Opposed） |
|---|---|---|
| Δ 符号 | >0 on (0, 2/3)；<0 on (2/3, 1) | 逐点反射：<0 on (0, 2/3)；>0 on (2/3, 1) |
| Δ 分段 | x/2 \| 1/3−x/2 \| 2/3−x \| x−1 | = −Δ_A |
| G | 3/2x \| x/2+1/3 \| 2/3 \| 2x−1 | x/2 \| 3/2x−1/3 \| 2x−2/3 \| 1 |
| density | 3/2 \| 1/2 \| 0 \| 2 | 1/2 \| 3/2 \| 2 \| 0 |
| max/min Δ | +1/6 @ x=1/3；−1/6 @ x=5/6 | 反射 |
| 断点 | 1/3、2/3、5/6（不取整到网格） | 相同 |

**Critical control property（已验证 112/112 + 测试）**：`Δ_O = −Δ_A` 逐点成立（dense + 全部 4 K 网格，max |Δ_O+Δ_A| = 1.1×10⁻¹⁶；G_O+G_A=2x）→ 同 breakpoints、同 regularity、同 piecewise complexity、同逐点 |Δ|、同 max amplitude、同初态、同控制器、同停止规则。**唯一操纵变量 = sign placement。**

共同配置：T=x、α=0.25、canonical q_near 初态（同 K 下逐位共享，runner 断言）、cumulative-excess selection、prefix removal、deterministic path、exact stop ≤1e-12、`q_j = G(j/K)−G((j−1)/K)` 解析 CDF 差分、τ_num=1e-12 仅验证用、H(K)=10K 或 exact stop（censoring 记录不延长）。

## 3. Execution matrix（冻结，全部执行）

XA @ K∈{50,100,200,400}；XO @ K∈{50,100,200,400}（固定矩阵，无顺序决策）。Per-round 诊断从第一轮就包含 audit 的 witness/risk 机制：`E_t(j)=M·Δ−γ`（repository `rebound_residual` helper）、raw + prefix-restricted argmax（tie 取小 j）、risk 处 Δ 符号分类、自身正支持区 in/out、R_t（M1C.4 语义）。gate：actual rebound witness 处 Δ 必须 positive（negative → halt；zero 仅在 numerical-boundary 轮容忍并记录）。

**静态验证 112/112；tests 全库 330 项通过；共享初态逐位一致；无任何 §29 correctness failure。**

## 4. Results

### 4.1 XA（Witness-Aligned）

全部 4 K：**t=6 exact stop、N_R=0**、max R_t 0.52–0.54、moved ≈0.76、D_max 单调降到机器零、numerical-boundary 0 例。risk（prefix-restricted）100% 在自身正支持区，x 中位 ≈0.32–0.33——即 XA 自己的 Δ 峰（1/3）附近，而非历史 witness 带中心（≈0.46）。

### 4.2 XO（Witness-Opposed）

全部 4 K：**10K horizon censoring**（Stage 3 首次 censoring，4/4）。

| K | N_R | f_R | S_R/K | max R_t | 类别 |
|---|---|---|---|---|---|
| 50 | 41 | 0.082 | 2.280 | 1.000 | Long-span recurrent |
| 100 | 581 | 0.581 | 9.930 | 1.021 | Long-span recurrent |
| 200 | 56 | 0.028 | 0.575 | 1.000 | Short recurrent |
| 400 | 2380 | 0.595 | 9.983 | 1.005 | Long-span recurrent |

动力学形态（finite-horizon 描述）：D_max 压到 **≈0.1344（K=50/100）/ ≈0.13778（K=200/400）的窄平台**并维持整个 horizon（K=400 p5–p95=[0.1378,0.1380]）；**边界只取 ~4 个离散值**（0.5、0.64、≈0.83、≈0.835），全部 rebound 在 a≈0.835（自身 Δ 峰 5/6 略右）；boundary switch freq ≈0；rebound 均为边际量（R max 1.00–1.02，inter-arrival 中位 2–3）；moved total 121.7–971.7（持续研磨）。与 canonical churn（边界跳转 sw_freq 0.92、D_max 带 0.003–0.08、f_R 0.36）是不同的动力学类型；不使用 attractor/周期轨道/non-convergence 词汇。

### 4.3 Theory verification

7524 active 轮：分解残差 max 5.33×10⁻¹⁵ ≤ τ_num；TP=2231、TN=4059；FP=407、FN=827 **全部 floating-boundary**（non-boundary FP=FN=0；max |worst|≈6.2×10⁻¹⁷；numerical-boundary 3492 轮全部在 XO）。

**Exact-tie riding（[Verified Numerical Analysis]）**：XO 平台动力学骑在判据 real-arithmetic 精确边际上——K=50 的 380 个 FP 轮中 377 个 actual d_D_max **逐位 = 0**（严格 `>` 为假；helper fp 求值 ±1.4×10⁻¹⁷ 摆动）。冻结规则处理：boundary 标记、语义不动、gate = non-boundary FP/FN=0。

### 4.4 Support migration（Gate 2 关键新诊断）

两个 kernel 的 criterion-relevant risk 都**迁移并锁定到自身正失配峰**（XA→1/3，XO→5/6；in-own-positive-support share：XA=1.000、XO=0.996–1.000）；全部 3058 个 actual rebound witness 处 Δ 严格 positive。**历史 canonical witness 带（audit [0.20,0.62]）对两个反事实 law 都没有预测力。**

## 5. Interpretation（frozen X1–X4）

**X4 成立（[Verified Simulation Result] / [Model-dependent Observation]）**：XO 比 XA 更复发——

> canonical historical witness distribution cannot be treated as a fixed intervention map for new kernels; witness geometry is endogenous to the law-induced trajectory.

附带（X1 倒置，不可 universal 化）：

> everywhere-positive mismatch is not necessary for long recurrent dynamics（正失配只占 1/3 定义域仍维持最长复发瞬态）；且"覆盖历史 witness 带 + 更大幅度的正失配"（XA）连一次 rebound 都没有。

X1/X2/X3 均与观察不符（观察到的恰好是 X4 的强证伪模式）。Gate 1 audit 的 static witness geometry 由本轮正式降级为 canonical 轨迹的 endogenous 产物。

## 6. Interpretation boundaries

- 全部结果只对：canonical initial state、α=0.25、K∈{50,100,200,400}、H=10K、deterministic path、这一对特定 piecewise-linear sign-changing profile 成立。
- XO 的"最终是否停止"在本轮 scope 内未知——只记录 pre-registered horizon censoring；禁止 non-convergent/attractor/chaos/ergodic/invariant-distribution 词汇；研究对象仍是 finite-time/pre-stopping transient dynamics。
- "边界锁定 ~4 值"是 horizon 内的描述性观察，不是周期轨道 claim。
- XA-vs-KC 不可用作 sign effect identification（profile/regularity/amplitude 全不同）；干净对照只有 XA-vs-XO。
- 类别 label 是 frozen reporting convention；每条附完整定量签名。

## 7. Unresolved questions（记录，不在本轮展开）

1. XO 的平台水平（≈0.1344 / ≈0.13778，近 K 无关）有没有解析来源？平台是否终结（何时）？
2. XO K=200 与 K=100/400 的类别跳变（Short vs Long-span）来源？
3. XA 零反弹的充分解释：正失配峰 1/3 处的 γ 结构是否系统性压制兑现（连接 Gate 1 ERC-1 的 Δ/γ 对齐问题）？
4. exact-tie 结构：哪些 piecewise-flat law/states 产生 real-arithmetic E≡0 边际骑乘？
5. gate return / episode 结构在 XO 平台动力学下是否存在（本轮未做 episode 审计）？

## 8. External Review Candidate Questions（本轮新增，登记 docs §11）

- **ERC-6**：什么 law-side 泛函（替代静态历史几何）预测哪些 sign placement 维持复发——正支持位置 × selection/export 几何的相互作用？
- **ERC-7**：XO exact-tie plateau 是否 admits 低维 reduced description（衔接 ERC-5）？piecewise-flat 结构产生精确边际 tie 的条件？

## 9. Explicit next step（冻结路线）

Gate 2 correctness **PASS** → 按 Stage 3 路线进入 **Stage 3 Closure（synthesis 文档）**，之后是 External Mathematical Review Gate。**不创建 Gate 3**；新问题只记录在 §7/§8。closure 文档本身等待 GPT + Owner review 后启动（本轮不做）。

## 10. 资产索引

| 交付 | 路径 |
|---|---|
| 冻结配置 | `experiments/m2_gate2_sign_placement/config.json` |
| Runner | `experiments/m2_gate2_sign_placement/run_experiment.py` |
| 实验报告 | `experiments/m2_gate2_sign_placement/M2_GATE2_REPORT.md` |
| 逐轮诊断（7524 行，含 witness/risk + sign + support 字段） | `results/per_round_diagnostics.csv` |
| 汇总表 / support migration / 静态验证 | `results/trajectory_summary.csv` / `support_migration.csv` / `static_kernel_validation.csv` |
| 图 A–D | `results/fig_gate2_*.png` |
| XA/XO kernel 构造器 | `src/sand_m0/model.py`（`witness_aligned_distribution`、`witness_opposed_distribution`） |
| 测试 | `tests/test_m2_gate2_sign_placement.py`（17 项；全库 330 项） |
