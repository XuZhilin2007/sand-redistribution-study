# M2 Gate 2 实验报告 — Sign-Placement Generality Experiment（XA vs XO）

> 日期：2026-09-20。运行方式：`python experiments/m2_gate2_sign_placement/run_experiment.py`（一次性重新生成全部 CSV / 图 / metadata）。
>
> 冻结配置：[config.json](config.json)。理论文档：[docs/M2_GATE2_SIGN_PLACEMENT.md](../../docs/M2_GATE2_SIGN_PLACEMENT.md)。
>
> 主要因果对照：**XA ↔ XO**（逐点符号反射对，`Δ_O = −Δ_A`，|Δ| 处处相同）。KC 仅为 contextual reference（committed Gate 1 results，未重跑）；XA-vs-KC 不是 sign effect identification。

## 1. 执行摘要

冻结 4-K 矩阵全部完成。**结果方向与 X4（强证伪分支）一致：把正失配放在历史 witness-support 区（XA）反而全程零反弹、6 轮 exact stop；把同一个 |Δ| profile 的符号完全翻转（XO，正失配只在 Far 侧）却在全部 4 个 K 上产生 Stage 3 至今最长的复发瞬态（10K horizon censoring）。**

| kernel | K=50 | K=100 | K=200 | K=400 |
|---|---|---|---|---|
| XA | t=6 exact stop，N_R=0 | t=6，N_R=0 | t=6，N_R=0 | t=6，N_R=0 |
| XO | **censored @10K**，N_R=41 | **censored @10K**，N_R=581 | **censored @10K**，N_R=56 | **censored @10K**，N_R=2380 |

理论验证：7524 active 轮，分解残差 max 5.33×10⁻¹⁵ ≤ τ_num；判据 TP=2231/TN=4059，**FP=407/FN=827 全部为 floating-boundary 情形（non-boundary FP=FN=0，max |worst| ≈6×10⁻¹⁷）**——XO 动力学骑在判据的精确边际 tie 上（详见 §6），严格不等式语义未动。

## 2. Gate 流程

| 阶段 | 内容 | 结果 |
|---|---|---|
| A2 bookkeeping | Gate 1 audit 的 KC Far share 核对：真实口径 **6/1958 = 0.31%**（per-K 0.19–0.85%）；audit 文档中"31/1958"为 prose 笔误 → 已最小更正（raw data 未动） | 完成 |
| A4 静态验证 | XA/XO 各 4 K：CDF 端点/单调/连续性（1/3、2/3、5/6 断点处分支公式求值）、归一化、密度非负、±1/6 极值（解析点 + 格点有界）、**符号反射 \|Δ_O+Δ_A\|≤τ 与 G_O+G_A=2x**、sign supports | **112/112 通过** |
| A5 tests | 新套件 17 项；全库 **330 项运行时断言通过** | 通过 |
| B 执行 | XA/XO × K∈{50,100,200,400} 固定矩阵（无中途增删）；共享 canonical initial state（逐位断言）；H(K)=10K | 完成 |
| C 验证/分析 | 理论验证、witness/risk + sign + support migration、表格、图 | 完成 |

## 3. XA Results（Witness-Aligned：正失配覆盖 (0,2/3)，含 canonical witness 带）

| K | active | t_stop | N_R | f_R | max R_t | moved total |
|---|---|---|---|---|---|---|
| 50 | 6 | 6 | 0 | 0.000 | 0.541 | 0.762 |
| 100 | 6 | 6 | 0 | 0.000 | 0.521 | 0.759 |
| 200 | 6 | 6 | 0 | 0.000 | 0.522 | 0.760 |
| 400 | 6 | 6 | 0 | 0.000 | 0.521 | 0.761 |

D_max 严格单调下降到机器零；numerical-boundary 0 例。witness/risk（prefix-restricted）：**100% 落在自身正支持区**，x 中位 ≈0.32–0.33（全部轮）/ ≈0.36（top-25% R_t 轮）——即 XA 自己的 Δ 峰（x=1/3）附近，而不是历史 canonical witness 带的中心（≈0.46）。

## 4. XO Results（Witness-Opposed：正失配只在 (2/3,1)）

| K | active | stop | N_R | f_R | S_R/K | max R_t | moved total | 类别 |
|---|---|---|---|---|---|---|---|---|
| 50 | 500 | **censored @10K** | 41 | 0.082 | 2.280 | 1.000 | 121.7 | Long-span recurrent |
| 100 | 1000 | **censored @10K** | 581 | 0.581 | 9.930 | 1.021 | 242.6 | Long-span recurrent |
| 200 | 2000 | **censored @10K** | 56 | 0.028 | 0.575 | 1.000 | 486.3 | Short recurrent |
| 400 | 4000 | **censored @10K** | 2380 | 0.595 | 9.983 | 1.005 | 971.7 | Long-span recurrent |

**XO 动力学形态（全部为 finite-horizon 描述，非 attractor/周期轨道 claim）**：

- D_max 快速压到 **≈0.1344（K=50/100）/ ≈0.13778（K=200/400）的窄平台**并在整个 10K horizon 上维持（K=400 的 p5–p95 = [0.1378, 0.1380]）；
- 边界只取 **~4 个离散值**（0.5、0.64、≈0.83、≈0.835），全部 rebound 轮 a≈0.835（即 XO 自身 Δ 峰 5/6≈0.833 略右）；boundary switch freq ≈ 0（|Δa| 中位 0–0.01）——与 canonical churn 的大幅边界跳转（sw_freq 0.92）完全不同；
- rebound 均为**边际量**（max R_t 1.00–1.02，inter-arrival 中位 2–3 轮）；K=50 出现 382 轮的最长无反弹连段（f_R 仅 0.082）；
- risk/witness（prefix-restricted）：**≈100% 钉在自身正支持区**（x_med ≈0.83–0.84 = 5/6 峰），sign 全部 positive（无 warning、无 gate 触发）。

这是 Stage 3 首次出现 **finite-horizon censoring**（4/4）。按冻结规则记录为 `not stopped within pre-registered Gate-2 horizon`，不延长、不作 non-convergence 解读。

## 5. XA vs XO 直接对照（本轮核心）

| 维度 | XA | XO |
|---|---|---|
| rebound / recurrence | 4 K 全部零反弹 | 4 K 全部复发（K=100/400 f_R≈0.58–0.60） |
| stop | 4 K t=6 exact stop | 4 K 10K censoring |
| S_R/K | — | 2.28 / 9.93 / 0.58 / 9.98 |
| 边界行为 | 正常移动（sw_freq 0.4） | 锁定 ~4 值、rebound 全在自身 Δ 峰处（sw_freq 0） |
| D_max 形态 | 单调降到机器零 | 窄平台 ≈0.134–0.138 维持整个 horizon |
| risk 位置 | 自身正峰 1/3（100% in-support） | 自身正峰 5/6（≈100% in-support） |
| moved total（K=100） | 0.76 | 242.6（>300×） |

唯一被操纵的变量是 sign placement（|Δ|、breakpoints、regularity、初态、控制器全部相同）。**符号放置完全决定了动力学命运。**

## 6. One-step theory verification

- 分解恒等式：7524 active 轮，max residual **5.33×10⁻¹⁵** ≤ τ_num。
- 判据：TP=2231、TN=4059；FP=407、FN=827 **全部 floating-boundary**（|worst| ≤ 1e-12；非边界 FP=FN=0；FP/FN 轮 max |worst| ≈ 6.2×10⁻¹⁷）。numerical-boundary 轮 3492（全部在 XO）。
- **精确边际结构（Verified Numerical Analysis）**：XO 平台动力学骑在判据的 real-arithmetic 精确 tie 上——K=50 的 380 个 FP 轮中 377 个 actual `d_D_max` **恰好 = 0.0**（状态更新把 D_max 逐位复现，严格 `>` 为假；helper 的 fp 求值在 ±1.4×10⁻¹⁷ 间摆动）。按冻结规则：boundary 标记、语义不动、以 non-boundary FP/FN=0 作为 correctness gate。
- 质量守恒 ≤ τ_num；共享初态逐位一致；witness sign gate：全部 3058 个 actual rebound witness 的 prefix-restricted argmax 处 Δ 严格 **positive**（0 个 zero、0 个 negative；negative 符号仅出现在 8 个非 rebound 轮的 risk 位置上）。

## 7. Support migration（Table 3 要点）

- **两个 kernel 的 criterion-relevant risk 都迁移并锁定到自身正失配峰**：XA → x≈1/3，XO → x≈5/6；in-own-positive-support share：XA=1.000（全部样本），XO=0.996–1.000（raw view 的 0.2–0.4% 负号轮全部是 suffix-trivial raw argmax；prefix view 的 zero/negative 仅在 boundary 轮）。
- **历史 canonical witness 带（audit 的 [0.20, 0.62] interior band）对两个反事实 law 都没有保留预测力**：XA 的正失配覆盖该带却没有产生任何 rebound；XO 的正失配完全在该带之外却产生最长复发瞬态。

## 8. Interpretation（frozen X1–X4 matrix）

| 分支 | 预测 | 观察 | 判定 |
|---|---|---|---|
| X1 | XA recurrent、XO 不复发 | **完全相反** | 不适用 |
| X2 | 两者都快速停止 | XA 停、XO censored | 不适用 |
| X3 | 两者类似复发 | 最大程度不同 | 不适用 |
| **X4** | **XO 更复发 → 强证伪** | **观察到** | **成立** |

**X4 支持（[Verified Simulation Result] / [Model-dependent Observation]）**：

> canonical historical witness distribution cannot be treated as a fixed intervention map for new kernels; witness geometry is endogenous to the law-induced trajectory.

附带支持（X1 的倒置形式，同样不可写成 universal claim）：

> everywhere-positive mismatch is not necessary for long recurrent dynamics——符号翻转后正失配只占 1/3 定义域（Far 侧）仍维持 Stage 3 至今最长的复发瞬态；而 100% 正覆盖（XA，且幅度 1/6 > KW 的 1/16）连一次 rebound 都没有。

**Gate 1 audit 的 static witness geometry 由此正式降级**：它只是 canonical 轨迹的 endogenous 产物，不是可迁移的干预地图。

## 9. Unexpected / Falsifying Results（完整保留）

1. **方向反转本身**：直觉（把正失配放 witness 带 → 更容易 rebound）被数据否定；真实行为由 law 自身诱发的 risk geometry 决定。
2. **XO 的平台锁定形态**是 Stage 3 新动力学类型（窄 D_max 平台 + ~4 值边界 + 边际 rebound + 零 boundary switching），与 canonical churn（hopping、宽 D_max 带）和 B/U（单调速停）都不同。
3. **exact-tie riding**：XO 大量轮次骑在判据 real-arithmetic 精确边界上（numerical-boundary 占其 active 轮的 20–99%），FP/FN 全部可由 floating boundary 解释。
4. XO K=200 与 K=100/400 的类别分辨率敏感性（Short vs Long-span；f_R 0.028 vs 0.58–0.60）——XO 内部也不是单一形态。

## 10. External Review Candidate Questions（本轮新增）

- **ERC-6**：什么 law-side 泛函（替代静态历史 witness 几何）能预测哪些 sign placement 维持复发？XO 的证据指向"正支持区位置 × selection/export 几何（平台边界处 1−G_O(a)）"的相互作用，而非 historical alignment。
- **ERC-7**：XO 的 exact-tie plateau（real-arithmetic E≡0 的边际平衡、边界锁定在 Δ 峰附近 ~4 个值）是否 admits 低维 reduced description（与 ERC-5 的 operator-level 问题衔接）？piecewise-flat G 结构在何种条件下产生精确边际 tie？

## 11. 完成度对照（§30 stop condition）

XA @ 4 K ✓、XO @ 4 K ✓、静态验证 ✓、理论验证 ✓、witness/risk 分析 ✓、support migration ✓、汇总指标 ✓、最小图（A/B/C/D）✓、tests（330 项官方口径）✓、canonical Gate 2 report（本文 + docs）✓、research log ✓、next-step update ✓ → **Gate 2 STOP 条件满足，PASS**。未加第三个 kernel、未调 breakpoint/amplitude/sign boundary、未开始 Gate 3 或 Stage 3 closure 文档。
