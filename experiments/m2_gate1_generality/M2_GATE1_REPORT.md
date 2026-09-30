# M2 Gate 1 实验报告 — Redistribution-Mismatch Generality Experiment

> 日期：2026-09-20。运行方式：`python experiments/m2_gate1_generality/run_experiment.py`（一次性重新生成全部 CSV / 图 / metadata）。
>
> 冻结配置：[config.json](config.json)。理论文档：[docs/M2_GATE1_GENERALITY_GATE.md](../../docs/M2_GATE1_GENERALITY_GATE.md)。
>
> 本文只报告可复现结果；解释边界与证据分级见 docs 文档。治理原则：repository reality > prompt > 聊天历史；所有 gate 未通过即停止科学解释（本轮全部 gate 通过，无中止）。

## 1. 执行摘要

Gate 1 冻结执行矩阵全部完成：**17 条轨迹**（5 kernels × 各自冻结 K 集合），**0 条触发 finite-horizon censoring**（全部远早于 H(K)=10K exact stop）。理论验证 **5510 个 active 轮 TP=1983 / TN=3527 / FP=0 / FN=0，numerical-boundary 0 例**，最大分解残差 7.55×10⁻¹⁵ ≤ τ_num=1e-12。

核心结果（K=100 主筛查，cross-K 见 §5）：

| kernel | 结果 | 解读分支（frozen contract） |
|---|---|---|
| K0（G=T） | t=6 exact stop，无反弹（全部 K） | integrity reference，复现 M1C.1 U |
| K-（G=x²） | t=3 exact stop，无反弹，严格单调 | 负控制全部通过（定理+模拟一致） |
| KW（λ=1/4） | t=8 exact stop，**零反弹**（全部 K） | **W0** |
| KC（canonical） | t=699，churn transient（reproduction） | integrity anchor，复现 Stage 2 |
| KS（flat-top） | t=26，**5 次反弹但短命**（K=400：10 次，t=29） | **S1** |

一句话：**单步 rebound/churn 判据完美泛化（零误分）；但 canonical 的 long-span churn-like transient 在两个 positive 对照上都没有出现——KW 连单次反弹都没有，KS 只有短促复发后快速 exact stop。**

## 2. Gate 流程与验收

| Gate | 内容 | 结果 |
|---|---|---|
| A2 静态验证 | 新 kernels（K-/KW/KS）在其全部执行 K 上：G(0)=0、G(1)=1、单调性、`q_j = G(j/K)−G((j−1)/K)` 精确 CDF 差分归一、mismatch 公式、KS 断点连续性 | **105/105 required 通过**（含 K0/KC 参照行 173/173） |
| A3 K0 integrity | 与已提交 M1C.1 Variant U rows 字符串比较（16 字段 × 4 K） | 通过（见合计） |
| A4 KC integrity | 与已提交 M1C.1 Variant A rows 比较 + t_stop 对 M1B.2 annihilation events | **合计 132/132 通过** |
| A5 K- 验收 | 零 genuine rebound、active 轮 D_max 严格下降（最小降幅 −0.130）、分解残差 ≤ τ、判据 FP=FN=0、质量守恒 ≤ τ | **5/5 通过** |

## 3. K- 负控制 @ K=100（定理后盾）

G_far(x)=x² ≤ x=T 触发 M1C.3 序定理 ⟹ 每个 active 轮 D_max 严格下降、永不反弹。实际：**t=3 exact stop**（比 K0 的 t=6 更快——远偏 law 的 net export `1−G(a)=1−a² > 1−a` 每轮带走更多前缀外质量），D_max 0.25 → ~2×10⁻¹⁶，total moved 0.579（KC 同 K 为 103.3）。定理与模拟完全一致；无任何需要解释的科学反弹。

## 4. KW 结果（amplitude ablation，λ=1/4）

全部 K∈{50,100,200,400}：**t=8 exact stop、N_R=0、零 rebound**。

| K | active | t_stop | N_R | f_R | R_t p5/p50/p95 | max R_t | moved total | max 残差 |
|---|---|---|---|---|---|---|---|---|
| 50 | 8 | 8 | 0 | 0.000 | 0.111 / 0.190 / 0.685 | 0.706 | 1.4222 | 4.8e-16 |
| 100 | 8 | 8 | 0 | 0.000 | 0.111 / 0.190 / 0.672 | 0.691 | 1.4221 | 1.1e-15 |
| 200 | 8 | 8 | 0 | 0.000 | 0.112 / 0.190 / 0.680 | 0.695 | 1.4223 | 8.3e-16 |
| 400 | 8 | 8 | 0 | 0.000 | 0.111 / 0.190 / 0.670 | 0.683 | 1.4209 | 1.8e-15 |

D_max 严格单调下降（0.25 → ~2×10⁻¹⁶）。定性签名：R_t 从 ~0.13 单调爬升到 ~0.69（K=100 的 t=7），**始终未穿过 1**；boundary 保持活跃（switch freq 0.571，|Δa| 中位 0.29）但没有反弹兑现。moved mass 1.42 vs K0 的 1.12——mismatch 使 exact stop 推迟 2 轮（6→8）但不产生任何反弹。

**解读（frozen contract W0）**：`positive mismatch existence alone is not sufficient to overcome the target-matched contraction margin along these trajectories`。（边界：不得写"positive mismatch never causes rebound"。）

## 5. KS 结果（spatial-shape control，flat-top）

| K | active | t_stop | N_R | f_R | S_R | S_R/K | rebound 轮 | max R_t | moved total | max 残差 |
|---|---|---|---|---|---|---|---|---|---|---|
| 50 | 26 | 26 | 5 | 0.192 | 11 | 0.220 | 13,15,21,22,24 | 2.558 | 4.129 | 1.0e-15 |
| 100 | 26 | 26 | 5 | 0.192 | 11 | 0.110 | 13,15,21,22,24 | 2.789 | 4.120 | 1.3e-15 |
| 200 | 26 | 26 | 5 | 0.192 | 11 | 0.055 | 13,15,21,22,24 | 2.601 | 4.120 | 4.1e-15 |
| 400 | 29 | 29 | 10 | 0.345 | 19 | 0.048 | 8,12,13,17,18,20,22,23,25,27 | 5.198 | 4.594 | 7.6e-15 |

定性签名（K=100）：D_max 0.25 → 0.049（t=13 首次反弹 +6.3×10⁻³）→ 0.056（t=15，+1.8×10⁻²）→ 回落 0.011（t=20）→ t=21 大反弹 +3.7×10⁻²（R=2.79）→ t=22、t=24 两次边际反弹（R=1.003、1.024）→ t=26 exact stop。rebound inter-arrival 中位 **2 轮**；longest contraction run 13；boundary switch freq 0.80；R_t p50≈0.75（KC 为 0.69）。**全部 25 个 rebound 轮的边界 a ≥ 0.66**。

**解读（frozen contract S1）**：`one-step rebound mechanism survives the shape change, while long recurrent dynamics remains shape-sensitive`。KS 证明 exact 单步判据在非抛物线 mismatch 下逐轮成立，但反弹只是短促复发（最长跨越 19 轮 / 0.05K），随后快速 exact stop——**没有出现 churn-like transient**。

**Scope boundary（冻结）**：KS 同时改变了 mismatch 正则性；KC vs KS 的差异只能表述为"对 mismatch 空间 profile / 正则性敏感"，不得说"flat top 导致差异"，不得说"shape 不重要"。

## 6. Cross-K 分析（§19 冻结问题）

> 第一轮只回答：phenomenon survives resolution changes or appears resolution-specific?

- **类别跨 K 稳定**：KW 全部 K = No rebound；KS 全部 K = Short recurrent transient；K0/K- 无反弹；KC 全部 K = Long-span churn-like。
- **无 censoring**：全部 17 条轨迹在 10K 之前 exact stop（最长 KC K=400 用 2876/4000）。
- **理论验证跨 K 稳定**：全部 active 轮 FP=FN=0，残差 ≤ 7.6×10⁻¹⁵。
- **定量细节有分辨率敏感性**：KS 的 rebound 轮次在 K≤200 完全相同（{13,15,21,22,24}），在 K=400 变为 10 次（t_first 8、max R 5.20）；S_R/K 随 K 单调下降（0.22→0.048，这是 S_R 有界而 K 增长的平凡标度+真实轻微增长共同作用）。KW 的 t_stop=8 与 R_t 分位在全部 K 近乎不变。

## 7. Table 1 — Kernel structure

| kernel | G | Δ(x)=G−T | sign | max Δ₊ | ∫Δ | 角色 |
|---|---|---|---|---|---|---|
| K0 | x | 0 | zero | 0 | 0 | target-matched integrity reference（M1C.1 U 复用） |
| K- | x² | −x(1−x) | ≤ 0（canonical 镜像） | 0 | −1/6 | sign control + 定理后盾 audit |
| KW | 5/4·x − 1/4·x² | x(1−x)/4 | ≥ 0 | 1/16 | 1/24 | amplitude ablation（λ=1/4 冻结） |
| KC | 2x−x² | x(1−x) | ≥ 0 | 1/4 | 1/6 | canonical anchor（integrity only） |
| KS | 分段 7/4·x \| x+1/4 \| 3/4+x/4 | 分段 flat-top | ≥ 0 | 1/4 | 1/6 | shape control（peak/integral 对齐，profile 不同） |

（max/integral 为解析值；grid 验证见 `results/static_kernel_validation.csv` 与测试套件。）

## 8. Table 2 — Dynamics

完整字段见 `results/rebound_summary.csv`；主列：

| kernel | K | active | stop | t_stop | N_R | f_R | S_R/K | sw freq | moved total | max resid | FP/FN | nb | category |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K0 | 50 | 6 | exact | 6 | 0 | 0.000 | — | 0.400 | 1.123 | 7.5e-16 | 0/0 | 0 | No rebound |
| K0 | 100 | 6 | exact | 6 | 0 | 0.000 | — | 0.400 | 1.123 | 2.1e-15 | 0/0 | 0 | No rebound |
| K0 | 200 | 6 | exact | 6 | 0 | 0.000 | — | 0.400 | 1.124 | 2.9e-15 | 0/0 | 0 | No rebound |
| K0 | 400 | 6 | exact | 6 | 0 | 0.000 | — | 0.400 | 1.122 | 5.4e-15 | 0/0 | 0 | No rebound |
| KC | 50 | 324 | exact | 324 | 117 | 0.361 | 6.320 | 0.923 | 49.02 | 7.9e-16 | 0/0 | 0 | Long-span churn-like |
| KC | 100 | 699 | exact | 699 | 251 | 0.359 | 6.860 | 0.930 | 103.30 | 1.1e-15 | 0/0 | 0 | Long-span churn-like |
| KC | 200 | 1445 | exact | 1445 | 531 | 0.367 | 7.175 | 0.941 | 210.59 | 1.6e-15 | 0/0 | 0 | Long-span churn-like |
| KC | 400 | 2876 | exact | 2876 | 1059 | 0.368 | 7.165 | 0.942 | 415.27 | 2.5e-15 | 0/0 | 0 | Long-span churn-like |
| K- | 100 | 3 | exact | 3 | 0 | 0.000 | — | 0.000 | 0.579 | 3.1e-16 | 0/0 | 0 | No rebound |
| KW | 100 | 8 | exact | 8 | 0 | 0.000 | — | 0.571 | 1.422 | 1.1e-15 | 0/0 | 0 | No rebound |
| KS | 100 | 26 | exact | 26 | 5 | 0.192 | 0.110 | 0.800 | 4.120 | 1.3e-15 | 0/0 | 0 | Short recurrent transient |
| KW | 50 | 8 | exact | 8 | 0 | 0.000 | — | 0.571 | 1.422 | 4.8e-16 | 0/0 | 0 | No rebound |
| KS | 50 | 26 | exact | 26 | 5 | 0.192 | 0.220 | 0.800 | 4.129 | 1.0e-15 | 0/0 | 0 | Short recurrent transient |
| KW | 200 | 8 | exact | 8 | 0 | 0.000 | — | 0.571 | 1.422 | 8.3e-16 | 0/0 | 0 | No rebound |
| KS | 200 | 26 | exact | 26 | 5 | 0.192 | 0.055 | 0.800 | 4.120 | 4.1e-15 | 0/0 | 0 | Short recurrent transient |
| KW | 400 | 8 | exact | 8 | 0 | 0.000 | — | 0.571 | 1.421 | 1.8e-15 | 0/0 | 0 | No rebound |
| KS | 400 | 29 | exact | 29 | 10 | 0.345 | 0.048 | 0.786 | 4.594 | 7.6e-15 | 0/0 | 0 | Short recurrent transient |

（nb = numerical-boundary 轮数；FP/FN 为 rebound criterion 误分；"Long-span churn-like" 是 [frozen 分类约定](config.json)下的 descriptive label，全签名见 CSV。）

## 9. 理论验证汇总（H）

对全部 17 条轨迹、**5510 个 active 轮**：

- **分解恒等式** `D_new = D_TM_new + M·Δ`：max |residual| = **7.55×10⁻¹⁵** ≤ τ_num=1e-12（发生在 KS K=400）。
- **精确 rebound 判据**（actual rebound ⟺ ∃j: M·Δ(j)>γ(j)）：**TP=1983、TN=3527、FP=0、FN=0**。
- **numerical-boundary cases**（|worst excess| ≤ 1e-12）：**0**。无需任何 boundary 解释，inequality 语义未做任何改动。
- 质量守恒：max |Σp−1| = 1.64×10⁻¹⁴ ≤ τ_num（KC K=100 的 fp 累积，canonical 语义未动）；min bin mass = 7.4×10⁻⁶ > 0。
- R_t（M1C.4 语义）与判据逐轮等价（R>1 ⟺ rebound），0 例外。

## 10. 图

- `fig_gate1_D_max_K100.png`（Figure A）：四 panel D_max(t)——KC 的 churn 平台 vs K-/KW/KS 的快速 exact stop。
- `fig_gate1_boundary_K100.png`（Figure B）：四 panel boundary a_t。
- `fig_gate1_rebound_R_K100.png`（Figure C）：R_t + rebound 标记——KW 的 R_t 逼近 1 但不穿；KS 的 5 次穿越。
- `fig_gate1_KC_KW_KS_comparison_K100.png`：KC vs KW vs KS（K- 参照）D_max 与 boundary 直接对比。

## 11. 与冻结预期对照

| 预期（config，运行前写定） | 结果 |
|---|---|
| K0：无反弹、快速停止（U t=6） | ✓ t=6 |
| K-：无反弹、严格下降、停止时间未知（可能 censoring） | ✓ t=3（无 censoring） |
| KW：判据允许反弹，dynamics unknown | t=8 无反弹 exact stop（W0） |
| KC：复现 324/699/1445/2876 | ✓ 132/132 |
| KS：dynamics unknown | 5–10 次短命反弹 + t=26/29 exact stop（S1） |

## 12. 完成度对照（§25 stop condition）

K0/KC integrity ✓、K- @K=100 ✓、KW/KS @ {50,100,200,400} ✓、theory verification ✓、summary metrics ✓、minimal plots ✓、tests（306 项运行时断言全通过，含本轮新增 18 项）✓、canonical report（本文 + docs/M2_GATE1_GENERALITY_GATE.md）✓、research-log update ✓、next-step checkpoint ✓ → **Gate 1 STOP 条件满足，PASS**。未追加任何 kernel、参数、扫描或 Gate 2 工作。
