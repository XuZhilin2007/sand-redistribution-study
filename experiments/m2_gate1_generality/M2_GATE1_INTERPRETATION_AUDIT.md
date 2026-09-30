# M2 Gate 1 Interpretation Audit — Rebound-Witness Spatial Localization

> 日期：2026-09-20。性质：**analysis / documentation-only 轮**（Owner + GPT 冻结 prompt）。无新 redistribution kernel、无 dynamics 修改、无 Gate 1 重跑、无 Gate 2 模拟。
>
> 唯一研究问题：**KC、KS 中真正触发 rebound 的 witness prefix 位于哪里？KW 在最接近 rebound 但始终未穿阈值时，最危险的位置在哪里？**
>
> 运行方式：`python experiments/m2_gate1_generality/run_interpretation_audit.py`。helper 测试：`tests/test_m2_gate1_audit.py`（全库 313 项运行时断言通过）。

## A. Repository Reality

- 起始/结束 HEAD：`7cf0fa3`（Gate 1 docs commit），branch `main`，working tree clean。
- Gate 1 资产在位：`M2_GATE1_REPORT.md`、`docs/M2_GATE1_GENERALITY_GATE.md`、`results/per_round_diagnostics.csv`（5510 行）等。
- 已冻结事实（本轮不改写）：KW 全部 4 K 零反弹、t=8、max R_t≈0.68–0.71（W0）；KS 全部 4 K 5–10 次反弹、t=26–29（S1）；KC long pre-stopping churn（本轮不复刻 M1C.4–6）。

## B. Audit Method

**数据来源与忠实性门禁**：committed per-round CSV 只保存标量，且其 `j_drive` 是 R_t **比值**的 argmax（M1C.4 语义），不是本轮定义的 E **差值** argmax；E_t(j) 依赖完整 D-profile，无法从标量恢复。因此 audit 从冻结 config 确定性再生轨迹（与 Gate 1 完全相同的 `run_m1a_deterministic` 调用），并在任何派生之前对**全部 5483 个 active 轮**逐字符串核对 committed 字段（D_max_before、j_t、a_t、moved_mass、rebound、worst_excess、R_t）——全部一致，零失配。再生核验不产生任何新 dynamics 结果，未触碰 Gate 1 raw results。

**Witness/risk 定义（冻结）**：

```text
E_t(j)      = M_t·Δ(j) − γ_t(j)        = M1C.3/M1C.5 rebound_residual helper（直接复用，零新语义）
γ_t(j)      = D_max(t) − D_TM_new,t(j)
j_witness   = argmax_j E_t(j)          （tie → 最小 j，np.argmax 语义，与 repository 边界 tie-breaking 一致）
x_witness   = j_witness / K
```

- KC/KS：仅对 **actual rebound rounds**（max_j E_t(j) > 0）提取 witness。
- KW：零反弹，argmax 是 **risk location**（max E ≤ 0），绝不称 witness。
- **结构性标记**：由 M1C.3 suffix-impossibility（E(j) ≤ 0 ∀ j > j\*；E(K) = −D_max 恒为平凡值），raw argmax 可能落在结构性不可能的 suffix 位置。audit 记录 `argmax_on_suffix` 标记，并补充 **prefix-restricted argmax**（j ≤ j\*，criterion-relevant）作为 supplementary 列——不改变冻结定义，只防止误读。

区域分类（冻结）：Near x<1/3，Middle 1/3≤x≤2/3，Far x>2/3。

## C. KC Witness Localization（actual rebound witnesses，1958 轮）

| K | n | x_witness median | p05–p95 | Near | Middle | Far | median \|x−a\| |
|---|---|---|---|---|---|---|---|
| 50 | 117 | 0.480 | 0.240–0.624 | 18.8% | 80.3% | 0.9% | 0.380 |
| 100 | 251 | 0.460 | 0.230–0.610 | 21.5% | 78.1% | 0.4% | 0.390 |
| 200 | 531 | 0.465 | 0.200–0.615 | 23.0% | 76.6% | 0.4% | 0.385 |
| 400 | 1059 | 0.455 | 0.203–0.613 | 24.6% | 75.2% | 0.2% | 0.393 |

**Verified Numerical Analysis / Model-dependent Observation**：

1. canonical rebound witnesses 高度集中在一个**稳定的 interior 带**：median ≈0.46，p05–p95 ≈ [0.20, 0.62]——正是 Δ_C(x)=x(1−x) 峰值（x=0.5）周围的 hump 区。跨 4 K 分布几乎不动（median 0.455–0.48，各区域 share 变化 ≤2.4 个百分点）。
2. witness 通常**远在扫掠边界左侧**：median |x−a| ≈ 0.38–0.39（rebound 轮的 a_t 典型 ≥0.6），不是贴近 boundary 的位置。
3. **64.0% 的 witness 在 x ≤ 1/2**（其中 [0.2, 0.5] 占 59.8%），35.6% 在 (1/2, 2/3]，**x > 2/3 仅 0.31%**（6/1958；per-K 0.19–0.85%）。
4. 0/1958 的 witness 落在 suffix（结构性不可能位置）——与 M1C.3 结构定理一致。

## D. KS Witness Localization（actual rebound witnesses，25 轮）

| K | n | x_witness median | p05–p95 | Near | Middle | Far | median \|x−a\| |
|---|---|---|---|---|---|---|---|
| 50 | 5 | 0.340 | 0.324–0.596 | 20% | 80% | 0% | 0.340 |
| 100 | 5 | 0.340 | 0.330–0.596 | 40% | 60% | 0% | 0.330 |
| 200 | 5 | 0.335 | 0.331–0.599 | 20% | 80% | 0% | 0.335 |
| 400 | 10 | 0.333 | 0.333–0.668 | 60% | 10% | 30% | 0.333 |

**核心观察（Model-dependent Observation）**：

1. **KS witness 不散布在 flat-top 内部，而是钉在 flat-top 的两个边界点 x≈1/3 与 x≈2/3 上**：25 个 witness 中 19 个在左边缘（≈1/3 侧）、6 个在右边缘（≈2/3 侧），**全部 25 个与最近 flat-top 边缘的距离 ≤ 0.0133**（多数 ≤ 1/(2K) 量级）。K=400 的 Near/Far share 波动是格点落在 1/3、2/3 哪一侧的网格伪影，不是真实的位置迁移。
2. 因此"frozen 三分区"（Near/Middle/Far）对 KS 的定位分辨力有限：witness 卡在 1/3 分类边界上，区域 share 随 K 的格点侧向翻转——比较 KC vs KS 应看连续位置（median ≈0.33 vs 0.46、双 edge 钉扎 vs 单 hump 带），不应只看三分区 share。
3. KS rebound witness 全部落在 a_t ≥ 0.66 的轮（同 Gate 1 记录），median |x−a| ≈ 0.33——同样远在边界左侧。
4. 0/25 witness 落在 suffix。

**问题回答（数值观察，无因果 claim）**：(1) 是——KS witness 集中在 flat-top 边缘点（≈1/3 主、≈2/3 次），不是整个 [1/3,2/3]；(2) KC 与 KS 明显错位（interior hump 带 vs 边缘钉扎）；(3) KS 的短 transient 对应 witness 在两个边缘点之间**快速来回迁移**（见 G），且 19–29 轮内消失（exact stop），不是某个位置的持续驻留；(4) KS witness **基本不在 flat-top 外部**——它们恰在 flat-top 的边界（斜率不连续点）上，即 Δ_S 从 ramp 转常数的转折处。

## E. KW Near-Rebound Risk Localization（risk location，非 witness）

KW 零反弹：32 个 active 轮（4 K × 8）的 max E 全部 ≤ 0（−gap 最小 −2.45×10⁻³）。

**Raw argmax（冻结定义）的警示**：24/32 轮的 raw argmax 落在 **suffix（结构性不可能位置）**——每 K 8 轮中 6 轮。closest round（4 K 一致为 t=7）的 raw argmax 在 x=1.0，E(K) = −D_max ≈ −2.5×10⁻³：这是"全场 E 都很负"时 argmax 退到平凡终点的产物，**不是**有意义的 rebound 候选位置。因此 KW 的 risk 定位必须读 prefix-restricted 列。

**Prefix-restricted（j ≤ j\*，criterion-relevant）**：

| 样本 | x median | p05–p95 | Near | Middle | Far |
|---|---|---|---|---|---|
| 全部 active 轮（pooled） | 0.223 | ≈[0.00, 0.66] | 59.4% | 37.5% | 3.1% |
| top 25% R_t 轮 | 0.119 | — | **100%** | 0% | 0% |
| 每 K closest round（t=7） | 0.0025–0.02 | — | 4/4 | 0 | 0 |

**回答（不把位置解释成因果）**：KW"接近但始终没穿 1"在 criterion-relevant 意义上主要发生在 **Near 区（x≲1/3，甚至贴近原点）**——relative 压力最高的轮（top R_t）与 margin 差距最小的轮的 prefix-restricted argmax 全部在 Near。注意与 R_t **比值** argmax（committed `j_drive`，KW t=7 K=100 时 x=0.55）是不同位置：差值 argmax 找"绝对 margin 最薄"处（小 j 处 target-matched 更新几乎不动 profile、γ 最小），比值 argmax 找"再注入/余量比"最大处。两个视角都保留在输出 CSV 中。

## F. KC vs KS Comparison

| 维度 | KC（1958 轮） | KS（25 轮） |
|---|---|---|
| witness 位置形态 | 单 interior hump 带，median ≈0.46，p05–p95 [0.20, 0.62] | 双 edge 钉扎：19/25 @≈1/3，6/25 @≈2/3，全部距 edge ≤0.013 |
| 与 Δ 形状的关系（描述性） | 覆盖 Δ_C 抛物 hump 内部（峰 0.5 周围） | 钉在 Δ_S 斜率不连续点（flat-top 边缘），不在幅度平台内部 |
| 跨 K 稳定性 | 分布几乎不动（4 K share 差 ≤2.4pp） | 同样两个 edge（跨 K 不动），仅轮数与格点侧翻转 |
| witness 生命周期 | churn 窗口内反复出现（S_R/K 6.3–7.2） | 19–29 轮内出现后随 exact stop 消失 |
| median \|x−a\| | 0.38–0.39 | 0.33 |

**两者的共同点**：witness 都远在扫掠边界左侧、都不在 Far（KC 0.3%、KS 0% 除 K=400 的 2/3 边缘格点侧）。**不同点**是位置形态：continuous hump 内部 vs 边缘点钉扎。这提示（仅 descriptive）：witness 位置跟踪的是 **Δ−γ 竞争的空间结构**，不是 Δ 幅度最大的点——KS 的 Δ 幅度在 [1/3,2/3] 内是常数，witness 却选斜率转折点。

## G. Rebound-to-Next-Round Movement

（已做：无需重构任何 pipeline——movement 从再生核验后的逐轮 argmax-E 序列直接读取。）

| kernel | n | median \|Δx\| | Δx p05–p95 |
|---|---|---|---|
| KC | 1958 | 0.300 | [−0.343, +0.710] |
| KS | 25 | 0.335 | [−0.334, +0.667] |

**没有固定 witness 驻留**：rebound 之后最危险位置发生大幅迁移（median ≈0.3），KC 的 p95=+0.71 显示可以跳到很远的前缀。KS 的位移分布两端 ±≈1/3、+2/3——与"两个 edge 点之间来回交替"一致（如 K=100：t=21 在 0.66 → t=22 在 0.34 → t=24 在 0.33）。这把"one-step rebound generation"与"multi-step state reconstruction"清楚分开：**单步的触发位置每轮重建，没有 persistent hot spot**（Model-dependent Observation；对 canonical 轨迹这与 M1C.6 的 reinjection-reshaping 结论一致，本轮未复刻其 pipeline）。

## H. Implication for Existing Gate 2 Candidate（descriptive assessment）

候选 `Δ_X(x) = x(1−x)(1−2x)`：x<1/2 为正（峰 +0.096 @ x≈0.211），x>1/2 为负（谷 −0.096 @ x≈0.789）。对照本轮定位结果，它会测试：

1. **对 canonical witness 带是 mixed intervention**：KC witness 的 64% 在 x≤1/2——Δ_X 在这半个带上**保持/加强**正失配（最多 +0.096）；(1/2, 2/3] 的 35.6% witness 区被**压制**（负瓣覆盖）；x>2/3 本来就几乎无 witness（0.3%）。也就是说 Δ_X 同时保留了大半 canonical witness-support region、只压制其上部——不是对 witness 带的干净"全抑制"或"全保留"实验。
2. **1/2 split 不是最有辨识力的 boundary（descriptive）**：KC witness 带的实际上界是 p95≈0.62，Far 侧近乎空；KS 的 witness 在 1/3 与 2/3——三组证据都不以 1/2 为自然分界。以 1/2 分正负更像是对称性驱动的选择，而不是 witness 定位驱动的选择。
3. **对 KS 是非对称边缘干预**：Δ_X 在 x=1/3 处 ≈ +0.074（保持左 edge 的正压），在 x=2/3 处 ≈ −0.074（抵消右 edge）——它会区分两个 KS witness site，但同时改变了 Δ 的幅度与正则性，读数不干净。
4. **幅度混杂**：|Δ_X|max ≈ 0.096，介于 KW（0.0625）与 KC 峰（0.25）之间——任何 sign-changing 结果都无法与"幅度降到 ~0.1"的效应区分（KW 已证明该量级附近可以零反弹）。
5. KW 侧的一个 observation（非设计建议）：KW 的 criterion-relevant risk 位置 100% 在 Near（top-R_t 轮 median x≈0.12）——若未来 Gate 2 想检验"压制 risk 区能否消灭反弹苗头"，负瓣的位置应参照 witness/risk 定位而非对称中点。

**结论（descriptive）：现有 candidate Δ_X = x(1−x)(1−2x) 的 1/2 split 与 witness 定位数据不匹配（mixed intervention + 幅度混杂），Gate 2 kernel should be redesigned。** 本轮不提出、不实现任何替代 kernel。

## I. Documentation Correction（措辞收紧，不改结果）

按冻结 prompt §13，已做最小修正（证据标签不变；实验数值零改动）：

1. `docs/M2_GATE1_GENERALITY_GATE.md` §8.1：原"canonical churn 在目前两个对照点上更像是 **full-amplitude smooth-parabolic profile 的特化行为**，不是该 kernel 类的普遍现象"（过强的特化性断言）→ 改为："当前证据表明，canonical long-span churn 依赖于比正 mismatch 的存在、峰值和积分强度更细的 mismatch 结构；amplitude 与 spatial profile/regularity 均会显著改变观察到的 dynamics，但维持 long transient 的具体多步机制仍未解决。"
2. 同处及 `docs/NEXT_STEPS.md`（Gate 1 追加条目）：否定对象的措辞统一限定为"**不支持 `positive mismatch alone broadly produces canonical-like churn` 这一强 generality claim**"，并显式标注"不得升级为 mismatch-driven churn 不是 model-class phenomenon"。
3. `experiments/m2_gate1_generality/M2_GATE1_REPORT.md` 检查后**无需修改**（其措辞本已限于"在两个 positive 对照上都没有出现"的现象陈述，无特化性断言）。
4. `RESEARCH_LOG.md` 为 append-only：Gate 1 旧 entry 不改写（其中含同类特化性短语），以本 audit entry 的收紧表述为准（见下）。

## J. External Review Candidate Questions（本轮新增）

- **ERC-5**：rebound witness localization 与 evolving contraction margin γ_t(j) 是否存在低维或 operator-level 表示，可预测 multi-step recurrence（witness 每轮重建、无 hot spot、S_R 与 stopping 的联合结构），而不仅是逐轮 rebound？（与 ERC-1 的逐轮对齐刻画互补：ERC-1 面向单步/静态泛函，ERC-5 面向轨迹级算子结构。）

（已登记进 [docs/M2_GATE1_GENERALITY_GATE.md §11](../../docs/M2_GATE1_GENERALITY_GATE.md)。）

## K. Gate 2 Readiness

```text
READY TO DESIGN GATE 2
```

依据：witness/risk 定位在全部冻结 K 上完成、跨 K 稳定、给出门控设计需要的三条约束（canonical witness 带 [0.20, 0.62] 双侧跨 1/2；KS witness 钉扎在 1/3、2/3 边缘；KW criterion-relevant risk 集中 Near）；现有候选 Δ_X 已给出 descriptive 评估（建议 redesign）与理由。**按 prompt 停止：不设计、不实现、不运行 Gate 2，等待 GPT + Owner。**

## 资产索引

| 交付 | 路径 |
|---|---|
| Audit 脚本（含再生核验门禁） | `experiments/m2_gate1_generality/run_interpretation_audit.py` |
| 逐轮 witness/risk 行（5483 轮，含 prefix-restricted 与 suffix 标记） | `results/witness_localization_rows.csv` |
| §9 汇总表（KC/KS witness；KW 三样本定义 × raw/prefix 两视图） | `results/witness_localization_summary.csv` |
| Rebound→next movement | `results/witness_next_movement.csv` |
| Figure 1（位置分布） | `results/witness_location_distribution.png` |
| Figure 2（rebound→next relocation） | `results/witness_next_movement.png` |
| Helper 测试 | `tests/test_m2_gate1_audit.py`（7 项；全库 313 项通过） |
