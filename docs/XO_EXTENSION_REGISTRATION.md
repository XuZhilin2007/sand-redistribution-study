# XO 补充成果正式登记 / XO supplement registration

登记版本：v0.2.0；日期：2026-10-07；作者：Zhilin Xu。
本登记将 A、T6 分别收录，原 [S3 登记](S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)、[证明](S3_XO_BOUNDARY_THEOREM_PROOF.md)和[审计](S3_XO_BOUNDARY_THEOREM_AUDIT.md)保持原范围。

## 共同模型 / Common model

一维有限维确定性质量反馈，单位总质量，canonical 初态 $`p_{0,i}=(2K-2i+1)/K^2`$，XO 固定喷撒律，最小指标 argmax 平局规则，精确算术，整数 $`K\ge6`$。更新前计算 $`d_t=\max_j(F_{t,j}-j/K)`$，在 $`d_t\le10^{-12}`$ 时停止。定义 $`m=\lceil5K/6\rceil`$、$`\ell=m-1`$、$`h=1/K`$、$`b=m/K`$、$`c=2b-b^2`$、$`a=b(1-b)`$、$`q=(6-K\bmod6)/(3K)`$。

The common assumptions are canonical initialization, the fixed XO law, smallest-index ties, exact arithmetic, integer K≥6, and the existing pre-update tolerance rule. Changing the initial state, law, or tie rule changes the theorem's scope.

## A：操作力度推广 / Removal-fraction extension

**前提：共同模型，且 $`0<\alpha<1`$。**
**陈述：** T1′–T5′ 成立：守恒与正下界；有限步进入边界区；首次进入后保持；精确边界坐标闭包；长期选择按 K mod 6 二分。进入时间满足

$$
t_*\le\left\lfloor\frac{81K}{80\alpha}\right\rfloor+1.
$$

进入后，令 $`x_t=p_{t,m}`$，则

$$
s_t=\begin{cases}m,&x_t>h,\\
\ell,&x_t\le h,\end{cases}
\qquad d_t=a+(h-x_t)^+,
\qquad
x_{t+1}=\begin{cases}(1-\alpha)x_t+\alpha cq,&x_t>h,\\
x_t+\alpha(c-x_t)q,&x_t\le h.\end{cases}
$$

余数 0、1、2 最终永远选择 m；余数 3、4、5 两边界被无限次选择。**不主张周期性、完整状态渐近或任意 α 的 H=10K 分类保证。**

证明：[A 正式收录稿](XO_ALPHA_EXTENSION_PROOF.md)。依赖：原 S3 中的 XO 几何/不变量方法及 A 自身的 α 推广、K=6 特例和长期标量分析；A 的证明没有以 T6 为前提。
验证：255 组完整状态 Fraction 矩阵（大 K 的实际时限为 300）、7 组 K=6 补充、255 组各 200000 步的有限浮点标量代理；这些计算不承担全称量词。实际配置和发布候选结果见[验收](XO_EXTENSION_VALIDATION.md)。

**A extends the long-term boundary theorem to every 0<α<1 under the common assumptions.** It proves finite entry, persistence, scalar closure of boundary observables, and the residue dichotomy. It does not prove periodicity or finite-window classifier reliability for arbitrary α.

## T6：冻结有限窗分类器 / Frozen finite-window classifier

**前提：共同模型，$`\alpha=1/4`$，$`H=10K`$。**
窗口恰有 H 个更新前选择项，时间为 t=0,…,H−1。使用原 [S2 的 classify_itinerary](../experiments/s4_w1_s2_boundary_microscope/run_s2.py)，其优先顺序为末尾 m 连续项数 R_m≥K 时 LOCK；否则首次进入后集合恰为 {ℓ,m} 时 CYCLE；否则 OTHER。

**陈述：所有整数 K≥6，余数 0、1、2 当且仅当分类为 LOCK；余数 3、4、5 当且仅当分类为 CYCLE；没有 OTHER。**

证明：[T6 正式收录稿](XO_T6_CLASSIFIER_PROOF.md)。证明只依赖原 S3 T1–T5 及自身驻留/计数界，**不依赖 A**。原 S2 分类函数与冻结材料未改。
验证：15 项精确条件、K=6,…,300 共 295 个实例，以及从仓内 S2 原文件提取函数的交叉核对；29 条独立整数完整状态轨迹另列在验收记录。实例数量不代替全 K 的证明。

**T6 establishes reliability of the unchanged classifier for every K≥6 at α=1/4 and H=10K.** Its mathematical dependency is original S3, not A. Its finite window contains the pre-update choices at t=0,…,H−1.

## 必须保留的范围外反例 / Scope counterexamples

| 参数 / Parameters | 冻结标签 / Frozen label | 与长期行为的关系 / Interpretation |
|---|---|---|
| K=9, α=1/8, H=90 | LOCK；t*=3，R_m=10，post-set={7,8} | A 仍给出无限双支切换；LOCK 优先顺序使该窗口不能保证长期分类。 / A still gives infinite switching; the finite label need not match it. |
| K=10, α=1/1000, H=100 | OTHER；窗口内未进入 | A 的有限进入时间未必落在 10K 内。 / Finite entry need not occur within 10K. |

## 审查与当前状态 / Review and status

2026-10-07 的独立只读技术审查在 A、T6 各自声明范围内未发现实质性证明漏洞；Owner 随后接受审查并授权 v0.2.0 入库发布。详细来源与计算边界见[单独的补充审查记录](XO_EXTENSION_REVIEW.md)。这不扩展原 S3 的审计来源，也不构成期刊同行评审、形式化证明或文献新颖性认定。

辅助 F、首次选择公式、限定探针、径向聚合和文献线索的定位见[辅助附录](XO_SUPPLEMENTARY_RESULTS.md)。新版本不以未完成的新研究终稿、周期性理论、一般噪声或二维研究为发布前提。
