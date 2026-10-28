# 2026-10-27 定向修正与验收记录

基准：[已完成评审](../review_2026-10-27/ASSESSMENT.md)及
[既有勘误 E1–E3](../review_2026-10-27/ERRATA.md)。本次仅落实这三项更正的
维护入口、引用范围与复验，不重新判断新颖性、不开展新数学研究。

**实际修改。** [READING_NOTES](../READING_NOTES.md) 顶部增加带日期的 E1–E3
索引，给出当前正确表述并列明历史原句及重复引用的位置。
[MODEL_ASSUMPTIONS](../MODEL_ASSUMPTIONS.md) 的 Phase 2A 条目明确：
`c_n>c_f` 是连续路径顶点判据；一般整数轮改善应按 E1 判断。
新增标准库核验入口 [verify_targeted_corrections.py](../../tools/verify_targeted_corrections.py)。

**保留原文的依据。** 既有 `ERRATA.md` 已给出两项 M0 修正的推导、精确反例
和 B 的斜率修正，无需重写。M0 主文、日期明确的报告/checkpoint/日志以及
已发布的 v0.3 B 证明均保留原字节，以明确的勘误索引限定其使用。这样也保留了
完成评审时记录的全部源文件哈希。维护中的模型假设台账则直接收缩其含糊表述。

| 已核对的引用位置 | 处理与当前含义 |
|---|---|
| `docs/M0_GENERAL_Q.md` 引言、§1.2、§1.4、§A.5–A.6 | 历史原文按 E1 读取。连续路径与顶点公式有效；其一般整数轮改善的充分条件已被更正。 |
| `experiments/m0_2_distribution_family/M0_2_REPORT.md` §2.1、§2.2、§7、§8 | 相同的过宽解释由 E1 限定；六个冻结案例的数值结果不变。 |
| `docs/RESEARCH_CHECKPOINT_2026-09-19.md`、`docs/RESEARCH_CHECKPOINT_2026-09-20.md`、`docs/RESEARCH_LOG.md` 的 M0 判据条目 | 日期明确的阶段材料，保留原文；阅读说明已明确列入 E1 的重复引用。 |
| `experiments/m0_2_distribution_family/run_experiment.py::predict` | 保留历史来源与输出模式。`predicted_interior_optimum` 检测连续顶点（含浮点容差），不承担一般整数轮改善定理；另算的 `predicted_best_round_discrete` 是记录窗口内的损失比较。 |
| `docs/M0_MECHANISM.md` §2.7、`experiments/m0_1_mechanism/run_experiment.py` 的 MC 均值比较及注释 | 历史精确期望声明按 E2 更正；旧脚本的固定均匀噪声参照是上界近似。未重生成或覆盖旧统计输出。 |
| `docs/V03_B_ROBUSTNESS_PROOF.md` §2(3) | 保留已发布证明；E3 和下面的推导明确替代错误上界，恢复原证明步骤。 |
| `src/sand_m0/model.py::noise_floor_l2` | 已明确只用于 iid 均匀样本，公式正确，无需改动。 |
| `tests/test_m0_2_family.py` | `Q_f*S_n>Q*S_f ⇔ z*<1` 的连续判据和有限案例检查正确，无需改动；不把其结果解释为一般整数轮证明。 |
| `src/sand_m0/plotting.py` 的相关图题 | “only when”给出必要条件；连续路径图的顶点标签也不等于一般整数轮充分条件，未作无关改动。 |

**E1 与 E2 的定向复核。** 对固定区域 M0 的二次损失，直接展开

\[
 V(\rho)-V(1)=A(\rho-1)(\rho+1-2z_*),\qquad A>0,\quad0<\rho<1.
\]

若首步不严格改善，则 \(\rho\le2z_*-1<z_*\)，以后所有 \(\rho^t\) 都远离顶点，
不能产生更优时刻。这验证 E1 的全时间逻辑。另检查等号时仅 t=0、1 并列，以及
α=0 静止情形。对 iid 直方图，从所有有限样本结果精确求期望，验证

\[
 \mathbb E V(\widehat p)=V(p)+\frac{1-\sum_i p_i^2}{N},
 \qquad V(p)+\frac{1-1/K}{N}-\mathbb E V(\widehat p)=\frac{V(p)}{N}.
\]

包含非均匀、均匀、零概率格及 N=1 的边界案例，并检查乘 K 后的归一化形式。

**E3 的推导核验。** 使用原 B 定理的假设 K≥7、α=1/4、
\(\max\{6/7,\Theta(K)\}<C<1\)。在选择 \(\ell=m-1\) 的一步，
\(0\le x=p_m\le h\)、\(w=F_\ell=C-x\)，因此

\[
 w\ge C-h>6/7-1/6=29/42,\qquad w\le C<1.
\]

对 \(H(y)=wG(y)-y\)，三段斜率分别为

\[
 w/2-1<-1/2<0,\quad 3w/2-1\ge1/28>0,\quad 2w-1\ge8/21>0.
\]

这里的非严格上界允许 x=0、w=C；既有 E3 中 `w<C` 的简写在该边界应读作
`w≤C`，但仍有严格的 w<1。完成评审的原文保留，这个通用边界写法已在维护
入口明确。第一个负号来自 w 的**上界**，其余正号来自 w 的下界。旧式
\(w/2-1<29/84-1\) 在 w=9/10 时被精确拒绝。

还需检查斜率符号能否完成后续论证。设 \(a_K=1-q-h/2\)，由 q≤2h、h≤1/7
可知 \(a_K\ge9/14>0\)。用原 \(\Theta(K)=h+(b-2h)/a_K\) 得

\[
 \begin{aligned}
 H(\ell h)-H(h)
 &=w a_K-b+2h\\
 &\ge(C-h)a_K-b+2h
 =a_K(C-\Theta(K))>0.
 \end{aligned}
\]

结合 H 的先降后升，得到全部内部格点 \(H(jh)<H(\ell h)\)。两者都被本轮前缀
更新扫过，故直接从累积更新展开

\[
 D'_\ell-D'_j=\frac34(D_\ell-D_j)+\frac14[H(\ell h)-H(jh)]>0.
\]

首次选择边界时内部平局已由最小指标规则排除。这就验证修正后的斜率论证、
端点因子及真实一步更新之间的连接；不是只替换一个符号。有限回归使用满足
原假设的完整质量状态，含 x=0 零质量边界及 x=h 平局，交叉检查上述恒等式及内部严格支配。
一般结论依赖这些解析不等式，不由有限案例数承担。

**实际执行结果。** 下列命令均以 Python 3.14.6、`-B` 执行，真实退出码均为 0。

| 检查 | 结果 | 本次证据 |
|---|---|---|
| `python -B tools/verify_targeted_corrections.py` | E1/E2 原反例、4 个 M0 边界案例、9 组 iid 精确枚举、B 的斜率/端点检查及 39 组完整状态更新全部通过 | [targeted-verification.json](targeted-verification.json) |
| `python -B tests/test_m0_mechanism.py` | 25 项通过 | [m0-regressions.json](m0-regressions.json) |
| `python -B tests/test_m0_2_family.py` | 54 项通过，冻结案例最佳轮次等既有锚点不变 | [m0-regressions.json](m0-regressions.json) |
| `python -B tools/audit_v040_review.py` | 18 个形式矩恒等式、93 条有限轨迹、21 组共轭、原反例及 E1/E2 通过 | [review-replay.json](review-replay.json) |
| `python -B tools/verify_v040.py` | 6/6 通过；37 个冻结输入前后未变；约 49 秒 | [v040-verification.json](v040-verification.json) |

保存报告时给命令加 `--report` 及一个新的工作区路径。定向核验器依赖同次评审
留下的 `tools/audit_v040_review.py`，仅用标准库；报告以排他创建模式写入，避免
覆盖既有证据。`targeted-checks*.json` 保留较早检查时的源哈希，最终维护文件以
`targeted-verification.json` 绑定的字节为准。M0 测试的实际命令、stdout 和退出码
完整保存在 `m0-regressions.json` 中。

**影响评估。** 历史 M0 的一般整数轮改善断言被限定为正确的条件；它的连续路径
公式和六例数值结果不变。MC 的错误精确等式改按 E2 读取，旧近似比较不能作为
该错误等式的验证。B 的原假设与驻留结论不变，引用该步骤的证明无需改变结论。
v0.4.0 主定理的模型、假设、科学声明和证明依赖均不变。

已检查实际 diff：仅修改阅读说明和模型假设台账的上述位置，新增本次定向核验器
及独立验收记录。原 M0 文稿、B 证明、实验来源/输出、完成评审、公开快照、
release 记录和所有 Git 标签/历史均保留；最终哈希和差异检查另见
[integrity.json](integrity.json)。本次三项修正没有遗留的数学修补项；不扩大既有
文献新颖性或研究方向等未决问题的评审范围。
