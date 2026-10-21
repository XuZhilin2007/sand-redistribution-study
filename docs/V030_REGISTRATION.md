# v0.3.0 成果登记 / Result registration

版本：v0.3.0。作者：Zhilin Xu。日期：2026-10-21。
登记状态：本页所列成果已通过本次发布验收，纳入 v0.3.0。
各证明形成时的“候选”标记按登记前的历史状态读取，证据等级仍按本页逐项区分。
Owner 已授权完成验收后合并并发布。本登记的发布条件为
[验收记录](V030_VALIDATION.md)列明的检查全部通过，远端版本以
`v0.3.0` 标签及对应 Release 为准。

本版以 v0.2.0 的 canonical XO 成果为基础，收录以下模型内数学结果。
S1/S2 的冻结有限窗口记录、原 S3 审计及 v0.2.0 登记保持各自原范围。
新结论由本版的独立证明与证据补充，不追溯改写早期证据等级。

## 1. 主成果：全部 K 的两轮最优停止

采用 canonical 初态
$`p_{0,i}=(2K-2i+1)/K^2`$、XO 喷撒律、$`\alpha=1/4`$、
精确算术、最小指标 argmax 平局规则及整数 $`K\ge6`$。设

```math
U(p)=K\sum_{i=1}^K\left(p_i-\frac1K\right)^2.
```

**全部整数 $`K\ge6`$、全部整数 $`t\ge0`$ 满足**

```math
U(p_t)\ge U(p_2).
```

因此执行恰两次更新给出一个全局最优停止时刻。结论不声称最优时刻唯一、
不声称后续 U 单调，也不适用于任意操作力度、初态、喷撒律或均匀度指标。
K=6 的首次边界选择为 t=1，最优停止仍可取 t=2；“首次边界选择时停”的
表述只在 K≥7 的 canonical 设置使用。

证明分为 [B-lock：LOCK 全类](V03_B_LOCK_ALLK.md)与
[E：CYCLE 全类](V03_E_CYCLE_ALLK.md)。B-lock 覆盖三个 LOCK 余数类；
E 的 96 子族配方证书覆盖三个 CYCLE 无限族，有限差集
{9,10,11,16,17,23} 由 [D](V03_D_STOPPING.md) 的全时间证书覆盖。
证据等级为**解析引理与精确计算辅助证明**，而非仅在大范围模拟中未见反例。

**For canonical XO, exact arithmetic, smallest-index ties, α=1/4 and every
integer K≥6, two updates attain a global minimum of U over all integer times.**
The proof combines analytic reduction with exact certificates for infinite
families and certified finite exceptions. It does not assert uniqueness,
time-monotonicity, or the same optimum for other models or metrics.

## 2. 同版登记的结果与适用范围

| 成果 | 本版登记范围与证据 |
|---|---|
| [T7：最终选择字与渐近状态](V03_T7_CYCLICITY_PROOF.md) | canonical α=1/4；每个 CYCLE K 的选择序列最终精确周期，质量状态渐近周期；显式字表、小 K 例外与尾部族覆盖。全域吸引形式分别覆盖 r=3 的 K≥195、r=4 的 K≥16、r=5 的 K≥65。 |
| [B：初态稳健性；B-II：有限回返](V03_B_ROBUSTNESS_PROOF.md) | B 要求原文的后缀支配、C∈(max{6/7,Θ(K)},1) 与正锚等条件；B-II 按原文全部合法可行初态范围，有限次内部回返后永久驻留两边界。C=1 和“首入不退”不纳入 B 的保证。 |
| [C：任意 α 的分类窗口](V03_C_ALPHA_WINDOW_PROOF.md) | canonical、0<α<1；新两测规则在声明窗口内正确；原冻结规则的 N_A^α<K 条件为充分条件，不写成 iff。 |
| [D：有限格点全时间证书与前缀分解](V03_D_STOPPING.md) | 203 个声明格点各自覆盖无限时间；102 条 CYCLE 记录有独立验证器。D-decomp 仅需两支驻留，不要求先锁定周期字。全 K 推论由 B-lock/E 承担。 |
| [F：旋转数与平台核对](V03_F_PLATEAU.md) | ρ 等于长期 B 频率；混合区短游程密度 d 满足 ρ=1/(10−d)。301 格点的平台/机械字核对保持有限认证等级。 |
| [G：一般 α 词典；G-ep：端点例外](V03_G_ALPHA_DICTIONARY.md) | 全 0<α<1 的 CYCLE 共轭词典；代数 α 可继承代数参数的有理旋转数结果。K=10 的明确代数 α* 实现无真正周期轨道的端点例外。 |
| [P：无限平台族](V03_P_PLATEAU_ALLK.md) | canonical 三个无限 CYCLE 区段，阈值 195/16/65；多项式符号证书与有限格点分类共同排除 canonical 右端点。 |
| [C-ep：端点精确表](V03_C_ENDPOINTS_ALPHA.md) | 13 个明确 (K,ρ) 目标的精确根数、隔离括号及极点排除；不升级为任意 K 的端点单调性定理。 |
| 文献继承定位 | 见[原文与继承说明](V030_RELATED_WORK.md)。有理旋转数、平台公式及适用情形下的周期吸引来自已有理论；本模型的归约、显式分类、量化证书和最优停止分别陈述。 |

LOCK/CYCLE 在旧 S2 中仍是其冻结有限窗口标签；本版的无限时间结论由新证明
给出。选择字最终周期与状态渐近周期始终区分。

## 3. 保留的反例与例外

- K=7、C=9/10 可在首入边界后回返内部：“首入不退”为假。
- K=9、C=1 是 B 的端点反例；其严格 C<1 假设不能删除。
- K=10 的代数 α* 是 G-ep 的真实端点例外，不能推广 canonical 的实际周期轨道结论。
- K=6：首次边界进入 t*=1，最优停止 t=2。
- K=10 的第二次选择 s₁=7，是通用前缀闭式的小 K 例外；E 使用 D 全时间证书。
- LP 固定 10⁻⁹ 间隙的不可行只排除达到该共同间隙的轨道；严格不可行的
  共同 ε 最大化仍是后续工作，B-II 的承载证明是解析推导。

## 4. 复现与后续

公开下载包首次运行前，先按 README 执行
`python -B tools/materialize_v030_certificates.py`，将 Git 中无损压缩的
两份大型证书恢复到原路径并核对原始 SHA-256。完整私有历史已有原 JSON，
可直接核验。

Windows PowerShell 先设 `$env:PYTHONUTF8="1"`，再运行
`python -B tests/count_tests.py`；macOS/Linux 可使用
`PYTHONUTF8=1 python -B tests/count_tests.py`。不要使用 `-O`。
固定依赖与完整命令见 [README](../README.md)。

现有结果附带生成器、JSON/CSV 证据、验证器与篡改拒绝回归。
有限轨迹核查与无限族证书的证据等级分开记录；“独立实现”不指形式化证明。

一般 α 的最优停止与完整相图、非 canonical 初态的停止问题、指标依赖性、
核/网格扰动稳健性、一般端点分类与混合字密度解释留作后续。
本版是确定性模型研究快照，不是经实测校准的沙粒物理模型，也不宣称已经
确立文献新颖性或完成期刊同行评审。
