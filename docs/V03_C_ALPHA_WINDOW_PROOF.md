# 定理 C — 任意移除比例下的有限窗口分类

> 研究分支 `research/v03-frontier`，2026-10-15。v0.3 frontier 成果之三。
> 依赖：定理 A（T2′–T5′，全 $`0<\alpha<1`$）与本文的 $`\alpha`$-驻留界。
> 精确验证入口 `experiments/v03_c_alpha_window/`。

## 1. 定义

对 $`K\ge6`$、$`0<\alpha<1`$（canonical 初态、XO 律、最小指标平局、精确算术），记
$`\lambda=c(6-r)/3`$。定义（均可由 A 的标量闭包直接导出）：

- 首次进入界 $`T_\alpha=\bigl\lfloor\tfrac{81K}{80\alpha}\bigr\rfloor+1`$（A, T2′）；
- 初始 B-游程界 $`N_B^\alpha=\min\{n:\ (1-\alpha q)^n c<c-h\}`$；
- 首个 A-游程界 $`U_1^\alpha=\min\{n:\ (1-\alpha)^n(c-L)\le h-L\}`$；
- 着陆后 A-游程界 $`N_A^\alpha=\min\{n:\ (1-\alpha)^n\bigl(h-(1-\alpha)L\bigr)\le h-L\}`$；
  显式 $`N_A^\alpha\le\bigl\lceil\ln 81/\ln\tfrac1{1-\alpha}\bigr\rceil`$（CYCLE 类
  $`\lambda\le\tfrac{80}{81}`$）；$`\alpha=\tfrac14`$ 时 $`N_A^\alpha\le11`$（T6 §3.3 同源）。

**两测分类器 $`\mathcal C_{2}`$**（窗口 $`t=0,\dots,H-1`$）：令 $`\tau_0=t_*+\max(N_B^\alpha,U_1^\alpha)`$，
输出 CYCLE 当且仅当存在 $`t>\tau_0`$ 且 $`t<H`$ 的 $`\ell`$-选择；否则 LOCK。

## 2. 陈述

**定理 C.** 取窗口

```math
H(\alpha,K)=T_\alpha+\max\bigl(N_B^\alpha,U_1^\alpha\bigr)+N_A^\alpha+2 .
```

则对一切 $`K\ge6`$、$`0<\alpha<1`$：**两测分类器** $`\mathcal C_2`$ 在该窗口上的输出恰为
$`r\in\{0,1,2\}\Rightarrow`$ LOCK，$`r\in\{3,4,5\}\Rightarrow`$ CYCLE。

**冻结优先分类器的充分条件（非当且仅当）。** 在更长的窗口
$`H'(\alpha,K)=T_\alpha+\max(N_B^\alpha,U_1^\alpha)+\max(N_A^\alpha+2,\,K)`$
上：若 $`N_A^\alpha<K`$，则冻结规则（尾段 $`m`$-连跑 $`\ge K`$ 判 LOCK）也完全正确。
当 $`N_A^\alpha\ge K`$ 时冻结规则**可能**在 CYCLE 真值上误判 LOCK（取决于窗口
终点是否落在长 A-游程深处）——不构成当且仅当，因为特定窗口终点也可能恰好
安全。

**对已登记反例的解释。** $`K=9,\ \alpha=\tfrac18`$：$`N_A^{1/8}(K=9)=18\ge9`$
（全 CYCLE 类的一致界为 $`\lceil\ln81/\ln\tfrac87\rceil=33`$），冻结规则在
$`H=90`$ 窗口输出 LOCK 而 $`A`$ 的无限时间结论为双支切换——正是 $`N_A^\alpha\ge K`$
预言的失效模式（`c_summary.json:counterexample` 精确复现：冻结=LOCK、两测=CYCLE）。

**LOCK 类的界约定。** $`r\in\{0,1,2\}`$ 时 $`L\ge h`$，A-游程永不结束，
$`N_A^\alpha,U_1^\alpha`$ 按一致解析界 $`\min\{n:(1-\alpha)^n\le\tfrac1{81}\}`$ 取值
（仅用于窗口公式的有限性；LOCK 侧的判定只需 $`T_\alpha`$ 与 $`N_B^\alpha`$）。

## 3. 证明

**LOCK 侧。** $`r\in\{0,1,2\}`$：$`L\ge h`$。进入后若 $`x_{t}>h`$，则 A-支
$`x'=L+(1-\alpha)(x-L)\ge\min(x,L)>h`$（$`L\ge h`$），从此恒选 $`m`$。进入若经 $`\ell`$
（$`x\le h`$），B-支以 $`(1-\alpha q)^n`$ 几何趋 $`c>h`$，至多 $`N_B^\alpha`$ 步穿越。
故 $`\tau_0`$ 之后永无 $`\ell`$-选择：$`\mathcal C_2`$ 输出 LOCK ✓；且尾段
$`m`$-连跑 $`\ge H-\tau_0\ge K`$，冻结规则亦正确。

**CYCLE 侧。** $`r\in\{3,4,5\}`$：$`L<h`$。首次 $`m`$-连跑始于 $`\tau_0`$ 之前
（$`t_*+\max(N_B^\alpha,U_1^\alpha)`$ 内必已出现 $`m`$-选择：或经 $`\ell`$ 进入后
$`\le N_B^\alpha`$ 步上穿，或直接经 $`m`$ 进入），该连跑内的 $`x`$-值
$`\le h+\alpha L`$（B-穿越过冲 $`\le\alpha(c-x)q\le\alpha L`$；经 $`m`$ 进入的首值
$`\le c`$，其首 A-游程 $`\le U_1^\alpha`$ 后着陆 $`\le h+\alpha L`$ 一类界），故
连跑长度 $`\le N_A^\alpha`$，于 $`\tau_0+N_A^\alpha<H`$ 内出现 $`\ell`$-选择：
$`\mathcal C_2`$ 输出 CYCLE ✓。

**冻结规则的充分条件。** 在窗口 $`H'`$ 上：LOCK 侧尾段
$`m`$-连跑 $`\ge H'-t_*-N_B^\alpha\ge K`$，LOCK 正确触发；CYCLE 侧当
$`N_A^\alpha<K`$ 时任何时刻的尾段 $`m`$-连跑 $`<N_A^\alpha<K`$，LOCK 不会触发，
且进入后选择集 $`=\{\ell,m\}`$（持久性），CYCLE 正确触发。当
$`N_A^\alpha\ge K`$ 时窗口终点可能落入长 A-游程使尾段 $`\ge K`$——误判可能但非
必然（反例给出其实例）。$`\square`$

## 4. 验证与范围

- `run_c.py`：$`(\alpha,K)`$ 网格（6 个 $`\alpha`$ 值 × 21 个 $`K`$）精确运行：
  驻留界全部成立（C1）、两测分类器全部正确（C2）、冻结分类器恰在
  $`N_A^\alpha<K`$ 处正确（C3）、LOCK 类 $`\tau_0`$ 后永无 $`\ell`$（C4）、反例精确
  复现。全部通过（`c_summary.json`）。
- 观测成本：$`H(\alpha,K)=\Theta(K/\alpha)`$——分类可靠窗口随 $`\alpha\to0`$
  线性增长；$`\alpha=\tfrac14`$ 时 $`H\le10K`$ 与 T6 一致（本定理不重复 T6 的全 $`K`$
  证明，T6 仍为 $`\alpha=\tfrac14`$ 的独立登记）。
- 范围：canonical 初态；$`\alpha\in(0,1)`$ 有理/无理皆可（界只依赖 $`\alpha`$ 的
  实值，实现验证用有理 $`\alpha`$）。不含非 canonical 初态（与 B 组合时窗口需
  以 B 的 $`\Delta_0`$ 界替换 $`T_\alpha`$）。
