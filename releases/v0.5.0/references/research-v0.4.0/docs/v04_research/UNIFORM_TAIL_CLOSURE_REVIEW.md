# 统一尾部的精确闭环：K≥204、整个 J、全部 t≥3

日期：2026-10-23。研究基线：`36d71a2`。复核者提供本证明，Zcode 已独立读证、回代并承接。最终组装与复验范围见 `../V040_REGISTRATION.md`。

**结论。** 对 canonical XO、精确算术、最小指标平局规则，令

```math
J=[1999/8000,2001/8000].
```

对全部整数 K≥204、全部 α∈J、全部整数 t≥3，有

```math
U(p_t)-U(p_2)>1/100.
```

证明依赖已复核的真实前缀选择引理、D₂ 进入引理及定理 A 的驻留，**没有要求以后固定同一个边界分支**。只用 48 个精确有理闭盒，覆盖连续参数域和无限 K。

结合全 K W0、t=2 恒等式及现有 K=6..203 的全时间区间证书，即得

```math
\forall K\ge6\quad\forall\alpha\in J\quad\forall t\in\mathbb Z_{\ge0}:\qquad
U(p_t)\ge U(p_2).
```

因此无需新增 K=204..2048 的逐 K 证书，无需扩大 J、优化误差常数或继续探索新窗口即可把当前板块闭合。

## 1. 真实两步状态与缩放

记 h=1/K、r=K mod 6、ρ=1−α、m=⌈5K/6⌉、ℓ=m−1，并置

```math
\delta=r/6,\quad\theta=(6-r)/3,\quad
b=m/K=5/6+\delta h,\quad L=\ell/K=b-h,\quad q=g_m=\theta h.
```

以 k(y) 表示 XO 的密度：在 [0,1/3]、[1/3,2/3]、[2/3,5/6] 分别为 1/2、3/2、2，之后为 0。记 f₀(y)=2−2y、P(y)=2y−y²。

由前缀引理，真实选择位置可写成

```math
u_0=1/2-(r\bmod2)h/2,\quad
\lambda=\alpha P(u_0),\quad
u_1=1/2+3\lambda/4+\eta h,\quad -1/2\le\eta\le1/2.
```

于是

```math
\nu=\alpha[P(u_1)-\lambda(4/3-3u_1/2)],\quad
\mu=\lambda+\nu,\quad\beta=\rho\lambda+\nu.
```

两步状态 p₂ 恰等于下列分段密度 f₂ 的格子积分，cutoff 都在格子边界：

| y 的范围 | f₂(y) |
|---|---|
| 0≤y≤u₀ | ρ² f₀(y)+β k(y) |
| u₀<y≤u₁ | ρ f₀(y)+β k(y) |
| u₁<y | f₀(y)+μ k(y) |

对 K≥204，两个 XO 跳点格分别位于 u₀ 之前、u₁ 之后。脚本也在每个盒上检查了强于所需的 u₀>1/3+h、u₁<2/3−h，因此不会误装跳点格所属分段。

令 c=P(b)、x₂=p₂,m。因为 m 未被前两次选择直接削减，

```math
x_2=hX,\qquad X=2-2b+h+\mu\theta,
\qquad\gamma=(c-x_2)/(1-q).
```

对 i<m 定义 eᵢ=p₂,i−γgᵢ。前 m 格总质量为 c，且 Σᵢ₍<m₎gᵢ=1−q，故

```math
E_1:=\sum_{i<m}e_i=(c-x_2)-\gamma(1-q)=0.
```

以下三个量均为 O(1)：

```math
C_g=K\sum_{i<m}g_i^2,\qquad
Z=K\sum_{i<m}g_ie_i,\qquad E=K\sum_{i<m}e_i^2.
```

## 2. 格子方差与协方差给出精确矩式

这是恒等式，不是渐近估计。令残差密度 v(y)=f₂(y)−γk(y)。在没有 XO 跳点的一段 [a,b]，写 v(y)=A f₀(y)+C k，k 为常数。设 d=b−a，则

```math
\overline v=A(2-a-b)+Ck,\quad
\int_a^b k v=d k\overline v,\quad
\int_a^b v^2=d\left(\overline v^2+A^2d^2/3\right).
```

把 [0,L] 分成五段，依次使用

| [a,b] | A | C | k |
|---|---|---|---|
| [0,1/3] | ρ² | β−γ | 1/2 |
| [1/3,u₀] | ρ² | β−γ | 3/2 |
| [u₀,u₁] | ρ | β−γ | 3/2 |
| [u₁,2/3] | 1 | μ−γ | 3/2 |
| [2/3,L] | 1 | μ−γ | 2 |

将上述积分分别求和，记为 I_Z=∫₀ᴸ kv、I_E=∫₀ᴸ v²。

两个内部跳点在所在格中的左侧比例为

```math
t_1=(r\bmod3)/3,\quad t_2=(2r\bmod3)/3,\quad v_j=t_j(1-t_j).
```

跳点恰在格边界时 tⱼ=0，修正项为零。对一个宽 h、跳跃 Δ、左侧比例 t 的格子，直接积分给出

```math
\int\operatorname{Var}_{\rm cell}(k)=h t(1-t)\Delta^2,
\qquad
\int\operatorname{Cov}_{\rm cell}(k,f_0)=-h^2t(1-t)\Delta.
```

此处第一式的左端表示 ∫(k−平均值)²，第二式表示 ∫(k−平均值)(f₀−平均值)。因此

```math
C_g=3/2-2\theta h-h(v_1+v_2/4),
```

```math
Z=I_Z+h^2(v_1\rho^2+v_2/2)
-h[v_1(\beta-\gamma)+v_2(\mu-\gamma)/4],
```

```math
V_{\rm lin}=\frac{h^2}{3}[\rho^4u_0+\rho^2(u_1-u_0)+L-u_1],
```

```math
V_{\rm jump}=h[v_1(\beta-\gamma)^2+v_2(\mu-\gamma)^2/4],
```

```math
V_{\rm cross}=-2h^2[v_1\rho^2(\beta-\gamma)+v_2(\mu-\gamma)/2],
\qquad E=I_E-V_{\rm lin}-V_{\rm jump}-V_{\rm cross}.
```

V_lin 是每个整格线性方差 A²h³/3 的总和；cutoff 位于格边界，故没有遗漏 cutoff 格修正。V_jump 与 V_cross 分别是两个 XO 跳点的平方项、与线性部分的交叉项。最后的 5/6 跳点位于被排除的第 m 格，不在 i<m 的矩中。

只需在真实 η 上这些等式与离散矩精确相同。盒里其他 η 可以不对应格点轨道；对恒等式的有理函数延拓作外包络仍覆盖全部真实参数。

## 3. 消去大数相减的尾部表达

从 t=2 起已边界驻留。设 x=x_t、z=ρ^(t−2)，有

```math
p_{t,i}=\frac{c-x_t}{1-q}g_i+z e_i\quad(i<m),
```

并且 i>m 的后缀冻结。该式对选择 m 和 ℓ 的两支分别归纳：γ′ 分别为 ργ+αc、(1−αq)γ；残差两支都乘 ρ。因此允许任意真实边界选择字。

展开 U 得

```math
U(p_t)=Ax^2+Bx+C+2z(dx+f)+Ez^2,
```

其中

```math
\widehat A=hA=1+\frac{hC_g}{(1-q)^2},\quad
B=-\frac{2cC_g}{(1-q)^2},\quad
d=-\frac Z{1-q},\quad f=\frac{cZ}{1-q}.
```

E₁=0 已用于 B、f 的简化；A>0。配方定义

```math
U(p_t)-U(p_2)=A\left(x+\frac{B+2dz}{2A}\right)^2+H(z).
```

Gram 不等式 Σg²Σe²≥(Σge)² 给出

```math
H_2=E-d^2/A\ge0.
```

利用 t=2 的精确恒等式消去 C−U₂，避免区间算术中大项相减，所需两式为

```math
H'(\rho)=2\left(E-\frac{d^2h}{\widehat A}\right)\rho
+2f-\frac{Bdh}{\widehat A},
```

```math
H(\rho)=E(\rho^2-1)+2f(\rho-1)
-hX(\widehat A X+B+2d)
-\frac{h(B+2d\rho)^2}{4\widehat A}.
```

这些式子在 h=0 处也有无奇点的连续延拓；计算不除以 h。

## 4. 48 个闭盒覆盖无限参数域

固定 r=0,…,5，取 α∈J，并用

```math
h\in[j/816,(j+1)/816]\quad(j=0,1,2,3),\qquad
\eta\in[-1/2,0]\ \text{或}\ [0,1/2].
```

这些 6×4×2=48 个闭盒的并集包含所有 K≥204 的真实参数。ρ、λ、ν 等派生量均从同一盒的 α、h、η 按上式计算；自然区间扩张只使包络变宽，不影响有效性。分母 1−q≥101/102，Â≥1；程序对每次除法拒绝包含零的分母。

标准有理区间四则运算、带零检查的平方，逐盒验证严格有理不等式

```math
H'(\rho)<-3/25,\qquad H(\rho)>1/100.
```

全部精确上下界保存在 `uniform_tail_independent_audit.json`，不使用浮点阈值或采样填补盒内。浮点近似只用于阅读：最大导数上界约 −0.12305，最小尾部值下界约 0.01167。

由于 H₂≥0，H′ 在 z 上递增；z∈[0,ρ] 时 H′(z)≤H′(ρ)<0，故 H(z)≥H(ρ)>1/100。对整数 t≥3，z=ρ^(t−2)∈[0,ρ]，平方项非负，所需尾部命题成立。

## 5. 全 K 组装与复验范围

| 参数/时间范围 | 承载证明 |
|---|---|
| K=6..203、α∈J、全部 t≥0 | 原 v3 见证，经独立 Bernstein 数学验收；不是 v4 版本字符串继承 |
| K≥204、α∈J、t=0,1 | `W0_UNIFORM_BOUND_REVIEW.md` 的全 K W0 |
| K≥204、α∈J、t=2 | 参考值恒等 |
| K≥204、α∈J、全部 t≥3 | 本证明；进入用 `D2_ENTRY_REVIEW.md`，驻留用定理 A |

两端 K=203/204 相接，没有遗漏、没有新增有限补集。端点 α 包含在闭盒和有限证书中；s₁ 平局均被真实 η 参数化覆盖。原见证 SHA-256：`c8c7acae690f939397b88084784df4162264ccef54e5f93c32b8d5e7c979b9d1`。

复验 `tools/audit_v04_uniform_tail.py`：标准库 Fraction 完成 48 盒证明；另由自含 CDF、全指标 argmax 和真实更新重放 177 个 (K,α) 点，逐位核对 Cg、Z、E、γ、x₂、A、B、d、f、H′、H，并在 t=2,3,4 回代二次型。有限重放负责发现矩装配错误；无限参数量词来自本节前的解析恒等式与闭盒外包络。

原生产引擎、v3 JSONL 与扫描数据均未修改。本结果为解析恒等式加精确计算辅助证明，范围为上述 canonical 模型与 J，不扩展至 J 以外、一般初态或一般喷撒律。
