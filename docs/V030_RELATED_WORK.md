# v0.3.0 文献与继承说明 / Related work and inherited results

日期：2026-10-21。本文给出发布版的结论边界与原文入口，不转载论文全文。
本项目尚未确立全面的文献新颖性，也不把程序验收等同于同行评审。

## 1. 已有理论与本模型的接口

Laurent–Nogueira 研究含两种斜率的区间映射。对本模型的 CYCLE 类，
令 $`h=1/K`$、$`q=g_m`$、$`L=cq`$，并定义
$`A(x)=(1-\alpha)x+\alpha L`$、
$`B(x)=(1-\alpha q)x+\alpha L`$。区间
$`U=(A(h),B(h)]`$ 一步不变。在
$`z=(B(h)-x)/(B(h)-A(h))\in[0,1)`$ 下，词典为

```math
\lambda=1-\alpha,\qquad
\mu=\frac{1-\alpha q}{1-\alpha},\qquad
\delta=1-(1-\alpha)\eta,\qquad
\eta=\frac{q(c-h)}{(1-q)h}.
```

因此 $`\lambda\mu=1-\alpha q<1`$。当 $`\alpha`$ 为代数数时，词典参数也是
代数数，可使用 Laurent–Nogueira 的有理旋转数结论。canonical
$`\alpha=1/4`$ 的全部 CYCLE 格点还由本项目证书排除了相关平台右端点，
从而可对接真正周期轨道及其吸引结论。完整状态为渐近周期；
选择字的最终精确周期性由 T7 的字许可与区间证书明确证明。

一般 $`\alpha`$ 下应区分平台内部与端点。G-ep 给出 K=10 的代数端点
$`\alpha_*`$，在那里原右连续映射没有真正周期轨道。不能把该端点的有限
$`\omega`$-极限集表述为原映射的实际周期轨道，也不能据此直接声称选择字最终周期。

The CYCLE reduction is conjugate to the two-slope interval family studied by
Laurent and Nogueira. Algebraic parameters imply a rational rotation number.
The attracting periodic-orbit statement also needs the appropriate non-endpoint
case. The canonical lattice is handled by this project's certificates; a
noncanonical algebraic removal fraction realizes an endpoint exception.

## 2. 继承与项目具体成果

| 内容 | 发布版定位 |
|---|---|
| 代数参数的有理旋转数、有理平台端点公式、适用情形下的周期轨道吸引 | 已有理论，经明确词典与假设核对使用 |
| 质量向量到边界标量、再到标准区间映射的归约 | 本模型的证明资产 |
| canonical K 格点上的显式字表、换挡阈值、平台与端点排除 | 本模型的显式分类与计算辅助证明 |
| 区间族覆盖、全称多项式证书、量化界与双实现核验 | 本项目的可复核证据资产 |
| 全 K 两轮最优停止、初态稳健性、LOCK 侧结构 | 独立陈述的模型结果；不由旋转数结论直接给出 |
| 一般 α 端点表、混合区密度解释、扰动稳健性 | 按各自有限验证或开放问题范围陈述 |

这些是对本项目成果的定位，不是“此前没有任何相似定理”的声明。
软件正确性、模型内数学结论、文献新颖性与真实沙坑适用性分别评估。

## 3. 原文入口

- [Michel Laurent and Arnaldo Nogueira, Dynamics of 2-interval piecewise affine
  maps and Hecke-Mahler series (arXiv:1907.08655)](https://arxiv.org/abs/1907.08655)：
  本项目使用的双斜率族、代数参数结果及有理平台端点情形。
- [José Pedro Gaivão, Michel Laurent and Arnaldo Nogueira, Rotation number of
  2-interval piecewise affine maps (arXiv:2212.10944)](https://arxiv.org/abs/2212.10944)：
  更一般参数族与阶梯结构；本项目一般 α 对接的背景。
- [Svante Janson and Anders Öberg, A piecewise contractive dynamical system
  and election methods (arXiv:1709.06398)](https://arxiv.org/abs/1709.06398)：
  分段收缩、旋转数与符号结构的相关模型。

原论文的范围以其正文为准；项目证明和有限计算的范围以
[正式登记](V030_REGISTRATION.md)及各证明文档为准。
