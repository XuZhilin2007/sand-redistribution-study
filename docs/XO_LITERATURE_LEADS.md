# S8：文献线索清理与范围说明

> v0.2.0 正式收录，2026-10-07。范围与证据等级见[补充成果登记](XO_EXTENSION_REGISTRATION.md)，
> 独立审查来源见[补充审查记录](XO_EXTENSION_REVIEW.md)。实验室旧稿与修订历史在私人仓库保留。


## 1. 可保留的框架与书目

模型在固定选择区域内为线性/仿射更新，选择由累计偏差比较决定。
归约后的标量映射在 h 处有跳跃：
f(h⁻)−f(h⁺)=αh(1−q)>0。
因此“不连续的分段收缩映射”比笼统的连续 border-collision normal form 更具体。
这是本模型结构的判断，不意味着已找到直接适用的已有定理。

- di Bernardo、Budd、Champneys、Kowalczyk，
  *Piecewise-smooth Dynamical Systems: Theory and Applications*，2008，
  Applied Mathematical Sciences 163。
  出版方：https://link.springer.com/book/10.1007/978-1-84628-708-4
- Ethan Akin，*The General Topology of Dynamical Systems*，1993，
  主题为拓扑动力学。作者资料：https://math.sci.ccny.cuny.edu/person/ethan-akin/
  可以纠正领域误配，不能据此推断 reviewer 原本想引用哪本书。
- Granados、Alsedà、Krupa，
  *The Period Adding and Incrementing Bifurcations: From Rotation Theory to Applications*，
  SIAM Review 59(2),225–292,2017，doi:10.1137/140996598。
  作者稿：https://arxiv.org/abs/1407.1895
  内容涉及不连续映射、符号序列、旋转数和 period-adding；
  对本模型的适用性需要另核对其具体假设。

## 2. 新颖性边界

原关键词检索未找到直接对应论文，只能保留为“检索未命中”。
撤销“无直接先验工作”“具体结果不属于任何标准推论”
以及未逐定理定位的“该机制已确定属于教科书内容”等强措辞。
本轮既未确认原创，也未判定已有文献已经覆盖本模型。
可保留的研究价值是严格归约、参数范围和实验判定的可解释性；
文献新颖性需要后续专门的内容对照。

原检索中的其他框架与书目仅作 archived leads，不进入本轮已核实引用清单。
本轮未开展新的查新、引文追踪或最终 related-work 写作。

## 仓内复现与依赖

完整验证入口与实际配置、结果、算术类型及源码校验值见[发布候选验收](XO_EXTENSION_VALIDATION.md)。原 [S3 证明](S3_XO_BOUNDARY_THEOREM_PROOF.md) 的陈述和历史范围保持原样。
