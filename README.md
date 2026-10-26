# 沙子再分配研究 / Sand Redistribution Study

**公开研究快照 · Public research snapshot — [v0.4.0](https://github.com/XuZhilin2007/sand-redistribution-study/releases/tag/v0.4.0) · 2026-10-26**

本项目研究由“回收过密前缀、再撒开、适时停止”启发的一维确定性质量反馈模型。
v0.4.0 的主结果是：在 canonical XO、精确算术、最小指标平局规则下，
**全部整数 K≥6、整个参数区间 J、全部整数时间 t≥0，两次更新达到平方损失的全局最小值**：

```math
J=[1999/8000,2001/8000]=[0.249875,0.250125],\qquad
U(p)=K\sum_i(p_i-1/K)^2,\qquad U(p_t)\ge U(p_2).
```

α 遍历 J 内全部实数。区间与 K 无关；v0.3.0 的 α=1/4 定理是其特例。
最优时刻不要求唯一，后续损失不要求单调。完整模型、假设与量词见
[成果登记](docs/V040_REGISTRATION.md)。

For canonical XO, exact arithmetic and smallest-index ties, **two updates
minimize squared nonuniformity for every integer K≥6, every real α in J,
and every integer time t≥0**. The parameter interval is uniform in K.
This is a deterministic model result, not a calibrated physical claim.

## v0.4.0 的完整证明 / Proof and evidence

| 部分 | 范围与证据 |
|---|---|
| [前缀选择](docs/v04_research/LARGE_K_PREFIX_LEMMA_REVIEW.md) | 真实选择位置的舍入公式；全部舍入误差和平局由 η 参数化 |
| [D₂ 一致进入](docs/v04_research/D2_ENTRY_REVIEW.md) | K≥33、α∈J 的边界支配；小 K 由有限证书补齐 |
| [驻留与不变量](docs/XO_ALPHA_EXTENSION_PROOF.md) | 已发布定理 A；首次边界选择后永久驻留，不要求固定选择字 |
| [全 K W0](docs/v04_research/W0_UNIFORM_BOUND_REVIEW.md) | U₀>U₂、U₁>U₂；L¹ 收缩、格子方差与精确多项式正性 |
| [统一尾部](docs/v04_research/UNIFORM_TAIL_CLOSURE_REVIEW.md) | K≥204、整个 J、全部 t≥3：U_t−U₂>1/100；48 个精确有理闭盒 |
| [有限区间证书验收](tools/audit_v04_interval_v3.py) | 目标集 203 个 K、207 分片、80,246 条符号义务；覆盖 K=6..203 并保留五个锚点 |

K=203/204 无缝相接。无限参数结论由解析恒等式与精确计算辅助证明承担，
有限轨迹重放作为装配交叉检查。α 端点和平局分片均纳入闭合证明。

Analytic identities, exact rational interval enclosures and independently
audited finite certificates carry the universal quantifiers. Full-vector
replays check the implementation; sampled trajectories do not replace the proof.
Independent implementations and AI-assisted checks are not formal proof-assistant certification.

## 快速复验 / Quick verification

从克隆或下载源码的仓库根目录运行：

```text
python -B tools/verify_v040.py
```

v0.4.0 主定理复验仅使用 **Python 标准库**，无需额外依赖。测试环境为
Python 3.14；程序顺序检查固定清单、五个数学审计器和受影响的回归，
约一分钟。默认不改写固定输入或旧报告。不要使用 `-O`。

The v0.4.0 verification path uses only the Python standard library and is
self-contained in the public source tree. It checks file hashes before and
after all six checks. No private repository or local absolute path is needed.

- [研究验收记录](docs/V040_VALIDATION.md)与[公开发布验收](docs/V040_PUBLICATION.md)
- [公开文件清单](docs/V040_MANIFEST.json)与[原研究冻结清单](docs/V040_FROZEN_MANIFEST.json)
- [阶段成果台账](docs/v04_research/RESULTS.md)：其他相图、反例、扰动与指标探针按各自形成时的范围保留；历史待办不作为当前主定理缺口

## v0.3.0 继承成果 / Retained results

[v0.3.0 登记](docs/V030_REGISTRATION.md)中的长期选择字、渐近质量状态、
初态结果、一般 α 词典、平台和端点证书保留各自原范围。选择字的最终精确
周期性与质量状态的渐近周期性始终区分。文献继承定位见
[V030_RELATED_WORK](docs/V030_RELATED_WORK.md)。

历史版本的完整回归使用 [requirements.txt](requirements.txt)；大型 B/E
证书以经过校验的 gzip 分片保存。执行旧版本全量验证前先运行

```text
python -B tools/materialize_v030_certificates.py
python -B tests/count_tests.py
```

原始证书从 `certificates/v0.3.0/` 离线恢复并核对 SHA-256。v0.4.0 的
一分钟主定理复验不依赖这一步。历史声明以各版本登记为准。
