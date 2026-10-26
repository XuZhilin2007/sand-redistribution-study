# v0.4.0 研究状态入口（STATUS）—— 主线闭合，本地成果冻结

> 主证明闭环：2026-10-23；成果与复验收敛：2026-10-24。
> 正式范围：[V040_REGISTRATION](../V040_REGISTRATION.md)；
> 实际验收：[V040_VALIDATION](../V040_VALIDATION.md)；
> 固定清单：[V040_MANIFEST](../V040_MANIFEST.json)。
> 分支：`codex/v04-research`（自 `main` @ `41219f0` 建立）。

## 主定理（本轮闭环，全部读证承接）

**canonical XO、精确算术、最小指标平局规则**：

```math
\forall K\in\mathbb Z_{\ge6}\quad
\forall\alpha\in J=[1999/8000,2001/8000]\quad
\forall t\in\mathbb Z_{\ge0}:\quad
U(p_t)\ge U(p_2),\qquad U(p)=K\textstyle\sum_i(p_i-1/K)^2 .
```

v0.3.0 的 α=1/4 全 K 定理是本定理在 α=1/4 处的推论；本定理给出
与 K 无关的双侧参数窗口。不涵盖一般初态、一般喷撒律或其他损失。

**证明链组装（无缺口）**：

| 参数/时间范围 | 承载证明 | 复验 |
|---|---|---|
| K=6..203、α∈J、全部 t≥0 | 有限区间证书（含 4 组真实平局分片、各分片真实 U₂），经**独立 Bernstein 数学验收**（203/203 K、207 分片、80,246 义务；SHA c8c7acae…） | `tools/audit_v04_interval_v3.py`（~21s） |
| K≥204、α∈J、t=0,1 | **全 K W0**（L¹ 收缩+格子方差；常数 11/(10K)、5/(2K)；极限余量 9/50、1/40；K≥145 拼接） | `tools/audit_v04_w0_uniform.py`（~4s） |
| K≥204、α∈J、t=2 | 参考值恒等 | — |
| K≥204、α∈J、t≥3 | **统一尾部**（48 精确有理闭盒：H′(ρ)<−3/25、H(ρ)>1/100；前置=前缀引理+D₂ 进入+定理 A 驻留） | `tools/audit_v04_uniform_tail.py`（~2.3s） |
| 进入（K≥33，α∈J） | **D₂ 一致进入**（μ∈[37/100,2/5]、中段<3/25、D₂(ℓ)≥247/1980、余量 47/9900；采用短证明） | `tools/audit_v04_d2_entry.py` |

K=203/204 无缝相接；α 端点含于闭盒与有限证书；s₁ 平局由真实 η
参数化覆盖；尾部不要求固定 s₂ 或未来字；二次型交叉项 2z(dx+f)。

**证据等级**：主定理 = 解析恒等式 + 精确计算辅助证明（闭盒/有限
证书）。推导、实现与交叉核验由 Zcode（GLM5.3）和 Codex 分别参与；
Owner 组织研究、设定范围与复核流程。独立实现检查不等于形式化证明。

## 快速数学复验路径（按序，全部真实退出码）

```powershell
$env:PYTHONUTF8="1"
python -B tools/audit_v04_interval_v3.py --scope common
python -B tools/audit_v04_w0_uniform.py
python -B tools/audit_v04_d2_entry.py
python -B tools/audit_v04_uniform_tail.py
python -B experiments/v04_w_interval/negative_tests_w.py
```

默认复验只输出结果。需要保存本次报告时，审计器支持 `--report` 指向
新的运行目录；不覆盖历史报告或原见证。完整命令与验收记录见 V040_VALIDATION。

## 停止状态

全 K、整个 J、全部整数时间的证明与有限补集均已闭合。当前范围没有
继续扩扫或补统一尾部的待办。其他方向按原证据范围保留在历史成果台账，
不作为本版本的未完成项。

## 历史状态（按日期，已全部被上表承接或按证据范围归档）

- 2026-10-22 第一/二波 + R1–R7 修复：相图、B-D 分解、C/D/E 探针、
  逐 K 区间证书 v3（I1/I2 修复、tie-split、Gram）。勘误台账见
  `RESULTS.md`。
- 2026-10-23 第三/四/五轮：弱非负内核重建、前缀引理承接、双分支
  连续极限（α*=16/63 双子序列修正）、W0 承接、D₂ 定理、工程门。
  曾计划的"聚合矩大模块、机械扩扫 204..2048、窗口扩张"**已被
  统一尾部证明取代**（2026-10-23 起为历史状态，不再待办）。
- 其余方向（最大窗口/完整相图、一般初态、多指标、EM 锐化、新
  无限平台分类、F/I）按各自证据范围归档于 `RESULTS.md` 与各方向
  文档，本轮不扩张。

## 来源（冻结提交所选文件与 Owner 复核链）

证明：`docs/v04_research/{W0_UNIFORM_BOUND_REVIEW,UNIFORM_TAIL_CLOSURE_REVIEW,D2_ENTRY_REVIEW,LARGE_K_PREFIX_LEMMA_REVIEW}.md`
核验：`tools/audit_v04_*.py`；报告：
`experiments/*/review_2026-10-23/*.json`；原始扫描与 v3 见证按原
SHA 保留。生成代码（interval_certs/symkernel/w0_sweep/d2_entry/
continuum_limit）与验证代码分离；具体哈希见 manifest。

## 历史成果速览（细节与勘误台账见 RESULTS.md 与各方向文档）

- 第一/二波 + R1–R7 修复包：25 锚点、相图扫描、边缘二分、
  B-D 分解/B-S 证书、C1/C2/D/E 探针（撤回项见台账）。
- 阶段二 v3：203/203 目标 K 逐 K 全时间区间证书（tie-split 含
  平局族；I1/I2 修复 + Gram 恒等式）；公共区间 [0.249875,0.250125]
  （后经外部独立验收链确认）。
- 精确反例族：K=9 α=1/5；K=6 α=7/24；K=9 α=1/3；K=9 P3 δ=1/10；
  K=9 倾斜 δ=1/10（详见 RESULTS.md CX 表）。

## 提交记录

41219f0（v0.3.0 main）→ …（各轮见 git log）→ 36d71a2
（上一复核基线）→ 本轮冻结提交（哈希见最终交付/manifest）。
