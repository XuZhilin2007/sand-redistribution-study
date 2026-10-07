# 沙子再分配研究 / Sand Redistribution Study

**公开研究快照 · Public research snapshot** — [v0.2.0](https://github.com/XuZhilin2007/sand-redistribution-study/releases/tag/v0.2.0) · 2026-10-07

本仓库研究一个由沙面修正启发的**一维确定性质量反馈模型**，保存代码、冻结实验和可核查证明。模型不是沙粒碰撞或真实沙坑仿真。v0.2.0 在原 canonical XO S3 成果上，正式收录 **A：操作力度推广**和 **T6：原冻结分类器可靠性**，另附限定辅助成果。问题、模型及结果见[双语研究报告](docs/RESEARCH_REPORT.md)。

This repository studies a **one-dimensional deterministic mass-feedback model** inspired by sand-surface correction. It contains code, frozen experiments, and checkable proofs, rather than a grain simulation. v0.2.0 adds two core supplements to canonical XO S3: **A, the removal-fraction extension**, and **T6, reliability of the original frozen classifier**, plus scoped auxiliary appendices. Start with the [bilingual report](docs/RESEARCH_REPORT.md).

## 当前结论 / Current results

共同前提为 canonical 初态、XO 固定喷撒律、最小指标平局、精确算术及整数 K≥6。

| 成果 / Result | 适用范围与结论 / Scope and conclusion |
|---|---|
| 原 S3 T1–T5 | α=1/4；有限进入、保持、边界坐标闭包与 K mod 6 的长期选择二分。原证明和审计范围保持原样。 / Original boundary theorem at α=1/4; its deposited proof and audit retain their scope. |
| [A](docs/XO_ALPHA_EXTENSION_PROOF.md) | 0<α<1；上述长期结构推广到任意该范围内的操作力度。 / The long-term boundary structure extends to every 0<α<1. |
| [T6](docs/XO_T6_CLASSIFIER_PROOF.md) | **仅 α=1/4、H=10K**；所有 K≥6，冻结分类器与长期余数类二分一致；证明依赖原 S3，不依赖 A。 / **Only α=1/4 and H=10K**: the unchanged classifier matches the dichotomy for every K≥6; the proof depends on S3, not A. |
| [限定附录](docs/XO_SUPPLEMENTARY_RESULTS.md) | 单调地板/条件守恒、首次选择公式、有限探针、径向聚合/半共轭和文献线索。 / Monotone floor, first-update formula, scoped probes, radial aggregation, and bibliographic leads. |

**A 与 T6 不能合并为任意 α 的有限窗口分类保证。** 例如 K=9、α=1/8、H=90，原分类器输出 LOCK，但长期两边界无限切换。K=10、α=1/1000、H=100 则可输出 OTHER。两例均进入回归。边界坐标闭包不恢复完整质量向量，双支无限选择不等于周期轨道；一般噪声、一般二维、其他初态/律和文献新颖性仍未解决。

**A does not extend T6's finite-window guarantee to arbitrary α.** At K=9, α=1/8, H=90, the original classifier returns LOCK despite infinite long-term switching. At K=10, α=1/1000, H=100, it returns OTHER before entry. Both are regression cases. Scalar closure does not reconstruct the full mass vector, and infinite switching does not prove periodicity. General noise/2D theory, other initial states/laws, and academic novelty remain open.

原 S3 为 **CLOSED — PASS**，其[审计来源](docs/S3_XO_BOUNDARY_THEOREM_AUDIT.md)保持不变。A/T6 有单独的[正式登记](docs/XO_EXTENSION_REGISTRATION.md)和[独立技术审查记录](docs/XO_EXTENSION_REVIEW.md)。这些记录不是期刊同行评审或形式化机器证明。

Original S3 remains **CLOSED — PASS** with its original [audit provenance](docs/S3_XO_BOUNDARY_THEOREM_AUDIT.md). A/T6 have separate [registration](docs/XO_EXTENSION_REGISTRATION.md) and [technical review provenance](docs/XO_EXTENSION_REVIEW.md). This is not a peer-reviewed paper or formal machine verification.

## 阅读路线 / Reading order

1. [双语研究报告](docs/RESEARCH_REPORT.md)与[中文导读](docs/XO_THEOREM_GUIDE_ZH.md)：问题、模型及当前范围。
2. [补充登记](docs/XO_EXTENSION_REGISTRATION.md)、[A](docs/XO_ALPHA_EXTENSION_PROOF.md)、[T6](docs/XO_T6_CLASSIFIER_PROOF.md)：两项核心陈述及各自依赖。
3. [原 S3 英文证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)、[中文完整译本](docs/S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md)与[原登记](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)：T6 的基础与原定理范围。
4. [补充验收](docs/XO_EXTENSION_VALIDATION.md)、[辅助附录](docs/XO_SUPPLEMENTARY_RESULTS.md)及[阅读说明](docs/READING_NOTES.md)：复现、证据等级与历史材料的读法。
5. [S1](docs/S4_W1_S1_XO_PRECISION.md)、[S2](docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md)、[Stage 2](docs/RESEARCH_CHECKPOINT_2026-09-20.md)和[Stage 3](docs/RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md)：保留形成时的范围；“尚未证明”由后续已登记成果补充，不是当前结论。

The current route is the bilingual report, separate A/T6 registration and proofs, original S3 sources, and the validation/appendix notes. Dated S1/S2 reports and checkpoints remain historical evidence. A/T6 supersede their earlier open questions only within the newly registered scopes. No manuscript, general-noise, or general-2D research was started for this release.

## 复现 / Reproduction

补充验证只使用 Python 标准库，在仓库根目录运行；输出目录必须为新目录或空目录：

```text
python -B tests/test_xo_extensions.py
python -B experiments/xo_extensions/validate.py --scope full --out-dir ../xo-v020-run
```

入口从本仓原 S2 文件读取 classify_itinerary，不需要实验室、base/、私人仓库或旧审查记录。配置描述精确 α、K、实际窗口与算术口径，输出 manifest.json 记录源码 SHA-256。完整矩阵是独立长入口；日常精简回归自动接入 tests/count_tests.py。

The supplements run with the standard library alone. The entry loads this repository's frozen S2 classifier, requires no lab or private checkout, and records configurations, actual horizons, arithmetic, and source hashes. The full matrices are separate from the fast routine regression.

既有模型测试使用 **Python 3.14.6** 和 [requirements.txt](requirements.txt) 的固定依赖。Windows PowerShell：

```text
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe -B tests\count_tests.py
```

macOS/Linux 使用 `.venv/bin/python`。原 v0.1.1 为 20 套、382 项运行时断言；本版新增独立补充套件。实际整合验收数量与脱离实验室复现记录见[验收页](docs/XO_EXTENSION_VALIDATION.md)。不要用 `-O` 关闭科学断言。

The existing suites use the pinned dependencies; on macOS/Linux use `.venv/bin/python`. v0.1.1 had 20 suites and 382 runtime assertions. This version adds a supplement suite; see the validation page for the actual integrated counts and isolated reproduction. Do not use Python's `-O` option.

## 许可与版本 / Licensing and versioning

软件、测试、实验脚本和配置采用 [Apache-2.0](LICENSE)；原创研究文稿、图表及结果采用 CC BY 4.0，具体边界见[内容许可](CONTENT_LICENSE.md)。作者署名为 Zhilin Xu；[CITATION.cff](CITATION.cff) 的项目 version 为 0.2.0，cff-version 1.2.0 是引用格式版本。

Software and configs use Apache-2.0; original research prose, figures, and results use CC BY 4.0 under the [license map](CONTENT_LICENSE.md). Cite Zhilin Xu and the project version in CITATION.cff. Its cff-version denotes the citation format, not the project version.

本仓库保持独立公开历史，不含私人 Git 历史、实验室 base/、旧错误稿、before/ 或本机管理记录。固定引用使用 [v0.2.0](https://github.com/XuZhilin2007/sand-redistribution-study/tree/v0.2.0)；旧 [v0.1.1](https://github.com/XuZhilin2007/sand-redistribution-study/tree/v0.1.1) 保留。变更见 [RELEASE_NOTES.md](RELEASE_NOTES.md)。

This independent public history excludes private Git history, lab base/, superseded drafts, and local management records. Cite the fixed v0.2.0 tag; v0.1.1 remains preserved. See the release notes for changes.
