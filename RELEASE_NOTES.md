# 公开版本记录 / Public release history

## v0.2.0 — A/T6 核心补充与限定附录 / A/T6 supplements and scoped appendices

> 2026-10-07。原 S3、S1/S2 冻结材料与 v0.1.1 标签保留；本版分别正式收录 A 与 T6。 / The original S3 and frozen S1/S2 record and v0.1.1 tag are preserved. A and T6 are registered separately.

- **A：** canonical 初态、XO、最小指标平局、精确算术、K≥6、0<α<1 下的长期边界结构推广。 / Long-term boundary extension to every 0<α<1 under the stated canonical XO assumptions.
- **T6：** 仅 α=1/4、H=10K；原冻结分类器对所有 K≥6 与长期二分一致，依赖原 S3，不依赖 A。 / Reliability of the unchanged classifier for every K≥6 at α=1/4 and H=10K; depends on S3, not A.
- 新增[独立补充登记](docs/XO_EXTENSION_REGISTRATION.md)与[审查来源](docs/XO_EXTENSION_REVIEW.md)，更新 README、双语报告、阅读说明和中文导读；原 S3 审计不扩写。 / Separate registration/review provenance and updated current reading routes; the original audit remains unchanged.
- 收录 F、首次选择公式、限定探针、径向聚合/半共轭和文献线索；S12 明确为标量 float 代理，C1 补齐参数和按被扫次数衰减。 / Scoped auxiliaries, explicitly labelled scalar float statistics, and parameterized radial aggregation.
- 验收迁移为仓内入口，从本仓 S2 文件读取冻结函数，记录配置、实际窗口、算术和源码 SHA-256；新增日常精简回归，两条任意 α 的有限窗反例保留。 / Repository-local validation, frozen configs and hashes, and routine exact regression with both scope counterexamples.
- 实际测试数量、完整矩阵与脱离实验室复现见[候选验收](docs/XO_EXTENSION_VALIDATION.md)。 / See integrated validation for actual suite counts, full matrices, and isolated reproduction.
- 修正 S3 中英文和导读的 GitHub 公式转义显示；等价 TeX 宏与换行调整不改变数学内容，18 个编号公式仍对应。 / Corrected GitHub display escaping with equivalent TeX macros and line breaks; all 18 tagged bilingual formulas still match.
- 项目引用版本为 0.2.0；cff-version 1.2.0 格式字段保持原样。 / Project citation version is 0.2.0; the CFF format version is unchanged.

本版不宣称任意 α 的 H=10K 分类保证、周期性、全状态渐近、一般噪声/二维理论或学术新颖性。实验室 base/、旧错误稿、私人 Git 历史和本机管理记录不进入公开版。 / No arbitrary-α finite-window guarantee, periodicity, full-state asymptotics, general noise/2D theory, or novelty claim is added. Private histories and lab management records are excluded.



## v0.1.1 — 文稿与引用维护 / Exposition and citation maintenance

> Released 2026-10-06. The canonical XO S3 theorem remains CLOSED — PASS within its existing scope. This patch adds no experiment, new theorem, or extension of the model.

### 更正 / Corrections

- M1C5 §0 第 3 条及 §5：将“h 恒为零”改为“h 非正、连续最大值为零”；§4 端点说明改用直接恒等式，原必要门槛结论不变。 / Corrected nonpositivity and the continuous maximum in M1C5, and replaced an inaccurate endpoint argument with a direct identity; the necessary-gate conclusion is unchanged.
- 新增[中英双语阅读说明](docs/READING_NOTES.md)，说明证据名称、容差停机、有限时限 LOCK/CYCLE、局部记号及历史更正。 / Added bilingual reading notes for evidence names, tolerance stopping, finite-horizon labels, local notation, and errata.
- S1 的 `m=m=` 笔误和测试说明中的旧文件名已修复；测试运行逻辑不变。 / Fixed the S1 notation typo and an obsolete path in a test docstring; test execution is unchanged.
- 定理登记明确区分私人来源提交与公开独立历史。公开证明自含，旧日志和私人评审 Pack 不是必要依赖。 / Clarified private-source provenance versus public history; the deposited proof is self-contained and does not require private logs or the review pack.
- 引用 metadata 更新为 v0.1.1、2026-10-06，署名 Zhilin Xu；内容许可说明补充项目版权署名，官方许可全文保持原样。 / Updated citation metadata and project attribution; official license texts are preserved.
- README 增加阅读说明、当前维护状态及固定版本入口；S8、S9 和其他研究分支当前暂不启动。 / Updated the README's reading route, maintenance status, and fixed-version link; S8, S9, and other research branches remain deferred.

### 验证 / Validation

- 20 套既有测试、382 项运行时断言全部通过（2026-10-06）；未新增测试或运行扩展实验。 / All 20 existing test suites and 382 runtime assertions passed; no new tests or expanded experiments were added.
- 相对 Markdown 链接目标缺失、排除目录文件、定向敏感信息形态检查命中均为 0。 / Zero missing relative Markdown targets, excluded-directory files, or targeted sensitive-pattern hits.
- 289 个软件、冻结配置及结果文件与本轮私人源文件逐字节一致；相对上一公开提交，`src/`、所有实验文件和冻结配置无变化，仅测试文件开头的说明路径更正。 / The 289 software, frozen-config, and result files match the updated private source byte for byte. Relative to the preceding public commit, source and all experiment files are unchanged; only a test's leading docstring path was corrected.
- S3 英文证明和中文译本均未改动；18 个带编号公式在剔除空白差异后完全一致。 / Both S3 proof documents are unchanged; all 18 tagged displays match after whitespace normalization.
- S3 定理范围和证据等级保持原样；LOCK/CYCLE 仍为有限时限分类标签，没有新增周期性或全状态约化主张。 / The theorem scope and evidence levels are unchanged; LOCK/CYCLE remain finite-horizon labels, with no added periodicity or full-state reduction claim.

版本与引用见 [CITATION.cff](CITATION.cff) 及 [v0.1.1 固定快照](https://github.com/XuZhilin2007/sand-redistribution-study/tree/v0.1.1)。 / See the citation metadata and the fixed v0.1.1 snapshot.

## v0.1.0 — 初始公开快照 / Initial public snapshot

> Prepared 2026-10-01, after the canonical XO S3 proof deposit. This is a curated research snapshot with an independent Git history, not a copy of the private research history.

### 范围 / Scope

本版收录模型与实验代码、测试、17 个当前实验目录的冻结配置和已提交输出、主要阶段性研究文档，以及新的中英双语[公开研究报告](docs/RESEARCH_REPORT.md)、[中文 XO 定理导读](docs/XO_THEOREM_GUIDE_ZH.md)和[完整中文证明译本](docs/S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md)。审计来源仍为[英文正式证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)；旧评审 Pack 的原 XO 不变量证明已完整重录于该证明 §2，因此公开版不依赖私人草稿。

The release includes the model and experiment code, tests, frozen configurations and committed outputs from 17 current experiment directories, principal scientific stage records, a bilingual [research report](docs/RESEARCH_REPORT.md), a [Chinese XO theorem guide](docs/XO_THEOREM_GUIDE_ZH.md), and a [full Chinese proof translation](docs/S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md). The [formal English S3 proof](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md) remains the audited theorem-level source. Its §2 reproduces the earlier XO invariant proof, removing any dependency on an unpublished review draft.

私人完整提交历史、内部评审草稿与导出件、原始聊天运行材料、bootstrap 归档和研究治理日志没有收入此版。日期较早的科研报告按形成时的证据层级保留；当前 XO 定理状态由本版 README、双语报告和正式证明说明。 / The private commit history, internal review drafts and exports, raw chat-run materials, bootstrap archive, and governance log are excluded. Dated scientific reports preserve their contemporary wording; use the README, bilingual report, and formal proof for the current XO status.

### 验证 / Validation

- 新建 Python 3.14.6 虚拟环境，安装 [requirements.txt](requirements.txt) 中的固定依赖：成功。 / Fresh environment and pinned dependency installation: passed.
- `python tests/count_tests.py`：20 suites，382 runtime assertions，全部通过。 / 20 suites and 382 runtime assertions passed.
- 相对 Markdown 链接目标缺失：0。 / Missing relative Markdown link targets: 0.
- 排除目录（`external_review`、`raw_chat_runs`、`archive`）的文件：0。 / Files from excluded private directories: 0.
- 对当前快照做凭据形态、私人作者邮箱域名及个人 Windows 路径的模式检查：0 命中。这是定向检查，不是对所有可能敏感信息的保证。 / Targeted credential, private-email-domain, and personal-path scan: 0 hits; this is not a guarantee against every possible sensitive item.
- 289 个软件、冻结配置和结果文件与私人研究源文件逐字节一致；原始实验结果未改。 / 289 software, frozen-config, and result files are byte-identical to the private research source; raw results were not changed.
- 中文证明译本与英文正式证明的 18 个带编号公式逐一核对；剔除空白差异后完全一致。译文叙述未进行新的独立数学审计。 / All 18 tagged displays in the Chinese proof translation match the English proof after whitespace normalization; the translated prose has not received a separate mathematical audit.

许可与引用见 [CONTENT_LICENSE.md](CONTENT_LICENSE.md) 及 [CITATION.cff](CITATION.cff)。这份仓库不是同行评审论文，也不宣称已确立学术新颖性。 / See the license map and citation metadata. This repository is not a peer-reviewed paper and does not claim established academic novelty.
