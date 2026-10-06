# 公开版本记录 / Public release history

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
