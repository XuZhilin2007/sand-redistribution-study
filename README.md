# 沙子再分配研究 / Sand Redistribution Study

**公开研究快照 · Public research snapshot** — 2026-10-01

本仓库保存一个由沙面修正问题启发的**一维确定性质量反馈模型**、相应代码与实验记录，以及 canonical XO 轨迹的边界约化证明。这里的“沙子”是问题来源；模型不是现实沙粒仿真，结论也不直接预测真实沙坑。主要研究结果已有中英双语[公开研究报告](docs/RESEARCH_REPORT.md)；核心定理同时提供[中文完整证明译本](docs/S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md)与[英文正式证明](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)。

This repository contains a **one-dimensional deterministic mass-feedback model** inspired by a sand-surface correction task, its code and computational record, and a boundary-reduction proof for one canonical XO trajectory. It is not a grain-level simulation or a prediction of real sand. Start with the bilingual [public research report](docs/RESEARCH_REPORT.md). The [formal English proof](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md) has a [full Chinese translation](docs/S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md) and a shorter [Chinese theorem guide](docs/XO_THEOREM_GUIDE_ZH.md).

## 目前的结论 / Current status

| 证据层级 / Evidence | 当前结果 / Result |
|---|---|
| 已有数学结果 / Existing math | 一般单步更新的分解与若干下降判据；canonical XO 的原不变量 T1。 / General one-step decomposition and descent statements; the earlier canonical XO invariant T1. |
| 新数学结果 / New mathematical result | S3 已证明 canonical XO 在明确前提下有限步进入并保持双边界区，边界观测量满足单变量精确递推，且无限时间选择行为按 \(K\bmod6\) 分两类（T2–T5）。 / Under stated assumptions, S3 proves finite boundary entry, persistence, scalar closure of boundary observables, and a residue-class selection dichotomy. |
| 计算核验 / Verified computation | S1 精确算术重放与 S2 的 51 个 \(K\) 的有限时限核验。LOCK/CYCLE 只表示冻结的 \(H=10K\) 实验分类。 / Exact S1 replay and the 51-case S2 finite-horizon checks; LOCK/CYCLE are classifier labels at \(H=10K\). |
| 未知 / Open | 学术新颖性、完整状态的渐近行为、其他规则和初态、二维及现实验证等。 / Literature novelty, full-state asymptotics, other laws and initial states, 2D, and real-world validation. |

S3 状态为 **CLOSED — PASS**，证明已入库。[审计记录](docs/S3_XO_BOUNDARY_THEOREM_AUDIT.md)中的独立 PASS 是项目所有者提供的审计结论；完整审计对话未随仓库发布。定理不宣称周期性、全状态一维化或一般沙子再分配定律。 / S3 is **CLOSED — PASS**, with its proof deposited. The independent PASS verdict in the [audit provenance](docs/S3_XO_BOUNDARY_THEOREM_AUDIT.md) was reported by the project owner; the full audit transcript is not included. The theorem does not claim periodicity, a one-dimensional reduction of the entire state, or a general law of sand redistribution.

## 阅读路线 / Reading order

1. [双语公开研究报告 / Bilingual research report](docs/RESEARCH_REPORT.md) — 问题、模型、证据层级、结果与边界 / question, model, evidence levels, results, limits.
2. [中文定理导读 / Chinese theorem guide](docs/XO_THEOREM_GUIDE_ZH.md)、[中文证明译本 / Full Chinese proof translation](docs/S3_XO_BOUNDARY_THEOREM_PROOF_ZH.md) 与 [英文正式证明 / Formal English proof](docs/S3_XO_BOUNDARY_THEOREM_PROOF.md)；[定理登记 / Registration](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md)给出精确适用范围。
3. [S1 精度报告](docs/S4_W1_S1_XO_PRECISION.md)、[S2 边界显微镜报告](docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md)与 `experiments/s4_w1_s1_xo_precision/`、`experiments/s4_w1_s2_boundary_microscope/` 内的冻结配置和结果。
4. [Stage 2](docs/RESEARCH_CHECKPOINT_2026-09-20.md) 与 [Stage 3](docs/RESEARCH_CHECKPOINT_2026-09-20_STAGE3.md) 历史 checkpoint：分别记录一般单步数学和跨规则对照。后续 S1–S3 结果会取代其中有关 XO 的早期开放问题或浮点解释。 / These dated checkpoints preserve earlier results and their then-current questions; later S1–S3 work supersedes some XO interpretations.

历史报告按写作时点保留，不应把其中的“下一步”或“尚未证明”当作当前状态。请以本页、双语公开报告及 S3 正式证明为当前入口。 / Dated reports are preserved as historical evidence; their old next steps and proof targets are not the live status. Use this page, the bilingual report, and the S3 proof as the current entry points.

## 复现 / Reproduction

已核对的环境为 **Python 3.14.6**、NumPy 2.5.3、Matplotlib 3.11.2、mpmath 1.4.1。创建虚拟环境并安装 [requirements.txt](requirements.txt) 后，在仓库根目录运行： / The checked environment uses **Python 3.14.6** with the package versions above. After creating a virtual environment and installing [requirements.txt](requirements.txt), run from the repository root (Windows PowerShell):

```text
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe tests\count_tests.py
```

macOS/Linux 将执行路径改为 `.venv/bin/python`。该测试入口逐个执行 `tests/test_*.py`。本快照在新建的 Python 3.14.6 虚拟环境中按上述固定版本安装依赖后，**20 套、382 项运行时断言全部通过**（2026-10-01）。/ On macOS/Linux, use `.venv/bin/python`. The runner executes each `tests/test_*.py`. In a fresh Python 3.14.6 environment with the pinned dependencies, this snapshot passed **20 suites and 382 runtime assertions** on 2026-10-01.

已提交的 CSV、JSON 与图像是证据记录；完整 S1/S2 运行可较耗时，不是阅读或运行测试的前提。冻结配置和原始输出应一同使用，不应把一个有限时限的分类当成无限时间定理。 / Committed CSV, JSON, and figures are evidence records. Full S1/S2 runs can be time-consuming and are not required for reading or running the test suite. Use frozen configurations together with their outputs; finite-horizon labels must not be read as asymptotic theorems.

## 仓库与许可 / Repository and licensing

`src/` 为模型实现；`tests/` 为回归测试；`experiments/` 保存冻结配置、执行脚本和结果；`docs/` 保存数学推导、阶段报告和当前双语导读。本仓库从私人研究记录制作**独立历史的公开快照**；内部评审草稿、原始聊天材料和旧归档不在此版本内。 / This is an independent-history public snapshot of the private research record. Internal review drafts, raw chat materials, and the old bootstrap archive are not included.

软件源代码、测试、实验脚本和运行配置采用 **Apache License 2.0**（[LICENSE](LICENSE)）；研究文字、图表和已发布的实验结果采用 **Creative Commons Attribution 4.0 International**（[CONTENT_LICENSE.md](CONTENT_LICENSE.md)）。具体文件边界以该说明为准。引用方式见 [CITATION.cff](CITATION.cff)。 / Software and configurations use **Apache-2.0**; research prose, figures, and published results use **CC BY 4.0**. The [content license note](CONTENT_LICENSE.md) defines the file boundary; see [CITATION.cff](CITATION.cff) for citation information.
