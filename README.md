# 沙子再分配研究 / Sand Redistribution Study

**公开研究快照 · Public research snapshot — [v0.3.0](https://github.com/XuZhilin2007/sand-redistribution-study/releases/tag/v0.3.0) · 2026-10-21**

本项目研究由“扫回过密区域、再撒开、适时停止”启发的一维确定性质量反馈模型。
v0.3.0 的核心结果是：在声明的 canonical XO 设置下，**全部整数 K≥6，
执行恰两轮使平方均匀度指标在全部时间中达到最优**。本版同时给出长期
选择字的显式周期分类、状态渐近结构、初态推广及平台/端点证书。

This project studies a one-dimensional deterministic mass-feedback model.
Under the canonical XO assumptions, **two updates minimize the squared
nonuniformity over all integer times, for every integer K≥6**.
Version 0.3.0 also deposits explicit eventual selection words, asymptotic
states, initial-state results, and exact platform/endpoint certificates.

## 当前成果 / Results

共同 canonical 前提：初态 p₀,ᵢ=(2K−2i+1)/K²、XO 喷撒律、最小指标平局、
精确算术。不同成果的 α、初态与范围分别读取，完整陈述见
[成果登记](docs/V030_REGISTRATION.md)。

| 结果 | 适用范围与实际保证 |
|---|---|
| [B-lock + E：两轮最优停止](docs/V03_E_CYCLE_ALLK.md) | α=1/4，全部 K≥6、t≥0，U(pₜ)≥U(p₂)，U=KΣ(pᵢ−1/K)²。最优不表示唯一，也不表示 U 随时间单调。 |
| [T7：长期结构](docs/V03_T7_CYCLICITY_PROOF.md) | α=1/4；选择序列最终精确周期，完整质量状态渐近周期。显式字表、例外、全域吸引适用区间见证明。 |
| [B/B-II：初态结果](docs/V03_B_ROBUSTNESS_PROOF.md) | 声明的安全邻域内稳健性，以及合法初态的有限内部回返后永久边界驻留；首入不退与 C=1 不在 B 的保证内。 |
| [C：分类窗口](docs/V03_C_ALPHA_WINDOW_PROOF.md) | 0<α<1 的新两测窗口规则；原冻结规则只在声明充分条件下保证正确。 |
| [F](docs/V03_F_PLATEAU.md)、[G](docs/V03_G_ALPHA_DICTIONARY.md)、[P](docs/V03_P_PLATEAU_ALLK.md) | 旋转数换算、一般 α 的标准形词典、canonical 无限平台族及端点排除；有限格点核对与全称证书分层陈述。 |
| [G-ep](docs/V03_G_ALPHA_DICTIONARY.md)、[C-ep](docs/V03_C_ENDPOINTS_ALPHA.md) | 真实代数 α 端点例外与 13 组精确端点表；不是一般 α 相图或全 K 端点分类。 |
| [原 S3](docs/S3_XO_BOUNDARY_THEOREM_REGISTRATION.md) 与 [v0.2.0 A/T6](docs/XO_EXTENSION_REGISTRATION.md) | 保持各自原范围与审查来源；v0.3.0 以新证明补充早期开放问题。 |

For the main theorem α=1/4. Eventual exact periodicity belongs to the selection
sequence; mass states converge asymptotically to periodic profiles.
The general-α dictionary, classifier windows, initial-state hypotheses and
finite endpoint tables have their own stated scopes.

## 阅读与证据 / Reading and evidence

1. [正式成果登记](docs/V030_REGISTRATION.md)：假设、量词、例外和后续问题。
2. [E](docs/V03_E_CYCLE_ALLK.md) 与 [B-lock](docs/V03_B_LOCK_ALLK.md)：全 K
   最优停止；[D](docs/V03_D_STOPPING.md) 保存有限格点的全时间证书与分解。
3. [T7](docs/V03_T7_CYCLICITY_PROOF.md)：最终选择字、状态渐近与族覆盖。
4. [发布验收](docs/V030_VALIDATION.md) 与 [文件摘要](docs/V030_FILE_MANIFEST.json)：
   实际测试、独立复验和精确字节来源。
5. [文献继承说明](docs/V030_RELATED_WORK.md)：已有分段仿射动力学理论与本模型
   显式结果之间的关系。周期性的一般机制与平台公式已有相关理论。

证明采用解析推导、精确有理数计算与无限族证书；有限轨迹交叉核查另列。
独立实现、AI 辅助审查与程序验收不等同于形式化机器证明或期刊同行评审。
学术新颖性仍需完整文献比较与外部评价。

Analytic reduction and exact computer-assisted certificates carry the
infinite claims. Sampled trajectory checks are recorded separately.
The project is a reproducible research record, not a peer-reviewed paper or
a formally verified proof development.

## 复现 / Reproduction

使用 Python 3.14.6，依赖固定于 [requirements.txt](requirements.txt)。
Windows PowerShell：

```text
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
$env:PYTHONUTF8="1"
.venv\Scripts\python.exe -B tools\materialize_v030_certificates.py
.venv\Scripts\python.exe -B tests\count_tests.py
```

macOS/Linux 先用 `.venv/bin/python -B tools/materialize_v030_certificates.py`
还原证书，再用 `.venv/bin/python` 并以
`PYTHONUTF8=1 .venv/bin/python -B tests/count_tests.py` 运行。
完整回归含精确子族证书复验，可能需要数十分钟。不要使用 `-O`。
科学文件按原始字节保存，避免平台换行转换破坏证书源码 SHA-256。
GitHub 的单文件限制使两个大 JSON 采用 gzip 分片存储于
`certificates/v0.3.0/`。首次运行上述还原脚本，用标准库离线解压并核验
原始 SHA-256；数学代码、原 JSON 字节和所有证明字段保持不变。
还原得到的两个 JSON 不再重复进入 Git，完整源代码下载包已包含所需分片。

进一步独立复验：

```text
python -B experiments/v03_t7_cyclicity/verify_independent.py
python -B experiments/v03_e_cycle_allk/verify_e_independent.py
python -B experiments/v03_d_stopping/verify_d_independent.py
```

The public checkout is self-contained. It needs no private repository,
research-session log or local absolute source path.

## 范围、许可与引用 / Scope, licensing, citation

该模型不是经真实沙坑数据校准的沙粒碰撞模拟；“两轮最优”不代表现实中
任何沙面都应操作两轮。其他 α、指标、初态、网格与喷撒律的最优停止和
扰动稳健性是后续研究。

软件、测试、配置使用 [Apache-2.0](LICENSE)；原创文稿、图表与结果依
[内容许可](CONTENT_LICENSE.md)采用 CC BY 4.0。署名 Zhilin Xu；
引用版本 0.3.0 见 [CITATION.cff](CITATION.cff)。
公开仓库保持独立历史，原 [v0.2.0](https://github.com/XuZhilin2007/sand-redistribution-study/tree/v0.2.0)
与 [v0.1.1](https://github.com/XuZhilin2007/sand-redistribution-study/tree/v0.1.1) 标签保留。
完整变更见 [RELEASE_NOTES.md](RELEASE_NOTES.md)。
