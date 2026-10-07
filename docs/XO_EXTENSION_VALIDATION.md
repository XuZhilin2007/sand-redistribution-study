# v0.2.0 补充成果验收 / Supplement validation

本页登记 v0.2.0 整合候选的实际运行结果。2026-10-07 已完成仓内完整验证、两仓常规回归的隔离执行、脱离实验室的完整复现和 GitHub Markdown 数学表达式保留检查。此前独立审查来源单独记录，未用其 PASS 替代本版验收。

## 仓内命令 / Repository-local commands

在仓库根目录运行，输出目录应为新目录或空目录；入口拒绝覆盖已有证据。补充脚本只使用 Python 标准库，无需实验室、base/ 或另一个仓库。

```text
python -B tests/test_xo_extensions.py
python -B experiments/xo_extensions/validate.py --scope core --out-dir ../xo-core-validation
python -B experiments/xo_extensions/validate.py --scope auxiliary --out-dir ../xo-aux-validation
python -B experiments/xo_extensions/validate.py --scope full --out-dir ../xo-full-validation
```

既有模型回归需安装根目录 requirements.txt 的固定依赖，再运行 `python -B tests/count_tests.py`。不能使用 `python -O`：这些程序以断言核验科学条件。

Run the supplement with the standard library alone. The entry reads this checkout's frozen S2 classifier and uses no external lab. Use a fresh output directory. The existing model suites additionally require the pinned project dependencies. Optimized Python execution is rejected because assertions are part of the checks.

## 冻结配置与口径 / Frozen configurations

| 核验 | 配置 | 算术和实际窗口 |
|---|---|---|
| A | [config](../experiments/xo_alpha_extension/config.json) | 255 组 Fraction 完整状态矩阵：K=6..53 为 10K，K=100/200/400 为 300；另 7 组 K=6 为 400。255 组标量 float 各 200000 步为独立有限代理。 |
| T6 | [config](../experiments/xo_t6_classifier/config.json) | α=1/4，295 个 K，均 H=10K；完整状态 Fraction 进入后用精确标量闭包；另外核对原 S2 分类函数。 |
| 辅助 | [config](../experiments/xo_supporting_results/config.json) | F/D/S6/C1 为精确有理或整数；S12 20 组各 200000 步为标量浮点代理。 |
| 独立交叉核对 | [config](../experiments/xo_extensions/config.json) | 29 条完整整数质量向量轨迹（含 K=400,H=4000）、10 组 K=6；特定强制规则完整状态各 33 步。 |

前三份配置描述原验证脚本的固定设置，统一入口复检配置与结果；第四份配置直接控制独立完整状态矩阵。所有 H 项均为更新前的选择状态 t=0,…,H−1。源码 SHA-256、配置、实际 α/K/H、Python 版本和算术类型随 manifest.json 记录；有限抽查不承担证明的全称量词。

## 整合候选结果 / Integrated candidate results


**实际结果：PASS。** [源码与配置 manifest](../experiments/xo_extensions/results/v0.2.0/manifest.json)记录 Python 3.14.6、完整 source_sha256、冻结设置和本轮全部结果；[整合验收](../experiments/xo_extensions/results/v0.2.0/integrated_validation.json)记录两仓回归与隔离复现。

- 两仓候选均为 **21 套、441 项运行时断言**：原 382 项加补充 59 项。既有 S2 测试会回写 CSV，故在与候选代码逐字节相同的隔离副本执行；正式仓库的旧实验结果未改。日志：[private](../experiments/xo_extensions/results/v0.2.0/routine_private.log)、[public](../experiments/xo_extensions/results/v0.2.0/routine_public.log)。
- A：[255 组精确矩阵与 7 组 K=6、255 组有限 float 代理](../experiments/xo_extensions/results/v0.2.0/A.json)，代理总步数 51000000，实际时限按逐例字段读取。
- T6：[295 个实例与 15 项条件](../experiments/xo_extensions/results/v0.2.0/T6.json)，[原冻结函数 295 个交叉分类](../experiments/xo_extensions/results/v0.2.0/frozen_classifier.json)。
- 独立核对：[29 条完整整数质量向量轨迹、10 组 K=6 参数边界、402 个有限代数配置及特定强制规则 exact 完整状态短程](../experiments/xo_extensions/results/v0.2.0/independent.json)。K=200/400 均跑完整 10K 窗口；强制规则 K=8/50/200 各 33 个观察状态。
- 辅助结果：[F](../experiments/xo_extensions/results/v0.2.0/F.json)、[首次选择 3995 个实例和 55 条第二步数据](../experiments/xo_extensions/results/v0.2.0/D.json)、[20 组 S12 标量代理与 12 个有效初态](../experiments/xo_extensions/results/v0.2.0/probes.json)、[C1 的 5 次活跃 exact 聚合](../experiments/xo_extensions/results/v0.2.0/C1.json)。
- 在没有实验室、base/、.git 或 venv 的新目录中，用 `python -S -B` 重跑完整入口：**PASS**。禁用第三方 site-packages，结果与最终候选逐项一致（计时除外），执行代码与配置 SHA-256 一致；没有借用私人仓库的运行模块。
- GitHub 官方 GFM 接口核对十份主要文稿，全部受保护的行内/显示公式在输出中保留。修正了既有换行接 \ell 和显示公式集合括号的 Markdown 转义问题；S3 中英文 18 个编号公式仍完全对应。仅使用等价 TeX 宏与换行调整，原数学内容及范围不变；这不是浏览器像素级检查。
- 公开相对链接、排除目录、定向凭据/个人路径/私人邮箱检查通过；原 S1/S2 冻结材料与结果哈希保留，实验室 31 个现行与 21 个旧稿 manifest 条目保持匹配。完整私人发布对应记录另存版本登记，不进入公开历史。

All six verifiers and both integrated candidate suites passed. The isolated standard-library reproduction matched all computational results apart from timings. The source/config hashes, exact full-state checks, float-proxy labels, and GFM math checks are recorded separately. Original S3 received equivalent formatting only; the scientific scope and frozen S1/S2 record are preserved.
