> 公开版本说明：以下保留 2026-10-24 研究冻结的登记/验收记录。
> v0.4.0 公开发布日期为 2026-10-26；公开快照验收见
> [V040_PUBLICATION](V040_PUBLICATION.md)。

# v0.4.0 成果验收与复验

日期：2026-10-24。被复核提交：`85a6d06`；主研究基线：`36d71a2`。
本轮将 v0.4.0 的定理、证明、证据和验收整理为独立的本地成果快照。
正式数学范围见 [V040_REGISTRATION](V040_REGISTRATION.md)，固定输入见
[V040_MANIFEST](V040_MANIFEST.json)。

## 1. 实际数学验收

环境：Windows、Python 3.14；下述快速复验仅使用 Python 标准库。
本轮六项顺序执行的退出码均为 0，合计约 49 秒。

| 检查 | 实际结果与范围 |
|---|---|
| 独立有限区间验收 common | 203/203 K、207 分片、80,246 条精确符号义务，整个 J、全时间 |
| 前缀引理 | 有理常数及 570 组精确交叉核验通过 |
| 全 K W0 | 解析常数、精确多项式、5535 点独立重放及 24 组完整向量核对通过 |
| D₂ 短证明 | μ、逐段上界、边界下界和小 K 衔接全部通过 |
| 统一尾部 | 48 个闭盒全覆盖；177 组完整向量回代通过，另核对缩放 Gram 恒等式 |
| 结构与预算回归 | 四项单字段错误被拒、恢复门排除坏记录、真实预算耗尽返回 OPEN；整套退出码 0 |

本轮运行输出保存于 `experiments/v04_all_k/review_2026-10-24/`：
`finite_common.json`、`prefix.json`、`w0.json`、`d2.json`、`tail.json` 和
`reverification.json`。最后一份记录完整命令、真实退出码、输出和输入不变检查。

原 v3 见证仍为 203 行；修复前后 SHA-256 均为

```
c8c7acae690f939397b88084784df4162264ccef54e5f93c32b8d5e7c979b9d1
```

重复运行的耗时、Git 上下文和新报告字节可以不同；数学输入、原见证与固定清单应保持一致。有限向量核对负责发现装配错误，主定理的无限参数量词由解析推导、闭盒外包络和全时间有限证书承担。

## 2. 一条命令完成复验

从仓库根目录运行：

```powershell
python -B tools/verify_v040.py
```

它先核对固定清单，再顺序执行五个数学审计器和受影响的回归，最后再次核对全部固定输入。任何缺文件、哈希不符或非零退出码均使总验收失败。已有虚拟环境也可使用 `.venv/Scripts/python.exe`。

默认只输出结果；保存一次新的运行记录可使用

```powershell
python -B tools/verify_v040.py --report scratch/v040-verification.json
```

该入口拒绝把报告写到任何固定输入或 manifest 自身。每个单独审计器也支持可选的 `--report`，默认不写回旧报告。

逐项复验命令：

```powershell
$env:PYTHONUTF8="1"
python -B tools/audit_v04_interval_v3.py --scope common
python -B tools/audit_v04_prefix_lemma.py
python -B tools/audit_v04_w0_uniform.py
python -B tools/audit_v04_d2_entry.py
python -B tools/audit_v04_uniform_tail.py
python -B experiments/v04_w_interval/negative_tests_w.py
```

命令中的 `python` 指向相应 Python 解释器；macOS/Linux 使用本机 Python 3 解释器运行相同的脚本路径。不要使用 `-O`。如需检查每 K 自己登记的较宽区间，可另运行 `audit_v04_interval_v3.py --scope recorded`；本版主定理只需 common 范围。

## 3. 冻结与传输修正

- W0、D₂、统一尾部和前缀核验器改为默认不写报告，报告输出显式指定，避免复验修改固定工件。
- 读取历史来源哈希时兼容 Windows 与 POSIX 路径分隔符；新报告使用 POSIX 相对路径。
- 固定输入按 `.gitattributes` 保留实际已验收字节，阻止自动 CRLF/LF 转换。原 v3 工件保留原 SHA，Git 存储字节应与 manifest 一致。
- 冻结清单覆盖主证明、审计代码、运行所需数据、生产/回归代码、历史小 K 接受记录和本轮报告；不把历史生成器版本号视为独立数学验收。
- 清单不包含自身哈希或要求把最终提交号写回自身，避免循环。最终本地冻结提交号在交付记录中报告。

隔离新检出复验使用冻结快照中的同一命令，要求全部固定文件哈希匹配、六项退出码为 0、运行前后固定输入不变。运行产物另存于检出目录之外。

## 4. 完成边界

全 K、整个 J、全部整数时间的证明及有限补集已闭合；结构/预算门已经验收；原见证按原哈希保留。当前范围没有继续扩扫、扩大参数窗口或追加证明单元的待办。
