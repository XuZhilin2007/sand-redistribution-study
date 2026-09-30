# S4 Wave 1 / S2 — run 产物说明（XO Boundary Switching Microscope）

canonical 报告：**[docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md](../../docs/S4_W1_S2_BOUNDARY_MICROSCOPE.md)**（本轮唯一权威叙述；本文只登记文件清单与复现方式）。

## 文件

| 文件 | 内容 |
|---|---|
| `config.json` | 冻结预注册（`e706a0d`，执行零改动；sha256 见 metadata） |
| `run_s2.py` | 执行器：Step-0 门（S0a–S0e）→ exact census（冻结分类器）→ C1–C7 约化核验 → H2 锚点逐字符串核验 → metadata |
| `results/step0_s1_csv_checks.csv` | 4 锚点 S0 检查（全部 PASS；S0e 残差 ~1e-63 ≪ 1e-50） |
| `results/per_k_census.csv` | 51 K 逐 K：t*、descent 行程、R_m、tie-tail、观测/预测类别、匹配、C1–C7 摘要 |
| `results/reduction_K{K}.csv` ×51 | 逐轮：t, s, branch, p_m, F_m, restriction/branch_eq/sel_branch/dpred ok, class_predicted/observed, tie_exact |
| `results/metadata.json` | git SHA、config hash、scientific events（0）、输出 sha256、runtime |

## 复现

```bash
.venv/Scripts/python.exe experiments/s4_w1_s2_boundary_microscope/run_s2.py   # ~21 min
.venv/Scripts/python.exe tests/test_s4_s2_microscope.py                        # 33 项断言
```

## 头条结果（详见 canonical 报告）

- mod-6 预测 **51/51 全中**（LOCK ⟺ K mod 6 ∈ {0,1,2}：25 K；CYCLE ⟺ {3,4,5}：26 K）；零 OTHER、零 scientific events；
- 约化核验 C1–C7 全部 51/51：**一维 p_m 映射与 full exact 轨道逐步 Fraction 恒等**；
- 三种 residue 类循环指纹（mod-3 慢 ~10.2 / mod-4 紧 ~2.5 / mod-5 B-重 ~3.3）均由映射精确生成。
