# S4 Wave 1 / S1 — run 产物说明（XO 高精度重放）

canonical 报告：**[docs/S4_W1_S1_XO_PRECISION.md](../../docs/S4_W1_S1_XO_PRECISION.md)**（本轮唯一权威叙述；本文只登记文件清单与复现方式）。

## 文件

| 文件 | 内容 |
|---|---|
| `config.json` | 轮前冻结的预注册 config（tiers、可行性门、比较规则、halt 条件；任何计算前写定） |
| `replay.py` | 三层重放引擎：exact Fraction / mpmath(50+15, 100+15 guard) / (float64 由 committed runner 提供)；内嵌逐轮自检（worst=d-增量恒等式、D_{t,m}=b(1−b)、s≤m、d>tol） |
| `run_s1.py` | 编排：Stage V（hash + 逐字符串再生核验）→ Stage R（阶梯）→ Stage C（比较/行程/元数据） |
| `results/stage_v_verification.csv` | 8 组再生核验（XO 4 K + KC 4 K），零失配 |
| `results/xo_{tier}_K{K}.csv`, `kc_*` | 各 tier 逐轮诊断（t, s, d, e, moved, worst, R, boundary, tie_exact, invariant, s_bound, rebound） |
| `results/comparison_pairs.csv` | tier 对比较（selection 零分歧、rebound 旗标失配轮、max\|Δd\|） |
| `results/summary_by_tier.csv` | 每 K×tier 汇总（N_R、ties、distinct s、\|e\| 带、dvals、maxR、wall） |
| `results/itinerary_xo.csv` | selection 行程（visits/runs/dwell/top transitions，全部 tier） |
| `results/comparison_kc.csv` | KC sanity（t_stop / N_R：F vs mp100 vs exact@K=50） |
| `results/metadata.json` | config 回显、integrity、环境、git SHA、输出 sha256 |

## 复现

```bash
.venv/Scripts/python.exe experiments/s4_w1_s1_xo_precision/run_s1.py   # ~19 min, 1152 s 记录
.venv/Scripts/python.exe tests/test_s4_s1_precision.py                 # 19 项断言
```

## 核验层（全部通过）

1. Stage V：committed CSV sha256 ↔ metadata 一致；XO 7500 + KC 5344 行逐字符串再生零失配；
2. exact tier 逐轮自检（四 K 全程 Fraction 级）；
3. repo helper（rebound_residual / target_matched_excess）交叉验证 ≤1e-15；
4. mp↔exact：selection 零分歧、d < 1e-38/1e-55（K=50 片段）；
5. 独立重实现（未读本引擎的 Fraction 复算）K=50 全程 500 轮全部吻合；
6. K=6 手算校验（D profile、不变量实例）。
