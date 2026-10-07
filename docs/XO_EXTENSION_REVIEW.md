# A/T6 独立技术审查记录 / Independent technical review record

日期：2026-10-07。对象为当日实验室修订稿中的 A、T6 及限定辅助材料，不是原 S3 审计记录的扩写。正式陈述以[补充登记](XO_EXTENSION_REGISTRATION.md)及两份收录证明为准。

## 来源与结论

项目 Owner 提供并接受了本次独立只读审查结论，随后明确授权修正发布准备问题、同步私人和公开仓库并发布 v0.2.0。审查由 Codex 在独立只读任务中重新阅读证明、原 S3 依赖、实现和结果后完成；不是期刊同行评审或人类专家委员会认定。完整本机审查记录、原始运行日志和路径信息保留在私人仓库，不是公开证明或复现的依赖。

审查结论：**A 和 T6 在各自声明范围内未发现实质性证明漏洞，可作为核心补充成果入库。** 该结论不推广有限窗口的 α 范围、不证明周期性、不确认学术新颖性。

The Owner accepted a separate, tool-assisted independent technical review of the deposited A/T6 manuscripts. The review found no substantive proof gap within their separate stated scopes. It is neither journal peer review nor formal machine verification. The private working records are not needed to read the proofs or run the public verifiers.

## 审查核实的证据

- 六个成果脚本重跑：A 的 255 组精确矩阵、T6 的 295 个分类实例和 15 项条件、首次选择公式的 3995 个实例及辅助脚本，结果与当时归档一致（计时除外）。
- 独立公共分母整数算法推进每一步完整质量向量，核对 29 条 α=1/4 轨迹，包括 K=200、400 的完整 10K 窗口；另核对 10 组 K=6 的 α 边界案例。
- 完整状态 exact 检查 K=8、50、200 的特定强制下支规则，各 33 个观察状态；与 200000 步标量浮点代理分开。
- 两条 α 范围外反例已复现；T6 的依赖仍为原 S3，不需要 A。
- 审查时的两个仓库各 20 套既有测试、382 项运行时断言通过；实验室 manifest 的 31 个现行条目及 21 个旧稿条目匹配。该项说明审查时基线，不能替代本版整合候选验收。

## 发布准备问题的处置

| 问题 | v0.2.0 收录方式 |
|---|---|
| 缺正式登记与审查来源 | A/T6 单独登记，本记录与原 S3 审计分开。 |
| 验收依赖实验室 base/ 和目录编号 | 新仓内入口只使用本仓冻结 S2、配置和脚本；输出指定到新目录。 |
| 当前入口仍称有限窗全 K 或其他 α 未覆盖 | 当前 README、双语报告、阅读说明和中文导读分别说明 A/T6；旧 S1/S2/S3 历史范围保留。 |
| S12 长程统计容易误读 | 字段改为 proxy_d_plateau_share、proxy_max_d_excess，并记录标量浮点口径；完整状态 exact 短程核对另列。 |
| C1 参数与名称 | 收录为径向聚合/半共轭，显式写明参数及停机范围，角向衰减按被扫次数。 |
| 核心脚本未纳入日常测试 | 独立精简回归接入 tests/count_tests.py；长矩阵保持另行入口。 |

实际迁移完成后的运行、脱离实验室复现和内容检查以[发布候选验收](XO_EXTENSION_VALIDATION.md)为准。审查的数学结论与发布操作的成功状态分别记录。
