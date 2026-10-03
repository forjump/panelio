# CLSC Benchmark Paper (panelio) — 章节结构说明

论文源码：`paper/panelio.tex`（LaTeX article 类，**未编译**，交付源码；图以占位文字说明，
表格用 `tabular`/`booktabs`，公式用数学环境）。正文所有数值均来自本次真实跑出的
panelio 基准结果，无编造数字；论文章节组织借鉴用户提供的 CLSC.docx（其全文提取见
`benchmark/results/CLSC_paper_extract.md`），方法与讨论表述参考其写法，数值一律为本轮实测。

## 章节结构

| 章节 | 内容 | 数据来源 |
|---|---|---|
| 摘要 | 框架定位 + 三项核心结论（Power/校准、消融、敏感性） | 全表汇总 |
| 1 引言 | 微生物组因果推断难点 → 现有方法不足 → CLSC 三轴分解 → 三点贡献 | — |
| 2 方法 | 2.1 问题设定与数据布局；2.2 成分预处理（滤过/零填充/CLR，式 1）；2.3 负二项观测模型与 softmax 映射（式 2–3）；2.4 带干预断点的时间平滑；2.5 四类因果估计器（DiD 式 4、Pooled OLS 式 5、Event study 式 6、simplex SCM 式 7–9）；2.6 推断层（placebo 置换 + BH-FDR、subject bootstrap、pre-fit 诊断式 10）；2.7 评估指标口径（Type I / Power / Bias / MAE / 反事实重建） | — |
| 3 仿真设计 | panelio A1–A12 × 30 rep 设计参数总表（表 1：T、t0、效应大小/特征/符号、innov_sd、深度、特设参数） | truth_table.csv、dataset_metadata.csv |
| 4 结果 | 4.1 四方法 Power（表 2）与 Type I（表 3）；4.2 符号命中与估计误差（表 4，softmax 压缩与符号翻转解释）；4.3 SCM 诊断（pre-RMSE≈3.6–4.0、SNR≈0.84–0.87、placebo 保守）与反事实重建 | summary_simulation.csv（72,000 行） |
| 5 消融 | 11 变体 × A1/A2/A6/A11/A12 的 Power/Type I/Bias/MAE（表 5）+ 四条结论（凸约束与 offset 关键、断点平滑权衡、滤过惰性、错误尺度灾难性） | summary_ablation.csv（82,500 行） |
| 6 敏感性 | 滤过阈值（表 6，三档全同）、均一化（表 7，仅 CLR/ILR/rCLR 校准）、零填充（表 8）、平滑带宽（表 9）、bootstrap 推断（表 10，A2 覆盖偏低原因） | summary_sensitivity_{filter,normalize,zero,band,bootstrap}.csv |
| 7 讨论 | 模块归因 + 五条局限/未来方向 + placebo 保守性说明 | — |
| 8 结论 | 收束 | — |
| 参考文献 | 24 条（沿用 CLSC.docx 文献列表，Harvard 风格转 bibitem） | CLSC_paper_extract.md |

## 数据与脚本

- 主实验：`benchmark/scripts/11_sim_shard.R`（4 分片并行）→ `results/sim_s{1..4}_perfeature.csv`
  + `sim_s{1..4}_diag.csv`（合计 72,000 行，16 列统一 schema，0 错误行）。
- 消融：`benchmark/scripts/12_ablation.R`（A1/A2/A6/A11/A12）→ `results/abl_*_perfeature.csv`
  （82,500 行）。
- 敏感性：`benchmark/scripts/13_sensitivity.R`（filter/normalize/zero/band/bootstrap 五类）→
  `results/sen_*_perfeature.csv`。
- 汇总：`benchmark/scripts/14_summarize.R` → `results/summary_*.csv`（论文全部表格数字的原始来源）。
- 真值口径：仅 `benchmark/truth/truth_table.csv` 的 effect_taxa / effect_signs / effect_size。

## 复现运行

```powershell
# 主实验（4 分片）
wsl -d Ubuntu-26.04 bash -lc "Rscript /mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/scripts/11_sim_shard.R A1 A2 A3 s1"  # 依次 s2=A4..A6, s3=A7..A9, s4=A10..A12
# 消融
wsl -d Ubuntu-26.04 bash -lc "Rscript /mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/scripts/12_ablation.R A1 A2 A6 A11 A12 abl_full"
# 敏感性
wsl -d Ubuntu-26.04 bash -lc "Rscript /mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/scripts/13_sensitivity.R <filter|normalize|zero|band|bootstrap> <A1 A11 A2...> <out_tag>"
# 汇总
wsl -d Ubuntu-26.04 bash -lc "Rscript /mnt/c/Users/Zhang/Desktop/clsc_r/benchmark/scripts/14_summarize.R <sim|abl|sen_*>" 
```

## 关键结论速览（真实值）

- **主实验**：大效应设计 pooled OLS Power 最高（A8 0.267、A12 0.272），SCM 0.11–0.23；
  A11（null）上 did/es Type I≈0.01/0.000、pooled 0.092、scm 0.090（placebo FDR=0，保守）；
  符号命中 0.42–0.63（softmax 闭包 + 干预点大创新噪声所致）；A12 bias≈−0.79–−0.81（trend_diff 混淆）。
- **消融**：nnls 失凸（Σw 1.29–3.01）；去 offset 后 A6 Power 0.233→0.028、A12 MAE 1.18→1.63；
  断点 Kalman 抹平效应（recon_null 0.052）；NB 模型在相对丰度上 Type I≈0.89–0.99。
- **敏感性**：滤过阈值三档无差异；仅 CLR/ILR/rCLR 在 null 上 Type I 达标（0.075–0.232），
  其余 6 种归一化 0.30–0.54；零填充各法 ±0.02；平滑带宽单调权衡；bootstrap 覆盖 A1 0.57–0.90、
  A2 0.27–0.43（大效应压缩所致，宜读作可重复性度量）。
