# 真值表（Ground Truth）

本目录单独记录 A1..A12 十二套基准数据集的**全部真值**，用于独立验证归档真伪，与归档本身分离，不依赖归档自述。

## 文件

| 文件 | 说明 |
| --- | --- |
| `truth_table.csv` | 每套数据集一行，共 12 行 32 列，机器可读真值表 |
| `build_truth.R` | 从源 `.rds` 生成真值表的脚本（可复现） |
| `verify_archive.R` | 验证任意归档真伪的脚本 |

## truth_table.csv 字段

- **身份**：`dataset`、`base_seed`、`N`、`n_treated`、`n_donor`、`D`、`T`、`t0`、`pre/post`、`n_rep`、`time`、`treatment`
- **效应真值**：`effect_size`、`n_effect_taxa`、`effect_taxa`（6 个效应菌群索引，分号分隔）、`effect_signs`（方向）；null 数据集（A11）为空
- **生成参数**：`baseline_sd`、`trend_sd`、`innov_sd`、`ref_dt`、`rho`、`phi_range`、`base_depth`、`depth_cv`、`structural_zero`、`trend_diff`、`effect_hetero`、`twin_trend`
- **指纹**：`relab_true_md5`、`Y_md5`（源 `.rds` 张量的 md5，用于核对源文件未被改动）
- **可复现**：`rep_seeds`（每套 30 个重复的随机种子）

有了参数 + 种子，可完整复现数据；有了指纹，可快速核对源 `.rds` 与真值表一致。

## 验证真伪

```bash
Rscript verify_archive.R ../panelio/A6/A6_rep12_v1.tar.gz
```

脚本对任意归档执行四类检查，全部通过输出 **AUTHENTIC**，任一失败输出 **CHECK FAILED**：

1. **结构身份**：N / D / T、`elapsed_time`、`t0` 干预断点、`label` 分组（treat/ctrl）与真值表一致；
2. **效应真值**：`effect_size`、`effect_taxa`、`effect_signs` 与真值表一致；
3. **源指纹**：源 `.rds` 的 `relab_true` 指纹与真值表一致（确认源本身未篡改）；
4. **回读保真**：`read_panel_table` 读回的丰度与真值 `relab_true` 张量逐位相等（max deviation = 0）。

## 重要说明

- `panelio` 的 `read_panel_table` 在**写入路径**强制组成和（和为 1）等校验，但**读入路径**只做列/类型/结构校验，**不**重新强制组成和与平衡性——因此仅靠读入函数**无法识别**丰度被篡改的归档。真伪必须通过本目录的真值表比对（即 `verify_archive.R`）确认。篡改测试已验证：把 `time_001.csv` 某个丰度改为 999 后，脚本正确判定 **CHECK FAILED**。
- 归档存的是组成型相对丰度 `relab_true`；原始计数 `Y`、深度 `depth`、效应真值由源 `.rds` 保存，真值表记录其指纹与复现信息。
