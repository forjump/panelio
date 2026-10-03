# panelio

balanced longitudinal panel 数据的归档、读写与对象转换 R 模块。

一个 dataset 被封装为一个 `tar.gz` 归档：内部按时间层存放多张宽表，外加一份机器可解析的 `manifest.txt` 承载全部观测元信息与时间层元信息。模块提供无损的读写闭环，以及 long_table 与 3D array 两种表示之间的双向等价契约。

**设计约束**：全部读/写/转换/预处理/因果路径为**向量化**实现，源码不含任何 `for` / `while` 原生循环；时序平滑与因果估计的迭代全部下沉到 `src/panelio_smoothing.f90` / `src/panelio_causal.f90` 的 Fortran 内核，R 层每个方法只做一次 `.Fortran` 派发，同样无 R 层循环；运行时**只依赖 base R**（`utils`、`stats`），除 C++/Fortran 外不引入多余依赖；所有可调参数集中为 `R/constants.R` 中的**命名常量**，无魔数、无硬编码路径、无外部配置与 bash 启动。

## 环境要求

- R >= 4.1.0（原生管道 `|>` 与简写函数 `\()` 需要 4.1+）
- GNU Fortran（`gfortran`）编译器（时序平滑的 Fortran 内核需要；包声明了 `NeedsCompilation: yes`）
- `tar` 命令（归档打包/解包依赖外部 `tar -C`，因此代码中不存在 `setwd()`）

## 安装

```r
# 开发安装（需要 pkgload）
pkgload::load_all(".")

# 或安装到库
R CMD INSTALL .
```

## 目录结构

```text
project_root/
├── R/
│   ├── constants.R          # 全部命名常量（无魔数），无外部配置
│   ├── utils_panel.R        # 内部校验、manifest 文法、结构转换（不导出）
│   ├── convert_panel.R      # as_panel_array / as_panel_table / as_*_meta
│   ├── read_panel.R         # read_panel_table
│   ├── write_panel.R        # write_panel_archive / write_manifest
│   ├── preprocess.R         # 预处理：filter_features / zero_fill_features / normalize_features / preprocess_array / preprocess_panel
│   ├── smooth.R             # 时序平滑：smooth_array / smooth_panel
│   ├── causal.R             # 因果估计总入口：panel_causal_estimate（did / pooled_ols / event_study / scm）
│   └── panelio-package.R    # 包级文档 + useDynLib 注册
├── src/
│   ├── panelio_smoothing.f90   # 平滑 Fortran 内核（kalman/locreg/ewma/spline）
│   ├── panelio_causal.f90      # 因果 Fortran 内核（demean2 / nnls_solve / scm_effects）
│   └── init.c                  # .Fortran 例程注册表
├── man/                     # roxygen2 自动生成的英文函数文档
├── tests/
│   ├── testthat.R
│   └── testthat/            # 单元测试：读写闭环 + 边界校验
├── DESCRIPTION
├── NAMESPACE                # roxygen2 自动生成
└── README.md
```

## 公开接口（15 个导出函数）

### I/O 与对象转换（7 个）

| 函数 | 方向 | 作用 |
| --- | --- | --- |
| `read_panel_table(archive_path)` | 读 | 读取 `tar.gz` 归档，返回长格式 `data.frame` |
| `write_panel_archive(x, dataset_name, version, out_dir)` | 写 | 长表或 `panel_array` → 归档，返回归档完整路径 |
| `write_manifest(observation_meta, time_meta, file_index, comments, outfile)` | 写 | 生成机器可读 `manifest.txt`，返回路径 |
| `as_panel_array(x)` | 转换 | 长表 → `panel_array` 对象（并附 S3 `dim`/`dimnames`/`print`） |
| `as_panel_table(x, abundance_array, observation_meta, time_meta)` | 转换 | 归档路径 → 长表；或数组+两套元信息 → 长表 |
| `as_observation_meta(x)` | 转换 | 长表 / `panel_array` / manifest 文本 → 观测元信息表 |
| `as_time_meta(x)` | 转换 | 长表 / `panel_array` / manifest 文本 → 时间层元信息表 |

其余全部为内部函数（`utils_panel.R`），不导出。

### 预处理（物种轴，5 个）

| 函数 | 作用 |
| --- | --- |
| `filter_features(abundance_array, ...)` | 按 prevalence / 丰度判据滤除特征（`none` = 跳过） |
| `zero_fill_features(abundance_array, method, ...)` | 零填充：`none` / `constant` / `half_min` / `sqrt_min` / `multiplicative` / `bayesian` |
| `normalize_features(abundance_array, method, ...)` | 均一化：`total` / `clr` / `ilr` / `rclr` / `z_score` / `min_max` / `css` / `rle` / `mclr` / `none` 共 10 种 |
| `preprocess_array(abundance_array, filter, zero_fill, normalize, ...)` | 三步管道，全部沿物种轴（维度 2），可自由组合 |
| `preprocess_panel(panel_array, ...)` | 在 `panel_array` 上执行 `preprocess_array` 并保留元信息 |

所有判据/方法/数值均可传 `NULL`（落到命名默认值）；每个样本逐列处理，观测轴与时间轴保持不变。

### 时序平滑（时间轴，2 个）

| 函数 | 作用 |
| --- | --- |
| `smooth_array(abundance_array, method, breakpoints, time_meta, ...)` | 沿时间轴（维度 3）平滑，R 层单次 `.Fortran` 派发 |
| `smooth_panel(panel_array, ...)` | 在 `panel_array` 上平滑并自动从 `time_meta` 派生断点 |

平滑方法（`SMOOTH_METHODS`，8 种）：

- **Kalman 族（默认，忽略断点）**：`ekf`（标量局部水平前向滤波，默认）、`rts`（Rauch-Tung-Striebel 后向平滑）。
- **benchmark 族（必须传断点，否则会把干预效应抹掉）**：`gaussian`（高斯核局部加权，band=σ）、`loess`（tricube 局部线性，span）、`savgol`（Savitzky-Golay，degree + halfwindow）、`moving_average`、`ewma`（递归加权平均，每段重置）、`spline`（二次差分惩罚的罚样条，lambda）。

断点语义：`breakpoints` 为 1-based 时间下标，标记**新分段起点**（干预后）；Kalman 族整条序列平滑、忽略断点；其余方法按断点切成段、每段独立平滑，从而保留干预效应。未显式传断点时，若 `time_meta` 含 `is_intervention_timepoint`，则以首个干预时刻作为断点。所有平滑的迭代都在 Fortran 内核内完成，R 层无任何 `for`/`while`。

**论文深度（Kalman 族观测模型）**：`smooth_array(..., observation_model, hierarchical, ...)` 支持 `observation_model = "nb"`，把零值直接送入**负二项观测方程**（不预先插补）：潜态为 log 相对丰度，观测 mean = lib·exp(z)，delta 方法方差 = mean + mean²/φ。`hierarchical = TRUE` 时用**层级共享参数的 EM**（约 80 次迭代）估计每特征创新方差 q 与离差 φ，`library_size` / `dispersion` / `em_iterations` 均可显式给定或置 NULL；`time_scaling = TRUE` 时过程噪声按采样间隔 Δt 缩放，适配不规则采样。`ekf` / `rts` 均忽略断点。

### 因果估计（1 个总入口）

| 函数 | 作用 |
| --- | --- |
| `panel_causal_estimate(formula, data, method, ...)` | 因果效应估计总入口；报告每特征点估计、标准误、t 统计量、p 值、Benjamini-Hochberg FDR 及逐时间点效应 |

- **公式**：fixest 风格三段式 `feature ~ label | subject + time`；`|` 后固定效应项接受 subject/observation/time 等别名。
- **数据**：`panel_array`、三维丰度数组（另需 `observation_meta` / `time_meta`）或长表（自动转 `panel_array`）。
- **方法**（`CAUSAL_METHODS`，4 种）：`did`（双重差分，默认）、`pooled_ols`（双向固定效应 OLS）、`event_study`（动态事件时间效应，基期固定 rel=-1）、`scm`（合成控制）。
- **scm 权重**：`scm_constraint = "nnls"`（默认，w≥0，NNLS KKT 最优）或 `"simplex"`（论文深度，w≥0 且 Σw=1，反事实落在供体凸包内，由**精确 active-set KKT 求解器**解出；`scm_offset = TRUE` 时吸收标量水平偏移）。`scm` 输出 `$weights`、`$weight_sum`（应恒为 1）与逐时间 gap。
- **实现**：双向固定效应去均值与合成控制求解在 `src/panelio_causal.f90` 内完成；R 层每个估计器只派发一次 `.Fortran`。
- **论文推断层**（全部可选、默认关）：`placebo = TRUE` 做 placebo 置换检验（每个对照单元轮流作拟处理、其余为供体池，每特征 p = |placebo|≥|real| 的比例，再 BH-FDR）；`bootstrap = TRUE` 做个体级重采样百分位 CI（`bootstrap_samples` 默认 500、`bootstrap_level` 默认 0.95）；`pre_fit = TRUE`（scm）报告 Abadie pre-RMSE 与 SNR。结果分别落在 `$effects` 的 `placebo_p_value`/`placebo_fdr`/`bootstrap_lower`/`bootstrap_upper` 及 `$placebo`/`$bootstrap`/`$pre_fit`。

## 调用示例

```r
# 读归档 → 长表
long_df <- read_panel_table("mydata_v1.tar.gz")

# 长表 → panel_array 对象
pa <- as_panel_array(long_df)

# 长表 → 两套元信息表
obs_meta <- as_observation_meta(long_df)
time_meta <- as_time_meta(long_df)

# 数组 + 元信息 → 长表（双向契约还原）
long_df2 <- as_panel_table(
  x = NULL,
  abundance_array = pa$abundance_array,
  observation_meta = obs_meta,
  time_meta = time_meta
)

# 对象 → 归档（无损闭环）
tar_path <- write_panel_archive(pa, dataset_name = "mydata", version = "v2")

# 预处理：先滤过，再零填充，再 CLR 均一化（全部沿物种轴）
pre <- preprocess_array(pa$abundance_array,
                        filter = list(min_prevalence = 0.2, min_abundance = 1e-6),
                        zero_fill = list(method = "half_min"),
                        normalize = list(method = "clr"))

# 时序平滑：默认 EKF；benchmark 方法带断点保护干预效应
sm_e <- smooth_array(pa$abundance_array, method = "ekf")
sm_g <- smooth_array(pa$abundance_array, method = "gaussian", breakpoints = 5L)
sm_p <- smooth_panel(pa, method = "loess")

# 因果估计：did（默认）对比真实干预特征
res <- panel_causal_estimate(feature ~ label | subject + time, data = pa, method = "did")
head(res$effects[, c("feature", "estimate", "p_value", "fdr")])
# 合成控制（输出权重）
res_scm <- panel_causal_estimate(feature ~ label | subject + time, data = pa, method = "scm")
dim(res_scm$weights)   # n_treated x n_control x n_feature
```

## 归档格式

归档命名：`{dataset_name}_{version}.tar.gz`

归档内部文件：

```text
manifest.txt
time_001.csv
time_002.csv
...
time_0N.csv
```

### `time_xxx.csv`（每个时间层一张宽表）

| 列 | 含义 |
| --- | --- |
| 第一列 `observation_id` | 样本唯一标识 |
| 其余列 `feature_0001` … | 特征名 |
| 单元格 | 丰度数值 |

约束：

- 全部 `time_*.csv` 的 `observation_id` 集合完全一致，且与 manifest 的 `observation_meta` 一致；
- 全部时间层特征列名与顺序完全一致；
- 丰度是**组成型**数据：每个观测在每个时间层内所有特征丰度和为 1（容差 `COMPOSITION_TOLERANCE`）；
- 本工具**仅支持平衡面板**，任何缺失单元格直接报错。

### `manifest.txt`（机器可解析）

键值对 `key = value`，同一行多字段以 `;` 分隔，`#` 开头为注释，区块由 `[name]` 界定，共三块：`[observation_meta]`、`[time_meta]`、`[file_index]`。`label` 以文本存放，物化为 `data.frame` 时转成 R factor。时间字段统一为 `elapsed_time`，单位由数据集自身定义。

### `delta_t`（时间间隔）

`delta_t` 是当前时间层与上一时间层的间隔：首层固定为 `0`，其余 `delta_t = elapsed_time[i] - elapsed_time[i-1]`，**允许不均匀**。由长表派生时按 `elapsed_time` 排序后再求差，因此打乱行序的长表与有序长表结果一致；该一致性在读取时被硬校验。

## 对象模型

### 1. 长格式 table（`data.frame`）

复合主键 `(observation, time_identifier, feature)` 三列联合唯一。

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `observation` | character | 观测标识（对应 manifest `observation_id`） |
| `time_identifier` | character | 时间层标识，形如 `time_001` |
| `elapsed_time` | numeric | 通用累积时间 |
| `feature` | character | 特征标识 |
| `abundance_value` | numeric | 丰度，组成型，观测 × 时间层内和为 1 |
| `label` | **factor** | 分组变量，水平数 ≤ 4 |
| `is_intervention` | logical | 是否处于干预期，全表至多一个断点 |

行顺序为规范顺序 `time_identifier × observation × feature`（feature 变化最快），与归档逐层拼接的物理布局一致。

### 2. S3 类 `panel_array`

```text
panel_array (S3 class)
├── abundance_array      # 三维 array[observation, feature, time_identifier]
├── observation_meta     # data.frame，列名 observation_id / label
└── time_meta            # data.frame，列名 time_identifier / elapsed_time / delta_t / is_intervention_timepoint
```

## 数据模型约定（硬校验，全部在 `R/utils_panel.R`）

1. `label` 是分组变量：必须是 R factor，水平数 ≤ `MAX_LABEL_GROUP_COUNT`（4）。
2. 至多一个干预断点：按 `elapsed_time` 升序，最早层必须为 `FALSE`，一旦变 `TRUE` 不得回退；全 `FALSE` 合法。
3. 丰度组成型：每 `观测 × 时间层` 丰度和为 `COMPOSITION_SUM_TARGET`（1），容差 `COMPOSITION_TOLERANCE`（1e-6）。
4. `delta_t` 允许不均匀：与 `elapsed_time` 逐层差值一致（容差 `DELTA_T_TOLERANCE`）。

## 无魔数与向量化

- 所有可调数值、前缀、分隔符等集中为 `R/constants.R` 中的命名常量（`MAX_LABEL_GROUP_COUNT`、`COMPOSITION_TOLERANCE`、`CSV_WRITE_DIGITS`、`TIME_IDENTIFIER_DIGITS` 等），源码中无未命名魔数、无外部配置文件、无 yaml 依赖。
- 读/写/转换/预处理路径全部向量化：`lapply` / `Map` / `vapply` / `aperm` / `rep` / `as.vector` / `t` / 命名派发表，不含 `for` / `while` 原生循环。
- 时序平滑的迭代（时间轴逐点滤波、局部加权、带状 Cholesky 等）全部在 `src/panelio_smoothing.f90` 内完成，R 层每个方法只做一次 `.Fortran` 派发；无 R 层循环、无 `switch` 分支（用命名派发表 `smoothing_handlers[[method]]`）。
- 平滑方法逐一与 R 参考实现对照验证：`ekf` / `rts` / `ewma` / `gaussian` / `moving_average` / `savgol` / `loess` / `spline` 最大误差 ≤ 6e-15；断点保留阶跃验证通过（带断点跳跃 1.709 vs 无断点被抹到 1.057）。
- 因果估计逐一与独立参考对照：`did` vs 闭式解（误差 0）、`pooled_ols` / `event_study` vs `lm` 双向固定效应（≤1e-6）、`scm`（nnls）权重用 **KKT 最优性**校验、`scm`（simplex 精确凸）权重 Σw=1 且活动集梯度 spread=0（严格最优），并验证多次调用完全确定。placebo / bootstrap / pre_fit 三个推断层均已向量化实现并通过功能验证。真实 A6 基准（40×50×9）四方法冒烟通过。
- 浮点精度：写出时数值单元格按 `CSV_WRITE_DIGITS`（17）位有效数字格式化，整数值补 `.0` 保证读回类型仍为 double。17 位有效数字足以唯一确定一个 IEEE-754 double，`1/3`、`1e-20`、`1e300`、`pi/4` 经归档后可逐位还原（往返误差为 0）。

## 测试与质量

```r
# 单元测试
Rscript -e 'testthat::test_local(".")'

# 重新生成 man/ 与 NAMESPACE
Rscript -e 'roxygen2::roxygenise(".")'

# 格式化
Rscript -e 'styler::style_dir(c("R", "tests"))'

# 静态检查（object_usage_linter 已因多文件包跨文件引用而在 .lintr 中关闭）
Rscript -e 'lintr::lint_dir()'
```

质量状态（实测）：`testthat` 全部通过（含读写闭环、预处理、断点保护、因果四方法对照与 scm KKT 校验）、`styler` 零改动、`lintr` 0 条、`R CMD check` **Status: OK**（0 ERROR / 0 WARNING / 0 NOTE）。

## 许可

MIT，声明见 `DESCRIPTION` 的 `License` 字段；`LICENSE` 为 R 标准的 DCF 许可存根。
