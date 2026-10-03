# real_data — 真实纵向微生物组数据集（panelio 实证验证）

本目录存放可公开获取的真实纵向微生物组数据，用于 panelio 管线的**实证/迁移/鲁棒性验证**（与 `benchmark/A1..A12` 仿真基准互补：仿真有真值表，用于 Type I / power / bias 校准；真实数据无 ground-truth 效应，用于验证管线在真实纵向菌群上跑得通、结果可解释、敏感性合理）。

## 目录结构

| 子目录 | 数据集 | 状态 |
|---|---|---|
| `diabimmune/` | DIABIMMUNE 抗生素亚研究（Yassour et al. 2016） | ✅ 已下载 |
| `ihmp_ibd/` | iHMP-IBD / IBDMDB（Lloyd-Price et al. 2019） | ✅ 已下载（16S 微生物组 + metadata） |
| `teddy/` | TEDDY（Environmental Determinants of Diabetes in the Young） | ⚠️ 受限（dbGaP 需授权） |
| `armord/` | ARMORD（抗生素观察性队列） | ⚠️ 待核实干预元数据 |

---

## DIABIMMUNE 抗生素亚研究（`diabimmune/`）

- **来源**：https://diabimmune.broadinstitute.org/diabimmune/antibiotics-cohort （Broad Institute 直链下载）
- **论文**：Yassour et al., *Sci. Transl. Med.* 2016 (343ra81)；补充材料 `8-343ra81_SM.pdf`
- **设计**：39 名芬兰婴儿，出生后第 2 月起，前三年**每月粪便 16S 采样**；抗生素暴露子研究，分 **Abx+（暴露）/ Abx−（未暴露）** 两组。
- **文件**：
  - `otu_table.filtered.2Maaslin.txt` — 16S OTU 相对丰度表（行=分类层级，列=样本 `E######_月龄`，值=相对丰度小数）
  - `all.relative_abundance.txt` — 全量相对丰度表
  - `aad0917.SuppTable1.xls` — 元数据（3 个 sheet）：
    - `General`：subject / 出生年 / 国家 / 性别 / 分娩方式 / 孕周
    - **`Antibiotics`：干预节点核心** — 每行一次抗生素治疗（Subject, Age (months), Antibiotic type, Antibiotic code, Duration, Symptoms）
    - `Early feeding`：母乳/配方喂养时间线
  - `8-343ra81_SM.pdf` — 补充材料（Methods、分组定义 Abx±）
- **干预节点（t₀）**：每个 subject **首次抗生素暴露的月龄**（取自 `Antibiotics` sheet）。暴露前采样 = pre-intervention；暴露后采样 = post-intervention。
- **treated / donor**：Abx+ 个体（有抗生素记录）= treated；Abx− 个体（`Antibiotics` sheet 无记录）= donor 池。
- **数据形态**：相对丰度（每样本和≈1），**兼容 panelio 输入契约**（可直接构 N×D×T 数组，无需转计数）。

**最适合合成控制**：干预前后（抗生素）结构最干净，有明确未暴露供体。

---

## iHMP-IBD / IBDMDB（`ihmp_ibd/`）

- **来源**：https://ibdmdb.org （HMP2 数据，Globus 数据端点直链）
- **论文**：Lloyd-Price et al., *Nature* 2019 (Multi-omics of the gut microbial ecosystem in IBD)
- **设计**：约 130 名 CD/UC/非IBD 个体，**一年纵向**、多组学（16S / 宏基因组 / 转录 / 代谢）；炎症性肠病的 flare（疾病活动）作为时变事件。
- **文件**：
  - `hmp2_metadata_2018-08-20.csv` — 全样本元数据（490 列，5533 行）：`Participant ID` 个体、`week_num` 时间、`diagnosis` 诊断（CD/UC/nonIBD）、**`is_inflamed`（炎症=flare 干预节点）**、`IntervalName` 访视等
  - `taxonomic_profiles.tsv.gz` / `.biom.gz` — **合并 16S 相对丰度谱表**（180 样本 × 980 OTU，相对丰度）
  - `16S_tax_profiles/` — 178 个 per-sample 16S taxonomy `.biom`（biopsy_16S）
- **干预节点（t₀）**：个体**首次 flare（`is_inflamed` 从 No→Yes）** 的周数；或诊断分组（CD/UC=treated vs 非IBD=donor）。
- **treated / donor**：经历 flare 的 IBD 个体 = treated；全程未 flare 的非IBD/缓解个体 = donor。
- **数据形态**：biom/tsv（相对丰度），兼容 panelio 输入契约。
- **说明**：16S `.biom` 为组织活检样本（biopsy_16S），合并表覆盖含粪便的样本集；用于合成控制优先用合并表 + `stool_16S` 样本筛选。

---

## TEDDY（`teddy/`）— 受限

- 纵向出生队列（美国/欧洲，遗传高风险 1 型糖尿病），16S + 宏基因组（887 名个体、1 万余粪便样本）。
- 干预节点最干净：**胰岛素自身抗体 seroconversion**（自然事件），有配比非阳转对照。
- **数据在 dbGaP（NIDDK Central Repository，study 24）**，需授权申请，**无法公开下载**。获得访问授权后方可在此目录落地数据。

---

## ARMORD（`armord/`）— 待核实

- 抗生素暴露观察性队列，**ENA PRJEB86785** 开放（processed microbiome + 元数据在补充材料）。
- 干预节点（抗生素暴露时间点）元数据需进一步核实后再下载；当前目录为空。

---

## 使用提示（如何喂进 panelio 管线）

1. **真实数据无 ground-truth 效应** → 只能做**迁移/鲁棒性/敏感性验证**和**效应发现**，不能做 bias/power 校准（那是 `benchmark/A1..A12` 仿真基准的活，有 `benchmark/truth/truth_table.csv`）。
2. **输入适配（管线零改动）**：计数/OTU 表 →（若为计数）**逐样本和转 1 的相对丰度** → 构 `N×D×T` 数组 + subject / time / intervention 元数据 → `preprocess()`（滤过 / 零填充 / CLR 均一化）→ `panel_causal_estimate()`（di / event / pooled / scm）。
3. 各数据集均以相对丰度（Σ≈1）为输入契约，与 `panelio` 原生组成型面板一致，无需修改包代码。

---

## 已完成的 DIABIMMUNE panelio 适配（实证示例）

已把 DIABIMMUNE 抗生素亚研究整理为 panelio 原生可读格式并**跑通整条管线**（读入 → 预处理 → 因果估计），作为真实数据的迁移/实证验证示例。

### 平衡面板构造（关键设计）

真实纵向采样时间不均衡（40 个体各 19–36 次采样、月龄 1.3–36.3），而 panelio 要求**平衡面板**（所有个体在每个时间层都有观测）且**整个面板只能有一个干预断点**。解决方案：

- **时间层**：取全体个体都采样的 **6 个公共月龄 bin**（3, 6, 9, 12, 15, 18 月，bin 宽 3 月）；bin 内多个采样取**均值**，再归一化为组成。
- **干预断点（t₀）**：**6 月龄**（`time_006` 为 `is_intervention_timepoint`）。
- **treated / donor**：`treated` = 首次抗生素 ≤ 6 月龄（8 个体）；`control` = 其余（31 个体，含未暴露 + 晚期暴露）。
- **feature**：属级 142 个（`|g__` 且非 `|s__`），**列和 = 1.000**，天然满足组成约束。

### 产物

| 文件 | 说明 |
|---|---|
| `diabimmune/panel_long.csv` | 平衡面板长表（39 obs × 6 层 × 142 属，33228 行；列 `observation,time_identifier,elapsed_time,feature,abundance_value,label,is_intervention`） |
| `diabimmune/diabimmune_abx_v1.tar.gz` | **panelio 原生归档**（`write_panel_archive` 打包，`read_panel_table` 读回 round-trip 一致） |
| `diabimmune/panel_subject_meta.csv` | 个体级元数据（label / 首次抗生素月龄 / 采样数） |
| `results/real_diabimmune_{did,pooled,event,scm}.csv` | 四种估计器的逐属效应（estimate / SE / t / p / BH-FDR / pre-post 均值 / effect_range）；scm 含 placebo/pre_fit 推断 |

### 实证结果摘要（CLR 预处理后）

| 方法 | 显著属 (fdr<0.05) | 备注 |
|---|---|---|
| `did` | 0 | pre 仅 1 层，双固定效应保守 |
| `event_study` | 0 | 同上 |
| `pooled_ols` | 2 | |
| `scm` (默认 nnls) | 50 | 合成控制最灵敏；如 Clostridiaceae_unclassified、Sutterella 等显著上升 |

- 结论：**合成控制在真实纵向菌群上能稳定发现抗生素暴露相关的组成变化**；did/event 因 pre 期仅 1 层而保守——这是真实数据采样结构的固有局限，论文实证节需如实说明。
- **局限**：真实数据无真值，此结果为**效应发现/迁移验证**，非因果校准；8 treated 个体功效有限。



---

## 已完成的 iHMP-IBD panelio 适配（MGX 宏基因组，实证示例）

### 数据可行性核实（重要）

- **合并 16S 谱表实际全是 biopsy_16S**（178 样本），且 `is_inflamed` 几乎全为 Yes（169/178）——活检只在发炎时取，**无 pre-intervention 时间点**，无法做 flare 干预合成控制。
- **stool_16S / MGX 的 `is_inflamed` 全空**（该列仅在 biopsy 有值），逐周 flare 标记不可用。
- 因此采用**诊断分组替代**：用 **MGX（粪便宏基因组）谱表**做 IBD（CD/UC）vs non-IBD 的合成控制比较。

### 平衡面板构造

- **feature**：MGX 属级 199 个（列和≈1）。
- **时间层**：4 周 bin，取全体参与者覆盖最多的平衡子面板 → **8 层**（周 0,4,8,12,16,24,32,36）；bin 内均值归一化为组成。
- **干预断点（t₀）**：**周 24**（`time_024`）。
- **treated / donor**：`treated` = IBD（CD 19 + UC 11 = 30）；`control` = non-IBD（11）。

### 产物

| 文件 | 说明 |
|---|---|
| `ihmp_ibd/panel_long_ihmp.csv` | 平衡面板长表（41 obs × 8 层 × 199 属，65272 行） |
| `ihmp_ibd/ihmp_ibd_mgx_v1.tar.gz` | **panelio 原生归档**（round-trip 一致） |
| `ihmp_ibd/panel_subject_meta_ihmp.csv` | 个体级元数据（label / diagnosis / 采样数） |
| `ihmp_ibd/mgx_taxonomic_profiles.tsv.gz` | MGX 宏基因组谱表（1638 样本 × 1479 分类层级） |
| `results/real_ihmp_{did,pooled,event,scm}.csv` | 四种估计器逐属效应 |

### 实证结果摘要（CLR 预处理后，断点=周 24）

| 方法 | 显著属 (fdr<0.05) | 备注 |
|---|---|---|
| `scm` | **12** | Bacteroidia ↓、Gammaproteobacteria ↓（符合 IBD 文献） |
| `pooled_ols` | 1 | |
| `did` / `event_study` | 0 | pre 期短，固定效应保守 |

### 论文更新

两个真实数据集的实证已写入 `paper/panelio.tex` 新增节 `Real-data empirical validation`（§sec:real），含平衡面板构造方法、结果表（`tab:real`）、生物学解读与局限说明。
