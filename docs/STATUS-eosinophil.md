# STATUS — eosinophil / 302025 + 302043（Rebecca 合作）

更新时间：2026-09-22 20:32 PDT（302043 部分为另一 Claude Code session 按 302043_git 仓库
commit 记录更新，[verified by run]；302025 部分本次由 adipose-eos-atlas 仓库内的 Claude Code
HPCC 会话按实际 SLURM job / commit 记录更新，同样 [verified by run]，见下）

> 302025（eosinophil scRNA-seq 重跑，原始 GSE311941 atlas cohort）+ 302043（eosinophil
> CITE-seq multiome 分析，后续 RR059 验证 cohort），均为与 Rebecca 合作的同一课题下的
> 两批数据——原始 atlas cohort 与后续验证 cohort，非同一批动物（RR059 内部 CITE-seq 与
> 部分 Multiome 反应共享 mouse pools，但与原 atlas 不是同一批），最终合并为一个课题。
> 原 `STATUS-multiomics.md` 已并入本文件（同一协作，无需两套状态）。

## 现在在哪一步

- 9/18–10月中全方面 focus sprint（Godzik 的命令）；10 月中 hard deadline——湿实验 postdoc 离开后数据不再产生。
- **302025（本 session 实际执行，[verified by run]，仓库 adipose-eos-atlas，RERUN_ROOT
  `/bigdata/nairlab/shared/XQiu/EOS_project/reruns/doublet_20260922_132819`）**：
  - Raw-data 分支：FastQC（SLURM job 29008756）运行中，`RUNNING`/`0:0`，已跑 6h53m，8 个
    FASTQ 文件中 3 个完成（RR1316 R1/R2、RR0104 R1），第 4 个（RR0104 R2）进行中 ~30%。
    4 个 `cellranger multi` job（29008758 RR1316、29008759 RR0104、29008760 RR0508、
    29008761 RR0912）已提交，`PENDING`，`afterok:29008756` 依赖已用 `scontrol` 核实。
    **尚未 COMPLETED，非完成证据**（sbatch 提交成功不等于完成，按 CLAUDE.md 规则）。
  - Cell Ranger 8.0.1 + mm10-2020-A reference + probe set v1.0.1：已下载、MD5 校验、
    `cellranger --version` 验证通过，PASS。RR1316 原始数据溯源问题（=RR001 改名，非
    lane-combine 或 RR0104 合并）已用 GEO submission package 校验和确认，commit
    `e6688ab`。`04_cellranger_multi.sh` 输出路径 + 去重复 highmem job 修复，commit
    `4d2cca6`。
  - R 环境分支：已完成，**PASS**（由已批准的 PASS_WITH_WARNING 关闭；两个 provenance
    caveat 保留：sccomp 1.7.16 是"1.7"的 best-effort patch 版本选择；scSHC@2671f6ee 是
    该仓库唯一 commit 的 best-effort reconstruction，非确认与原始分析完全一致）。renv
    隔离环境（R 4.5.0，非
    drift 的 shared library）：Seurat 5.1.0、CellChat 1.6.1、monocle3 1.2.9、harmony
    1.2.1 均为精确 pin；全部 68 条 pipeline 脚本已加运行时 guard（`RERUN_R_LIB` 未设时
    行为不变向后兼容；设了但 shared library 会抢先则硬失败，PASS/FAIL 两种情况均已
    smoke test 验证）。Commits：`d3c1530`/`70a0af3`/`2e7441c`/`c8e1ebc`/`2103893`/
    `77a7a3a`/`6db1252`。全部已 push 到 origin，HPCC checkout 已 pull 同步（本地/origin/
    HPCC 三方 commit 一致，已核实）。
- **302043（本 session 实际执行，[verified by run]，仓库 302043_git，commit 见下）**：
  - Phase 1（sex-differential ATAC regulon test，register 1.1-1.5）：closed。5 个已发表
    sex regulon 里只有 Stat5a 通过全部预设标准（rank 21/816, p=0.025，naive+infected
    两组都 replicate）；chrY/Xist 阳性对照通过（chrY 2.12x p=7e-145, Xist 4.71x
    p=3.6e-114）。commit `af041f6` 等。
  - C0/MHC-II genome-wide motif discovery（register 3.3）：closed，正结果。CITE-seq→ATAC
    label transfer 通过了预设 complexity gate，未再观察到原先 8-10x 的复杂度差异；
    genome-wide chromVAR + B=500
    permutation null 后，4 个预设候选中 Cebpa/Mafb 通过 rank-stability 判据；
    genome-wide top-20 motif 全部通过 max-T FWER 校正，其中 13 个同时达到预设的
    rock-solid rank-stability tier，另外 7 个属于仍然 FWER 显著但区分度较弱的
    moderate tier。commit `39d7e38`。
  - Register 4.1/4.2：图谱论文 "Male, KO-female" C0 富集 claim 只有部分支持——KO-female
    方向在最干净的对比里是反的。commit `bc3c43d`。
  - Register 4.3（男性样本降采样检验）：closed。无男性专属测序批次；女性信号在
    3/4 配对里能在功效匹配后存活（0.92-1.00），KO_F_Inf 是例外（0.28）。verdict：原来
    男性不显著更可能是功效不足，不是没有生物学信号。commit `3519dd0`。
  - Register 2.0（Multiome ATAC donor demux）：**进行中，未完成**。cellSNP-lite
    de novo 作业（29012023-29012026）已提交运行，截至最后一次确认无报错、无类似
    CITE-seq KO_F_Inf 那种卡死迹象；vireo sweep + 性别阳性对照检验还没跑。commit
    `3519dd0`（基础设施部分）。
  - 可复现性审计（reproducibility audit）：所有 critical 已修复，`make audit` 返回 0
    critical；剩余约 62 个 major 已分类，主要为当前 HPCC 可运行但影响跨环境可移植性的
    硬编码路径，因此不阻塞本次 HPCC 分析，但不是已修复项。commit `af276c8`, `4da97a1`。
  - Eosinophil vs macrophage identity check（为图谱论文 rebuttal 找独立证据）：4 轮后
    关闭，结果混合——SiglecF/CCR3 蛋白阳性率（79-100%）干净可引用；macrophage-marker
    共表达信号 4 轮方法论迭代后仍未完全解决（v4 只在功效最强的 1/6 样本里干净地消失）。
  - **本 session 尚未 push 到 origin**（用户明确要求先不 push，等另一个并行 session 的
    `archive/` 重组落地后一起 review 再 push）。
- 两线最终合并为一个课题；合并的时间点和形式待定。

## 下一步是什么

1. Chase Rebecca：拿全 bench records / metadata（9/18 邮件已发出，若无回复要在 postdoc 走之前 chase）。
2. 在湿实验 postdoc 10 月中离开前谈定 authorship / credit。
3. **302025（[verified by run] 的具体下一步）**：等 FastQC（job 29008756）跑完 → 4 个
   cellranger multi 跑完（COMPLETED/0:0 + web-summary QC gate，非仅 sbatch 提交成功）→
   提交 `00a_doublet_detection_per_sample.py`（Python，cellrank_env）→ R 分支从
   `00_decontx_per_sample.R` 开始，须带 `RERUN_R_LIB=$RERUN_ROOT/r_env/renv_library`、
   `CMDSTAN=$RERUN_ROOT/r_env/cmdstan/cmdstan-2.40.0`、`module load gcc/12.2.0 R/4.5.0`。
   当前无阻塞，纯等计算（FastQC 卡在大 FLEX 文件 I/O，非故障——已用日志确认稳定推进）。
4. **302043（[verified by run] 的具体下一步）**：
   a. 等 register 2.0 的 cellSNP-lite 作业跑完，跑 vireo sweep + 性别阳性对照，closed
      2.0。
   b. Review 另一个并行 session 未提交的 `archive/` 仓库重组，确认不冲突。
   c. 上述都完成后，push 302043_git 到 origin。
5. 交付后先跟 Godzik 通气并获同意，10 月下旬再约 Lonardi（备 Elencho preprint 初稿 + CCTL v2 10-seed 结果 + Replogle verdict）。

## 卡点 / 在等谁

- Rebecca 回复（bench records / metadata）。
- Postdoc 10 月中离开是硬 deadline。
- 302025：无阻塞，等计算（FastQC + cellranger multi，完成时间待实际进度确认——不同
  FASTQ 文件大小/压缩吞吐差异大，无可靠估计依据；FastQC 6h53m 时 3/8 文件完成，
  进度稳定）。
- 302043 的 register 2.0 在等 HPCC cellSNP-lite 作业跑完，非阻塞性等待。

## 上次动作

- [2026-09-22 20:32 PDT] [verified by run] adipose-eos-atlas（302025）：audited HPCC
  rerun 进行中，RERUN_ROOT `doublet_20260922_132819`。Raw-data 分支：FastQC job
  29008756 运行中（6h53m，3/8 文件完成），4 个 cellranger multi job 已提交并
  dependency-gated，尚未跑。R 环境分支：PASS（两个 provenance caveat 保留），Cindy 已批准 PASS_WITH_WARNING 关闭为 PASS，68 脚本
  运行时 guard 已实现并 smoke-test 验证（commit `6db1252`），本地/origin/HPCC 三方
  commit 同步核实一致 → 下一步：等 FastQC/cellranger multi 完成后进入
  `00a_doublet_detection_per_sample.py` → `00_decontx_per_sample.R`。
- [2026-09-22 20:27 PDT] [verified by run] 302043_git：register 2.0（Multiome donor
  demux）cellSNP-lite 作业提交运行中，尚未完成 → 基础设施已建好，commit `3519dd0`
  的部分 → 下一步：等作业跑完，跑 vireo + 性别对照。
- [2026-09-22] [verified by run] 302043_git：register 4.3（男性降采样）closed，女性
  MHC-II bimodality 信号在功效匹配后 3/4 配对存活 → commit `3519dd0` → 下一步：无
  （该项已关闭）。
- [2026-09-22] [verified by run] 302043_git：C0/MHC-II genome-wide motif discovery
  permutation null（B=500）跑完，4 个预设候选中 Cebpa/Mafb 通过 rank-stability 判据，
  genome-wide top-20 全部通过 max-T FWER 校正（13 个另达 rock-solid rank-stability
  tier，7 个为 FWER 显著但区分度较弱的 moderate tier）→ commit `39d7e38` → 下一步：
  已关闭，manuscript text 已同步更新。
- [2026-09-22] [verified by run] 302043_git：reproducibility audit 所有 critical
  已修复，`make audit` 返回 0 critical → commit `af276c8`, `4da97a1` → 下一步：剩余
  ~62 个 major 已分类，为当前 HPCC 可运行但影响跨环境可移植性的硬编码路径，不阻塞
  本次 HPCC 分析，但不是已修复项。
- [2026-09-22] Cindy 确认 302025 / 302043 均为 Rebecca 合作同一批数据；原 `STATUS-multiomics.md` 并入本文件。
  **[2026-09-22 修正]**：审核后澄清，"同一批数据"表述不准确——302025 是原始 GSE311941
  atlas cohort，302043/RR059 是后续验证 cohort，RR059 内部 CITE-seq 与部分 Multiome
  反应共享 mouse pools，但与原 atlas 非同一批动物。准确表述见文件开头 blockquote。
- [2026-09-18] 给 Rebecca 的 bench-records 邮件已发出（Cindy 确认）。
- [2026-09-22] 防编造条款（FOR-claude-code.md）同样适用于本项目所有开工会话。

## Someday / out of scope

- （合并完成后，在此记录合并决议：时间、合并后的统一命名）
