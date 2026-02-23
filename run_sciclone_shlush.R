#!/usr/bin/env Rscript
################################################################################
# Shlush 队列 SciClone 克隆演化分析脚本
# 
# 主要功能：
#   - 批量遍历 Shlush_analysis 目录中的所有患者样本
#   - 直接读取已经预处理好的 SNV / CNV 表格（CSV）
#   - 调用本仓库中的 sciClone 核心函数完成克隆聚类
#   - 为每个患者单独建立结果文件夹，输出所有 SciClone 相关结果与可视化图
# 重要约定：
#   - 不再在本脚本中做原始 VCF / CNV 的预处理，仅做“读取 + 调用 sciClone”
#   - 处理状态全部记录在 Shlush_analysis/processing_status.csv 中，可断点续跑
################################################################################

# 加载依赖的 R 包（仅用于数据操作与绘图，克隆聚类逻辑来自本仓库 R/ 目录）
cat("正在加载所需 R 包...\n")
suppressPackageStartupMessages({
  library(reshape2)  # 用于数据重塑（melt/cast 操作）
  library(limma)     # 提供 strsplit2 函数（与示例脚本保持一致）
  library(ggplot2)   # 用于绘制 2D 散点图等结果图像
})

# 固定工作目录到项目根目录，避免从不同路径运行脚本导致相对路径失效
project_root <- "/Users/yanjiechen/Documents/Github/sciclone"
setwd(project_root)

# 直接加载本仓库中的 sciClone 源文件
# 注意：这里不使用 library(sciClone)，而是显式 source，本仓库即是一个 R 包的源码结构
cat("从本仓库加载 sciClone 函数...\n")
source(file.path(project_root, "R", "object.R"))      # 先加载 object.R，定义 initScClass 函数
initScClass()                                          # 显式初始化 scObject 类（相当于包加载时的 .onLoad）
source(file.path(project_root, "R", "sciClone.R"))
source(file.path(project_root, "R", "clustering.R"))
source(file.path(project_root, "R", "plots.R"))

# 定义数据与结果目录（均位于 Shlush_analysis 子目录下）
# mutation_dir: 已预处理好的 SNV 表格（每行一个突变）
# cn_dir      : 已预处理好的 CNV 表格（每行一个 CN 区段）
# results_dir : 每个样本一个子文件夹，存放 SciClone 输出的所有结果
# data_dir    : Shlush 分析相关的通用文件（如驱动基因列表）
mutation_dir <- file.path(project_root, "Shlush_analysis", "mutation_table")
cn_dir <- file.path(project_root, "Shlush_analysis", "copy_number_table")
results_dir <- file.path(project_root, "Shlush_analysis", "results")
data_dir <- file.path(project_root, "Shlush_analysis")

# 基础路径检查：如果输入目录不存在，直接报错退出
if (!dir.exists(mutation_dir)) {
  stop(sprintf("Mutation table directory does not exist: %s", mutation_dir))
}
if (!dir.exists(cn_dir)) {
  stop(sprintf("Copy number table directory does not exist: %s", cn_dir))
}

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 定义克隆簇的配色方案（与示例脚本保持一致，最多 9 个簇）
COLORS <- c(
  "#ABD9F6", "#A8DC69", "#D9BDE3", "#FFA50A", "#FC6F6F",
  "#3B97CB", "#575757", "#30AE29", "#FF0000"
)

################################################################################
# 函数：clonalEvolution
# 作用：对某一位患者（诊断期 + 复发期）进行完整的克隆演化分析，并写出所有 SciClone 结果文件
#
# 参数说明：
#   Tumor_mut   : 诊断期突变表（数据框），列格式必须为：
#                  contig, position, norm_read, tum_read, tum_vaf, gene_name, sample, variant
#                ——前 5 列会按 sciClone 约定重命名为：
#                  chr, st, ref, var, vaf（vaf 为 0–100 的百分比）
#   Relapse_mut : 复发期突变表，列格式与 Tumor_mut 相同
#   Tumor_cnv   : 诊断期 CNV 表，列格式：
#                  contig, start, stop, cn, sample, type
#                ——前 4 列会传入 sciClone 作为 copyNumberCalls
#   Relapse_cnv : 复发期 CNV 表，列格式与 Tumor_cnv 相同
#   SampName    : 患者 ID（例如 "Patient-4"），仅用于输出文件命名
#   MIN_DEPTH   : 最小测序深度阈值（sciClone 的过滤参数，默认 50）
#   VER         : 输出版本号，用于区分不同参数设置产生的结果
#   output_dir  : 本患者所有输出文件所在的目录（通常为 Shlush_analysis/results/<SampleID>/）
################################################################################
clonalEvolution <- function(Tumor_mut,
                            Relapse_mut,
                            Tumor_cnv,
                            Relapse_cnv,
                            SampName,
                            MIN_DEPTH = 50,
                            VER = 1,
                            output_dir = ".") {
  factorNames <- c("Diagnosis", "Relapse")

  # 提取 LOH（Loss of Heterozygosity，杂合性缺失）区域
  # 在 LOH 区域中，等位基因频率会由于拷贝数不平衡而被严重扭曲，
  # 若不排除这些区域，克隆聚类会受到很大干扰，因此在 sciClone 中作为排除区域传入
  loh_regions_tum <- Tumor_cnv[Tumor_cnv$type == "LOH", 1:3]
  loh_regions_rel <- Relapse_cnv[Relapse_cnv$type == "LOH", 1:3]

  # ========== 使用 sciClone 进行克隆聚类 ==========
  # sciClone 会自动：
  #   1）按 chr / pos 合并多个时间点的 SNV
  #   2）在 CN=2 且测序深度足够的位点上进行聚类
  #   3）输出每个突变的聚类编号与后验概率
  # vafs / copyNumberCalls 的列顺序非常关键，这里直接子集前几列以匹配 sciClone 约定
  sc <- sciClone(
    vafs = list(Tumor_mut[, 1:5], Relapse_mut[, 1:5]),              # 两个时间点的 VAF 数据
    copyNumberCalls = list(Tumor_cnv[, 1:4], Relapse_cnv[, 1:4]),   # 两个时间点的 CNV 数据
    regionsToExclude = list(loh_regions_tum, loh_regions_rel),      # 需要排除的 LOH 区域
    annotation = Tumor_mut[Tumor_mut$gene_name != "", c(1, 2, 6)],  # (chr, pos, gene_name) 注释信息
    sampleNames = factorNames,          # 样本名称：诊断期 / 复发期
    useSexChrs = TRUE,                  # 是否使用性染色体上的突变
    minimumDepth = MIN_DEPTH,           # 最小测序深度
    verbose = 0,                        # 关闭详细日志输出（由本脚本控制日志）
    doClusteringAlongMargins = TRUE,    # 同时对边缘（低 VAF 区域）做聚类
    plotIntermediateResults = 0,        # 不保存中间调试图
    cnCallsAreLog2 = FALSE,             # CNV 输入已经是绝对拷贝数，不是 log2 比值
    maximumClusters = 10                 # 最大允许的克隆簇数量
  )

  # 写出完整的聚类结果表
  # 表中包含：
  #   - 各时间点的 VAF / read depth / CN 等信息
  #   - 每个位点的聚类编号（cluster）与聚类概率（cluster.prob）
  #   - 若提供 annotation，则会附带基因名等注释
  writeClusterTable(
    sc,
    file.path(output_dir, paste0(SampName, "_dp", MIN_DEPTH, "_clusters_full_v", VER, ".txt"))
  )

  # 计算并写出“按簇汇总”的简化结果表
  # 每个簇包含：
  #   - Diagnosis 列：该簇在诊断期样本中的平均 VAF
  #   - Relapse   列：该簇在复发期样本中的平均 VAF
  #   - Nmut      列：该簇包含的突变数量
  Nclust <- max(sc@vafs.merged$cluster, na.rm = TRUE)
  ClustStats <- array(0, dim = c(Nclust, length(factorNames) + 1))
  
  for (k in seq_len(Nclust)) {
    # Find all mutations belonging to cluster k
    idx <- sc@vafs.merged$cluster == k & !is.na(sc@vafs.merged$cluster)
    
    # Calculate average VAF at both time points and total number of mutations for this cluster
    ClustStats[k, ] <- c(
      colMeans(
        sc@vafs.merged[idx, paste0(factorNames, ".vaf")],
        na.rm = TRUE
      ),
      Nmut = sum(idx)
    )
  }
  
  colnames(ClustStats) <- c(factorNames, "Nmut")
  ClustStatsDF <- data.frame(cluster = seq_len(Nclust), round(ClustStats, 3))
  write.table(
    ClustStatsDF,
    file = file.path(output_dir, paste0(SampName, "_dp", MIN_DEPTH, "_clusters_v", VER, ".txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  # 写出 cluster summary 表：
  #   *_cluster_summary_vX.means : 各簇在每个样本中的 VAF 均值
  #   *_cluster_summary_vX.lower : VAF 均值的下置信区间
  #   *_cluster_summary_vX.upper : VAF 均值的上置信区间
  writeClusterSummaryTable(
    sc,
    file.path(output_dir, paste0(SampName, "_dp", MIN_DEPTH, "_cluster_summary_v", VER))
  )

  # 提取 sciClone 合并后的长表（所有样本 / 簇信息已经 merge 完成）
  DTA <- sc@vafs.merged
  
  # 标记哪些突变位点落在癌症驱动基因上
  # 这些突变在后续可视化中会被额外高亮和加上文本标签
  DTA$iscancergene <- FALSE

  # ========== 绘制定制的 2D VAF 散点图 ==========
  # 横轴：诊断期 VAF（Diagnosis.vaf）
  # 纵轴：复发期 VAF（Relapse.vaf）
  # 颜色：克隆簇编号
  # 额外：对于驱动基因所在的突变，叠加黑色描边 + 文本标签
  maxX <- max(DTA$Diagnosis.vaf, na.rm = TRUE)
  maxY <- max(DTA$Relapse.vaf, na.rm = TRUE)
  DTA$cluster <- factor(DTA$cluster, levels = as.character(1:10))
  
  # Filter mutations in cancer driver genes
  DTASel <- DTA[DTA$iscancergene & !is.na(DTA$cluster) & DTA$cluster != "0", ]
  
  # 散点图主体：所有聚类成功的突变位点
  p <- ggplot(DTA[!is.na(DTA$cluster), ],
              aes(Diagnosis.vaf, Relapse.vaf, color = cluster)) +
    geom_point(alpha = 0.3) +  # 所有点设置一定透明度，方便观察重叠
    theme_bw() +
    scale_color_manual(values = COLORS) +  # 使用预定义的配色方案
    geom_point(
      data = DTASel,
      aes(Diagnosis.vaf, Relapse.vaf),
      colour = "black",
      size = 3
    ) +  # 驱动基因突变：先画一层黑色实心点作为描边
    geom_point(
      data = DTASel,
      aes(x = Diagnosis.vaf, y = Relapse.vaf, color = cluster),
      alpha = 0.5
    ) +  # 再叠加一层按簇着色的点，使其既有颜色又有描边
    geom_text(
      data = DTASel,
      aes(
        Diagnosis.vaf + maxX / 18,
        Relapse.vaf + maxY / 35,
        label = gene_name
      ),
      color = "black",
      angle = 15,
      size = 2.5
    ) +  # 添加驱动基因名称标签（略微平移避免遮挡点）
    xlab("Diagnosis VAF") +
    ylab("Relapse VAF")
  
  ggsave(
    p,
    file = file.path(output_dir, paste0(SampName, "_dp", MIN_DEPTH, "_clusters.2d_custom_v", VER, ".png")),
    scale = 2,
    width = 2.5,
    height = 2.2,
    dpi = 300
  )

  # 生成标准 SciClone 自带的可视化结果
  # 1D 图：展示各拷贝数状态下 VAF 的分布及簇结构
  tryCatch({
    plot_file <- file.path(output_dir, "plot_1d.pdf")
    sc.plot1d(sc, outputFile = plot_file)
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to generate 1D plot: %s\n", e$message))
  })

  # 2D 图：标准的二维克隆结构图（与 README 中示例一致）
  tryCatch({
    plot_2d_file <- file.path(output_dir, "plot_2d_comparison.pdf")
    sc.plot2d(sc, outputFile = plot_2d_file)
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to generate 2D plot: %s\n", e$message))
  })

  # 2D + margins 图：带边缘密度 / 聚类的二维图（若 marginalClust 信息可用）
  tryCatch({
    # Check if marginalClust exists and is valid
    if (is.null(sc@marginalClust) || length(sc@marginalClust) == 0) {
      cat("WARNING: marginalClust is NULL or empty, skipping plot_2d_with_margins\n")
    } else {
      # Check if all required elements exist
      valid_marginal <- TRUE
      for (d in seq_along(sc@marginalClust)) {
        if (is.null(sc@marginalClust[[d]]) || 
            !is.list(sc@marginalClust[[d]]) ||
            is.null(sc@marginalClust[[d]]$fit.x) ||
            is.null(sc@marginalClust[[d]]$fit.y)) {
          valid_marginal <- FALSE
          break
        }
      }
      
      if (valid_marginal) {
        plot_2d_margins_file <- file.path(output_dir, "plot_2d_with_margins.pdf")
        sc.plot2dWithMargins(sc, outputFile = plot_2d_margins_file)
      } else {
        cat("WARNING: marginalClust structure is invalid, skipping plot_2d_with_margins\n")
      }
    }
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to generate 2D plot with margins: %s\n", e$message))
  })

  return(sc)
}

# 工具函数：从突变表目录中推断所有样本 ID
# 约定的文件命名格式（新）：
#   <SampleID>_Dx_mutation_table.csv   （诊断期）
#   <SampleID>_Rx_mutation_table.csv   （复发期）
# 例如：
#   Patient-4_Dx_mutation_table.csv  => SampleID = "Patient-4"
################################################################################
find_all_samples <- function(mutation_dir) {
  # 筛选出所有“诊断期突变表”的文件名（Dx = Diagnosis）
  mut_files <- list.files(mutation_dir, pattern = "_Dx_mutation_table\\.csv$",
                          full.names = FALSE)

  # 根据前缀部分提取样本 ID（去掉统一的后缀）
  sample_ids <- gsub("_Dx_mutation_table\\.csv$", "", mut_files)

  return(sort(sample_ids))
}

# 函数：process_sample
# 作用：对单个样本（Patient-X）执行完整的“读取数据 + 克隆分析 + 结果输出”流程，
#       并返回给主程序一个简洁的状态列表，便于写入 processing_status.csv
#
# 返回值：list(
#   processed      : 是否尝试过处理该样本（逻辑值），本脚本固定为 TRUE
#   status         : "success" / "failed"
#   error_message  : 若失败，则包含具体的错误信息字符串；成功则为空字符串
# )
################################################################################
process_sample <- function(sample_id, mutation_dir, cn_dir, results_dir) {
  cat("\n")
  cat("========================================\n")
  cat(sprintf("Processing sample: %s\n", sample_id))
  cat("========================================\n\n")
  
  # 为该样本拼接所有输入文件路径（新命名：Dx/Rx + mutation_table / copy_number_table）
  mut_diagnosis_file <- file.path(mutation_dir,
                                   paste0(sample_id, "_Dx_mutation_table.csv"))
  mut_relapse_file <- file.path(mutation_dir,
                                 paste0(sample_id, "_Rx_mutation_table.csv"))
  cn_dx_file <- file.path(cn_dir, paste0(sample_id, "_Dx_copy_number_table.csv"))
  cn_rx_file <- file.path(cn_dir, paste0(sample_id, "_Rx_copy_number_table.csv"))
  
  # 检查四个关键输入文件是否齐全
  required_files <- list(
    "Diagnosis mutation table" = mut_diagnosis_file,
    "Relapse mutation table" = mut_relapse_file,
    "Dx copy number table" = cn_dx_file,
    "Rx copy number table" = cn_rx_file
  )
  
  missing_files <- c()
  for (file_type in names(required_files)) {
    if (!file.exists(required_files[[file_type]])) {
      missing_files <- c(missing_files, 
                        sprintf("%s: %s", file_type, required_files[[file_type]]))
    }
  }
  
  if (length(missing_files) > 0) {
    error_msg <- sprintf("Missing required files for sample %s:\n  %s", 
                        sample_id, paste(missing_files, collapse = "\n  "))
    cat(sprintf("ERROR: %s\n", error_msg))
    return(list(processed = TRUE, status = "failed", error_message = error_msg))
  }
  
  # 为该样本创建独立的结果目录（若不存在则新建）
  output_dir <- file.path(results_dir, sample_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 读取诊断期 / 复发期突变表
  # 注意：此处假定表格已经在外部脚本中预处理好，满足 sciClone 所需的列结构
  cat("\n--- Step 1: Reading mutation tables ---\n")
  Tumor_mut <- NULL
  Relapse_mut <- NULL
  tryCatch({
    Tumor_mut <- read.csv(mut_diagnosis_file, stringsAsFactors = FALSE)
    Relapse_mut <- read.csv(mut_relapse_file, stringsAsFactors = FALSE)
    
    # 去掉 contig 列中的 "chr" 前缀（与 sciClone 示例保持一致）
    Tumor_mut$contig <- gsub("^chr", "", Tumor_mut$contig)
    Relapse_mut$contig <- gsub("^chr", "", Relapse_mut$contig)
    
    cat(sprintf("  Diagnosis: %d mutations\n", nrow(Tumor_mut)))
    cat(sprintf("  Relapse: %d mutations\n", nrow(Relapse_mut)))
  }, error = function(e) {
    error_msg <- sprintf("Failed to read mutation tables: %s", e$message)
    cat(sprintf("ERROR: %s\n", error_msg))
    return(list(processed = TRUE, status = "failed", error_message = error_msg))
  })
  
  # 读取诊断期 / 复发期 CNV 表
  cat("\n--- Step 2: Reading copy number tables ---\n")
  Tumor_cnv <- NULL
  Relapse_cnv <- NULL
  tryCatch({
    Tumor_cnv <- read.csv(cn_dx_file, stringsAsFactors = FALSE)
    Relapse_cnv <- read.csv(cn_rx_file, stringsAsFactors = FALSE)
    
    # 同样去掉 CNV 中 contig 的 "chr" 前缀
    Tumor_cnv$contig <- gsub("^chr", "", Tumor_cnv$contig)
    Relapse_cnv$contig <- gsub("^chr", "", Relapse_cnv$contig)
    
    cat(sprintf("  Diagnosis: %d CNV segments\n", nrow(Tumor_cnv)))
    cat(sprintf("  Relapse: %d CNV segments\n", nrow(Relapse_cnv)))
  }, error = function(e) {
    error_msg <- sprintf("Failed to read copy number tables: %s", e$message)
    cat(sprintf("ERROR: %s\n", error_msg))
    return(list(processed = TRUE, status = "failed", error_message = error_msg))
  })
  
  # 调用 clonalEvolution 函数执行克隆演化分析（核心逻辑与示例脚本对齐）
  cat("\n--- Step 3: Running clonal evolution analysis ---\n")
  cat("Running sciClone with the following parameters:\n")
  cat("  - Number of samples: 2\n")
  cat("  - Sample names: Diagnosis, Relapse\n")
  cat("  - Minimum depth: 25\n")
  cat("  - Maximum clusters: 5\n\n")
  
  sc_result <- tryCatch({
    sc <- clonalEvolution(
      Tumor_mut,
      Relapse_mut,
      Tumor_cnv,
      Relapse_cnv,
      sample_id,
      MIN_DEPTH = 25,
      VER = 5,
      output_dir = output_dir
    )
    list(sc = sc, error = NULL)
  }, error = function(e) {
    error_msg <- sprintf("Clonal evolution analysis failed: %s", e$message)
    cat(sprintf("\n!!! ERROR during clonal evolution analysis for sample %s !!!\n", sample_id))
    cat(sprintf("Error message: %s\n", e$message))
    return(list(sc = NULL, error = error_msg))
  })
  
  # 检查 clonalEvolution 是否返回错误信息
  if (!is.null(sc_result$error)) {
    cat(sprintf("Failed to run clonal evolution analysis for sample %s\n", sample_id))
    return(list(processed = TRUE, status = "failed", error_message = sc_result$error))
  }
  
  sc <- sc_result$sc
  if (is.null(sc)) {
    error_msg <- "Clonal evolution analysis returned NULL (unknown error)"
    cat(sprintf("Failed to run clonal evolution analysis for sample %s\n", sample_id))
    return(list(processed = TRUE, status = "failed", error_message = error_msg))
  }
  
  # 确认返回对象是否为合法的 sciClone scObject
  if (!inherits(sc, "scObject")) {
    error_msg <- sprintf("Analysis returned invalid object type: %s", class(sc)[1])
    cat(sprintf("Failed to run clonal evolution analysis for sample %s\n", sample_id))
    return(list(processed = TRUE, status = "failed", error_message = error_msg))
  }
  
  cat("\n========================================\n")
  cat(sprintf("Analysis completed successfully for sample %s!\n", sample_id))
  cat("========================================\n\n")
  
  # 打印本样本的简要统计信息，方便在命令行快速浏览结果
  cat("Summary of results:\n")
  tryCatch({
    cat(sprintf("  - Total variants analyzed: %d\n", nrow(sc@vafs.merged)))
    
    if (!is.null(sc@clust) && is.list(sc@clust) && !is.null(sc@clust$cluster.assignments)) {
      n_clusters <- max(sc@clust$cluster.assignments, na.rm = TRUE)
      cat(sprintf("  - Number of clusters identified: %d\n", n_clusters))
      
      # 打印每个簇中包含的突变数量（忽略 outlier 簇 0）
      cluster_counts <- table(sc@clust$cluster.assignments)
      cat("\n  Cluster sizes:\n")
      for (clust_id in names(cluster_counts)) {
        if (clust_id != "0") {  # 0 means outlier
          cat(sprintf("    Cluster %s: %d variants\n", clust_id, cluster_counts[clust_id]))
        }
      }
      
      outliers <- sum(sc@clust$cluster.assignments == 0, na.rm = TRUE)
      if (outliers > 0) {
        cat(sprintf("    Outliers: %d variants\n", outliers))
      }
    } else {
      cat("  - No clusters identified\n")
    }
  }, error = function(e) {
    cat(sprintf("  WARNING: Failed to print summary statistics: %s\n", e$message))
  })
  
  cat(sprintf("\nAll outputs saved to: %s/\n", output_dir))
  
  return(list(processed = TRUE, status = "success", error_message = ""))
}

################################################################################
# 主程序入口：批量遍历所有样本并调度 process_sample
################################################################################

cat("\n========================================\n")
cat("SciClone Analysis for Shlush (批量克隆演化分析)\n")
cat("========================================\n\n")

# 1）扫描突变表目录，自动推断所有待分析的样本 ID
cat(sprintf("扫描突变表目录以获取样本列表：%s\n", mutation_dir))
sample_ids <- find_all_samples(mutation_dir)
cat(sprintf("\n共发现 %d 个样本待处理：%s\n", 
            length(sample_ids), paste(sample_ids, collapse = ", ")))

# 2）准备 / 读取 processing_status 文件，实现断点续跑
#    - Sample_ID     : 患者 ID
#    - processed     : 是否已经跑过（TRUE / FALSE）
#    - status        : "success" / "failed"
#    - error_message : 失败原因（若有）
#    - timestamp     : 最后一次更新该样本状态的时间
status_file <- file.path(project_root, "Shlush_analysis", "processing_status.csv")

# 若已存在历史的 processing_status 文件，则读取并跳过已经完成的样本
if (file.exists(status_file)) {
  cat(sprintf("\nReading processing status file: %s\n", status_file))
  status_df <- read.csv(status_file, stringsAsFactors = FALSE)
  
  # 找出 processed == TRUE 的样本，视为已经完成，不再重复运行
  processed_samples <- unique(status_df$Sample_ID[status_df$processed == TRUE | 
                                                   status_df$processed == "TRUE" | 
                                                   status_df$processed == "true"])
  if (length(processed_samples) > 0) {
    cat(sprintf("检测到 %d 个样本已处理完成，跳过：%s\n", 
                length(processed_samples), paste(processed_samples, collapse = ", ")))
    sample_ids <- setdiff(sample_ids, processed_samples)
    cat(sprintf("剩余待处理样本数量：%d\n", length(sample_ids)))
  }
} else {
  # 若不存在，则新建一个空的 processing_status 文件
  cat("未发现 processing_status.csv，创建新的状态跟踪文件\n")
  status_df <- data.frame(
    Sample_ID = character(0),
    processed = logical(0),
    status = character(0),
    error_message = character(0),
    timestamp = character(0),
    stringsAsFactors = FALSE
  )
  write.csv(status_df, status_file, row.names = FALSE)
}

if (length(sample_ids) == 0) {
  cat("\n所有样本均已处理完成，无需再次运行。\n")
  quit(status = 0)
}

# 3）逐个样本调用 process_sample，并在每次完成后更新 processing_status 文件
success_count <- 0
fail_count <- 0
failed_samples <- c()

for (sample_id in sample_ids) {
  result <- process_sample(sample_id, mutation_dir, cn_dir, results_dir)
  
  # Update status file
  # 将本次运行的结果写回 processing_status：
  #   - 若样本之前已有记录，则覆盖状态
  #   - 若是新样本，则追加一行
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Check if this sample already exists in status_df
  existing_idx <- which(status_df$Sample_ID == sample_id)
  
  if (length(existing_idx) > 0) {
    # Update existing record
    status_df$processed[existing_idx] <- result$processed
    status_df$status[existing_idx] <- result$status
    status_df$error_message[existing_idx] <- result$error_message
    status_df$timestamp[existing_idx] <- timestamp
  } else {
    # Add new record
    new_row <- data.frame(
      Sample_ID = sample_id,
      processed = result$processed,
      status = result$status,
      error_message = result$error_message,
      timestamp = timestamp,
      stringsAsFactors = FALSE
    )
    status_df <- rbind(status_df, new_row)
  }
  
  # Save status file after each sample
  write.csv(status_df, status_file, row.names = FALSE)
  
  if (result$status == "success") {
    success_count <- success_count + 1
    cat(sprintf("\n状态文件已更新：%s 标记为 success\n", sample_id))
  } else {
    fail_count <- fail_count + 1
    failed_samples <- c(failed_samples, sample_id)
    cat(sprintf("\n状态文件已更新：%s 标记为 failed，原因：%s\n", 
                sample_id, result$error_message))
  }
}

# 4）最终汇总打印
cat("\n")
cat("========================================\n")
cat("运行总结\n")
cat("========================================\n")
cat(sprintf("本次运行共处理样本数: %d\n", length(sample_ids)))
cat(sprintf("成功样本数          : %d\n", success_count))
cat(sprintf("失败样本数          : %d\n", fail_count))

if (length(failed_samples) > 0) {
  cat(sprintf("\n失败样本列表: %s\n", paste(failed_samples, collapse = ", ")))
}

cat("\n=== 脚本运行结束 ===\n")
