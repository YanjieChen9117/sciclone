#!/usr/bin/env Rscript
# 从多样本 sciClone mutation assignment 表格和 cluster VAF 汇总表
# 计算并输出 Dx / Rx 的 cellular prevalence 表格
# 用法：在项目根目录运行，或设置 results_dir 为 Shlush_analysis/results 的路径

# ---- 参数 ----
# mutation 数量少于该比例（相对该时间点总 mutation 数）的 cluster 将被筛去
mutation_fraction_threshold <- 0.01  # 1%

# results 目录（脚本在项目根时使用 Shlush_analysis/results）
results_dir <- file.path("Shlush_analysis", "results")
if (!dir.exists(results_dir)) {
  results_dir <- "results"
}
if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

# 查找所有 patient 的 mutation assignment 文件
assignment_files <- list.files(
  results_dir,
  pattern = "^Patient-[0-9]+_sciClone_mutation_assignment\\.csv$",
  full.names = TRUE
)

if (length(assignment_files) == 0) {
  stop("No *_sciClone_mutation_assignment.csv found in ", results_dir)
}

# 从文件名提取 patient ID，例如 Patient-1_sciClone_mutation_assignment.csv -> Patient-1
patient_ids <- sub("_sciClone_mutation_assignment\\.csv$", "", basename(assignment_files))

for (i in seq_along(assignment_files)) {
  patient_id <- patient_ids[i]
  assignment_path <- assignment_files[i]
  # 对应 patient 子目录下的 means 文件，例如 results/Patient-1/Patient-1_dp25_cluster_summary_v5.means
  means_path <- file.path(results_dir, patient_id, paste0(patient_id, "_dp25_cluster_summary_v5.means"))

  if (!file.exists(means_path)) {
    message("Skip ", patient_id, ": means file not found: ", means_path)
    next
  }

  # 读取 mutation assignment
  assign_dt <- read.csv(assignment_path, stringsAsFactors = FALSE)
  if (!"mutation_assignment" %in% names(assign_dt)) {
    message("Skip ", patient_id, ": no column mutation_assignment")
    next
  }
  if (!"exist_Dx" %in% names(assign_dt) || !"exist_Rx" %in% names(assign_dt)) {
    message("Skip ", patient_id, ": missing exist_Dx or exist_Rx")
    next
  }

  # 读取 cluster VAF（tab 分隔，第一列为 cluster 名，列为 Diagnosis / Relapse）
  means_dt <- read.delim(means_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Diagnosis" %in% names(means_dt) || !"Relapse" %in% names(means_dt)) {
    message("Skip ", patient_id, ": means file missing Diagnosis or Relapse column")
    next
  }

  # 从第一列（cluster1, cluster2, ...）解析 cluster index（列名可能为空，用列索引）
  cluster_labels <- means_dt[, 1]
  means_dt$cluster_index <- as.integer(sub("^cluster", "", cluster_labels))
  means_dt <- means_dt[order(means_dt$cluster_index), ]

  # ---- Dx: exist_Dx == TRUE 的 mutation 按 cluster 计数，合并 Diagnosis VAF ----
  dx_sub <- assign_dt[toupper(as.character(assign_dt$exist_Dx)) == "TRUE", ]
  dx_counts <- as.data.frame(table(cluster_index = dx_sub$mutation_assignment), stringsAsFactors = FALSE)
  dx_counts$cluster_index <- as.integer(as.character(dx_counts$cluster_index))
  dx_prev <- merge(
    data.frame(cluster_index = means_dt$cluster_index),
    dx_counts,
    by = "cluster_index",
    all.x = TRUE
  )
  dx_prev$Freq[is.na(dx_prev$Freq)] <- 0
  names(dx_prev)[names(dx_prev) == "Freq"] <- "mutation_count"
  dx_prev$mean_VAF <- means_dt$Diagnosis[match(dx_prev$cluster_index, means_dt$cluster_index)]
  dx_prev$`2*mean_VAF` <- 2 * dx_prev$mean_VAF
  total_dx_mutations <- nrow(dx_sub)
  dx_threshold <- max(1L, ceiling(mutation_fraction_threshold * total_dx_mutations))
  # 筛去 mutation 数量少于该时间点总 mutation 数 1% 的 cluster
  dx_out <- dx_prev[dx_prev$mutation_count >= dx_threshold, , drop = FALSE]
  dx_out <- dx_out[order(dx_out$`2*mean_VAF`, decreasing = TRUE), , drop = FALSE]
  dx_out$rank <- seq_len(nrow(dx_out))
  dx_out <- dx_out[, c("rank", "cluster_index", "mutation_count", "2*mean_VAF")]
  dx_path <- file.path(results_dir, paste0(patient_id, "_Dx_cellular_prevalence.csv"))
  write.csv(dx_out, dx_path, row.names = FALSE)
  message("Written: ", dx_path)

  # ---- Rx: exist_Rx == TRUE 的 mutation 按 cluster 计数，合并 Relapse VAF ----
  rx_sub <- assign_dt[toupper(as.character(assign_dt$exist_Rx)) == "TRUE", ]
  rx_counts <- as.data.frame(table(cluster_index = rx_sub$mutation_assignment), stringsAsFactors = FALSE)
  rx_counts$cluster_index <- as.integer(as.character(rx_counts$cluster_index))
  rx_prev <- merge(
    data.frame(cluster_index = means_dt$cluster_index),
    rx_counts,
    by = "cluster_index",
    all.x = TRUE
  )
  rx_prev$Freq[is.na(rx_prev$Freq)] <- 0
  names(rx_prev)[names(rx_prev) == "Freq"] <- "mutation_count"
  rx_prev$mean_VAF <- means_dt$Relapse[match(rx_prev$cluster_index, means_dt$cluster_index)]
  rx_prev$`2*mean_VAF` <- 2 * rx_prev$mean_VAF
  total_rx_mutations <- nrow(rx_sub)
  rx_threshold <- max(1L, ceiling(mutation_fraction_threshold * total_rx_mutations))
  # 筛去 mutation 数量少于该时间点总 mutation 数 1% 的 cluster
  rx_out <- rx_prev[rx_prev$mutation_count >= rx_threshold, , drop = FALSE]
  rx_out <- rx_out[order(rx_out$`2*mean_VAF`, decreasing = TRUE), , drop = FALSE]
  rx_out$rank <- seq_len(nrow(rx_out))
  rx_out <- rx_out[, c("rank", "cluster_index", "mutation_count", "2*mean_VAF")]
  rx_path <- file.path(results_dir, paste0(patient_id, "_Rx_cellular_prevalence.csv"))
  write.csv(rx_out, rx_path, row.names = FALSE)
  message("Written: ", rx_path)
}

message("Done. Outputs are in: ", normalizePath(results_dir))
