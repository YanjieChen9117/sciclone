#!/usr/bin/env Rscript
# 从单样本 sciClone 结果计算 Dx / Rx 的 cellular prevalence 表格
# 单样本模式：Dx 与 Rx 分别来自两次独立的 sciClone 运行，各有各自的 assignment 与 cluster summary
# 输出写入 Shlush_analysis/results_single，与 compute_cellular_prevalence.R（双样本 results）对应

# ---- 参数 ----
# mutation 数量少于该比例（相对该时间点总 mutation 数）的 cluster 将被筛去
mutation_fraction_threshold <- 0.01  # 1%

results_dir <- file.path("Shlush_analysis", "results_single")
if (!dir.exists(results_dir)) {
  results_dir <- "results_single"
}
if (!dir.exists(results_dir)) {
  stop("Results directory not found: ", results_dir)
}

# 从单样本 assignment 文件名提取 patient ID
# 例如 Patient-1_diagnosis_sciClone_mutation_assignment.csv -> Patient-1
assignment_files_dx <- list.files(
  results_dir,
  pattern = "^Patient-[0-9]+_diagnosis_sciClone_mutation_assignment\\.csv$",
  full.names = TRUE
)
assignment_files_rx <- list.files(
  results_dir,
  pattern = "^Patient-[0-9]+_relapse_sciClone_mutation_assignment\\.csv$",
  full.names = TRUE
)
patient_ids_dx <- sub("_diagnosis_sciClone_mutation_assignment\\.csv$", "", basename(assignment_files_dx))
patient_ids_rx <- sub("_relapse_sciClone_mutation_assignment\\.csv$", "", basename(assignment_files_rx))
patient_ids <- sort(unique(c(patient_ids_dx, patient_ids_rx)))

if (length(patient_ids) == 0) {
  stop("No single-sample assignment files found in ", results_dir)
}

# 单时间点 CP：用该时间点的 assignment + 同时间点的 means 表
compute_cp_one <- function(assignment_path, means_path, out_path, timepoint_label) {
  if (!file.exists(assignment_path)) return(invisible(NULL))
  if (!file.exists(means_path)) {
    message("Skip ", timepoint_label, ": means not found: ", means_path)
    return(invisible(NULL))
  }

  assign_dt <- read.csv(assignment_path, stringsAsFactors = FALSE)
  if (!"mutation_assignment" %in% names(assign_dt)) {
    message("Skip ", timepoint_label, ": no column mutation_assignment in ", assignment_path)
    return(invisible(NULL))
  }

  # 单样本 assignment 中无 exist_Dx/exist_Rx，所有行为该样本的突变
  # 可选：若有 exist 列则只保留 exist==TRUE
  if ("exist" %in% names(assign_dt)) {
    assign_dt <- assign_dt[toupper(as.character(assign_dt$exist)) == "TRUE", , drop = FALSE]
  }

  means_dt <- read.delim(means_path, check.names = FALSE, stringsAsFactors = FALSE)
  # 单样本 means: 第 1 列为 cluster 标签 (cluster1, cluster2, ...)，第 2 列为该样本的 mean VAF
  if (ncol(means_dt) < 2) {
    message("Skip ", timepoint_label, ": means file has < 2 columns: ", means_path)
    return(invisible(NULL))
  }
  cluster_labels <- means_dt[, 1]
  mean_vaf <- as.numeric(means_dt[, 2])
  means_dt$cluster_index <- as.integer(sub("^cluster", "", cluster_labels))
  means_dt$mean_VAF <- mean_vaf
  means_dt <- means_dt[order(means_dt$cluster_index), ]

  # 按 cluster 计数
  counts <- as.data.frame(table(cluster_index = assign_dt$mutation_assignment), stringsAsFactors = FALSE)
  counts$cluster_index <- as.integer(as.character(counts$cluster_index))
  prev <- merge(
    data.frame(cluster_index = means_dt$cluster_index),
    counts,
    by = "cluster_index",
    all.x = TRUE
  )
  prev$Freq[is.na(prev$Freq)] <- 0
  names(prev)[names(prev) == "Freq"] <- "mutation_count"
  prev$mean_VAF <- means_dt$mean_VAF[match(prev$cluster_index, means_dt$cluster_index)]
  prev$`2*mean_VAF` <- 2 * prev$mean_VAF

  total_mutations <- nrow(assign_dt)
  threshold <- max(1L, ceiling(mutation_fraction_threshold * total_mutations))
  prev <- prev[prev$mutation_count >= threshold, , drop = FALSE]
  prev <- prev[order(prev$`2*mean_VAF`, decreasing = TRUE), , drop = FALSE]
  prev$rank <- seq_len(nrow(prev))
  prev <- prev[, c("rank", "cluster_index", "mutation_count", "2*mean_VAF")]

  write.csv(prev, out_path, row.names = FALSE)
  message("Written: ", out_path)
  invisible(NULL)
}

for (patient_id in patient_ids) {
  # Dx: 来自单样本 Diagnosis sciClone 结果
  assignment_dx <- file.path(results_dir, paste0(patient_id, "_diagnosis_sciClone_mutation_assignment.csv"))
  means_dx     <- file.path(results_dir, patient_id, paste0(patient_id, "_Diagnosis_dp25_cluster_summary_v5.means"))
  out_dx       <- file.path(results_dir, paste0(patient_id, "_Dx_cellular_prevalence.csv"))
  compute_cp_one(assignment_dx, means_dx, out_dx, paste0(patient_id, " Dx"))

  # Rx: 来自单样本 Relapse sciClone 结果
  assignment_rx <- file.path(results_dir, paste0(patient_id, "_relapse_sciClone_mutation_assignment.csv"))
  means_rx     <- file.path(results_dir, patient_id, paste0(patient_id, "_Relapse_dp25_cluster_summary_v5.means"))
  out_rx       <- file.path(results_dir, paste0(patient_id, "_Rx_cellular_prevalence.csv"))
  compute_cp_one(assignment_rx, means_rx, out_rx, paste0(patient_id, " Rx"))
}

message("Done. Outputs are in: ", normalizePath(results_dir))
