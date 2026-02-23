#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  has_data_table <- requireNamespace("data.table", quietly = TRUE)
}))

read_clusters_full <- function(path) {
  if (has_data_table) {
    return(data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, showProgress = FALSE))
  }
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

stop_with <- function(...) stop(paste0(...), call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
project_root <- getwd()
results_dir <- file.path(project_root, "Shlush_analysis", "results")

if (length(args) >= 1 && nzchar(args[[1]])) {
  results_dir <- args[[1]]
}

if (!dir.exists(results_dir)) {
  stop_with("results_dir 不存在: ", results_dir, "\n",
            "用法: Rscript extract_sciclone_mutation_assignment.R [results_dir]\n",
            "例如: Rscript extract_sciclone_mutation_assignment.R Shlush_analysis/results")
}

patient_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
patient_dirs <- patient_dirs[file.info(patient_dirs)$isdir %in% TRUE]

if (length(patient_dirs) == 0) {
  stop_with("在目录下未找到子文件夹: ", results_dir)
}

pick_latest_clusters_full <- function(patient_dir) {
  files <- list.files(patient_dir, pattern = "_clusters_full_.*\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files[which.max(file.info(files)$mtime)]
}

required_cols <- c("chr", "st", "Diagnosis.var", "Relapse.var", "cluster")

processed <- 0L
skipped <- 0L

for (patient_dir in patient_dirs) {
  patient_id <- basename(patient_dir)
  in_path <- pick_latest_clusters_full(patient_dir)

  if (is.na(in_path) || !file.exists(in_path)) {
    message("[跳过] 未找到 *_clusters_full_*.txt: ", patient_id)
    skipped <- skipped + 1L
    next
  }

  dt <- read_clusters_full(in_path)

  missing_required <- setdiff(required_cols, colnames(dt))
  if (length(missing_required) > 0) {
    message("[跳过] 缺少必要列: ", patient_id, " | ", paste(missing_required, collapse = ", "), " | 文件: ", basename(in_path))
    skipped <- skipped + 1L
    next
  }

  prob_cols <- grep("^cluster\\.prob\\.", colnames(dt), value = TRUE)
  # 让 prob 列按数字顺序排列（cluster.prob.1, 2, 10 ...）
  if (length(prob_cols) > 0) {
    prob_idx <- suppressWarnings(as.integer(sub("^cluster\\.prob\\.", "", prob_cols)))
    if (all(!is.na(prob_idx))) {
      prob_cols <- prob_cols[order(prob_idx)]
    }
  }

  out <- data.frame(
    Chromosome = as.character(dt[["chr"]]),
    Position = dt[["st"]],
    exist_Dx = !is.na(dt[["Diagnosis.var"]]) & (dt[["Diagnosis.var"]] != 0),
    exist_Rx = !is.na(dt[["Relapse.var"]]) & (dt[["Relapse.var"]] != 0),
    mutation_assignment = dt[["cluster"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (length(prob_cols) > 0) {
    for (pc in prob_cols) {
      k <- sub("^cluster\\.prob\\.", "", pc)
      out[[paste0("prob_cluster_", k)]] <- dt[[pc]]
    }
  }

  # 去掉 mutation_assignment = NA 的数据
  out <- out[!is.na(out[["mutation_assignment"]]), , drop = FALSE]

  # 最后检查：去掉任何包含 NA 的行
  out <- out[complete.cases(out), , drop = FALSE]

  out_path <- file.path(results_dir, paste0(patient_id, "_sciClone_mutation_assignment.csv"))
  # 不给字符串加双引号（包括 Chromosome）
  utils::write.csv(out, out_path, row.names = FALSE, quote = FALSE)
  message("[完成] ", patient_id, " -> ", out_path, " | 输入: ", basename(in_path))
  processed <- processed + 1L
}

message("全部完成: processed=", processed, ", skipped=", skipped, ", results_dir=", results_dir)

