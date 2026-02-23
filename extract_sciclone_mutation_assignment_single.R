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
results_dir <- file.path(project_root, "Shlush_analysis", "results_single")

if (length(args) >= 1 && nzchar(args[[1]])) {
  results_dir <- args[[1]]
}

if (!dir.exists(results_dir)) {
  stop_with(
    "results_dir 不存在: ", results_dir, "\n",
    "用法: Rscript extract_sciclone_mutation_assignment_single.R [results_dir]\n",
    "例如: Rscript extract_sciclone_mutation_assignment_single.R Shlush_analysis/results_single"
  )
}

patient_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)
patient_dirs <- patient_dirs[file.info(patient_dirs)$isdir %in% TRUE]

if (length(patient_dirs) == 0) {
  stop_with("在目录下未找到子文件夹: ", results_dir)
}

pick_latest_clusters_full <- function(patient_dir, timepoint) {
  # 例如: Patient-12_Diagnosis_dp25_clusters_full_v5.txt
  patt <- paste0("_", timepoint, "_.*_clusters_full_.*\\.txt$")
  files <- list.files(patient_dir, pattern = patt, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files[which.max(file.info(files)$mtime)]
}

find_timepoint_var_col <- function(df, timepoint) {
  # 例如: Patient-12_Diagnosis.var / Patient-12_Relapse.var
  patt <- paste0("_", timepoint, "\\.var$")
  hits <- grep(patt, colnames(df), value = TRUE)
  if (length(hits) == 0) return(NA_character_)
  if (length(hits) > 1) {
    # 理论上只应有一个；取第一个并提示
    message("[警告] 发现多个 ", timepoint, ".var 列，使用: ", hits[[1]], " | all: ", paste(hits, collapse = ", "))
  }
  hits[[1]]
}

extract_probs <- function(df) {
  prob_cols <- grep("^cluster\\.prob\\.", colnames(df), value = TRUE)
  if (length(prob_cols) == 0) return(character(0))
  prob_idx <- suppressWarnings(as.integer(sub("^cluster\\.prob\\.", "", prob_cols)))
  if (all(!is.na(prob_idx))) {
    prob_cols <- prob_cols[order(prob_idx)]
  }
  prob_cols
}

processed <- 0L
skipped <- 0L

timepoints <- c("Diagnosis", "Relapse")

for (patient_dir in patient_dirs) {
  patient_id <- basename(patient_dir)

  for (tp in timepoints) {
    in_path <- pick_latest_clusters_full(patient_dir, tp)
    if (is.na(in_path) || !file.exists(in_path)) {
      message("[跳过] 未找到 ", tp, " 的 *_clusters_full_*.txt: ", patient_id)
      skipped <- skipped + 1L
      next
    }

    dt <- read_clusters_full(in_path)

    required_cols <- c("chr", "st", "cluster")
    missing_required <- setdiff(required_cols, colnames(dt))
    if (length(missing_required) > 0) {
      message("[跳过] 缺少必要列: ", patient_id, " | ", tp, " | ", paste(missing_required, collapse = ", "), " | 文件: ", basename(in_path))
      skipped <- skipped + 1L
      next
    }

    var_col <- find_timepoint_var_col(dt, tp)
    if (is.na(var_col)) {
      message("[跳过] 未找到 ", tp, ".var 列: ", patient_id, " | 文件: ", basename(in_path))
      skipped <- skipped + 1L
      next
    }

    prob_cols <- extract_probs(dt)

    out <- data.frame(
      Chromosome = as.character(dt[["chr"]]),
      Position = dt[["st"]],
      exist = !is.na(dt[[var_col]]) & (dt[[var_col]] != 0),
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

    out_path <- file.path(results_dir, paste0(patient_id, "_", tolower(tp), "_sciClone_mutation_assignment.csv"))
    # 不给字符串加双引号（包括 Chromosome）
    utils::write.csv(out, out_path, row.names = FALSE, quote = FALSE)
    message("[完成] ", patient_id, " | ", tp, " -> ", out_path, " | 输入: ", basename(in_path))
    processed <- processed + 1L
  }
}

message("全部完成: processed=", processed, ", skipped=", skipped, ", results_dir=", results_dir)

