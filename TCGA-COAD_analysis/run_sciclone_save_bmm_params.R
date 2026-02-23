#!/usr/bin/env Rscript
################################################################################
# 单样本 SciClone 分析 + Beta 混合模型（BMM）参数保存
#
# 基于 run_sciclone_shlush_single.R，在运行单样本 SciClone 的同时，
# 将 sc_object@clust 中所有 Beta 分布相关参数保存，便于后续复现/自定义绘制
# “模型拟合线”（各成分的后验预测密度及总拟合）。
#
# 数据结构：TCGA-COAD，无 Dx/Rx 后缀，无 copy number 表（已做二倍体筛选）。
#
# 保存内容：
#   - <patient>_bmm_params.rds：完整 R 列表，含 mu, alpha, nu, beta, pi,
#     cluster.means/lower/upper, fit.x, fit.y, individual.fits.y 等，用于重绘拟合曲线。
#   - <patient>_bmm_params.csv：一个 CSV 包含该患者全部 BMM 参数（mu, alpha, nu, beta, pi），
#     行标签为 mu_1, alpha_1, ... , pi，列为 cluster_1, cluster_2, ...
################################################################################

##========================= 配置参数 =========================##

FILTER_ZERO_TUM_READ <- TRUE
MIN_DEPTH <- 25
VER <- 1

# 是否写出 BMM 参数为单个 CSV（每个 patient 一个文件）
SAVE_BMM_CSV <- TRUE

# 子采样：当突变数超过此值时，不放回抽取 n 个 mutation 再跑 SciClone，以加速大样本聚类。
# 设为 NULL 表示不子采样，使用全部突变。
SUBSAMPLE_N <- 20000L

##========================= 加载依赖 =========================##

cat("正在加载所需 R 包...\n")
suppressPackageStartupMessages({
  library(reshape2)
  library(limma)
  library(ggplot2)
})

##========================= 设置工作目录 & 加载 sciClone 源码 =========================##

project_root <- "/Users/yanjiechen/Documents/Github/sciclone"
setwd(project_root)

cat("从本仓库加载 sciClone 函数...\n")
source(file.path(project_root, "R", "object.R"))
initScClass()
source(file.path(project_root, "R", "sciClone.R"))
source(file.path(project_root, "R", "clustering.R"))
source(file.path(project_root, "R", "plots.R"))

##========================= 目录定义（TCGA-COAD_analysis）=========================##

base_dir       <- file.path(project_root, "TCGA-COAD_analysis")
mutation_dir   <- file.path(base_dir, "mutation_table")
results_dir    <- file.path(base_dir, "results_single")
bmm_params_dir <- file.path(base_dir, "bmm_params")

if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
if (!dir.exists(mutation_dir)) {
  dir.create(mutation_dir, recursive = TRUE)
  cat(sprintf("已创建目录: %s\n", mutation_dir))
}
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(bmm_params_dir)) {
  dir.create(bmm_params_dir, recursive = TRUE)
  cat(sprintf("已创建 BMM 参数目录: %s\n", bmm_params_dir))
}

##========================= 保存 BMM 参数 =========================##

# 从 scObject 的 clust 中提取并保存：RDS（完整列表）+ 一个 CSV（该 patient 全部参数）
# CSV 结构：行 = mu_1, mu_2, ..., alpha_1, ..., nu_1, ..., beta_1, ..., pi；列 = row_id, cluster_1, cluster_2, ...
save_bmm_params <- function(sc, sample_label, output_dir, save_csv = TRUE) {
  if (is.null(sc) || !inherits(sc, "scObject")) {
    warning("save_bmm_params: sc 为空或非 scObject，跳过保存。")
    return(invisible(NULL))
  }
  clust <- sc@clust
  if (is.null(clust) || is.null(clust[[1L]]) || !identical(clust$cluster.method, "bmm")) {
    cat(sprintf("  [BMM] 样本 %s 未使用 BMM 或无聚类结果，跳过 BMM 参数保存。\n", sample_label))
    return(invisible(NULL))
  }

  bmm_list <- list(
    sample_label    = sample_label,
    sampleNames     = sc@sampleNames,
    dimensions      = sc@dimensions,
    cluster.method  = clust$cluster.method,
    mu = clust$mu, alpha = clust$alpha, nu = clust$nu, beta = clust$beta, pi = clust$pi,
    cluster.means = clust$cluster.means, cluster.lower = clust$cluster.lower, cluster.upper = clust$cluster.upper,
    fit.x = clust$fit.x, fit.y = clust$fit.y, individual.fits.y = clust$individual.fits.y,
    cluster.assignments = clust$cluster.assignments, cluster.probabilities = clust$cluster.probabilities
  )

  prefix <- file.path(output_dir, paste0(sample_label, "_bmm"))
  saveRDS(bmm_list, paste0(prefix, "_params.rds"))
  cat(sprintf("  [BMM] 已保存: %s_params.rds\n", prefix))

  if (save_csv) {
    K <- ncol(clust$mu)
    n_dim <- nrow(clust$mu)
    col_names <- c("row_id", paste0("cluster_", seq_len(K)))
    rows <- list()
    for (i in seq_len(n_dim)) {
      rows[[length(rows) + 1L]] <- c(paste0("mu_", i), clust$mu[i, ])
      rows[[length(rows) + 1L]] <- c(paste0("alpha_", i), clust$alpha[i, ])
      rows[[length(rows) + 1L]] <- c(paste0("nu_", i), clust$nu[i, ])
      rows[[length(rows) + 1L]] <- c(paste0("beta_", i), clust$beta[i, ])
    }
    rows[[length(rows) + 1L]] <- c("pi", clust$pi)
    tab <- do.call(rbind, lapply(rows, function(r) c(r[1], as.character(r[-1]))))
    df <- as.data.frame(tab, stringsAsFactors = FALSE)
    names(df) <- col_names
    csv_file <- paste0(prefix, "_params.csv")
    write.csv(df, csv_file, row.names = FALSE)
    cat(sprintf("  [BMM] 已保存: %s\n", csv_file))
  }

  invisible(bmm_list)
}

##========================= 单样本聚类函数（含 BMM 保存）=========================##

singleSampleClonalAnalysis <- function(mutation_df,
                                       sample_label,
                                       MIN_DEPTH = 25,
                                       VER = 1,
                                       output_dir = ".",
                                       bmm_output_dir = NULL,
                                       SAVE_BMM_CSV = TRUE) {
  factorName <- sample_label

  # TCGA-COAD：mutation 已做二倍体筛选，不传入 copy number
  sc <- sciClone(
    vafs = mutation_df[, 1:5],
    copyNumberCalls = NULL,
    regionsToExclude = NULL,
    annotation = if ("gene_name" %in% colnames(mutation_df) && any(mutation_df$gene_name != "")) {
      mutation_df[mutation_df$gene_name != "", c(1, 2, 6), drop = FALSE]
    } else {
      NULL
    },
    sampleNames = factorName,
    useSexChrs = TRUE,
    minimumDepth = MIN_DEPTH,
    verbose = 0,
    doClusteringAlongMargins = FALSE,
    plotIntermediateResults = 0,
    cnCallsAreLog2 = FALSE,
    maximumClusters = 10
  )

  # 写出完整突变级结果表
  writeClusterTable(
    sc,
    file.path(output_dir, paste0(sample_label, "_dp", MIN_DEPTH, "_clusters_full_v", VER, ".txt"))
  )

  DTA <- sc@vafs.merged
  if (!"cluster" %in% colnames(DTA) || all(is.na(DTA$cluster))) {
    warning("No clusters identified for sample: ", sample_label)
  } else {
    Nclust <- max(DTA$cluster, na.rm = TRUE)
    ClustStats <- array(0, dim = c(Nclust, 2))
    for (k in seq_len(Nclust)) {
      idx <- DTA$cluster == k & !is.na(DTA$cluster)
      ClustStats[k, ] <- c(
        mean(DTA[idx, paste0(factorName, ".vaf")], na.rm = TRUE),
        Nmut = sum(idx)
      )
    }
    colnames(ClustStats) <- c(factorName, "Nmut")
    ClustStatsDF <- data.frame(cluster = seq_len(Nclust), round(ClustStats, 3))
    write.table(
      ClustStatsDF,
      file = file.path(output_dir, paste0(sample_label, "_dp", MIN_DEPTH, "_clusters_v", VER, ".txt")),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
    writeClusterSummaryTable(
      sc,
      file.path(output_dir, paste0(sample_label, "_dp", MIN_DEPTH, "_cluster_summary_v", VER))
    )
  }

  # 生成 1D 克隆图
  tryCatch({
    sc.plot1d(sc, outputFile = file.path(output_dir, paste0(sample_label, "_plot_1d.pdf")))
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to generate 1D plot for %s: %s\n", sample_label, e$message))
  })

  # 保存 BMM 参数（便于复现拟合线）
  save_bmm_params(
    sc,
    sample_label,
    output_dir = if (!is.null(bmm_output_dir)) bmm_output_dir else output_dir,
    save_csv = SAVE_BMM_CSV
  )

  return(sc)
}

##========================= 样本发现（无 Dx/Rx 后缀）=========================##

find_all_samples <- function(mutation_dir) {
  mut_files <- list.files(
    mutation_dir,
    pattern = "_mutation_table\\.csv$",
    full.names = FALSE
  )
  # 只保留 <Patient>_mutation_table.csv，排除已带 Dx/Rx 的（若有）
  mut_files <- mut_files[!grepl("_(Dx|Rx)_mutation_table\\.csv$", mut_files)]
  sample_ids <- gsub("_mutation_table\\.csv$", "", mut_files)
  sort(unique(sample_ids))
}

##========================= 单患者处理（无 Dx/Rx，无 CN）=========================##

process_patient_single <- function(sample_id,
                                   mutation_dir,
                                   results_dir,
                                   bmm_params_dir,
                                   FILTER_ZERO_TUM_READ = TRUE,
                                   MIN_DEPTH = 25,
                                   VER = 1,
                                   SAVE_BMM_CSV = TRUE,
                                   SUBSAMPLE_N = NULL) {
  cat("\n========================================\n")
  cat(sprintf("Processing patient: %s\n", sample_id))
  cat("========================================\n\n")

  mut_file <- file.path(mutation_dir, paste0(sample_id, "_mutation_table.csv"))
  if (!file.exists(mut_file)) {
    cat("ERROR: No mutation table for patient ", sample_id, ".\n")
    return(list(processed = TRUE, status = "failed", error_message = "No mutation table"))
  }

  output_dir <- file.path(results_dir, sample_id)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(bmm_params_dir)) dir.create(bmm_params_dir, recursive = TRUE)

  read_mutation_table <- function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if ("tum_read" %in% colnames(df) && FILTER_ZERO_TUM_READ) {
      df <- df[df$tum_read != 0, , drop = FALSE]
    }
    if ("contig" %in% colnames(df)) {
      df$contig <- gsub("^chr", "", df$contig)
    }
    df
  }

  res <- tryCatch({
    mut_df <- read_mutation_table(mut_file)
    n_orig <- nrow(mut_df)
    if (!is.null(SUBSAMPLE_N) && n_orig > SUBSAMPLE_N) {
      set.seed(42L)
      idx <- sample.int(n_orig, size = SUBSAMPLE_N, replace = FALSE)
      mut_df <- mut_df[idx, , drop = FALSE]
      cat(sprintf("  [子采样] 突变数 %d -> %d（SUBSAMPLE_N = %d）\n", n_orig, nrow(mut_df), SUBSAMPLE_N))
    }
    sc <- singleSampleClonalAnalysis(
      mutation_df = mut_df,
      sample_label = sample_id,
      MIN_DEPTH = MIN_DEPTH,
      VER = VER,
      output_dir = output_dir,
      bmm_output_dir = bmm_params_dir,
      SAVE_BMM_CSV = SAVE_BMM_CSV
    )
    list(ok = !is.null(sc) && inherits(sc, "scObject"), msg = "")
  }, error = function(e) {
    list(ok = FALSE, msg = e$message)
  })

  list(
    processed = TRUE,
    status = if (res$ok) "success" else "failed",
    error_message = res$msg
  )
}

##========================= 主程序 =========================##

cat("\n========================================\n")
cat("SciClone 单样本分析 + BMM 参数保存 (TCGA-COAD_analysis)\n")
cat("========================================\n\n")

cat(sprintf("突变表目录: %s\n", mutation_dir))
sample_ids <- find_all_samples(mutation_dir)
if (length(sample_ids) == 0) {
  cat("未发现 *_mutation_table.csv，请将突变表放入 mutation_table 目录。\n")
  quit(status = 0)
}
cat(sprintf("共 %d 个患者: %s\n", length(sample_ids), paste(sample_ids, collapse = ", ")))
cat(sprintf("BMM 参数将保存到: %s（每患者一个 RDS + 一个 CSV）\n", bmm_params_dir))
if (!is.null(SUBSAMPLE_N)) {
  cat(sprintf("子采样: 突变数 > %d 时将不放回抽取 %d 个再跑 SciClone\n", SUBSAMPLE_N, SUBSAMPLE_N))
} else {
  cat("子采样: 关闭（使用全部突变）\n")
}

status_file <- file.path(base_dir, "processing_status_single.csv")
if (file.exists(status_file)) {
  status_df <- read.csv(status_file, stringsAsFactors = FALSE)
  processed <- unique(status_df$Sample_ID[
    status_df$processed == TRUE | status_df$processed == "TRUE" | status_df$processed == "true"
  ])
  if (length(processed) > 0) {
    sample_ids <- setdiff(sample_ids, processed)
    cat(sprintf("已处理（跳过）: %s；剩余: %d\n", paste(processed, collapse = ", "), length(sample_ids)))
  }
} else {
  status_df <- data.frame(
    Sample_ID = character(0), processed = logical(0), status = character(0),
    error_message = character(0), timestamp = character(0), stringsAsFactors = FALSE
  )
  write.csv(status_df, status_file, row.names = FALSE)
}

if (length(sample_ids) == 0) {
  cat("无待处理样本。\n")
  quit(status = 0)
}

for (sample_id in sample_ids) {
  result <- process_patient_single(
    sample_id,
    mutation_dir,
    results_dir,
    bmm_params_dir,
    FILTER_ZERO_TUM_READ = FILTER_ZERO_TUM_READ,
    MIN_DEPTH = MIN_DEPTH,
    VER = VER,
    SAVE_BMM_CSV = SAVE_BMM_CSV,
    SUBSAMPLE_N = SUBSAMPLE_N
  )
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  idx <- which(status_df$Sample_ID == sample_id)
  if (length(idx) > 0) {
    status_df$processed[idx] <- result$processed
    status_df$status[idx] <- result$status
    status_df$error_message[idx] <- result$error_message
    status_df$timestamp[idx] <- timestamp
  } else {
    status_df <- rbind(status_df, data.frame(
      Sample_ID = sample_id, processed = result$processed, status = result$status,
      error_message = result$error_message, timestamp = timestamp, stringsAsFactors = FALSE
    ))
  }
  write.csv(status_df, status_file, row.names = FALSE)
}

cat("\n=== 完成。BMM 参数已保存至 bmm_params 目录（每患者一个 .rds + 一个 .csv）。===\n")
cat("重绘拟合线示例：\n")
cat("  bmm <- readRDS(\"TCGA-COAD_analysis/bmm_params/<patient>_bmm_params.rds\")\n")
cat("  plot(bmm$fit.x, bmm$fit.y[1,], type=\"l\", xlab=\"VAF\", ylab=\"Density\")\n")
cat("  for (k in seq_along(bmm$individual.fits.y)) lines(bmm$fit.x, bmm$individual.fits.y[[k]][1,], lty=2)\n")
