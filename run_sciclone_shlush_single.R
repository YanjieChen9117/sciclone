#!/usr/bin/env Rscript
################################################################################
# Shlush 队列 SciClone 单样本克隆分析脚本
#
# 主要功能：
#   - 批量遍历 Shlush_analysis 目录中的所有患者样本
#   - 对每个 patient，分别对 Diagnosis 与 Relapse 两个样本独立运行一次 sciClone
#   - 输入使用已经预处理好的单样本 mutation table（CSV）及对应时间点的 copy number 表（若存在则传入，否则假定 CN=2）
#   - 写出 sciClone 的表格结果与 1D/2D 图（单样本不画 2D 对比图和 margins 图）
#   - 在 Shlush_analysis/processing_status_single.csv 中记录处理状态，可断点续跑
#   - 将 BMM（Beta 混合模型）参数保存到 results_single/<sample_id>/，供后续复现/绘制拟合线
#   - 支持手动跳过列表 SKIP_MANUAL：指定 (Sample_ID, Dx/Rx) 不运行 sciClone
#
# 与 run_sciclone_shlush.R 的区别：
#   - 这里是“单样本模式”：一次 sciClone 只处理一个样本（一个时间点）
#   - 同一个 patient 运行两次 sciClone：<patient>_Diagnosis 与 <patient>_Relapse
#   - 读入 mutation table 后必须先过滤掉 tum_read == 0 的行（仅本脚本需要，由 FILTER_ZERO_TUM_READ 控制，默认 TRUE）
################################################################################

##========================= 配置参数 =========================##

# 是否在读入 mutation table 后过滤掉 tum_read == 0 的行（单样本模式建议保持 TRUE）
FILTER_ZERO_TUM_READ <- TRUE

# 最小测序深度阈值（传给 sciClone 的 minimumDepth）
MIN_DEPTH <- 25

# 版本号，用于输出文件名后缀
VER <- 5

# 是否写出 BMM 参数为 CSV（除 .rds 外再写一份表格，便于查看）
SAVE_BMM_CSV <- TRUE

# 手动跳过列表：列出要跳过的 (Sample_ID, Dx/Rx)
# 每项为 c("Sample_ID", "Dx") 或 c("Sample_ID", "Rx")，例如 c("Patient-1", "Dx") 表示跳过该患者的 Diagnosis
SKIP_MANUAL <- list(
  c("Patient-4", "Rx")
  # c("Patient-2", "Rx")
)

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

##========================= 目录定义 =========================##

mutation_dir <- file.path(project_root, "Shlush_analysis", "mutation_table")
cn_dir      <- file.path(project_root, "Shlush_analysis", "copy_number_table")
results_dir <- file.path(project_root, "Shlush_analysis", "results_single")

if (!dir.exists(mutation_dir)) {
  stop(sprintf("Mutation table directory does not exist: %s", mutation_dir))
}
if (!dir.exists(cn_dir)) {
  stop(sprintf("Copy number table directory does not exist: %s", cn_dir))
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 判断 (sample_id, Dx/Rx) 是否在手动跳过列表中
is_skipped_manual <- function(sample_id, dx_or_rx) {
  if (length(SKIP_MANUAL) == 0) return(FALSE)
  any(vapply(SKIP_MANUAL, function(x) {
    length(x) == 2L && identical(x[1], sample_id) && identical(x[2], dx_or_rx)
  }, logical(1)))
}

##========================= 保存 BMM 参数 =========================##

# 从 scObject 的 clust 中提取并保存：RDS（完整列表）+ 可选 CSV
# 结果写入 output_dir（即 results_single/<sample_id>/），便于后续复现/绘制拟合线
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

##========================= 单样本聚类函数 =========================##

# 对单个样本（单一时间点）运行 sciClone，并写出结果
# 若提供 cnv_df（列格式：contig, start, stop, cn, sample, type），则传入 copy number 与 LOH 排除区域；否则假定 CN=2
singleSampleClonalAnalysis <- function(mutation_df,
                                       sample_label,   # 例如 "Patient-1_Diagnosis"
                                       cnv_df = NULL,  # 可选，该时间点 CN 表；NULL 时假定所有位点 CN=2
                                       MIN_DEPTH = 25,
                                       VER = 1,
                                       output_dir = ".",
                                       save_bmm_csv = TRUE) {
  factorName <- sample_label

  # 单样本时 vafs 为 data.frame，sciClone 内部会执行 copyNumberCalls=list(copyNumberCalls)，
  # 故这里必须传 data.frame 而非 list(data.frame)，否则会变成 list(list(...)) 导致 "incorrect number of dimensions"
  if (!is.null(cnv_df) && nrow(cnv_df) > 0) {
    loh_regions <- cnv_df[cnv_df$type == "LOH", 1:3, drop = FALSE]
    copyNumberCalls <- cnv_df[, 1:4]
    regionsToExclude <- list(loh_regions)
  } else {
    copyNumberCalls <- NULL
    regionsToExclude <- NULL
  }

  sc <- sciClone(
    vafs = mutation_df[, 1:5],
    copyNumberCalls = copyNumberCalls,
    regionsToExclude = regionsToExclude,
    annotation = mutation_df[mutation_df$gene_name != "", c(1, 2, 6), drop = FALSE],
    sampleNames = factorName,
    useSexChrs = TRUE,
    minimumDepth = MIN_DEPTH,
    verbose = 0,
    doClusteringAlongMargins = FALSE,   # 单样本不做边缘 1D clustering
    plotIntermediateResults = 0,
    cnCallsAreLog2 = FALSE,
    maximumClusters = 10
  )

  # 写出完整突变级结果表
  writeClusterTable(
    sc,
    file.path(output_dir, paste0(sample_label, "_dp", MIN_DEPTH, "_clusters_full_v", VER, ".txt"))
  )

  # 写出簇级汇总表（本质与双样本版本类似，但只有 1 列平均 VAF）
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

    # 写出 cluster summary（三个文件：means/lower/upper）
    writeClusterSummaryTable(
      sc,
      file.path(output_dir, paste0(sample_label, "_dp", MIN_DEPTH, "_cluster_summary_v", VER))
    )
  }

  # 生成 1D 克隆图（单样本最有用的图）
  tryCatch({
    plot_file_1d <- file.path(output_dir, paste0(sample_label, "_plot_1d.pdf"))
    sc.plot1d(sc, outputFile = plot_file_1d)
  }, error = function(e) {
    cat(sprintf("WARNING: Failed to generate 1D plot for %s: %s\n", sample_label, e$message))
  })

  # 对单样本而言，标准 sc.plot2d / sc.plot2dWithMargins 需要至少两维，这里不调用

  # 保存 BMM 参数到 output_dir（results_single/<sample_id>/），便于后续复现/绘制拟合线
  save_bmm_params(sc, sample_label, output_dir, save_csv = save_bmm_csv)

  return(sc)
}

##========================= 样本发现函数 =========================##

# 扫描 mutation_table 目录，推断所有 patient ID
# 与双样本脚本一致：基于 *_Dx_mutation_table.csv（新命名）
find_all_samples <- function(mutation_dir) {
  mut_files <- list.files(
    mutation_dir,
    pattern = "_Dx_mutation_table\\.csv$",
    full.names = FALSE
  )
  sample_ids <- gsub("_Dx_mutation_table\\.csv$", "", mut_files)
  return(sort(sample_ids))
}

##========================= 单个 patient 的处理函数 =========================##

# 对单个 patient，同步完成两次单样本 sciClone：
#   - <id>_Diagnosis
#   - <id>_Relapse
# 若存在对应时间点的 CN 表则读入并传入 sciClone，否则假定 CN=2
process_patient_single <- function(sample_id,
                                   mutation_dir,
                                   cn_dir,
                                   results_dir,
                                   FILTER_ZERO_TUM_READ = TRUE) {
  cat("\n")
  cat("========================================\n")
  cat(sprintf("Processing patient (single-sample mode): %s\n", sample_id))
  cat("========================================\n\n")

  # 文件路径（mutation table + copy number table，新命名：Dx/Rx）
  mut_diagnosis_file <- file.path(
    mutation_dir,
    paste0(sample_id, "_Dx_mutation_table.csv")
  )
  mut_relapse_file <- file.path(
    mutation_dir,
    paste0(sample_id, "_Rx_mutation_table.csv")
  )
  cn_dx_file <- file.path(cn_dir, paste0(sample_id, "_Dx_copy_number_table.csv"))
  cn_rx_file <- file.path(cn_dir, paste0(sample_id, "_Rx_copy_number_table.csv"))

  # 至少需要某个时间点的 mutation 表存在即可
  has_dx <- file.exists(mut_diagnosis_file)
  has_rx <- file.exists(mut_relapse_file)

  if (!has_dx && !has_rx) {
    error_msg <- sprintf(
      "No single-sample mutation table (Dx or Rx) for patient %s.",
      sample_id
    )
    cat("ERROR: ", error_msg, "\n")
    return(list(
      processed = TRUE,
      status = "failed",
      error_message = error_msg
    ))
  }

  # 为该患者创建结果目录
  output_dir <- file.path(results_dir, sample_id)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # 内部小工具：读入 + 可选过滤 + 规范化 contig
  read_mutation_table <- function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if (!("tum_read" %in% colnames(df))) {
      warning("Column 'tum_read' not found in file: ", path, "; FILTER_ZERO_TUM_READ will be ignored for this file.")
    } else if (FILTER_ZERO_TUM_READ) {
      before_n <- nrow(df)
      df <- df[df$tum_read != 0, , drop = FALSE]
      after_n <- nrow(df)
      cat(sprintf("  Filtering tum_read == 0 in %s: %d -> %d rows\n",
                  basename(path), before_n, after_n))
    }
    # 去掉 "chr" 前缀
    if ("contig" %in% colnames(df)) {
      df$contig <- gsub("^chr", "", df$contig)
    }
    df
  }

  # 读入 CN 表并规范化 contig（与双样本脚本一致）
  read_cn_table <- function(path) {
    if (!file.exists(path)) return(NULL)
    df <- read.csv(path, stringsAsFactors = FALSE)
    if ("contig" %in% colnames(df)) {
      df$contig <- gsub("^chr", "", df$contig)
    }
    df
  }

  # 记录每个时间点的运行结果（方便组合成患者级状态）
  dx_status <- NULL
  rx_status <- NULL

  ##----------------- Diagnosis 单样本 -----------------##
  if (has_dx) {
    if (is_skipped_manual(sample_id, "Dx")) {
      cat("\n--- [手动跳过] Diagnosis ---\n")
      dx_status <- list(ok = TRUE, msg = "skipped (manual)")
    } else {
      cat("\n--- Running single-sample sciClone for Diagnosis ---\n")
      cat(sprintf("  Mutation file : %s\n", mut_diagnosis_file))
      Tumor_cnv <- read_cn_table(cn_dx_file)
      if (!is.null(Tumor_cnv)) {
        cat(sprintf("  CN file (Dx)  : %s (%d segments)\n", cn_dx_file, nrow(Tumor_cnv)))
      } else {
        cat("  CN file (Dx)  : not found, assuming CN=2\n")
      }

      dx_res <- tryCatch({
        Tumor_mut <- read_mutation_table(mut_diagnosis_file)

        cat(sprintf("  Diagnosis mutations: %d\n", nrow(Tumor_mut)))

        sc_dx <- singleSampleClonalAnalysis(
          mutation_df = Tumor_mut,
          sample_label = paste0(sample_id, "_Diagnosis"),
          cnv_df = Tumor_cnv,
          MIN_DEPTH = MIN_DEPTH,
          VER = VER,
          output_dir = output_dir,
          save_bmm_csv = SAVE_BMM_CSV
        )

        list(sc = sc_dx, error = NULL)
      }, error = function(e) {
      msg <- sprintf("Diagnosis single-sample analysis failed: %s", e$message)
      cat(msg, "\n")
      list(sc = NULL, error = msg)
    })

      if (is.null(dx_res$sc)) {
        dx_status <- list(ok = FALSE, msg = dx_res$error)
      } else if (!inherits(dx_res$sc, "scObject")) {
        msg <- sprintf("Diagnosis analysis returned invalid object type: %s",
                       class(dx_res$sc)[1])
        dx_status <- list(ok = FALSE, msg = msg)
      } else {
        dx_status <- list(ok = TRUE, msg = "")
      }
    }
  }

  ##----------------- Relapse 单样本 -----------------##
  if (has_rx) {
    if (is_skipped_manual(sample_id, "Rx")) {
      cat("\n--- [手动跳过] Relapse ---\n")
      rx_status <- list(ok = TRUE, msg = "skipped (manual)")
    } else {
      cat("\n--- Running single-sample sciClone for Relapse ---\n")
      cat(sprintf("  Mutation file : %s\n", mut_relapse_file))
      Relapse_cnv <- read_cn_table(cn_rx_file)
      if (!is.null(Relapse_cnv)) {
        cat(sprintf("  CN file (Rx)  : %s (%d segments)\n", cn_rx_file, nrow(Relapse_cnv)))
      } else {
        cat("  CN file (Rx)  : not found, assuming CN=2\n")
      }

      rx_res <- tryCatch({
        Relapse_mut <- read_mutation_table(mut_relapse_file)

        cat(sprintf("  Relapse mutations: %d\n", nrow(Relapse_mut)))

        sc_rx <- singleSampleClonalAnalysis(
          mutation_df = Relapse_mut,
          sample_label = paste0(sample_id, "_Relapse"),
          cnv_df = Relapse_cnv,
          MIN_DEPTH = MIN_DEPTH,
          VER = VER,
          output_dir = output_dir,
          save_bmm_csv = SAVE_BMM_CSV
        )

        list(sc = sc_rx, error = NULL)
      }, error = function(e) {
      msg <- sprintf("Relapse single-sample analysis failed: %s", e$message)
      cat(msg, "\n")
      list(sc = NULL, error = msg)
    })

      if (is.null(rx_res$sc)) {
        rx_status <- list(ok = FALSE, msg = rx_res$error)
      } else if (!inherits(rx_res$sc, "scObject")) {
        msg <- sprintf("Relapse analysis returned invalid object type: %s",
                       class(rx_res$sc)[1])
        rx_status <- list(ok = FALSE, msg = msg)
      } else {
        rx_status <- list(ok = TRUE, msg = "")
      }
    }
  }

  ##----------------- 汇总患者级状态 -----------------##
  msgs <- c()
  ok_vec <- c()
  if (!is.null(dx_status)) {
    ok_vec <- c(ok_vec, dx_status$ok)
    if (!dx_status$ok) msgs <- c(msgs, paste0("Diagnosis: ", dx_status$msg))
  }
  if (!is.null(rx_status)) {
    ok_vec <- c(ok_vec, rx_status$ok)
    if (!rx_status$ok) msgs <- c(msgs, paste0("Relapse: ", rx_status$msg))
  }

  all_ok <- all(ok_vec)
  if (all_ok) {
    cat(sprintf("\nPatient %s single-sample analyses completed successfully.\n", sample_id))
    return(list(processed = TRUE, status = "success", error_message = ""))
  } else {
    err_msg <- paste(msgs, collapse = " | ")
    cat(sprintf("\nPatient %s single-sample analyses finished with errors: %s\n",
                sample_id, err_msg))
    return(list(processed = TRUE, status = "failed", error_message = err_msg))
  }
}

##========================= 主程序：批量调度 =========================##

cat("\n========================================\n")
cat("SciClone Analysis for Shlush (单样本批量分析)\n")
cat("========================================\n\n")

cat(sprintf("扫描突变表目录以获取样本列表：%s\n", mutation_dir))
sample_ids <- find_all_samples(mutation_dir)
cat(sprintf("\n共发现 %d 个患者待处理（基于 Diagnosis 突变表）：%s\n",
            length(sample_ids), paste(sample_ids, collapse = ", ")))

if (length(SKIP_MANUAL) > 0) {
  cat(sprintf("手动跳过列表（共 %d 项）：\n", length(SKIP_MANUAL)))
  for (x in SKIP_MANUAL) {
    if (length(x) >= 2) cat(sprintf("  - %s / %s\n", x[1], x[2]))
  }
}

# 单样本模式的状态文件
status_file <- file.path(project_root, "Shlush_analysis", "processing_status_single.csv")

if (file.exists(status_file)) {
  cat(sprintf("\nReading single-sample processing status file: %s\n", status_file))
  status_df <- read.csv(status_file, stringsAsFactors = FALSE)

  processed_samples <- unique(status_df$Sample_ID[
    status_df$processed == TRUE |
      status_df$processed == "TRUE" |
      status_df$processed == "true"
  ])
  if (length(processed_samples) > 0) {
    cat(sprintf("检测到 %d 个患者（单样本模式）已处理完成，跳过：%s\n",
                length(processed_samples), paste(processed_samples, collapse = ", ")))
    sample_ids <- setdiff(sample_ids, processed_samples)
    cat(sprintf("剩余待处理患者数量：%d\n", length(sample_ids)))
  }
} else {
  cat("未发现 processing_status_single.csv，创建新的状态跟踪文件\n")
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
  cat("\n所有患者（单样本模式）均已处理完成，无需再次运行。\n")
  quit(status = 0)
}

success_count <- 0
fail_count <- 0
failed_samples <- c()

for (sample_id in sample_ids) {
  result <- process_patient_single(
    sample_id,
    mutation_dir,
    cn_dir,
    results_dir,
    FILTER_ZERO_TUM_READ = FILTER_ZERO_TUM_READ
  )

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  existing_idx <- which(status_df$Sample_ID == sample_id)

  if (length(existing_idx) > 0) {
    status_df$processed[existing_idx] <- result$processed
    status_df$status[existing_idx] <- result$status
    status_df$error_message[existing_idx] <- result$error_message
    status_df$timestamp[existing_idx] <- timestamp
  } else {
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

  write.csv(status_df, status_file, row.names = FALSE)

  if (result$status == "success") {
    success_count <- success_count + 1
    cat(sprintf("\n状态文件已更新：%s (single) 标记为 success\n", sample_id))
  } else {
    fail_count <- fail_count + 1
    failed_samples <- c(failed_samples, sample_id)
    cat(sprintf("\n状态文件已更新：%s (single) 标记为 failed，原因：%s\n",
                sample_id, result$error_message))
  }
}

cat("\n")
cat("========================================\n")
cat("单样本模式运行总结\n")
cat("========================================\n")
cat(sprintf("本次运行共处理患者数: %d\n", length(sample_ids)))
cat(sprintf("成功患者数          : %d\n", success_count))
cat(sprintf("失败患者数          : %d\n", fail_count))

if (length(failed_samples) > 0) {
  cat(sprintf("\n失败患者列表: %s\n", paste(failed_samples, collapse = ", ")))
}

cat("\n=== 单样本脚本运行结束 ===\n")

