#!/usr/bin/env Rscript
################################################################################
# Shlush：从 BMM 参数 CSV 绘制“拟合线”堆积图，并叠加样本原始 VAF 直方图
#
# 针对 Shlush_analysis/results_single 目录结构：
#   results_single/<Patient-ID>/<Patient-ID>_Diagnosis_bmm_params.csv
#   results_single/<Patient-ID>/<Patient-ID>_Relapse_bmm_params.csv
# 突变表：mutation_table/<Patient-ID>_Dx_mutation_table.csv / _Rx_mutation_table.csv
#
# 用法：
#   source("Shlush_analysis/plot_bmm_fit.R")
#   plot_bmm_fit_from_csv("Shlush_analysis/results_single/Patient-1/Patient-1_Diagnosis_bmm_params.csv",
#                         output_file = "Shlush_analysis/plots/Patient-1_Diagnosis_bmm_fit.png")
################################################################################

# 与 run_sciclone_shlush_single.R 中 VER 一致，用于定位 cluster_summary .means 文件
SHLUSH_CLUSTER_SUMMARY_VER <- 5L
SHLUSH_MIN_DEPTH <- 25L

# 计算密度需用 bmm 包（与 sciClone 一致）
if (!requireNamespace("bmm", quietly = TRUE)) {
  stop("请先安装 bmm: install.packages('bmm') 或 devtools::install_github('genome/bmm')")
}
suppressPackageStartupMessages(library(bmm))

##========================= 原始数据直方图（mutation table）=========================##

#' 从样本 mutation table 计算 VAF 密度直方图（1% bin）
#' 过滤与 run_sciclone_shlush_single.R 一致：norm_read + tum_read >= 25，且 tum_read != 0（去掉 VAF=0）
compute_data_histogram <- function(mutation_table_path) {
  if (!file.exists(mutation_table_path)) return(NULL)
  d <- read.csv(mutation_table_path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("norm_read", "tum_read", "tum_vaf")
  if (!all(need %in% names(d))) return(NULL)
  d <- d[(d$norm_read + d$tum_read) >= 25, , drop = FALSE]
  if ("tum_read" %in% names(d)) d <- d[d$tum_read != 0, , drop = FALSE]
  n <- nrow(d)
  if (n == 0) return(list(breaks = 0:100, density = rep(0, 100), n = 0))
  vaf_pct <- as.numeric(d$tum_vaf)
  vaf_pct <- vaf_pct[!is.na(vaf_pct) & vaf_pct > 0 & vaf_pct <= 100]
  n <- length(vaf_pct)
  if (n == 0) return(list(breaks = 0:100, density = rep(0, 100), n = 0))
  bin_idx <- pmin(floor(vaf_pct), 99) + 1L
  counts <- tabulate(bin_idx, nbins = 100)
  density <- counts / n
  list(breaks = 0:100, density = density, n = n)
}

##========================= 解析 BMM 参数 CSV =========================##

#' 从 BMM 参数 CSV 读入并解析为矩阵/向量
read_bmm_params_csv <- function(csv_path) {
  d <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"row_id" %in% names(d)) {
    if (ncol(d) >= 2 && names(d)[1] == "row_id") { } else {
      stop("CSV 需包含列 row_id 及 cluster_1, cluster_2, ...")
    }
  }
  rid <- d[, 1L]
  cluster_cols <- setdiff(names(d), "row_id")
  cluster_cols <- cluster_cols[grepl("^cluster_", cluster_cols)]
  K <- length(cluster_cols)
  if (K == 0) stop("未找到 cluster_* 列")

  n_dim <- sum(grepl("^mu_[0-9]+$", rid))
  if (n_dim == 0) n_dim <- 1L

  mu    <- matrix(as.numeric(NA), nrow = n_dim, ncol = K)
  alpha <- matrix(as.numeric(NA), nrow = n_dim, ncol = K)
  nu    <- matrix(as.numeric(NA), nrow = n_dim, ncol = K)
  beta  <- matrix(as.numeric(NA), nrow = n_dim, ncol = K)

  for (i in seq_len(n_dim)) {
    r_mu    <- which(rid == paste0("mu_", i))
    r_alpha <- which(rid == paste0("alpha_", i))
    r_nu    <- which(rid == paste0("nu_", i))
    r_beta  <- which(rid == paste0("beta_", i))
    if (length(r_mu) == 1)    mu[i, ]    <- as.numeric(d[r_mu,    cluster_cols])
    if (length(r_alpha) == 1) alpha[i, ] <- as.numeric(d[r_alpha, cluster_cols])
    if (length(r_nu) == 1)    nu[i, ]    <- as.numeric(d[r_nu,    cluster_cols])
    if (length(r_beta) == 1)  beta[i, ]  <- as.numeric(d[r_beta,  cluster_cols])
  }
  r_pi <- which(rid == "pi")
  if (length(r_pi) != 1) stop("CSV 需包含一行 row_id == 'pi'")
  pi <- as.numeric(d[r_pi, cluster_cols])

  list(mu = mu, alpha = alpha, nu = nu, beta = beta, pi = pi,
       n_dim = n_dim, K = K, cluster_cols = cluster_cols)
}

##========================= 由参数计算密度曲线（调用 bmm）=========================##

#' 给定 BMM 参数，在 x 网格上计算各 cluster 的后验预测密度
compute_bmm_density_curves <- function(params, n_grid = 999, num_samples_bmm = 100) {
  x <- seq(0, 1, length.out = n_grid + 2)[-c(1, n_grid + 2)]
  n <- length(x)
  K <- params$K
  n_dim <- params$n_dim
  mu    <- params$mu
  alpha <- params$alpha
  nu    <- params$nu
  beta  <- params$beta
  pi    <- params$pi

  y_per_cluster <- vector("list", K)
  for (k in seq_len(K)) {
    y_per_cluster[[k]] <- matrix(0, nrow = n_dim, ncol = n)
    for (dim in seq_len(n_dim)) {
      for (i in seq_len(n)) {
        y_per_cluster[[k]][dim, i] <- bmm.component.posterior.predictive.density(
          x[i], mu[dim, k], alpha[dim, k], nu[dim, k], beta[dim, k], pi[k],
          num.samples = num_samples_bmm
        )
      }
    }
  }
  y_total <- matrix(0, nrow = n_dim, ncol = n)
  for (k in seq_len(K)) y_total <- y_total + y_per_cluster[[k]]

  list(
    x_vaf = x * 100,
    y_per_cluster = y_per_cluster,
    y_total = y_total,
    n_dim = n_dim,
    K = K
  )
}

##========================= 堆积图绘制 =========================##

cluster_colors <- function(K, alpha = 0.6) {
  palette <- c(
    "Cluster 1" = "#D55E00",
    "Cluster 2" = "#0072B2",
    "Cluster 3" = "#009E73",
    "Cluster 4" = "#CC79A7",
    "Cluster 5" = "#E69F00",
    "Cluster 6" = "#56B4E9"
  )
  base <- palette[seq_len(min(K, length(palette)))]
  if (K > length(palette)) {
    extra <- rep(palette, ceiling(K / length(palette)))[(length(palette) + 1):K]
    base <- c(base, extra)
  }
  cols <- adjustcolor(base, alpha.f = alpha)
  names(cols) <- NULL
  cols
}

#' 绘制单维度 BMM 拟合线堆积图，可选叠加原始数据直方图
plot_bmm_stacked_one_dim <- function(curves, dim = 1L, output_file = NULL,
                                    xlim = c(0, 100),
                                    sample_id = NULL, n_mutations = NULL, truncal_mean_vaf = NULL,
                                    hist_data = NULL,
                                    xlab = "Variant Allele Frequency", ylab = "Density",
                                    legend = FALSE) {
  x <- curves$x_vaf
  K <- curves$K
  y_list <- lapply(curves$y_per_cluster, function(m) m[dim, ])
  y_total <- curves$y_total[dim, ]
  y_total <- y_total / 100
  y_list <- lapply(y_list, function(y) y / 100)

  y_max_data <- 0
  if (!is.null(hist_data) && hist_data$n > 0) {
    y_max_data <- 100 * max(hist_data$density, na.rm = TRUE)
  }
  y_max_bmm <- max(y_total, na.rm = TRUE)
  y_max <- max(y_max_data, y_max_bmm, 1e-6)
  ylim <- c(0, y_max * 1.05)

  cols <- cluster_colors(K)
  if (!is.null(output_file)) png(output_file, width = 6, height = 3, units = "in", res = 300, bg = "white")
  par(mar = c(2.5, 2.5, 1.5, 1.2), cex.lab = 0.85, cex.axis = 0.85)
  plot(NA, type = "n", xlim = xlim, ylim = ylim,
       xaxs = "i", yaxs = "i",
       xaxt = "n", yaxt = "n",
       bty = "n",
       main = "", xlab = xlab, ylab = ylab)
  if (!is.null(hist_data) && hist_data$n > 0) {
    br <- hist_data$breaks
    d <- hist_data$density
    for (i in seq_len(length(br) - 1)) {
      rect(br[i], 0, br[i + 1], 100 * d[i], col = "#999999", border = NA)
    }
  }
  y_bottom <- rep(0, length(x))
  for (k in seq_len(K)) {
    y_top <- y_bottom + y_list[[k]]
    polygon(c(x, rev(x)), c(y_bottom, rev(y_top)), col = cols[k], border = NA)
    y_bottom <- y_top
  }
  lines(x, y_total, col = "grey20", lwd = 1.5)
  if (!is.null(output_file)) dev.off()
  invisible(list(x = x, y_per_cluster = y_list, y_total = y_total, colors = cols))
}

##========================= Shlush 路径约定 =========================##

#' 从 Shlush BMM CSV 路径推导同一样本的 mutation table 路径
#' CSV 如 results_single/Patient-1/Patient-1_Diagnosis_bmm_params.csv
#' -> mutation_table/Patient-1_Dx_mutation_table.csv 或 _Rx_
mutation_table_path_from_bmm_csv_shlush <- function(csv_path, shlush_base_dir = NULL) {
  if (is.null(shlush_base_dir)) {
    # csv_path: .../Shlush_analysis/results_single/Patient-1/Patient-1_Diagnosis_bmm_params.csv
    parts <- strsplit(normalizePath(csv_path, mustWork = FALSE), .Platform$file.sep)[[1]]
    idx <- which(parts == "Shlush_analysis")
    if (length(idx) > 0) {
      shlush_base_dir <- do.call(file.path, as.list(parts[seq_len(idx)]))
    } else {
      # results_single/<id>/file.csv -> 上两级为 results_single，再上一级为 Shlush_analysis
      shlush_base_dir <- dirname(dirname(dirname(csv_path)))
    }
  }
  sample_label <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  if (grepl("_Diagnosis$", sample_label)) {
    patient_id <- sub("_Diagnosis$", "", sample_label)
    suffix <- "Dx"
  } else if (grepl("_Relapse$", sample_label)) {
    patient_id <- sub("_Relapse$", "", sample_label)
    suffix <- "Rx"
  } else {
    return(NULL)
  }
  file.path(shlush_base_dir, "mutation_table", paste0(patient_id, "_", suffix, "_mutation_table.csv"))
}

#' 从 BMM CSV 路径读取 truncal mean VAF（同目录下 cluster_summary .means）
read_truncal_mean_vaf_shlush <- function(csv_path, cluster_summary_ver = SHLUSH_CLUSTER_SUMMARY_VER) {
  csv_dir <- dirname(csv_path)
  sample_label <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  means_file <- paste0(sample_label, "_dp", SHLUSH_MIN_DEPTH, "_cluster_summary_v", cluster_summary_ver, ".means")
  path <- file.path(csv_dir, means_file)
  if (!file.exists(path)) return(NULL)
  d <- tryCatch(
    read.delim(path, header = FALSE, stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) NULL
  )
  if (is.null(d) || ncol(d) < 2 || nrow(d) < 1) return(NULL)
  v <- suppressWarnings(as.numeric(d[[2L]]))
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NULL)
  max(v)
}

##========================= 主入口：从 CSV 绘制 =========================##

#' 从 Shlush BMM 参数 CSV 绘制拟合线堆积图
#' @param csv_path  BMM 参数 CSV 路径（如 results_single/Patient-1/Patient-1_Diagnosis_bmm_params.csv）
#' @param output_file  输出 PNG 路径，NULL 则仅屏幕显示
#' @param mutation_table_path  可选，mutation table 路径；NULL 则从 csv_path 推导
#' @param shlush_base_dir  可选，Shlush_analysis 目录；NULL 则从 csv_path 推导
plot_bmm_fit_from_csv <- function(csv_path, output_file = NULL, mutation_table_path = NULL,
                                  shlush_base_dir = NULL,
                                  dim = 1L, n_grid = 999) {
  if (!file.exists(csv_path)) stop("文件不存在: ", csv_path)
  params <- read_bmm_params_csv(csv_path)
  curves <- compute_bmm_density_curves(params, n_grid = n_grid)
  sample_label <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  if (is.null(mutation_table_path)) {
    mutation_table_path <- mutation_table_path_from_bmm_csv_shlush(csv_path, shlush_base_dir)
  }
  hist_data <- compute_data_histogram(mutation_table_path)
  n_mutations <- if (!is.null(hist_data)) hist_data$n else NULL
  truncal_mean_vaf <- read_truncal_mean_vaf_shlush(csv_path)
  plot_bmm_stacked_one_dim(curves, dim = dim, output_file = output_file,
                           sample_id = sample_label, n_mutations = n_mutations,
                           truncal_mean_vaf = truncal_mean_vaf,
                           hist_data = hist_data)
}

##========================= 命令行 / 测试 =========================##

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1L) {
  csv_path <- .args[1L]
  output_file <- if (length(.args) >= 2L) .args[2L] else {
    sample_label <- sub("_bmm_params\\.csv$", "", basename(csv_path))
    shlush_base <- dirname(dirname(dirname(csv_path)))
    plot_dir <- file.path(shlush_base, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    file.path(plot_dir, paste0(sample_label, "_bmm_fit.png"))
  }
  plot_bmm_fit_from_csv(csv_path, output_file = output_file)
  cat("已写出:", output_file, "\n")
}
