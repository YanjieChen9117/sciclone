#!/usr/bin/env Rscript
################################################################################
# 从 BMM 参数 CSV 绘制“拟合线”堆积图，并叠加样本原始 VAF 直方图
#
# 参考 R/plots.R 中 1D 模型叠加方式，根据 CSV 中的 mu/alpha/nu/beta/pi
# 用 bmm 包计算各 cluster 的后验预测密度，绘制堆积图；背景为 mutation table
# 过滤（norm_read + tum_read >= 25）后按 1% bin 的密度直方图。
#
# 用法：
#   source("TCGA-COAD_analysis/plot_bmm_fit.R")
#   plot_bmm_fit_from_csv("TCGA-COAD_analysis/bmm_params/TCGA-A6-6141_bmm_params.csv",
#                         output_file = "TCGA-A6-6141_bmm_fit.png")
################################################################################

# 计算密度需用 bmm 包（与 sciClone 一致）
if (!requireNamespace("bmm", quietly = TRUE)) {
  stop("请先安装 bmm: install.packages('bmm') 或 devtools::install_github('genome/bmm')")
}
suppressPackageStartupMessages(library(bmm))

##========================= 原始数据直方图（mutation table）=========================##

#' 从样本 mutation table 计算 VAF 密度直方图（1% bin）
#' 过滤条件：norm_read + tum_read >= 25
#' bin：tum_vaf 按 1% 分箱（0–1%, 1–2%, …, 99–100%），密度 = 每 bin 比例（面积和为 1）
#' @param mutation_table_path  mutation_table CSV 路径（含 contig, norm_read, tum_read, tum_vaf）
#' @return list(breaks = 0:100, density = 长度 100 的向量, n = 过滤后突变数)，无文件时返回 NULL
compute_data_histogram <- function(mutation_table_path) {
  if (!file.exists(mutation_table_path)) return(NULL)
  d <- read.csv(mutation_table_path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("norm_read", "tum_read", "tum_vaf")
  if (!all(need %in% names(d))) return(NULL)
  d <- d[ (d$norm_read + d$tum_read) >= 25, , drop = FALSE ]
  n <- nrow(d)
  if (n == 0) return(list(breaks = 0:100, density = rep(0, 100), n = 0))
  vaf_pct <- as.numeric(d$tum_vaf)
  vaf_pct <- vaf_pct[!is.na(vaf_pct) & vaf_pct >= 0 & vaf_pct <= 100]
  n <- length(vaf_pct)
  if (n == 0) return(list(breaks = 0:100, density = rep(0, 100), n = 0))
  # 1% bins: [0,1), [1,2), ..., [99,100]; bin index 0..99
  bin_idx <- pmin(floor(vaf_pct), 99) + 1L   # 1-based, 1..100
  counts <- tabulate(bin_idx, nbins = 100)
  density <- counts / n   # 每 1% bin 的比例 = 密度（bin width = 1%）
  list(breaks = 0:100, density = density, n = n)
}

##========================= 解析 BMM 参数 CSV =========================##

#' 从 BMM 参数 CSV 读入并解析为矩阵/向量
#' CSV 格式：row_id = mu_1, alpha_1, nu_1, beta_1, pi；列为 cluster_1, cluster_2, ...
#' 返回 list(mu, alpha, nu, beta, pi)，其中 mu/alpha/nu/beta 为 维度 x 聚类数 矩阵，pi 为向量
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

#' 给定 BMM 参数，在 x 网格上计算各 cluster 的后验预测密度（与 clustering.R 一致）
#' x 为 [0,1] 区间；返回 list(x_vaf = x*100, y_per_cluster = list of length K, y_total = vector)
compute_bmm_density_curves <- function(params, n_grid = 999, num_samples_bmm = 100) {
  x <- seq(0, 1, length.out = n_grid + 2)[-c(1, n_grid + 2)]  # 避开 0 和 1
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

#' 获取 cluster 颜色（半透明），最多 6 个 cluster 使用指定色板
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
#' curves: list from compute_bmm_density_curves; dim: 使用第几维（通常 1）
#' hist_data: list(breaks, density, n) from compute_data_histogram，NULL 则不画直方图
#' 标题：main = sample_id；副标题：n = ..., truncal frequency = pi[1]
#' Y 轴为密度（与直方图一致）；BMM 密度从 [0,1] 转为 per-percent 与直方图同尺度
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
  # BMM 密度是 [0,1] 上的密度，x 用 0–100 时需除以 100 得到 per-percent 密度
  y_total <- y_total / 100
  y_list <- lapply(y_list, function(y) y / 100)

  # 数据直方图是 proportion（0–1），模型是 per-% 密度，需将直方图放大 100 倍与模型同尺度
  y_max_data <- 0
  if (!is.null(hist_data) && hist_data$n > 0) {
    y_max_data <- 100 * max(hist_data$density, na.rm = TRUE)
  }
  y_max_bmm <- max(y_total, na.rm = TRUE)
  y_max <- max(y_max_data, y_max_bmm, 1e-6)
  ylim <- c(0, y_max * 1.05)

  cols <- cluster_colors(K)
  # 2:1 长宽比（宽:高），去掉标题和副标题
  if (!is.null(output_file)) png(output_file, width = 6, height = 3, units = "in", res = 300, bg = "white")
  par(mar = c(2.5, 2.5, 1.5, 1.2), cex.lab = 0.85, cex.axis = 0.85)
  plot(NA, type = "n", xlim = xlim, ylim = ylim,
       xaxs = "i", yaxs = "i",
       xaxt = "n", yaxt = "n",  # 去掉 X/Y 轴数值刻度
       bty = "n",              # 去掉绘图边框
       main = "", xlab = xlab, ylab = ylab)
  # 背景：原始数据直方图（灰色），高度 ×100 与模型（per-%）同尺度
  if (!is.null(hist_data) && hist_data$n > 0) {
    br <- hist_data$breaks
    d <- hist_data$density
    for (i in seq_len(length(br) - 1)) {
      rect(br[i], 0, br[i + 1], 100 * d[i], col = "#999999", border = NA)
    }
  }
  # 堆积多边形：从下往上 cluster_1, cluster_2, ...
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

##========================= 主入口：从 CSV 绘制 =========================##

#' 从 BMM 参数 CSV 路径推导同一样本的 mutation_table 路径
#' 例如 bmm_params/TCGA-A6-6141_bmm_params.csv -> mutation_table/TCGA-A6-6141_mutation_table.csv
mutation_table_path_from_bmm_csv <- function(csv_path) {
  dir_bmm <- dirname(csv_path)
  base_dir <- dirname(dir_bmm)
  sample_id <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  file.path(base_dir, "mutation_table", paste0(sample_id, "_mutation_table.csv"))
}

#' 从样本 results_single 下的 cluster_summary .means 文件读取 truncal mean VAF
#' 即第二列（各 cluster 的 mean VAF）的最大值；文件不存在或无法解析时返回 NULL
#' 路径约定：results_single/<sample_id>/<sample_id>_dp25_cluster_summary_v1.means
read_truncal_mean_vaf <- function(csv_path, means_path = NULL) {
  if (!is.null(means_path)) {
    path <- means_path
  } else {
    dir_bmm <- dirname(csv_path)
    base_dir <- dirname(dir_bmm)
    sample_id <- sub("_bmm_params\\.csv$", "", basename(csv_path))
    path <- file.path(base_dir, "results_single", sample_id,
                      paste0(sample_id, "_dp25_cluster_summary_v1.means"))
  }
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

#' 从 BMM 参数 CSV 绘制拟合线堆积图，并叠加该样本的原始 VAF 直方图
#' @param csv_path  BMM 参数 CSV 路径（如 TCGA-A6-6141_bmm_params.csv）
#' @param output_file  输出 PNG 路径，NULL 则仅屏幕显示
#' @param mutation_table_path  可选，mutation table CSV 路径；NULL 则从 csv_path 推导
#' @param dim  使用第几维（单样本为 1）
#' @param n_grid  密度网格点数
plot_bmm_fit_from_csv <- function(csv_path, output_file = NULL, mutation_table_path = NULL,
                                  dim = 1L, n_grid = 999) {
  if (!file.exists(csv_path)) stop("文件不存在: ", csv_path)
  params <- read_bmm_params_csv(csv_path)
  curves <- compute_bmm_density_curves(params, n_grid = n_grid)
  sample_id <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  if (is.null(mutation_table_path)) {
    mutation_table_path <- mutation_table_path_from_bmm_csv(csv_path)
  }
  hist_data <- compute_data_histogram(mutation_table_path)
  n_mutations <- if (!is.null(hist_data)) hist_data$n else NULL
  truncal_mean_vaf <- read_truncal_mean_vaf(csv_path)
  plot_bmm_stacked_one_dim(curves, dim = dim, output_file = output_file,
                           sample_id = sample_id, n_mutations = n_mutations,
                           truncal_mean_vaf = truncal_mean_vaf,
                           hist_data = hist_data)
}

##========================= 命令行 / 测试 =========================##

# 若传入参数：第 1 个为 CSV 路径，第 2 个可选为输出 PNG 路径
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1L) {
  csv_path <- .args[1L]
  output_file <- if (length(.args) >= 2L) .args[2L] else sub("\\.csv$", "_bmm_fit.png", csv_path)
  plot_bmm_fit_from_csv(csv_path, output_file = output_file)
  cat("已写出:", output_file, "\n")
} else if (identical(Sys.getenv("RUN_PLOT_BMM_TEST"), "1")) {
  project_root <- "/Users/yanjiechen/Documents/Github/sciclone"
  csv_path <- file.path(project_root, "TCGA-COAD_analysis", "bmm_params", "TCGA-A6-6141_bmm_params.csv")
  if (file.exists(csv_path)) {
    out_png <- file.path(project_root, "TCGA-COAD_analysis", "bmm_params", "TCGA-A6-6141_bmm_fit.png")
    plot_bmm_fit_from_csv(csv_path, output_file = out_png)
    cat("已写出:", out_png, "\n")
  } else {
    cat("测试文件不存在:", csv_path, "\n")
  }
}
