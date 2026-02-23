#!/usr/bin/env Rscript
################################################################################
# Shlush：为 results_single 下所有样本的 BMM 参数 CSV 批量绘制拟合图（PNG）
#
# 扫描 Shlush_analysis/results_single 及其子目录（按患者分目录），
# 为每个 *_bmm_params.csv 在 Shlush_analysis/plots 下生成 <sample_label>_bmm_fit.png。
#
# 用法（在项目根目录下）：
#   Rscript Shlush_analysis/plot_all_bmm_fit.R
# 或在 Shlush_analysis 目录下：
#   Rscript plot_all_bmm_fit.R
################################################################################

args <- commandArgs(trailingOnly = FALSE)
script_match <- grep("^--file=", args)
if (length(script_match) > 0) {
  script_path <- sub("^--file=", "", args[script_match])
  base_dir <- normalizePath(dirname(script_path), mustWork = FALSE)
} else {
  base_dir <- normalizePath("Shlush_analysis", mustWork = FALSE)
  if (!dir.exists(base_dir)) base_dir <- getwd()
}

results_dir <- file.path(base_dir, "results_single")
plots_dir   <- file.path(base_dir, "plots")
plot_script <- file.path(base_dir, "plot_bmm_fit.R")

if (!file.exists(plot_script)) {
  stop("未找到 plot_bmm_fit.R: ", plot_script)
}
if (!dir.exists(results_dir)) {
  stop("未找到 results_single 目录: ", results_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
  cat("已创建目录:", plots_dir, "\n")
}

source(plot_script)

csv_files <- list.files(results_dir, pattern = "_bmm_params\\.csv$", full.names = TRUE, recursive = TRUE)
if (length(csv_files) == 0) {
  cat("results_single 下没有 *_bmm_params.csv 文件。\n")
  quit(save = "no", status = 0)
}

cat("共 ", length(csv_files), " 个样本，开始绘图...\n", sep = "")
for (i in seq_along(csv_files)) {
  csv_path <- csv_files[i]
  sample_label <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  output_file <- file.path(plots_dir, paste0(sample_label, "_bmm_fit.png"))
  tryCatch({
    plot_bmm_fit_from_csv(csv_path, output_file = output_file, shlush_base_dir = base_dir)
    cat("[", i, "/", length(csv_files), "] ", sample_label, " -> ", output_file, "\n", sep = "")
  }, error = function(e) {
    cat("[", i, "/", length(csv_files), "] ", sample_label, " 失败: ", conditionMessage(e), "\n", sep = "")
  })
}
cat("完成。\n")
