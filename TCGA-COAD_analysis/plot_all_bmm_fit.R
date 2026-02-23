#!/usr/bin/env Rscript
################################################################################
# 为 bmm_params 目录下所有样本的 BMM 参数 CSV 批量绘制拟合图（PNG）
#
# 会反复调用 plot_bmm_fit.R 中的 plot_bmm_fit_from_csv，为每个 *_bmm_params.csv
# 在 plots 目录下生成对应的 *_bmm_fit.png。
#
# 用法（在项目根目录下）：
#   Rscript TCGA-COAD_analysis/plot_all_bmm_fit.R
# 或在 TCGA-COAD_analysis 目录下：
#   Rscript plot_all_bmm_fit.R
################################################################################

# 确定脚本所在目录（即 TCGA-COAD_analysis）
args <- commandArgs(trailingOnly = FALSE)
script_match <- grep("^--file=", args)
if (length(script_match) > 0) {
  script_path <- sub("^--file=", "", args[script_match])
  base_dir <- normalizePath(dirname(script_path), mustWork = FALSE)
} else {
  base_dir <- normalizePath("TCGA-COAD_analysis", mustWork = FALSE)
  if (!dir.exists(base_dir)) base_dir <- getwd()
}

bmm_dir <- file.path(base_dir, "bmm_params")
plots_dir <- file.path(base_dir, "plots")
plot_script <- file.path(base_dir, "plot_bmm_fit.R")

if (!file.exists(plot_script)) {
  stop("未找到 plot_bmm_fit.R: ", plot_script)
}
if (!dir.exists(bmm_dir)) {
  stop("未找到 bmm_params 目录: ", bmm_dir)
}
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
  cat("已创建目录:", plots_dir, "\n")
}

source(plot_script)

csv_files <- list.files(bmm_dir, pattern = "_bmm_params\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  cat("bmm_params 下没有 *_bmm_params.csv 文件。\n")
  quit(save = "no", status = 0)
}

cat("共", length(csv_files), "个样本，开始绘图...\n")
for (i in seq_along(csv_files)) {
  csv_path <- csv_files[i]
  sample_id <- sub("_bmm_params\\.csv$", "", basename(csv_path))
  output_file <- file.path(plots_dir, paste0(sample_id, "_bmm_fit.png"))
  tryCatch({
    plot_bmm_fit_from_csv(csv_path, output_file = output_file)
    cat("[", i, "/", length(csv_files), "] ", sample_id, " -> ", output_file, "\n", sep = "")
  }, error = function(e) {
    cat("[", i, "/", length(csv_files), "] ", sample_id, " 失败: ", conditionMessage(e), "\n", sep = "")
  })
}
cat("完成。\n")
