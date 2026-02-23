#!/usr/bin/env Rscript
################################################################################
# 从 TCGA-COAD raw_data 提取 mutation_table
#
# 1. 读取 sample_table.txt，取 Patient 列得到样本 ID
# 2. 对每个样本使用二倍体突变表：raw_data/<Patient>_1_1.txt
# 3. 按 Shlush 结构写出 mutation_table，列对应关系：
#    Chr --> contig (加 chr 前缀)
#    Pos --> position
#    ref_counts --> norm_read
#    alt_counts --> tum_read
#    alt_counts/(ref_counts+alt_counts)*100 --> tum_vaf
#    SYMBOL --> gene_name
#    Patient --> sample
#    Ref/Alt --> variant
#
# 输出：
#   <Patient>_mutation_table.csv
################################################################################

project_root <- "/Users/yanjiechen/Documents/Github/sciclone"
base_dir     <- file.path(project_root, "TCGA-COAD_analysis")
raw_dir      <- file.path(base_dir, "raw_data")
sample_file  <- file.path(raw_dir, "sample_table.txt")
out_dir      <- file.path(base_dir, "mutation_table")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("已创建 mutation_table 目录:", out_dir, "\n")
}

# 读取样本列表
sample_df <- read.delim(sample_file, stringsAsFactors = FALSE, check.names = FALSE)
patient_ids <- unique(sample_df$Patient)
patient_ids <- patient_ids[nzchar(patient_ids)]
cat(sprintf("样本表中共 %d 个患者 ID\n", length(patient_ids)))

n_ok <- 0
n_missing <- 0
n_empty <- 0

for (pid in patient_ids) {
  f <- file.path(raw_dir, paste0(pid, "_1_1.txt"))
  if (!file.exists(f)) {
    cat("  跳过（无二倍体表）:", pid, "\n")
    n_missing <- n_missing + 1
    next
  }

  # 读入二倍体突变表（tab 分隔，带引号）
  d <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"")

  if (nrow(d) == 0) {
    cat("  跳过（无行）:", pid, "\n")
    n_empty <- n_empty + 1
    next
  }

  # 必需列
  need <- c("Chr", "Pos", "Ref", "Alt", "ref_counts", "alt_counts", "Patient")
  miss <- setdiff(need, colnames(d))
  if (length(miss) > 0) {
    cat("  跳过", pid, "：缺少列 ", paste(miss, collapse = ", "), "\n")
    n_missing <- n_missing + 1
    next
  }

  # 数值列
  d$ref_counts <- as.numeric(d$ref_counts)
  d$alt_counts <- as.numeric(d$alt_counts)
  d$Pos        <- as.numeric(d$Pos)
  d[is.na(d$ref_counts), "ref_counts"] <- 0
  d[is.na(d$alt_counts), "alt_counts"] <- 0

  # 构建 mutation_table（与 Shlush 结构一致）
  total <- d$ref_counts + d$alt_counts
  tum_vaf <- ifelse(total > 0, 100 * d$alt_counts / total, 0)

  gene_col <- if ("SYMBOL" %in% colnames(d)) "SYMBOL" else "ID"
  gene_name <- d[[gene_col]]
  if (is.null(gene_name)) gene_name <- ""
  gene_name[is.na(gene_name) | !nzchar(gene_name)] <- ""

  out <- data.frame(
    contig    = paste0("chr", as.character(d$Chr)),
    position  = as.integer(d$Pos),
    norm_read = as.integer(d$ref_counts),
    tum_read  = as.integer(d$alt_counts),
    tum_vaf   = round(tum_vaf, 1),
    gene_name = as.character(gene_name),
    sample    = as.character(d$Patient),
    variant   = paste(d$Ref, d$Alt, sep = "/"),
    stringsAsFactors = FALSE
  )

  out_path <- file.path(out_dir, paste0(pid, "_mutation_table.csv"))
  write.csv(out, out_path, row.names = FALSE)
  n_ok <- n_ok + 1
  cat(sprintf("  已写 %s (%d 行)\n", basename(out_path), nrow(out)))
}

cat("\n完成: 成功 ", n_ok, "，无文件 ", n_missing, "，空表 ", n_empty, "\n")
cat("输出目录:", out_dir, "\n")
