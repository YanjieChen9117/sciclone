## 配置根目录（请根据你的实际路径修改或直接 setwd 到 repo 根目录）
root_dir <- "/Users/yanjiechen/Documents/Github/sciclone"
results_dir <- file.path(root_dir, "TCGA-COAD_analysis", "results_single")
output_csv <- file.path(root_dir, "TCGA-COAD_analysis", "purity_estimates.csv")

## 列出所有 patient 子目录
patient_dirs <- list.dirs(results_dir, full.names = TRUE, recursive = FALSE)

## 用来存结果的列表
res_list <- list()

for (pd in patient_dirs) {
  ## patient id = 子目录名
  patient_id <- basename(pd)
  
  ## 查找 *_dp25_cluster_summary_v1.means 文件
  means_files <- list.files(
    pd,
    pattern = "_dp25_cluster_summary_v1\\.means$",
    full.names = TRUE
  )
  
  if (length(means_files) == 0) {
    ## 没有找到文件，记录 NA
    res_list[[length(res_list) + 1]] <- data.frame(
      patient = patient_id,
      truncal_purity = NA_real_
    )
    next
  }
  
  ## 如果有多个匹配文件，这里取第一个（如有需要可改为其它逻辑）
  f <- means_files[1]
  
  ## 更稳健地逐行解析 .means 文件
  lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) {
    res_list[[length(res_list) + 1]] <- data.frame(
      patient = patient_id,
      truncal_purity = NA_real_
    )
    next
  }
  
  ## 去掉空行和首尾空白
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  ## 只保留以 "cluster" 开头的行（排除第一行 patient id）
  cluster_lines <- grep("^cluster", lines, value = TRUE)
  
  if (length(cluster_lines) == 0) {
    res_list[[length(res_list) + 1]] <- data.frame(
      patient = patient_id,
      truncal_purity = NA_real_
    )
    next
  }
  
  ## 从每一行取出最后一个字段作为 cellular prevalence
  cp <- sapply(strsplit(cluster_lines, "\\s+"), function(x) {
    suppressWarnings(as.numeric(tail(x, 1)))
  })
  cp <- as.numeric(cp)
  cp <- cp[!is.na(cp)]
  
  if (length(cp) == 0) {
    res_list[[length(res_list) + 1]] <- data.frame(
      patient = patient_id,
      truncal_purity = NA_real_
    )
    next
  }
  
  ## 取 <= 0.5 中的最大值
  cp_le_05 <- cp[cp <= 0.5]
  if (length(cp_le_05) == 0) {
    truncal_purity <- NA_real_
  } else {
    truncal_cp <- max(cp_le_05, na.rm = TRUE)
    truncal_purity <- truncal_cp * 2
  }
  
  res_list[[length(res_list) + 1]] <- data.frame(
    patient = patient_id,
    truncal_purity = truncal_purity
  )
}

## 合并并输出 CSV
result_df <- do.call(rbind, res_list)

write.csv(result_df, file = output_csv, row.names = FALSE)

cat("完成：写出文件到", output_csv, "\n")