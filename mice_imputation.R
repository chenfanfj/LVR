############################################################
## 完整优化版: 多重插补 (MI) + 逆概率加权 (IPW)
## 版本: 2.0 - 已修复所有诊断问题
## 日期: 2024
############################################################

## ---------------------------
## 环境设置
## ---------------------------
rm(list = ls())
gc()

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║   多重插补 (MI) + IPW 完整流程 - 优化版 v2.0       ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## 加载包
packages <- c(
  "here", "readxl", "dplyr", "tidyr", "tibble", "stringr",
  "mice", "naniar", "VIM",
  "broom", "purrr", "pscl",
  "ggplot2", "corrplot", "gridExtra",
  "cobalt", "survey",
  "flextable", "officer", "gtsummary"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("安装包:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("✓ 所有包已加载\n\n")

## 创建输出目录
outputs_dir <- here::here("outputs")
plots_dir <- here::here("plots")
dir.create(outputs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

## ═══════════════════════════════════════════════════════
## 步骤 1: 数据读取
## ═══════════════════════════════════════════════════════

cat("【步骤 1】数据读取\n")

data_path <- here::here("data", "AMI_PCI_metal_clean3.xlsx")
dat_raw <- readxl::read_excel(
  path = data_path,
  col_names = TRUE,
  na = c("", "NA", "N/A", "#N/A", "NULL", "null", "-999", "-9999")
)

cat("  原始数据:", nrow(dat_raw), "行 ×", ncol(dat_raw), "列\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 数据清洗（增强版）
## ═══════════════════════════════════════════════════════

cat("【步骤 2】数据清洗（增强版）\n")

dat <- dat_raw

## ---------------------------
## 2.1 标准化缺失值
## ---------------------------
cat("  2.1 标准化缺失值\n")

for (col in names(dat)) {
  if (is.character(dat[[col]])) {
    dat[[col]][trimws(dat[[col]]) == ""] <- NA
  }
  if (is.numeric(dat[[col]])) {
    dat[[col]][dat[[col]] %in% c(-999, -9999, Inf, -Inf)] <- NA
  }
}

## ---------------------------
## 2.2 处理含不等号的变量（完整列表）
## ---------------------------
cat("  2.2 处理含不等号的变量\n")

## 扩展列表：所有可能含不等号的变量
vars_with_symbols <- c(
  ## 心肌标志物
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "NTproBNP_baseline", "NTproBNP_peak",
  "CK", "CKMB",
  
  ## 凝血功能
  "DD", "APTT", "PT_sec", "TT", "INR", "FIB",
  
  ## 甲状腺功能
  "FT3", "FT4", "S_TSH",
  
  ## 其他可能的实验室指标
  "IL_6", "CRP", "h_CT"
)

cleaning_log <- data.frame(
  variable = character(),
  n_with_symbols = integer(),
  examples = character(),
  stringsAsFactors = FALSE
)

for (var in intersect(vars_with_symbols, names(dat))) {
  
  if (is.character(dat[[var]])) {
    ## 检测含符号的值
    has_symbol <- grepl("＜|<|﹤|＞|>|﹥|≤|≥", dat[[var]], ignore.case = TRUE)
    n_with_symbols <- sum(has_symbol, na.rm = TRUE)
    
    if (n_with_symbols > 0) {
      examples <- paste(head(unique(dat[[var]][has_symbol]), 3), collapse = "; ")
      
      cleaning_log <- rbind(cleaning_log, data.frame(
        variable = var,
        n_with_symbols = n_with_symbols,
        examples = examples
      ))
      
      cat("    清洗:", var, "-", n_with_symbols, "个含符号\n")
      
      ## 移除所有类型的不等号
      dat[[var]] <- gsub("＜|<|﹤|≤", "", dat[[var]])
      dat[[var]] <- gsub("＞|>|﹥|≥", "", dat[[var]])
      dat[[var]] <- gsub(",", "", dat[[var]])  # 移除千位分隔符
      dat[[var]] <- gsub("\\s+", "", dat[[var]])  # 移除空格
      dat[[var]] <- trimws(dat[[var]])
      
      ## 转换为数值
      dat[[var]] <- suppressWarnings(as.numeric(dat[[var]]))
    }
  }
}

if (nrow(cleaning_log) > 0) {
  cat("\n    字符型数值清洗摘要:\n")
  print(cleaning_log, row.names = FALSE)
  write.csv(cleaning_log, 
            file.path(outputs_dir, "symbol_cleaning_log.csv"),
            row.names = FALSE)
}

## ---------------------------
## 2.3 创建目标变量
## ---------------------------
cat("\n  2.3 创建目标变量\n")

dat <- dat %>%
  mutate(
    has_fu_echo = factor(
      if_else(! is.na(LVEDV_fu), 1, 0),
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

fu_summary <- table(dat$has_fu_echo)
cat("    有随访:", fu_summary["Yes"], "例 (", 
    round(fu_summary["Yes"]/sum(fu_summary)*100, 1), "%)\n")
cat("    无随访:", fu_summary["No"], "例 (", 
    round(fu_summary["No"]/sum(fu_summary)*100, 1), "%)\n")

## ---------------------------
## 2.4 极端值处理（Winsorization）
## ---------------------------
cat("\n  2.4 极端值处理\n")

continuous_vars_for_winsorize <- c(
  "age", "height", "weight", "BMI", "HR", "SBP", "DBP",
  "WBC", "HGB", "PLT", "RBC", "HCT",
  "SCR", "EGFR", "BUN",
  "cTnI_baseline", "cTnIpeak", "NTproBNP_baseline", "NTproBNP_peak",
  "CK", "CKMB", "CRP", "IL_6",
  "EF_baseline", "LVEDV_baseline", "LVESV_baseline",
  "EF_fu", "LVEDV_fu", "LVESV_fu",
  "Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
  "DD", "APTT", "FIB"
)

winsorize_count <- 0

for (var in intersect(continuous_vars_for_winsorize, names(dat))) {
  if (is.numeric(dat[[var]])) {
    valid_vals <- dat[[var]][!is.na(dat[[var]]) & is.finite(dat[[var]])]
    
    if (length(valid_vals) > 10) {
      q1 <- quantile(valid_vals, 0.25)
      q3 <- quantile(valid_vals, 0.75)
      iqr <- q3 - q1
      
      lower_fence <- q1 - 3 * iqr
      upper_fence <- q3 + 3 * iqr
      
      n_outliers <- sum(valid_vals < lower_fence | valid_vals > upper_fence)
      
      if (n_outliers > 0 && n_outliers / length(valid_vals) < 0.05) {
        dat[[var]][dat[[var]] < lower_fence & ! is.na(dat[[var]])] <- lower_fence
        dat[[var]][dat[[var]] > upper_fence & !is.na(dat[[var]])] <- upper_fence
        winsorize_count <- winsorize_count + 1
      }
    }
  }
}

cat("    处理了", winsorize_count, "个变量的极端值\n")

cat("\n  ✓ 数据清洗完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 变量选择
## ═══════════════════════════════════════════════════════

cat("【步骤 3】插补变量选择\n")

## ---------------------------
## 3.1 定义变量分组
## ---------------------------

## 结局变量
outcome_vars <- c(
  "LVEDV_fu", "LVESV_fu", "EF_fu"
)

## 暴露变量（金属）
exposure_vars <- c(
  "Cu", "Zn", "Fe", "Se", "Pb",
  "Al", "As", "Cr", "Mn", "Ni", "Mo", 
  "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li"
)

## IPW 模型协变量（基于您的研究）
ipw_covariates <- c(
  ## 人口学
  "age", "gender", "resident", "Career",
  
  ## 合并症
  "DM", "hypertension", "smoking", "Cancer", "his_stroke", "AF",
  
  ## 心梗特征
  "STEMI", "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",
  "ST_dev", "ST_dep",
  
  ## 治疗
  "pPCI", "Stent_no", "Lesion_no",
  
  ## 并发症
  "Cardio_shock", "VF", "Cardiac_arrest_in",
  
  ## 严重程度
  "IN_killip", "OUT_killip", "GRACE_in", "GRACE_in_str",
  
  ## 基线超声（仅EF和LVEDV，移除LVESV避免共线性）
  "EF_baseline", "LVEDV_baseline",
  
  ## 心肌标志物
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "NTproBNP_baseline", "NTproBNP_peak",
  "CK", "CKMB",
  
  ## 血常规
  "WBC", "HGB", "PLT", "RBC", "HCT", "NE", "LY", "MO",
  
  ## 肾功能
  "SCR", "EGFR", "BUN",
  
  ## 炎症
  "CRP", "IL_6"
)

## 辅助变量
auxiliary_vars <- c(
  "height", "weight",  # 用于预测BMI
  "HR", "SBP", "DBP",
  "ALT", "AST", "TBIL",
  "CHOL", "TG", "LDL", "HDL",
  "GLU", "HbAlc",
  "DD", "FIB", "APTT",
  "K", "Na", "Ca", "Mg",
  "LAD", "LCX", "RCA", "LM",
  "TL_LAD", "TL_LCX", "TL_RCA",
  "Aspirin", "Statin", "ACEIorARB", "β_block"
)

## ---------------------------
## 3.2 合并并筛选
## ---------------------------

all_imputation_vars <- unique(c(
  "has_fu_echo",
  outcome_vars,
  exposure_vars,
  ipw_covariates,
  auxiliary_vars
))

## 明确排除的变量
exclude_vars <- c(
  ## 派生变量（稍后被动插补）
  "BMI", "ΔLVEDV", "AST_ALT", "BUN_CREA", "APROB_APROA",
  
  ## 避免共线性
  "LVESV_baseline",  # 与EF_baseline高度相关
  
  ## ID和其他非插补变量
  "ID"
)

imputation_vars <- setdiff(
  intersect(all_imputation_vars, names(dat)),
  exclude_vars
)

## 准备插补数据
dat_for_mi <- dat %>%
  select(ID, all_of(imputation_vars))

cat("  最终插补变量数:", length(imputation_vars), "\n")
cat("    - 结局变量  :", length(intersect(outcome_vars, imputation_vars)), "个\n")
cat("    - 暴露变量  :", length(intersect(exposure_vars, imputation_vars)), "个\n")
cat("    - IPW协变量 :", length(intersect(ipw_covariates, imputation_vars)), "个\n")
cat("    - 辅助变量  :", length(intersect(auxiliary_vars, imputation_vars)), "个\n\n")

## ---------------------------
## 3.3 缺失模式分析
## ---------------------------
cat("  缺失模式分析.. .\n")

miss_summary <- data.frame(
  variable = names(dat_for_mi),
  n_missing = colSums(is.na(dat_for_mi)),
  pct_missing = round(colMeans(is.na(dat_for_mi)) * 100, 2)
) %>%
  filter(n_missing > 0) %>%
  arrange(desc(pct_missing))

cat("    缺失变量数:", nrow(miss_summary), "/", ncol(dat_for_mi), "\n")

high_miss <- miss_summary %>% filter(pct_missing > 40)
if (nrow(high_miss) > 0) {
  cat("    ⚠ 缺失率 > 40% 的变量:\n")
  for (i in 1:min(5, nrow(high_miss))) {
    cat("      ", high_miss$variable[i], ":", high_miss$pct_missing[i], "%\n")
  }
}

write.csv(miss_summary, 
          file.path(outputs_dir, "missing_summary_for_mi.csv"),
          row.names = FALSE)

## 缺失模式可视化
png(file.path(plots_dir, "missing_pattern_before_mi.png"),
    width = 14, height = 10, units = "in", res = 300)

vis_miss(dat_for_mi %>% select(-ID), 
         sort_miss = TRUE, 
         warn_large_data = FALSE) +
  labs(title = "插补前的缺失模式")

dev.off()

cat("    ✓ 缺失分析完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 变量质量诊断（增强版 v2）
## ═══════════════════════════════════════════════════════

cat("【步骤 4】变量质量诊断（增强版）\n")

## ---------------------------
## 4.1 第一轮：检查明显的常数变量
## ---------------------------
cat("  4.1 检查明显常数变量\n")

constant_vars_round1 <- character()

for (var in names(dat_for_mi)) {
  if (var == "ID" || var == "has_fu_echo") next
  
  if (is.numeric(dat_for_mi[[var]])) {
    non_missing <- na.omit(dat_for_mi[[var]])
    
    if (length(non_missing) > 0) {
      n_unique <- length(unique(non_missing))
      sd_val <- sd(non_missing, na.rm = TRUE)
      
      ## 更严格的标准
      if (n_unique == 1 || sd_val == 0 || sd_val < 1e-10) {
        cat("    ✗ 移除:", var, 
            "(唯一值:", n_unique, 
            ", SD:", round(sd_val, 10), ")\n")
        constant_vars_round1 <- c(constant_vars_round1, var)
      }
    }
  } else if (is.factor(dat_for_mi[[var]])) {
    if (nlevels(dat_for_mi[[var]]) <= 1) {
      cat("    ✗ 移除因子:", var, 
          "(水平数:", nlevels(dat_for_mi[[var]]), ")\n")
      constant_vars_round1 <- c(constant_vars_round1, var)
    }
  }
}

if (length(constant_vars_round1) > 0) {
  dat_for_mi <- dat_for_mi %>% select(-all_of(constant_vars_round1))
  imputation_vars <- setdiff(imputation_vars, constant_vars_round1)
  cat("    第一轮移除:", length(constant_vars_round1), "个变量\n")
} else {
  cat("    ✓ 第一轮未发现常数变量\n")
}

## ---------------------------
## 4.2 第二轮：检查在完整案例中变为常数的变量
## ---------------------------
cat("\n  4.2 检查完整案例中的常数变量\n")

constant_vars_round2 <- character()

numeric_cols <- sapply(dat_for_mi, is.numeric)
numeric_vars <- names(dat_for_mi)[numeric_cols]

## 排除 ID
numeric_vars <- setdiff(numeric_vars, "ID")

for (var in numeric_vars) {
  ## 检查在完整案例中的标准差
  complete_cases_data <- dat_for_mi[[var]][complete.cases(dat_for_mi[, numeric_vars])]
  
  if (length(complete_cases_data) > 0) {
    sd_complete <- sd(complete_cases_data, na.rm = TRUE)
    n_unique_complete <- length(unique(complete_cases_data))
    
    ## 在完整案例中变为常数
    if (! is.na(sd_complete) && (sd_complete == 0 || sd_complete < 1e-10 || n_unique_complete == 1)) {
      cat("    ✗ 移除:", var, 
          "(完整案例 SD:", round(sd_complete, 10), 
          ", 唯一值:", n_unique_complete, ")\n")
      constant_vars_round2 <- c(constant_vars_round2, var)
    }
  }
}

if (length(constant_vars_round2) > 0) {
  dat_for_mi <- dat_for_mi %>% select(-all_of(constant_vars_round2))
  imputation_vars <- setdiff(imputation_vars, constant_vars_round2)
  cat("    第二轮移除:", length(constant_vars_round2), "个变量\n")
} else {
  cat("    ✓ 第二轮未发现常数变量\n")
}

## ---------------------------
## 4. 3 检查极低方差变量
## ---------------------------
cat("\n  4.3 检查极低方差变量\n")

low_variance_vars <- character()
low_variance_report <- data.frame(
  variable = character(),
  sd = numeric(),
  cv = numeric(),
  n_unique = integer(),
  stringsAsFactors = FALSE
)

for (var in numeric_vars) {
  if (var %in% c(constant_vars_round1, constant_vars_round2)) next
  
  if (var %in% names(dat_for_mi) && is.numeric(dat_for_mi[[var]])) {
    non_missing <- na.omit(dat_for_mi[[var]])
    
    if (length(non_missing) > 10) {
      sd_val <- sd(non_missing)
      mean_val <- mean(non_missing)
      cv <- if (mean_val != 0) abs(sd_val / mean_val) else Inf
      n_unique <- length(unique(non_missing))
      
      low_variance_report <- rbind(low_variance_report, data.frame(
        variable = var,
        sd = round(sd_val, 6),
        cv = round(cv, 6),
        n_unique = n_unique
      ))
      
      ## 标记极低方差（变异系数 < 0.01 且唯一值 <= 3）
      if (cv < 0.01 && n_unique <= 3) {
        cat("    ⚠ 极低方差:", var, 
            "(CV:", round(cv, 4), 
            ", 唯一值:", n_unique, ")\n")
        low_variance_vars <- c(low_variance_vars, var)
      }
    }
  }
}

if (nrow(low_variance_report) > 0) {
  low_variance_report <- low_variance_report %>% arrange(cv)
  write.csv(low_variance_report, 
            file.path(outputs_dir, "variable_variance_detailed.csv"),
            row.names = FALSE)
  
  cat("    变异系数最小的10个变量:\n")
  print(head(low_variance_report, 10), row.names = FALSE)
}

if (length(low_variance_vars) > 0) {
  cat("\n    建议移除极低方差变量:", paste(low_variance_vars, collapse = ", "), "\n")
  cat("    是否移除？（建议保留，除非确认为数据错误）\n")
  
  ## 保存警告列表
  write.csv(data.frame(variable = low_variance_vars, reason = "极低方差"),
            file.path(outputs_dir, "low_variance_warning.csv"),
            row.names = FALSE)
}

## ---------------------------
## 4.4 安全的相关性检查（修复版）
## ---------------------------
cat("\n  4.4 检查完美共线性\n")

numeric_cols_final <- sapply(dat_for_mi, is.numeric)
numeric_vars_final <- names(dat_for_mi)[numeric_cols_final]
numeric_vars_final <- setdiff(numeric_vars_final, 
                              c("ID", constant_vars_round1, constant_vars_round2))

if (length(numeric_vars_final) > 1) {
  numeric_data_final <- dat_for_mi[, numeric_vars_final, drop = FALSE]
  
  ## 预先检查每个变量的标准差
  cat("    预检查各变量标准差.. .\n")
  
  sd_check_results <- data.frame(
    variable = character(),
    n_total = integer(),
    n_missing = integer(),
    n_valid = integer(),
    sd = numeric(),
    stringsAsFactors = FALSE
  )
  
  zero_sd_vars <- character()
  
  for (var in names(numeric_data_final)) {
    non_na <- na.omit(numeric_data_final[[var]])
    n_total <- length(numeric_data_final[[var]])
    n_missing <- sum(is.na(numeric_data_final[[var]]))
    n_valid <- length(non_na)
    
    if (n_valid > 0) {
      sd_val <- sd(non_na)
      
      sd_check_results <- rbind(sd_check_results, data.frame(
        variable = var,
        n_total = n_total,
        n_missing = n_missing,
        n_valid = n_valid,
        sd = sd_val
      ))
      
      ## 检测零标准差
      if (is.na(sd_val) || sd_val < 1e-10) {
        cat("      ✗", var, ": SD =", 
            if (is.na(sd_val)) "NA" else sprintf("%.2e", sd_val), 
            "(n_valid =", n_valid, ")\n")
        zero_sd_vars <- c(zero_sd_vars, var)
      }
    } else {
      cat("      ⊗", var, ": 无有效值\n")
      zero_sd_vars <- c(zero_sd_vars, var)
    }
  }
  
  ## 保存标准差检查结果
  if (nrow(sd_check_results) > 0) {
    sd_check_results <- sd_check_results %>% arrange(sd)
    write.csv(sd_check_results, 
              file.path(outputs_dir, "sd_check_before_correlation.csv"),
              row.names = FALSE)
  }
  
  ## 移除零标准差变量
  if (length(zero_sd_vars) > 0) {
    cat("\n    移除零标准差变量 (", length(zero_sd_vars), "个):\n")
    cat("     ", paste(zero_sd_vars, collapse = ", "), "\n")
    
    dat_for_mi <- dat_for_mi %>% select(-all_of(zero_sd_vars))
    imputation_vars <- setdiff(imputation_vars, zero_sd_vars)
    
    ## 更新数值数据
    numeric_data_final <- dat_for_mi[, setdiff(numeric_vars_final, zero_sd_vars), drop = FALSE]
  }
  
  ## 现在安全地计算相关性
  if (ncol(numeric_data_final) > 1) {
    cat("\n    计算相关性矩阵 (", ncol(numeric_data_final), "个变量).. .\n")
    
    tryCatch({
      ## 正确的参数：无空格！
      cor_matrix <- cor(numeric_data_final, 
                        use = "pairwise.complete.obs")  # ← 修复：移除空格
      
      cat("      ✓ 相关性矩阵计算成功\n")
      
      ## 检查完美相关
      perfect_pairs <- character()
      high_cor_pairs <- list()
      
      for (i in 1:(ncol(cor_matrix) - 1)) {
        for (j in (i + 1):ncol(cor_matrix)) {
          r <- cor_matrix[i, j]
          
          if (! is.na(r) && abs(r) > 0.995) {
            var1 <- rownames(cor_matrix)[i]
            var2 <- colnames(cor_matrix)[j]
            
            cat("      ⚠ 高度相关:", var1, "与", var2, 
                "(r =", round(r, 4), ")\n")
            
            ## 选择保留缺失较少的
            miss1 <- mean(is.na(dat_for_mi[[var1]]))
            miss2 <- mean(is.na(dat_for_mi[[var2]]))
            
            var_to_remove <- if (miss1 > miss2) var1 else var2
            var_to_keep <- if (miss1 > miss2) var2 else var1
            
            cat("        → 保留:", var_to_keep, 
                "(缺失:", round(min(miss1, miss2)*100, 1), "%)\n")
            cat("        → 移除:", var_to_remove, 
                "(缺失:", round(max(miss1, miss2)*100, 1), "%)\n")
            
            perfect_pairs <- c(perfect_pairs, var_to_remove)
            
            high_cor_pairs[[length(high_cor_pairs) + 1]] <- list(
              var1 = var1,
              var2 = var2,
              correlation = r,
              removed = var_to_remove
            )
          }
        }
      }
      
      ## 移除高度相关变量
      if (length(perfect_pairs) > 0) {
        perfect_pairs <- unique(perfect_pairs)
        
        cat("\n      移除高度相关变量 (", length(perfect_pairs), "个):\n")
        cat("       ", paste(perfect_pairs, collapse = ", "), "\n")
        
        dat_for_mi <- dat_for_mi %>% select(-all_of(perfect_pairs))
        imputation_vars <- setdiff(imputation_vars, perfect_pairs)
        
        ## 保存相关性报告
        if (length(high_cor_pairs) > 0) {
          cor_report <- do.call(rbind, lapply(high_cor_pairs, function(x) {
            data.frame(
              var1 = x$var1,
              var2 = x$var2,
              correlation = round(x$correlation, 4),
              removed = x$removed,
              stringsAsFactors = FALSE
            )
          }))
          
          write.csv(cor_report, 
                    file.path(outputs_dir, "high_correlation_pairs. csv"),
                    row.names = FALSE)
        }
      } else {
        cat("      ✓ 未发现完美共线性 (|r| > 0.995)\n")
      }
      
    }, warning = function(w) {
      cat("      ⚠ 相关性计算警告:", conditionMessage(w), "\n")
      
      ## 尝试逐个变量检查
      cat("      → 逐个检查变量标准差...\n")
      
      for (var in names(numeric_data_final)) {
        vals <- na.omit(numeric_data_final[[var]])
        if (length(vals) > 0) {
          sd_val <- sd(vals)
          if (is.na(sd_val) || sd_val < 1e-8) {
            cat("        ✗ 发现低方差:", var, 
                "SD =", sprintf("%.2e", sd_val), "\n")
            
            dat_for_mi <- dat_for_mi %>% select(-all_of(var))
            imputation_vars <- setdiff(imputation_vars, var)
          }
        }
      }
      
    }, error = function(e) {
      cat("      ✗ 相关性计算失败:", conditionMessage(e), "\n")
      cat("      → 检查参数和数据类型\n")
      
      ## 诊断信息
      cat("\n      诊断信息:\n")
      cat("        数据维度 :", nrow(numeric_data_final), "×", ncol(numeric_data_final), "\n")
      cat("        数据类型 :", class(numeric_data_final), "\n")
      cat("        变量类型 :", paste(unique(sapply(numeric_data_final, class)), collapse = ", "), "\n")
      
      ## 检查是否所有列都是数值型
      non_numeric <- names(numeric_data_final)[! sapply(numeric_data_final, is.numeric)]
      if (length(non_numeric) > 0) {
        cat("        ⚠ 非数值列:", paste(non_numeric, collapse = ", "), "\n")
      }
    })
    
  } else {
    cat("    ⊙ 数值变量不足2个，跳过相关性检查\n")
  }
  
} else {
  cat("    ⊙ 数值变量不足，跳过共线性检查\n")
}

cat("\n")

## ---------------------------
## 4.5 生成诊断摘要
## ---------------------------
cat("\n  4.5 变量清理摘要\n")

total_removed <- length(c(constant_vars_round1, constant_vars_round2, zero_sd_vars))

cat("    移除变量总数   :", total_removed, "\n")
cat("      - 第一轮常数 :", length(constant_vars_round1), "\n")
cat("      - 第二轮常数 :", length(constant_vars_round2), "\n")
if (exists("zero_sd_vars")) {
  cat("      - 零标准差   :", length(zero_sd_vars), "\n")
}
if (exists("perfect_pairs")) {
  cat("      - 高度相关   :", length(unique(perfect_pairs)), "\n")
}

cat("    保留变量数     :", ncol(dat_for_mi) - 1, "个（排除ID）\n")

## 保存清理日志
all_removed_vars <- unique(c(constant_vars_round1, constant_vars_round2, 
                             zero_sd_vars, perfect_pairs))

if (length(all_removed_vars) > 0) {
  removed_log <- data.frame(
    variable = all_removed_vars,
    reason = sapply(all_removed_vars, function(v) {
      if (v %in% constant_vars_round1) return("第一轮常数")
      if (v %in% constant_vars_round2) return("完整案例中常数")
      if (exists("zero_sd_vars") && v %in% zero_sd_vars) return("零标准差")
      if (exists("perfect_pairs") && v %in% perfect_pairs) return("高度相关")
      return("其他")
    })
  )
  
  write.csv(removed_log, 
            file.path(outputs_dir, "removed_variables_log.csv"),
            row.names = FALSE)
  
  cat("\n    移除变量清单已保存: removed_variables_log.csv\n")
}

cat("\n  ✓ 变量质量诊断完成\n\n")

## 保存最终清理后的数据
saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_cleaned.rds"))


## ═══════════════════════════════════════════════════════
## 步骤 5: 多重插补配置（完全修复版）
## ═══════════════════════════════════════════════════════

cat("【步骤 5】配置多重插补（完全修复版）\n")

## ---------------------------
## 5.1 初始化
## ---------------------------
cat("  5.1 初始化 MICE\n")

ini <- mice(dat_for_mi, maxit = 0, print = FALSE)

if (! is.null(ini$loggedEvents) && nrow(ini$loggedEvents) > 0) {
  cat("    ⚠ 初始化警告:\n")
  print(ini$loggedEvents)
}

meth <- ini$method
pred <- ini$predictorMatrix

cat("    ✓ 初始化完成\n")

## ---------------------------
## 5. 2 变量类型转换 + 插补方法
## ---------------------------
cat("\n  5.2 变量类型转换与插补方法配置\n")

## ID 和目标变量不插补
meth["ID"] <- ""
pred[, "ID"] <- 0
pred["ID", ] <- 0
meth["has_fu_echo"] <- ""

## 二元变量转换为 factor
binary_vars <- c(
  "gender", "DM", "hypertension", "smoking", "Cancer", 
  "his_stroke", "AF", "pPCI", "VF", "Cardio_shock",
  "STEMI", "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",
  "ST_dev", "ST_dep", "Cardiac_arrest_in",
  "LAD", "LCX", "RCA", "LM", "TL_LAD", "TL_LCX", "TL_RCA",
  "resident", "Career",
  "Aspirin", "Statin", "ACEIorARB", "β_block"
)

cat("    转换二元变量:\n")
n_converted_binary <- 0

for (var in intersect(binary_vars, names(dat_for_mi))) {
  if (is.numeric(dat_for_mi[[var]])) {
    dat_for_mi[[var]] <- factor(dat_for_mi[[var]], 
                                levels = c(0, 1), 
                                labels = c("No", "Yes"))
    n_converted_binary <- n_converted_binary + 1
  }
  if (var %in% names(meth) && meth[var] != "") {
    meth[var] <- "logreg"
  }
}

cat("      转换了", n_converted_binary, "个二元变量\n")

## 有序变量转换
ordinal_vars <- c("IN_killip", "OUT_killip", "GRACE_in_str")

cat("    转换有序变量:\n")
n_converted_ordinal <- 0

for (var in intersect(ordinal_vars, names(dat_for_mi))) {
  if (is.numeric(dat_for_mi[[var]])) {
    unique_vals <- sort(unique(na.omit(dat_for_mi[[var]])))
    dat_for_mi[[var]] <- factor(dat_for_mi[[var]], 
                                levels = unique_vals, 
                                ordered = TRUE)
    n_converted_ordinal <- n_converted_ordinal + 1
  }
  if (var %in% names(meth) && meth[var] != "") {
    meth[var] <- "polr"
  }
}

cat("      转换了", n_converted_ordinal, "个有序变量\n")

## 连续变量用 PMM
continuous_vars <- setdiff(names(meth)[meth != ""], 
                           c(binary_vars, ordinal_vars))

for (var in continuous_vars) {
  if (meth[var] != "") {
    meth[var] <- "pmm"
  }
}

cat("\n    插补方法统计:\n")
cat("      PMM    :", sum(meth == "pmm"), "\n")
cat("      logreg :", sum(meth == "logreg"), "\n")
cat("      polr   :", sum(meth == "polr"), "\n")

## ---------------------------
## 5.3 被动插补（确保数值型）
## ---------------------------
cat("\n  5.3 配置被动插补\n")

## BMI
if (all(c("height", "weight") %in% names(dat_for_mi))) {
  
  ## 强制数值型
  dat_for_mi$height <- as.numeric(as.character(dat_for_mi$height))
  dat_for_mi$weight <- as.numeric(as.character(dat_for_mi$weight))
  
  if (! "BMI" %in% names(dat_for_mi)) {
    dat_for_mi$BMI <- NA_real_  ## 初始化为数值型 NA
  }
  
  meth["BMI"] <- "~I(weight / (height/100)^2)"
  cat("    ✓ BMI 被动插补\n")
}

## ΔLVEDV
if (all(c("LVEDV_fu", "LVEDV_baseline") %in% names(dat_for_mi))) {
  
  ## 强制数值型
  dat_for_mi$LVEDV_fu <- as.numeric(as.character(dat_for_mi$LVEDV_fu))
  dat_for_mi$LVEDV_baseline <- as.numeric(as.character(dat_for_mi$LVEDV_baseline))
  
  if (! "ΔLVEDV" %in% names(dat_for_mi)) {
    dat_for_mi$ΔLVEDV <- NA_real_
  }
  
  meth["ΔLVEDV"] <- "~I((LVEDV_fu - LVEDV_baseline) / LVEDV_baseline * 100)"
  cat("    ✓ ΔLVEDV 被动插补\n")
}

## 被动变量不预测
passive_vars <- names(meth)[grepl("^~", meth)]
for (var in passive_vars) {
  if (var %in% rownames(pred)) {
    pred[var, ] <- 0
  }
}

## ---------------------------
## 5.4 重新同步（关键步骤）
## ---------------------------
cat("\n  5. 4 同步配置\n")

## 因为添加了 BMI/ΔLVEDV，需要重新初始化
ini_final <- mice(dat_for_mi, maxit = 0, print = FALSE)

## 获取新的维度
meth_final <- ini_final$method
pred_final <- ini_final$predictorMatrix

## 恢复所有自定义设置
meth_final["ID"] <- ""
pred_final[, "ID"] <- 0
pred_final["ID", ] <- 0
meth_final["has_fu_echo"] <- ""

for (var in intersect(binary_vars, names(meth_final))) {
  if (meth_final[var] != "") meth_final[var] <- "logreg"
}

for (var in intersect(ordinal_vars, names(meth_final))) {
  if (meth_final[var] != "") meth_final[var] <- "polr"
}

if ("BMI" %in% names(meth_final)) {
  meth_final["BMI"] <- "~I(weight / (height/100)^2)"
  pred_final["BMI", ] <- 0
}

if ("ΔLVEDV" %in% names(meth_final)) {
  meth_final["ΔLVEDV"] <- "~I((LVEDV_fu - LVEDV_baseline) / LVEDV_baseline * 100)"
  pred_final["ΔLVEDV", ] <- 0
}

meth <- meth_final
pred <- pred_final

cat("    同步后: 数据", ncol(dat_for_mi), "列, meth", length(meth), "个\n")

## ---------------------------
## 5.5 检查并处理罕见事件
## ---------------------------
cat("\n  5.5 检查罕见事件变量\n")

for (var in names(dat_for_mi)) {
  if (is.factor(dat_for_mi[[var]]) && 
      nlevels(dat_for_mi[[var]]) == 2 && 
      var %in% names(meth) &&
      meth[var] == "logreg") {
    
    tab <- table(na.omit(dat_for_mi[[var]]))
    
    if (min(tab) < 10) {
      cat("    ⚠", var, ": 最少事件", min(tab), "→ 改用 PMM\n")
      
      ## 转回数值使用 PMM
      dat_for_mi[[var]] <- as.numeric(dat_for_mi[[var]]) - 1
      meth[var] <- "pmm"
    }
  }
}

## ═══════════════════════════════════════════════════════
## 修复插补配置（针对日志事件）
## ═══════════════════════════════════════════════════════

cat("\n【修复】根据日志事件调整插补配置\n")

## ---------------------------
## 修复 1: Lesion_no
## ---------------------------
cat("\n  修复 1: Lesion_no\n")

if ("Lesion_no" %in% names(dat_for_mi)) {
  
  ## 检查当前数据类型
  cat("    当前类型:", class(dat_for_mi$Lesion_no), "\n")
  cat("    唯一值  :", paste(sort(unique(na.omit(dat_for_mi$Lesion_no))), collapse = ", "), "\n")
  
  ## 选项 1: 改用 PMM（最稳健，推荐）
  if (is.factor(dat_for_mi$Lesion_no)) {
    dat_for_mi$Lesion_no <- as.numeric(as.character(dat_for_mi$Lesion_no))
    cat("    → 转换为数值型\n")
  }
  
  meth["Lesion_no"] <- "pmm"
  cat("    → 插补方法改为 PMM\n")
  
  ## 选项 2: 使用有序逻辑回归（如果明确是有序变量）
  # dat_for_mi$Lesion_no <- factor(dat_for_mi$Lesion_no, 
  #                                levels = c(0, 1, 2, 3),
  #                                ordered = TRUE)
  # meth["Lesion_no"] <- "polr"
  # cat("    → 转换为有序因子，使用 polr\n")
}

## ---------------------------
## 修复 2: has_fu_echo（关键！）
## ---------------------------
cat("\n  修复 2: has_fu_echo\n")

if ("has_fu_echo" %in% names(dat_for_mi)) {
  
  ## 检查缺失值
  n_missing <- sum(is.na(dat_for_mi$has_fu_echo))
  
  if (n_missing > 0) {
    cat("    ⚠ 检测到", n_missing, "个缺失值，重新创建\n")
    
    ## 基于 LVEDV_fu 重新定义
    dat_for_mi$has_fu_echo <- factor(
      if_else(!  is.na(dat_for_mi$LVEDV_fu), 1, 0),
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
    
    cat("    → 已重新创建（基于 LVEDV_fu）\n")
  }
  
  ## 强制设置为不插补
  meth["has_fu_echo"] <- ""
  
  ## 从预测矩阵中移除（不应作为预测因子）
  pred["has_fu_echo", ] <- 0
  pred[, "has_fu_echo"] <- 0
  
  cat("    ✓ has_fu_echo 已设置为不插补\n")
  cat("    ✓ 已从预测矩阵中移除\n")
}

## ---------------------------
## 验证修复
## ---------------------------
cat("\n  验证修复结果:\n")

cat("    Lesion_no 插补方法  :", meth["Lesion_no"], "\n")
cat("    Lesion_no 数据类型  :", class(dat_for_mi$Lesion_no), "\n")
cat("    has_fu_echo 插补方法:", meth["has_fu_echo"], "(应为空字符串)\n")
cat("    has_fu_echo 缺失数  :", sum(is.na(dat_for_mi$has_fu_echo)), "(应为 0)\n")

## 保存修复后的配置
saveRDS(list(
  data = dat_for_mi,
  method = meth,
  predictorMatrix = pred
), file. path(outputs_dir, "mice_config_after_fix.rds"))

cat("\n  ✓ 配置修复完成\n\n")

## ---------------------------
## 5.6 执行插补
## ---------------------------
cat("\n  5.6 执行多重插补\n")

## 最终验证
if (ncol(dat_for_mi) != length(meth)) {
  stop("维度不匹配: data=", ncol(dat_for_mi), ", meth=", length(meth))
}

m <- 20
maxit <- 15
seed <- 20240101

set.seed(seed)

cat("    参数: m=", m, ", maxit=", maxit, "\n")
cat("    开始插补..  .\n")

imp <- mice(
  dat_for_mi,
  m = m,
  maxit = maxit,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE,  # 显示进度
  seed = seed
)

cat("    ✓ 插补完成\n\n")

## 保存
saveRDS(imp, file.path(outputs_dir, "mice_imputation_object.  rds"))
saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_with_factors. rds"))

cat("  ✓ 多重插补完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 5. 7: 检查插补日志事件（新增）
## ═══════════════════════════════════════════════════════

cat("【步骤 5.7】检查插补日志事件\n")

if (! is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0) {
  
  cat("  检测到", nrow(imp$loggedEvents), "个日志事件\n\n")
  
  ## 分析日志事件类型
  event_summary <- imp$loggedEvents %>%
    group_by(out) %>%
    summarise(
      n_events = n(),
      meth = unique(meth),
      .groups = "drop"
    ) %>%
    arrange(desc(n_events))
  
  cat("  日志事件摘要（按变量）:\n")
  print(head(event_summary, 20), row.names = FALSE)
  
  ## 保存完整日志
  write.csv(imp$loggedEvents, 
            file.path(outputs_dir, "mice_logged_events_full.csv"),
            row.names = FALSE)
  
  ## 识别问题变量
  problem_vars <- event_summary %>%
    filter(n_events > 10) %>%  # 事件数 > 10 的变量
    pull(out)
  
  if (length(problem_vars) > 0) {
    cat("\n  问题变量 (事件数 > 10):\n")
    for (pv in problem_vars) {
      n_ev <- sum(imp$loggedEvents$out == pv)
      meth_pv <- unique(imp$loggedEvents$meth[imp$loggedEvents$out == pv])
      cat("    -", pv, ":", n_ev, "次 (方法:", meth_pv, ")\n")
    }
    
    ## 建议处理
    cat("\n  建议处理方案:\n")
    cat("    1. 将这些变量改用 PMM 插补\n")
    cat("    2.  检查这些变量的缺失模式和预测因子\n")
    cat("    3. 考虑移除事件数过多的变量\n")
  }
  
} else {
  cat("  ✓ 无日志事件\n")
}

cat("\n")
## ═══════════════════════════════════════════════════════
## 步骤 5.8: 诊断插补结果（新增）
## ═══════════════════════════════════════════════════════

cat("【步骤 5.8】诊断插补结果\n")

## 获取第一个插补数据集
imp_data_1 <- complete(imp, 1)

## 检查每个变量的插补情况
imputation_diagnostic <- data.frame(
  variable = character(),
  original_missing = integer(),
  after_imputation_missing = integer(),
  imputed_count = integer(),
  imputation_rate = numeric(),
  method = character(),
  stringsAsFactors = FALSE
)

for (var in names(dat_for_mi)) {
  
  if (var == "ID") next
  
  original_na <- sum(is.na(dat_for_mi[[var]]))
  after_na <- sum(is.na(imp_data_1[[var]]))
  imputed <- original_na - after_na
  imp_rate <- if (original_na > 0) imputed / original_na * 100 else 100
  
  imputation_diagnostic <- rbind(imputation_diagnostic, data.frame(
    variable = var,
    original_missing = original_na,
    after_imputation_missing = after_na,
    imputed_count = imputed,
    imputation_rate = round(imp_rate, 1),
    method = if (var %in% names(imp$method)) imp$method[var] else ""
  ))
}

## 排序：插补率最低的在前
imputation_diagnostic <- imputation_diagnostic %>%
  arrange(imputation_rate, desc(original_missing))

cat("\n  插补效果摘要:\n")
cat("    完全插补变量:", sum(imputation_diagnostic$imputation_rate == 100), "个\n")
cat("    部分插补变量:", sum(imputation_diagnostic$imputation_rate > 0 & 
                         imputation_diagnostic$imputation_rate < 100), "个\n")
cat("    未插补变量  :", sum(imputation_diagnostic$imputation_rate == 0 & 
                          imputation_diagnostic$original_missing > 0), "个\n")

## 显示插补失败或不完全的变量
failed_imputation <- imputation_diagnostic %>%
  filter(original_missing > 0, imputation_rate < 100)

if (nrow(failed_imputation) > 0) {
  cat("\n  插补不完全的变量:\n")
  print(failed_imputation, row.names = FALSE)
  
  ## 保存
  write.csv(failed_imputation, 
            file.path(outputs_dir, "failed_imputation_variables.csv"),
            row.names = FALSE)
  
  cat("\n  ⚠ 这些变量可能导致密度图错误\n")
}

## 保存完整诊断
write.csv(imputation_diagnostic, 
          file.path(outputs_dir, "imputation_diagnostic_report.csv"),
          row.names = FALSE)

cat("\n  ✓ 插补诊断完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 插补质量诊断
## ═══════════════════════════════════════════════════════

cat("【步骤 6】插补质量诊断\n")

## 收敛性
png(file.path(plots_dir, "mi_convergence.png"),
    width = 16, height = 14, units = "in", res = 300)
plot(imp, layout = c(5, 5))
dev.off()
cat("  ✓ 收敛性图已保存\n")

## 密度图
key_vars <- c("age", "BMI", "EF_baseline", "LVEDV_baseline", 
              "Cu", "Zn", "NTproBNP_peak", "CKMB", "DD", "APTT")
key_vars <- intersect(key_vars, names(dat_for_mi))

if (length(key_vars) > 0) {
  png(file.path(plots_dir, "mi_density_check.png"),
      width = 14, height = 12, units = "in", res = 300)
  
  densityplot(imp, 
              as.formula(paste("~", paste(key_vars, collapse = " + "))),
              layout = c(3, ceiling(length(key_vars)/3)))
  
  dev.off()
  cat("  ✓ 密度图已保存\n\n")
}


## ═══════════════════════════════════════════════════════
## 步骤 7: 生成综合报告
## ═══════════════════════════════════════════════════════

cat("【步骤 7】生成综合报告\n")

report_path <- file.path(outputs_dir, "MI_Report.txt")

sink(report_path)

cat("═══════════════════════════════════════════════════════\n")
cat("    多重插补 (MI)分析完整报告\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. 数据清洗\n")
cat("─────────────────────────────────────\n")
cat("原始样本数        :", nrow(dat_raw), "\n")
cat("清洗后样本数      :", nrow(dat), "\n")
cat("处理含符号变量    :", nrow(cleaning_log), "个\n")
cat("Winsorize 变量数  :", winsorize_count, "个\n\n")

cat("2. 变量选择\n")
cat("─────────────────────────────────────\n")
cat("插补变量数        :", length(imputation_vars), "\n")
cat("  - 结局变量      :", length(intersect(outcome_vars, imputation_vars)), "个\n")
cat("  - 暴露变量      :", length(intersect(exposure_vars, imputation_vars)), "个\n")
cat("  - IPW协变量     :", length(intersect(ipw_covariates, imputation_vars)), "个\n")
cat("  - 辅助变量      :", length(intersect(auxiliary_vars, imputation_vars)), "个\n\n")

cat("3. 缺失数据\n")
cat("─────────────────────────────────────\n")
cat("缺失变量数        :", nrow(miss_summary), "\n")
cat("平均缺失率        :", round(mean(miss_summary$pct_missing), 1), "%\n")
cat("最高缺失率        :", max(miss_summary$pct_missing), "% (", 
    miss_summary$variable[which. max(miss_summary$pct_missing)], ")\n\n")

cat("4. 多重插补\n")
cat("─────────────────────────────────────\n")
cat("插补数据集数      : m =", m, "\n")
cat("最大迭代数        : maxit =", maxit, "\n")
cat("插补方法          :\n")
cat("  - PMM           :", sum(meth == "pmm"), "个\n")
cat("  - logreg        :", sum(meth == "logreg"), "个\n")
cat("  - polr          :", sum(meth == "polr"), "个\n")
cat("  - 被动插补      :", sum(grepl("^~", meth)), "个\n\n")

cat("5. 输出文件\n")
cat("─────────────────────────────────────\n")
cat("插补对象:\n")
cat("  - mice_imputation_object.rds\n")
cat("结果文件:\n")
cat("  - ps_model_coefficients_pooled.csv\n")
cat("  - missing_summary_for_mi.csv\n")
cat("  - symbol_cleaning_log.csv\n\n")
cat("图形文件:\n")
cat("  - missing_pattern_before_mi.png\n")
cat("  - mi_convergence.png\n")
cat("  - mi_density_check.png\n\n")

sink()

cat("  ✓ 综合报告已保存\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════