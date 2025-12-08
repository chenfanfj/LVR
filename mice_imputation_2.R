############################################################
## 完整优化版: 多重插补 (MI) 
############################################################

## ---------------------------
## 环境设置
## ---------------------------
## 清理环境
rm(list = ls())
gc()

## 强制释放内存
invisible(gc())
memory_limit <- memory.size()
cat("当前内存限制:", memory_limit, "MB\n")

if (memory_limit < 1e4) {
  cat("警告: 内存可能不足！\n")
}

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║   多重插补 (MI) 流程    ║\n")
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
## 步骤 2: 数据清洗（保留原始步骤）
## ═══════════════════════════════════════════════════════

cat("【步骤 2】数据清洗\n")

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
## 2.2 处理含不等号的变量
## ---------------------------
cat("  2.2 处理含不等号的变量\n")

vars_with_symbols <- c(
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "NTproBNP_baseline", "NTproBNP_peak",
  "CK", "CKMB",
  "DD", "APTT", "PT_sec", "TT", "INR", "FIB",
  "FT3", "FT4", "S_TSH",
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
      
      dat[[var]] <- gsub("＜|<|﹤|≤", "", dat[[var]])
      dat[[var]] <- gsub("＞|>|﹥|≥", "", dat[[var]])
      dat[[var]] <- gsub(",", "", dat[[var]])
      dat[[var]] <- gsub("\\s+", "", dat[[var]])
      dat[[var]] <- trimws(dat[[var]])
      dat[[var]] <- suppressWarnings(as.numeric(dat[[var]]))
    }
  }
}

if (nrow(cleaning_log) > 0) {
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

## 结局变量（关键修改：LVEDV_fu 和 ΔLVEDV 不插补）
## 仅 LVESV_fu 和 EF_fu 需要插补（用于计算派生指标）
outcome_vars <- c(
  "LVESV_fu", "EF_fu"
  ## 注意：LVEDV_fu 和 ΔLVEDV 不包含在插补变量中
)

## 暴露变量（金属）
exposure_vars <- c(
  "Cu", "Zn", "Fe", "Se", "Pb",
  "Al", "As", "Cr", "Mn", "Ni", "Mo", 
  "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li"
)

## IPW 模型协变量
ipw_covariates <- c(
  "age", "gender", "resident", "Career",
  "DM", "hypertension", "smoking", "Cancer", "his_stroke", "AF",
  "STEMI", "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",
  "ST_dev", "ST_dep",
  "pPCI", "Stent_no", "Lesion_no",
  "Cardio_shock", "VF", "Cardiac_arrest_in",
  "IN_killip", "OUT_killip", "GRACE_in", "GRACE_in_str",
  "EF_baseline", "LVEDV_baseline",
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "NTproBNP_baseline", "NTproBNP_peak",
  "CK", "CKMB",
  "WBC", "HGB", "PLT", "RBC", "HCT", "NE", "LY", "MO",
  "SCR", "EGFR", "BUN",
  "CRP", "IL_6"
)

## 辅助变量
auxiliary_vars <- c(
  "height", "weight",
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
  
  ## 结局变量（不应插补）
  "LVEDV_fu",  # ← 关键：左室舒张末容积随访值不插补
  
  ## 避免共线性
  "LVESV_baseline",
  
  ## ID
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

## ---------------------------
## 3.4 应用缺失率阈值
## ---------------------------
cat("\n  3.4 应用缺失率阈值\n")

MISSINGNESS_THRESHOLD <- 40  ## 严格阈值: 40%

high_miss_to_remove <- miss_summary %>%
  filter(pct_missing > MISSINGNESS_THRESHOLD, 
         !variable %in% c("has_fu_echo", "LVESV_fu", "EF_fu")) %>%  ## 保护关键结局
  pull(variable)

if (length(high_miss_to_remove) > 0) {
  cat("    移除高缺失变量 (>", MISSINGNESS_THRESHOLD, "%):", 
      length(high_miss_to_remove), "个\n")
  
  ## 保存移除清单
  write.csv(
    data.frame(variable = high_miss_to_remove, 
               reason = paste0("缺失率 >", MISSINGNESS_THRESHOLD, "%")),
    file.path(outputs_dir, "removed_high_missing_vars.csv"),
    row.names = FALSE
  )
  
  dat_for_mi <- dat_for_mi %>% select(-all_of(high_miss_to_remove))
  imputation_vars <- setdiff(imputation_vars, high_miss_to_remove)
}

cat("    保留变量数:", ncol(dat_for_mi) - 1, "个\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 变量质量诊断（完全重写版）
## ═══════════════════════════════════════════════════════

cat("【步骤 4】变量质量诊断\n")

## ---------------------------
## 4.1 检查并移除常数变量
## ---------------------------
cat("  4.1 检查常数变量\n")

vars_to_remove <- character()

for (var in names(dat_for_mi)) {
  if (var %in% c("ID", "has_fu_echo")) next
  
  if (is.numeric(dat_for_mi[[var]])) {
    vals <- na.omit(dat_for_mi[[var]])
    
    if (length(vals) > 0) {
      n_unique <- length(unique(vals))
      sd_val <- sd(vals, na.rm = TRUE)
      
      if (n_unique == 1 || is.na(sd_val) || sd_val < 1e-10) {
        cat("    ✗ 移除:", var, "(唯一值:", n_unique, ", SD:", 
            if (is.na(sd_val)) "NA" else sprintf("%. 2e", sd_val), ")\n")
        vars_to_remove <- c(vars_to_remove, var)
      }
    } else {
      cat("    ✗ 移除:", var, "(无有效值)\n")
      vars_to_remove <- c(vars_to_remove, var)
    }
    
  } else if (is.factor(dat_for_mi[[var]])) {
    if (nlevels(dat_for_mi[[var]]) <= 1) {
      cat("    ✗ 移除因子:", var, "(水平数:", nlevels(dat_for_mi[[var]]), ")\n")
      vars_to_remove <- c(vars_to_remove, var)
    }
  }
}

if (length(vars_to_remove) > 0) {
  dat_for_mi <- dat_for_mi %>% select(-all_of(vars_to_remove))
  imputation_vars <- setdiff(imputation_vars, vars_to_remove)
  cat("    第一轮移除:", length(vars_to_remove), "个变量\n")
} else {
  cat("    ✓ 第一轮未发现常数变量\n")
}

## ---------------------------
## 4.2 安全的共线性检查
## ---------------------------
cat("\n  4.2 检查共线性（安全模式）\n")

## 仅对数值变量进行检查
numeric_cols <- sapply(dat_for_mi, is.numeric)
numeric_vars <- names(dat_for_mi)[numeric_cols]
numeric_vars <- setdiff(numeric_vars, "ID")

cat("    待检查的数值变量:", length(numeric_vars), "个\n")

if (length(numeric_vars) > 1) {
  
  ## 第一步：再次验证每个变量的标准差
  valid_vars <- character()
  
  for (var in numeric_vars) {
    vals <- na.omit(dat_for_mi[[var]])
    
    if (length(vals) >= 2) {
      sd_val <- sd(vals)
      
      if (! is.na(sd_val) && sd_val >= 1e-8) {
        valid_vars <- c(valid_vars, var)
      } else {
        cat("    ⊗ 跳过:", var, "(标准差过小)\n")
      }
    }
  }
  
  cat("    有效变量数:", length(valid_vars), "\n")
  
  ## 第二步：仅对有效变量计算相关性
  if (length(valid_vars) > 1) {
    
    valid_data <- dat_for_mi[, valid_vars, drop = FALSE]
    
    ## 使用 tryCatch 完全捕获错误
    cor_result <- tryCatch({
      cor(valid_data, use = "pairwise.complete.obs")
    }, warning = function(w) {
      cat("    ⚠ 警告:", conditionMessage(w), "\n")
      return(NULL)
    }, error = function(e) {
      cat("    ✗ 错误:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (! is.null(cor_result)) {
      cat("    ✓ 相关性矩阵计算成功 (", nrow(cor_result), "×", ncol(cor_result), ")\n")
      
      ## 查找高相关
      high_cor_vars <- character()
      
      for (i in 1:(nrow(cor_result) - 1)) {
        for (j in (i + 1):ncol(cor_result)) {
          r <- cor_result[i, j]
          
          if (!is.na(r) && abs(r) > 0.995) {
            var1 <- rownames(cor_result)[i]
            var2 <- colnames(cor_result)[j]
            
            cat("      ⚠ 高相关:", var1, "~", var2, "(r =", round(r, 4), ")\n")
            
            ## 移除缺失更多的
            miss1 <- mean(is.na(dat_for_mi[[var1]]))
            miss2 <- mean(is.na(dat_for_mi[[var2]]))
            
            var_to_remove <- if (miss1 > miss2) var1 else var2
            high_cor_vars <- c(high_cor_vars, var_to_remove)
            
            cat("        → 移除:", var_to_remove, "\n")
          }
        }
      }
      
      if (length(high_cor_vars) > 0) {
        high_cor_vars <- unique(high_cor_vars)
        dat_for_mi <- dat_for_mi %>% select(-all_of(high_cor_vars))
        imputation_vars <- setdiff(imputation_vars, high_cor_vars)
        cat("\n    移除高相关变量:", length(high_cor_vars), "个\n")
      } else {
        cat("    ✓ 未发现完美共线性\n")
      }
      
    } else {
      cat("    ⊙ 相关性计算失败，跳过共线性检查\n")
    }
    
  } else {
    cat("    ⊙ 有效变量不足 2 个，跳过相关性检查\n")
  }
  
} else {
  cat("    ⊙ 数值变量不足，跳过共线性检查\n")
}
## 在共线性检查失败后添加
if (is.null(cor_result)) {
  cat("    ⊙ 启用备用清理策略\n")
  
  ## 备用方案: 逐个检查变量是否可用于回归
  problematic_vars <- character()
  
  for (var in numeric_vars) {
    test_data <- dat_for_mi[, c(var, sample(setdiff(numeric_vars, var), min(5, length(numeric_vars)-1)))]
    test_data <- test_data[complete.cases(test_data), ]
    
    if (nrow(test_data) < 20) {
      problematic_vars <- c(problematic_vars, var)
      next
    }
    
    ## 尝试拟合简单线性模型
    tryCatch({
      formula <- as.formula(paste(var, "~ . "))
      lm(formula, data = test_data)
    }, error = function(e) {
      problematic_vars <<- c(problematic_vars, var)
    })
  }
  
  if (length(problematic_vars) > 0) {
    cat("      移除无法建模的变量:", length(problematic_vars), "个\n")
    dat_for_mi <- dat_for_mi %>% select(-all_of(problematic_vars))
    imputation_vars <- setdiff(imputation_vars, problematic_vars)
  }
}
## ---------------------------
## 4.3 生成诊断摘要
## ---------------------------
cat("\n  4.3 诊断摘要\n")

cat("    最终保留变量数:", ncol(dat_for_mi) - 1, "个（排除 ID）\n")
cat("    总移除变量数  :", length(vars_to_remove) + 
      ifelse(exists("high_cor_vars"), length(high_cor_vars), 0), "个\n")

## 保存诊断结果
if (length(vars_to_remove) > 0 || (exists("high_cor_vars") && length(high_cor_vars) > 0)) {
  
  all_removed <- unique(c(vars_to_remove, 
                          if (exists("high_cor_vars")) high_cor_vars else character()))
  
  removed_log <- data.frame(
    variable = all_removed,
    reason = sapply(all_removed, function(v) {
      if (v %in% vars_to_remove) return("常数或零标准差")
      if (exists("high_cor_vars") && v %in% high_cor_vars) return("高度相关")
      return("其他")
    })
  )
  
  write.csv(removed_log, 
            file.path(outputs_dir, "removed_variables_log.csv"),
            row.names = FALSE)
  
  cat("    移除变量清单已保存\n")
}

cat("\n  ✓ 变量质量诊断完成\n\n")

saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_cleaned.rds"))


## ═══════════════════════════════════════════════════════
## 步骤 5: 多重插补配置
## ═══════════════════════════════════════════════════════

cat("【步骤 5】配置多重插补\n")

## ---------------------------
## 5.1 初始化
## ---------------------------
cat("  5.1 初始化 MICE\n")

ini <- mice(dat_for_mi, maxit = 0, print = FALSE)

meth <- ini$method
pred <- ini$predictorMatrix

cat("    ✓ 初始化完成\n")

## ---------------------------
## 5.2 设置不插补的变量（关键修改）
## ---------------------------
cat("\n  5.2 设置不插补的变量\n")

## ID 不插补
meth["ID"] <- ""
pred[, "ID"] <- 0
pred["ID", ] <- 0

## has_fu_echo 不插补（目标变量）
meth["has_fu_echo"] <- ""
pred["has_fu_echo", ] <- 0
pred[, "has_fu_echo"] <- 0
cat("    ✓ has_fu_echo 设置为不插补\n")

## LVEDV_fu 不插补（结局变量，不应插补）
if ("LVEDV_fu" %in% names(meth)) {
  meth["LVEDV_fu"] <- ""
  cat("    ✓ LVEDV_fu 设置为不插补（结局变量）\n")
}

## ΔLVEDV 不插补（派生的结局变量，不应插补）
if ("ΔLVEDV" %in% names(meth)) {
  meth["ΔLVEDV"] <- ""
  cat("    ✓ ΔLVEDV 设置为不插补（派生结局变量）\n")
}

## ---------------------------
## 5.3 变量类型转换 + 插补方法
## ---------------------------
cat("\n  5.3 变量类型转换与插补方法配置\n")

## 二元变量
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
n_converted <- 0

for (var in intersect(binary_vars, names(dat_for_mi))) {
  if (is.numeric(dat_for_mi[[var]])) {
    dat_for_mi[[var]] <- factor(dat_for_mi[[var]], 
                                levels = c(0, 1), 
                                labels = c("No", "Yes"))
    n_converted <- n_converted + 1
  }
  
  ## 检查事件数，决定插补方法
  if (var %in% names(meth) && meth[var] != "") {
    if (is.factor(dat_for_mi[[var]]) && nlevels(dat_for_mi[[var]]) == 2) {
      tab <- table(na.omit(dat_for_mi[[var]]))
      if (min(tab) < 10) {
        ## 事件数太少，改用 PMM
        dat_for_mi[[var]] <- as.numeric(dat_for_mi[[var]]) - 1
        meth[var] <- "pmm"
      } else {
        meth[var] <- "logreg"
      }
    }
  }
}

cat("      转换了", n_converted, "个二元变量\n")

## 有序变量
ordinal_vars <- c("IN_killip", "OUT_killip", "GRACE_in_str")

cat("    转换有序变量:\n")
n_ordinal <- 0

for (var in intersect(ordinal_vars, names(dat_for_mi))) {
  if (is.numeric(dat_for_mi[[var]])) {
    unique_vals <- sort(unique(na.omit(dat_for_mi[[var]])))
    dat_for_mi[[var]] <- factor(dat_for_mi[[var]], 
                                levels = unique_vals, 
                                ordered = TRUE)
    n_ordinal <- n_ordinal + 1
  }
  if (var %in% names(meth) && meth[var] != "") {
    meth[var] <- "polr"
  }
}

cat("      转换了", n_ordinal, "个有序变量\n")

## Lesion_no 特殊处理（关键修改：改用 PMM）
cat("\n    特殊处理 Lesion_no:\n")

if ("Lesion_no" %in% names(dat_for_mi)) {
  cat("      当前类型:", class(dat_for_mi$Lesion_no), "\n")
  
  if (is.factor(dat_for_mi$Lesion_no)) {
    dat_for_mi$Lesion_no <- as.numeric(as.character(dat_for_mi$Lesion_no))
    cat("      → 转换为数值型\n")
  }
  
  meth["Lesion_no"] <- "pmm"
  cat("      ✓ Lesion_no 插补方法设为 PMM\n")
}

## 连续变量用 PMM
continuous_vars <- setdiff(names(meth)[meth != ""], 
                           c(binary_vars, ordinal_vars))

for (var in continuous_vars) {
  if (meth[var] != "") {
    meth[var] <- "pmm"
  }
}

cat("\n    插补方法统计:\n")
cat("      PMM    :", sum(meth == "pmm"), "个\n")
cat("      logreg :", sum(meth == "logreg"), "个\n")
cat("      polr   :", sum(meth == "polr"), "个\n")

## ---------------------------
## 5.4 被动插补
## ---------------------------
cat("\n  5.4 配置被动插补\n")

## BMI
if (all(c("height", "weight") %in% names(dat_for_mi))) {
  dat_for_mi$height <- as.numeric(as.character(dat_for_mi$height))
  dat_for_mi$weight <- as.numeric(as.character(dat_for_mi$weight))
  
  if (!"BMI" %in% names(dat_for_mi)) {
    dat_for_mi$BMI <- NA_real_
  }
  
  meth["BMI"] <- "~I(weight / (height/100)^2)"
  cat("    ✓ BMI 被动插补已配置\n")
}

## 重要说明：ΔLVEDV 不进行被动插补
## 原因：LVEDV_fu 是结局变量，不应插补
## 因此 ΔLVEDV 也不应通过插补计算
cat("    ⊙ ΔLVEDV 不进行被动插补（因 LVEDV_fu 为结局变量）\n")

## 被动变量不预测其他变量
passive_vars <- names(meth)[grepl("^~", meth)]
for (var in passive_vars) {
  if (var %in% rownames(pred)) {
    pred[var, ] <- 0
  }
}

## ---------------------------
## 5.5 重新同步
## ---------------------------
cat("\n  5.5 同步配置\n")

## 因为添加了 BMI，需要重新初始化
ini_final <- mice(dat_for_mi, maxit = 0, print = FALSE)

meth_final <- ini_final$method
pred_final <- ini_final$predictorMatrix

## 恢复所有自定义设置
meth_final["ID"] <- ""
pred_final[, "ID"] <- 0
pred_final["ID", ] <- 0

meth_final["has_fu_echo"] <- ""
pred_final["has_fu_echo", ] <- 0
pred_final[, "has_fu_echo"] <- 0

if ("LVEDV_fu" %in% names(meth_final)) {
  meth_final["LVEDV_fu"] <- ""
}

if ("ΔLVEDV" %in% names(meth_final)) {
  meth_final["ΔLVEDV"] <- ""
}

## 恢复二元变量设置
for (var in intersect(binary_vars, names(meth_final))) {
  if (var %in% names(dat_for_mi) && is.factor(dat_for_mi[[var]])) {
    if (nlevels(dat_for_mi[[var]]) == 2) {
      tab <- table(na.omit(dat_for_mi[[var]]))
      if (min(tab) < 10) {
        meth_final[var] <- "pmm"
      } else {
        meth_final[var] <- "logreg"
      }
    }
  } else if (var %in% names(dat_for_mi) && is.numeric(dat_for_mi[[var]])) {
    meth_final[var] <- "pmm"
  }
}

## 恢复有序变量
for (var in intersect(ordinal_vars, names(meth_final))) {
  if (var %in% names(dat_for_mi) && is.ordered(dat_for_mi[[var]])) {
    meth_final[var] <- "polr"
  }
}

## 恢复 Lesion_no
if ("Lesion_no" %in% names(meth_final)) {
  meth_final["Lesion_no"] <- "pmm"
}

## 恢复 BMI
if ("BMI" %in% names(meth_final)) {
  meth_final["BMI"] <- "~I(weight / (height/100)^2)"
  pred_final["BMI", ] <- 0
}

meth <- meth_final
pred <- pred_final

cat("    同步后: 数据", ncol(dat_for_mi), "列, meth", length(meth), "个\n")

cat("\n  ✓ 配置完成\n\n")

## ---------------------------
## 5.6 优化预测矩阵
## ---------------------------
cat("\n  5.6 优化预测矩阵\n")

## 使用 quickpred 生成稀疏预测矩阵
pred <- quickpred(
  dat_for_mi, 
  mincor = 0.2,      ## 只保留相关性 >0.2 的预测变量
  minpuc = 0.3,      ## 预测变量至少在 30% 样本中可用
  include = c("age", "gender", "EF_baseline", "LVEDV_baseline"),  ## 强制包含核心变量
  exclude = c("ID", "has_fu_echo")
)

## 手动调整关键变量的预测设置
outcome_cols <- c("LVESV_fu", "EF_fu")
for (var in outcome_cols) {
  if (var %in% rownames(pred)) {
    ## 结局变量可被所有变量预测
    pred[var, ] <- ifelse(colnames(pred) %in% c("ID", "has_fu_echo", var), 0, 1)
  }
}

## 验证预测矩阵密度
pred_density <- sum(pred) / (nrow(pred) * ncol(pred))
cat("    预测矩阵密度:", round(pred_density * 100, 1), "% (推荐 <20%)\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 执行多重插补
## ═══════════════════════════════════════════════════════

cat("【步骤 6】执行多重插补\n")

## 最终验证
if (ncol(dat_for_mi) != length(meth)) {
  stop("维度不匹配: data=", ncol(dat_for_mi), ", meth=", length(meth))
}

## 验证关键变量不被插补
cat("  关键变量验证:\n")
cat("    has_fu_echo 方法  :", meth["has_fu_echo"], "(应为空)\n")
if ("LVEDV_fu" %in% names(meth)) {
  cat("    LVEDV_fu 方法     :", meth["LVEDV_fu"], "(应为空)\n")
}
if ("ΔLVEDV" %in% names(meth)) {
  cat("    ΔLVEDV 方法       :", meth["ΔLVEDV"], "(应为空)\n")
}
cat("    Lesion_no 方法    :", meth["Lesion_no"], "(应为 pmm)\n")

## 第一阶段: 快速测试
m_test <- 2
maxit_test <- 3

cat("\n  阶段 1: 小规模测试 (m=", m_test, ", maxit=", maxit_test, ")\n")

set.seed(seed)
imp_test <- mice(
  dat_for_mi,
  m = m_test,
  maxit = maxit_test,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE,   ## 显示进度
  seed = seed
)

## 检查测试结果
if (class(imp_test)[1] != "mids") {
  stop("测试插补失败！请检查数据和方法设置。")
}

cat("\n  ✓ 测试插补成功\n")
cat("    日志事件数:", ifelse(is.null(imp_test$loggedEvents), 0, nrow(imp_test$loggedEvents)), "\n")

## 第二阶段: 完整插补
m <- 20
maxit <- 10  ## 从 15 降低到 10
seed <- 20240101
cat("\n  阶段 2: 完整插补 (m=", m, ", maxit=", maxit, ")\n")

set.seed(seed)
imp <- mice(
  dat_for_mi,
  m = m,
  maxit = maxit,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE,
  seed = seed
)
cat("\n  ✓ 插补完成\n")

## 检查日志事件
if (! is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0) {
  cat("\n  日志事件数:", nrow(imp$loggedEvents), "\n")
  
  if (nrow(imp$loggedEvents) < 100) {
    cat("  ✓ 日志事件数可接受\n")
  } else {
    cat("  ⚠ 日志事件较多\n")
    event_summary <- imp$loggedEvents %>%
      count(out, meth) %>%
      arrange(desc(n))
    cat("\n  前5个问题变量:\n")
    print(head(event_summary, 5), row.names = FALSE)
  }
} else {
  cat("  ✓ 无日志事件\n")
}

## 保存
saveRDS(imp, file.path(outputs_dir, "mice_imputation_final.rds"))
saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_final.rds"))

cat("\n  ✓ 插补对象已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 7: 插补质量诊断
## ═══════════════════════════════════════════════════════

cat("【步骤 7】插补质量诊断\n")

## 7.1 插补完成率
imp_data_1 <- complete(imp, 1)

imputation_diagnostic <- data.frame(
  variable = names(dat_for_mi),
  original_missing = colSums(is.na(dat_for_mi)),
  after_imputation_missing = colSums(is.na(imp_data_1))
) %>%
  filter(original_missing > 0) %>%
  mutate(
    imputed_count = original_missing - after_imputation_missing,
    imputation_rate = round((1 - after_imputation_missing / original_missing) * 100, 1)
  ) %>%
  arrange(imputation_rate, desc(original_missing))

cat("\n  插补效果摘要:\n")
cat("    完全插补变量:", sum(imputation_diagnostic$imputation_rate == 100), "个\n")
cat("    部分插补变量:", sum(imputation_diagnostic$imputation_rate > 0 & 
                         imputation_diagnostic$imputation_rate < 100), "个\n")
cat("    未插补变量  :", sum(imputation_diagnostic$imputation_rate == 0), "个\n")

## 显示未插补的变量
not_imputed <- imputation_diagnostic %>% filter(imputation_rate == 0)
if (nrow(not_imputed) > 0) {
  cat("\n  未插补的变量（应为结局变量）:\n")
  print(not_imputed[, c("variable", "original_missing")], row.names = FALSE)
}

write.csv(imputation_diagnostic, 
          file.path(outputs_dir, "imputation_diagnostic_report.csv"),
          row.names = FALSE)

## 7.2 收敛性图
cat("\n  生成诊断图形...\n")

imputed_vars <- names(meth)[meth != "" & ! grepl("^~", meth) & names(meth) != "ID"]
plot_vars <- head(imputed_vars, 20)

png(file.path(plots_dir, "mi_convergence. png"),
    width = 14, height = 12, res = 300)
plot(imp, plot_vars, layout = c(4, 5))
dev.off()

cat("    ✓ 收敛性图已保存\n")

## 7.3 密度图
density_vars <- character()

for (var in imputed_vars) {
  if (var %in% names(imp_data_1) && is.numeric(imp_data_1[[var]])) {
    vals <- na.omit(imp_data_1[[var]])
    if (length(vals) >= 10 && length(unique(vals)) >= 2) {
      density_vars <- c(density_vars, var)
    }
  }
}

if (length(density_vars) > 0) {
  density_vars_plot <- head(density_vars, 9)
  
  png(file.path(plots_dir, "mi_density.png"),
      width = 12, height = 9, res = 300)
  
  densityplot(imp, 
              as.formula(paste("~", paste(density_vars_plot, collapse = " + "))),
              layout = c(3, 3))
  
  dev.off()
  
  cat("    ✓ 密度图已保存\n")
}

cat("\n  ✓ 质量诊断完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 8: 生成综合报告
## ═══════════════════════════════════════════════════════

cat("【步骤 8】生成综合报告\n")

report_path <- file.path(outputs_dir, "MI_Final_Report.txt")

sink(report_path)

cat("═══════════════════════════════════════════════════════\n")
cat("         多重插补 (MI) 分析最终报告\n")
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

cat("3. 关键设置\n")
cat("─────────────────────────────────────\n")
cat("不插补的变量:\n")
cat("  - has_fu_echo   : 目标变量\n")
cat("  - LVEDV_fu      : 结局变量\n")
cat("  - ΔLVEDV        : 派生结局变量\n\n")

cat("特殊处理:\n")
cat("  - Lesion_no     : 改用 PMM（多分类变量）\n\n")

cat("4. 多重插补\n")
cat("─────────────────────────────────────\n")
cat("插补数据集数      : m =", m, "\n")
cat("最大迭代数        : maxit =", maxit, "\n")
cat("插补方法          :\n")
cat("  - PMM           :", sum(meth == "pmm"), "个\n")
cat("  - logreg        :", sum(meth == "logreg"), "个\n")
cat("  - polr          :", sum(meth == "polr"), "个\n")
cat("  - 被动插补      :", sum(grepl("^~", meth)), "个\n\n")

cat("5. 插补质量\n")
cat("─────────────────────────────────────\n")
cat("完全插补变量      :", sum(imputation_diagnostic$imputation_rate == 100), "个\n")
cat("部分插补变量      :", sum(imputation_diagnostic$imputation_rate > 0 & 
                           imputation_diagnostic$imputation_rate < 100), "个\n")
cat("未插补变量        :", sum(imputation_diagnostic$imputation_rate == 0), "个\n\n")

cat("6. 输出文件\n")
cat("─────────────────────────────────────\n")
cat("插补对象:\n")
cat("  - mice_imputation_final.rds\n")
cat("  - data_for_mi_final.rds\n\n")
cat("诊断文件:\n")
cat("  - imputation_diagnostic_report.csv\n")
cat("  - missing_summary_for_mi.csv\n")
cat("  - symbol_cleaning_log.csv\n\n")
cat("图形文件:\n")
cat("  - missing_pattern_before_mi.png\n")
cat("  - mi_convergence.png\n")
cat("  - mi_density.png\n\n")

cat("═══════════════════════════════════════════════════════\n")

sink()

cat("  ✓ 综合报告已保存\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║            多重插补流程全部完成！                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("关键输出:\n")
cat("  1. mice_imputation_final.rds - 插补对象（用于 IPW）\n")
cat("  2. data_for_mi_final.rds - 插补前数据\n")
cat("  3.  MI_Final_Report.txt - 综合报告\n")
cat("  4. mi_convergence.png - 收敛性诊断\n")
cat("  5.  mi_density.png - 密度图\n\n")

cat("重要说明:\n")
cat("  ✓ has_fu_echo 未被插补（目标变量）\n")
cat("  ✓ LVEDV_fu 未被插补（结局变量）\n")
cat("  ✓ ΔLVEDV 未被插补（派生结局变量）\n")
cat("  ✓ Lesion_no 使用 PMM 插补（多分类变量）\n\n")

cat("下一步:\n")
cat("  使用插补数据构建 IPW 权重并进行结局分析\n\n")

############################################################
## END
############################################################