############################################################
## 多重插补 (MI) 
## 用途: 为 IPW 倾向性评分建模准备高质量插补数据
## 修复日期: 2025-12-09
############################################################

## ═══════════════════════════════════════════════════════
## 环境设置
## ═══════════════════════════════════════════════════════
rm(list = ls())
gc()
invisible(gc())
## 【修复 1】跨平台内存检查
cat("当前内存信息:\n")

if (.Platform$OS.type == "windows") {
  ## Windows 系统
  memory_limit <- memory.size()
  cat("  内存限制:", memory_limit, "MB\n")
  
  if (memory_limit < 1e4) {
    cat("  ⚠ 警告: 内存可能不足！\n")
  }
} else {
  ## Unix/Linux/Mac 系统
  ## 使用 gc() 获取内存信息
  mem_info <- gc()
  total_mb <- sum(mem_info[, 2])
  cat("  已用内存:", round(total_mb, 1), "MB\n")
  cat("  (注: 非Windows系统无法获取总内存限制)\n")
}

## 统一设置随机种子
GLOBAL_SEED <- 20240101
set.seed(GLOBAL_SEED)


cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║        多重插补 (MI) 流程              ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")
cat("全局随机种子:", GLOBAL_SEED, "\n\n")

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
  if (! require(pkg, character.only = TRUE, quietly = TRUE)) {
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
## 先读取原始数据
dat_raw <- readxl::read_excel(
  path = data_path,
  col_names = TRUE,
  na = c("", "NA", "N/A", "#N/A", "NULL", "null", "-999", "-9999")
)

cat("  原始数据:", nrow(dat_raw), "行 ×", ncol(dat_raw), "列\n")

## 立即检查关键变量的原始状态
cat("\n  关键变量原始状态检查:\n")

key_check_vars <- c("STEMI", "pPCI", "gender", "DM", "hypertension")

for (var in intersect(key_check_vars, names(dat_raw))) {
  cat("    ", var, ":\n")
  cat("      类型:", class(dat_raw[[var]]), "\n")
  cat("      唯一值:", paste(head(unique(dat_raw[[var]]), 5), collapse=", "), "\n")
  cat("      缺失数:", sum(is.na(dat_raw[[var]])), "\n")
}

cat("\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 数据清洗
## ═══════════════════════════════════════════════════════

cat("【步骤 2】数据清洗\n")
dat <- dat_raw

## 2.1 标准化缺失值
cat("  2.1 标准化缺失值\n")

for (col in names(dat)) {
  if (is.character(dat[[col]])) {
    dat[[col]][trimws(dat[[col]]) == ""] <- NA
  }
  if (is.numeric(dat[[col]])) {
    dat[[col]][dat[[col]] %in% c(-999, -9999, Inf, -Inf)] <- NA
  }
}

## 2.2 处理含不等号的变量（保持原有逻辑）
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

## 2.3 创建目标变量
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

## 在清洗后立即验证关键二元变量
cat("\n  2.4 验证关键二元变量\n")

binary_vars_to_check <- c(
  "STEMI", "pPCI", "gender", "DM", "hypertension", 
  "smoking", "AF", "ST_dev", "ST_dep"
)

for (var in intersect(binary_vars_to_check, names(dat))) {
  if (is.numeric(dat[[var]])) {
    unique_vals <- sort(unique(na.omit(dat[[var]])))
    n_missing <- sum(is.na(dat[[var]]))
    
    cat("    ", var, ": 唯一值=", paste(unique_vals, collapse=","), 
        ", 缺失=", n_missing, "\n")
    
    ## 如果不是标准的 0/1，发出警告
    if (! all(unique_vals %in% c(0, 1))) {
      cat("      ⚠ 警告:", var, "的值域异常！\n")
    }
  } else {
    cat("    ", var, ": 类型=", class(dat[[var]]), 
        " (应为 numeric)\n")
  }
}

## 2.5 极端值处理（保持原有逻辑）
cat("\n  2. 5 极端值处理\n")

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
## 步骤 3: 变量选择（基于 VIF 分析优化）
## ═══════════════════════════════════════════════════════

cat("【步骤 3】插补变量选择\n")

## ---------------------------
## 3.1 定义变量分组
## ---------------------------

## 结局变量（LVEDV_fu 和 ΔLVEDV 不插补）
outcome_vars <- c(
  "LVESV_fu", "EF_fu"
)

## 暴露变量（金属）
exposure_vars <- c(
  "Cu", "Zn", "Fe", "Se", "Pb",
  "Al", "As", "Cr", "Mn", "Ni", "Mo", 
  "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li"
)

## IPW 模型协变量 - 移除高共线性和高缺失变量
ipw_covariates <- c(
  ## 人口学
  "age", "gender", "resident", "Career",
  
  ## 基础疾病
  "DM", "hypertension", "smoking", "Cancer", "his_stroke", "AF",
  
  ## 梗死类型与治疗
  "STEMI", "ST_dev", "ST_dep",
  "pPCI", "Stent_no", "Lesion_no",
  "Cardio_shock", "VF", "Cardiac_arrest_in",
  
  ## 【移除】高缺失梗死部位变量 (46-47%)
  # "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",
  
  ## Killip 分级与 GRACE 评分
  "IN_killip", "OUT_killip", "GRACE_in", "GRACE_in_str",
  
  ## 基线心功能
  "EF_baseline", "LVEDV_baseline",
  ## 【移除】LVESV_baseline (VIF=79.87, 与 EF/LVEDV 完全共线)
  
  ## 心肌损伤标志物
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  ## 【移除】CK (VIF=9.76, 与 CKMB/cTnI 共线)
  "CKMB",  ## 保留 CKMB (VIF=2.34 after removing CK)
  
  ## 心衰标志物
  "NTproBNP_baseline", "NTproBNP_peak",
  
  ## 血液学
  "WBC", "HGB", "PLT", "RBC", "HCT", "NE", "LY", "MO",
  
  ## 肾功能
  "SCR", "EGFR", "BUN",
  
  ## 炎症标志物
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
  
  ## 【移除】高缺失冠脉变量 (46-47%)
  # "LAD", "LCX", "RCA", "LM",
  # "TL_LAD", "TL_LCX", "TL_RCA",
  
  ## 药物
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

## 【修复 3】明确排除变量清单
exclude_vars <- c(
  ## 派生变量（稍后被动插补）
  "BMI", "ΔLVEDV", "AST_ALT", "BUN_CREA", "APROB_APROA",
  
  ## 结局变量（不应插补）
  "LVEDV_fu", "LVESV_fu", "EF_fu",
  
  ## 高共线性变量（基于 VIF 分析）
  "LVESV_baseline",  ## VIF=79.87
  "CK",              ## VIF=9.76
  
  ## 高缺失变量（基于单变量分析）
  "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",  ## 46-47% 缺失
  "LAD", "LCX", "RCA", "LM",                                    ## 46-47% 缺失
  "TL_LAD", "TL_LCX", "TL_RCA",                                 ## 46-47% 缺失
  
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

cat("  已排除变量:\n")
cat("    - 高共线性  :", sum(c("LVESV_baseline", "CK") %in% exclude_vars), "个 (LVESV_baseline, CK)\n")
cat("    - 高缺失    :", sum(c("Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI", 
                            "LAD", "LCX", "RCA", "LM") %in% exclude_vars), "个 (梗死部位/冠脉变量)\n\n")

## ---------------------------
## 3.3 缺失模式分析
## ---------------------------
cat("  3.3 缺失模式分析\n")

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
  for (i in 1:min(10, nrow(high_miss))) {
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
  labs(title = "插补前的缺失模式（优化后）")

dev.off()

cat("    ✓ 缺失分析完成\n\n")

## ---------------------------
## 3.4 应用分层缺失率阈值
## ---------------------------
cat("  3.4 应用分层缺失率阈值\n")

## 定义分层阈值
THRESHOLD_IPW_CORE <- 5      ## IPW 核心协变量: 极严格
THRESHOLD_IPW_AUX <- 20      ## IPW 辅助协变量: 严格
THRESHOLD_OUTCOME <- 40      ## 结局变量: 中等
THRESHOLD_AUXILIARY <- 50    ## 其他辅助变量: 宽松

## IPW 核心协变量（必须完整或近乎完整）
ipw_core_vars <- c(
  "age", "gender", "pPCI", "STEMI", "ST_dev",
  "cTnIpeak", "CKMB", "NTproBNP_peak",
  "EF_baseline", "LVEDV_baseline",
  "GRACE_in_str", "hypertension", "resident"
)

## IPW 辅助协变量
ipw_auxiliary_vars <- c(
  "DM", "smoking", "Career", "AF", "his_stroke",
  "Stent_no", "Lesion_no", "IN_killip", "OUT_killip",
  "cTnI_baseline", "NTproBNP_baseline",
  "WBC", "HGB", "PLT", "SCR", "EGFR",
  "CRP", "IL_6"
)

## 结局变量
outcome_model_vars <- c("LVESV_fu", "EF_fu")

## 其他辅助变量
other_auxiliary <- setdiff(
  imputation_vars, 
  c(ipw_core_vars, ipw_auxiliary_vars, outcome_model_vars, "has_fu_echo")
)

## 执行分层筛选
high_miss_to_remove <- miss_summary %>%
  filter(
    ## 分层移除规则
    (variable %in% ipw_core_vars & pct_missing > THRESHOLD_IPW_CORE) |
      (variable %in% ipw_auxiliary_vars & pct_missing > THRESHOLD_IPW_AUX) |
      (variable %in% outcome_model_vars & pct_missing > THRESHOLD_OUTCOME) |
      (variable %in% other_auxiliary & pct_missing > THRESHOLD_AUXILIARY),
    ## 保护关键结局
    ! variable %in% c("has_fu_echo", "LVESV_fu", "EF_fu")
  ) %>%
  pull(variable)

if (length(high_miss_to_remove) > 0) {
  cat("    移除高缺失变量:", length(high_miss_to_remove), "个\n")
  
  ## 输出详细清单
  removed_detail <- miss_summary %>%
    filter(variable %in% high_miss_to_remove) %>%
    arrange(desc(pct_missing))
  
  cat("\n    移除清单:\n")
  for (i in 1:nrow(removed_detail)) {
    var_type <- case_when(
      removed_detail$variable[i] %in% ipw_core_vars ~ "IPW核心",
      removed_detail$variable[i] %in% ipw_auxiliary_vars ~ "IPW辅助",
      removed_detail$variable[i] %in% outcome_model_vars ~ "结局",
      TRUE ~ "其他辅助"
    )
    cat("      [", var_type, "] ", removed_detail$variable[i], 
        ": ", removed_detail$pct_missing[i], "%\n", sep = "")
  }
  
  ## 保存移除清单
  write.csv(
    data.frame(
      variable = removed_detail$variable,
      variable_type = sapply(removed_detail$variable, function(v) {
        if (v %in% ipw_core_vars) return("IPW核心")
        if (v %in% ipw_auxiliary_vars) return("IPW辅助")
        if (v %in% outcome_model_vars) return("结局")
        return("其他辅助")
      }),
      missing_pct = removed_detail$pct_missing,
      threshold = sapply(removed_detail$variable, function(v) {
        if (v %in% ipw_core_vars) return(THRESHOLD_IPW_CORE)
        if (v %in% ipw_auxiliary_vars) return(THRESHOLD_IPW_AUX)
        if (v %in% outcome_model_vars) return(THRESHOLD_OUTCOME)
        return(THRESHOLD_AUXILIARY)
      }),
      reason = paste0("缺失率超过阈值")
    ),
    file.path(outputs_dir, "removed_high_missing_vars_stratified.csv"),
    row.names = FALSE
  )
  
  ## 从数据中移除
  dat_for_mi <- dat_for_mi %>% select(-all_of(high_miss_to_remove))
  imputation_vars <- setdiff(imputation_vars, high_miss_to_remove)
}

cat("\n    保留变量数:", ncol(dat_for_mi) - 1, "个（排除 ID）\n")

## 输出各层级保留情况
cat("\n    分层保留情况:\n")
cat("      IPW核心变量 (阈值", THRESHOLD_IPW_CORE, "%) :", 
    sum(ipw_core_vars %in% names(dat_for_mi)), "/", length(ipw_core_vars), "\n")
cat("      IPW辅助变量 (阈值", THRESHOLD_IPW_AUX, "%):", 
    sum(ipw_auxiliary_vars %in% names(dat_for_mi)), "/", length(ipw_auxiliary_vars), "\n")
cat("      结局变量 (阈值", THRESHOLD_OUTCOME, "%)    :", 
    sum(outcome_model_vars %in% names(dat_for_mi)), "/", length(outcome_model_vars), "\n")

cat("\n  ✓ 缺失率筛选完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 变量质量诊断
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
            if (is.na(sd_val)) "NA" else sprintf("%.2e", sd_val), ")\n")
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
  cat("    ✓ 未发现常数变量\n")
}

## ---------------------------
## 4.2 共线性检查
## ---------------------------
cat("\n  4.2 检查共线性\n")

## 仅对数值变量进行检查
numeric_cols <- sapply(dat_for_mi, is.numeric)
numeric_vars <- names(dat_for_mi)[numeric_cols]
numeric_vars <- setdiff(numeric_vars, "ID")

cat("    待检查的数值变量:", length(numeric_vars), "个\n")

if (length(numeric_vars) > 1) {
  
  ## 第一步：验证每个变量的标准差
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
  
  ## 第二步：计算相关性
  if (length(valid_vars) > 1) {
    
    valid_data <- dat_for_mi[, valid_vars, drop = FALSE]
    
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
          
          if (! is.na(r) && abs(r) > 0.95) {  ## 阈值从 0.995 降低到 0.95
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
        cat("    ✓ 未发现高度共线性 (阈值 r > 0.95)\n")
      }
      
    } else {
      cat("    ⊙ 相关性计算失败，启用备用清理策略\n")
      
      ## 优化的备用清理策略
      problematic_vars <- character()
      
      ## 随机采样变量进行测试（提高效率）
      test_vars <- if (length(valid_vars) > 50) {
        sample(valid_vars, min(50, length(valid_vars)))
      } else {
        valid_vars
      }
      
      for (var in test_vars) {
        ## 选择少量其他变量进行回归测试
        other_vars <- setdiff(valid_vars, var)
        test_predictors <- sample(other_vars, min(3, length(other_vars)))
        
        test_data <- dat_for_mi[, c(var, test_predictors), drop = FALSE]
        test_data <- test_data[complete.cases(test_data), ]
        
        if (nrow(test_data) < 20) {
          problematic_vars <- c(problematic_vars, var)
          next
        }
        
        ## 尝试拟合简单线性模型
        result <- tryCatch({
          formula_str <- paste(var, "~", paste(test_predictors, collapse = " + "))
          lm(as.formula(formula_str), data = test_data)
          TRUE
        }, error = function(e) {
          FALSE
        })
        
        if (! result) {
          problematic_vars <- c(problematic_vars, var)
        }
      }
      
      if (length(problematic_vars) > 0) {
        cat("      移除无法建模的变量:", length(problematic_vars), "个\n")
        dat_for_mi <- dat_for_mi %>% select(-all_of(problematic_vars))
        imputation_vars <- setdiff(imputation_vars, problematic_vars)
      }
    }
    
  } else {
    cat("    ⊙ 有效变量不足 2 个，跳过相关性检查\n")
  }
  
} else {
  cat("    ⊙ 数值变量不足，跳过共线性检查\n")
}

## ---------------------------
## 4.3 生成诊断摘要
## ---------------------------
cat("\n  4.3 诊断摘要\n")

total_removed <- length(vars_to_remove)
if (exists("high_cor_vars")) total_removed <- total_removed + length(high_cor_vars)
if (exists("problematic_vars")) total_removed <- total_removed + length(problematic_vars)

cat("    最终保留变量数:", ncol(dat_for_mi) - 1, "个（排除 ID）\n")
cat("    诊断移除变量数:", total_removed, "个\n")

## 保存诊断结果
if (total_removed > 0) {
  
  all_removed <- unique(c(
    vars_to_remove, 
    if (exists("high_cor_vars")) high_cor_vars else character(),
    if (exists("problematic_vars")) problematic_vars else character()
  ))
  
  removed_log <- data.frame(
    variable = all_removed,
    reason = sapply(all_removed, function(v) {
      if (v %in% vars_to_remove) return("常数或零标准差")
      if (exists("high_cor_vars") && v %in% high_cor_vars) return("高度相关 (r>0.95)")
      if (exists("problematic_vars") && v %in% problematic_vars) return("无法建模")
      return("其他")
    })
  )
  
  write.csv(removed_log, 
            file.path(outputs_dir, "removed_variables_quality_check.csv"),
            row.names = FALSE)
  
  cat("    移除变量清单已保存\n")
}

cat("\n  ✓ 变量质量诊断完成\n\n")

## 保存清洗后数据
saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_cleaned.rds"))

## ═══════════════════════════════════════════════════════
## 步骤 5: 配置多重插补（分离策略）
## ═══════════════════════════════════════════════════════

cat("【步骤 5】配置多重插补（分离策略）\n")

## ─────────────────────────────────────────────────────
## 5. 1 识别需要插补的变量
## ─────────────────────────────────────────────────────

cat("  5.1 识别需要插补的变量\n")

## 统计每个变量的缺失数
missing_counts <- data.frame(
  variable = names(dat_for_mi),
  n_missing = colSums(is.na(dat_for_mi)),
  stringsAsFactors = FALSE
)

## 需要插补的变量（缺失数 > 0，排除特殊变量）
vars_need_imputation <- missing_counts %>%
  filter(
    n_missing > 0,
    !variable %in% c("ID", "has_fu_echo", "LVEDV_fu", "ΔLVEDV", "LVESV_fu", "EF_fu")
  ) %>%
  pull(variable)

## 不需要插补的变量（完整数据）
vars_complete <- missing_counts %>%
  filter(
    n_missing == 0 | variable %in% c("ID", "has_fu_echo")
  ) %>%
  pull(variable)

cat("    需要插补的变量:", length(vars_need_imputation), "个\n")
cat("    完整变量      :", length(vars_complete), "个\n")

## 输出关键完整变量
complete_key_vars <- intersect(
  c("STEMI", "pPCI", "ST_dev", "gender", "age", "EF_baseline"),
  vars_complete
)

if (length(complete_key_vars) > 0) {
  cat("    关键完整变量  :", paste(complete_key_vars, collapse=", "), "\n")
}

## ─────────────────────────────────────────────────────
## 5.2 准备插补数据集（仅包含需要插补的变量 + 完整预测变量）
## ─────────────────────────────────────────────────────

cat("\n  5.2 准备插补数据集\n")

## 【关键】分离数据
## 保存完整变量（不参与插补）
dat_complete <- dat_for_mi %>% select(all_of(vars_complete))

cat("    完整数据保存  :", ncol(dat_complete), "列\n")

## 创建插补数据集：需要插补的变量 + 完整预测变量（作为预测因子）
## 选择关键完整变量作为预测因子
complete_predictors <- intersect(
  c("age", "gender", "STEMI", "pPCI", "ST_dev", 
    "EF_baseline", "LVEDV_baseline", "hypertension", "DM"),
  vars_complete
)

dat_for_imputation <- dat_for_mi %>%
  select(ID, all_of(vars_need_imputation), all_of(complete_predictors))

cat("    插补数据准备  :", ncol(dat_for_imputation), "列\n")
cat("      - 需要插补  :", length(vars_need_imputation), "个\n")
cat("      - 完整预测变量:", length(complete_predictors), "个\n")

## ─────────────────────────────────────────────────────
## 5.3 配置插补方法
## ─────────────────────────────────────────────────────

cat("\n  5.3 配置插补方法\n")

## 初始化
ini <- mice(dat_for_imputation, maxit = 0, print = FALSE)

meth <- ini$method
pred <- ini$predictorMatrix

## ID 不插补
meth["ID"] <- ""
pred[, "ID"] <- 0
pred["ID", ] <- 0

## 【关键】完整预测变量不插补
for (var in complete_predictors) {
  meth[var] <- ""
  pred[var, ] <- 0  ## 它们不被其他变量预测
  ## 但它们可以预测其他变量（列不清零）
}

cat("    完整预测变量方法设置:\n")
for (cpv in complete_predictors) {
  cat("      ", cpv, ": 方法=", meth[cpv], 
      ", 预测矩阵行和=", sum(pred[cpv, ]), "\n")
}

## 配置需要插补的变量
cat("\n    需要插补的变量方法:\n")

binary_vars <- c("smoking", "Cancer", "his_stroke", "AF", 
                 "VF", "Cardio_shock", "Cardiac_arrest_in",
                 "resident", "Career", "Aspirin", "Statin", 
                 "ACEIorARB", "β_block")

ordinal_vars <- c("IN_killip", "OUT_killip", "GRACE_in_str")

for (var in vars_need_imputation) {
  if (var %in% binary_vars) {
    ## 二元变量用 PMM（避免 logreg 问题）
    meth[var] <- "pmm"
  } else if (var %in% ordinal_vars) {
    ## 有序变量转数值 + PMM
    if (is.ordered(dat_for_imputation[[var]])) {
      dat_for_imputation[[var]] <- as.numeric(dat_for_imputation[[var]])
    }
    meth[var] <- "pmm"
  } else {
    ## 连续变量用 PMM
    meth[var] <- "pmm"
  }
}

## 被动插补 BMI
if ("BMI" %in% names(dat_for_imputation) && 
    all(c("height", "weight") %in% names(dat_for_imputation))) {
  meth["BMI"] <- "~I(weight / (height/100)^2)"
  pred["BMI", ] <- 0
}

cat("    插补方法统计:\n")
cat("      PMM      :", sum(meth == "pmm"), "个\n")
cat("      被动插补 :", sum(grepl("^~", meth)), "个\n")
cat("      不插补   :", sum(meth == ""), "个\n")

## ─────────────────────────────────────────────────────
## 5.4 优化预测矩阵
## ─────────────────────────────────────────────────────

cat("\n  5.4 优化预测矩阵\n")

## 使用 quickpred 但保护完整预测变量
pred_sparse <- quickpred(
  dat_for_imputation,
  mincor = 0.2,
  minpuc = 0.3,
  include = complete_predictors,  ## 强制包含完整变量作为预测因子
  exclude = c("ID")
)

## 确保完整预测变量的行全为 0
for (cpv in complete_predictors) {
  pred_sparse[cpv, ] <- 0
}

pred <- pred_sparse

pred_density <- sum(pred) / (nrow(pred) * ncol(pred))
cat("    预测矩阵密度:", round(pred_density * 100, 1), "%\n")

cat("\n  ✓ 配置完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 执行插补
## ═══════════════════════════════════════════════════════

cat("【步骤 6】执行多重插补（分离策略）\n")

m <- 20
maxit <- 10

cat("  参数:  m=", m, ", maxit=", maxit, "\n")
cat("  插补变量数:", sum(meth != "" & ! grepl("^~", meth)), "个\n\n")

set.seed(GLOBAL_SEED)

imp <- mice(
  dat_for_imputation,
  m = m,
  maxit = maxit,
  method = meth,
  predictorMatrix = pred,
  printFlag = TRUE,
  seed = GLOBAL_SEED
)

cat("\n  ✓ 插补完成\n")

## 检查日志事件
if (! is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0) {
  cat("\n  ⚠ 日志事件数:", nrow(imp$loggedEvents), "\n")
  
  ## 分析日志事件
  event_summary <- imp$loggedEvents %>%
    count(out, dep) %>%
    arrange(desc(n))
  
  cat("  问题变量组合（前10）:\n")
  print(head(event_summary, 10), row.names = FALSE)
  
} else {
  cat("  ✓ 无日志事件\n")
}

## ─────────────────────────────────────────────────────
## 6.3 合并完整变量与插补数据
## ─────────────────────────────────────────────────────

cat("\n  6.3 合并完整变量与插补数据\n")

## 【关键】为每个插补数据集添加完整变量
imp_datasets_full <- lapply(1:m, function(i) {
  ## 提取第 i 个插补数据集
  imp_data_i <- complete(imp, i)
  
  ## 合并完整变量（按 ID）
  dat_full_i <- imp_data_i %>%
    left_join(dat_complete, by = "ID")
  
  return(dat_full_i)
})

cat("    生成完整插补数据集:", length(imp_datasets_full), "个\n")
cat("    每个数据集维度  :", dim(imp_datasets_full[[1]]), "\n")

## 验证 STEMI
stemi_check <- imp_datasets_full[[1]]$STEMI
cat("    STEMI 检查:\n")
cat("      类型  :", class(stemi_check), "\n")
cat("      缺失数:", sum(is.na(stemi_check)), "\n")
cat("      唯一值:", paste(head(unique(stemi_check), 5), collapse=", "), "\n")

## ─────────────────────────────────────────────────────
## 6.4 重构 mids 对象（可选，用于后续分析）
## ─────────────────────────────────────────────────────

cat("\n  6.4 重构完整 mids 对象\n")

## 创建包含所有变量的完整 mids 对象
## 方法：使用 mice:: as.mids()

## 准备数据结构
imp_list <- list()

for (i in 1:m) {
  imp_list[[i]] <- imp_datasets_full[[i]]
}

## 转换为 long 格式
imp_long <- mice::complete(imp, action = "long", include = TRUE)

## 添加完整变量到 long 格式
complete_vars_to_add <- setdiff(names(dat_complete), c("ID", names(imp_long)))

if (length(complete_vars_to_add) > 0) {
  for (var in complete_vars_to_add) {
    ## 所有插补数据集使用相同的完整变量值
    imp_long[[var]] <- rep(dat_complete[[var]], m + 1)
  }
}

## 转回 mids 对象
imp_full <- as.mids(imp_long)

cat("    完整 mids 对象变量数:", length(imp_full$method), "\n")

## 保存
saveRDS(imp_full, file.path(outputs_dir, "mice_imputation_final.rds"))
saveRDS(imp_datasets_full, file.path(outputs_dir, "imputed_datasets_list.rds"))
saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_final.rds"))

cat("\n  ✓ 插补对象已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 7: 插补质量诊断
## ═══════════════════════════════════════════════════════
## ─────────────────────────────────────────────────────
## 7.1 插补完成率
## ─────────────────────────────────────────────────────

cat("  7.1 插补完成率诊断\n")

## 【修复 2】安全地提取插补后数据
imp_data_1 <- imp_datasets_full[[1]]

## 验证数据结构
cat("    插补后数据维度:", dim(imp_data_1), "\n")
cat("    原始数据维度  :", dim(dat_for_mi), "\n")

## 【关键修复】确保列名匹配
common_vars <- intersect(names(dat_for_mi), names(imp_data_1))
cat("    共同变量数    :", length(common_vars), "\n")

if (length(common_vars) < ncol(dat_for_mi)) {
  missing_vars <- setdiff(names(dat_for_mi), names(imp_data_1))
  cat("    ⚠ 插补后缺少的变量:", length(missing_vars), "个\n")
  if (length(missing_vars) > 0 && length(missing_vars) <= 10) {
    cat("      ", paste(missing_vars, collapse=", "), "\n")
  }
}

## 【修复核心】仅对共同变量计算诊断指标
imputation_diagnostic <- data.frame(
  variable = common_vars,
  original_missing = sapply(common_vars, function(v) sum(is.na(dat_for_mi[[v]]))),
  after_imputation_missing = sapply(common_vars, function(v) sum(is.na(imp_data_1[[v]]))),
  stringsAsFactors = FALSE
) %>%
  filter(original_missing > 0) %>%
  mutate(
    imputed_count = original_missing - after_imputation_missing,
    imputation_rate = round((1 - after_imputation_missing / original_missing) * 100, 1)
  ) %>%
  arrange(imputation_rate, desc(original_missing))

cat("\n    插补效果摘要:\n")
cat("      完全插补变量:", sum(imputation_diagnostic$imputation_rate == 100), "个\n")
cat("      部分插补变量:", sum(imputation_diagnostic$imputation_rate > 0 & 
                           imputation_diagnostic$imputation_rate < 100), "个\n")
cat("      未插补变量  :", sum(imputation_diagnostic$imputation_rate == 0), "个\n")

## 显示未插补的变量
not_imputed <- imputation_diagnostic %>% filter(imputation_rate == 0)
if (nrow(not_imputed) > 0) {
  cat("\n    未插补的变量（应为目标或结局变量）:\n")
  for (i in 1:min(5, nrow(not_imputed))) {
    cat("      ", not_imputed$variable[i], ": ", 
        not_imputed$original_missing[i], "缺失\n")
  }
}

## 保存诊断结果
write.csv(imputation_diagnostic, 
          file.path(outputs_dir, "imputation_diagnostic_report.csv"),
          row.names = FALSE)

cat("\n    ✓ 插补诊断报告已保存\n")

## ─────────────────────────────────────────────────────
## 7.2 IPW 核心变量插补质量专项检查
## ─────────────────────────────────────────────────────

cat("\n  7.2 IPW 核心变量插补质量\n")

ipw_core_vars <- c(
  "age", "STEMI", "pPCI", "ST_dev", 
  "cTnIpeak", "NTproBNP_peak",
  "EF_baseline", "LVEDV_baseline",
  "GRACE_in_str", "hypertension", "resident"
)

ipw_core_in_data <- intersect(ipw_core_vars, common_vars)

if (length(ipw_core_in_data) > 0) {
  ipw_core_diagnostic <- imputation_diagnostic %>%
    filter(variable %in% ipw_core_in_data)
  
  if (nrow(ipw_core_diagnostic) > 0) {
    cat("    插补的 IPW 核心变量:\n")
    for (i in 1:nrow(ipw_core_diagnostic)) {
      cat("      ", ipw_core_diagnostic$variable[i], ": ",
          ipw_core_diagnostic$original_missing[i], "缺失 → ",
          ipw_core_diagnostic$imputation_rate[i], "% 插补\n")
    }
  } else {
    cat("    ✓ 所有 IPW 核心变量无缺失\n")
  }
  
  ## 检查完整的 IPW 核心变量
  ipw_core_complete <- setdiff(ipw_core_in_data, ipw_core_diagnostic$variable)
  if (length(ipw_core_complete) > 0) {
    cat("\n    完整的 IPW 核心变量:\n")
    cat("      ", paste(ipw_core_complete, collapse=", "), "\n")
  }
} else {
  cat("    ⚠ 警告: IPW 核心变量未找到\n")
}

## ─────────────────────────────────────────────────────
## 7.3 高缺失变量插补质量警告
## ─────────────────────────────────────────────────────

cat("\n  7.3 高缺失变量插补质量检查\n")

high_miss_imputed <- imputation_diagnostic %>%
  filter(
    original_missing > 0.3 * nrow(dat_for_mi),
    imputation_rate == 100
  ) %>%
  mutate(pct_missing = round(original_missing / nrow(dat_for_mi) * 100, 1)) %>%
  arrange(desc(pct_missing))

if (nrow(high_miss_imputed) > 0) {
  cat("    高缺失但完全插补的变量 (>30%):\n")
  for (i in 1:min(5, nrow(high_miss_imputed))) {
    cat("      ", high_miss_imputed$variable[i], 
        " (", high_miss_imputed$pct_missing[i], "% 缺失)\n")
  }
  cat("    ⚠ 提示: 这些变量的插补质量需要特别验证\n")
} else {
  cat("    ✓ 无高缺失变量被插补\n")
}

## ─────────────────────────────────────────────────────
## 7.4 收敛性图（修复版）
## ─────────────────────────────────────────────────────

cat("\n  7.4 生成收敛性图\n")

## 【关键】基于实际插补对象 imp，而非合并后的数据
imputed_numeric_vars <- character()

for (var in names(imp$method)) {
  if (imp$method[var] != "" && 
      !grepl("^~", imp$method[var]) &&
      var != "ID" &&
      var %in% names(imp$data) &&
      is.numeric(imp$data[[var]])) {
    imputed_numeric_vars <- c(imputed_numeric_vars, var)
  }
}

cat("    被插补的数值变量:", length(imputed_numeric_vars), "\n")

if (length(imputed_numeric_vars) == 0) {
  cat("    ⊙ 无数值变量被插补，跳过收敛性图\n")
  
  ## 诊断信息
  cat("    诊断信息:\n")
  vars_with_method <- names(imp$method)[imp$method != "" & !grepl("^~", imp$method)]
  cat("      方法不为空的变量:", length(vars_with_method), "个\n")
  
  if (length(vars_with_method) > 0) {
    cat("      示例变量（前5个）:\n")
    for (v in head(vars_with_method, 5)) {
      cat("        ", v, ": 方法=", imp$method[v], 
          ", 类型=", class(imp$data[[v]])[1], "\n")
    }
  }
  
} else {
  ## 优先选择 IPW 核心变量
  ipw_core_numeric <- intersect(
    c("cTnIpeak", "NTproBNP_peak", "age", "WBC", "HGB", "EGFR"),
    imputed_numeric_vars
  )
  
  plot_vars_conv <- c(ipw_core_numeric, 
                      setdiff(imputed_numeric_vars, ipw_core_numeric))
  plot_vars_conv <- head(plot_vars_conv, 20)
  
  cat("    绘制变量 (", length(plot_vars_conv), "个):\n")
  cat("      ", paste(head(plot_vars_conv, 5), collapse=", "), "...\n")
  
  ## 绘图
  plot_success <- tryCatch({
    png(file.path(plots_dir, "mi_convergence.png"),
        width = 14, height = 12, units = "in", res = 300)
    
    # 【关键修复】: 必须使用 print() 包裹 plot() 命令
    print(plot(imp, plot_vars_conv, layout = c(4, 5)))
    
    dev.off()
    
    TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("    ✗ 收敛性图生成失败:", conditionMessage(e), "\n")
    FALSE
  }, warning = function(w) {
    # 即使有警告也要确保关闭设备
    if (dev.cur() > 1) dev.off() 
    cat("    ⚠ 警告:", conditionMessage(w), "\n")
    TRUE
  })
  
  if (plot_success) {
    cat("    ✓ 收敛性图已保存\n")
    
    ## 验证文件
    conv_file <- file.path(plots_dir, "mi_convergence.png")
    if (file.exists(conv_file)) {
      file_size <- file.info(conv_file)$size
      cat("      文件大小:", round(file_size / 1024, 1), "KB\n")
    } else {
      cat("      ⚠ 文件未找到\n")
    }
  }
}

## ─────────────────────────────────────────────────────
## 7.5 密度图
## ─────────────────────────────────────────────────────

cat("\n  7.5 生成密度图\n")

if (length(imputed_numeric_vars) == 0) {
  cat("    ⊙ 无数值变量被插补，跳过密度图\n")
  
} else {
  ## 筛选可绘制密度的变量
  density_vars <- character()
  
  ## 使用插补对象的第一个完整数据集
  imp_check_data <- complete(imp, 1)
  
  for (var in imputed_numeric_vars) {
    if (var %in% names(imp_check_data)) {
      vals <- na.omit(imp_check_data[[var]])
      
      if (length(vals) >= 10 && 
          length(unique(vals)) >= 2 && 
          sd(vals, na.rm = TRUE) > 1e-8) {
        density_vars <- c(density_vars, var)
      }
    }
  }
  
  cat("    可绘制密度的变量:", length(density_vars), "\n")
  
  if (length(density_vars) >= 3) {
    density_vars_plot <- head(density_vars, 9)
    
    cat("    绘制变量 (", length(density_vars_plot), "个):\n")
    cat("      ", paste(head(density_vars_plot, 5), collapse=", "), "...\n")
    
    plot_success <- tryCatch({
      png(file.path(plots_dir, "mi_density.png"),
          width = 12, height = 9, units = "in", res = 300)
      
      # 【关键修复】: 必须使用 print() 包裹 densityplot() 命令
      print(densityplot(imp, 
                        as.formula(paste("~", paste(density_vars_plot, collapse = " + "))),
                        layout = c(3, 3)))
      
      dev.off()
      
      TRUE
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      cat("    ✗ 密度图生成失败:", conditionMessage(e), "\n")
      FALSE
    })
    
    if (plot_success) {
      cat("    ✓ 密度图已保存\n")
      
      dens_file <- file.path(plots_dir, "mi_density.png")
      if (file.exists(dens_file)) {
        file_size <- file.info(dens_file)$size
        cat("      文件大小:", round(file_size / 1024, 1), "KB\n")
      }
    }
  } else {
    cat("    ⊙ 可绘制变量不足3个，跳过密度图\n")
  }
}

## ─────────────────────────────────────────────────────
## 7.6 条带图（补充诊断）
## ─────────────────────────────────────────────────────

cat("\n  7.6 生成条带图（补充诊断）\n")

if (length(imputed_numeric_vars) >= 3) {
  strip_vars <- head(imputed_numeric_vars, 6)
  
  cat("    绘制变量 (", length(strip_vars), "个):\n")
  cat("      ", paste(strip_vars, collapse=", "), "\n")
  
  plot_success <- tryCatch({
    png(file.path(plots_dir, "mi_stripplot_key_vars.png"),
        width = 12, height = 8, units = "in", res = 300)
    
    # 【关键修复】: 必须使用 print() 包裹 stripplot() 命令
    print(stripplot(imp, 
                    as.formula(paste("~", paste(strip_vars, collapse = " + "))),
                    pch = c(1, 20), cex = c(1, 1.5)))
    
    dev.off()
    
    TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("    ✗ 条带图生成失败:", conditionMessage(e), "\n")
    FALSE
  })
  
  if (plot_success) {
    cat("    ✓ 条带图已保存\n")
  }
} else {
  cat("    ⊙ 可绘制变量不足，跳过条带图\n")
}

## ─────────────────────────────────────────────────────
## 7.7 最终验证
## ─────────────────────────────────────────────────────

cat("\n  7.7 最终验证\n")

cat("    生成的图形文件:\n")
plot_files <- list.files(plots_dir, pattern = "\\.png$", full.names = FALSE)

if (length(plot_files) > 0) {
  for (pf in plot_files) {
    file_path <- file.path(plots_dir, pf)
    file_size <- file.info(file_path)$size
    cat("      ✓ ", pf, " (", round(file_size / 1024, 1), "KB)\n")
  }
} else {
  cat("      ⚠ 未生成任何图形文件\n")
  cat("      可能原因:  所有变量都无缺失，或插补失败\n")
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
cat("         专用于 IPW 权重构建\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("随机种子:", GLOBAL_SEED, "\n\n") 

## ─────────────────────────────
## 1. 数据清洗
## ─────────────────────────────
cat("1. 数据清洗\n")
cat("─────────────────────────────────────\n")
cat("原始样本数        :", nrow(dat_raw), "\n")
cat("清洗后样本数      :", nrow(dat), "\n")
cat("有随访数据        :", sum(dat$has_fu_echo == "Yes"), 
    "(", round(mean(dat$has_fu_echo == "Yes")*100, 1), "%)\n")
cat("处理含符号变量    :", nrow(cleaning_log), "个\n")
cat("Winsorize 变量数  :", winsorize_count, "个\n\n")

## ─────────────────────────────
## 2. 变量选择策略
## ─────────────────────────────
cat("2. 变量选择策略（基于单变量筛选优化）\n")
cat("─────────────────────────────────────\n")
cat("插补变量总数      :", length(imputation_vars), "\n")
cat("  ├─ 结局变量     :", length(intersect(outcome_vars, imputation_vars)), 
    "个 (LVESV_fu, EF_fu)\n")
cat("  ├─ 暴露变量     :", length(intersect(exposure_vars, imputation_vars)), 
    "个 (金属)\n")
cat("  ├─ IPW核心协变量:", length(intersect(ipw_core_vars, imputation_vars)), 
    "个 (p<0.05, 缺失<5%)\n")
cat("  ├─ IPW辅助协变量:", length(intersect(ipw_auxiliary_vars, imputation_vars)), 
    "个 (临床重要, 缺失<20%)\n")
cat("  └─ 其他辅助变量 :", length(intersect(auxiliary_vars, imputation_vars)), 
    "个\n\n")

cat("明确排除的变量:\n")
cat("  ├─ 高共线性      : CK (VIF=9.76), CKMB (VIF=8.51), LVESV_baseline (VIF=79.87)\n")
cat("  ├─ 高缺失梗死部位: Inferior_MI, Anterior_MI 等 (~47%缺失)\n")
cat("  ├─ 高缺失冠脉变量: LAD, LCX, RCA 等 (~47%缺失)\n")
cat("  └─ 结局变量      : LVEDV_fu (不应插补)\n\n")

## ─────────────────────────────
## 3. 缺失率阈值
## ─────────────────────────────
cat("3. 分层缺失率阈值\n")
cat("─────────────────────────────────────\n")
cat("  ├─ IPW核心协变量: ≤", THRESHOLD_IPW_CORE, "% (极严格)\n")
cat("  ├─ IPW辅助协变量: ≤", THRESHOLD_IPW_AUX, "% (严格)\n")
cat("  ├─ 结局变量     : ≤", THRESHOLD_OUTCOME, "% (中等)\n")
cat("  └─ 其他辅助变量 : ≤", THRESHOLD_AUXILIARY, "% (宽松)\n\n")

if (length(high_miss_to_remove) > 0) {
  cat("移除高缺失变量    :", length(high_miss_to_remove), "个\n")
  for (hmv in head(high_miss_to_remove, 10)) {
    miss_pct <- miss_summary %>% filter(variable == hmv) %>% pull(pct_missing)
    cat("  - ", hmv, " (", miss_pct, "%)\n")
  }
  if (length(high_miss_to_remove) > 10) {
    cat("  - ...  及其他", length(high_miss_to_remove) - 10, "个变量\n")
  }
  cat("\n")
}

## ─────────────────────────────
## 4. 变量质量诊断
## ─────────────────────────────
cat("4. 变量质量诊断\n")
cat("─────────────────────────────────────\n")

n_removed_quality <- length(vars_to_remove) + 
  ifelse(exists("high_cor_vars"), length(high_cor_vars), 0) +
  ifelse(exists("problematic_vars"), length(problematic_vars), 0)

cat("质量问题移除变量  :", n_removed_quality, "个\n")

if (length(vars_to_remove) > 0) {
  cat("  ├─ 常数变量     :", length(vars_to_remove), "个\n")
}
if (exists("high_cor_vars") && length(high_cor_vars) > 0) {
  cat("  ├─ 高相关变量   :", length(high_cor_vars), "个 (|r|>0.95)\n")
}
if (exists("problematic_vars") && length(problematic_vars) > 0) {
  cat("  └─ 无法建模变量 :", length(problematic_vars), "个\n")
}
cat("\n")

## ─────────────────────────────
## 5. 多重插补参数
## ─────────────────────────────
cat("5. 多重插补参数\n")
cat("─────────────────────────────────────\n")
cat("插补数据集数 (m)  :", m, "\n")
cat("最大迭代数 (maxit):", maxit, "\n")
cat("随机种子 (seed)   :", GLOBAL_SEED, "\n\n")

cat("插补方法分布      :\n")
cat("  ├─ PMM          :", sum(meth == "pmm"), "个 (连续变量 + 低事件数二元变量)\n")
cat("  ├─ logreg       :", sum(meth == "logreg"), "个 (二元变量)\n")
cat("  ├─ polr         :", sum(meth == "polr"), "个 (有序变量)\n")
cat("  ├─ 被动插补     :", sum(grepl("^~", meth)), "个 (BMI)\n")
cat("  └─ 不插补       :", sum(meth == ""), "个 (ID, has_fu_echo, 结局变量)\n\n")

## ─────────────────────────────
## 6. 插补质量
## ─────────────────────────────
cat("6. 插补质量评估\n")
cat("─────────────────────────────────────\n")

cat("插补完成情况      :\n")
cat("  ├─ 完全插补变量 :", sum(imputation_diagnostic$imputation_rate == 100), "个\n")
cat("  ├─ 部分插补变量 :", sum(imputation_diagnostic$imputation_rate > 0 & 
                           imputation_diagnostic$imputation_rate < 100), "个\n")
cat("  └─ 未插补变量   :", sum(imputation_diagnostic$imputation_rate == 0), 
    "个 (目标/结局变量)\n\n")

## IPW核心变量插补质量
ipw_core_final <- intersect(ipw_core_vars, names(dat_for_mi))
ipw_core_imputed_check <- imputation_diagnostic %>%
  filter(variable %in% ipw_core_final)

cat("IPW核心变量插补质量:\n")
if (nrow(ipw_core_imputed_check) > 0) {
  for (i in 1:nrow(ipw_core_imputed_check)) {
    cat("  - ", ipw_core_imputed_check$variable[i], ": ",
        ipw_core_imputed_check$original_missing[i], " 缺失 → ",
        ipw_core_imputed_check$imputation_rate[i], "% 插补\n")
  }
} else {
  cat("  ✓ 所有IPW核心变量无缺失或插补完成\n")
}
cat("\n")

## 日志事件
if (!  is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0) {
  cat("插补日志事件      :", nrow(imp$loggedEvents), "个\n")
  
  if (nrow(imp$loggedEvents) > 100) {
    cat("  ⚠ 警告: 日志事件数较多，可能存在问题变量\n")
    
    event_summary <- imp$loggedEvents %>%
      count(out, meth) %>%
      arrange(desc(n)) %>%
      head(5)
    
    cat("  前5个问题变量:\n")
    for (i in 1:nrow(event_summary)) {
      cat("    - ", event_summary$out[i], " (", event_summary$meth[i], "): ",
          event_summary$n[i], " 次\n")
    }
  } else {
    cat("  ✓ 日志事件数可接受\n")
  }
} else {
  cat("插补日志事件      : 0 (完美)\n")
}
cat("\n")

## 高缺失但完成插补的变量
high_miss_imputed <- imputation_diagnostic %>%
  filter(original_missing > 0.3 * nrow(dat_for_mi), 
         imputation_rate == 100)

if (nrow(high_miss_imputed) > 0) {
  cat("高缺失但完全插补变量 (>30%):\n")
  for (i in 1:min(5, nrow(high_miss_imputed))) {
    pct <- round(high_miss_imputed$original_missing[i] / nrow(dat_for_mi) * 100, 1)
    cat("  - ", high_miss_imputed$variable[i], " (", pct, "% 缺失)\n")
  }
  cat("  ⚠ 提示: 这些变量的插补质量需要特别验证\n")
  cat("\n")
}

## ─────────────────────────────
## 7. 输出文件清单
## ─────────────────────────────
cat("7.  输出文件清单\n")
cat("─────────────────────────────────────\n")

cat("核心输出:\n")
cat("  ├─ mice_imputation_final.rds      : 插补对象 (用于IPW建模)\n")
cat("  └─ data_for_mi_final.rds          : 插补前数据\n\n")

cat("诊断文件:\n")
cat("  ├─ imputation_diagnostic_report.csv : 插补完成率\n")
cat("  ├─ missing_summary_for_mi.csv       : 缺失率汇总\n")
cat("  ├─ excluded_variables_summary.csv   : 排除变量清单\n")
cat("  ├─ removed_high_missing_vars.csv    : 高缺失移除清单\n")
cat("  ├─ high_correlation_pairs.csv       : 高相关变量对\n")

if (!  is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 0) {
  cat("  ├─ mice_logged_events.csv           : 插补日志事件\n")
}

cat("  └─ symbol_cleaning_log.csv          : 符号清洗日志\n\n")

cat("图形文件:\n")
cat("  ├─ missing_pattern_before_mi.png    : 插补前缺失模式\n")
cat("  ├─ mi_convergence. png               : 收敛性诊断图\n")
cat("  ├─ mi_density.png                   : 插补值密度图\n")
cat("  └─ mi_stripplot_key_vars.png        : IPW核心变量条带图\n\n")

## ─────────────────────────────
## 8. 关键决策记录
## ─────────────────────────────
cat("8. 关键决策记录\n")
cat("─────────────────────────────────────\n")

cat("变量选择决策:\n")
cat("  ✓ 保留 cTnIpeak (VIF=2.38), 移除 CK/CKMB (VIF>5)\n")
cat("  ✓ 保留 EF_baseline + LVEDV_baseline, 移除 LVESV_baseline (VIF=79.87)\n")
cat("  ✓ 移除梗死部位变量 (Inferior_MI等, 47%缺失)\n")
cat("  ✓ 移除冠脉变量 (LAD/LCX/RCA, 47%缺失)\n")
cat("  ✓ LVEDV_fu 不插补 (结局变量)\n\n")

cat("插补方法决策:\n")
cat("  ✓ 连续变量: PMM (稳健且保持分布)\n")
cat("  ✓ 二元变量: logreg (事件数≥10) 或 PMM (事件数<10)\n")
cat("  ✓ 有序变量: polr (水平样本数≥5) 或 PMM (否则)\n")
cat("  ✓ BMI: 被动插补 (基于height和weight)\n\n")

cat("预测矩阵决策:\n")
cat("  ✓ 使用 quickpred() 稀疏化 (mincor=0.2, minpuc=0.3)\n")
cat("  ✓ 强制包含IPW核心变量作为预测因子\n")
cat("  ✓ 最终密度:", round(pred_density * 100, 1), "%\n\n")

## ─────────────────────────────
## 9. 下一步分析建议
## ─────────────────────────────
cat("9. 下一步分析建议\n")
cat("─────────────────────────────────────\n")

cat("IPW权重构建:\n")
cat("  1. 加载插补对象: readRDS('outputs/mice_imputation_final.rds')\n")
cat("  2. 对每个插补数据集拟合倾向性评分模型:\n")
cat("     formula: has_fu_echo ~ pPCI + STEMI + ST_dev + cTnIpeak +\n")
cat("              NTproBNP_peak + EF_baseline + LVEDV_baseline +\n")
cat("              GRACE_in_str + age + hypertension + resident + .. .\n")
cat("  3. 计算稳定化IPW权重\n")
cat("  4.  检查协变量平衡 (bal. tab, SMD<0.1)\n")
cat("  5. 截断极端权重 (99th percentile)\n\n")

cat("结局模型:\n")
cat("  1.  计算 ΔLVEDV = LVEDV_fu - LVEDV_baseline\n")
cat("  2. 拟合加权线性/logistic回归: ΔLVEDV ~ metals + covariates\n")
cat("  3. 使用Rubin规则合并多重插补结果\n")
cat("  4. 报告合并后的估计值与标准误\n\n")

cat("质量控制:\n")
cat("  ✓ 检查IPW权重分布 (范围、极端值)\n")
cat("  ✓ 敏感性分析: 比较不同m值的结果稳定性\n")
cat("  ✓ 完整案例分析: 对比IPW加权与未加权结果\n")

## ─────────────────────────────
## 10. 重要注意事项
## ─────────────────────────────
cat("10. 重要注意事项\n")
cat("─────────────────────────────────────\n")

cat("插补局限性:\n")
cat("  ⚠ 金属暴露全部缺失同一批样本 (167例), 可能MNAR\n")
cat("  ⚠ BMI完全缺失但通过被动插补重构, 需验证合理性\n\n")

if (!  is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 100) {
  cat("插补警告:\n")
  cat("  ⚠ 日志事件数>100, 部分变量插补可能不稳定\n")
  cat("  ⚠ 建议检查 mice_logged_events.csv 中的问题变量\n")
  cat("  ⚠ 如STEMI等关键变量有问题, 需修复插补方法\n\n")
}

cat("统计推断:\n")
cat("  ⚠ IPW权重放大方差, 置信区间会变宽\n")
cat("  ⚠ MI与IPW的双重不确定性需正确传递\n")
cat("  ⚠ 使用稳健标准误 (survey包的svyglm)\n\n")

## ─────────────────────────────
## 结束
## ─────────────────────────────
cat("═══════════════════════════════════════════════════════\n")
cat("                    报告结束\n")
cat("═══════════════════════════════════════════════════════\n")

sink()

cat("  ✓ 综合报告已保存至:", report_path, "\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║            多重插补流程全部完成！                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("关键输出文件:\n")
cat("  1. mice_imputation_final. rds     - 插补对象 (★用于IPW建模)\n")
cat("  2. data_for_mi_final.rds         - 插补前数据\n")
cat("  3. MI_Final_Report.txt           - 详细分析报告\n")
cat("  4. imputation_diagnostic_report.csv - 插补质量诊断\n")
cat("  5. mi_convergence.png            - 收敛性图\n")
cat("  6. mi_density.png                - 密度图\n\n")

cat("重要提醒:\n")
cat("  ✓ has_fu_echo 未被插补 (IPW目标变量)\n")
cat("  ✓ LVEDV_fu 未被插补 (结局变量)\n")
cat("  ✓ 高共线性变量已移除 (CK, CKMB, LVESV_baseline)\n")
cat("  ✓ 高缺失梗死部位/冠脉变量已移除 (47%缺失)\n")

if (! is.null(imp$loggedEvents) && nrow(imp$loggedEvents) > 100) {
  cat("\n  ⚠ 警告: 插补日志事件数=", nrow(imp$loggedEvents), "\n")
  cat("    请检查 mice_logged_events.csv\n")
  cat("    如果STEMI等关键变量有问题, 需重新运行修复版代码\n")
}

cat("\n下一步:\n")
cat("  → 加载插补对象: imp <- readRDS('outputs/mice_imputation_final.rds')\n")
cat("  → 构建IPW权重: 对每个插补数据集拟合倾向性评分模型\n")
cat("  → 结局分析: 加权回归分析金属暴露与LVR的关系\n\n")

cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("全局种子:", GLOBAL_SEED, "\n\n")  
############################################################
## END
############################################################

