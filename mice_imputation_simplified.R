############################################################
## 多重插补 (MI) - 精简版
## 用途: 为 IPW 倾向性评分建模准备插补数据
############################################################

## ═══════════════════════════════════════════════════════
## 环境设置
## ═══════════════════════════════════════════════════════
rm(list = ls())
gc()

GLOBAL_SEED <- 20240101
set.seed(GLOBAL_SEED)

cat("多重插补流程启动 | 种子:", GLOBAL_SEED, "\n\n")

## 加载包
packages <- c(
  "here", "readxl", "dplyr", "tidyr", "mice", "naniar", "VIM",
  "ggplot2", "corrplot"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

## 创建输出目录
outputs_dir <- here::here("outputs")
plots_dir <- here::here("plots")
dir.create(outputs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

## ═══════════════════════════════════════════════════════
## 步骤 1: 数据读取与清洗
## ═══════════════════════════════════════════════════════
cat("【步骤 1】数据读取与清洗\n")

data_path <- here::here("data", "AMI_PCI_metal_clean3.xlsx")
dat_raw <- readxl::read_excel(
  path = data_path,
  col_names = TRUE,
  na = c("", "NA", "N/A", "#N/A", "NULL", "null", "-999", "-9999")
)

cat("  原始数据:", nrow(dat_raw), "行 ×", ncol(dat_raw), "列\n")

dat <- dat_raw

## 标准化缺失值
for (col in names(dat)) {
  if (is.character(dat[[col]])) {
    dat[[col]][trimws(dat[[col]]) == ""] <- NA
  }
  if (is.numeric(dat[[col]])) {
    dat[[col]][dat[[col]] %in% c(-999, -9999, Inf, -Inf)] <- NA
  }
}

## 处理含不等号的变量
vars_with_symbols <- c(
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "NTproBNP_baseline", "NTproBNP_peak",
  "CK", "CKMB", "DD", "APTT", "PT_sec", "TT", "INR", "FIB",
  "FT3", "FT4", "S_TSH", "IL_6", "CRP", "h_CT"
)

for (var in intersect(vars_with_symbols, names(dat))) {
  if (is.character(dat[[var]])) {
    dat[[var]] <- gsub("＜|<|﹤|≤|＞|>|﹥|≥", "", dat[[var]])
    dat[[var]] <- gsub(",", "", dat[[var]])
    dat[[var]] <- trimws(dat[[var]])
    dat[[var]] <- suppressWarnings(as.numeric(dat[[var]]))
  }
}

## 创建目标变量
dat <- dat %>%
  mutate(
    has_fu_echo = factor(
      if_else(!is.na(LVEDV_fu), 1, 0),
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

## Winsorize 极端值
continuous_vars <- c(
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

for (var in intersect(continuous_vars, names(dat))) {
  if (is.numeric(dat[[var]])) {
    valid_vals <- dat[[var]][!is.na(dat[[var]]) & is.finite(dat[[var]])]
    if (length(valid_vals) > 10) {
      q1 <- quantile(valid_vals, 0.25)
      q3 <- quantile(valid_vals, 0.75)
      iqr <- q3 - q1
      lower_fence <- q1 - 3 * iqr
      upper_fence <- q3 + 3 * iqr
      dat[[var]][dat[[var]] < lower_fence & !is.na(dat[[var]])] <- lower_fence
      dat[[var]][dat[[var]] > upper_fence & !is.na(dat[[var]])] <- upper_fence
    }
  }
}

cat("  ✓ 数据清洗完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 变量选择
## ═══════════════════════════════════════════════════════
cat("【步骤 2】变量选择\n")

outcome_vars <- c("LVESV_fu", "EF_fu")

exposure_vars <- c(
  "Cu", "Zn", "Fe", "Se", "Pb",
  "Al", "As", "Cr", "Mn", "Ni", "Mo", 
  "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li"
)

ipw_covariates <- c(
  "age", "gender", "resident", "Career",
  "DM", "hypertension", "smoking", "Cancer", "his_stroke", "AF",
  "STEMI", "ST_dev", "ST_dep",
  "pPCI", "Stent_no", "Lesion_no",
  "Cardio_shock", "VF", "Cardiac_arrest_in",
  "IN_killip", "OUT_killip", "GRACE_in", "GRACE_in_str",
  "EF_baseline", "LVEDV_baseline",
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "CKMB",
  "NTproBNP_baseline", "NTproBNP_peak",
  "WBC", "HGB", "PLT", "RBC", "HCT", "NE", "LY", "MO",
  "SCR", "EGFR", "BUN",
  "CRP", "IL_6"
)

auxiliary_vars <- c(
  "height", "weight",
  "HR", "SBP", "DBP",
  "ALT", "AST", "TBIL",
  "CHOL", "TG", "LDL", "HDL",
  "GLU", "HbAlc",
  "DD", "FIB", "APTT",
  "K", "Na", "Ca", "Mg",
  "Aspirin", "Statin", "ACEIorARB", "β_block"
)

all_imputation_vars <- unique(c(
  "has_fu_echo",
  outcome_vars,
  exposure_vars,
  ipw_covariates,
  auxiliary_vars
))

exclude_vars <- c(
  "BMI", "delta_LVEDV", "AST_ALT", "BUN_CREA", "APROB_APROA",
  "LVEDV_fu",
  "LVESV_baseline", "CK",
  "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",
  "LAD", "LCX", "RCA", "LM",
  "TL_LAD", "TL_LCX", "TL_RCA",
  "ID"
)

imputation_vars <- setdiff(
  intersect(all_imputation_vars, names(dat)),
  exclude_vars
)

dat_for_mi <- dat %>%
  select(ID, LVEDV_fu, LVESV_fu, EF_fu, all_of(imputation_vars))

cat("  插补变量数:", length(imputation_vars), "\n")

## 缺失模式分析
miss_summary <- data.frame(
  variable = names(dat_for_mi),
  n_missing = colSums(is.na(dat_for_mi)),
  pct_missing = round(colMeans(is.na(dat_for_mi)) * 100, 2)
) %>%
  filter(n_missing > 0) %>%
  arrange(desc(pct_missing))

write.csv(miss_summary, 
          file.path(outputs_dir, "missing_summary.csv"),
          row.names = FALSE)

## 应用分层缺失率阈值
THRESHOLD_IPW_CORE <- 5
THRESHOLD_IPW_AUX <- 20
THRESHOLD_OUTCOME <- 40
THRESHOLD_AUXILIARY <- 50

ipw_core_vars <- c(
  "age", "gender", "pPCI", "STEMI", "ST_dev",
  "cTnIpeak", "CKMB", "NTproBNP_peak",
  "EF_baseline", "LVEDV_baseline",
  "GRACE_in_str", "hypertension", "resident"
)

ipw_auxiliary_vars <- c(
  "DM", "smoking", "Career", "AF", "his_stroke",
  "Stent_no", "Lesion_no", "IN_killip", "OUT_killip",
  "cTnI_baseline", "NTproBNP_baseline",
  "WBC", "HGB", "PLT", "SCR", "EGFR",
  "CRP", "IL_6"
)

high_miss_to_remove <- miss_summary %>%
  filter(
    (variable %in% ipw_core_vars & pct_missing > THRESHOLD_IPW_CORE) |
      (variable %in% ipw_auxiliary_vars & pct_missing > THRESHOLD_IPW_AUX) |
      (variable %in% outcome_vars & pct_missing > THRESHOLD_OUTCOME) |
      (pct_missing > THRESHOLD_AUXILIARY),
    !variable %in% c("has_fu_echo", "LVESV_fu", "EF_fu", "LVEDV_fu")
  ) %>%
  pull(variable)

if (length(high_miss_to_remove) > 0) {
  dat_for_mi <- dat_for_mi %>% select(-all_of(high_miss_to_remove))
  imputation_vars <- setdiff(imputation_vars, high_miss_to_remove)
}

cat("  最终保留:", ncol(dat_for_mi) - 1, "个变量\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 变量质量诊断
## ═══════════════════════════════════════════════════════
cat("【步骤 3】变量质量诊断\n")

## 检查常数变量
vars_to_remove <- character()

for (var in names(dat_for_mi)) {
  if (var %in% c("ID", "has_fu_echo")) next
  
  if (is.numeric(dat_for_mi[[var]])) {
    vals <- na.omit(dat_for_mi[[var]])
    if (length(vals) > 0) {
      n_unique <- length(unique(vals))
      sd_val <- sd(vals, na.rm = TRUE)
      if (n_unique == 1 || is.na(sd_val) || sd_val < 1e-10) {
        vars_to_remove <- c(vars_to_remove, var)
      }
    }
  }
}

if (length(vars_to_remove) > 0) {
  dat_for_mi <- dat_for_mi %>% select(-all_of(vars_to_remove))
}

## 共线性检查
numeric_vars <- names(dat_for_mi)[sapply(dat_for_mi, is.numeric)]
numeric_vars <- setdiff(numeric_vars, "ID")

if (length(numeric_vars) > 1) {
  cor_result <- tryCatch(
    cor(dat_for_mi[, numeric_vars], use = "pairwise.complete.obs"),
    error = function(e) NULL
  )
  
  if (!is.null(cor_result)) {
    high_cor_vars <- character()
    for (i in 1:(nrow(cor_result) - 1)) {
      for (j in (i + 1):ncol(cor_result)) {
        r <- cor_result[i, j]
        if (!is.na(r) && abs(r) > 0.95) {
          var1 <- rownames(cor_result)[i]
          var2 <- colnames(cor_result)[j]
          miss1 <- mean(is.na(dat_for_mi[[var1]]))
          miss2 <- mean(is.na(dat_for_mi[[var2]]))
          var_to_remove <- if (miss1 > miss2) var1 else var2
          high_cor_vars <- c(high_cor_vars, var_to_remove)
        }
      }
    }
    
    if (length(high_cor_vars) > 0) {
      high_cor_vars <- unique(high_cor_vars)
      dat_for_mi <- dat_for_mi %>% select(-all_of(high_cor_vars))
    }
  }
}

cat("  ✓ 质量诊断完成\n\n")

saveRDS(dat_for_mi, file.path(outputs_dir, "data_for_mi_cleaned.rds"))

## ═══════════════════════════════════════════════════════
## 步骤 4: 配置多重插补
## ═══════════════════════════════════════════════════════
cat("【步骤 4】配置多重插补\n")

## 识别需要插补的变量
missing_counts <- data.frame(
  variable = names(dat_for_mi),
  n_missing = colSums(is.na(dat_for_mi)),
  stringsAsFactors = FALSE
)

vars_need_imputation <- missing_counts %>%
  filter(
    n_missing > 0,
    !variable %in% c("ID", "has_fu_echo", "delta_LVEDV")
  ) %>%
  pull(variable)

vars_complete <- missing_counts %>%
  filter(n_missing == 0 | variable %in% c("ID", "has_fu_echo")) %>%
  pull(variable)

## 分离数据
dat_complete <- dat_for_mi %>% select(all_of(vars_complete))

complete_predictors <- intersect(
  c("age", "gender", "STEMI", "pPCI", "ST_dev", 
    "EF_baseline", "LVEDV_baseline", "hypertension", "DM"),
  vars_complete
)

dat_for_imputation <- dat_for_mi %>%
  select(ID, all_of(vars_need_imputation), all_of(complete_predictors))

## 配置插补方法
ini <- mice(dat_for_imputation, maxit = 0, print = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

meth["ID"] <- ""
pred[, "ID"] <- 0
pred["ID", ] <- 0

for (var in complete_predictors) {
  meth[var] <- ""
  pred[var, ] <- 0
}

## 结局变量不插补
outcomes_to_skip <- c("LVEDV_fu", "LVESV_fu", "EF_fu")
for (var in outcomes_to_skip) {
  if (var %in% names(meth)) {
    meth[var] <- ""
    pred[var, ] <- 0
    pred[, var] <- 0
  }
}

## 其他变量使用 PMM
binary_vars <- c("smoking", "Cancer", "his_stroke", "AF", 
                 "VF", "Cardio_shock", "Cardiac_arrest_in",
                 "resident", "Career", "Aspirin", "Statin", 
                 "ACEIorARB", "β_block")

for (var in vars_need_imputation) {
  if (var %in% binary_vars) {
    meth[var] <- "pmm"
  } else {
    meth[var] <- "pmm"
  }
}

## 被动插补 BMI
if ("BMI" %in% names(dat_for_imputation) && 
    all(c("height", "weight") %in% names(dat_for_imputation))) {
  meth["BMI"] <- "~I(weight / (height/100)^2)"
  pred["BMI", ] <- 0
}

## 优化预测矩阵
pred <- quickpred(
  dat_for_imputation,
  mincor = 0.2,
  minpuc = 0.3,
  include = complete_predictors,
  exclude = c("ID")
)

for (cpv in complete_predictors) {
  pred[cpv, ] <- 0
}

cat("  ✓ 配置完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 5: 执行插补
## ═══════════════════════════════════════════════════════
cat("【步骤 5】执行多重插补\n")

m <- 20
maxit <- 10

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

cat("  ✓ 插补完成\n\n")

## 合并完整变量与插补数据
imp_datasets_full <- lapply(1:m, function(i) {
  imp_data_i <- complete(imp, i)
  dat_full_i <- imp_data_i %>%
    left_join(dat_complete, by = "ID")
  return(dat_full_i)
})

## 重构 mids 对象
imp_long <- mice::complete(imp, action = "long", include = TRUE)

complete_vars_to_add <- setdiff(names(dat_complete), c("ID", names(imp_long)))

if (length(complete_vars_to_add) > 0) {
  for (var in complete_vars_to_add) {
    imp_long[[var]] <- rep(dat_complete[[var]], m + 1)
  }
}

imp_full <- as.mids(imp_long)

## 保存
saveRDS(imp_full, file.path(outputs_dir, "mice_imputation_final.rds"))
saveRDS(imp_datasets_full, file.path(outputs_dir, "imputed_datasets_list.rds"))

cat("  ✓ 插补对象已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 插补质量诊断
## ═══════════════════════════════════════════════════════
cat("【步骤 6】插补质量诊断\n")

imp_data_1 <- imp_datasets_full[[1]]

imputation_diagnostic <- data.frame(
  variable = intersect(names(dat_for_mi), names(imp_data_1)),
  original_missing = sapply(intersect(names(dat_for_mi), names(imp_data_1)), 
                           function(v) sum(is.na(dat_for_mi[[v]]))),
  after_imputation_missing = sapply(intersect(names(dat_for_mi), names(imp_data_1)), 
                                   function(v) sum(is.na(imp_data_1[[v]])))
) %>%
  filter(original_missing > 0) %>%
  mutate(
    imputation_rate = round((1 - after_imputation_missing / original_missing) * 100, 1)
  )

write.csv(imputation_diagnostic, 
          file.path(outputs_dir, "imputation_diagnostic.csv"),
          row.names = FALSE)

## 生成诊断图
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

if (length(imputed_numeric_vars) > 0) {
  plot_vars <- head(imputed_numeric_vars, 20)
  
  png(file.path(plots_dir, "mi_convergence.png"),
      width = 14, height = 12, units = "in", res = 300)
  print(plot(imp, plot_vars, layout = c(4, 5)))
  dev.off()
  
  if (length(imputed_numeric_vars) >= 9) {
    density_vars <- head(imputed_numeric_vars, 9)
    png(file.path(plots_dir, "mi_density.png"),
        width = 12, height = 9, units = "in", res = 300)
    print(densityplot(imp, 
                      as.formula(paste("~", paste(density_vars, collapse = " + "))),
                      layout = c(3, 3)))
    dev.off()
  }
}

cat("  ✓ 诊断完成\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════
cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║            多重插补流程完成！                        ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("核心输出:\n")
cat("  • mice_imputation_final.rds - 插补对象\n")
cat("  • imputed_datasets_list.rds - 插补数据集列表\n")
cat("  • imputation_diagnostic.csv - 插补质量报告\n\n")

cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
