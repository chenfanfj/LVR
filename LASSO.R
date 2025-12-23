# ==============================================================================
# 脚本名称: lasso_metabolite_selection_MICE.R
# 目的: 基于多重插补数据进行 LASSO 变量筛选
# 核心逻辑: 在每个插补数据集运行 LASSO，统计变量选中频率
# ==============================================================================
# ==============================================================================
# 脚本名称: lasso_metabolite_selection_MICE.R
# ==============================================================================

rm(list = ls())
gc()
library(glmnet)
library(dplyr)
library(mice)
library(readxl)

# ------------------------------------------------------------------------------
# 步骤 1: 加载数据与结局定义
# ------------------------------------------------------------------------------
cat("===== 步骤 1: 数据准备与结局定义 =====\n")

# 1.1 读取插补数据
mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")

# 1.2 在 mids 对象中定义 LVR_EF
long_data <- complete(mids_obj, action = "long", include = TRUE)

long_data <- long_data %>%
  mutate(
    LVR_EF = case_when(
      is.na(EF_baseline) | is.na(EF_fu) ~ NA_real_,
      EF_baseline >= 50 & EF_fu < 50 ~ 1,
      EF_baseline >= 50 & EF_fu >= 50 ~ 0,
      TRUE ~ NA_real_ 
    )
  )

mids_obj <- as.mids(long_data)
cat("- 结局变量 LVR_EF 已定义\n")

# 1.3: 细化变量池
check_data <- complete(mids_obj, 1)

# 获取所有列名用于后续交叉校验 (修复 all_vars 缺失问题)
all_column_names <- names(check_data)

# 定义必须排除的非代谢物列名
exclude_cols <- c("ID", "NO", "MNO", "name", "age", "gender", "sex", 
                  "EF_baseline", "EF_fu", "LVEDV_baseline", "cTnIpeak", 
                  "NTproBNP_peak", "pPCI", "STEMI", "hypertension", "DM", 
                  "diabetes", "smoking", "GRACE_in", "LVR_EF", ".imp", ".id")

# 自动筛选：必须是数值型列，且不在排除名单内
metabo_pool <- all_column_names[sapply(check_data, is.numeric)]
metabo_pool <- setdiff(metabo_pool, exclude_cols)

cat("- 最终确认的数值型代谢物候选池总数:", length(metabo_pool), "\n")

# 1.4 定义强制保留的协变量 (修复 all_vars 报错)
fixed_covs <- c("age", "gender", "EF_baseline", "cTnIpeak", "GRACE_in")
# 将原先错误的 all_vars 替换为定义过的 all_column_names
fixed_covs <- intersect(fixed_covs, all_column_names) 

# ------------------------------------------------------------------------------
# 步骤 2: 运行 LASSO 稳定性筛选
# ------------------------------------------------------------------------------
cat("\n===== 步骤 2: 运行 LASSO 稳定性筛选 =====\n")

M <- mids_obj$m
lasso_selected_list <- list()



for(i in 1:M) {
  cat(sprintf("\r正在处理第 %d/%d 个数据集...", i, M))
  
  data_i <- complete(mids_obj, i) %>%
    filter(!is.na(LVR_EF))
  
  # log2 转化 (仅针对 metabo_pool 中的数值列)
  data_i <- data_i %>% mutate(across(all_of(metabo_pool), ~log2(.x + 1)))
  
  # 构建矩阵
  X <- as.matrix(data_i[, c(fixed_covs, metabo_pool)])
  Y <- data_i$LVR_EF
  
  # 惩罚因子：协变量权重为 0 (不惩罚)，代谢物为 1 (惩罚)
  p_factor <- c(rep(0, length(fixed_covs)), rep(1, length(metabo_pool)))
  
  set.seed(123 + i)
  try({
    cv_fit <- cv.glmnet(X, Y, family = "binomial", alpha = 1, 
                        penalty.factor = p_factor, nfolds = 10, type.measure = "auc")
    
    coeffs <- coef(cv_fit, s = "lambda.min")
    non_zero_vars <- rownames(coeffs)[which(coeffs != 0)]
    
    lasso_selected_list[[i]] <- setdiff(non_zero_vars, c("(Intercept)", fixed_covs))
  }, silent = TRUE)
}

# ------------------------------------------------------------------------------
# 步骤 3: 汇总筛选频率 (增加空结果检查)
# ------------------------------------------------------------------------------
cat("\n\n===== 步骤 3: 汇总筛选频率 =====\n")

# 检查是否有代谢物被选中
if (length(unlist(lasso_selected_list)) == 0) {
  stop("❌ 警告：所有插补数据集中的 LASSO 均未选出任何代谢物。请检查：
       1. 结局 LVR_EF 是否分布极度不均？
       2. 代谢物数据是否存在大量缺失或零值？
       3. 考虑将结局更换为连续变量 (ΔEF) 进行筛选。")
}

selection_summary <- unlist(lasso_selected_list) %>%
  table() %>%
  as.data.frame(stringsAsFactors = FALSE) 

# 安全地重命名列
colnames(selection_summary) <- c("Metabolite", "Frequency")

selection_summary <- selection_summary %>%
  mutate(Selection_Rate = (Frequency / M) * 100) %>%
  arrange(desc(Frequency))

# 筛选入选频率阈值
threshold_rate <- 50 # 如果结果太少，可以降至 30
final_candidates <- selection_summary %>% filter(Selection_Rate >= threshold_rate)

# 创建目录并保存
output_dir <- "outputs/EF_based/lasso_selection/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(selection_summary, paste0(output_dir, "metabolite_selection_frequency.csv"), row.names = FALSE)
saveRDS(as.character(final_candidates$Metabolite), paste0(output_dir, "lasso_final_metabolites.rds"))

cat("- 最终入选代谢物数 (频率 >= ", threshold_rate, "%):", nrow(final_candidates), "\n")
print(head(selection_summary, 15))

cat("\n✅ 筛选完成。最终候选代谢物列表已保存至 outputs/EF_metabolite/lasso_selection/\n")