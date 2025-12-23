# ==============================================================================
# 脚本名称: mwas_binary_lvr_selection.R
# 目的: 基于多重插补数据的全代谢组关联分析 (二元结局 LVR_EF)
# 策略: 单变量逻辑回归 + 协变量校正 + 结果池化
# ==============================================================================

rm(list = ls())
gc()
library(dplyr)
library(broom)
library(mice)

# ------------------------------------------------------------------------------
# 步骤 1: 加载数据与结局定义
# ------------------------------------------------------------------------------
cat("===== 步骤 1: 数据准备与结局定义 =====\n")

mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")

# 在 mids 对象中定义 LVR_EF (确保排除基线 EF < 50 的人群)
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

# 确定数值型代谢物变量池
check_data <- complete(mids_obj, 1)
exclude_cols <- c("ID", "NO", "MNO", "name", "age", "gender", "sex", 
                  "EF_baseline", "EF_fu", "LVEDV_baseline", "cTnIpeak", 
                  "NTproBNP_peak", "pPCI", "STEMI", "hypertension", "DM", 
                  "diabetes", "smoking", "GRACE_in", "LVR_EF", ".imp", ".id")

metabo_pool <- names(check_data)[sapply(check_data, is.numeric)]
metabo_pool <- setdiff(metabo_pool, exclude_cols)

# 核心校正协变量
fixed_covs <- c("age", "gender", "EF_baseline", "cTnIpeak", "GRACE_in")
fixed_covs <- intersect(fixed_covs, names(check_data))

cat("- 代谢物池总数:", length(metabo_pool), "\n")
cat("- 协变量已校正:", paste(fixed_covs, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 步骤 2: 运行单变量 MWAS (Logistic Regression)
# ------------------------------------------------------------------------------
cat("\n===== 步骤 2: 运行全代谢组关联分析 (Logistic MWAS) =====\n")

M <- mids_obj$m
mwas_results_list <- list()

for(i in 1:M) {
  cat(sprintf("\r正在处理第 %d/%d 个插补集...", i, M))
  
  data_i <- complete(mids_obj, i) %>%
    filter(!is.na(LVR_EF)) %>%
    # 对代谢物进行 log2 转化以符合逻辑回归假设
    mutate(across(all_of(metabo_pool), ~log2(.x + 1)))
  
  # 对每个代谢物运行逻辑回归
  res_i <- lapply(metabo_pool, function(m) {
    # 构建公式: LVR_EF ~ 代谢物 + 临床协变量
    form <- as.formula(paste0("LVR_EF ~ ", m, " + ", paste(fixed_covs, collapse = " + ")))
    
    tryCatch({
      fit <- glm(form, data = data_i, family = binomial())
      # 提取该代谢物的系数、标准误和 P 值
      tidy(fit) %>% 
        filter(term == m) %>% 
        select(term, estimate, std.error, p.value)
    }, error = function(e) return(NULL))
  })
  
  mwas_results_list[[i]] <- bind_rows(res_i)
}

# ------------------------------------------------------------------------------
# 步骤 3: 结果池化与筛选
# ------------------------------------------------------------------------------
cat("\n\n===== 步骤 3: 汇总池化结果 =====\n")

# 使用 Rubin's Rules 的简化版进行池化 (取中位数 P 值或计算池化 P)
final_mwas <- bind_rows(mwas_results_list) %>%
  group_by(term) %>%
  summarise(
    Pooled_Beta = mean(estimate),
    Pooled_SE = sqrt(mean(std.error^2) + (1 + 1/M) * var(estimate)), # Rubin's Rules
    Median_P = median(p.value)
  ) %>%
  mutate(
    Pooled_P = 2 * pt(-abs(Pooled_Beta / Pooled_SE), df = 100), # 近似自由度
    OR = exp(Pooled_Beta)
  ) %>%
  arrange(Pooled_P)

# ------------------------------------------------------------------------------
# 步骤 4: 保存候选名单
# ------------------------------------------------------------------------------
# 筛选标准：取 Pooled_P 前 50 名 (或者 Pooled_P < 0.05)
final_candidates <- final_mwas %>% 
  filter(Pooled_P < 0.1) %>% # 先放宽到 0.1 以确保有足够的数量
  head(50)

output_dir <- "outputs/EF_metabolite/lasso_selection/" # 复用目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(final_mwas, paste0(output_dir, "mwas_binary_full_results.csv"), row.names = FALSE)
saveRDS(as.character(final_candidates$term), paste0(output_dir, "lasso_final_metabolites.rds"))

cat("- MWAS 分析完成。已选取前", nrow(final_candidates), "个代谢物进入下一步回归分析。\n")
print(head(final_mwas, 15))