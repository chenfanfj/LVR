# ==============================================================================
# 脚本名称: mwas_continuous_delta_ef_selection.R
# 目的: 基于多重插补数据的全代谢组关联分析 (连续变量结局: ΔEF)
# 策略: 单变量线性回归 + 协变量校正 + 结果池化 (Rubin's Rules)
# ==============================================================================

rm(list = ls())
gc()
library(dplyr)
library(broom)
library(mice)

# ------------------------------------------------------------------------------
# 步骤 1: 加载数据与连续结局定义
# ------------------------------------------------------------------------------
cat("===== 步骤 1: 数据准备与 ΔEF 定义 =====\n")

# 1.1 读取插补数据
mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")

# 1.2 在 mids 对象中定义 ΔEF (随访 EF - 基线 EF)
long_data <- complete(mids_obj, action = "long", include = TRUE)
long_data <- long_data %>%
  mutate(
    # 严格遵循定义：仅分析基线 EF 正常的患者 (>= 50)
    delta_EF = ifelse(EF_baseline >= 50, EF_fu - EF_baseline, NA_real_)
  )
mids_obj <- as.mids(long_data)

# 1.3 确定数值型代谢物变量池 (排除非数值列和临床变量)
check_data <- complete(mids_obj, 1)
exclude_cols <- c("ID", "NO", "MNO", "name", "age", "gender", "sex", 
                  "EF_baseline", "EF_fu", "LVEDV_baseline", "cTnIpeak", 
                  "NTproBNP_peak", "pPCI", "STEMI", "hypertension", "DM", 
                  "diabetes", "smoking", "GRACE_in", "LVR_EF", "delta_EF", ".imp", ".id")

metabo_pool <- names(check_data)[sapply(check_data, is.numeric)]
metabo_pool <- setdiff(metabo_pool, exclude_cols)

# 核心校正协变量
fixed_covs <- c("age", "gender", "EF_baseline", "cTnIpeak", "GRACE_in")
fixed_covs <- intersect(fixed_covs, names(check_data))

cat("- 样本总数:", nrow(check_data), "\n")
cat("- 代谢物池总数:", length(metabo_pool), "\n")
cat("- 结局变量: ΔEF (随访 - 基线)\n")

# ------------------------------------------------------------------------------
# 步骤 2: 运行单变量 MWAS (Linear Regression)
# ------------------------------------------------------------------------------
cat("\n===== 步骤 2: 运行全代谢组关联分析 (Linear MWAS) =====\n")

M <- mids_obj$m
mwas_results_list <- list()

for(i in 1:M) {
  cat(sprintf("\r正在处理第 %d/%d 个插补集...", i, M))
  
  data_i <- complete(mids_obj, i) %>%
    filter(!is.na(delta_EF)) %>%
    # 对代谢物进行 log2 转化，使回归系数 Beta 代表“浓度翻倍”后的 EF 变化量
    mutate(across(all_of(metabo_pool), ~log2(.x + 1)))
  
  # 对每个代谢物运行线性回归
  res_i <- lapply(metabo_pool, function(m) {
    form <- as.formula(paste0("delta_EF ~ ", m, " + ", paste(fixed_covs, collapse = " + ")))
    
    tryCatch({
      fit <- lm(form, data = data_i)
      tidy(fit) %>% 
        filter(term == m) %>% 
        select(term, estimate, std.error)
    }, error = function(e) return(NULL))
  })
  
  mwas_results_list[[i]] <- bind_rows(res_i)
}

# ------------------------------------------------------------------------------
# 步骤 3: 使用 Rubin's Rules 池化结果
# ------------------------------------------------------------------------------
cat("\n\n===== 步骤 3: 结果池化与筛选 =====\n")

final_mwas <- bind_rows(mwas_results_list) %>%
  group_by(term) %>%
  summarise(
    Pooled_Beta = mean(estimate),
    # 计算总方差 (组内方差 + 组间方差修正)
    Within_Var = mean(std.error^2),
    Between_Var = var(estimate),
    Total_Var = Within_Var + (1 + 1/M) * Between_Var,
    Pooled_SE = sqrt(Total_Var)
  ) %>%
  mutate(
    T_stat = Pooled_Beta / Pooled_SE,
    # 计算近似 P 值 (基于 t 分布)
    Pooled_P = 2 * pt(-abs(T_stat), df = nrow(check_data) - length(fixed_covs) - 2),
    # FDR 校正
    FDR = p.adjust(Pooled_P, method = "fdr")
  ) %>%
  arrange(Pooled_P)

# ------------------------------------------------------------------------------
# 步骤 4: 保存与输出
# ------------------------------------------------------------------------------
output_dir <- "outputs/EF_metabolite/delta_ef_selection/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 筛选前 50 名代谢物作为候选名单
final_candidates <- final_mwas %>% head(50)

write.csv(final_mwas, paste0(output_dir, "mwas_delta_ef_full_results.csv"), row.names = FALSE)
saveRDS(as.character(final_candidates$term), paste0(output_dir, "delta_ef_final_metabolites.rds"))

cat("- MWAS 分析完成。已选取 ΔEF 关联最显著的前", nrow(final_candidates), "个代谢物。\n")
print(head(final_mwas, 15))