# ==============================================================================
# 脚本名称: final_verification_IPW_LVR.R
# 目的: 对 delta_EF 筛选出的潜力代谢物进行二元 LVR 结局的最终验证
# 统计方法: IPW加权 + 逻辑回归 + 多重插补池化 (Rubin's Rules)
# ==============================================================================

rm(list = ls())
gc()
library(dplyr)
library(broom)
library(mice)
library(survey)
library(readxl)
library(ggplot2)

# ------------------------------------------------------------------------------
# 步骤 1: 加载数据与候选名单提取
# ------------------------------------------------------------------------------
cat("===== 步骤 1: 加载数据与候选名单 =====\n")

# 1.1 加载 MWAS 结果并提取 P < 0.1 的代谢物 (确保路径正确)
mwas_path <- "outputs/EF_metabolite/delta_ef_selection/mwas_delta_ef_full_results.csv"
if(!file.exists(mwas_path)) stop("❌ 找不到 MWAS 结果文件，请检查路径。")

mwas_results <- read.csv(mwas_path)

# 排除非代谢物列 (临床指标)
clinical_exclude <- c("LVEDV_fu", "DD", "LVEDV_baseline", "NTproBNP_baseline", "ST_dep", "Zn")
sig_metabolites <- mwas_results %>%
  filter(Pooled_P < 0.1 & !(term %in% clinical_exclude)) %>%
  pull(term)

cat("- 提取到 P < 0.1 的潜力代谢物:", length(sig_metabolites), "个\n")

# 1.2 加载插补数据与代谢物名称映射
mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")
metabo_mapping <- read_excel("data/metabolism.xlsx", sheet = "original") %>%
  select(MNO, name) %>% rename(Metabolite = MNO)

# ------------------------------------------------------------------------------
# 步骤 2: 结局定义与数据预处理
# ------------------------------------------------------------------------------
cat("\n===== 步骤 2: 结局定义与数据预处理 =====\n")

long_data <- complete(mids_obj, action = "long", include = TRUE)
long_data <- long_data %>%
  mutate(
    LVR_EF = case_when(
      is.na(EF_baseline) | is.na(EF_fu) ~ NA_real_,
      EF_baseline >= 50 & EF_fu < 50 ~ 1,
      EF_baseline >= 50 & EF_fu >= 50 ~ 0,
      TRUE ~ NA_real_
    ),
    # 进行 log2 转化
    across(all_of(sig_metabolites), ~log2(.x + 1))
  )
mids_obj <- as.mids(long_data)

# --- 权重变量识别 (修复找不到 sw_trunc 的关键步骤) ---
check_data <- complete(mids_obj, 1)
# 自动识别权重变量：按优先级寻找常见命名
possible_weights <- c("sw_trunc", "ipw", "weight", "weights", "iptw", "sw")
weight_var <- intersect(possible_weights, names(check_data))[1]

if(is.na(weight_var)) {
  cat("⚠️ 未在数据中找到权重变量，将进行【未加权】回归分析。\n")
  use_weights <- FALSE
} else {
  cat("✅ 成功识别权重变量:", weight_var, "\n")
  use_weights <- TRUE
}

# ------------------------------------------------------------------------------
# 步骤 3: 运行 IPW 逻辑回归验证
# ------------------------------------------------------------------------------
cat("\n===== 步骤 3: 运行逻辑回归验证 =====\n")

fixed_covs <- c("age", "gender", "EF_baseline", "cTnIpeak", "GRACE_in")
fixed_covs <- intersect(fixed_covs, names(check_data))
M <- mids_obj$m

run_final_model <- function(metabo_name) {
  formula_str <- paste0("LVR_EF ~ ", metabo_name, " + ", paste(fixed_covs, collapse = " + "))
  
  coef_list <- vector("numeric", M)
  se_list <- vector("numeric", M)
  
  for(i in 1:M) {
    data_i <- complete(mids_obj, i) %>% filter(!is.na(LVR_EF))
    
    if(use_weights) {
      # 排除权重缺失样本
      data_i <- data_i %>% filter(!is.na(!!sym(weight_var)))
      # 权重截断 (0.99分位数)
      w_upper <- quantile(data_i[[weight_var]], 0.99, na.rm = TRUE)
      data_i[[weight_var]] <- pmin(data_i[[weight_var]], w_upper)
      
      design_i <- svydesign(ids = ~1, weights = as.formula(paste0("~", weight_var)), data = data_i)
      fit_i <- svyglm(as.formula(formula_str), design = design_i, family = quasibinomial())
    } else {
      fit_i <- glm(as.formula(formula_str), data = data_i, family = binomial())
    }
    
    res <- tidy(fit_i) %>% filter(term == metabo_name)
    coef_list[i] <- res$estimate
    se_list[i] <- res$std.error
  }
  
  # Rubin's Rules 池化
  pooled_beta <- mean(coef_list)
  within_var <- mean(se_list^2)
  between_var <- var(coef_list)
  total_var <- within_var + (1 + 1/M) * between_var
  pooled_se <- sqrt(total_var)
  
  # 修正自由度 (Barnard-Rubin)
  df_obs <- nrow(data_i) - length(fixed_covs) - 2
  lambda <- (between_var + between_var/M) / total_var
  df_adj <- ((M-1)/lambda^2 * df_obs) / ((M-1)/lambda^2 + df_obs)
  
  p_val <- 2 * pt(-abs(pooled_beta/pooled_se), df = df_adj)
  
  return(data.frame(
    Metabolite = metabo_name, Beta = pooled_beta, SE = pooled_se,
    OR = exp(pooled_beta), OR_lower = exp(pooled_beta - 1.96*pooled_se),
    OR_upper = exp(pooled_beta + 1.96*pooled_se), P_value = p_val
  ))
}

final_results <- bind_rows(lapply(sig_metabolites, run_final_model))

# ------------------------------------------------------------------------------
# 步骤 4: 结果输出与森林图绘制
# ------------------------------------------------------------------------------
cat("\n===== 步骤 4: 输出结果与绘制森林图 =====\n")

final_results <- final_results %>%
  left_join(metabo_mapping, by = "Metabolite") %>%
  mutate(
    Display_Name = ifelse(!is.na(name), name, Metabolite),
    Sig_Label = case_when(P_value < 0.05 ~ "*", P_value < 0.1 ~ "†", TRUE ~ ""),
    Display_Name = factor(Display_Name, levels = rev(Display_Name[order(OR)]))
  )

write.csv(final_results, "outputs/EF_metabolite/final_verification_LVR_results.csv", row.names = FALSE)

# 绘制森林图

p_forest <- ggplot(final_results, aes(x = OR, y = Display_Name)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2) +
  geom_point(aes(size = -log10(P_value), color = (P_value < 0.05))) +
  geom_text(aes(label = Sig_Label), vjust = -0.5, size = 5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(title = "Final Verification: Metabolites & LVR", x = "Odds Ratio (95% CI)", y = "") +
  theme_minimal()

dir.create("plots/EF_metabolite", recursive = TRUE, showWarnings = FALSE)
ggsave("plots/EF_metabolite/final_verification_forest_plot.png", p_forest, width = 10, height = 8)

cat("✅ 分析完成，结果已保存至 outputs/ 和 plots/ 文件夹。\n")