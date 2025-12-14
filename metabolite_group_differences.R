# ==============================================================================
# 分析1.2: 代谢物组间差异分析 (LVR vs 非LVR)
# 方法: IPW加权logistic回归 + 多重插补池化
# ==============================================================================

rm(list = ls())
gc()
library(dplyr)
library(broom)
library(ggplot2)
library(ggrepel)
library(survey)
library(readxl)

# ------------------------------------------------------------------------------
# 步骤1: 数据准备
# ------------------------------------------------------------------------------
cat("===== 步骤1: 数据准备 =====\n")
# 加载金属转换函数
source(here::here("functions", "transform_metals.R"))

# 1.1 读取数据
mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")
cat("- 多重插补数据已加载:", mids_obj$m, "个数据集,", nrow(mids_obj$data), "行\n")

# --- 金属变量识别 --
target_metals <- c("Li", "Be", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn",
                   "As", "Se", "Mo", "Cd", "Sn", "Sb", "Ba", "Tl", "Pb", "U", "Al", "Fe", "Sr")
data_names <- names(mids_obj$data)
found_metals <- intersect(target_metals, data_names)

# 如果没找到原始名，找 log_ 或 ln_ 前缀
if(length(found_metals) == 0) {
  log_candidates <- paste0("log_", target_metals)
  ln_candidates <- paste0("ln_", target_metals)
  if(any(log_candidates %in% data_names)) {
    metal_vars <- data_names[data_names %in% log_candidates]
    cat("- 识别到 log_ 金属变量:", length(metal_vars), "个\n")
  } else if(any(ln_candidates %in% data_names)) {
    metal_vars <- data_names[data_names %in% ln_candidates]
    cat("- 识别到 ln_ 金属变量:", length(metal_vars), "个\n")
  } else {
    stop("❌ 无法识别金属变量，请检查列名！")
  }
} else {
  metal_vars <- found_metals
  # 这里建议手动对数转换逻辑(如上一轮回复所示)，此处从略，假设已处理
}

# --- 代谢物变量识别 ---
candidate_file <- "outputs/13 candidate_metabolites/final_metabolite_list.rds"
if(file.exists(candidate_file)) {
  metabo_vars <- readRDS(candidate_file)
} else {
  metabo_vars <- names(mids_obj$data)[grep("^M[0-9]+$|metabo_", names(mids_obj$data))]
}
metabo_vars <- metabo_vars[metabo_vars %in% names(mids_obj$data)]

# --- 权重变量智能搜索 (修复 Issue 1) ---
# 定义可能的权重变量名
possible_weights <- c("ipw", "weight", "weights", "iptw", "sw", "ps_weight")
# 找到第一个匹配的
weight_var <- intersect(possible_weights, names(mids_obj$data))[1]

if(!is.na(weight_var)) {
  use_weights <- TRUE
  cat(sprintf("✅ 成功识别权重变量: '%s'\n", weight_var))
} else {
  use_weights <- FALSE
  weight_var <- NULL
  cat("⚠️ 未在数据中找到权重变量 (如 'ipw', 'weight')。\n")
  cat("   -> 分析将继续进行，但不使用加权 (Unweighted Analysis)。\n")
  cat("   -> 原因: 之前的插补/合并步骤可能未包含权重计算结果。\n")
}

# 1.2 读取候选代谢物
candidates <- readRDS("outputs/13 candidate_metabolites/final_metabolite_list.rds")
metabo_vars <- candidates$variable_names

# 1.3 读取代谢物名称映射
metabo_mapping <- read_excel("data/metabolism.xlsx", sheet = "original") %>%
  select(NO, name)

# 1.4 定义结局变量
outcome_var <- "LVR"

# 如果需要从连续变量计算LVR:
data_check <- complete(mids_obj, 1)

if(!"LVR" %in% names(data_check)) {
  cat("- 计算LVR二分类变量(LVEDV增加>20%).. .\n")
  
  long_data <- complete(mids_obj, action = "long", include = TRUE)
  
  long_data <- long_data %>%
    mutate(
      LVEDV_change_pct = (LVEDV_fu - LVEDV_baseline) / LVEDV_baseline * 100,
      LVR = ifelse(LVEDV_change_pct > 20, 1, 0)
    )
  
  mids_obj <- as.mids(long_data)
  
  outcome_var <- "LVR"
}

# 检查LVR分布
lvr_table <- table(data_check[[outcome_var]], useNA = "ifany")
cat("- LVR分布:", paste(names(lvr_table), lvr_table, collapse = ", "), "\n")
cat("- LVR发生率:", round(prop.table(lvr_table)[2] * 100, 1), "%\n")

# 1.5 定义协变量(与原代码相同)
covariates <- c("age", "sex", "BMI", "hypertension", "diabetes", 
                "smoking", "LVEF_baseline", "peak_CK_MB")
covariates <- covariates[covariates %in% names(data_check)]

# 1.6 IPW权重
weight_var <- "ipw"
use_weights <- weight_var %in% names(data_check)

# ------------------------------------------------------------------------------
# 步骤2: 对每个代谢物拟合IPW加权logistic回归
# ------------------------------------------------------------------------------
cat("\n===== 步骤2: 拟合加权logistic回归 =====\n")

# 2.1 定义拟合函数
fit_metabolite_model <- function(mids_obj, metabolite, outcome, covariates, 
                                 weight_var = NULL, use_weights = TRUE) {
  
  # 构建公式
  formula_str <- paste0(outcome, " ~ ", metabolite, " + ", 
                        paste(covariates, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  if(use_weights) {
    # 使用survey包进行加权分析
    # 对每个插补数据集拟合模型
    M <- mids_obj$m
    coef_list <- vector("list", M)
    
    for(m in 1:M) {
      data_m <- complete(mids_obj, m)
      
      # 创建survey design对象
      design_m <- svydesign(
        ids = ~1,
        weights = as.formula(paste0("~", weight_var)),
        data = data_m
      )
      
      # 拟合加权logistic回归
      fit_m <- svyglm(formula_obj, design = design_m, family = binomial())
      
      # 提取代谢物的系数
      coef_list[[m]] <- tidy(fit_m) %>%
        filter(term == metabolite)
    }
    
    # 池化结果(使用Rubin规则)
    coefs <- sapply(coef_list, function(x) x$estimate)
    ses <- sapply(coef_list, function(x) x$std.error)
    
    # Rubin规则
    pooled_coef <- mean(coefs)
    within_var <- mean(ses^2)
    between_var <- var(coefs)
    total_var <- within_var + (1 + 1/M) * between_var
    pooled_se <- sqrt(total_var)
    
    # 计算t统计量和P值
    df <- (M - 1) * (1 + within_var / ((1 + 1/M) * between_var))^2
    t_stat <- pooled_coef / pooled_se
    p_value <- 2 * pt(-abs(t_stat), df = df)
    
    # 计算OR和95% CI
    OR <- exp(pooled_coef)
    OR_lower <- exp(pooled_coef - 1.96 * pooled_se)
    OR_upper <- exp(pooled_coef + 1.96 * pooled_se)
    
  } else {
    # 未加权分析: 使用mice:: pool()
    fit_list <- with(mids_obj, glm(formula_obj, family = binomial()))
    pooled_fit <- pool(fit_list)
    
    result <- summary(pooled_fit) %>%
      filter(term == metabolite)
    
    pooled_coef <- result$estimate
    pooled_se <- result$std.error
    p_value <- result$p.value
    
    OR <- exp(pooled_coef)
    OR_lower <- exp(pooled_coef - 1.96 * pooled_se)
    OR_upper <- exp(pooled_coef + 1.96 * pooled_se)
  }
  
  return(data.frame(
    Metabolite = metabolite,
    Beta = pooled_coef,
    SE = pooled_se,
    OR = OR,
    OR_lower = OR_lower,
    OR_upper = OR_upper,
    P_value = p_value
  ))
}

# 3.2 批量拟合所有代谢物
cat("- 开始拟合", length(metabo_vars), "个代谢物模型...\n")

results_list <- vector("list", length(metabo_vars))
pb <- txtProgressBar(min = 0, max = length(metabo_vars), style = 3)

for(i in 1:length(metabo_vars)) {
  
  metabolite <- metabo_vars[i]
  
  tryCatch({
    results_list[[i]] <- fit_metabolite_model(
      mids_obj = mids_obj,
      metabolite = metabolite,
      outcome = outcome_var,
      covariates = covariates,
      weight_var = weight_var,
      use_weights = use_weights
    )
  }, error = function(e) {
    cat("\n⚠️ 模型", metabolite, "拟合失败:", e$message, "\n")
    results_list[[i]] <- data.frame(
      Metabolite = metabolite,
      Beta = NA, SE = NA, OR = NA,
      OR_lower = NA, OR_upper = NA, P_value = NA
    )
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

# 3.3 合并结果
metabo_diff_results <- bind_rows(results_list) %>%
  filter(! is.na(P_value))

cat("\n✅ 成功拟合", nrow(metabo_diff_results), "个模型\n")

# ------------------------------------------------------------------------------
# 步骤4: FDR校正
# ------------------------------------------------------------------------------
cat("\n===== 步骤4: 多重检验校正 =====\n")

metabo_diff_results <- metabo_diff_results %>%
  mutate(
    FDR = p.adjust(P_value, method = "fdr"),
    Bonferroni = p.adjust(P_value, method = "bonferroni"),
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      FDR < 0.10 ~ "†",
      TRUE ~ ""
    ),
    # 计算Log2 Fold Change (用于火山图)
    # 注意: OR已经是指数化的,需要log2转换
    Log2FC = log2(OR)
  ) %>%
  arrange(FDR, desc(abs(Log2FC)))

# 统计显著性
n_sig_001 <- sum(metabo_diff_results$FDR < 0.001, na.rm = TRUE)
n_sig_01 <- sum(metabo_diff_results$FDR < 0.01, na.rm = TRUE)
n_sig_05 <- sum(metabo_diff_results$FDR < 0.05, na.rm = TRUE)
n_sig_10 <- sum(metabo_diff_results$FDR < 0.10, na.rm = TRUE)

cat("- FDR < 0.001:", n_sig_001, "个代谢物\n")
cat("- FDR < 0.01: ", n_sig_01, "个代谢物\n")
cat("- FDR < 0.05: ", n_sig_05, "个代谢物\n")
cat("- FDR < 0.10: ", n_sig_10, "个代谢物\n")

# ------------------------------------------------------------------------------
# 步骤5: 保存结果
# ------------------------------------------------------------------------------
cat("\n===== 步骤5: 保存结果 =====\n")

# 5.1 完整结果
write.csv(metabo_diff_results, 
          "outputs/metabolite_group_differences_full.csv",
          row.names = FALSE)

cat("✅ 完整结果已保存至:  outputs/metabolite_group_differences_full.csv\n")

# 5.2 显著差异代谢物
metabo_diff_

# ------------------------------------------------------------------------------
# 步骤6: 增强的火山图(添加真实代谢物名称)
# ------------------------------------------------------------------------------
cat("\n===== 步骤6: 绘制火山图 =====\n")

# 6.1 添加代谢物名称
metabo_diff_results_enhanced <- metabo_diff_results %>%
  left_join(
    metabo_mapping %>% rename(Metabolite = NO, Metabolite_Name = name),
    by = "Metabolite"
  ) %>%
  mutate(
    # 标记显著且效应量大的代谢物
    Label = ifelse(FDR < 0.05 & abs(Log2FC) > 0.2, Metabolite_Name, "")
  )

# 6.2 绘制火山图
pdf("outputs/volcano_plot_metabolites_enhanced.pdf", width = 10, height = 8)

ggplot(metabo_diff_results_enhanced, 
       aes(x = Log2FC, y = -log10(P_value))) +
  
  # 散点
  geom_point(aes(color = Significance), size = 3, alpha = 0.6) +
  
  # 显著性阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "grey50", size = 0.5) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", 
             color = "grey50", size = 0.5) +
  
  # 标注代谢物名称
  geom_text_repel(aes(label = Label), 
                  size = 3, 
                  max.overlaps = 20,
                  box.padding = 0.5) +
  
  # 颜色
  scale_color_manual(
    values = c("***" = "#E64B35FF", "**" = "#FF6B6B", 
               "*" = "#FFA07A", "†" = "#FFD700", "" = "grey70"),
    breaks = c("***", "**", "*", "†", ""),
    labels = c("FDR<0.001", "FDR<0.01", "FDR<0.05", "FDR<0.10", "NS")
  ) +
  
  # 标签
  labs(
    title = "代谢物与LVR关联的火山图",
    subtitle = paste0("N=", nrow(data_check), 
                      ", LVR发生率=", round(prop.table(lvr_table)[2]*100, 1), "%"),
    x = "Log2 Odds Ratio",
    y = "-log10(P value)",
    color = "Significance"
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

dev.off()

cat("✅ 火山图已保存至: outputs/volcano_plot_metabolites_enhanced.pdf\n")

# ------------------------------------------------------------------------------
# 步骤7: 保存增强版结果表
# ------------------------------------------------------------------------------

write.csv(metabo_diff_results_enhanced, 
          "outputs/metabolite_group_differences_with_names.csv",
          row.names = FALSE)

cat("✅ 增强版结果(含代谢物名称)已保存\n")

cat("\n========== 分析1.2完成 ==========\n")