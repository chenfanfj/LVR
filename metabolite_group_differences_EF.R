# ==============================================================================
# 分析1.2_EF版本: 代谢物组间差异分析 (基于EF定义的LVR)
# LVR定义: EF从≥50%降至<50%
# 方法: IPW加权logistic回归 + 多重插补池化
# ==============================================================================

rm(list = ls())
gc()
library(dplyr)
library(broom)
library(mice)
library(ggplot2)
library(ggrepel)
library(survey)
library(readxl)

# ------------------------------------------------------------------------------
# 步骤1: 数据准备
# ------------------------------------------------------------------------------
cat("===== 步骤1: 数据准备 (修正版) =====\n")

# 1.1 读取多重插补数据
mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")
cat("- 多重插补数据已加载\n")

# 1.2 读取候选代谢物 (修正 atomic vector 报错)
candidates <- readRDS("outputs/13 candidate_metabolites/final_metabolite_list.rds")
if (is.list(candidates) && "variable_names" %in% names(candidates)) {
  metabo_vars <- candidates$variable_names
} else {
  metabo_vars <- as.character(candidates) # 如果是向量，直接转换
}
cat("- 候选代谢物数:", length(metabo_vars), "\n")

# 1.3 读取代谢物名称映射 (根据 csv 修正列名为 MNO)
metabo_mapping <- read_excel("data/metabolism.xlsx", sheet = "original") %>%
  # 根据 file_difference.csv，original 表使用的是 MNO 和 name 列
  select(MNO, name) %>% 
  rename(NO = MNO) # 重命名为 NO 以匹配后续代码

# ------------------------------------------------------------------------------
# 代谢物数据 log2 转化
# ------------------------------------------------------------------------------
cat("\n===== 预处理: 代谢物数据 log2 转化 =====\n")
long_data <- complete(mids_obj, action = "long", include = TRUE)
long_data <- long_data %>%
  mutate(across(all_of(metabo_vars), ~log2(.x + 1)))
mids_obj <- as.mids(long_data)
cat("- 代谢物已完成 log2(x+1) 转化\n")


# ------------------------------------------------------------------------------
# 步骤2: 计算基于EF的LVR结局
# ------------------------------------------------------------------------------
cat("\n===== 步骤2: 计算LVR结局 (EF定义) =====\n")

# 检查第一个数据集
data_check <- complete(mids_obj, 1)

# 验证必需变量
required_vars <- c("EF_baseline", "EF_fu")
missing_vars <- setdiff(required_vars, names(data_check))

if(length(missing_vars) > 0) {
  stop("❌ 缺少必需变量: ", paste(missing_vars, collapse = ", "))
}

# 计算LVR (EF从≥50%降至<50%)
long_data <- complete(mids_obj, action = "long", include = TRUE)

long_data <- long_data %>%
  mutate(
    LVR_EF = case_when(
      is.na(EF_baseline) | is.na(EF_fu) ~ NA_real_,
      EF_baseline < 50 ~ NA_real_, # 排除基线不符合条件的样本
      EF_baseline >= 50 & EF_fu < 50 ~ 1,
      EF_baseline >= 50 & EF_fu >= 50 ~ 0,
      TRUE ~ 0
    )
  )

# 转回mids对象
mids_obj <- as.mids(long_data)

# 检查LVR分布
data_check <- complete(mids_obj, 1)
lvr_table <- table(data_check$LVR_EF, useNA = "ifany")
cat("- LVR_EF分布:\n")
print(lvr_table)
cat("- LVR发生率:", round(mean(data_check$LVR_EF, na.rm=TRUE) * 100, 1), "%\n")

if(sum(data_check$LVR_EF == 1, na.rm=TRUE) < 10) {
  stop("❌ LVR事件数过少 (<10), 无法进行logistic回归")
}

outcome_var <- "LVR_EF"

# ------------------------------------------------------------------------------
# 步骤3: 定义协变量 (针对EF结局优化)
# ------------------------------------------------------------------------------
cat("\n===== 步骤3: 定义协变量 =====\n")

# 核心协变量 (与EF下降相关)
covariates <- c(
  "age",              # 年龄
  "gender",           # 性别 (或 sex)
  "EF_baseline",      # 基线EF (关键!)
  "LVEDV_baseline",   # 基线左室容积
  "cTnIpeak",         # 峰值肌钙蛋白
  "NTproBNP_peak",    # 峰值NT-proBNP
  "pPCI",             # 急诊PCI
  "STEMI",            # ST段抬高型心梗
  "hypertension",     # 高血压
  "DM",               # 糖尿病 (或 diabetes)
  "smoking",          # 吸烟
  "GRACE_in"          # GRACE评分
)

# 检查变量存在性并替换
if("sex" %in% names(data_check) & !"gender" %in% names(data_check)) {
  covariates[covariates == "gender"] <- "sex"
}

if("diabetes" %in% names(data_check) & !"DM" %in% names(data_check)) {
  covariates[covariates == "DM"] <- "diabetes"
}

covariates <- covariates[covariates %in% names(data_check)]
cat("- 实际使用协变量:", length(covariates), "个\n")
cat("  ", paste(covariates, collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 步骤4: 识别IPW权重变量
# ------------------------------------------------------------------------------
cat("\n===== 步骤4: 识别权重变量 =====\n")

possible_weights <- c("sw_trunc", "ipw", "weight", "weights", "iptw", "sw", "ps_weight")
weight_var <- intersect(possible_weights, names(data_check))[1]

if(!is.na(weight_var)) {
  use_weights <- TRUE
  cat(sprintf("✅ 成功识别权重变量: '%s'\n", weight_var))
} else {
  use_weights <- FALSE
  weight_var <- NULL
  cat("⚠️ 未找到权重变量,将进行未加权分析\n")
}

# ------------------------------------------------------------------------------
# 共线性检查 (VIF)
# ------------------------------------------------------------------------------
cat("\n===== 诊断: 协变量共线性检查 (VIF) =====\n")
library(car)
test_data <- complete(mids_obj, 1)
vif_formula <- as.formula(paste0(outcome_var, " ~ ", paste(covariates, collapse = " + ")))
vif_fit <- glm(vif_formula, data = test_data, family = binomial())
print(vif(vif_fit)) 

# ------------------------------------------------------------------------------
# 步骤 5: 拟合函数 
# ------------------------------------------------------------------------------
fit_metabolite_model <- function(mids_obj, metabolite, outcome, covariates, 
                                 weight_var = NULL, use_weights = TRUE) {
  formula_obj <- as.formula(paste0(outcome, " ~ ", metabolite, " + ", 
                                   paste(covariates, collapse = " + ")))
  
  M <- mids_obj$m
  coef_list <- vector("list", M)
  se_list <- vector("list", M)
  n_obs_list <- vector("list", M)
  
  for(m in 1:M) {
    data_m <- complete(mids_obj, m)
    data_m <- data_m[complete.cases(data_m[, c(outcome, metabolite, covariates)]), ]
    
    # --- 建议 4: 权重截断 (Trimming) ---
    if(use_weights && !is.null(weight_var)) {
      w_limits <- quantile(data_m[[weight_var]], c(0.01, 0.99), na.rm = TRUE)
      data_m[[weight_var]] <- pmax(pmin(data_m[[weight_var]], w_limits[2]), w_limits[1])
      
      design_m <- svydesign(ids = ~1, weights = as.formula(paste0("~", weight_var)), data = data_m)
      fit_m <- svyglm(formula_obj, design = design_m, family = quasibinomial())
      res <- tidy(fit_m) %>% filter(term == metabolite)
      coef_list[[m]] <- res$estimate
      se_list[[m]] <- res$std.error
    } else {
      fit_m <- glm(formula_obj, data = data_m, family = binomial())
      res <- tidy(fit_m) %>% filter(term == metabolite)
      coef_list[[m]] <- res$estimate
      se_list[[m]] <- res$std.error
    }
    n_obs_list[[m]] <- nrow(data_m)
  }
  
  # --- 建议 3: Barnard-Rubin 调整后的池化 ---
  b_coefs <- unlist(coef_list)
  b_ses <- unlist(se_list)
  
  pooled_coef <- mean(b_coefs)
  within_var <- mean(b_ses^2)
  between_var <- var(b_coefs)
  total_var <- within_var + (1 + 1/M) * between_var
  pooled_se <- sqrt(total_var)
  
  # 修正自由度
  lambda <- (between_var + between_var/M) / total_var
  df_old <- (M - 1) / (lambda^2)
  df_obs <- mean(unlist(n_obs_list)) - length(covariates) - 2
  df_adj <- (df_old * df_obs) / (df_old + df_obs)
  
  p_val <- 2 * pt(-abs(pooled_coef / pooled_se), df = df_adj)
  
  return(data.frame(
    Metabolite = metabolite, Beta = pooled_coef, SE = pooled_se,
    OR = exp(pooled_coef), OR_lower = exp(pooled_coef - 1.96*pooled_se),
    OR_upper = exp(pooled_coef + 1.96*pooled_se), P_value = p_val
  ))
}

# 5.2 批量拟合
cat("- 开始拟合", length(metabo_vars), "个代谢物模型.. .\n")

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
    results_list[[i]] <- data.frame(
      Metabolite = metabolite,
      Beta = NA, SE = NA, OR = NA,
      OR_lower = NA, OR_upper = NA, P_value = NA, N = NA
    )
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

# 5.3 合并结果
metabo_diff_results <- bind_rows(results_list) %>%
  filter(! is.na(P_value))

cat("\n✅ 成功拟合", nrow(metabo_diff_results), "个模型\n")

# ------------------------------------------------------------------------------
# 步骤6: FDR校正 & 效应量计算
# ------------------------------------------------------------------------------
cat("\n===== 步骤6: 多重检验校正 =====\n")

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
    Log2FC = log2(OR)  # Log2 Fold Change
  ) %>%
  arrange(FDR, desc(abs(Log2FC)))

# 统计显著性
cat("- FDR < 0.001:", sum(metabo_diff_results$FDR < 0.001, na.rm=TRUE), "个\n")
cat("- FDR < 0.01: ", sum(metabo_diff_results$FDR < 0.01, na.rm=TRUE), "个\n")
cat("- FDR < 0.05: ", sum(metabo_diff_results$FDR < 0.05, na.rm=TRUE), "个\n")
cat("- FDR < 0.10: ", sum(metabo_diff_results$FDR < 0.10, na.rm=TRUE), "个\n")

# ------------------------------------------------------------------------------
# 步骤7: 保存结果
# ------------------------------------------------------------------------------
cat("\n===== 步骤7: 保存结果 =====\n")

# 创建输出目录
dir.create("outputs/EF_metabolite", showWarnings = FALSE, recursive = TRUE)
dir.create("plots/EF_metabolite", showWarnings = FALSE, recursive = TRUE)

# 7.1 完整结果
write.csv(metabo_diff_results, 
          "outputs/EF_metabolite/metabolite_group_differences_EF_full.csv",
          row.names = FALSE)

# 7.2 显著差异代谢物
metabo_sig <- metabo_diff_results %>%
  filter(FDR < 0.10) %>%
  arrange(FDR)

write.csv(metabo_sig, 
          "outputs/EF_metabolite/metabolite_group_differences_EF_significant.csv",
          row.names = FALSE)

cat("✅ 结果已保存至:  outputs/EF_metabolite\n")
cat("   显著代谢物数:", nrow(metabo_sig), "\n")

# ------------------------------------------------------------------------------
# 步骤8: 增强火山图
# ------------------------------------------------------------------------------
cat("\n===== 步骤8: 绘制火山图 =====\n")

# 8.1 添加代谢物名称
metabo_diff_enhanced <- metabo_diff_results %>%
  left_join(
    metabo_mapping %>% rename(Metabolite = NO, Metabolite_Name = name),
    by = "Metabolite"
  ) %>%
  mutate(
    # 优先显示名称,否则显示编号
    Display_Name = ifelse(! is.na(Metabolite_Name) & Metabolite_Name != "", 
                          Metabolite_Name, Metabolite),
    # 标记显著且效应量大的代谢物
    Label = ifelse(FDR < 0.10 & abs(Log2FC) > 0.2, Display_Name, ""),
    # 分类
    Category = case_when(
      FDR < 0.05 & Log2FC > 0.2 ~ "Up in LVR (FDR<0.05)",
      FDR < 0.05 & Log2FC < -0.2 ~ "Down in LVR (FDR<0.05)",
      FDR < 0.10 & abs(Log2FC) > 0.2 ~ "Trend (FDR<0.10)",
      TRUE ~ "NS"
    )
  )

# 8.2 绘制火山图
p_volcano <- ggplot(metabo_diff_enhanced, 
                    aes(x = Log2FC, y = -log10(P_value))) +
  
  # 散点
  geom_point(aes(color = Category), size = 2.5, alpha = 0.7) +
  
  # 阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", 
             color = "grey40", linewidth = 0.5) +
  
  # 标注
  geom_text_repel(aes(label = Label), 
                  size = 3, 
                  max.overlaps = 30,
                  box.padding = 0.5,
                  segment.color = "grey50") +
  
  # 颜色
  scale_color_manual(
    values = c(
      "Up in LVR (FDR<0.05)" = "#E64B35FF",
      "Down in LVR (FDR<0.05)" = "#4DBBD5FF",
      "Trend (FDR<0.10)" = "#FFD700",
      "NS" = "grey70"
    )
  ) +
  
  # 标签
  labs(
    title = "Volcano Plot: Metabolites Associated with LVR",
    subtitle = paste0("LVR Definition: EF ≥50% → <50%  |  N=", 
                      nrow(data_check), 
                      ", LVR Rate=", 
                      round(mean(data_check$LVR_EF, na.rm=TRUE)*100, 1), "%"),
    x = "Log2 Odds Ratio",
    y = "-log10(P value)",
    color = ""
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )

ggsave("plots/EF_metabolite/volcano_plot_metabolites_EF.pdf", 
       p_volcano, width = 10, height = 8, dpi = 300)
ggsave("plots/EF_metabolite/volcano_plot_metabolites_EF.png", 
       p_volcano, width = 10, height = 8, dpi = 300)

cat("✅ 火山图已保存至: plots/EF_metabolite/\n")

# ------------------------------------------------------------------------------
# 步骤9: TOP代谢物条形图
# ------------------------------------------------------------------------------

if(nrow(metabo_sig) > 0) {
  
  top_n <- min(20, nrow(metabo_sig))
  
  metabo_top <- metabo_diff_enhanced %>%
    filter(FDR < 0.10) %>%
    arrange(FDR) %>%
    slice(1:top_n) %>%
    mutate(
      Direction = ifelse(Log2FC > 0, "↑ LVR", "↓ LVR"),
      Display_Name = factor(Display_Name, levels = Display_Name)
    )
  
  p_bar <- ggplot(metabo_top, aes(x = Log2FC, y = reorder(Display_Name, Log2FC))) +
    geom_bar(stat = "identity", aes(fill = Direction), width = 0.7) +
    geom_errorbarh(aes(xmin = log2(OR_lower), xmax = log2(OR_upper)), 
                   height = 0.3, color = "grey30") +
    scale_fill_manual(values = c("↑ LVR" = "#E64B35FF", "↓ LVR" = "#4DBBD5FF")) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    labs(
      title = paste0("Top ", top_n, " Metabolites Associated with LVR"),
      subtitle = "Sorted by FDR, Error bars = 95% CI",
      x = "Log2 Odds Ratio",
      y = "",
      fill = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30")
    )
  
  ggsave("plots/EF_metabolite/top_metabolites_barplot_EF.pdf", 
         p_bar, width = 10, height = max(6, top_n*0.4), dpi = 300)
  ggsave("plots/EF_metabolite/top_metabolites_barplot_EF.png", 
         p_bar, width = 10, height = max(6, top_n*0.4), dpi = 300)
  
  cat("✅ TOP代谢物条形图已保存\n")
}

# ------------------------------------------------------------------------------
# 步骤10: 保存增强版结果
# ------------------------------------------------------------------------------

write.csv(metabo_diff_enhanced, 
          "outputs/EF_metabolite/metabolite_group_differences_EF_with_names.csv",
          row.names = FALSE)

cat("✅ 增强版结果(含代谢物名称)已保存\n")

# 打印TOP10结果
cat("\n========== TOP 10 代谢物 (按FDR排序) ==========\n")
print(metabo_diff_enhanced %>%
        select(Display_Name, OR, OR_lower, OR_upper, P_value, FDR, Significance) %>%
        head(10))

cat("\n========== 分析1. 2 (EF定义) 完成 ==========\n")