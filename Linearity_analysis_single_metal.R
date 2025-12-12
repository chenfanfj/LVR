############################################################
## 单金属主效应分析 - 步骤1: 线性关联
## 输入:  outcome_analysis_data_with_weights_extended.rds
## 输出: 19种金属的关联结果 + 森林图
############################################################

rm(list = ls())
gc()

library(dplyr)
library(survey)
library(mitools)
library(ggplot2)
library(here)
library(broom)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║     单金属主效应分析 - 线性关联                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 1: 加载数据
## ═══════════════════════════════════════════════════════

cat("【步骤 1】加载结局分析数据\n")

outcome_data <- readRDS(here:: here("outputs", "4 IPW2_complete", "outcome_analysis_data_with_weights_extended.rds"))

cat("  数据维度:", dim(outcome_data), "\n")
cat("  插补数据集数:", length(unique(outcome_data$. imp)), "\n")
cat("  分析样本数（每个插补）:", nrow(outcome_data) / length(unique(outcome_data$ .imp)), "\n\n")

# 验证结局变量
cat("  结局变量检查:\n")
cat("    delta_LVEDV:  ", sum(! is.na(outcome_data$delta_LVEDV)), "个非缺失\n")
cat("    范围: [", round(min(outcome_data$delta_LVEDV, na.rm=T), 1), ", ",
    round(max(outcome_data$delta_LVEDV, na.rm=T), 1), "]\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 定义金属列表和协变量
## ═══════════════════════════════════════════════════════

cat("【步骤 2】定义分析变量\n")

# 19种金属
metals <- c("Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
            "Mo", "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li")

cat("  金属数量:", length(metals), "\n")
cat("   ", paste(head(metals, 10), collapse=", "), "...\n\n")

# 协变量（基于文献和临床重要性）
covariates_core <- c(
  # 核心人口学
  "age", "gender",
  
  # 核心心功能
  "EF_baseline", "LVEDV_baseline",
  
  # 核心心肌损伤
  "cTnIpeak",
  
  # 核心治疗
  "pPCI", "STEMI"
)

covariates_extended <- c(
  covariates_core,
  
  # 扩展协变量
  "smoking", "DM", "hypertension",
  "NTproBNP_peak", "GRACE_in",
  "WBC", "HGB", "CRP",
  "CHOL", "LDL", "AST"
)

cat("  协变量设置:\n")
cat("    核心模型（Model 1）:", length(covariates_core), "个变量\n")
cat("    扩展模型（Model 2）:", length(covariates_extended), "个变量\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 对数转换验证
## ═══════════════════════════════════════════════════════

cat("【步骤 3】验证对数转换\n")

# 选择第一个插补数据集
dat_check <- outcome_data %>% filter(.imp == 1)

# 检查对数转换后的偏度
cat("  对数转换后的偏度（前5个金属）:\n")

for(metal in head(metals, 5)) {
  log_metal <- paste0("log_", metal)
  
  if(log_metal %in% names(dat_check)) {
    vals <- na.omit(dat_check[[log_metal]])
    skew <- e1071::skewness(vals)
    
    cat("    ", metal, ": 偏度 =", round(skew, 2))
    
    if(abs(skew) < 1) {
      cat(" ✓ (轻度偏态)\n")
    } else if(abs(skew) < 2) {
      cat(" (中度偏态)\n")
    } else {
      cat(" ⚠ (严重偏态)\n")
    }
  }
}

cat("\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 拟合单金属模型
## ═══════════════════════════════════════════════════════

cat("【步骤 4】拟合单金属模型\n")

# 准备存储结果
results_model1 <- data.frame()
results_model2 <- data.frame()

# 进度条
cat("  进度:  ")

for(i in 1:length(metals)) {
  metal <- metals[i]
  log_metal <- paste0("log_", metal)
  
  # 进度指示
  if(i %% 5 == 0) cat(i, "..  ")
  
  # 检查变量存在性
  if(! log_metal %in% names(outcome_data)) {
    cat("\n  ⚠ 跳过", metal, ": log变量不存在\n")
    next
  }
  
  ## ─────────────────────────────
  ## Model 1: 核心协变量
  ## ─────────────────────────────
  
  formula_m1 <- as.formula(paste0(
    "delta_LVEDV ~ ", log_metal, " + ",
    paste(covariates_core, collapse = " + ")
  ))
  
  # 对每个插补数据集拟合
  fit_list_m1 <- lapply(1:20, function(imp_i) {
    dat_i <- outcome_data %>% filter(.imp == imp_i)
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    tryCatch({
      svyglm(formula_m1, design = design_i, family = gaussian())
    }, error = function(e) {
      NULL
    })
  })
  
  # 移除NULL
  fit_list_m1 <- Filter(Negate(is.null), fit_list_m1)
  
  if(length(fit_list_m1) > 0) {
    # 合并
    pooled_m1 <- MIcombine(fit_list_m1)
    summary_m1 <- summary(pooled_m1)
    
    # 提取金属系数
    if(log_metal %in% rownames(summary_m1)) {
      beta_m1 <- summary_m1[log_metal, "results"]
      se_m1 <- summary_m1[log_metal, "se"]
      
      results_model1 <- rbind(results_model1, data.frame(
        Metal = metal,
        Model = "Model 1 (核心)",
        Beta = beta_m1,
        SE = se_m1,
        Lower95 = beta_m1 - 1.96 * se_m1,
        Upper95 = beta_m1 + 1.96 * se_m1,
        P_value = 2 * (1 - pnorm(abs(beta_m1 / se_m1))),
        N_imputations = length(fit_list_m1)
      ))
    }
  }
  
  ## ─────────────────────────────
  ## Model 2: 扩展协变量
  ## ─────────────────────────────
  
  formula_m2 <- as.formula(paste0(
    "delta_LVEDV ~ ", log_metal, " + ",
    paste(covariates_extended, collapse = " + ")
  ))
  
  fit_list_m2 <- lapply(1:20, function(imp_i) {
    dat_i <- outcome_data %>% filter(.imp == imp_i)
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    tryCatch({
      svyglm(formula_m2, design = design_i, family = gaussian())
    }, error = function(e) {
      NULL
    })
  })
  
  fit_list_m2 <- Filter(Negate(is.null), fit_list_m2)
  
  if(length(fit_list_m2) > 0) {
    pooled_m2 <- MIcombine(fit_list_m2)
    summary_m2 <- summary(pooled_m2)
    
    if(log_metal %in% rownames(summary_m2)) {
      beta_m2 <- summary_m2[log_metal, "results"]
      se_m2 <- summary_m2[log_metal, "se"]
      
      results_model2 <- rbind(results_model2, data.frame(
        Metal = metal,
        Model = "Model 2 (扩展)",
        Beta = beta_m2,
        SE = se_m2,
        Lower95 = beta_m2 - 1.96 * se_m2,
        Upper95 = beta_m2 + 1.96 * se_m2,
        P_value = 2 * (1 - pnorm(abs(beta_m2 / se_m2))),
        N_imputations = length(fit_list_m2)
      ))
    }
  }
}

cat("完成!\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 5: FDR校正
## ═══════════════════════════════════════════════════════

cat("【步骤 5】多重检验校正\n")

# 合并两个模型的结果
results_all <- rbind(results_model1, results_model2)

# FDR校正（分别对两个模型）
results_all <- results_all %>%
  group_by(Model) %>%
  mutate(
    P_FDR = p.adjust(P_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(Model, P_value)

cat("  Model 1: ", sum(results_model1$P_value < 0.05, na.rm=T), 
    "个显著 (P<0.05),",
    sum(results_model1$P_FDR < 0.10, na.rm=T), "个显著 (FDR<0.10)\n")

cat("  Model 2: ", sum(results_model2$P_value < 0.05, na.rm=T), 
    "个显著 (P<0.05),",
    sum(results_model2$P_FDR < 0.10, na.rm=T), "个显著 (FDR<0.10)\n\n")

# 保存结果
write.csv(results_all,
          here::here("outputs", "Single_Metal_Linear_Associations_Both_Models.csv"),
          row.names = FALSE)

cat("  ✓ 结果已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 生成森林图
## ═══════════════════════════════════════════════════════

cat("【步骤 6】生成森林图\n")

# 仅使用Model 2（扩展模型）绘图
results_plot <- results_all %>%
  filter(Model == "Model 2 (扩展)") %>%
  mutate(
    Significant = case_when(
      P_FDR < 0.05 ~ "FDR < 0.05",
      P_FDR < 0.10 ~ "FDR < 0.10",
      P_value < 0.05 ~ "P < 0.05",
      TRUE ~ "不显著"
    ),
    Significant = factor(Significant, 
                         levels = c("FDR < 0.05", "FDR < 0.10", 
                                    "P < 0.05", "不显著"))
  )

# 森林图
p_forest <- ggplot(results_plot, aes(x = Beta, y = reorder(Metal, Beta))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
  geom_errorbarh(aes(xmin = Lower95, xmax = Upper95, color = Significant),
                 height = 0.3, size = 0.9) +
  geom_point(aes(color = Significant, shape = Significant), 
             size = 4, alpha = 0.9) +
  scale_color_manual(
    values = c("FDR < 0.05" = "#D73027", 
               "FDR < 0.10" = "#FC8D59",
               "P < 0.05" = "#91BFDB",
               "不显著" = "gray60"),
    name = "显著性"
  ) +
  scale_shape_manual(
    values = c("FDR < 0.05" = 18, 
               "FDR < 0.10" = 17,
               "P < 0.05" = 16,
               "不显著" = 16),
    name = "显著性"
  ) +
  labs(
    x = "ΔLVEDV变化 (%) per log单位金属增加",
    y = "",
    title = "19种血清金属与左室重构的关联",
    subtitle = paste0("调整", length(covariates_extended), "个协变量 | ",
                      "IPW加权多重插补分析 (m=20)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11)
  )

ggsave(here::here("plots", "Forest_Plot_Single_Metal_Extended.png"), 
       p_forest,
       width = 12, height = 10, dpi = 300, bg = "white")

cat("  ✓ 森林图已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 7: 生成汇总表
## ═══════════════════════════════════════════════════════

cat("【步骤 7】生成汇总表\n")

# 识别显著金属
sig_metals_fdr <- results_all %>%
  filter(Model == "Model 2 (扩展)", P_FDR < 0.10) %>%
  pull(Metal)

sig_metals_p <- results_all %>%
  filter(Model == "Model 2 (扩展)", P_value < 0.05, ! Metal %in% sig_metals_fdr) %>%
  pull(Metal)

if(length(sig_metals_fdr) > 0) {
  cat("\n  显著金属 (FDR < 0.10):\n")
  for(metal in sig_metals_fdr) {
    res <- results_all %>% 
      filter(Metal == metal, Model == "Model 2 (扩展)")
    
    cat("    ", metal, ": Beta =", round(res$Beta, 3), 
        "(95% CI:", round(res$Lower95, 3), "-", round(res$Upper95, 3),
        "), P =", format.pval(res$P_value, digits=3),
        ", FDR =", round(res$P_FDR, 3), "\n")
  }
} else {
  cat("\n  无金属达到FDR < 0.10\n")
}

if(length(sig_metals_p) > 0) {
  cat("\n  边缘显著金属 (P < 0.05, FDR ≥ 0.10):\n")
  for(metal in sig_metals_p) {
    res <- results_all %>% 
      filter(Metal == metal, Model == "Model 2 (扩展)")
    
    cat("    ", metal, ": Beta =", round(res$Beta, 3), 
        "(", round(res$Lower95, 3), "-", round(res$Upper95, 3),
        "), P =", format.pval(res$P_value, digits=3), "\n")
  }
}

# 保存显著金属清单
sig_metals_all <- unique(c(sig_metals_fdr, sig_metals_p))

if(length(sig_metals_all) > 0) {
  write.csv(data.frame(Metal = sig_metals_all, 
                       Priority = c(rep("High", length(sig_metals_fdr)),
                                    rep("Medium", length(sig_metals_p)))),
            here::here("outputs", "Significant_Metals_for_RCS. csv"),
            row.names = FALSE)
  
  cat("\n  ✓ 显著金属清单已保存（用于RCS分析）\n")
}

cat("\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║         单金属线性关联分析完成！                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("输出文件:\n")
cat("  1. Single_Metal_Linear_Associations_Both_Models.csv\n")
cat("  2. Forest_Plot_Single_Metal_Extended.png\n")
cat("  3. Significant_Metals_for_RCS.csv\n\n")

cat("下一步:\n")
cat("  → 运行RCS分析代码（针对显著金属）\n")
cat("  → 运行分层分析代码（探索效应修饰）\n\n")

############################################################
## END
############################################################