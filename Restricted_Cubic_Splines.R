############################################################
##  A1: 非线性关系分析 (Restricted Cubic Splines)
##  输入:  imputed_data_with_ipw_weights_extended.rds
##  输出: RCS曲线图 + 非线性检验结果
############################################################

rm(list = ls())
gc()

library(mice)
library(survey)
library(mitools)
library(rms)
library(aod)
library(ggplot2)
library(dplyr)
library(patchwork)
library(here)

# 创建输出目录
dir.create(here("outputs", "EF"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("plots", "EF"), showWarnings = FALSE, recursive = TRUE)

## ═══════════════════════════════════════════════════════
## 步骤 1: 数据准备 (保持之前的严谨性)
## ═══════════════════════════════════════════════════════

imputed_full <- readRDS(here("outputs", "4 IPW2_complete", "imputed_data_with_ipw_weights_extended.rds"))

core_metals <- c("Pb", "Zn", "Fe", "Se", "Cu", "As")
covariates  <- c("age", "gender", "EF_baseline", "LVEDV_baseline", 
                 "cTnIpeak", "pPCI", "STEMI", "smoking", "DM", 
                 "hypertension", "NTproBNP_peak", "GRACE_in", 
                 "WBC", "HGB", "CRP", "CHOL", "LDL", "AST")

imputed_processed <- imputed_full %>%
  filter(has_fu_echo == "Yes", !is.na(LVEDV_fu)) %>%
  mutate(
    LVR = as.numeric(EF_baseline >= 50 & EF_fu < 50),
    across(c(gender, pPCI, STEMI, smoking, DM, hypertension), as.factor)
  )

imp_list <- split(imputed_processed, imputed_processed$.imp)
imp_list <- imp_list[names(imp_list) != "0"]

# 设置 Data Dist
df_for_dd <- imp_list[[1]] %>% select(all_of(covariates), all_of(paste0("log_", core_metals)), LVR)
dd <- datadist(df_for_dd)
options(datadist = "dd")

## ═══════════════════════════════════════════════════════
## 步骤 2: 核心分析与优化绘图
## ═══════════════════════════════════════════════════════

summary_results <- data.frame()
plot_list <- list()

for (metal in core_metals) {
  cat("\n>>> 分析中:", metal)
  log_metal <- paste0("log_", metal)
  
  # A. 拟合 RCS (svyglm)
  fits_rcs <- lapply(imp_list, function(df) {
    design <- svydesign(ids = ~1, weights = ~sw_trunc, data = df)
    formula_rcs <- as.formula(paste("LVR ~ rcs(", log_metal, ", 4) +", paste(covariates, collapse = " + ")))
    svyglm(formula_rcs, design = design, family = quasibinomial())
  })
  
  # B. 合并结果与 D1 非线性检验
  combined_fit <- MIcombine(fits_rcs)
  coef_names <- names(coef(combined_fit))
  nonlin_indices <- grep(paste0("rcs\\(", log_metal, ".*\\'"), coef_names)
  w_test <- wald.test(b = coef(combined_fit), Sigma = vcov(combined_fit), Terms = nonlin_indices)
  p_nonlin <- w_test$result$chi2["P"]
  
  # C. 汇总线性效应
  fits_linear <- lapply(imp_list, function(df) {
    design <- svydesign(ids = ~1, weights = ~sw_trunc, data = df)
    formula_lin <- as.formula(paste("LVR ~", log_metal, "+", paste(covariates, collapse = " + ")))
    svyglm(formula_lin, design = design, family = quasibinomial())
  })
  combined_lin <- MIcombine(fits_linear)
  lin_or <- exp(coef(combined_lin)[log_metal])
  lin_ci <- exp(confint(combined_lin)[log_metal, ])
  
  knots_pos <- attr(rcs(imp_list[[1]][[log_metal]], 4), "parms")
  summary_results <- rbind(summary_results, data.frame(
    Metal = metal,
    Linear_OR_95CI = sprintf("%.2f (%.2f-%.2f)", lin_or, lin_ci[1], lin_ci[2]),
    P_Nonlinearity = p_nonlin,
    Knots_Original = paste(round(exp(knots_pos), 2), collapse = ", ")
  ))
  
  # D. 绘图数据准备 (深度聚焦优化)
  df_plot <- imp_list[[1]]
  test_data <- data.frame(
    x_val = seq(quantile(df_plot[[log_metal]], 0.05, na.rm=T), 
                quantile(df_plot[[log_metal]], 0.95, na.rm=T), length.out = 100)
  )
  names(test_data) <- log_metal
  for(cv in covariates) {
    if(is.numeric(df_plot[[cv]])) test_data[[cv]] <- median(df_plot[[cv]], na.rm=T)
    else test_data[[cv]] <- levels(df_plot[[cv]])[1]
  }
  
  # 提取预测值
  pred_obj <- predict(fits_rcs[[1]], newdata = test_data, type = "link")
  test_data$fit_link <- as.numeric(pred_obj)
  test_data$se_link  <- as.numeric(SE(pred_obj))
  
  test_data$prob  <- plogis(test_data$fit_link)
  test_data$lower <- plogis(test_data$fit_link - 1.96 * test_data$se_link)
  test_data$upper <- plogis(test_data$fit_link + 1.96 * test_data$se_link)
  test_data$metal_orig <- exp(test_data[[log_metal]])
  
  # 计算动态 y 轴范围 (深度聚焦：移除 0 限制，缩小缓冲至 5%)
  y_min <- min(test_data$lower, na.rm=T)
  y_max <- max(test_data$upper, na.rm=T)
  y_range <- y_max - y_min
  
  # E. 优化后的 ggplot 绘图 (超大字体 + 深度变焦版)
  p <- ggplot(test_data, aes(x = metal_orig, y = prob)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.15) +
    geom_line(color = "steelblue", linewidth = 2.5) + # 进一步加粗趋势线至 2.5
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) + # 提高百分比精度
    # 核心修改：允许 y 轴根据数据范围自由起始
    coord_cartesian(ylim = c(y_min - 0.05 * y_range, y_max + 0.05 * y_range)) + 
    labs(title = paste("金属:", metal), 
         subtitle = paste0("非线性P (D1) = ", sprintf("%.4f", p_nonlin)), # 提高P值精度
         x = paste(metal, "浓度 (μg/L)"), 
         y = "预测 LVR 概率 (%)") +
    theme_minimal(base_size = 48) + # 超大基础字体
    theme(
      plot.title = element_text(size = 60, face = "bold", hjust = 0.5, margin = margin(b = 20)), 
      plot.subtitle = element_text(size = 50, hjust = 0.5, color = "darkred", face = "italic"), 
      axis.title.x = element_text(size = 52, face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(size = 52, face = "bold", margin = margin(r = 20)),
      axis.text = element_text(size = 42, color = "black", face = "bold"), 
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 2), # 加粗边框
      plot.margin = margin(30, 30, 30, 30) # 增加留白防止文字切断
    )
  
  plot_list[[metal]] <- p
}

## ═══════════════════════════════════════════════════════
## 步骤 3: 结果保存
## ═══════════════════════════════════════════════════════

write.csv(summary_results, here("outputs", "EF", "A1_RCS_Summary.csv"), row.names = FALSE)

# 组合图：调整列数与间距
combined_plot <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(title = "图1 金属与LVR风险之间的剂量-反应关系",
                  theme = theme(plot.title = element_text(size = 64, face = "bold", hjust = 0.5)))

ggsave(here("plots", "EF", "A1_RCS_Summary.png"), 
       combined_plot, width = 16, height = 20, dpi = 300, bg = "white")

cat("\n✅ 优化版 RCS 分析完成！图表已保存。\n")