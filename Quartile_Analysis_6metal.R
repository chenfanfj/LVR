############################################################
## 主要金属四分位数分析 
############################################################

rm(list = ls())
library(dplyr)
library(survey)
library(mitools)
library(flextable)
library(ggplot2)
library(here)
library(officer)

# 1. 数据加载与对齐
outcome_data <- readRDS(here::here("outputs", "4 IPW2_complete", "outcome_analysis_data_with_weights_extended.rds"))
outcome_data_filtered <- outcome_data %>%
  filter(EF_baseline >= 50) %>%
  mutate(LVR = ifelse(EF_fu < 50, 1, 0))

# 挑选 6 种核心金属
target_metals <- c("Fe", "Li", "Sr", "Zn", "Pb", "V")
covariates_all <- c("age", "gender", "EF_baseline", "LVEDV_baseline", "cTnIpeak", 
                    "pPCI", "STEMI", "smoking", "DM", "hypertension", 
                    "NTproBNP_peak", "GRACE_in", "WBC", "HGB", "CRP", 
                    "CHOL", "LDL", "AST")

# 建立映射 (键为代码变量名, 值为显示名称)
metal_map <- c("Fe"="铁 (Fe)", "Li"="锂 (Li)", "Sr"="锶 (Sr)", 
               "Zn"="锌 (Zn)", "Pb"="铅 (Pb)", "V"="钒 (V)")

# 2. 循环分析生成结果
all_results_list <- list()
for(metal in target_metals) {
  log_metal_name <- paste0("log_", metal)
  cat("\n正在处理金属:", metal)
  
  fits_q <- list(); fits_trend <- list()
  n_imps <- length(unique(outcome_data_filtered$.imp))
  
  for(i in 1:n_imps) {
    dat_i <- outcome_data_filtered %>% filter(.imp == i)
    
    # 计算四分位数
    qs <- quantile(dat_i[[log_metal_name]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    dat_i$Quartile <- cut(dat_i[[log_metal_name]], breaks = qs, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
    dat_i$Q_trend <- as.numeric(dat_i$Quartile)
    
    # 重新计算 PS 评分 (确保模型收敛)
    ps_form <- as.formula(paste0(log_metal_name, " ~ ", paste(covariates_all, collapse = " + ")))
    dat_i$ps_score <- predict(lm(ps_form, data = dat_i))
    
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    fits_q[[i]] <- svyglm(LVR ~ Quartile + ps_score, design = design_i, family = quasibinomial())
    fits_trend[[i]] <- svyglm(LVR ~ Q_trend + ps_score, design = design_i, family = quasibinomial())
  }
  
  # 合并结果
  pooled_q <- MIcombine(fits_q)
  ci_mat <- confint(pooled_q)
  betas <- pooled_q$coefficients
  p_vals <- 2 * (1 - pnorm(abs(betas / sqrt(diag(pooled_q$variance)))))
  q_rows <- c("QuartileQ2", "QuartileQ3", "QuartileQ4")
  
  pooled_t <- MIcombine(fits_trend)
  p_trend_val <- 2 * (1 - pnorm(abs(pooled_t$coefficients["Q_trend"] / sqrt(pooled_t$variance["Q_trend", "Q_trend"]))))
  
  # 存储结果数据
  all_results_list[[metal]] <- data.frame(
    Metal_ID = metal,
    Metal_Label = unname(metal_map[metal]),
    Group = c("Q1 (Reference)", "Q2", "Q3", "Q4"),
    OR = unname(c(1.00, exp(betas[q_rows]))),
    Lower = unname(c(1.00, exp(ci_mat[q_rows, 1]))),
    Upper = unname(c(1.00, exp(ci_mat[q_rows, 2]))),
    P_val = unname(c("-", format.pval(p_vals[q_rows], digits = 3))),
    P_trend = unname(c(format.pval(p_trend_val, digits = 3), "", "", ""))
  )
}

# 3. 准备绘图数据 (定义 plot_data 并修复铁曲线)
final_df <- bind_rows(all_results_list)

plot_data <- final_df %>%
  mutate(
    Group_Factor = factor(Group, levels = c("Q1 (Reference)", "Q2", "Q3", "Q4")),
    # 核心步骤：截断铁(Fe)等极端置信区间上限至 5，使 Y 轴比例恢复正常
    Upper_Capped = ifelse(Upper > 5, 5, Upper),
    Lower_Capped = ifelse(Lower < 0, 0, Lower)
  )

# 4. 可视化：大字体组合图
p_main <- ggplot(plot_data, aes(x = Group_Factor, y = OR, color = Metal_Label, group = Metal_Label)) +
  geom_line(linewidth = 1.8) + 
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = Lower_Capped, ymax = Upper_Capped), width = 0.2, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", linewidth = 1) +
  facet_wrap(~Metal_Label, scales = "free_y", ncol = 3) +
  
  # 样式美化与字体放大
  theme_minimal(base_size = 36) + 
  labs(
    title = "核心血清金属暴露与 LVR 风险的趋势关联",
    subtitle = "注：Q1为参考组；铁(Fe)等金属的极端置信区间上限已截断以优化视觉展示",
    x = "暴露水平 (Quartiles)", 
    y = "Odds Ratio (95% CI)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 36, face = "bold", margin = margin(b=15)),
    plot.subtitle = element_text(hjust = 0.5, size = 28, color = "grey40"),
    axis.title = element_text(face = "bold", size = 28),
    axis.text = element_text(color = "black", size = 24),
    strip.background = element_rect(fill = "#ECEFF1", color = NA),
    strip.text = element_text(size = 28, face = "bold"),
    panel.spacing = unit(2, "lines"),
    panel.grid.minor = element_blank()
  )

# 保存高清晰度图片
ggsave(here::here("plots", "EF", "Main_Figure_6_Metals.png"), 
       p_main, width = 15, height = 11, dpi = 300)

cat("\n✅ 分析与绘图已完成。")