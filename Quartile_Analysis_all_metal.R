############################################################
## 最终整合版：金属四分位数分析 (解决报错、大字体、铁曲线优化)
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

target_metals <- c("Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
                   "Mo", "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li")
covariates_all <- c("age", "gender", "EF_baseline", "LVEDV_baseline", "cTnIpeak", 
                    "pPCI", "STEMI", "smoking", "DM", "hypertension", 
                    "NTproBNP_peak", "GRACE_in", "WBC", "HGB", "CRP", 
                    "CHOL", "LDL", "AST")
metal_map <- c("Cu"="铜 (Cu)", "Zn"="锌 (Zn)", "Fe"="铁 (Fe)", "Se"="硒 (Se)", "Pb"="铅 (Pb)", 
               "Al"="铝 (Al)", "As"="砷 (As)", "Cr"="铬 (Cr)", "Mn"="锰 (Mn)", "Ni"="镍 (Ni)", 
               "Mo"="钼 (Mo)", "Rb"="铷 (Rb)", "Sb"="锑 (Sb)", "Sn"="锡 (Sn)", "Sr"="锶 (Sr)", 
               "V"="钒 (V)", "Ba"="钡 (Ba)", "B"="硼 (B)", "Li"="锂 (Li)")

# 2. 循环分析生成结果
all_results_list <- list()
for(metal in target_metals) {
  log_metal_name <- paste0("log_", metal)
  fits_q <- list(); fits_trend <- list()
  n_imps <- length(unique(outcome_data_filtered$.imp))
  
  for(i in 1:n_imps) {
    dat_i <- outcome_data_filtered %>% filter(.imp == i)
    qs <- quantile(dat_i[[log_metal_name]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    dat_i$Quartile <- cut(dat_i[[log_metal_name]], breaks = qs, include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
    dat_i$Q_trend <- as.numeric(dat_i$Quartile)
    ps_form <- as.formula(paste0(log_metal_name, " ~ ", paste(covariates_all, collapse = " + ")))
    dat_i$ps_score <- predict(lm(ps_form, data = dat_i))
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    fits_q[[i]] <- svyglm(LVR ~ Quartile + ps_score, design = design_i, family = quasibinomial())
    fits_trend[[i]] <- svyglm(LVR ~ Q_trend + ps_score, design = design_i, family = quasibinomial())
  }
  
  pooled_q <- MIcombine(fits_q)
  ci_mat <- confint(pooled_q); betas <- pooled_q$coefficients
  ses <- sqrt(diag(pooled_q$variance)); p_vals <- 2 * (1 - pnorm(abs(betas / ses)))
  q_rows <- c("QuartileQ2", "QuartileQ3", "QuartileQ4")
  pooled_t <- MIcombine(fits_trend)
  p_trend_val <- 2 * (1 - pnorm(abs(pooled_t$coefficients["Q_trend"] / sqrt(pooled_t$variance["Q_trend", "Q_trend"]))))
  
  # 存储数据 (unname 防止警告)
  all_results_list[[metal]] <- data.frame(
    Metal_Name = unname(metal_map[metal]),
    Group = c("Q1 (Reference)", "Q2", "Q3", "Q4"),
    OR = unname(c(1.00, exp(betas[q_rows]))),
    Lower = unname(c(1.00, exp(ci_mat[q_rows, 1]))),
    Upper = unname(c(1.00, exp(ci_mat[q_rows, 2]))),
    P_val = unname(c("-", format.pval(p_vals[q_rows], digits = 3))),
    P_trend = unname(c(format.pval(p_trend_val, digits = 3), "", "", ""))
  )
}

# 3. 汇总数据
final_df <- bind_rows(all_results_list)
## ═══════════════════════════════════════════════════════
## 2. 导出 Word 三线表 
## ═══════════════════════════════════════════════════════

final_df_table <- final_df %>%
  mutate(OR_CI = sprintf("%.2f (%.2f, %.2f)", OR, Lower, Upper)) %>%
  mutate(OR_CI = ifelse(Group == "Q1 (Reference)", "1.00", OR_CI)) %>%
  select(Metal_Name, Group, OR_CI, P_val, P_trend) %>%
  rename(金属 = Metal_Name, 暴露分组 = Group, `OR (95% CI)` = OR_CI, P值 = P_val, 趋势P值 = P_trend)

ft <- flextable(final_df_table) %>%
  add_header_lines("表. 血清金属浓度四分位数与 LVR 风险的 Logistic 回归关联") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 2)) %>%
  hline_bottom(part = "all", border = fp_border(width = 2)) %>%
  merge_v(j = c("金属", "趋势P值")) %>%
  align(align = "center", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

save_as_docx(ft, path = here::here("outputs", "EF", "LVR_Quartile_Table_all_metal.docx"))

## ═══════════════════════════════════════════════════════
## 3. 可视化趋势图 
## ═══════════════════════════════════════════════════════

# --- 核心修改：准备绘图数据 (定义 plot_data 并解决 Fe 的直线问题) ---
plot_data <- final_df %>%
  mutate(
    Group_Factor = factor(Group, levels = c("Q1 (Reference)", "Q2", "Q3", "Q4")),
    # 截断过大的置信区间上限至 5，确保 Y 轴刻度协调
    Upper_Capped = ifelse(Upper > 5, 5, Upper),
    Lower_Capped = ifelse(Lower < 0, 0, Lower)
  )

# 4. 可视化 (大字体 + 协调比例)
p_trend <- ggplot(plot_data, aes(x = Group_Factor, y = OR, color = Metal_Name, group = Metal_Name)) +
  # 加粗线条和点
  geom_line(linewidth = 1.5) + 
  geom_point(size = 3) +
  # 绘制截断后的误差线
  geom_errorbar(aes(ymin = Lower_Capped, ymax = Upper_Capped), width = 0.2, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", linewidth = 1) +
  facet_wrap(~Metal_Name, scales = "free_y") +
  
  # --- 极致大字体设置 ---
  theme_minimal(base_size = 28) + 
  labs(
    title = "金属浓度四分位数与 LVR 风险的趋势关联",
    subtitle = "注：铁 (Fe) 的上限已做截断处理以看清趋势",
    x = "暴露分组 (Quartiles)", y = "Odds Ratio (95% CI)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 44, face = "bold", margin = margin(b=15)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 24),
    strip.background = element_rect(fill = "#ECEFF1"),
    strip.text = element_text(size = 28, face = "bold"),
    panel.spacing = unit(1, "lines") # 增加子图间距
  )

# 保存
ggsave(here::here("plots", "EF", "Metal_Quartile_all_metal.png"), p_trend, width = 12, height = 9, dpi = 300)

cat("✅ 分析与绘图已完成，请查看 plots/EF 目录。\n")
