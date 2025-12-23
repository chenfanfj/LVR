############################################################
##  A3.1: 多金属同时纳入logistic模型
##  输入:  imputed_data_with_ipw_weights. rds
##  输出: 多金属模型系数表 + 森林图
############################################################

rm(list = ls())
gc()

library(mice)
library(survey)
library(ggplot2)
library(dplyr)
library(car)  # VIF检查
library(here)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║       A3.1: 多金属联合模型 (同时纳入)              ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤1: 读取数据
## ═══════════════════════════════════════════════════════

cat("【步骤1】读取数据\n")

imputed_data <- readRDS(here("outputs", "4 IPW2_complete", "imputed_data_with_ipw_weights_extended.rds"))

outcome_data <- imputed_data %>% 
  filter(has_fu_echo == "Yes", !is.na(LVEDV_fu)) %>%
  mutate(LVR = as.numeric(EF_baseline >= 50 & EF_fu < 50))

cat("  样本数:", nrow(outcome_data), "\n\n")

## ═══════════════════════════════════════════════════════
## 步骤2: 定义金属与协变量
## ═══════════════════════════════════════════════════════

# 核心金属 (文献相关 + 本研究显著)
metals_core <- c("Pb", "Zn", "Fe", "Cu", "Se", "As")
log_metals <- paste0("log_", metals_core)

# 协变量
covariates <- c("age", "gender", "EF_baseline", "LVEDV_baseline", 
                "cTnIpeak", "pPCI", "STEMI", "smoking", "DM", 
                "hypertension", "NTproBNP_peak", "GRACE_in", 
                "WBC", "HGB", "CRP", "CHOL", "LDL", "AST")

## ═══════════════════════════════════════════════════════
## 步骤3: 检查多重共线性 (VIF)
## ═══════════════════════════════════════════════════════

cat("【步骤3】检查金属间多重共线性\n")

# 使用第1个插补
dat_1 <- outcome_data %>% filter(.imp == 1)

# 移除缺失
dat_1_complete <- dat_1 %>%
  select(LVR, all_of(log_metals), all_of(covariates), sw_trunc) %>%
  na.omit()

# 拟合完整模型
formula_full <- as.formula(paste(
  "LVR ~", paste(c(log_metals, covariates), collapse = " + ")
))

svy_design_1 <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_1_complete)
fit_full <- svyglm(formula_full, design = svy_design_1, family = quasibinomial)

# 计算VIF (使用car包)
# 注意:  svyglm对象可能不支持vif()，需提取系数矩阵手动计算
# 简化方法: 用未加权的glm计算VIF (仅作参考)
fit_unweighted <- glm(formula_full, data = dat_1_complete, family = binomial)
vif_values <- vif(fit_unweighted)

cat("\n  金属VIF值:\n")
print(vif_values[log_metals])

# 判断共线性
high_vif <- names(vif_values[vif_values > 5])
if (length(high_vif) > 0) {
  cat("\n  ⚠️ 高VIF变量 (>5):", paste(high_vif, collapse = ", "), "\n")
  cat("  建议:  考虑使用正则化方法或删除高VIF金属\n\n")
} else {
  cat("\n  ✓ 无严重共线性问题\n\n")
}

## ═══════════════════════════════════════════════════════
## 步骤4: 对20个插补拟合多金属模型
## ═══════════════════════════════════════════════════════

cat("【步骤4】拟合多金属模型\n")

multi_metal_results <- list()

for (i in 1:20) {
  
  cat("  插补", i, "\n")
  
  dat_i <- outcome_data %>% filter(.imp == i)
  
  # 移除缺失
  dat_i_complete <- dat_i %>%
    select(LVR, all_of(log_metals), all_of(covariates), sw_trunc) %>%
    na.omit()
  
  # 加权设计
  svy_design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i_complete)
  
  # 拟合
  fit_i <- svyglm(formula_full, design = svy_design_i, family = quasibinomial)
  
  # 提取金属系数
  coef_metals <- coef(fit_i)[log_metals]
  se_metals <- summary(fit_i)$coefficients[log_metals, "Std. Error"]
  
  multi_metal_results[[i]] <- data.frame(
    .imp = i,
    Metal = metals_core,
    log_OR = coef_metals,
    SE = se_metals
  )
}

## ═══════════════════════════════════════════════════════
## 步骤5: Rubin规则合并
## ═══════════════════════════════════════════════════════

cat("\n【步骤5】合并结果 (Rubin规则)\n")

multi_metal_pooled <- bind_rows(multi_metal_results) %>%
  group_by(Metal) %>%
  summarise(
    log_OR_pooled = mean(log_OR, na.rm = TRUE),
    within_var = mean(SE^2, na.rm = TRUE),
    between_var = var(log_OR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_var = within_var + (1 + 1/20) * between_var,
    SE_pooled = sqrt(total_var),
    OR_pooled = exp(log_OR_pooled),
    OR_lower = exp(log_OR_pooled - 1.96 * SE_pooled),
    OR_upper = exp(log_OR_pooled + 1.96 * SE_pooled),
    Z = log_OR_pooled / SE_pooled,
    P_value = 2 * (1 - pnorm(abs(Z)))
  )

cat("\n多金属模型结果:\n")
print(multi_metal_pooled %>% select(Metal, OR_pooled, OR_lower, OR_upper, P_value), row.names = FALSE)

## ═══════════════════════════════════════════════════════
## 步骤6: 保存结果
## ═══════════════════════════════════════════════════════

write.csv(multi_metal_pooled,
          here("outputs", "A3_Multi_Metal_Model_Results.csv"),
          row.names = FALSE)

cat("\n  ✓ 结果已保存: outputs/EF/A3_Multi_Metal_Model_Results.csv\n")


## ═══════════════════════════════════════════════════════
## 步骤7: 森林图
## ═══════════════════════════════════════════════════════
# 根据 A3_Multi_Metal_Model_Results.csv 实际数值设置截断值
# 最小值 0.16，最大值 32.42。为了兼顾 Fe/Zn (0.16) 和 Pb (3.85) 的展示：
lower_limit_mm <- 0.1 
upper_limit_mm <- 5.0 # Cu (32.42) 和 Se (11.58) 将被截断显示为箭头

# 1. 准备绘图表格数据
plot_data_mm <- multi_metal_pooled %>%
  arrange(desc(OR_pooled)) %>%
  mutate(
    OR_fmt = sprintf("%.2f (%.2f-%.2f)", OR_pooled, OR_lower, OR_upper),
    P_fmt  = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value)),
    is_summary = FALSE
  )

header_mm <- data.frame(
  Metal = "金属",
  OR_fmt = "OR (95% CI)",
  P_fmt = "P值",
  OR_pooled = as.numeric(NA), 
  OR_lower = as.numeric(NA), 
  OR_upper = as.numeric(NA),
  is_summary = TRUE
)

final_plot_df <- bind_rows(header_mm, plot_data_mm)

table_text_mm <- list(
  final_plot_df$Metal,
  final_plot_df$OR_fmt,
  final_plot_df$P_fmt
)

# 2. 执行绘图
output_path_mm <- here("plots", "EF", "A3_Multi_Metal_Forest_Plot.png")

png(filename = output_path_mm, width = 2800, height = 1800, res = 300)

p_mm <- forestplot(
  labeltext  = table_text_mm,
  mean       = final_plot_df$OR_pooled,
  lower      = final_plot_df$OR_lower,
  upper      = final_plot_df$OR_upper,
  is.summary = final_plot_df$is_summary,
  
  graphwidth = unit(70, "mm"), # 进一步收窄中心图以突出大文字
  graph.pos  = 2, 
  
  xlog       = TRUE,
  zero       = 1.0, 
  # 调整刻度，使其在 0.1-5.0 之间分布均匀
  xticks     = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
  clip       = c(lower_limit_mm, upper_limit_mm),
  
  boxsize    = 0.20,
  lineheight = unit(15, 'mm'), 
  colgap     = unit(12, 'mm'), 
  lwd.ci     = 2.2,
  col        = fpColors(box = '#458B00', lines = 'black', zero = '#7AC5CD'),
  
  xlab       = "Odds Ratio (OR) on Log Scale",
  title      = "多金属联合模型：与LVR风险的关联",
  
  txt_gp     = fpTxtGp(
    label   = gpar(cex = 2.2), # 大文字
    summary = gpar(cex = 2.2, fontface = "bold"),
    ticks   = gpar(cex = 1.8),
    xlab    = gpar(cex = 2.2),
    title   = gpar(cex = 3.0, fontface = "bold")
  )
)

print(p_mm)

# 底部说明
grid.text("注: 所有金属同时纳入模型 (互相调整); 调整18个协变量; 坐标轴根据数据分布截断 (0.1-5.0); * P < 0.05", 
          x = unit(0.05, "npc"), y = unit(0.04, "npc"), 
          just = "left", gp = gpar(cex = 2.0, fontface = "italic"))

dev.off()

cat("  ✓ 多金属森林图已成功保存:", output_path_mm, "\n")

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║            A3.1 多金属模型完成！                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")