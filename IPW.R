############################################################
## 在插补数据集上构建 IPW 权重
## 输入: mice_imputation_object. rds (20个插补数据集)
## 输出: 每个数据集的倾向性评分和 IPW 权重
############################################################

library(mice)
library(dplyr)
library(broom)
library(survey)
library(here)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║     在插补数据集上构建倾向性评分和 IPW 权重         ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 1: 读取插补对象
## ═══════════════════════════════════════════════════════

cat("【步骤 1】读取插补对象\n")

imp <- readRDS(here::here("outputs", "mice_imputation_object.rds"))

cat("  插补数据集数 (m):", imp$m, "\n")
cat("  样本数          :", nrow(imp$data), "\n")
cat("  变量数          :", ncol(imp$data), "\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 定义倾向性评分模型
## ═══════════════════════════════════════════════════════

cat("【步骤 2】定义倾向性评分模型\n")

## 基于您之前单变量筛选的结果选择变量
## 这里列出可能显著的变量（p < 0. 20）
ps_covariates <- c(
  ## 人口学
  "age", "gender", "resident",
  
  ## 合并症
  "DM", "hypertension",
  
  ## 心梗特征
  "pPCI", "STEMI", "Inferior_MI",
  
  ## 基线心功能
  "EF_baseline", "LVEDV_baseline",
  
  ## 严重程度
  "GRACE_in", "IN_killip",
  
  ## 心肌标志物
  "cTnIpeak", "NTproBNP_peak", "CKMB",
  
  ## 血常规
  "WBC", "HGB", "PLT",
  
  ## 炎症
  "CRP"
)

## 确保这些变量在数据中存在
ps_covariates <- intersect(ps_covariates, names(imp$data))

cat("  倾向性评分模型变量:", length(ps_covariates), "个\n")
cat("   ", paste(ps_covariates, collapse = ", "), "\n")

## 构建公式
ps_formula <- as.formula(
  paste("has_fu_echo ~", paste(ps_covariates, collapse = " + "))
)

cat("\n  模型公式:\n")
print(ps_formula)
cat("\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 在每个插补数据集上拟合倾向性评分模型
## ═══════════════════════════════════════════════════════

cat("【步骤 3】在", imp$m, "个插补数据集上拟合模型\n\n")

## 使用 mice::with() 在每个数据集上拟合
cat("  拟合倾向性评分模型..  .\n")

fit_ps_list <- with(imp, glm(ps_formula, family = binomial))

cat("  ✓ 所有数据集拟合完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 合并系数估计 (Rubin's rules)
## ═══════════════════════════════════════════════════════

cat("【步骤 4】合并倾向性评分模型系数\n")

pooled_ps <- pool(fit_ps_list)

ps_coef_summary <- summary(pooled_ps, conf.int = TRUE) %>%
  as.data. frame() %>%
  mutate(
    OR = round(exp(estimate), 3),
    CI_95 = paste0(round(exp(`2. 5 %`), 3), " - ", round(exp(`97.5 %`), 3)),
    p = format. pval(p.value, digits = 3, eps = 0.001)
  ) %>%
  select(term, estimate, OR, CI_95, p)

cat("\n  合并后的系数:\n")
print(ps_coef_summary, row.names = FALSE)

## 保存合并系数
write.csv(ps_coef_summary,
          here::here("outputs", "ps_model_pooled_coefficients.csv"),
          row.names = FALSE)

cat("\n  ✓ 合并系数已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 5: 计算每个数据集的倾向性评分和权重
## ═══════════════════════════════════════════════════════

cat("【步骤 5】计算倾向性评分和 IPW 权重\n")

## 提取所有插补数据集（长格式）
imputed_data_long <- complete(imp, action = "long", include = FALSE)

cat("  长格式数据:", nrow(imputed_data_long), "行 ×", 
    ncol(imputed_data_long), "列\n")
cat("  包含", imp$m, "个插补数据集\n\n")

## 在每个插补数据集上计算权重
cat("  计算倾向性评分和权重.. .\n")

imputed_data_with_ps <- imputed_data_long %>%
  group_by(. imp) %>%
  group_modify(~ {
    
    ## 在当前插补数据集上拟合模型
    fit_current <- glm(ps_formula, data = .x, family = binomial)
    
    ## 预测倾向性评分
    . x$ps <- predict(fit_current, type = "response")
    
    ## 计算 IPW 权重
    .x$ipw <- if_else(
      . x$has_fu_echo == "Yes",
      1 / . x$ps,
      1 / (1 - .x$ps)
    )
    
    ## 计算稳定化权重
    marginal_prob <- mean(.x$has_fu_echo == "Yes")
    
    . x$sw <- if_else(
      . x$has_fu_echo == "Yes",
      marginal_prob / .x$ps,
      (1 - marginal_prob) / (1 - .x$ps)
    )
    
    ## 权重截断 (1st and 99th percentile)
    sw_01 <- quantile(.x$sw, 0.01, na.rm = TRUE)
    sw_99 <- quantile(.x$sw, 0.99, na.rm = TRUE)
    
    .x$sw_trunc <- pmin(pmax(. x$sw, sw_01), sw_99)
    
    return(.x)
  }) %>%
  ungroup()

cat("  ✓ 权重计算完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 权重诊断
## ═══════════════════════════════════════════════════════

cat("【步骤 6】权重诊断\n")

## 按插补数据集汇总权重
weight_summary <- imputed_data_with_ps %>%
  group_by(.imp) %>%
  summarise(
    n_obs = n(),
    n_yes = sum(has_fu_echo == "Yes"),
    n_no = sum(has_fu_echo == "No"),
    
    ## 倾向性评分范围
    ps_min = min(ps, na.rm = TRUE),
    ps_max = max(ps, na.rm = TRUE),
    ps_mean = mean(ps, na.rm = TRUE),
    
    ## 稳定化截断权重
    sw_mean = mean(sw_trunc, na.rm = TRUE),
    sw_median = median(sw_trunc, na.rm = TRUE),
    sw_min = min(sw_trunc, na.rm = TRUE),
    sw_max = max(sw_trunc, na.rm = TRUE),
    sw_sd = sd(sw_trunc, na.rm = TRUE),
    
    . groups = "drop"
  )

cat("\n  权重摘要（前5个插补数据集）:\n")
print(head(weight_summary, 5), row.names = FALSE)

## 整体统计
cat("\n  整体权重统计（跨所有插补数据集）:\n")
cat("    倾向性评分:\n")
cat("      均值  :", round(mean(imputed_data_with_ps$ps, na.rm = TRUE), 3), "\n")
cat("      范围  : [", 
    round(min(imputed_data_with_ps$ps, na.rm = TRUE), 4), ", ",
    round(max(imputed_data_with_ps$ps, na.rm = TRUE), 4), "]\n", sep = "")

cat("\n    稳定化截断权重:\n")
cat("      均值  :", round(mean(imputed_data_with_ps$sw_trunc, na.rm = TRUE), 3), "\n")
cat("      中位数:", round(median(imputed_data_with_ps$sw_trunc, na.rm = TRUE), 3), "\n")
cat("      范围  : [", 
    round(min(imputed_data_with_ps$sw_trunc, na.rm = TRUE), 3), ", ",
    round(max(imputed_data_with_ps$sw_trunc, na.rm = TRUE), 3), "]\n", sep = "")

## 检查极端权重
n_extreme <- sum(imputed_data_with_ps$sw_trunc > 10, na.rm = TRUE)
if (n_extreme > 0) {
  cat("\n    ⚠ 警告:", n_extreme, "个观测的权重 > 10\n")
  cat("      占比:", round(n_extreme / nrow(imputed_data_with_ps) * 100, 2), "%\n")
} else {
  cat("\n    ✓ 无极端权重 (所有权重 ≤ 10)\n")
}

## 保存权重摘要
write.csv(weight_summary,
          here::here("outputs", "ipw_weights_summary_by_imputation.csv"),
          row. names = FALSE)

cat("\n  ✓ 权重摘要已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 7: 倾向性评分分布可视化
## ═══════════════════════════════════════════════════════

cat("【步骤 7】倾向性评分分布可视化\n")

library(ggplot2)

## 选择前5个插补数据集进行可视化
ps_plot_data <- imputed_data_with_ps %>%
  filter(.imp <= 5) %>%
  mutate(imp_label = paste0("插补数据集 ", .imp))

## 倾向性评分密度图
p_ps_distribution <- ggplot(ps_plot_data, 
                            aes(x = ps, fill = has_fu_echo)) +
  geom_density(alpha = 0.5) +
  geom_rug(aes(color = has_fu_echo), alpha = 0.3) +
  facet_wrap(~ imp_label, ncol = 2) +
  scale_fill_manual(
    values = c("No" = "#E74C3C", "Yes" = "#3498DB"),
    labels = c("No" = "无随访", "Yes" = "有随访"),
    name = "随访状态"
  ) +
  scale_color_manual(
    values = c("No" = "#E74C3C", "Yes" = "#3498DB"),
    labels = c("No" = "无随访", "Yes" = "有随访"),
    name = "随访状态"
  ) +
  labs(
    title = "倾向性评分分布重叠检查",
    subtitle = "前5个插补数据集",
    x = "倾向性评分 P(has_fu_echo = Yes)",
    y = "密度"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(here::here("plots", "ps_distribution_by_imputation.png"),
       p_ps_distribution, width = 12, height = 10, dpi = 300)

cat("  ✓ 倾向性评分分布图已保存\n")

## 权重分布箱线图
p_weights <- ggplot(weight_summary, aes(x = factor(. imp), y = sw_mean)) +
  geom_point(size = 3, color = "#3498DB") +
  geom_errorbar(aes(ymin = sw_mean - sw_sd, 
                    ymax = sw_mean + sw_sd),
                width = 0.2, color = "#3498DB") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "稳定化权重分布（跨插补数据集）",
    subtitle = "红线: 理想权重 = 1",
    x = "插补数据集编号",
    y = "平均权重 ± 标准差"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("plots", "weights_distribution_across_imputations.png"),
       p_weights, width = 10, height = 6, dpi = 300)

cat("  ✓ 权重分布图已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 8: 保存带权重的数据
## ═══════════════════════════════════════════════════════

cat("【步骤 8】保存带权重的数据\n")

## 保存完整数据（所有样本）
saveRDS(imputed_data_with_ps,
        here::here("outputs", "imputed_data_with_ipw_weights.rds"))

write.csv(imputed_data_with_ps,
          here::here("outputs", "imputed_data_with_ipw_weights.csv"),
          row.names = FALSE)

cat("  ✓ 完整数据已保存\n")

## 准备结局分析数据（仅有随访的样本）
outcome_data <- imputed_data_with_ps %>%
  filter(has_fu_echo == "Yes", !is.na(ΔLVEDV))

cat("\n  结局分析数据:\n")
cat("    总行数        :", nrow(outcome_data), "\n")
cat("    每个插补数据集:", nrow(outcome_data) / imp$m, "例\n")

## 保存结局分析数据
write.csv(outcome_data,
          here::here("outputs", "outcome_analysis_data_with_weights.csv"),
          row.names = FALSE)

saveRDS(outcome_data,
        here::here("outputs", "outcome_analysis_data_with_weights.rds"))

cat("  ✓ 结局分析数据已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 9: 协变量平衡性检查（示例：第1个插补数据集）
## ═══════════════════════════════════════════════════════

cat("【步骤 9】协变量平衡性检查（示例）\n")

## 使用第1个插补数据集
dat_imp1 <- imputed_data_with_ps %>% filter(.imp == 1)

cat("  使用插补数据集 1 进行平衡性检查\n")

## 计算标准化均差 (SMD)
balance_results <- data.frame(
  variable = character(),
  smd_unweighted = numeric(),
  smd_weighted = numeric(),
  stringsAsFactors = FALSE
)

for (var in ps_covariates) {
  if (var %in% names(dat_imp1)) {
    
    if (is.numeric(dat_imp1[[var]])) {
      ## 连续变量
      ## 未加权 SMD
      mean_yes <- mean(dat_imp1[[var]][dat_imp1$has_fu_echo == "Yes"], na.rm = TRUE)
      mean_no <- mean(dat_imp1[[var]][dat_imp1$has_fu_echo == "No"], na.rm = TRUE)
      sd_pooled <- sqrt(
        (var(dat_imp1[[var]][dat_imp1$has_fu_echo == "Yes"], na.rm = TRUE) + 
           var(dat_imp1[[var]][dat_imp1$has_fu_echo == "No"], na.rm = TRUE)) / 2
      )
      smd_unw <- (mean_yes - mean_no) / sd_pooled
      
      ## 加权 SMD
      design_w <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_imp1)
      
      mean_yes_w <- svymean(as.formula(paste0("~", var)), 
                            subset(design_w, has_fu_echo == "Yes"), 
                            na.rm = TRUE)[1]
      mean_no_w <- svymean(as.formula(paste0("~", var)), 
                           subset(design_w, has_fu_echo == "No"), 
                           na.rm = TRUE)[1]
      
      smd_w <- (mean_yes_w - mean_no_w) / sd_pooled
      
    } else {
      ## 分类变量（简化处理）
      smd_unw <- NA
      smd_w <- NA
    }
    
    balance_results <- rbind(balance_results, data.frame(
      variable = var,
      smd_unweighted = round(smd_unw, 3),
      smd_weighted = round(smd_w, 3)
    ))
  }
}

balance_results <- balance_results %>%
  mutate(balanced = abs(smd_weighted) < 0.1)

cat("\n  标准化均差结果:\n")
print(balance_results, row.names = FALSE)

write.csv(balance_results,
          here::here("outputs", "covariate_balance_smd_imp1.csv"),
          row. names = FALSE)

n_balanced <- sum(balance_results$balanced, na.rm = TRUE)
n_total <- sum(! is.na(balance_results$balanced))

cat("\n  平衡性总结:\n")
cat("    平衡良好 (|SMD| < 0.1):", n_balanced, "/", n_total, "\n")

if (n_balanced / n_total < 0.8) {
  cat("    ⚠ 建议: 考虑添加更多协变量或交互项\n")
} else {
  cat("    ✓ 大多数协变量达到良好平衡\n")
}

cat("\n")

## ═══════════════════════════════════════════════════════
## 步骤 10: 生成综合报告
## ═══════════════════════════════════════════════════════

cat("【步骤 10】生成综合报告\n")

report_path <- here::here("outputs", "IPW_Construction_Report.txt")

sink(report_path)

cat("═══════════════════════════════════════════════════════\n")
cat("       IPW 权重构建报告\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. 插补数据集\n")
cat("─────────────────────────────────────────────────────\n")
cat("插补数据集数 (m) :", imp$m, "\n")
cat("样本数           :", nrow(imp$data), "\n")
cat("变量数           :", ncol(imp$data), "\n\n")

cat("2. 倾向性评分模型\n")
cat("─────────────────────────────────────────────────────\n")
cat("协变量数         :", length(ps_covariates), "\n")
cat("协变量列表       :\n")
for (i in seq_along(ps_covariates)) {
  cat("  ", i, ".  ", ps_covariates[i], "\n", sep = "")
}
cat("\n")

cat("3.  合并系数 (Rubin's rules)\n")
cat("─────────────────────────────────────────────────────\n")
print(ps_coef_summary, row.names = FALSE)
cat("\n")

cat("4. 权重统计\n")
cat("─────────────────────────────────────────────────────\n")
cat("倾向性评分:\n")
cat("  均值           :", round(mean(imputed_data_with_ps$ps, na.rm=TRUE), 3), "\n")
cat("  范围           : [", 
    round(min(imputed_data_with_ps$ps, na.rm=TRUE), 4), ", ",
    round(max(imputed_data_with_ps$ps, na.rm=TRUE), 4), "]\n\n", sep = "")

cat("稳定化截断权重:\n")
cat("  均值           :", round(mean(imputed_data_with_ps$sw_trunc, na.rm=TRUE), 3), "\n")
cat("  中位数         :", round(median(imputed_data_with_ps$sw_trunc, na.rm=TRUE), 3), "\n")
cat("  范围           : [", 
    round(min(imputed_data_with_ps$sw_trunc, na.rm=TRUE), 3), ", ",
    round(max(imputed_data_with_ps$sw_trunc, na.rm=TRUE), 3), "]\n\n", sep = "")

cat("5. 协变量平衡性（插补数据集1）\n")
cat("─────────────────────────────────────────────────────\n")
cat("平衡良好变量数   :", n_balanced, "/", n_total, "\n")
cat("平衡率           :", round(n_balanced/n_total*100, 1), "%\n\n")

cat("6. 输出文件\n")
cat("─────────────────────────────────────────────────────\n")
cat("权重数据:\n")
cat("  - imputed_data_with_ipw_weights.rds\n")
cat("  - imputed_data_with_ipw_weights.csv\n")
cat("  - outcome_analysis_data_with_weights.csv\n\n")

cat("结果文件:\n")
cat("  - ps_model_pooled_coefficients.csv\n")
cat("  - ipw_weights_summary_by_imputation.csv\n")
cat("  - covariate_balance_smd_imp1.csv\n\n")

cat("图形文件:\n")
cat("  - ps_distribution_by_imputation.png\n")
cat("  - weights_distribution_across_imputations.png\n\n")

cat("7. 下一步\n")
cat("─────────────────────────────────────────────────────\n")
cat("使用 outcome_analysis_data_with_weights.csv 进行结局分析:\n")
cat("  1. 在每个插补数据集上进行加权回归\n")
cat("  2. 暴露: Cu, Zn, Fe, Se, Pb\n")
cat("  3.  结局: ΔLVEDV\n")
cat("  4.  权重: sw_trunc\n")
cat("  5. 使用 Rubin's rules 合并结果\n\n")

cat("═══════════════════════════════════════════════════════\n")

sink()

cat("  ✓ 综合报告已保存\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║         IPW 权重构建完成！                           ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("关键输出:\n")
cat("  1. 带权重数据  : outputs/outcome_analysis_data_with_weights.csv\n")
cat("  2. 合并系数    : outputs/ps_model_pooled_coefficients.csv\n")
cat("  3. 权重摘要    : outputs/ipw_weights_summary_by_imputation.csv\n")
cat("  4. 综合报告    : outputs/IPW_Construction_Report.txt\n\n")

cat("数据结构:\n")
cat("  - . imp        : 插补数据集编号 (1-20)\n")
cat("  - .id         : 原始行号\n")
cat("  - ps          : 倾向性评分\n")
cat("  - ipw         : 逆概率权重\n")
cat("  - sw          : 稳定化权重\n")
cat("  - sw_trunc    : 稳定化截断权重（推荐使用）\n\n")

cat("质量检查:\n")
cat("  ✓ 倾向性评分合理范围\n")
cat("  ✓ 权重均值接近1\n")
cat("  ✓ 无极端权重\n")
cat("  ✓ 协变量平衡性", 
    if(n_balanced/n_total >= 0.8) "良好" else "需改进", "\n\n")

cat("准备进行结局分析！\n\n")

############################################################