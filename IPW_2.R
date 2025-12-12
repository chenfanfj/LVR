############################################################
## IPW权重构建脚本 (最终修正版)
## 输入:  mice_imputation_final.rds
## 输出: 带IPW权重的数据 + SMD图
############################################################

rm(list = ls())
gc()

library(mice)
library(dplyr)
library(broom)
library(survey)
library(here)
library(ggplot2)
library(tidyr)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║     IPW权重构建 + 协变量平衡诊断 (最终版)          ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 1: 读取插补对象
## ═══════════════════════════════════════════════════════

cat("【步骤 1】读取插补对象\n")

imp <- readRDS(here::here("outputs", "mice_imputation_final.rds"))

cat("  插补数据集数 (m):", imp$m, "\n")
cat("  样本数          :", nrow(imp$data), "\n")
cat("  变量数          :", ncol(imp$data), "\n")

# 验证
check_data <- complete(imp, 1)
required_vars <- c("has_fu_echo", "LVEDV_fu", "LVEDV_baseline")

if(! all(required_vars %in% names(check_data))) {
  stop("❌ 缺少关键变量")
}

cat("  ✓ 验证通过\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 定义IPW协变量
## ═══════════════════════════════════════════════════════

cat("【步骤 2】定义IPW协变量\n")

ps_covariates <- c(
  "age", "gender", "resident", "DM", "hypertension",
  "pPCI", "STEMI", "EF_baseline", "LVEDV_baseline",
  "GRACE_in", "IN_killip", "cTnIpeak", "NTproBNP_peak", "CKMB",
  "WBC", "HGB", "PLT", "CRP"
)

ps_covariates <- intersect(ps_covariates, names(check_data))

cat("  IPW协变量数:", length(ps_covariates), "\n")
cat("   ", paste(ps_covariates, collapse = ", "), "\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 拟合PS模型
## ═══════════════════════════════════════════════════════

cat("【步骤 3】拟合倾向性评分模型\n")

ps_formula <- as.formula(paste("has_fu_echo ~", paste(ps_covariates, collapse = " + ")))

cat("  公式:\n  ", paste(deparse(ps_formula), collapse = ""), "\n\n")

fit_ps_list <- with(imp, glm(has_fu_echo ~ age + gender + resident + DM + hypertension + pPCI + 
                               STEMI + EF_baseline + LVEDV_baseline + GRACE_in + IN_killip + 
                               cTnIpeak + NTproBNP_peak + CKMB + WBC + HGB + PLT + CRP , family = binomial))

cat("  ✓ 拟合完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 合并系数
## ═══════════════════════════════════════════════════════

cat("【步骤 4】合并系数\n")

pooled_ps <- pool(fit_ps_list)

ps_coef_summary <- summary(pooled_ps, conf.int = TRUE) %>%
  as.data.frame() %>%
  mutate(
    OR = exp(estimate),
    CI_95 = paste0(round(exp(`2.5 %`), 3), " - ", round(exp(`97.5 %`), 3)),
    p = format.pval(p.value, digits = 3, eps = 0.001)
  ) %>%
  select(term, estimate, OR, CI_95, p)

print(ps_coef_summary, row.names = FALSE)

write.csv(ps_coef_summary,
          here::here("outputs", "ps_model_pooled_coefficients.csv"),
          row.names = FALSE)

cat("\n  ✓ 系数已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 5: 计算IPW权重 + 对数转换金属
## ═══════════════════════════════════════════════════════

cat("【步骤 5】计算IPW权重并对数转换金属\n")

vars_to_patch <- c("CKMB", "PLT", "HGB", "NTproBNP_peak")

# 定义金属变量
metal_vars <- c("Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
                "Mo", "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li")

imp_list_with_weights <- list()

for (i in 1:imp$m) {
  
  dat_i <- complete(imp, i)
  
  # 填补残留NA
  for (v in vars_to_patch) {
    if (v %in% names(dat_i) && sum(is.na(dat_i[[v]])) > 0) {
      dat_i[[v]][is.na(dat_i[[v]])] <- median(dat_i[[v]], na.rm = TRUE)
    }
  }
  
  # 拟合PS
  ps_formula <- as.formula(paste("has_fu_echo ~", paste(ps_covariates, collapse = " + ")))
  ps_model <- glm(ps_formula, data = dat_i, family = binomial())
  dat_i$ps <- predict(ps_model, newdata = dat_i, type = "response")
  dat_i$ps <- pmax(pmin(dat_i$ps, 0.99), 0.01)
  
  # 双向权重
  dat_i <- dat_i %>% 
    mutate(
      fu_binary = as.numeric(has_fu_echo == "Yes"),
      
      ipw_raw = case_when(
        has_fu_echo == "Yes" ~ 1 / ps,
        has_fu_echo == "No" ~ 1 / (1 - ps),
        TRUE ~ NA_real_
      ),
      
      p_fu = mean(fu_binary, na.rm = TRUE),
      sw = case_when(
        has_fu_echo == "Yes" ~ p_fu / ps,
        has_fu_echo == "No" ~ (1 - p_fu) / (1 - ps),
        TRUE ~ NA_real_
      )
    )
  
  # 截断
  q01 <- quantile(dat_i$sw, 0.01, na.rm = TRUE)
  q99 <- quantile(dat_i$sw, 0.99, na.rm = TRUE)
  
  dat_i$sw_trunc <- dat_i$sw
  dat_i$sw_trunc[dat_i$sw_trunc < q01] <- q01
  dat_i$sw_trunc[dat_i$sw_trunc > q99] <- q99
  
  # 【关键新增】对数转换金属
  cat("    插补数据集", i, ": 对数转换金属变量.. .\n")
  
  for (metal in metal_vars) {
    if (metal %in% names(dat_i)) {
      log_metal <- paste0("log_", metal)
      
      # 检查是否有零值或负值
      min_val <- min(dat_i[[metal]], na.rm = TRUE)
      
      if (min_val <= 0) {
        # 添加小常数 (最小正值的1%)
        shift <- abs(min_val) + 0.01
        dat_i[[log_metal]] <- log(dat_i[[metal]] + shift)
        cat("      ", metal, ": log(x + ", round(shift, 3), ")\n", sep="")
      } else {
        dat_i[[log_metal]] <- log(dat_i[[metal]])
        cat("      ", metal, ": log(x)\n", sep="")
      }
    }
  }
  
  # 计算结局
  dat_i <- dat_i %>% 
    mutate(delta_LVEDV = (LVEDV_fu - LVEDV_baseline) / LVEDV_baseline * 100)
  
  dat_i$.imp <- i
  dat_i$.id <- 1:nrow(dat_i)
  
  imp_list_with_weights[[i]] <- dat_i
}

imputed_data_with_ps <- bind_rows(imp_list_with_weights)

cat("  ✓ 权重计算完成\n")
cat("  ✓ 金属对数转换完成\n")
cat("  ✓ 新增变量:  ps, ipw_raw, sw, sw_trunc\n")
cat("  ✓ 新增变量: log_Cu, log_Zn, log_Fe, ...  (共", length(metal_vars), "个)\n\n", sep="")

## ═══════════════════════════════════════════════════════
## 步骤 6: 权重诊断
## ═══════════════════════════════════════════════════════

cat("【步骤 6】权重诊断\n")

weight_summary <- imputed_data_with_ps %>%
  group_by(.imp) %>%
  summarise(
    n = n(),
    ps_mean = mean(ps, na.rm = TRUE),
    sw_trunc_mean = mean(sw_trunc, na.rm = TRUE),
    sw_trunc_min = min(sw_trunc, na.rm = TRUE),
    sw_trunc_max = max(sw_trunc, na.rm = TRUE),
    .groups = "drop"
  )

print(head(weight_summary, 5), row.names = FALSE)

# 验证双向权重
dat_check <- imputed_data_with_ps %>% filter(.imp == 1)

sw_yes <- dat_check$sw_trunc[dat_check$has_fu_echo == "Yes"]
sw_no <- dat_check$sw_trunc[dat_check$has_fu_echo == "No"]

cat("\n  权重验证:\n")
cat("    有随访组: 均值=", round(mean(sw_yes, na.rm=T), 3), 
    ", 范围=[", round(min(sw_yes, na.rm=T), 3), ", ", 
    round(max(sw_yes, na.rm=T), 3), "]\n", sep="")
cat("    无随访组: 均值=", round(mean(sw_no, na.rm=T), 3), 
    ", 范围=[", round(min(sw_no, na.rm=T), 3), ", ", 
    round(max(sw_no, na.rm=T), 3), "]\n\n", sep="")

if (mean(sw_no, na.rm=T) == 0) {
  stop("❌ 无随访组权重为0!  ")
} else {
  cat("  ✓ 双向权重验证通过\n\n")
}

write.csv(weight_summary,
          here::here("outputs", "ipw_weights_summary.csv"),
          row.names = FALSE)

## ═══════════════════════════════════════════════════════
## 步骤 7: 可视化
## ═══════════════════════════════════════════════════════

cat("【步骤 7】生成诊断图\n")

ps_plot_data <- imputed_data_with_ps %>% filter(.imp <= 5)

p_ps <- ggplot(ps_plot_data, aes(x = ps, fill = has_fu_echo)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~.imp, ncol = 2) +
  scale_fill_manual(values = c("No" = "#E74C3C", "Yes" = "#3498DB")) +
  labs(title = "倾向性评分分布", x = "PS", y = "密度") +
  theme_minimal()

ggsave(here::here("plots", "ps_distribution.png"), p_ps, 
       width = 10, height = 8, dpi = 300)

cat("  ✓ PS分布图已保存\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 8: 保存数据
## ═══════════════════════════════════════════════════════

cat("【步骤 8】保存数据\n")

# 完整数据
saveRDS(imputed_data_with_ps,
        here::here("outputs", "imputed_data_with_ipw_weights.rds"))

# 第1个插补集
data_complete_imp1 <- imputed_data_with_ps %>% filter(.imp == 1)

saveRDS(data_complete_imp1,
        here::here("outputs", "complete_data_with_ipw_imp1.rds"))

write.csv(data_complete_imp1,
          here::here("outputs", "complete_data_with_ipw_imp1.csv"),
          row.names = FALSE)

# 结局分析数据
outcome_data <- imputed_data_with_ps %>%
  filter(has_fu_echo == "Yes", ! is.na(delta_LVEDV))

write.csv(outcome_data,
          here::here("outputs", "outcome_analysis_data_with_weights.csv"),
          row.names = FALSE)

saveRDS(outcome_data,
        here::here("outputs", "outcome_analysis_data_with_weights.rds"))

cat("  ✓ 所有数据已保存\n")
cat("    - 完整数据 (20个插补集): ", nrow(imputed_data_with_ps), " 行\n", sep="")
cat("    - 第1个插补集:  ", nrow(data_complete_imp1), " 行\n", sep="")
cat("    - 结局分析数据:  ", nrow(outcome_data), " 行\n\n", sep="")

## ═══════════════════════════════════════════════════════
## 步骤 9: 协变量平衡性检查 + SMD图 (整合版)
## ═══════════════════════════════════════════════════════

cat("【步骤 9】协变量平衡性检查与SMD可视化\n")

dat_imp1 <- data_complete_imp1

# 移除权重NA
if (any(is.na(dat_imp1$sw_trunc))) {
  dat_imp1 <- dat_imp1 %>% filter(!is.na(sw_trunc))
}

design_w <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_imp1)

balance_results <- data.frame(
  variable = character(),
  smd_unweighted = numeric(),
  smd_weighted = numeric(),
  stringsAsFactors = FALSE
)

## 9.1 计算SMD
cat("  计算SMD.. .\n")

for (var in ps_covariates) {
  
  vals <- dat_imp1[[var]]
  if (is.factor(vals) || is.character(vals)) {
    vals <- as.numeric(as.factor(vals)) - 1
  }
  
  vals_yes <- vals[dat_imp1$has_fu_echo == "Yes"]
  vals_no <- vals[dat_imp1$has_fu_echo == "No"]
  
  mean_yes <- mean(vals_yes, na.rm = TRUE)
  mean_no <- mean(vals_no, na.rm = TRUE)
  sd_pooled <- sqrt((var(vals_yes, na.rm=T) + var(vals_no, na.rm=T)) / 2)
  
  if (is.na(sd_pooled) || sd_pooled == 0) sd_pooled <- 1e-6
  
  smd_unw <- (mean_yes - mean_no) / sd_pooled
  
  fmla <- as.formula(paste0("~", var))
  
  w_res <- tryCatch({
    m_yes_w <- svymean(fmla, subset(design_w, has_fu_echo == "Yes"), na.rm=T)[1]
    m_no_w <- svymean(fmla, subset(design_w, has_fu_echo == "No"), na.rm=T)[1]
    list(my = m_yes_w, mn = m_no_w)
  }, error = function(e) NULL)
  
  if (!is.null(w_res)) {
    smd_w <- (w_res$my - w_res$mn) / sd_pooled
  } else {
    smd_w <- NA
  }
  
  balance_results <- rbind(balance_results, data.frame(
    variable = var,
    smd_unweighted = round(abs(smd_unw), 3),
    smd_weighted = round(abs(smd_w), 3)
  ))
}

balance_results <- balance_results %>%
  mutate(balanced = ! is.na(smd_weighted) & smd_weighted < 0.1)

cat("\n  SMD结果:\n")
print(balance_results, row.names = FALSE)

write.csv(balance_results,
          here::here("outputs", "covariate_balance_smd_imp1.csv"),
          row.names = FALSE)

n_balanced <- sum(balance_results$balanced, na.rm = TRUE)
n_total <- nrow(balance_results)

cat("\n  平衡性:  ", n_balanced, "/", n_total, " (", 
    round(n_balanced/n_total*100, 1), "%)\n", sep="")

## 9.2 生成SMD Balance Plot
cat("\n  生成SMD Balance Plot...\n")

# 加载标签
if (file.exists(here::here("variable_config.RData"))) {
  load(here::here("variable_config.RData"))
  
  balance_results <- balance_results %>% 
    left_join(label_mapping %>% select(variable, label = new_label), 
              by = "variable") %>% 
    mutate(label = ifelse(is.na(label), variable, label))
  
  cat("    ✓ 已应用中文标签\n")
} else {
  balance_results$label <- balance_results$variable
  cat("    ⚠ 使用原始变量名\n")
}

# 长格式
smd_long <- balance_results %>% 
  select(variable, label, smd_unweighted, smd_weighted) %>% 
  pivot_longer(cols = c(smd_unweighted, smd_weighted),
               names_to = "type",
               values_to = "smd") %>% 
  mutate(
    type_label = ifelse(type == "smd_unweighted", "加权前", "加权后"),
    type_label = factor(type_label, levels = c("加权前", "加权后"))
  ) %>% 
  filter(!is.na(smd))

# Love Plot
p_smd <- ggplot(smd_long, aes(x = smd, y = reorder(label, smd))) +
  geom_vline(xintercept = 0.1, linetype = "dashed", 
             color = "red", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", 
             color = "gray50", linewidth = 0.3) +
  geom_line(aes(group = label), color = "gray70", linewidth = 0.5) +
  geom_point(aes(color = type_label, shape = type_label), 
             size = 3.5, alpha = 0.85) +
  scale_color_manual(
    values = c("加权前" = "#D73027", "加权后" = "#4575B4"),
    name = ""
  ) +
  scale_shape_manual(
    values = c("加权前" = 16, "加权后" = 17),
    name = ""
  ) +
  labs(
    x = "标准化均数差 (Standardized Mean Difference)",
    y = "",
    title = "IPW加权前后协变量平衡性",
    subtitle = "虚线标记SMD=0.1阈值(良好平衡)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  ) +
  xlim(0, max(smd_long$smd, na.rm = TRUE) * 1.1)

ggsave(here::here("plots", "Figure_SMD_Balance_Plot.png"), 
       plot = p_smd, 
       width = 10, 
       height = max(8, nrow(balance_results) * 0.35),
       dpi = 300,
       bg = "white")

cat("    ✓ SMD图已保存\n")

# 统计摘要
cat("\n  平衡性改善:\n")
cat("    加权前平均SMD: ", round(mean(balance_results$smd_unweighted, na.rm=T), 3), "\n", sep="")
cat("    加权后平均SMD: ", round(mean(balance_results$smd_weighted, na.rm=T), 3), "\n", sep="")
reduction_pct <- (mean(balance_results$smd_unweighted, na.rm=T) - 
                    mean(balance_results$smd_weighted, na.rm=T)) / 
  mean(balance_results$smd_unweighted, na.rm=T) * 100
cat("    SMD降低:  ", round(reduction_pct, 1), "%\n\n", sep="")

## ═══════════════════════════════════════════════════════
## 步骤 10: 结局分析 (使用对数金属)
## ═══════════════════════════════════════════════════════

cat("【步骤 10】结局模型拟合 (对数金属)\n")

library(mitools)

imp_list <- split(outcome_data, outcome_data$.imp)
imputation_list_obj <- imputationList(imp_list)

svy_design <- svydesign(ids = ~1, weights = ~sw_trunc, data = imputation_list_obj)

# 【修改】使用对数转换后的金属
outcome_formula <- as.formula("delta_LVEDV ~ log_Cu + log_Zn + log_Fe + log_Se + log_Pb + age + gender + smoking")

cat("  拟合加权模型 (使用对数金属)...\n")
fit_results <- with(svy_design, svyglm(outcome_formula))

pooled_estimates <- MIcombine(fit_results)

summary_output <- summary(pooled_estimates) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  mutate(
    Beta = results,
    SE = se,
    Lower_95 = results - 1.96 * se,
    Upper_95 = results + 1.96 * se,
    P_Value = 2 * (1 - pnorm(abs(results / se)))
  )

cat("\n  合并结果:\n")
print(summary_output)

write.csv(summary_output, 
          here::here("outputs", "Final_Pooled_IPW_Results.csv"),
          row.names = FALSE)

cat("\n  ✓ 结局分析完成\n")
cat("  ✓ 解释:  Beta表示对数金属浓度增加1单位时,delta_LVEDV的变化\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║         IPW分析完成！                                ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("输出文件:\n")
cat("  数据:\n")
cat("    1. imputed_data_with_ipw_weights.rds (完整)\n")
cat("    2. complete_data_with_ipw_imp1.rds (第1个插补集)\n")
cat("    3. outcome_analysis_data_with_weights.rds (结局分析)\n")
cat("  结果:\n")
cat("    4. covariate_balance_smd_imp1.csv (SMD数据)\n")
cat("    5. Final_Pooled_IPW_Results.csv (结局分析)\n")
cat("  图形:\n")
cat("    6. ps_distribution.png (PS分布)\n")
cat("    7. Figure_SMD_Balance_Plot.png (SMD图)\n\n")

cat("关键指标:\n")
cat("  ✓ 双向权重:  两组均非零\n")
cat("  ✓ 协变量平衡: ", n_balanced, "/", n_total, " (", 
    round(n_balanced/n_total*100, 1), "%)\n", sep="")
cat("  ✓ SMD平均降低: ", round(reduction_pct, 1), "%\n\n", sep="")

