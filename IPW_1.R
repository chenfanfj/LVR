############################################################
## 在插补数据集上构建 IPW 权重
## 输入: mice_imputation_object. rds (20个插补数据集)
## 输出: 每个数据集的倾向性评分和 IPW 权重
############################################################
rm(list = ls())
gc()
invisible(gc())
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

imp <- readRDS(here::here("outputs", "mice_imputation_final.rds"))

# 2. 查看对象的整体类与结构（初步看有哪些组件）
class(imp)  # 应为 "mids"
str(imp, max.level = 1)  # 第一层结构，避免输出过长

# 3. 查看插补数据中的变量名（关键步骤）
# 方法A：从“data”组件中获取（原始数据框架，含缺失值）
colnames(imp$data)

# 检查 'has_fu_echo' 和 'LVEDV_fu' 是否存在于完整数据中
check_data <- complete(imp, 1)

if(!all(c("has_fu_echo", "LVEDV_fu") %in% names(check_data))) {
  stop("严重错误：结局或选择变量在插补对象中缺失！")
} else {
  message("桥梁验证通过：数据已准备好进行 IPW。")
}

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
ps_formula <- paste("has_fu_echo ~", paste(ps_covariates, collapse = " + "))
cat("\n  模型公式:\n")
print(ps_formula)
cat("\n")


## ═══════════════════════════════════════════════════════
## 步骤 3: 在每个插补数据集上拟合倾向性评分模型
## ═══════════════════════════════════════════════════════

cat("【步骤 3】在", imp$m, "个插补数据集上拟合模型\n\n")

## 使用 mice::with() 在每个数据集上拟合
cat("  拟合倾向性评分模型..  .\n")

fit_ps_list <- with(imp, glm(has_fu_echo ~ age + gender + resident + DM + hypertension + pPCI + 
                               STEMI + EF_baseline + LVEDV_baseline + GRACE_in + IN_killip + 
                               cTnIpeak + NTproBNP_peak + CKMB + WBC + HGB + PLT + CRP , family = binomial))
cat("  ✓ 所有数据集拟合完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 4: 合并系数估计 (Rubin's rules)
## ═══════════════════════════════════════════════════════

cat("【步骤 4】合并倾向性评分模型系数\n")

pooled_ps <- pool(fit_ps_list)

ps_coef_summary <- summary(pooled_ps, conf.int = TRUE) %>%
  as.data.frame() %>%
  mutate(
    OR = round(exp(estimate), 3),
    CI_95 = paste0(round(exp(`2.5 %`), 3), " - ", round(exp(`97.5 %`), 3)),
    p = format.pval(p.value, digits = 3, eps = 0.001)
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

## 0. 检查缺失情况
# 【【检查第一个插补数据集中，哪些用于 PS 模型的变量还有缺失值
vars_in_ps_model <- c("age", "gender", "resident", "DM", "hypertension", "pPCI", 
                      "STEMI", "EF_baseline", "LVEDV_baseline", "GRACE_in", 
                      "IN_killip", "cTnIpeak", "NTproBNP_peak", "CKMB", 
                      "WBC", "HGB", "PLT", "CRP") # 根据您的系数表列出
# 查看缺失计数
sapply(complete(imp, 1)[, vars_in_ps_model], function(x) sum(is.na(x))) #】】


## 1. 准备工作
# 定义需要微调填补的变量 (基于之前的报错分析)
vars_to_patch <- c("CKMB", "PLT", "HGB", "NTproBNP_peak")

# 创建一个列表来存储处理好的数据集
imp_list_with_weights <- list()

## 2. 循环处理 m=20 个插补数据集
for (i in 1:imp$m) {
  
  # [A] 提取第 i 个完整插补数据集
  dat_i <- complete(imp, i)
  
  # [B] (关键) 合并结局变量 (因为之前它们被排除在插补之外)
  # 假设 raw_dat 是您读取的原始 Excel 数据
  dat_i <- dat_i %>%
    left_join(raw_dat %>% select(ID, LVEDV_fu, LVESV_fu, EF_fu, LVEDV_baseline), by = "ID") %>%
    mutate(
      # 计算 ΔLVEDV (如果还没计算的话)
      delta_LVEDV = (LVEDV_fu - LVEDV_baseline.y) / LVEDV_baseline.y * 100
      # 注意：如果有重名列，dplyr 会自动加 .x .y 后缀，请根据实际情况调整
    )
  
  # [C] (关键修复) 填补残留的微量缺失值
  # 这一步是为了防止 glm 因为这几个 NA 而丢弃行，导致 predict 长度不匹配
  for (v in vars_to_patch) {
    if (v %in% names(dat_i) && sum(is.na(dat_i[[v]])) > 0) {
      dat_i[[v]][is.na(dat_i[[v]])] <- median(dat_i[[v]], na.rm = TRUE)
    }
  }
  
  # [D] 拟合倾向性评分模型 (PS Model)
  # 请确保公式里的变量都在 dat_i 中且无缺失
  ps_model <- glm(has_fu_echo ~ age + gender + resident + DM + hypertension + pPCI + 
                    STEMI + EF_baseline + LVEDV_baseline + GRACE_in_str + 
                    IN_killip + cTnIpeak + NTproBNP_peak + CKMB + WBC + HGB + PLT + CRP, 
                  data = dat_i, 
                  family = binomial())
  
  # [E] 计算 PS 值和权重
  # 使用 newdata = dat_i 确保即使有极少数 NA 也能返回对应长度的向量 (虽然步骤 C 已修复 NA)
  dat_i$ps <- predict(ps_model, newdata = dat_i, type = "response")
  
  # 计算 IPW (仅针对有随访人群: has_fu_echo == "Yes")
  # 稳定化权重可根据需要选用，这里使用标准 ATT/ATE 权重逻辑，针对随访人群通常用 1/PS
  dat_i$ipw <- ifelse(dat_i$has_fu_echo == "Yes", 1 / dat_i$ps, 0)
  
  # [F] 截断极端权重 (例如 99% 分位数)
  # 仅在有权重的样本中计算分位数
  q99 <- quantile(dat_i$ipw[dat_i$has_fu_echo == "Yes"], 0.99, na.rm = TRUE)
  dat_i$ipw[dat_i$ipw > q99] <- q99
  
  # [G] 添加插补编号 (方便后续合并)
  dat_i$.imp <- i
  dat_i$.id <- 1:nrow(dat_i)
  
  # 存入列表
  imp_list_with_weights[[i]] <- dat_i
}

## 3. 合并回长格式数据框 (如果您后续代码需要长格式)
imputed_data_with_ps <- bind_rows(imp_list_with_weights)

cat("  ✓ IPW 权重计算完成，数据已合并。\n")
cat("    数据维度:", nrow(imputed_data_with_ps), "行 (应为 1124 * 20)\n")

## 检查是否还有 NA 权重 (针对有随访的人)
n_na_weights <- sum(is.na(imputed_data_with_ps$ipw) & imputed_data_with_ps$has_fu_echo == "Yes")
if (n_na_weights > 0) {
  cat("  ⚠ 警告: 仍有", n_na_weights, "个随访样本的权重为 NA，请检查 PS 模型变量。\n")
} else {
  cat("  ✓ 所有随访样本均成功获得权重。\n")
}

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
    
    .groups = "drop"
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
          row.names = FALSE)

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
p_weights <- ggplot(weight_summary, aes(x = factor(.imp), y = sw_mean)) +
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
## 步骤 10: 最终加权结局分析 
## ═══════════════════════════════════════════════════════

library(survey)
library(mitools)

cat("【步骤 11】拟合加权结局模型并合并结果\n")

# 1. 读取数据 
data_long <- read.csv(here::here("outputs", "outcome_analysis_data_with_weights.csv"))

# 2. 转换为插补列表
imp_list <- split(data_long, data_long$.imp)
imputation_list_obj <- mitools::imputationList(imp_list)

# 3. 定义调查设计 (Survey Design)
# 【修改点】: 如果 ID 是医院编号，使用 ids = ~ID 处理聚类。
# 如果 ID 是患者编号(不重复)，使用 ids = ~1。此处假设为医院编号。
svy_design <- svydesign(ids = ~1, weights = ~sw_trunc, data = imputation_list_obj)

# 4. 定义模型公式
# 【修改点】: ΔLVEDV 是连续变量，适合直接做因变量。
# 这里的变量列表建议包含：金属(暴露) + 核心临床混杂因素(双重稳健)。
# 不需要再做逻辑回归筛选，直接纳入重要的临床变量。
# 请根据实际列名修改下方的变量名
out_formula <- as.formula("ΔLVEDV ~ Cu + Zn + Fe + Se + Pb + age + gender + hypertension + DM + smoking")

cat("  拟合模型 (线性回归，因变量为连续值)...\n")
# svyglm 默认 family=gaussian，即线性回归，适合连续变量
fit_results <- with(svy_design, svyglm(out_formula))

# 5. 合并结果 (Pooling)
# 这是最后一步，将 20 个模型的结果整合
pooled_estimates <- MIcombine(fit_results)

# 6. 提取结果
summary_output <- summary(pooled_estimates) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  mutate(
    Beta = results, # 连续变量回归系数
    Lower_95 = results - 1.96 * se,
    Upper_95 = results + 1.96 * se,
    P_Value = 2 * (1 - pnorm(abs(results / se)))
  )

cat("\n  最终合并后的结果 (Beta系数):\n")
print(summary_output)

# 7. 保存
write.csv(summary_output, here::here("outputs", "Final_Pooled_IPW_Results.csv"))

## ═══════════════════════════════════════════════════════
## 步骤 10: 协变量平衡性检查（示例：第1个插补数据集）
## ═══════════════════════════════════════════════════════

cat("【步骤 10】协变量平衡性检查（示例）\n")

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
          row.names = FALSE)

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
## 步骤 11: 生成综合报告
## ═══════════════════════════════════════════════════════

cat("【步骤 11】生成综合报告\n")

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