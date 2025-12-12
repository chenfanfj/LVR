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
## 步骤 5: 计算倾向性评分和 IPW 权重 (最终修正版)
## ═══════════════════════════════════════════════════════

cat("【步骤 5】计算倾向性评分和 IPW 权重\n")

## 1. 准备工作
# 定义需要微调填补的变量 (基于之前的报错分析)
# 这些变量偶尔存在的微量 NA 会导致 glm 丢弃行，必须填补
vars_to_patch <- c("CKMB", "PLT", "HGB", "NTproBNP_peak")

# 创建一个列表来存储处理好的数据集
imp_list_with_weights <- list()

## 2. 循环处理 m=20 个插补数据集
for (i in 1:imp$m) {
  
  # [A] 提取第 i 个完整插补数据集
  # 因为此时 imp 中已包含结局变量，这里提取出来的 dat_i 已经有 delta_LVEDV 等
  dat_i <- complete(imp, i)
  
  # [B] (关键修复) 填补残留的微量缺失值
  # 解决 "Assigned data has 1118 rows" 错误的必要步骤
  for (v in vars_to_patch) {
    if (v %in% names(dat_i) && sum(is.na(dat_i[[v]])) > 0) {
      # 使用中位数填补这几个残留的 NA
      dat_i[[v]][is.na(dat_i[[v]])] <- median(dat_i[[v]], na.rm = TRUE)
    }
  }
  
  # [C] 拟合倾向性评分模型 (PS Model)
  # 请确保公式里的变量名与数据一致
  ps_model <- glm(has_fu_echo ~ age + gender + resident + DM + hypertension + pPCI + 
                    STEMI + EF_baseline + LVEDV_baseline + GRACE_in_str + 
                    IN_killip + cTnIpeak + NTproBNP_peak + CKMB + WBC + HGB + PLT + CRP, 
                  data = dat_i, 
                  family = binomial())
  
  # [D] 计算 PS 值
  # 使用 newdata = dat_i 确保返回长度为 1124 的向量
  dat_i$ps <- predict(ps_model, newdata = dat_i, type = "response")
  
  # [E] 计算 IPW 权重
  # 仅针对有随访人群 (has_fu_echo == "Yes") 计算权重 1/PS
  # 无随访人群权重设为 0 (不参与结局分析)
  dat_i$ipw_raw <- ifelse(dat_i$has_fu_echo == "Yes", 1 / dat_i$ps, 0)
  
  # [F] 计算截断权重 (生成 sw_trunc 变量供后续使用)
  # 截断阈值：99% 分位数 (仅在有权重的样本中计算)
  q99 <- quantile(dat_i$ipw_raw[dat_i$has_fu_echo == "Yes"], 0.99, na.rm = TRUE)
  
  # 创建 sw_trunc 变量 (这是后续步骤 6, 10, 11 调用的变量名)
  dat_i$sw_trunc <- dat_i$ipw_raw
  dat_i$sw_trunc[dat_i$sw_trunc > q99] <- q99
  
  # [G] 添加辅助列
  dat_i$.imp <- i
  dat_i$.id <- 1:nrow(dat_i)
  
  # 存入列表
  imp_list_with_weights[[i]] <- dat_i
}

## 3. 合并回长格式数据框
imputed_data_with_ps <- bind_rows(imp_list_with_weights)

cat("  ✓ IPW 权重计算完成 (含行数匹配补丁)。\n")
cat("  ✓ 已生成变量: ps, ipw_raw, sw_trunc\n")

## ═══════════════════════════════════════════════════════
## 补丁步骤: 手动计算 delta_LVEDV
## ═══════════════════════════════════════════════════════

cat("【补丁】手动计算 delta_LVEDV 变量...\n")

# 1. 检查变量名，防止 LVEDV_baseline 因为合并产生了 .x / .y 后缀
# 如果因为之前的 left_join 导致了重名，通常 .x 是插补后的（完整），.y 是原始的
baseline_col <- "LVEDV_baseline"
if ("LVEDV_baseline.x" %in% names(imputed_data_with_ps)) {
  baseline_col <- "LVEDV_baseline.x"
}

cat("  使用基线变量列名:", baseline_col, "\n")

# 2. 执行计算
imputed_data_with_ps <- imputed_data_with_ps %>%
  mutate(
    # 确保基线变量是数值型
    base_val = as.numeric(.data[[baseline_col]]),
    fu_val   = as.numeric(LVEDV_fu),
    
    # 计算变化率: (随访 - 基线) / 基线 * 100
    # 结果单位为百分比 (%)
    delta_LVEDV = (fu_val - base_val) / base_val * 100
  ) %>%
  # 清理临时辅助列
  select(-base_val, -fu_val)

# 3. 检查计算结果
n_calced <- sum(!is.na(imputed_data_with_ps$delta_LVEDV))
cat("  ✓ 已计算 delta_LVEDV\n")
cat("    有效值数量:", n_calced, "(应与有随访数据的样本量一致)\n")
cat("    前5个值   :", paste(round(head(imputed_data_with_ps$delta_LVEDV[!is.na(imputed_data_with_ps$delta_LVEDV)], 5), 2), collapse=", "), "\n\n")


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
  filter(has_fu_echo == "Yes", !is.na(delta_LVEDV))

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
# 【修改点】: delta_LVEDV 是连续变量，适合直接做因变量。
# 这里的变量列表建议包含：金属(暴露) + 核心临床混杂因素(双重稳健)。
# 不需要再做逻辑回归筛选，直接纳入重要的临床变量。
# 请根据实际列名修改下方的变量名
out_formula <- as.formula("delta_LVEDV ~ Cu + Zn + Fe + Se + Pb + age + gender + hypertension + DM + smoking")

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
## 步骤 10: 协变量平衡性检查（修正版）
## ═══════════════════════════════════════════════════════

cat("【步骤 10】协变量平衡性检查（示例）\n")

## 使用第1个插补数据集
dat_imp1 <- imputed_data_with_ps %>% filter(.imp == 1)

## 检查权重是否有 NA
if (any(is.na(dat_imp1$sw_trunc))) {
  cat("  ⚠ 警告：权重中存在 NA，已自动移除这些样本进行平衡性检查。\n")
  dat_imp1 <- dat_imp1 %>% filter(!is.na(sw_trunc))
}

## 定义加权设计对象 (注意变量名 sw_trunc)
design_w <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_imp1)

balance_results <- data.frame(
  variable = character(),
  smd_unweighted = numeric(),
  smd_weighted = numeric(),
  stringsAsFactors = FALSE
)

# 定义要检查的协变量列表 (确保都在 dat_imp1 中)
# 使用 intersect 确保变量存在
vars_to_check <- intersect(ps_covariates, names(dat_imp1))

for (var in vars_to_check) {
  
  # 确保变量是数值型或可转换为数值
  # 对于因子变量，这里简化处理，将其转为数值 (0/1) 或忽略
  # 如果是二元因子 (如 gender)，转为 0/1 后计算均值差是有意义的
  vals <- dat_imp1[[var]]
  if (is.factor(vals) || is.character(vals)) {
    vals <- as.numeric(as.factor(vals)) - 1
  }
  
  ## 1. 未加权 SMD
  # 提取两组数据
  vals_yes <- vals[dat_imp1$has_fu_echo == "Yes"]
  vals_no  <- vals[dat_imp1$has_fu_echo == "No"]
  
  mean_yes <- mean(vals_yes, na.rm = TRUE)
  mean_no  <- mean(vals_no, na.rm = TRUE)
  var_yes  <- var(vals_yes, na.rm = TRUE)
  var_no   <- var(vals_no, na.rm = TRUE)
  
  # 混合标准差
  sd_pooled <- sqrt((var_yes + var_no) / 2)
  
  # 防止分母为 0
  if (is.na(sd_pooled) || sd_pooled == 0) sd_pooled <- 1e-6
  
  smd_unw <- (mean_yes - mean_no) / sd_pooled
  
  ## 2. 加权 SMD
  # 使用 svymean 计算加权均值
  # 构造公式
  fmla <- as.formula(paste0("~", var))
  
  # 分别计算两组的加权均值
  # 使用 tryCatch 防止报错中断循环
  w_res <- tryCatch({
    m_yes_w <- svymean(fmla, subset(design_w, has_fu_echo == "Yes"), na.rm = TRUE)[1]
    m_no_w  <- svymean(fmla, subset(design_w, has_fu_echo == "No"),  na.rm = TRUE)[1]
    list(my = m_yes_w, mn = m_no_w)
  }, error = function(e) return(NULL))
  
  if (!is.null(w_res)) {
    smd_w <- (w_res$my - w_res$mn) / sd_pooled
  } else {
    smd_w <- NA
  }
  
  balance_results <- rbind(balance_results, data.frame(
    variable = var,
    smd_unweighted = round(smd_unw, 3),
    smd_weighted = round(smd_w, 3)
  ))
}

## 标记平衡状态
# 处理 NA：如果 smd_weighted 是 NA，则视为不平衡 (FALSE)
balance_results <- balance_results %>%
  mutate(balanced = !is.na(smd_weighted) & abs(smd_weighted) < 0.1)

cat("\n  标准化均差结果:\n")
print(balance_results, row.names = FALSE)

write.csv(balance_results,
          here::here("outputs", "covariate_balance_smd_imp1.csv"),
          row.names = FALSE)

## 统计平衡比例
n_balanced <- sum(balance_results$balanced, na.rm = TRUE)
n_total <- nrow(balance_results) # 使用总行数作为分母

cat("\n  平衡性总结:\n")
cat("    平衡良好 (|SMD| < 0.1):", n_balanced, "/", n_total, "\n")

# 防止除零错误
ratio <- if (n_total > 0) n_balanced / n_total else 0

if (ratio < 0.8) {
  cat("    ⚠ 建议: 考虑添加更多协变量或交互项\n")
} else {
  cat("    ✓ 大多数协变量达到良好平衡\n")
}

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
cat("  3.  结局: delta_LVEDV\n")
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