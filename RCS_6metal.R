############################################################
## RCS非线性分析 - 探索剂量-反应关系
## 针对：Se, Pb, Zn, Cu, Fe, As
## 选择金属（基于生物学合理性，而非P值）：
## 优先分析的金属（6个）：

## Se（硒）：抗氧化，文献最多，可能U型
## Pb（铅）：基线失衡最大（SMD=0.347），已知心血管毒性
## Zn（锌）：Model 2中Beta最大（-2.04），P=0.23接近0.2
## Cu（铜）：与心衰相关（Frustaci研究）
## Fe（铁）：氧化应激，心肌铁过载
## As（砷）：环境毒物，心血管风险明确
############################################################

rm(list = ls())
gc()
library(dplyr)
library(survey)
library(mitools)
library(rms)
library(ggplot2)
library(gridExtra)
library(here)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║     RCS非线性分析                                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 1: 加载数据
## ═══════════════════════════════════════════════════════

outcome_data <- readRDS(here:: here("outputs", "4 IPW2_complete", "outcome_analysis_data_with_weights_extended.rds"))

cat("【步骤 1】数据加载完成\n")
cat("  样本数:", nrow(outcome_data) / 20, "（每个插补）\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 定义节点（在合并数据上，确保跨插补一致）
## ═══════════════════════════════════════════════════════

cat("【步骤 2】定义RCS节点\n")

metals_rcs <- c("Se", "Pb", "Zn", "Cu", "Fe", "As")

# 在所有插补的合并数据上计算节点
knots_list <- list()

for(metal in metals_rcs) {
  log_metal <- paste0("log_", metal)
  
  if(log_metal %in% names(outcome_data)) {
    vals <- na.omit(outcome_data[[log_metal]])
    
    # 3节点（10th, 50th, 90th百分位数）
    knots <- quantile(vals, probs = c(0.1, 0.5, 0.9))
    knots_list[[metal]] <- knots
    
    cat("  ", metal, ": 节点 =", round(knots, 2), "\n")
  }
}

cat("\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 对每个金属拟合RCS模型
## ═══════════════════════════════════════════════════════

cat("【步骤 3】拟合RCS模型\n")

# 协变量
covariates <- c("age", "gender", "EF_baseline", "LVEDV_baseline", "cTnIpeak",
                "pPCI", "STEMI", "smoking", "DM", "hypertension",
                "NTproBNP_peak", "GRACE_in", "WBC", "CRP", "CHOL", "LDL", "AST")

# 存储结果
rcs_results <- list()
lrt_pvalues <- data.frame()

for(metal in metals_rcs) {
  log_metal <- paste0("log_", metal)
  
  cat("\n  分析", metal, ".. .\n")
  
  if(!  log_metal %in% names(outcome_data)) {
    cat("    ⊙ 跳过（变量不存在）\n")
    next
  }
  
  knots <- knots_list[[metal]]
  
  ## ─────────────────────────────
  ## 3.1 拟合线性模型
  ## ─────────────────────────────
  
  formula_linear <- as.formula(paste0(
    "delta_LVEDV ~ ", log_metal, " + ",
    paste(covariates, collapse = " + ")
  ))
  
  fit_linear_list <- lapply(1:20, function(i) {
    dat_i <- outcome_data %>% filter(.imp == i)
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    tryCatch({
      svyglm(formula_linear, design = design_i)
    }, error = function(e) NULL)
  })
  
  fit_linear_list <- Filter(Negate(is.null), fit_linear_list)
  
  ## ─────────────────────────────
  ## 3.2 拟合RCS模型
  ## ─────────────────────────────
  
  fit_rcs_list <- list()
  
  for(i in 1:20) {
    dat_i <- outcome_data %>% filter(.imp == i)
    
    # 创建RCS项（手动）
    # rcs()函数需要datadist，这里用rcspline. eval()替代
    rcs_terms <- rcspline.eval(dat_i[[log_metal]], knots = knots, 
                               inclx = TRUE)
    
    # 添加RCS项到数据
    dat_i$rcs1 <- rcs_terms[, 1]
    dat_i$rcs2 <- rcs_terms[, 2]
    
    # 拟合
    formula_rcs <- as.formula(paste0(
      "delta_LVEDV ~ rcs1 + rcs2 + ",
      paste(covariates, collapse = " + ")
    ))
    
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    fit_rcs_i <- tryCatch({
      svyglm(formula_rcs, design = design_i)
    }, error = function(e) NULL)
    
    if(! is.null(fit_rcs_i)) {
      fit_rcs_list[[i]] <- fit_rcs_i
    }
  }
  
  ## ─────────────────────────────
  ## 3.3 LRT检验非线性
  ## ─────────────────────────────
  
  if(length(fit_rcs_list) > 0 && length(fit_linear_list) > 0) {
    
    lrt_p_vec <- numeric(length(fit_rcs_list))
    
    for(i in 1:length(fit_rcs_list)) {
      # 计算似然比（近似）
      # 使用残差平方和比较
      rss_linear <- sum(residuals(fit_linear_list[[i]], type="response")^2)
      rss_rcs <- sum(residuals(fit_rcs_list[[i]], type="response")^2)
      
      # 卡方统计量（自由度 = RCS额外参数数 = 1）
      n <- nobs(fit_rcs_list[[i]])
      chi_sq <- n * log(rss_linear / rss_rcs)
      
      lrt_p_vec[i] <- 1 - pchisq(chi_sq, df = 1)
    }
    
    # 合并p值（使用Fisher方法）
    if(requireNamespace("metap", quietly = TRUE)) {
      library(metap)
      pooled_lrt_p <- sumz(lrt_p_vec)$p
    } else {
      # 简单平均（保守）
      pooled_lrt_p <- mean(lrt_p_vec)
    }
    
    lrt_pvalues <- rbind(lrt_pvalues, data.frame(
      Metal = metal,
      LRT_P_pooled = pooled_lrt_p,
      Nonlinear = ifelse(pooled_lrt_p < 0.05, "是", "否")
    ))
    
    cat("    LRT p值（合并）:", round(pooled_lrt_p, 4))
    
    if(pooled_lrt_p < 0.05) {
      cat(" ✓ 存在非线性\n")
    } else {
      cat(" (线性足够)\n")
    }
    
    # 保存拟合对象（用于后续预测）
    rcs_results[[metal]] <- list(
      fit_list = fit_rcs_list,
      knots = knots,
      lrt_p = pooled_lrt_p
    )
  }
}

cat("\n  ✓ RCS拟合完成\n\n")

# 保存LRT结果
write.csv(lrt_pvalues,
          here::here("outputs", "RCS_Nonlinearity_Test_Results.csv"),
          row.names = FALSE)

## ═══════════════════════════════════════════════════════
## 步骤 4: 生成剂量-反应曲线
## ═══════════════════════════════════════════════════════

cat("【步骤 4】生成剂量-反应曲线\n")

plot_list <- list()

for(metal in metals_rcs) {
  
  if(! metal %in% names(rcs_results)) next
  
  cat("  绘制", metal, ".. .\n")
  
  log_metal <- paste0("log_", metal)
  knots <- rcs_results[[metal]]$knots
  fit_list <- rcs_results[[metal]]$fit_list
  
  # 生成预测网格
  log_metal_range <- range(outcome_data[[log_metal]], na.rm = TRUE)
  pred_grid_log <- seq(log_metal_range[1], log_metal_range[2], length.out = 100)
  
  # 对每个插补预测
  pred_list <- list()
  
  for(i in 1:length(fit_list)) {
    dat_i <- outcome_data %>% filter(.imp == i)
    
    # 创建预测数据（协变量取中位数/众数）
    newdata_i <- data.frame(
      rcs1 = rcspline.eval(pred_grid_log, knots = knots, inclx = TRUE)[, 1],
      rcs2 = rcspline.eval(pred_grid_log, knots = knots, inclx = TRUE)[, 2],
      age = median(dat_i$age, na.rm = TRUE),
      gender = 1,
      EF_baseline = median(dat_i$EF_baseline, na.rm = TRUE),
      LVEDV_baseline = median(dat_i$LVEDV_baseline, na.rm = TRUE),
      cTnIpeak = median(dat_i$cTnIpeak, na.rm = TRUE),
      pPCI = 1,
      STEMI = 1,
      smoking = 0,
      DM = 0,
      hypertension = 1,
      NTproBNP_peak = median(dat_i$NTproBNP_peak, na.rm = TRUE),
      GRACE_in = median(dat_i$GRACE_in, na.rm = TRUE),
      WBC = median(dat_i$WBC, na.rm = TRUE),
      CRP = median(dat_i$CRP, na.rm = TRUE),
      CHOL = median(dat_i$CHOL, na.rm = TRUE),
      LDL = median(dat_i$LDL, na.rm = TRUE),
      AST = median(dat_i$AST, na.rm = TRUE)
    )
    
    # 预测
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    pred_i <- predict(fit_list[[i]], newdata = newdata_i, se.fit = TRUE)
    
    pred_list[[i]] <- data.frame(
      log_metal_val = pred_grid_log,
      fit = pred_i$fit,
      se = pred_i$se.fit,
      .imp = i
    )
  }
  
  # 合并预测（Rubin规则）
  pred_all <- bind_rows(pred_list)
  
  pred_pooled <- pred_all %>%
    group_by(log_metal_val) %>%
    summarise(
      fit_mean = mean(fit),
      se_within = sqrt(mean(se^2)),
      se_between = sd(fit),
      se_total = sqrt(se_within^2 + (1 + 1/20) * se_between^2),
      lower = fit_mean - 1.96 * se_total,
      upper = fit_mean + 1.96 * se_total,
      .groups = "drop"
    )
  
  # 转换回原始尺度（如果需要）
  # 这里保持log尺度，横轴标注原始浓度
  
  # 计算原始金属浓度（用于横轴标签）
  original_metal_range <- exp(log_metal_range)
  
  # 绘图
  p <- ggplot(pred_pooled, aes(x = exp(log_metal_val), y = fit_mean)) +
    geom_line(size = 1.2, color = "#2E86AB") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#2E86AB") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
    
    # 添加参考线（中位数）
    geom_vline(xintercept = exp(knots[2]), linetype = "dotted", 
               color = "gray40", size = 0.5) +
    
    # 标注节点
    annotate("text", x = exp(knots[2]), y = max(pred_pooled$upper) * 0.9,
             label = "中位数", size = 3, color = "gray40", hjust = -0.1) +
    
    labs(
      x = paste0(metal, " 血清浓度 (μg/L)"),
      y = "ΔLVEDV预测值 (%)",
      title = paste0(metal, " 与左室重构的剂量-反应关系"),
      subtitle = paste0("RCS (3节点) | LRT p = ", 
                        round(rcs_results[[metal]]$lrt_p, 3))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40")
    )
  
  plot_list[[metal]] <- p
  
  # 保存单独图片
  ggsave(here:: here("plots", paste0("RCS_", metal, "_DoseResponse.png")),
         p, width = 8, height = 6, dpi = 300, bg = "white")
}

# 组合图（6张拼接）
p_combined <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)

ggsave(here::here("plots", "RCS_All_Metals_Combined.png"),
       p_combined, width = 14, height = 18, dpi = 300, bg = "white")

cat("  ✓ 剂量-反应曲线已保存\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║         RCS非线性分析完成！                          ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("输出文件:\n")
cat("  1. RCS_Nonlinearity_Test_Results.csv（LRT检验结果）\n")
cat("  2. RCS_[Metal]_DoseResponse.png（单个金属）\n")
cat("  3. RCS_All_Metals_Combined.png（组合图）\n\n")

cat("解释提示:\n")
cat("  - 如果LRT p < 0.05：存在非线性，解释曲线形状\n")
cat("  - 如果LRT p > 0.05：线性关系足够，报告线性Beta\n")
cat("  - 关注：U型、倒U型、阈值效应\n\n")

############################################################
## END
############################################################