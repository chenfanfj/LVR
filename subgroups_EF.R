############################################################
##  A2: 亚组分析 (效应修饰检验)
##  输入: imputed_data_with_ipw_weights.rds
##  输出: 亚组森林图 + 交互作用检验表
############################################################

rm(list = ls())
gc()

library(mice)
library(survey)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forestplot)
library(here)
library(grid)


cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║        A2: 金属与LVR的亚组分析 (效应修饰)          ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤1: 读取数据
## ═══════════════════════════════════════════════════════

cat("【步骤1】读取数据\n")

imputed_data <- readRDS(here("outputs", "4 IPW2_complete", "imputed_data_with_ipw_weights_extended.rds"))

outcome_data <- imputed_data %>% 
  filter(has_fu_echo == "Yes", ! is.na(LVEDV_fu)) %>%
  mutate(LVR = as.numeric(EF_baseline >= 50 & EF_fu < 50))

cat("  样本数:", nrow(outcome_data), "\n\n")

## ═══════════════════════════════════════════════════════
## 步骤2: 定义亚组与核心金属
## ═══════════════════════════════════════════════════════

# 核心金属 (仅分析显著/边缘显著的)
metals_to_analyze <- c("Pb", "Zn", "Fe")

# 定义亚组变量
subgroups <- list(
  
  # 1. 基线EF分层
  EF_baseline_group = list(
    var = "EF_baseline",
    cut = 55,
    label_low = "低EF (<50%)",
    label_high = "正常EF (≥50%)"
  ),
  
  # 2. 肾功能分层
  EGFR_group = list(
    var = "EGFR",
    cut = 60,
    label_low = "肾功能不全 (<60)",
    label_high = "肾功能正常 (≥60)"
  ),
  
  # 3. 糖尿病
  DM_group = list(
    var = "DM",
    label_0 = "无糖尿病",
    label_1 = "糖尿病"
  ),
  
  # 4. 年龄
  Age_group = list(
    var = "age",
    cut = 65,
    label_low = "< 65岁",
    label_high = "≥ 65岁"
  ),
  
  # 5. 急诊PCI
  pPCI_group = list(
    var = "pPCI",
    label_0 = "非急诊PCI",
    label_1 = "急诊PCI"
  )
)

# 协变量
covariates <- c("age", "gender", "EF_baseline", "LVEDV_baseline", 
                "cTnIpeak", "pPCI", "STEMI", "smoking", "DM", 
                "hypertension", "NTproBNP_peak", "GRACE_in", 
                "WBC", "HGB", "CRP", "CHOL", "LDL", "AST")

## ═══════════════════════════════════════════════════════
## 步骤3: 创建亚组变量
## ═══════════════════════════════════════════════════════

cat("【步骤3】创建亚组变量\n")

outcome_data <- outcome_data %>%
  mutate(
    # EF分层
    EF_baseline_group = ifelse(EF_baseline < 50, "低EF (<50%)", "正常EF (≥50%)"),
    
    # EGFR分层
    EGFR_group = ifelse(EGFR < 60, "肾功能不全 (<60)", "肾功能正常 (≥60)"),
    
    # 糖尿病
    DM_group = ifelse(DM == 1, "糖尿病", "无糖尿病"),
    
    # 年龄
    Age_group = ifelse(age < 65, "< 65岁", "≥ 65岁"),
    
    # pPCI
    pPCI_group = ifelse(pPCI == 1, "急诊PCI", "非急诊PCI")
  )

cat("  ✓ 亚组变量已创建\n\n")

## ═══════════════════════════════════════════════════════
## 步骤4: 亚组分析函数
## ═══════════════════════════════════════════════════════

run_subgroup_analysis <- function(data, metal, subgroup_var, covariates) {
  
  log_metal <- paste0("log_", metal)
  
  # 移除协变量中的亚组变量 (避免共线性)
  covars_adjusted <- setdiff(covariates, c(subgroup_var, gsub("_group", "", subgroup_var)))
  
  # 获取亚组水平
  subgroup_levels <- unique(data[[subgroup_var]])
  subgroup_levels <- subgroup_levels[! is.na(subgroup_levels)]
  
  results_list <- list()
  
  # 对每个亚组分别拟合
  for (level in subgroup_levels) {
    
    # 筛选子集
    data_sub <- data %>% filter(!! sym(subgroup_var) == level)
    
    # 移除缺失
    data_sub_complete <- data_sub %>%
      select(LVR, all_of(log_metal), all_of(covars_adjusted), sw_trunc) %>%
      na.omit()
    
    if (nrow(data_sub_complete) < 50) {
      # 样本量太小，跳过
      next
    }
    
    # 加权设计
    svy_design_sub <- svydesign(ids = ~1, weights = ~sw_trunc, data = data_sub_complete)
    
    # 拟合模型
    formula_sub <- as.formula(paste(
      "LVR ~", log_metal, "+",
      paste(covars_adjusted, collapse = " + ")
    ))
    
    fit_sub <- svyglm(formula_sub, design = svy_design_sub, family = quasibinomial)
    
    # 提取系数
    coef_metal <- coef(fit_sub)[log_metal]
    se_metal <- summary(fit_sub)$coefficients[log_metal, "Std. Error"]
    
    # OR和95%CI
    OR <- exp(coef_metal)
    OR_lower <- exp(coef_metal - 1.96 * se_metal)
    OR_upper <- exp(coef_metal + 1.96 * se_metal)
    P_value <- summary(fit_sub)$coefficients[log_metal, "Pr(>|t|)"]
    
    # 存储
    results_list[[level]] <- data.frame(
      Subgroup_level = level,
      N = nrow(data_sub_complete),
      OR = OR,
      OR_lower = OR_lower,
      OR_upper = OR_upper,
      P_value = P_value
    )
  }
  
  # 合并结果
  results_df <- bind_rows(results_list)
  
  # 交互作用检验 (完整数据)
  data_complete <- data %>%
    select(LVR, all_of(log_metal), all_of(subgroup_var), all_of(covars_adjusted), sw_trunc) %>%
    na.omit()
  
  svy_design_full <- svydesign(ids = ~1, weights = ~sw_trunc, data = data_complete)
  
  # 交互模型
  formula_interaction <- as.formula(paste(
    "LVR ~", log_metal, "*", subgroup_var, "+",
    paste(covars_adjusted, collapse = " + ")
  ))
  
  fit_interaction <- svyglm(formula_interaction, design = svy_design_full, family = quasibinomial)
  
  # 提取交互项P值
  interaction_term <- paste0(log_metal, ":", subgroup_var)
  
  # 检查交互项是否存在
  if (interaction_term %in% names(coef(fit_interaction))) {
    P_interaction <- summary(fit_interaction)$coefficients[interaction_term, "Pr(>|t|)"]
  } else {
    # 如果是多分类，取第一个交互项
    interaction_terms <- grep(paste0(log_metal, ":"), names(coef(fit_interaction)), value = TRUE)
    if (length(interaction_terms) > 0) {
      P_interaction <- summary(fit_interaction)$coefficients[interaction_terms[1], "Pr(>|t|)"]
    } else {
      P_interaction <- NA
    }
  }
  
  results_df$P_interaction <- P_interaction
  
  return(results_df)
}

## ═══════════════════════════════════════════════════════
## 步骤5: 对20个插补分别分析 + 合并
## ═══════════════════════════════════════════════════════

cat("【步骤5】亚组分析拟合\n")

subgroup_results_all <- list()

for (metal in metals_to_analyze) {
  
  cat("\n  分析金属:", metal, "\n")
  
  for (subgroup_name in names(subgroups)) {
    
    cat("    亚组:", subgroup_name, "\n")
    
    # 对每个插补
    results_by_imp <- list()
    
    for (i in 1:20) {
      
      dat_i <- outcome_data %>% filter(.imp == i)
      
      tryCatch({
        res <- run_subgroup_analysis(dat_i, metal, subgroup_name, covariates)
        res$.imp <- i
        results_by_imp[[i]] <- res
      }, error = function(e) {
        cat("      插补", i, "失败:", e$message, "\n")
      })
    }
    
    # 合并20个插补
    results_combined <- bind_rows(results_by_imp)
    
    # 按亚组水平汇总 (取均值)
    results_pooled <- results_combined %>%
      group_by(Subgroup_level) %>%
      summarise(
        N = round(mean(N, na.rm = TRUE)),
        OR_pooled = exp(mean(log(OR), na.rm = TRUE)),
        OR_lower_pooled = exp(mean(log(OR_lower), na.rm = TRUE)),
        OR_upper_pooled = exp(mean(log(OR_upper), na.rm = TRUE)),
        P_pooled = median(P_value, na.rm = TRUE),
        P_interaction_pooled = median(P_interaction, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Metal = metal,
        Subgroup = subgroup_name
      )
    
    subgroup_results_all[[paste(metal, subgroup_name, sep = "_")]] <- results_pooled
  }
}

## ═══════════════════════════════════════════════════════
## 步骤6: 保存结果
## ═══════════════════════════════════════════════════════

cat("\n【步骤6】保存亚组分析结果\n")

subgroup_results_df <- bind_rows(subgroup_results_all)

write.csv(subgroup_results_df,
          here("outputs", "EF", "A2_Subgroup_Analysis_Results.csv"),
          row.names = FALSE)

cat("  ✓ 结果已保存: outputs/EF/A2_Subgroup_Analysis_Results.csv\n")

## ═══════════════════════════════════════════════════════
## 步骤7: 修正版亚组森林图（确保文件可打开 & 动态截断）
## ═══════════════════════════════════════════════════════

cat("\n【步骤7】正在安全生成森林图...\n")

# 强制关闭之前可能遗留的所有图形设备，确保文件可写入
while (!is.null(dev.list())) dev.off()

# 定义截断界限 (OR 尺度)
lower_limit <- 0.1
upper_limit <- 5.0 

subgroup_labels <- c(
  "EF_baseline_group" = "EF基线",
  "EGFR_group"        = "肾小球滤过率 (eGFR)",
  "DM_group"          = "糖尿病",
  "Age_group"         = "年龄",
  "pPCI_group"        = "急诊PCI"
)

for (metal in metals_to_analyze) {
  
  cat("  正在处理金属数据并写入文件:", metal, "\n")
  
  # 初始化绘图表数据 (确保 mean, lower, upper 为 numeric)
  table_data <- data.frame(
    label = "亚组", n = "样本量", or_ci = "OR (95% CI)",
    p_val = "P值", p_int = "交互P值",
    mean = as.numeric(NA), lower = as.numeric(NA), upper = as.numeric(NA),
    is_summary = TRUE, stringsAsFactors = FALSE
  )
  
  metal_results <- subgroup_results_df %>% filter(Metal == metal)
  
  for (sg_var in names(subgroup_labels)) {
    sg_res <- metal_results %>% filter(Subgroup == sg_var)
    if (nrow(sg_res) == 0) next
    
    # 插入标题行
    p_int_val <- if(!is.na(sg_res$P_interaction_pooled[1])) sprintf("%.3f", sg_res$P_interaction_pooled[1]) else ""
    table_data <- table_data %>% add_row(
      label = subgroup_labels[sg_var], p_int = p_int_val, is_summary = TRUE,
      mean = as.numeric(NA), lower = as.numeric(NA), upper = as.numeric(NA)
    )
    
    # 插入数据行
    for (i in 1:nrow(sg_res)) {
      orig_mean  <- sg_res$OR_pooled[i]
      orig_lower <- sg_res$OR_lower_pooled[i]
      orig_upper <- sg_res$OR_upper_pooled[i]
      
      # 处理绘图截断数值 (仅用于绘图，不影响文本显示)
      plot_mean  <- orig_mean
      plot_lower <- max(orig_lower, lower_limit * 0.9) # 略微出界以触发箭头显示
      plot_upper <- min(orig_upper, upper_limit * 1.1)
      
      table_data <- table_data %>% add_row(
        label = paste0("    ", sg_res$Subgroup_level[i]),
        n     = as.character(sg_res$N[i]),
        or_ci = sprintf("%.2f (%.2f-%.2f)", orig_mean, orig_lower, orig_upper),
        p_val = ifelse(sg_res$P_pooled[i] < 0.001, "<0.001", sprintf("%.3f", sg_res$P_pooled[i])),
        p_int = "",
        mean  = plot_mean, lower = plot_lower, upper = plot_upper,
        is_summary = FALSE
      )
    }
  }
  
  # 【核心修正点】安全清洗文本列，不触碰数值列
  table_text_data <- table_data %>%
    mutate(across(c(label, n, or_ci, p_val, p_int), ~replace_na(as.character(.), "")))
  
  table_text_list <- list(
    table_text_data$label, table_text_data$n, 
    table_text_data$or_ci, table_text_data$p_val, table_text_data$p_int
  )
  
  # 设置输出路径
  output_path <- here("plots", "EF", paste0("A2_Forest_Final_", metal, ".png"))
  
  # 打开图形设备
  png(filename = output_path, width = 3600, height = 2400, res = 300)
  
  # 绘图逻辑
  tryCatch({
    p <- forestplot(
      labeltext   = table_text_list,
      mean        = table_data$mean,
      lower       = table_data$lower,
      upper       = table_data$upper,
      is.summary = table_data$is_summary,# 这里控制哪些行加粗
      
      # --- 修改 1：进一步压缩中间图形区域宽度，使其更窄 ---
      graphwidth = unit(75, "mm"), # 原为 25mm，调整为 20mm 以减小占比
      graph.pos  = 3,              # 确保图形放置在 N 和 OR (95% CI) 列之间
      
      xlog       = TRUE,
      zero       = 1.0,
      
      # --- 修改 2：动态设置横轴刻度 (xticks) 匹配实际数值范围，优化下方布局 ---
      # 显式指定刻度，使其在窄图形中分布更均匀，符合实际截断值 lower_limit 和 upper_limit
      xticks     = c(0.1, 0.5, 1.0, 2.0, 5.0), 
      clip       = c(lower_limit, upper_limit),
      
      boxsize    = 0.15,
      lineheight = unit(10, 'mm'),
      colgap     = unit(8, 'mm'),
      lwd.ci     = 2.0,
      col        = fpColors(box = '#458B00', lines = 'black', zero = '#7AC5CD'),
      
      xlab       = "Odds Ratio (OR) on Log Scale",
      title      = paste0("亚组分析: ", metal, " & LVR"),
      
      # --- 修改 3：大幅调大文字大小以增强可读性 ---
      txt_gp     = fpTxtGp(
        label = gpar(cex = 2.3),  # 调大左侧和右侧表格文字
        ticks = gpar(cex = 2.0),  # 调大下方刻度数字
        xlab  = gpar(cex = 2.3),  # 调大横轴标题
        title = gpar(cex = 3.0, fontface = "bold") # 调大顶部标题
      )
    )
    print(p)
    
    # 底部说明
    grid.text("注: 调整18个协变量; 坐标轴已动态截断 (范围: 0.1 - 5.0); * P < 0.05", 
              x = unit(0.05, "npc"), y = unit(0.03, "npc"), 
              just = "left", gp = gpar(cex = 1.8, fontface = "italic"))
    
  }, error = function(e) {
    cat("    绘制过程中出错:", e$message, "\n")
  })
  
  # 强制关闭设备，确保写入磁盘
  dev.off()
  cat("  ✓ 图片已保存:", output_path, "\n")
}
cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║              A2 亚组分析任务全部完成！              ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")