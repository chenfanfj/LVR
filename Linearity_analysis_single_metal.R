############################################################
## 单金属主效应分析 - 步骤1: Logistic回归 (新 LVR 定义 + PS 调整)
## 输入:  outcome_analysis_data_with_weights_extended.rds
## 输出: 19种金属的关联结果 (OR, AUC) + 森林图
## 修改: 
##   1. LVR 定义: 基线 EF >= 50 且 随访 EF_fu < 50
##   2. 采用倾向性评分 (PS) 处理小样本 (Event=32) 协变量调整
##   3. 新增 C-statistic (AUC) 输出
############################################################

rm(list = ls())
gc()

library(dplyr)
library(survey)
library(mitools)
library(ggplot2)
library(here)
library(broom)
library(pROC) # 用于计算 C-statistic (AUC)

# 绘图包
if(!require(forestploter)) install.packages("forestploter"); library(forestploter)
if(!require(grid)) library(grid)
if(!require(showtext)) install.packages("showtext"); library(showtext)

showtext_auto()

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║     单金属主效应分析 - Logistic回归 (PS调整 & AUC)     ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 1: 加载数据、过滤样本并定义 LVR
## ═══════════════════════════════════════════════════════

cat("【步骤 1】数据预处理\n")

# 加载数据 (假设路径不变)
outcome_data <- readRDS(here::here("outputs", "4 IPW2_complete", "outcome_analysis_data_with_weights_extended.rds"))

# 【关键修改】
# 1. 过滤基线 EF 正常的患者 (EF >= 50)
# 2. 定义二分类结局 LVR: 随访 EF_fu < 50
outcome_data_filtered <- outcome_data %>%
  filter(EF_baseline >= 50) %>%
  mutate(
    LVR = ifelse(EF_fu < 50, 1, 0)
  )

n_imps <- length(unique(outcome_data_filtered$.imp))

# 打印新定义下的样本量情况
cat("  基线正常样本量 (EF >= 50):", nrow(outcome_data_filtered %>% filter(.imp == 1)), "\n")
cat("  LVR 事件数 (EF_fu < 50):", sum((outcome_data_filtered %>% filter(.imp == 1))$LVR), "\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 2: 定义金属列表和用于 PS 的协变量
## ═══════════════════════════════════════════════════════

metals <- c("Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
            "Mo", "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li")

# 待压缩进入 PS 的全量协变量
covariates_to_ps <- c("age", "gender", "EF_baseline", "LVEDV_baseline", "cTnIpeak", 
                      "pPCI", "STEMI", "smoking", "DM", "hypertension", 
                      "NTproBNP_peak", "GRACE_in", "WBC", "HGB", "CRP", 
                      "CHOL", "LDL", "AST")

## ═══════════════════════════════════════════════════════
## 步骤 4: 拟合 Logistic 模型 (PS 压缩法 + AUC 计算)
## ═══════════════════════════════════════════════════════

cat("【步骤 4】拟合模型并计算 C-statistic\n")

results_final <- data.frame()

# 开始循环分析 19 种金属
for(i in 1:length(metals)) {
  metal <- metals[i]
  log_metal <- paste0("log_", metal)
  
  if(i %% 5 == 0) cat(i, "..  ")
  
  fit_list <- list()
  auc_vector <- c()
  
  for(imp_i in 1:n_imps) {
    # 提取当前插补集
    dat_i <- outcome_data_filtered %>% filter(.imp == imp_i)
    
    # --- 4.1 计算 PS (将协变量压缩为单一得分) ---
    ps_formula <- as.formula(paste0(log_metal, " ~ ", paste(covariates_to_ps, collapse = " + ")))
    ps_model <- lm(ps_formula, data = dat_i)
    dat_i$ps_score <- predict(ps_model)
    
    # --- 4.2 建立加权设计 ---
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    # --- 4.3 拟合 Logistic 回归 (仅调整金属和 PS) ---
    fit_i <- tryCatch({
      svyglm(as.formula(paste0("LVR ~ ", log_metal, " + ps_score")), 
             design = design_i, family = quasibinomial())
    }, error = function(e) NULL)
    
    if(!is.null(fit_i)) {
      fit_list[[length(fit_list) + 1]] <- fit_i
      
      # --- 4.4 计算 C-statistic (AUC) ---
      probs <- as.numeric(predict(fit_i, type = "response"))
      roc_obj <- roc(dat_i$LVR, probs, weights = dat_i$sw_trunc, quiet = TRUE)
      auc_vector <- c(auc_vector, as.numeric(auc(roc_obj)))
    }
  }
  
  # --- 4.5 汇总 MI 结果 ---
  if(length(fit_list) > 0) {
    pooled_fit <- MIcombine(fit_list)
    summary_fit <- summary(pooled_fit)
    
    if(log_metal %in% rownames(summary_fit)) {
      beta <- summary_fit[log_metal, "results"]
      se <- summary_fit[log_metal, "se"]
      
      results_final <- rbind(results_final, data.frame(
        Metal = metal,
        OR = exp(beta),
        Lower95 = exp(beta - 1.96 * se),
        Upper95 = exp(beta + 1.96 * se),
        P_value = 2 * (1 - pnorm(abs(beta / se))),
        C_statistic = mean(auc_vector) # 输出插补集平均 AUC
      ))
    }
  }
}

cat("完成!\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 5: 保存结果 (含 AUC)
## ═══════════════════════════════════════════════════════

results_final <- results_final %>%
  mutate(P_FDR = p.adjust(P_value, method = "BH")) %>%
  arrange(P_value)

write.csv(results_final,
          here::here("outputs", "EF", "Single_Metal_PS_Adjusted_LVR_with_AUC.csv"),
          row.names = FALSE)


## ═══════════════════════════════════════════════════════
## 步骤 6: 最终视觉优化版 - 图形减半，字体加粗加大
## ═══════════════════════════════════════════════════════

if(!require(forestplot)) install.packages("forestplot"); library(forestplot)

cat("【步骤 6】生成图形压缩、字体放大版森林图...\n")

# 1. 路径准备
target_dir <- here::here("plots", "EF")
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
file_path <- here::here("plots", "EF", "Forest_Plot_Optimized_Visual.png")
metal_map <- c(
  "Cu" = "铜 (Cu)", "Zn" = "锌 (Zn)", "Fe" = "铁 (Fe)", "Se" = "硒 (Se)", 
  "Pb" = "铅 (Pb)", "Al" = "铝 (Al)", "As" = "砷 (As)", "Cr" = "铬 (Cr)", 
  "Mn" = "锰 (Mn)", "Ni" = "镍 (Ni)", "Mo" = "钼 (Mo)", "Rb" = "铷 (Rb)", 
  "Sb" = "锑 (Sb)", "Sn" = "锡 (Sn)", "Sr" = "锶 (Sr)", "V"  = "钒 (V)", 
  "Ba" = "钡 (Ba)", "B"  = "硼 (B)",  "Li" = "锂 (Li)"
)
# 2. 数据准备 (按照 P 值排序)
fp_data <- results_final %>%
  arrange(P_value) %>%
  mutate(
    Metal_Lab = metal_map[as.character(Metal)],
    OR_CI = sprintf("%.2f (%.2f, %.2f)", OR, Lower95, Upper95),
    AUC_val = sprintf("%.3f", C_statistic),
    P_val = format.pval(P_value, digits = 3, eps = 0.001)
  )

# 3. 构造文本矩阵
table_text <- rbind(
  c("金属", "OR (95% CI)", "C-统计量", "P值"),
  as.matrix(fp_data %>% select(Metal_Lab, OR_CI, AUC_val, P_val))
)

# 4. 构造数值向量
mean_vec  <- c(NA, fp_data$OR)
lower_vec <- c(NA, fp_data$Lower95)
upper_vec <- c(NA, fp_data$Upper95)

# 5. 绘图与保存
# 画布宽度稍微加宽到 7.5 英寸，以容纳更大的字体，高度保持紧凑
png(file_path, width = 7.5, height = 6, units = "in", res = 300)

forestplot(labeltext = table_text,
           mean  = mean_vec,
           lower = lower_vec,
           upper = upper_vec,
           
           # --- 关键修改 1：压缩中间图形区域 ---
           # 将图形宽度设为 25-30mm，能显著减小图形在全图中的占比
           graphwidth = unit(28, "mm"), 
           graph.pos = 2,
           
           # --- 关键修改 2：全面放大字体 (cex) ---
           txt_gp = fpTxtGp(
             label = gpar(cex = 1.2),     # 侧边表格文字显著放大
             ticks = gpar(cex = 1.0),     # 刻度数值放大
             xlab  = gpar(cex = 1.1),     # X 轴标签放大
             title = gpar(cex = 1.3, fontface = "bold") # 标题加粗加大
           ),
           
           # --- 样式与间距优化 ---
           is.summary = c(TRUE, rep(FALSE, nrow(fp_data))),
           zero = 1,
           boxsize = 0.15,
           lineheight = unit(6, 'mm'),    # 配合大字体，稍微增加行高
           colgap = unit(5, 'mm'),        # 增加列间距，使大字体不拥挤
           lwd.zero = 1.5,
           lwd.ci = 1.8, 
           col = fpColors(box = '#458B00', lines = 'black', zero = '#7AC5CD'),
           
           # --- 辅助设置 ---
           xlab = "Odds Ratio (OR)",
           clip = c(0.1, 3.5), 
           lty.ci = "solid",
           title = "金属暴露与 LVR 风险关联",
           line.margin = 0.05,
           fn.get_labels = "注: 调整18个协变量; * P < 0.05")

dev.off()

cat(paste0("  ✓ 最终优化版森林图已保存: ", file_path, "\n"))
