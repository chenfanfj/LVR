############################################################
##  A3.3: Bayesian Kernel Machine Regression (BKMR)
##  输入: imputed_data_with_ipw_weights_extended.rds
##  输出: 单金属效应、联合效应、交互作用图
##  ⚠️ 计算密集，运行时间长 (可能数小时)
############################################################

rm(list = ls())
gc()

library(bkmr)  # BKMR包
library(ggplot2)
library(dplyr)
library(here)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║           A3.3: BKMR模型 (多金属交互)              ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("⚠️  警告:  BKMR计算量大，建议使用高性能计算集群\n")
cat("   本示例使用较少迭代次数 (iter=5000) 以加速\n")
cat("   正式分析建议:  iter=50000, burnin=25000\n\n")

## ═══════════════════════════════════════════════════════
## 步骤1: 读取数据
## ═══════════════════════════════════════════════════════

cat("【步骤1】读取数据\n")

imputed_data <- readRDS(here("outputs", "4 IPW2_complete", "imputed_data_with_ipw_weights_extended.rds"))

outcome_data <- imputed_data %>% 
  filter(has_fu_echo == "Yes", !is.na(LVEDV_fu), .imp == 1) %>%
  mutate(LVR = as.numeric(EF_baseline >= 50 & EF_fu < 50))

cat("  样本数:", nrow(outcome_data), "\n\n")

## ═══════════════════════════════════════════════════════
## 步骤2: 准备BKMR数据
## ═══════════════════════════════════════════════════════

metals_bkmr <- c("log_Pb", "log_Zn", "log_Fe", "log_Cu", "log_Se", "log_As")
covariates <- c("age", "gender", "EF_baseline", "LVEDV_baseline", "cTnIpeak")

bkmr_data <- outcome_data %>%
  select(LVR, all_of(metals_bkmr), all_of(covariates)) %>%
  na.omit()

# 构建矩阵
Y <- bkmr_data$LVR
Z <- as.matrix(bkmr_data[, metals_bkmr])
X <- as.matrix(bkmr_data[, covariates])

cat("【步骤3】拟合BKMR模型\n")
cat("  迭代次数:  5000 (演示), 建议50000+\n")
cat("  预计耗时: 10-30分钟 (取决于CPU)\n\n")

## ═══════════════════════════════════════════════════════
## 步骤3: 拟合BKMR
## ═══════════════════════════════════════════════════════

cat("【步骤3】拟合BKMR模型 \n")

set.seed(2024)

bkmr_fit <- kmbayes(
  y = Y,
  Z = Z,
  X = X,
  iter = 50000,      # 总迭代次数
  verbose = TRUE,
  varsel = TRUE,     # 既然增加到5万次，可以尝试保持变量选择开启
  family = "binomial"
)

cat("\n  ✓ BKMR拟合完成\n\n")

# 定义选择范围：从 25001 到 50000 次
sel_idx <- seq(25001, 50000, by = 1)


## ═══════════════════════════════════════════════════════
## 步骤4: 后验包含概率 (PIP)
## ═══════════════════════════════════════════════════════
cat("【步骤4】变量重要性 (PIP)\n")

# 1. 获取 PIP 对象
pip_raw <- ExtractPIPs(bkmr_fit, sel = sel_idx)

# 2. 定向提取数值部分
# 在 varsel=TRUE 时，bkmr 通常将具体变量的 PIP 放在第二个元素或名为 'selPips' 的部分
if (is.list(pip_raw)) {
  # 优先尝试通过名字提取，如果名字不存在则取最后一个元素（通常是 selPips）
  if ("selPips" %in% names(pip_raw)) {
    pip_values <- pip_raw$selPips
  } else {
    pip_values <- pip_raw[[length(pip_raw)]] 
  }
} else {
  pip_values <- pip_raw
}

# 3. 确保提取出来的是纯数值向量
pip_values <- as.numeric(pip_values)

# 4. 构建数据框并排序
pip_df <- data.frame(
  Metal = gsub("log_", "", colnames(Z)),
  PIP = pip_values
) %>% 
  arrange(desc(PIP))

cat("\n后验包含概率 (PIP):\n")
print(pip_df)

cat("\n  解释:  PIP>0.5表示该金属对LVR有重要影响\n\n")

## ═══════════════════════════════════════════════════════
## 步骤5: 单金属剂量-反应曲线
## ═══════════════════════════════════════════════════════
cat("【步骤5】绘制单金属效应曲线 - 动态位置索引版\n")
# 在步骤 5 开始处添加以下代码
# 选择 PIP > 0.3 的金属进行绘图，若无则绘制全部
metals_to_plot <- pip_df$Metal[pip_df$PIP > 0.3]
if (length(metals_to_plot) == 0) metals_to_plot <- pip_df$Metal

for (metal in metals_to_plot) {
  target_col <- paste0("log_", metal)
  cat("  正在处理:", metal, "\n")
  
  pred_response <- PredictorResponseUnivar(
    fit = bkmr_fit,
    which.z = which(colnames(Z) == target_col),
    sel = sel_idx # 关键
  )
  
  # 1. 强制转为纯数值矩阵，避免因子干扰
  res_mat <- as.matrix(pred_response)
  
  # 2. 识别真实数据列
  # 在大多数版本中：
  # 第1列是变量名(字符)，需跳过
  # 第2列是 z (浓度分位数)
  # 第3列是 est (效应值)
  # 倒数第2列和最后1列是 lower 和 upper
  n_total <- ncol(res_mat)
  
  plot_df <- data.frame(
    m_val = as.numeric(res_mat[, 2]),          # 尝试第2列作为浓度
    est   = as.numeric(res_mat[, 3]),          # 尝试第3列作为效应
    lower = as.numeric(res_mat[, n_total - 1]),# 倒数第2列
    upper = as.numeric(res_mat[, n_total])    # 最后1列
  )
  
  # 3. 容错检查：如果第2列全是NA，说明数据可能从第1列就开始了
  if(all(is.na(plot_df$m_val))) {
    plot_df$m_val <- as.numeric(res_mat[, 1])
    plot_df$est   <- as.numeric(res_mat[, 2])
  }
  
  # 4. 还原浓度并清洗
  plot_df <- plot_df[complete.cases(plot_df), ]
  plot_df$metal_orig <- exp(plot_df$m_val)
  
  # 5. 绘图 (保持之前的 ggplot 代码)
  p <- ggplot(plot_df, aes(x = metal_orig, y = est)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "darkred") +
    geom_line(color = "darkred", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_log10() +
    labs(title = paste(metal, "单金属效应 (BKMR)"),
         x = paste0(metal, " 浓度 (μg/L)"), y = "h(浓度) - 相对风险") +
    theme_minimal(base_size = 24)
  
  ggsave(filename = here("plots", "EF", paste0("A3_BKMR_UnivarEffect_", metal, ".png")),
         plot = p, width = 8, height = 6)
}
cat("  ✓ 单金属效应图已保存\n")

## ═══════════════════════════════════════════════════════
## 步骤5.1: 单金属剂量-反应曲线
## ═══════════════════════════════════════════════════════
cat("\n【补充】检查模型收敛性 (Trace Plots)\n")

# 确保目录存在
if(!dir.exists(here("plots", "EF", "Diagnostics"))) dir.create(here("plots", "EF", "Diagnostics"), recursive = TRUE)

# 1. 检查误差方差 sigsq.eps (反映整体拟合稳定性)
png(here("plots", "EF", "Diagnostics", "Trace_Sigma.png"), width = 800, height = 600)
TracePlot(fit = bkmr_fit, par = "sigsq.eps")
dev.off()

# 2. 检查协变量系数 beta (反映 age, gender 等校正是否稳定)
png(here("plots", "EF", "Diagnostics", "Trace_Beta.png"), width = 800, height = 600)
TracePlot(fit = bkmr_fit, par = "beta")
dev.off()

# 3. 检查变量选择参数 r (反映各金属的重要性选择过程)
png(here("plots", "EF", "Diagnostics", "Trace_r.png"), width = 800, height = 600)
TracePlot(fit = bkmr_fit, par = "r")
dev.off()

cat("  ✓ 收敛性诊断图已保存至: plots/EF/Diagnostics/\n")

## ═══════════════════════════════════════════════════════
## 步骤6: 多金属联合效应
## ═══════════════════════════════════════════════════════
cat("\n【步骤6】计算多金属联合效应 - 强力重命名版\n")

overall_effect <- OverallRiskSummaries(
  fit = bkmr_fit,
  qs = seq(0.25, 0.75, by = 0.05),
  sel = sel_idx
)

overall_df <- as.data.frame(as.matrix(overall_effect))
n_cols <- ncol(overall_df)

# 根据列数执行强制重命名
if (n_cols == 5) {
  colnames(overall_df) <- c("quantile", "est", "sd", "lower", "upper")
} else if (n_cols == 4) {
  colnames(overall_df) <- c("quantile", "est", "lower", "upper")
} else {
  colnames(overall_df)[1:2] <- c("quantile", "est")
}

# 关键：彻底数值化
overall_df[] <- lapply(overall_df, function(x) as.numeric(as.character(x)))

if ("lower" %in% colnames(overall_df)) {
  p_overall <- ggplot(overall_df, aes(x = quantile, y = est)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgreen") +
    geom_line(color = "darkgreen", linewidth = 1.2) +
    geom_point(color = "darkgreen", size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "多金属联合效应 (95% CI)", x = "金属共同分位数", y = "h(z)") +
    theme_minimal(base_size = 24)
} else {
  p_overall <- ggplot(overall_df, aes(x = quantile, y = est)) +
    geom_line(color = "darkgreen", linewidth = 1.2, linetype = "dotted") +
    geom_point(color = "darkgreen", size = 3) +
    labs(title = "多金属联合效应 (无CI趋势图)") +
    theme_minimal()
}

ggsave(here("plots", "EF", "A3_BKMR_Overall_Effect.png"), p_overall, width = 8, height = 6)
cat("  ✓ 联合效应图已保存\n")


## ═══════════════════════════════════════════════════════
## 步骤 7.1: 双金属交互作用等高线图 (方案 A)
## ═══════════════════════════════════════════════════════
cat("\n【步骤 7.1】绘制双金属交互作用地形图 - 轻量化版\n")

inter_pairs <- list(c("Zn", "Cu"), c("Zn", "Pb"), c("Cu", "Pb"))

for (pair in inter_pairs) {
  m1 <- pair[1]; m2 <- pair[2]
  cat("  正在处理:", m1, "&", m2, "\n")
  
  inter_grid <- PredictorResponseBivar(
    fit = bkmr_fit,
    z.ptr = c(which(colnames(Z) == paste0("log_", m1)), 
              which(colnames(Z) == paste0("log_", m2))),
    gridsize = 20, # 稍微降低网格密度以节省内存
    sel = sel_idx
  )
  
  plot_data <- as.data.frame(as.matrix(inter_grid))
  # 强制根据位置提取：3:z1, 4:z2, 5:est
  plot_data <- data.frame(
    x_orig = exp(as.numeric(plot_data[, 3])),
    y_orig = exp(as.numeric(plot_data[, 4])),
    z_val  = as.numeric(plot_data[, 5])
  )
  
  # 剔除异常值
  plot_data <- plot_data[is.finite(plot_data$x_orig) & is.finite(plot_data$y_orig), ]
  
  # 绘图：使用 interpolate = TRUE 让颜色平滑
  p_contour <- ggplot(plot_data, aes(x = x_orig, y = y_orig)) +
    geom_tile(aes(fill = z_val)) + 
    # 显著加粗等高线，并减少条数（bins=5）使其更清晰
    geom_contour(aes(z = z_val), color = "white", linewidth = 1.2, bins = 5) +
    scale_fill_viridis_c(option = "magma", name = "h(z)") + # 换一个对比度更高的配色
    scale_x_log10() + scale_y_log10() +
    labs(title = paste(m1, "与", m2, "交互作用地形图"),
         subtitle = "白色粗线为风险等高线",
         x = paste(m1, "浓度"), y = paste(m2, "浓度")) +
    theme_bw(base_size = 18)
  
  ggsave(filename = here("plots", "EF", paste0("A4_BKMR_Contour_", m1, "_", m2, ".png")),
         plot = p_contour, width = 8, height = 7, dpi = 300)
}


## ═══════════════════════════════════════════════════════
## 步骤 7.2: 分组单效应图 (方案 B)
## ═══════════════════════════════════════════════════════
cat("\n【步骤 7.2】绘制分组单效应图 - 强化显示版\n")

modifier_pairs <- list(
  list(exposure = "Zn", modifier = "Pb"),
  list(exposure = "Cu", modifier = "Pb"),
  list(exposure = "Zn", modifier = "Cu")
)

for (comb in modifier_pairs) {
  exp_m <- comb$exposure; mod_m <- comb$modifier
  cat("  分析:", exp_m, "在不同", mod_m, "分位数下的效应\n")
  
  mod_col <- paste0("log_", mod_m)
  q_vals <- quantile(Z[, mod_col], probs = c(0.10, 0.50, 0.90))
  group_results <- data.frame()
  
  for (i in 1:3) {
    z_fixed <- apply(Z, 2, median)
    z_fixed[mod_col] <- q_vals[i]
    
    res <- PredictorResponseUnivar(
      fit = bkmr_fit,
      which.z = which(colnames(Z) == paste0("log_", exp_m)),
      z.fixed = z_fixed,
      sel = sel_idx
    )
    
    tmp_mat <- as.matrix(res)
    tmp_clean <- data.frame(
      z_orig = exp(as.numeric(tmp_mat[, 2])),
      est    = as.numeric(tmp_mat[, 3]),
      lower  = as.numeric(tmp_mat[, ncol(tmp_mat)-1]),
      upper  = as.numeric(tmp_mat[, ncol(tmp_mat)]),
      Group  = factor(c("10th", "50th", "90th")[i], levels = c("10th", "50th", "90th"))
    )
    group_results <- rbind(group_results, tmp_clean)
  }
  
  # 彻底剔除 NA，防止 ggplot 丢弃整层数据
  group_results <- group_results[complete.cases(group_results), ]
  
  p_group <- ggplot(group_results, aes(x = z_orig, y = est, group = Group, color = Group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Group), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_log10() +
    scale_color_manual(values = c("10th" = "#1f78b4", "50th" = "#33a02c", "90th" = "#e31a1c")) +
    scale_fill_manual(values = c("10th" = "#1f78b4", "50th" = "#33a02c", "90th" = "#e31a1c")) +
    labs(title = paste(exp_m, "剂量-反应曲线 (受", mod_m, "调节)"),
         x = paste(exp_m, "浓度 (μg/L)"), y = "h(z)") +
    theme_minimal(base_size = 36) + 
    theme(legend.position = "bottom")
  
  ggsave(filename = here("plots", "EF", paste0("A4_BKMR_GroupEffect_", exp_m, "_by_", mod_m, ".png")),
         plot = p_group, width = 9, height = 7, dpi = 300, bg = "white")
}

## ═══════════════════════════════════════════════════════
## 步骤8: 保存BKMR对象
## ═══════════════════════════════════════════════════════

saveRDS(bkmr_fit, here("outputs", "EF", "A3_BKMR_Fit_Object.rds"))

write.csv(pip_df, here("outputs", "EF", "A3_BKMR_PIP. csv"), row.names = FALSE)

cat("\n  ✓ BKMR对象已保存:  outputs/EF/A3_BKMR_Fit_Object.rds\n")
cat("  ✓ PIP结果已保存: outputs/EF/A3_BKMR_PIP.csv\n")

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║              A3.3 BKMR模型完成！                    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("\n📊 BKMR结果解释:\n")
cat("  • PIP (后验包含概率): 金属被选入模型的概率\n")
cat("  • 单金属效应: 固定其他金属在中位数时的剂量-反应\n")