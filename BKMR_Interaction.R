rm(list = ls())
gc()

library(bkmr)
library(ggplot2)
library(dplyr)
library(here)

# 1. 加载模型对象
cat("正在读取 BKMR 模型对象...\n")
bkmr_fit <- readRDS(here("outputs", "EF", "A3_BKMR_Fit_Object.rds"))

# 2. 准备基础数据 (必须从原始模型中提取 Z 矩阵以获取分位数)
Z <- bkmr_fit$Z 
sel_idx <- seq(25001, 50000, by = 1) # 保持与之前一致的 Burn-in

# 定义分析组合
inter_pairs <- list(c("Zn", "Cu"), c("Zn", "Pb"), c("Cu", "Pb"))

## ═══════════════════════════════════════════════════════
## 步骤 7.1: 双金属交互作用地形图
## ═══════════════════════════════════════════════════════
for (pair in inter_pairs) {
  m1 <- pair[1]; m2 <- pair[2]
  cat("正在处理地形图:", m1, "&", m2, "\n")
  
  inter_grid <- PredictorResponseBivar(
    fit = bkmr_fit,
    z.ptr = c(which(colnames(Z) == paste0("log_", m1)), 
              which(colnames(Z) == paste0("log_", m2))),
    gridsize = 20, 
    sel = sel_idx
  )
  
  plot_data <- as.data.frame(as.matrix(inter_grid))
  plot_data <- data.frame(
    x_orig = exp(as.numeric(plot_data[, 3])),
    y_orig = exp(as.numeric(plot_data[, 4])),
    z_val  = as.numeric(plot_data[, 5])
  )
  
  p_contour <- ggplot(plot_data, aes(x = x_orig, y = y_orig)) +
    geom_tile(aes(fill = z_val)) + 
    geom_contour(aes(z = z_val), color = "white", linewidth = 1.0, bins = 5) +
    scale_fill_viridis_c(option = "magma", name = "h(z)") +
    scale_x_log10() + scale_y_log10() +
    labs(title = paste(m1, "与", m2, "交互作用地形图"),
         x = paste(m1, "浓度"), y = paste(m2, "浓度")) +
    theme_bw()
  
  ggsave(here("plots", "EF", paste0("A4_BKMR_Contour_", m1, "_", m2, ".png")), p_contour)
  
  # 及时清理大对象释放内存
  rm(inter_grid, plot_data, p_contour); gc() 
}

## ═══════════════════════════════════════════════════════
## 步骤 7.2: 分组单效应图
## ═══════════════════════════════════════════════════════
cat("\n正在处理分组单效应图...\n")

modifier_pairs <- list(
  list(exposure = "Zn", modifier = "Pb"),
  list(exposure = "Cu", modifier = "Pb"),
  list(exposure = "Zn", modifier = "Cu")
)

for (comb in modifier_pairs) {
  exp_m <- comb$exposure; mod_m <- comb$modifier
  mod_col <- paste0("log_", mod_m)
  q_vals <- quantile(Z[, mod_col], probs = c(0.10, 0.50, 0.90))
  group_results <- data.frame()
  
  for (i in 1:3) {
    z_fixed <- apply(Z, 2, median); z_fixed[mod_col] <- q_vals[i]
    res <- PredictorResponseUnivar(
      fit = bkmr_fit,
      which.z = which(colnames(Z) == paste0("log_", exp_m)),
      z.fixed = z_fixed, sel = sel_idx
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
  
  # 修正后的分组效应绘图代码（含地毯图与独立密度引用）
  p_group <- ggplot(group_results, aes(x = z_orig, y = est, group = Group, color = Group)) +
    # 1. 置信区间阴影（最底层）
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Group), alpha = 0.05, color = NA) +
    
    # 2. 地毯图：显式引用原始观测数据 Z，以反映真实密度
    # 注意：inherit.aes = FALSE 必须加上，防止它去寻找 Group 变量
    geom_rug(data = data.frame(x_rug = exp(Z[, paste0("log_", exp_m)])), 
             aes(x = x_rug), inherit.aes = FALSE, 
             sides = "b", alpha = 0.2, color = "grey30") +
    
    # 3. 预测线条（中间层）
    # 使用 linetype 区分，即使线条重合也能看出是由三组构成的
    geom_line(aes(linetype = Group), linewidth = 1.2) +
    
    # 4. 辅助线与坐标轴
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_log10() +
    
    # 5. 手动配置颜色与线型
    scale_color_manual(values = c("10th" = "#377eb8", "50th" = "#4daf4a", "90th" = "#e41a1c")) +
    scale_fill_manual(values = c("10th" = "#377eb8", "50th" = "#4daf4a", "90th" = "#e41a1c")) +
    scale_linetype_manual(values = c("10th" = "dotted", "50th" = "solid", "90th" = "dashed")) +
    
    # 6. 标签与主题
    labs(title = paste(exp_m, "在不同", mod_m, "分位数下的效应曲线"),
         subtitle = "底部刻度代表原始样本分布；阴影随样本稀疏而变宽是正常现象",
         x = paste(exp_m, "浓度 (μg/L)"), 
         y = "h(z) - 风险贡献",
         linetype = "调节因子水平") +
    theme_bw(base_size = 12) + 
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank()) # 关闭细网格让画面更干净
  
  ggsave(here("plots", "EF", paste0("A4_BKMR_GroupEffect_", exp_m, "_by_", mod_m, ".png")), p_group)
  rm(group_results, p_group); gc()
}


## ═══════════════════════════════════════════════════════
## 步骤 7.2.2: 分组单效应图
## ═══════════════════════════════════════════════════════
library(patchwork)

cat("\n【步骤 7.2.2】生成分组效应汇总组合图\n")

# 1. 初始化一个列表用于存储图形
all_plots <- list()

for (comb in modifier_pairs) {
  exp_m <- comb$exposure
  mod_m <- comb$modifier
  cat("  正在处理组合:", exp_m, "by", mod_m, "\n")
  
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
  
  # 绘制单个图形
  p_group <- ggplot(group_results, aes(x = z_orig, y = est, group = Group, color = Group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Group), alpha = 0.05, color = NA) +
    geom_rug(data = data.frame(x_rug = exp(Z[, paste0("log_", exp_m)])), 
             aes(x = x_rug), inherit.aes = FALSE, 
             sides = "b", alpha = 0.2, color = "grey30") +
    geom_line(aes(linetype = Group), linewidth = 1.0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_x_log10() +
    scale_color_manual(values = c("10th" = "#377eb8", "50th" = "#4daf4a", "90th" = "#e41a1c")) +
    scale_fill_manual(values = c("10th" = "#377eb8", "50th" = "#4daf4a", "90th" = "#e41a1c")) +
    scale_linetype_manual(values = c("10th" = "dotted", "50th" = "solid", "90th" = "dashed")) +
    labs(title = paste(exp_m, "(受", mod_m, "调节)"),
         x = paste(exp_m, "(μg/L)"), 
         y = "h(z) - 风险贡献") +
    theme_bw(base_size = 10) + 
    theme(legend.position = "none", # 在子图中先隐藏图例，最后统一显示
          panel.grid.minor = element_blank())
  
  # 将图形存入列表
  all_plots[[paste0(exp_m, "_", mod_m)]] <- p_group
}

# 2. 使用 patchwork 组合图形
# ncol = 2 表示两列排列，guides = "collect" 表示合并公共图例
combined_plot <- wrap_plots(all_plots, ncol = 2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "多金属交互作用分组效应汇总分析图",
    subtitle = "各曲线重合且底部地毯图反映样本分布；阴影代表 95% 置信区间",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5))
  ) & theme(legend.position = "bottom")

# 3. 保存最终组合大图
ggsave(here("plots", "EF", "A4_BKMR_Combined_GroupEffects0.png"), 
       combined_plot, width = 12, height = 10, dpi = 300, bg = "white")

cat("\n✓ 组合图已保存至: plots/EF/A4_BKMR_Combined_GroupEffects.png\n")
## ═══════════════════════════════════════════════════════
## 步骤 7.3: 交互作用数值验证 
## ═══════════════════════════════════════════════════════
cat("\n【补充】交互作用数值验证 (判断曲线是否重合)\n")

# 定义需要验证的组合
check_pairs <- list(
  list(exposure = "Zn", modifier = "Pb"),
  list(exposure = "Cu", modifier = "Pb"),
  list(exposure = "Zn", modifier = "Cu")
)

# 准备结果容器
interaction_summary <- data.frame()

for (comb in check_pairs) {
  exp_m <- comb$exposure; mod_m <- comb$modifier
  cat("  正在验证:", exp_m, "是否受", mod_m, "调节...\n")
  
  mod_col <- paste0("log_", mod_m)
  q_vals <- quantile(Z[, mod_col], probs = c(0.10, 0.90))
  
  # 1. 计算 10th 分位数下的曲线值
  z_fixed_10 <- apply(Z, 2, median); z_fixed_10[mod_col] <- q_vals[1]
  res_10 <- PredictorResponseUnivar(
    fit = bkmr_fit, 
    which.z = which(colnames(Z) == paste0("log_", exp_m)), 
    z.fixed = z_fixed_10, 
    sel = sel_idx
  )
  
  # 2. 计算 90th 分位数下的曲线值
  z_fixed_90 <- apply(Z, 2, median); z_fixed_90[mod_col] <- q_vals[2]
  res_90 <- PredictorResponseUnivar(
    fit = bkmr_fit, 
    which.z = which(colnames(Z) == paste0("log_", exp_m)), 
    z.fixed = z_fixed_90, 
    sel = sel_idx
  )
  
  # 3. 提取 est 列进行数值对比
  est_10 <- as.numeric(as.matrix(res_10)[, 3])
  est_90 <- as.numeric(as.matrix(res_90)[, 3])
  
  # 4. 计算差异指标
  diff_vec <- est_10 - est_90
  max_diff <- max(abs(diff_vec), na.rm = TRUE)
  mean_diff <- mean(abs(diff_vec), na.rm = TRUE)
  
  # 5. 保存结果
  interaction_summary <- rbind(interaction_summary, data.frame(
    Exposure = exp_m,
    Modifier = mod_m,
    Max_Abs_Diff = max_diff,
    Mean_Abs_Diff = mean_diff,
    Conclusion = ifelse(max_diff < 0.01, "No Interaction (Overlapped)", "Potential Interaction")
  ))
}

# 打印最终对比表
cat("\n--- 交互作用数值对比总结 ---\n")
print(interaction_summary)

# 保存数值结果
write.csv(interaction_summary, here("outputs", "EF", "A4_BKMR_Interaction_Numerical_Check.csv"), row.names = FALSE)

# 检查 r 参数的汇总统计
#在 BKMR 中，参数r_m决定了金属m的影响强度。如果两个金属的r值都较大且模型收敛，它们共同进入核函数。
#虽然这不直接等于交互 PIP，但 $r$ 的波动情况能反映交互的潜力。
r_summary <- summary(bkmr_fit$r[sel_idx, ])
print(r_summary)

# 计算每个变量的 PIP (即 r > 0 的比例)
r_pips <- colMeans(bkmr_fit$r[sel_idx, ] > 0)
names(r_pips) <- colnames(Z)
print(r_pips)
