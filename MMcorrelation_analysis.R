# ==============================================================================
# 分析1.1: 金属-代谢物相关矩阵 (IPW加权 + 多重插补池化)
# ==============================================================================

rm(list = ls())
gc()

# ------------------------------------------------------------------------------
# 步骤1: 环境准备
# ------------------------------------------------------------------------------
library(mice)
library(dplyr)
library(tidyr)
library(psych)        # 提供 corr.test() 用于相关检验
library(weights)      # 提供 wtd.cor() 用于加权相关
library(pheatmap)     # 热图绘制
library(ggplot2)
library(RColorBrewer)
library(here)         # 确保路径管理

# 创建保存目录(如果不存在)
dir.create("plots/14 MMcorrelation_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/14 MMcorrelation_analysis", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 修正后的步骤2: 智能识别变量 (金属、代谢物、权重)
# ==============================================================================
cat("===== 步骤2: 数据读取与变量智能识别 =====\n")

# 2.1 读取多重插补数据
mids_obj <- readRDS("outputs/mice_imputation_with_metabo.rds")
cat("- 多重插补数据已加载:", mids_obj$m, "个数据集,", nrow(mids_obj$data), "行\n")

# --- 2.2 金属变量识别 (保持之前的修正) ---
target_metals <- c("Li", "Be", "V", "Cr", "Mn", "Co", "Ni", "Cu", "Zn",
                   "As", "Se", "Mo", "Cd", "Sn", "Sb", "Ba", "Tl", "Pb", "U", "Al", "Fe", "Sr")
data_names <- names(mids_obj$data)
found_metals <- intersect(target_metals, data_names)

# 如果没找到原始名，找 log_ 或 ln_ 前缀
if(length(found_metals) == 0) {
  log_candidates <- paste0("log_", target_metals)
  ln_candidates <- paste0("ln_", target_metals)
  if(any(log_candidates %in% data_names)) {
    metal_vars <- data_names[data_names %in% log_candidates]
    cat("- 识别到 log_ 金属变量:", length(metal_vars), "个\n")
  } else if(any(ln_candidates %in% data_names)) {
    metal_vars <- data_names[data_names %in% ln_candidates]
    cat("- 识别到 ln_ 金属变量:", length(metal_vars), "个\n")
  } else {
    stop("❌ 无法识别金属变量，请检查列名！")
  }
} else {
  metal_vars <- found_metals
  # 这里建议手动对数转换逻辑(如上一轮回复所示)，此处从略，假设已处理
}

# --- 2.3 代谢物变量识别 ---
candidate_file <- "outputs/13 candidate_metabolites/final_metabolite_list.rds"
if(file.exists(candidate_file)) {
  metabo_vars <- readRDS(candidate_file)
} else {
  metabo_vars <- names(mids_obj$data)[grep("^M[0-9]+$|metabo_", names(mids_obj$data))]
}
metabo_vars <- metabo_vars[metabo_vars %in% names(mids_obj$data)]

# --- 2.4 权重变量智能搜索 (修复 Issue 1) ---
# 定义可能的权重变量名
possible_weights <- c("ipw", "weight", "weights", "iptw", "sw", "ps_weight")
# 找到第一个匹配的
weight_var <- intersect(possible_weights, names(mids_obj$data))[1]

if(!is.na(weight_var)) {
  use_weights <- TRUE
  cat(sprintf("✅ 成功识别权重变量: '%s'\n", weight_var))
} else {
  use_weights <- FALSE
  weight_var <- NULL
  cat("⚠️ 未在数据中找到权重变量 (如 'ipw', 'weight')。\n")
  cat("   -> 分析将继续进行，但不使用加权 (Unweighted Analysis)。\n")
  cat("   -> 原因: 之前的插补/合并步骤可能未包含权重计算结果。\n")
}

# ------------------------------------------------------------------------------
# 步骤3: 在每个插补数据集中计算加权Spearman相关
# ------------------------------------------------------------------------------
cat("\n===== 步骤3: 计算加权相关矩阵 =====\n")

# 3.1 定义加权Spearman相关函数
weighted_spearman_cor <- function(x, y, weights = NULL) {
  if(is.null(weights)) {
    return(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
  }
  
  complete_cases <- complete.cases(x, y, weights)
  x <- x[complete_cases]
  y <- y[complete_cases]
  weights <- weights[complete_cases]
  
  # 在秩次上计算加权相关
  rank_x <- rank(x)
  rank_y <- rank(y)
  
  # wtd.cor 返回相关矩阵，取 [1,1] 即可
  # 注意 weights 包的 wtd.cor 对单值会报错，需增加 tryCatch
  tryCatch({
    w_cor <- wtd.cor(rank_x, rank_y, weight = weights)[1, 1]
    return(w_cor)
  }, error = function(e) {
    return(NA)
  })
}

# 3.2 计算循环
M <- mids_obj$m
cor_matrices_list <- vector("list", M)
n_samples_list <- vector("list", M)

cat("开始计算...\n")
for(m in 1:M) {
  if(m %% 5 == 0) cat(sprintf("- 处理插补数据集 %d/%d.. .\n", m, M))
  
  data_m <- complete(mids_obj, m)
  
  cor_matrix <- matrix(NA, nrow = length(metal_vars), ncol = length(metabo_vars),
                       dimnames = list(metal_vars, metabo_vars))
  n_matrix <- matrix(NA, nrow = length(metal_vars), ncol = length(metabo_vars),
                     dimnames = list(metal_vars, metabo_vars))
  
  for(i in 1:length(metal_vars)) {
    for(j in 1:length(metabo_vars)) {
      metal <- metal_vars[i]
      metabo <- metabo_vars[j]
      
      x <- data_m[[metal]]
      y <- data_m[[metabo]]
      w <- if(use_weights) data_m[[weight_var]] else NULL
      
      if(use_weights) {
        cor_matrix[i, j] <- weighted_spearman_cor(x, y, w)
      } else {
        cor_matrix[i, j] <- cor(x, y, method = "spearman", use = "pairwise.complete.obs")
      }
      n_matrix[i, j] <- sum(complete.cases(x, y))
    }
  }
  cor_matrices_list[[m]] <- cor_matrix
  n_samples_list[[m]] <- n_matrix
}
cat("✅ 相关矩阵计算完成\n")

# ------------------------------------------------------------------------------
# 步骤4: 使用Fisher's Z转换池化相关系数
# ------------------------------------------------------------------------------
cat("\n===== 步骤4: 池化相关系数 =====\n")

fisher_z <- function(r) { 0.5 * log((1 + r) / (1 - r)) }
fisher_z_inverse <- function(z) { (exp(2 * z) - 1) / (exp(2 * z) + 1) }

# 转换 -> 平均 -> 逆转换
z_matrices <- lapply(cor_matrices_list, function(r_mat) fisher_z(r_mat))
z_pooled <- Reduce("+", z_matrices) / M
r_pooled <- fisher_z_inverse(z_pooled)

# 计算P值
n_avg <- Reduce("+", n_samples_list) / M
se_z <- 1 / sqrt(n_avg - 3)
z_stat <- z_pooled / se_z
p_values <- 2 * pnorm(-abs(z_stat))

cat("- 池化完成。平均样本量:", round(mean(n_avg), 1), "\n")

# ------------------------------------------------------------------------------
# 步骤5: FDR校正
# ------------------------------------------------------------------------------
cat("\n===== 步骤5: 多重检验校正 =====\n")

p_values_vector <- as.vector(p_values)
p_adj_fdr <- p.adjust(p_values_vector, method = "fdr")

p_adj_matrix <- matrix(p_adj_fdr, nrow = nrow(p_values), ncol = ncol(p_values),
                       dimnames = dimnames(p_values))

# 生成标记
sig_matrix <- matrix("", nrow = nrow(p_adj_matrix), ncol = ncol(p_adj_matrix),
                     dimnames = dimnames(p_adj_matrix))
sig_matrix[p_adj_matrix < 0.001] <- "***"
sig_matrix[p_adj_matrix >= 0.001 & p_adj_matrix < 0.01] <- "**"
sig_matrix[p_adj_matrix >= 0.01 & p_adj_matrix < 0.05] <- "*"
sig_matrix[p_adj_matrix >= 0.05 & p_adj_matrix < 0.10] <- "†"

cat("- 显著相关 (FDR < 0.05): ", sum(p_adj_matrix < 0.05), "对\n")

# ==============================================================================
# 步骤6: 生成摘要表 (双版本：ID版 + 名称版) [修正列名匹配错误]
# ==============================================================================
cat("\n===== 步骤6: 生成摘要表 =====\n")

# --- 6.1 [原始流程] 生成并保存 ID 版 ---
cor_summary <- data.frame(
  Metal = rep(rownames(r_pooled), times = ncol(r_pooled)),
  Metabolite = rep(colnames(r_pooled), each = nrow(r_pooled)),
  Spearman_r = as.vector(r_pooled),
  P_value = as.vector(p_values),
  FDR = as.vector(p_adj_matrix),
  N = as.vector(n_avg),
  Significance = as.vector(sig_matrix)
) %>%
  arrange(FDR, desc(abs(Spearman_r)))

# 保存 ID 版全表
write.csv(cor_summary, 
          "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_full.csv",
          row.names = FALSE)

# 保存 ID 版显著表
cor_sig <- cor_summary %>% filter(FDR < 0.05)
write.csv(cor_sig, 
          "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_significant.csv",
          row.names = FALSE)

cat("✅ [ID版] 结果表已保存 (outputs/..._full.csv & ..._significant.csv)\n")

# --- 6.2 [新增流程] 生成并保存 名称版 (指定 MNO 和 name) ---
if(file.exists("data/metabolism.xlsx")) {
  library(readxl)
  
  # 读取原始表
  metabo_raw <- read_excel("data/metabolism.xlsx", sheet = "original")
  
  # 【关键修正】：直接指定从 MNO 和 name 列提取
  # 优先检查列名是否存在，防止报错
  if("MNO" %in% names(metabo_raw) & "name" %in% names(metabo_raw)) {
    metabo_mapping <- metabo_raw %>%
      dplyr::select(NO = MNO, Metabolite_Name = name)
  } else {
    # 如果找不到 name，尝试找 NAME 或 Name
    name_col <- intersect(c("Name", "NAME", "name"), names(metabo_raw))[1]
    id_col <- intersect(c("MNO", "NO", "ID"), names(metabo_raw))[1]
    
    if(is.na(name_col) | is.na(id_col)) {
      warning("⚠️ 无法在 Excel 中找到指定的 'MNO' 或 'name' 列，名称版将使用 ID。")
      metabo_mapping <- data.frame(NO = unique(cor_summary$Metabolite), Metabolite_Name = unique(cor_summary$Metabolite))
    } else {
      metabo_mapping <- metabo_raw %>%
        dplyr::select(NO = all_of(id_col), Metabolite_Name = all_of(name_col))
    }
  }
} else {
  warning("⚠️ 未找到 metabolism.xlsx，名称版将继续使用 ID")
  metabo_mapping <- data.frame(NO = unique(cor_summary$Metabolite), Metabolite_Name = unique(cor_summary$Metabolite))
}

# 确保 ID 是字符型，防止匹配失败
metabo_mapping$NO <- as.character(metabo_mapping$NO)
metabo_mapping$Metabolite_Name <- as.character(metabo_mapping$Metabolite_Name)

# 合并名称 (Left Join)
cor_summary_named <- cor_summary %>%
  left_join(metabo_mapping, by = c("Metabolite" = "NO")) %>%
  # 【修正】：此时 Metabolite_Name 列肯定存在了（即使是NA）
  mutate(Metabolite_Name = ifelse(is.na(Metabolite_Name), Metabolite, Metabolite_Name)) %>%
  dplyr::select(Metal, Metabolite, Metabolite_Name, everything()) %>%
  arrange(FDR, desc(abs(Spearman_r)))

# 保存名称版全表
write.csv(cor_summary_named, 
          "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_full_with_names.csv",
          row.names = FALSE)

# 保存名称版显著表
cor_sig_named <- cor_summary_named %>% filter(FDR < 0.05)
write.csv(cor_sig_named, 
          "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_significant_with_names.csv",
          row.names = FALSE)

cat("✅ [名称版] 结果表已保存 (outputs/..._with_names.csv)\n")


# ==============================================================================
# 步骤7: 可视化 - 热图 (双版本)
# ==============================================================================
cat("\n===== 步骤7: 绘制热图 (双版本) =====\n")

# 7.1 准备绘图数据
has_sig_metal <- apply(p_adj_matrix < 0.05, 1, any)
has_sig_metabo <- apply(p_adj_matrix < 0.05, 2, any)

if(sum(has_sig_metal) > 0 & sum(has_sig_metabo) > 0) {
  metals_to_plot <- rownames(r_pooled)[has_sig_metal]
  metabo_to_plot <- colnames(r_pooled)[has_sig_metabo]
} else {
  metals_to_plot <- rownames(r_pooled)
  metabo_to_plot <- colnames(r_pooled)
}

r_for_plot <- r_pooled[metals_to_plot, metabo_to_plot, drop = FALSE]
sig_for_plot <- sig_matrix[metals_to_plot, metabo_to_plot, drop = FALSE]

# 动态计算尺寸
plot_width <- max(10, ncol(r_for_plot) * 0.4 + 2)
plot_height <- max(8, nrow(r_for_plot) * 0.3 + 2)

# --- 7.2 [原始流程] 绘制 ID 版热图 ---
png_id <- "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_heatmap.png"

tryCatch({
  png(png_id, width = plot_width, height = plot_height, units = "in", res = 300)
  pheatmap(
    r_for_plot,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    display_numbers = sig_for_plot,
    number_color = "black",
    fontsize_number = 8,
    cluster_rows = TRUE, cluster_cols = TRUE,
    main = "Metal-Metabolite Correlation (ID)",
    fontsize = 10
  )
  dev.off()
  cat("✅ [ID版] 热图已保存: ..._heatmap.png\n")
}, error = function(e) { try(dev.off(), silent=TRUE); cat("❌ [ID版] 绘图失败:", e$message, "\n") })

# --- 7.3 [新增流程] 绘制 名称版热图 ---
# 准备标签 (从 metabo_mapping 获取)
plot_col_labels <- metabo_mapping$Metabolite_Name[match(colnames(r_for_plot), metabo_mapping$NO)]
# 填补 NA
plot_col_labels[is.na(plot_col_labels)] <- colnames(r_for_plot)[is.na(plot_col_labels)]

png_named <- "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_heatmap_named.png"

tryCatch({
  png(png_named, width = plot_width + 1, height = plot_height + 1, units = "in", res = 300)
  pheatmap(
    r_for_plot,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks = seq(-1, 1, length.out = 101),
    display_numbers = sig_for_plot,
    number_color = "black",
    fontsize_number = 8,
    labels_col = plot_col_labels, # 使用真实名称
    cluster_rows = TRUE, cluster_cols = TRUE,
    main = "Metal-Metabolite Correlation (Named)",
    fontsize = 10,
    angle_col = 45
  )
  dev.off()
  cat("✅ [名称版] 热图已保存: ..._heatmap_named.png\n")
}, error = function(e) { try(dev.off(), silent=TRUE); cat("❌ [名称版] 绘图失败:", e$message, "\n") })

# --- 7.4 [新增流程] 绘制 分类注释热图 ---
annot_df <- NULL
hybrid_file <- "outputs/13 candidate_metabolites/candidate_metabolites_final_hybrid.csv"

if(file.exists(hybrid_file)) {
  temp_df <- read.csv(hybrid_file)
  if("Source" %in% names(temp_df)) {
    annot_df <- temp_df %>% dplyr::select(MNO, Category = Source)
  }
} 

if(!is.null(annot_df)) {
  annotation_col <- data.frame(row.names = colnames(r_for_plot))
  annotation_col$Category <- annot_df$Category[match(rownames(annotation_col), annot_df$MNO)]
  annotation_col$Category[is.na(annotation_col$Category)] <- "Unknown"
  
  png_annot <- "outputs/14 MMcorrelation_analysis/metal_metabolite_correlation_heatmap_annotated_named.png"
  
  tryCatch({
    png(png_annot, width = plot_width + 2, height = plot_height + 1, units = "in", res = 300)
    pheatmap(
      r_for_plot,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-1, 1, length.out = 101),
      display_numbers = sig_for_plot,
      number_color = "black",
      fontsize_number = 8,
      annotation_col = annotation_col, 
      labels_col = plot_col_labels,   
      cluster_rows = TRUE, cluster_cols = TRUE,
      main = "Metal-Metabolite Correlation (Annotated)",
      fontsize = 10,
      angle_col = 45
    )
    dev.off()
    cat("✅ [分类版] 热图已保存: ..._heatmap_annotated_named.png\n")
  }, error = function(e) { try(dev.off(), silent=TRUE); cat("❌ [分类版] 绘图失败:", e$message, "\n") })
}

# ==============================================================================
# 步骤8: 可视化 - 网络图 (双版本)
# ==============================================================================
cat("\n===== 步骤8: 绘制网络图 (双版本) =====\n")

if(exists("cor_sig") && nrow(cor_sig) > 0) {
  library(igraph)
  library(ggraph)
  
  # --- 8.1 数据构建 ---
  edges <- cor_sig %>%
    dplyr::select(from = Metal, to = Metabolite, r_val = Spearman_r, FDR) %>%
    mutate(
      edge_type = ifelse(r_val > 0, "Positive", "Negative"),
      layout_weight = abs(r_val),
      plot_width = abs(r_val),
      edge_alpha = -log10(FDR)
    )
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  V(g)$type <- ifelse(V(g)$name %in% metal_vars, "Metal", "Metabolite")
  V(g)$color <- ifelse(V(g)$type == "Metal", "#4DBBD5FF", "#E64B35FF")
  V(g)$degree <- degree(g)
  V(g)$size <- log1p(V(g)$degree) * 2 + 3 
  E(g)$layout_weight <- pmax(E(g)$layout_weight, 0.001) # 防止0权重
  
  # --- 8.2 [ID 版] ---
  set.seed(123)
  p_net_id <- ggraph(g, layout = "fr", weights = layout_weight) +  
    geom_edge_link(aes(color = edge_type, width = plot_width, alpha = edge_alpha), show.legend = TRUE) +
    geom_node_point(aes(color = type, size = size), show.legend = c(color=TRUE, size=FALSE)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "bold") +
    scale_edge_color_manual(values = c("Positive" = "#D6604D", "Negative" = "#4393C3")) +
    scale_edge_width(range = c(0.5, 2.0), name = "|r|") +
    scale_edge_alpha(range = c(0.3, 0.8), guide = "none") +
    scale_color_manual(values = c("Metal" = "#4DBBD5FF", "Metabolite" = "#E64B35FF")) +
    labs(title = "Network (ID Version)", subtitle = paste("FDR < 0.05, N =", nrow(cor_sig))) +
    theme_void()
  
  ggsave("outputs/14 MMcorrelation_analysis/metal_metabolite_network.png", 
         plot = p_net_id, width = 12, height = 10, dpi = 300)
  cat("✅ [ID版] 网络图已保存: ..._network.png\n")
  
  # --- 8.3 [名称版] ---
  g_named <- g
  # 构建 ID -> Name 字典
  id_to_name <- setNames(metabo_mapping$Metabolite_Name, metabo_mapping$NO)
  
  # 替换 Label
  current_labels <- V(g_named)$name
  new_labels <- sapply(current_labels, function(x) {
    if(x %in% names(id_to_name)) return(id_to_name[x]) else return(x)
  })
  V(g_named)$label <- new_labels
  
  set.seed(123)
  p_net_named <- ggraph(g_named, layout = "fr", weights = layout_weight) +  
    geom_edge_link(aes(color = edge_type, width = plot_width, alpha = edge_alpha), show.legend = TRUE) +
    geom_node_point(aes(color = type, size = size), show.legend = c(color=TRUE, size=FALSE)) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3.0, fontface = "bold") + 
    scale_edge_color_manual(values = c("Positive" = "#D6604D", "Negative" = "#4393C3")) +
    scale_edge_width(range = c(0.5, 2.0), name = "|r|") +
    scale_edge_alpha(range = c(0.3, 0.8), guide = "none") +
    scale_color_manual(values = c("Metal" = "#4DBBD5FF", "Metabolite" = "#E64B35FF")) +
    labs(title = "Network (Named Version)", subtitle = paste("FDR < 0.05, N =", nrow(cor_sig))) +
    theme_void()
  
  ggsave("outputs/14 MMcorrelation_analysis/metal_metabolite_network_named.png", 
         plot = p_net_named, width = 14, height = 12, dpi = 300)
  cat("✅ [名称版] 网络图已保存: ..._network_named.png\n")
  
} else {
  cat("⚠️ 无显著相关，跳过网络图绘制\n")
}

# ==============================================================================
# 步骤9: 生成分析报告
# ==============================================================================
cat("\n===== 步骤9: 生成摘要报告 =====\n")

report_file <- "outputs/14 MMcorrelation_analysis/correlation_analysis_report.txt"
sink(report_file)

cat("=============================================================\n")
cat("           金属-代谢物相关分析摘要报告\n")
cat("=============================================================\n\n")
cat("分析日期:", format(Sys.Date(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("显著相关 (FDR < 0.05): ", sum(p_adj_matrix < 0.05), "对\n\n")

if(exists("cor_sig_named") && nrow(cor_sig_named) > 0) {
  cat("Top 10 最强相关 (详细名单见CSV):\n")
  # 打印时指定列名
  print(head(cor_sig_named %>% dplyr::select(Metal, Metabolite_Name, Spearman_r, FDR), 10), row.names=FALSE)
} else {
  cat("无显著相关。\n")
}

cat("\n输出文件说明:\n")
cat("1. ID版: ..._full.csv, ..._heatmap.png, ..._network.png\n")
cat("2. 名称版: ..._with_names.csv, ..._heatmap_named.png, ..._network_named.png\n")
sink()
cat("✅ 报告已生成: correlation_analysis_report.txt\n")

