# ==============================================================================
# 分析1.2: 代谢物组间差异分析 (LVR vs 非LVR) 
# 方法: IPW加权logistic回归 + 多重插补池化 (Doubly Robust Estimation)
# 修正人: Senior Statistical Consultant
# ==============================================================================

rm(list = ls())
gc()

# 加载必要的包
# 检查并自动加载包的辅助函数
ensure_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(paste("需要安装包:", pkg))
  library(pkg, character.only = TRUE)
}

ensure_library("dplyr")
ensure_library("mice")
ensure_library("broom")
ensure_library("ggplot2")
ensure_library("ggrepel")
ensure_library("survey")
ensure_library("readxl")
ensure_library("here") # 推荐使用here处理路径

# 创建输出目录
output_dir <- "outputs/15 metabolite_group_differences" # 修改编号以保持顺序清晰
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  代谢物组间差异分析 (LVR vs 非LVR) [优化版]    ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")


# ==============================================================================
# 步骤1: 数据准备 (修正版 V4 - 智能识别文件类型，彻底解决兼容性)
# ==============================================================================
cat("===== 步骤1: 数据准备 (智能兼容模式) =====\n")

library(dplyr)
library(mice)
library(readxl)

# 1.1 定义文件路径
metabo_file_path <- "outputs/mice_imputation_with_metabo.rds"
weight_file_path <- "outputs/4 IPW2_complete/imputed_data_with_ipw_weights_extended.rds"

# 1.2 加载数据
if(!file.exists(metabo_file_path)) stop("❌ 找不到代谢物文件")
if(!file.exists(weight_file_path)) stop("❌ 找不到权重文件")

raw_metabo <- readRDS(metabo_file_path)
raw_weights <- readRDS(weight_file_path)

cat("- 文件加载成功\n")
cat(sprintf("  - 代谢物文件类型: %s\n", paste(class(raw_metabo), collapse=", ")))
cat(sprintf("  - 权重文件类型:   %s\n", paste(class(raw_weights), collapse=", ")))

# --- 核心函数：智能转 Long 格式 ---
# 无论输入是 mids 对象还是 data.frame，都能统统搞定
smart_to_long <- function(obj, obj_name) {
  cat(sprintf("\n  > 正在处理 %s ...\n", obj_name))
  
  # 情况 A: 输入是 data.frame (已经是表格了)
  if (inherits(obj, "data.frame")) {
    cat("    识别为 data.frame，无需解包，直接整理。\n")
    df <- obj
    
    # 检查是否有 .imp 列 (如果没有，可能是一个单独的数据集)
    if (!".imp" %in% names(df)) {
      warning(paste("⚠️", obj_name, "缺少 .imp 列，假设这是原始数据 (.imp=0)"))
      df$.imp <- 0
    }
    # 检查是否有 .id 列
    if (!".id" %in% names(df)) {
      df$.id <- 1:nrow(df)
    }
    return(df)
  }
  
  # 情况 B: 输入是 mids 对象 (mice 插补对象)
  if (inherits(obj, "mids")) {
    cat("    识别为 mids 对象，开始手动硬解包...\n")
    
    # 获取参数
    m <- obj$m
    if (is.null(m)) m <- 0 # 防止 m 为空
    
    # 如果 m=0 或数据为空，直接返回原始数据
    if (m == 0) {
      df <- obj$data
      df$.imp <- 0
      df$.id <- 1:nrow(df)
      return(df)
    }
    
    # 硬解包逻辑 (之前的 V3 逻辑，但加了 m 检查)
    out_list <- vector("list", m + 1)
    
    # 原始数据
    temp_0 <- obj$data
    temp_0$.imp <- 0
    temp_0$.id <- 1:nrow(temp_0)
    out_list[[1]] <- temp_0
    
    # 循环提取
    vars_with_na <- names(obj$imp)
    # 进度条
    for(i in 1:m) {
      temp_i <- obj$data
      for(var in vars_with_na) {
        imp_mat <- obj$imp[[var]]
        if(!is.null(imp_mat) && ncol(imp_mat) >= i) {
          na_idx <- which(is.na(obj$data[[var]]))
          if(length(na_idx) > 0) {
            temp_i[na_idx, var] <- imp_mat[, i]
          }
        }
      }
      temp_i$.imp <- i
      temp_i$.id <- 1:nrow(temp_i)
      out_list[[i + 1]] <- temp_i
    }
    
    return(bind_rows(out_list))
  }
  
  stop(paste("❌ 未知的数据类型:", paste(class(obj), collapse=", ")))
}

# 1.3 执行智能转换
long_metabo <- smart_to_long(raw_metabo, "代谢物数据")
long_weights <- smart_to_long(raw_weights, "权重数据")

# 1.4 统一 ID 列名
standardize_id <- function(df, df_name) {
  # 移除可能的 rowname 列干扰
  if("rownames" %in% names(df)) df <- df %>% select(-rownames)
  
  possible_ids <- c("ID", "NO", "MNO", "patient_id", "PatientID")
  found_id <- intersect(possible_ids, names(df))[1]
  
  if(is.na(found_id)) {
    # 尝试找第一列如果是整数
    first_col <- names(df)[1]
    warning(paste("⚠️", df_name, "未找到标准ID，尝试使用第一列:", first_col))
    found_id <- first_col
  }
  
  if(found_id != "ID") {
    cat(sprintf("    - 将 '%s' 重命名为 'ID'\n", found_id))
    df <- df %>% rename(ID = all_of(found_id))
  }
  # 强制 ID 为字符型以确保匹配
  df$ID <- as.character(df$ID)
  return(df)
}

long_metabo <- standardize_id(long_metabo, "代谢物数据")
long_weights <- standardize_id(long_weights, "权重数据")

# 1.5 合并权重
cat("\n- 正在合并...\n")
possible_weights <- c("sw_trunc", "ipw", "sw", "weight", "weights", "iptw")
weight_var_name <- intersect(possible_weights, names(long_weights))[1]

if(is.na(weight_var_name)) {
  # 打印一下列名帮助调试
  cat("  可用列名:", paste(names(long_weights), collapse=", "), "\n")
  stop("❌ 无法识别权重变量名！")
} else {
  cat(sprintf("  ✓ 权重变量: '%s'\n", weight_var_name))
}

# 准备合并子集
weights_subset <- long_weights %>%
  select(.imp, .id, ID, all_of(weight_var_name))

# 执行合并
long_merged <- long_metabo %>%
  left_join(weights_subset, by = c(".imp", "ID")) 
# 注意：这里去掉了 .id 作为合并键，因为如果两个文件来源不同，行号(.id)可能对不上，
# 但 .imp (插补集序号) 和 ID (病人ID) 是一定能对上的。

# 检查合并结果
if(nrow(long_merged) != nrow(long_metabo)) {
  warning("⚠️ 合并后行数发生变化，请检查 ID 是否唯一。")
}

# 1.6 重建 Mids 对象
mids_obj <- as.mids(long_merged)
cat("  ✓ 数据准备完成 (V4)！\n\n")

# --- 加载辅助数据 ---
candidate_file <- "outputs/13 candidate_metabolites/final_metabolite_list.rds"
if(file.exists(candidate_file)) {
  metabo_vars <- readRDS(candidate_file)
} else {
  metabo_vars <- names(mids_obj$data)[grep("^M[0-9]+$", names(mids_obj$data))]
}
metabo_vars <- metabo_vars[metabo_vars %in% names(mids_obj$data)]

metabo_mapping <- read_excel("data/metabolism.xlsx", sheet = "original")
if("MNO" %in% names(metabo_mapping)) {
  metabo_mapping <- metabo_mapping %>% select(NO = MNO, Metabolite_Name = name)
} else {
  metabo_mapping <- metabo_mapping %>% select(NO, Metabolite_Name = name)
}
cat("- 辅助文件加载完成。\n")


# ==============================================================================
# 步骤2: 定义结局变量 (LVR) - 确保在插补对象中正确生成
# ==============================================================================
cat("===== 步骤2: 定义结局变量 (LVR) =====\n")

data_check <- complete(mids_obj, 1)

# 检查是否需要重新计算LVR
need_calc_lvr <- !"LVR" %in% names(data_check)

if(need_calc_lvr) {
  cat("- 未找到LVR变量，开始基于插补值计算...\n")
  
  if(!all(c("LVEDV_fu", "LVEDV_baseline") %in% names(data_check))) {
    stop("❌ 缺少 LVEDV_fu 或 LVEDV_baseline 变量，无法计算 LVR。")
  }
  
  # [顾问修正]: 使用 complete(..., include=TRUE) 确保包含原始数据，这对于 as.mids 至关重要
  long_data <- complete(mids_obj, action = "long", include = TRUE)
  
  long_data <- long_data %>%
    mutate(
      delta_LVEDV = (LVEDV_fu - LVEDV_baseline) / LVEDV_baseline * 100,
      # 定义 LVR: delta >= 20% 为 1, 否则为 0
      LVR = ifelse(!is.na(delta_LVEDV) & delta_LVEDV >= 20, 1, 0)
    )
  
  # 重建 mids 对象
  mids_obj <- as.mids(long_data)
  data_check <- complete(mids_obj, 1) # 更新检查数据
  cat("  ✓ LVR变量已创建并回写至mids对象\n")
}

# 检查LVR分布
lvr_counts <- table(data_check$LVR, useNA = "ifany")
lvr_rate <- prop.table(table(data_check$LVR))[2] * 100

cat("- LVR分布 (基于第1个插补集):\n")
print(lvr_counts)
cat("  发生率:", round(lvr_rate, 2), "%\n\n")

if(is.na(lvr_rate) || lvr_rate < 1) warning("⚠️ LVR事件极少，模型可能不收敛！")

# ==============================================================================
# 步骤3: 定义协变量 (双重稳健估计配置)
# ==============================================================================
cat("===== 步骤3: 定义协变量 =====\n")

# [顾问提示]: 在IPW分析中再次加入协变量是为了实现“双重稳健”。
# 即使权重模型有误，结果模型正确，或者反之，估计量仍是一致的。

# 原始定义列表
covariates_ideal <- c(
  "age", "gender", "resident", "DM", "hypertension", "pPCI", "STEMI",
  "EF_baseline", "LVEDV_baseline", "GRACE_in", "IN_killip",
  "cTnIpeak", "NTproBNP_peak", "CKMB", "WBC", "HGB", "PLT",
  "CRP", "CHOL", "LDL", "AST", "Statin", "Lesion_no"
)

# 变量名兼容性处理 (Map old names to new names)
name_map <- c(
  "sex" = "gender", "LVEF_baseline" = "EF_baseline", 
  "peak_CK_MB" = "CKMB", "diabetes" = "DM"
)

# 检查并重命名mids对象中的变量
data_names <- names(mids_obj$data)
for(old in names(name_map)) {
  new <- name_map[old]
  if(old %in% data_names && !(new %in% data_names)) {
    # 如果旧名存在且新名不存在，需要重命名
    # 注意：直接修改mids$data的列名比较危险，建议重新构建
    # 这里采用简化方案，仅修改data副本，实际分析时需注意
    cat(sprintf("  - 检测到旧变量名 '%s'，将在分析中使用\n", old))
    # 更新covariates_ideal中的名称匹配实际数据
    covariates_ideal[covariates_ideal == new] <- old
  }
}

# 最终筛选可用协变量
covariates <- intersect(covariates_ideal, names(data_check))
missing_covs <- setdiff(covariates_ideal, names(data_check))

cat("- 纳入模型的协变量:", length(covariates), "个\n")
if(length(missing_covs) > 0) cat("⚠️ 缺失协变量:", paste(missing_covs, collapse=", "), "\n")

# ==============================================================================
# 步骤4: 识别IPW权重变量
# ==============================================================================
cat("===== 步骤4: 识别IPW权重 =====\n")

weight_vars <- c("sw_trunc", "ipw", "sw", "weights", "iptw")
weight_var <- intersect(weight_vars, names(data_check))[1]

use_weights <- FALSE
if(!is.na(weight_var)) {
  use_weights <- TRUE
  cat(sprintf("✅ 锁定权重变量: '%s'\n", weight_var))
  cat("   分析模式: IPW加权 Logistic 回归 (Doubly Robust)\n")
} else {
  cat("🛑 警告: 未在数据中找到权重变量！\n")
  cat("   分析模式: 普通 Logistic 回归 (Unweighted)\n")
  cat("   [顾问建议]: 如果本意是做IPW，请检查输入数据是否包含权重列。\n")
}
cat("\n")

# ==============================================================================
# 步骤5: 拟合模型 (核心修正部分)
# ==============================================================================
cat("===== 步骤5: 拟合模型与池化 (Rubin's Rules) =====\n")

# 定义核心拟合函数
fit_pooled_model <- function(mids_data, metab, covars, outcome="LVR", 
                             w_var=NULL, weighted=TRUE) {
  
  # 构建公式
  # 注意：代谢物通常需要标准化(scale)以便比较OR，或者保持原样看单位变化
  # 这里保持原样
  form <- as.formula(paste(outcome, "~", metab, "+", paste(covars, collapse = " + ")))
  
  if(weighted && !is.null(w_var)) {
    # --- 加权分析路径 (Manual Pooling) ---
    M <- mids_data$m
    coefs <- numeric(M)
    vars  <- numeric(M) # 存储方差 (SE^2)
    
    # 遍历每个插补数据集
    for(i in 1:M) {
      # 提取单个完整数据集
      dat_i <- complete(mids_data, i)
      
      # [重要] 移除权重缺失或为0的行，防止svydesign报错
      dat_i <- dat_i[!is.na(dat_i[[w_var]]) & dat_i[[w_var]] > 0, ]
      
      # 定义设计对象
      des <- svydesign(ids = ~1, weights = as.formula(paste("~", w_var)), data = dat_i)
      
      # 拟合模型
      # quasibinomial 用于避免非整数权重的警告，但在二分类下系数与binomial一致
      mod <- svyglm(form, design = des, family = quasibinomial())
      
      # 提取代谢物的系数和标准误
      summ <- summary(mod)$coefficients
      if(metab %in% rownames(summ)) {
        coefs[i] <- summ[metab, "Estimate"]
        vars[i]  <- summ[metab, "Std. Error"]^2
      } else {
        # 如果模型因共线性剔除了变量
        coefs[i] <- NA
        vars[i] <- NA
      }
    }
    
    # 如果有任何插补集失败，返回NA
    if(any(is.na(coefs))) return(NULL)
    
    # Rubin's Rules 池化
    pool_q <- mean(coefs)                 # Combined estimate
    pool_u <- mean(vars)                  # Within variance
    pool_b <- var(coefs)                  # Between variance
    pool_t <- pool_u + (1 + 1/M) * pool_b # Total variance
    
    pool_se <- sqrt(pool_t)
    
    # 计算自由度 (Barnard-Rubin adjustment for small samples is better, but traditional is ok)
    r <- (1 + 1/M) * pool_b / pool_u
    v_old <- (M - 1) * (1 + 1/r)^2
    # 我们可以直接使用 v_old 作为 df
    
    p_val <- 2 * (1 - pt(abs(pool_q / pool_se), df = v_old))
    
    res <- data.frame(
      Metabolite = metab,
      Beta = pool_q,
      SE = pool_se,
      OR = exp(pool_q),
      OR_lower = exp(pool_q - 1.96 * pool_se),
      OR_upper = exp(pool_q + 1.96 * pool_se),
      P_value = p_val
    )
    
  } else {
    # --- 未加权分析路径 (Standard MICE) ---
    # 使用 with() 和 pool()
    mod_list <- with(mids_data, glm(as.formula(paste(outcome, "~", metab, "+", paste(covars, collapse = " + "))), 
                                    family = binomial()))
    pooled <- pool(mod_list)
    res_summ <- summary(pooled, conf.int = TRUE)
    
    # 提取目标行
    target <- res_summ[res_summ$term == metab, ]
    
    if(nrow(target) == 0) return(NULL)
    
    res <- data.frame(
      Metabolite = metab,
      Beta = target$estimate,
      SE = target$std.error,
      OR = exp(target$estimate),
      OR_lower = exp(target$`2.5 %`), # 注意broom版本差异，可能是 2.5 %
      OR_upper = exp(target$`97.5 %`),
      P_value = target$p.value
    )
  }
  return(res)
}

# 批量运行
cat("- 开始拟合", length(metabo_vars), "个模型...\n")
results_list <- list()
pb <- txtProgressBar(min = 0, max = length(metabo_vars), style = 3)

for(i in seq_along(metabo_vars)) {
  metab <- metabo_vars[i]
  
  # 使用 tryCatch 防止单个模型报错中断循环
  out <- tryCatch({
    fit_pooled_model(mids_obj, metab, covariates, 
                     w_var = weight_var, weighted = use_weights)
  }, error = function(e) {
    return(NULL) # 忽略错误
  })
  
  if(!is.null(out)) results_list[[i]] <- out
  setTxtProgressBar(pb, i)
}
close(pb)

final_results <- bind_rows(results_list)

cat("\n✅ 模型拟合完成。成功:", nrow(final_results), "/", length(metabo_vars), "\n")

# ==============================================================================
# 步骤6: 多重检验校正与注释
# ==============================================================================
cat("===== 步骤6: 结果整理与校正 =====\n")

final_results <- final_results %>%
  mutate(
    FDR = p.adjust(P_value, method = "fdr"),
    Bonferroni = p.adjust(P_value, method = "bonferroni"),
    Log2FC = log2(OR), # Log2 Odds Ratio 作为替代 FC
    Significance = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      P_value < 0.05 ~ "P < 0.05",
      TRUE ~ "NS"
    )
  ) %>%
  # 关联名称
  left_join(metabo_mapping, by = c("Metabolite" = "NO")) %>%
  mutate(
    # 如果没有映射名，使用ID
    Metabolite_Name = ifelse(is.na(Metabolite_Name), Metabolite, Metabolite_Name),
    Label = ifelse(P_value < 0.01 & abs(Log2FC) > 0.5, Metabolite_Name, NA) # 仅标注显著点
  ) %>%
  arrange(P_value)

# 打印Top结果
cat("Top 5 显著代谢物 (按P值):\n")
print(head(final_results %>% select(Metabolite_Name, OR, P_value, FDR), 5))

# ==============================================================================
# 步骤7: 可视化 (火山图)
# ==============================================================================
cat("===== 步骤7: 绘图 =====\n")

p_vol <- ggplot(final_results, aes(x = Log2FC, y = -log10(P_value))) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("FDR < 0.05" = "red", "P < 0.05" = "orange", "NS" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey") +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(
    title = paste0("Volcano Plot: Metabolites associated with LVR"),
    subtitle = ifelse(use_weights, "Method: IPW-weighted Logistic Regression (Doubly Robust)", "Method: Unweighted Logistic Regression"),
    x = "Log2 Odds Ratio (Log2FC)",
    y = "-Log10 P-value"
  )

ggsave(file.path(output_dir, "volcano_plot.png"), p_vol, width = 8, height = 6)

# ==============================================================================
# 步骤8: 保存文件
# ==============================================================================
write.csv(final_results, file.path(output_dir, "metabolite_LVR_results.csv"), row.names = FALSE)
cat("✅ 结果已保存至:", output_dir, "\n")
