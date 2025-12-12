# =============================================================================
# 脚本名称: 01_Generate_Table1.R
# 功能: 生成Table 1 - IPW核心变量基线特征
# 输入: outputs/imputed_data_with_ipw_weights.rds, variable_config.RData
# 输出: outputs/Table1_Baseline_Core. docx
# =============================================================================
rm(list = ls())
gc()
invisible(gc())
# 1. 环境准备
# -----------------------------------------------------------------------------
if(! require(pacman)) install.packages("pacman")
pacman::p_load(flextable, officer, dplyr, tibble)

if(!dir.exists("outputs")) dir.create("outputs")

# 2. 加载数据和配置
# -----------------------------------------------------------------------------
data_all <- readRDS("outputs/imputed_data_with_ipw_weights.rds")
load("variable_config.RData")

# 3. 提取第1个插补数据集
# -----------------------------------------------------------------------------
data_complete <- data_all %>% 
  filter(.imp == 1) %>% 
  select(-.imp, -.id)

cat("✅ 数据提取完成:  N =", nrow(data_complete), "\n")

# 4. 【修复】智能识别并划分队列
# -----------------------------------------------------------------------------
cat("\n📊 识别 has_fu_echo 编码.. .\n")

# 转换为字符型进行统一处理
has_fu_char <- as.character(data_complete$has_fu_echo)
unique_vals <- unique(has_fu_char[! is.na(has_fu_char)])

cat("  唯一值:", paste(unique_vals, collapse = ", "), "\n")

# 频数表
freq_table <- table(has_fu_char, useNA = "always")
print(freq_table)

# 智能匹配"有随访"的标识
positive_patterns <- c("Yes", "yes", "Y", "1", "有", "是", "TRUE", "True", "true")
has_fu_indicator <- NULL

for (pattern in positive_patterns) {
  if (pattern %in% unique_vals) {
    has_fu_indicator <- pattern
    break
  }
}

# 如果没找到,使用频数较小的那个(假设有随访的人较少)
if (is.null(has_fu_indicator)) {
  freq_no_na <- freq_table[names(freq_table) != "<NA>"]
  has_fu_indicator <- names(freq_no_na)[which.min(freq_no_na)]
  warning("无法自动识别,使用频数较小的值: ", has_fu_indicator)
}

no_fu_indicator <- setdiff(unique_vals, has_fu_indicator)[1]

cat(sprintf("✅ 识别结果:  '%s' = 有随访, '%s' = 无随访\n\n", 
            has_fu_indicator, no_fu_indicator))

# 创建标准化的分组变量
data_complete <- data_complete %>% 
  mutate(
    fu_status = case_when(
      as.character(has_fu_echo) == has_fu_indicator ~ "有随访",
      as.character(has_fu_echo) == no_fu_indicator ~ "无随访",
      TRUE ~ NA_character_
    )
  )

# 划分队列
cohort_total <- data_complete
cohort_analysis <- data_complete %>% filter(fu_status == "有随访")
cohort_missing <- data_complete %>% filter(fu_status == "无随访")

# 输出结果
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("📋 队列划分结果:\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat(sprintf("总人群:    N = %d (100.0%%)\n", nrow(cohort_total)))
cat(sprintf("分析人群:   N = %d (%.1f%%)\n", 
            nrow(cohort_analysis), 
            100 * nrow(cohort_analysis) / nrow(cohort_total)))
cat(sprintf("剔除人群:   N = %d (%.1f%%)\n", 
            nrow(cohort_missing),
            100 * nrow(cohort_missing) / nrow(cohort_total)))
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# 安全检查
stopifnot("分析人群为空!" = nrow(cohort_analysis) > 0)
stopifnot("剔除人群为空!" = nrow(cohort_missing) > 0)

# 5. 定义Table 1核心变量
# -----------------------------------------------------------------------------
ipw_covariates <- c(
  "age", "gender", "resident", "DM", "hypertension",
  "pPCI", "STEMI", "EF_baseline", "LVEDV_baseline",
  "GRACE_in", "IN_killip", "cTnIpeak", "NTproBNP_peak",
  "CKMB", "WBC", "HGB", "PLT", "CRP"
)

exposure_vars <- c("Cu", "Zn", "Fe", "Se", "Pb")
outcome_vars <- c("LVEDV_fu", "LVESV_fu", "EF_fu")

table1_vars <- c(ipw_covariates, exposure_vars, outcome_vars)
table1_vars <- intersect(table1_vars, names(data_complete))

cat(sprintf("📊 Table 1 包含变量: %d 个\n\n", length(table1_vars)))

# 6. 识别分类变量
# -----------------------------------------------------------------------------
categorical_vars_table1 <- intersect(categorical_vars, table1_vars)

# 7. 统一因子水平
# -----------------------------------------------------------------------------
for (var in categorical_vars_table1) {
  if (var %in% names(data_complete)) {
    all_levels <- unique(c(
      as.character(cohort_total[[var]]),
      as.character(cohort_analysis[[var]]),
      as.character(cohort_missing[[var]])
    ))
    all_levels <- all_levels[!is.na(all_levels)]
    all_levels <- sort(all_levels)
    
    cohort_total[[var]] <- factor(cohort_total[[var]], levels = all_levels)
    cohort_analysis[[var]] <- factor(cohort_analysis[[var]], levels = all_levels)
    cohort_missing[[var]] <- factor(cohort_missing[[var]], levels = all_levels)
  }
}

# 8. 计算统计量
# -----------------------------------------------------------------------------
calc_stats <- function(data, vars, cat_vars, all_levels_list) {
  results <- tibble(variable = character(), label = character(), 
                    statistics = character(), row_type = character())
  
  for (var in vars) {
    if (! var %in% names(data)) next
    
    var_label <- label_mapping %>% 
      filter(variable == var) %>% 
      pull(new_label)
    if (length(var_label) == 0) var_label <- var
    
    if (var %in% cat_vars) {
      levels_to_use <- all_levels_list[[var]]
      if (is.null(levels_to_use)) levels_to_use <- levels(data[[var]])
      
      tbl <- table(factor(data[[var]], levels = levels_to_use), useNA = "no")
      total <- sum(tbl)
      
      results <- add_row(results, variable = var, label = var_label, 
                         statistics = "", row_type = "label")
      
      for (level_name in levels_to_use) {
        n <- tbl[level_name]
        if (is.na(n)) n <- 0
        pct <- ifelse(total > 0, 100 * n / total, 0)
        results <- add_row(results, variable = var, 
                           label = paste0("  ", level_name),
                           statistics = sprintf("%d (%.1f)", n, pct),
                           row_type = "level")
      }
    } else {
      vals <- data[[var]][!is.na(data[[var]])]
      if (length(vals) > 0) {
        stat_str <- sprintf("%. 1f ± %.1f", mean(vals), sd(vals))
      } else {
        stat_str <- "—"
      }
      results <- add_row(results, variable = var, label = var_label,
                         statistics = stat_str, row_type = "label")
    }
  }
  return(results)
}

all_levels_list <- list()
for (var in categorical_vars_table1) {
  if (var %in% names(data_complete)) {
    all_levels_list[[var]] <- levels(cohort_total[[var]])
  }
}

stats_total <- calc_stats(cohort_total, table1_vars, 
                          categorical_vars_table1, all_levels_list)
stats_analysis <- calc_stats(cohort_analysis, table1_vars, 
                             categorical_vars_table1, all_levels_list)
stats_missing <- calc_stats(cohort_missing, table1_vars, 
                            categorical_vars_table1, all_levels_list)

# 9. 合并统计量
# -----------------------------------------------------------------------------
table_data <- stats_total %>% 
  select(variable, label, row_type) %>% 
  mutate(
    col1 = stats_total$statistics,
    col2 = stats_analysis$statistics,
    col3 = stats_missing$statistics
  )

colnames(table_data) <- c(
  "variable", "特征", "row_type",
  sprintf("总人群\n(N=%d)", nrow(cohort_total)),
  sprintf("分析人群\n(N=%d)", nrow(cohort_analysis)),
  sprintf("剔除人群\n(N=%d)", nrow(cohort_missing))
)

# 10. 插入分组标题
# -----------------------------------------------------------------------------
core_groups <- list(
  list(
    title = "人口学特征与危险因素",
    vars = c("age", "gender", "resident", "DM", "hypertension")
  ),
  list(
    title = "急性心梗特征与治疗",
    vars = c("pPCI", "STEMI", "GRACE_in", "IN_killip")
  ),
  list(
    title = "基线实验室检查",
    vars = c("cTnIpeak", "NTproBNP_peak", "CKMB", "WBC", "HGB", "PLT", "CRP")
  ),
  list(
    title = "基线与随访超声心动图",
    vars = c("EF_baseline", "LVEDV_baseline", "EF_fu", "LVEDV_fu", "LVESV_fu")
  ),
  list(
    title = "血清金属元素浓度(μg/L)",
    vars = exposure_vars
  )
)

insert_positions <- list()
for (i in seq_along(core_groups)) {
  group <- core_groups[[i]]
  group_vars <- intersect(group$vars, table1_vars)
  if (length(group_vars) == 0) next
  
  first_row <- which(table_data$variable == group_vars[1] & 
                       table_data$row_type == "label")[1]
  if (! is.na(first_row)) {
    insert_positions[[length(insert_positions) + 1]] <- 
      list(position = first_row, title = group$title)
  }
}

insert_positions <- insert_positions[order(
  sapply(insert_positions, function(x) x$position), decreasing = TRUE)]

for (pos_info in insert_positions) {
  header_row <- table_data[pos_info$position, ]
  header_row$特征 <- pos_info$title
  header_row$row_type <- "group_header"
  header_row[, 4: 6] <- ""
  
  table_data <- bind_rows(
    table_data[1:(pos_info$position - 1), ],
    header_row,
    table_data[pos_info$position: nrow(table_data), ]
  )
}

table_data_final <- select(table_data, -variable, -row_type)

# 11. 生成Word三线表
# -----------------------------------------------------------------------------
ft <- flextable(table_data_final)
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- font(ft, fontname = "SimSun", part = "body", j = 1)
ft <- font(ft, fontname = "SimSun", part = "header")
ft <- fontsize(ft, size = 10, part = "all")
ft <- align(ft, align = "left", part = "body", j = 1)
ft <- align(ft, align = "center", part = "body", j = 2: 4)
ft <- align(ft, align = "center", part = "header")
ft <- width(ft, j = 1, width = 2.8)
ft <- width(ft, j = 2:4, width = 1.4)

group_rows <- which(table_data$row_type == "group_header")
if (length(group_rows) > 0) ft <- bold(ft, i = group_rows, j = 1)

ft <- border_remove(ft)
ft <- hline_top(ft, border = fp_border(color = "black", width = 1.5), part = "all")
ft <- hline(ft, border = fp_border(color = "black", width = 1.5), part = "header")
ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1.5), part = "body")
if (length(group_rows) > 0) {
  ft <- hline(ft, i = group_rows, 
              border = fp_border(color = "black", width = 0.5), part = "body")
}

ft <- add_footer_lines(ft, 
                       paste0("注: 连续变量表示为均值±标准差;分类变量表示为频数(百分比)。",
                              "分析人群指有完整随访超声数据的患者;剔除人群指缺失随访数据的患者。",
                              "表中变量为逆概率加权(IPW)模型中的协变量、主要暴露变量和结局变量。"))
ft <- font(ft, fontname = "SimSun", part = "footer")
ft <- fontsize(ft, size = 9, part = "footer")
ft <- align(ft, align = "left", part = "footer")

doc <- read_docx()
doc <- body_add_par(doc, "Table 1. 研究人群基线特征(IPW核心变量)", 
                    style = "heading 1")
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "outputs/Table1_Baseline_Core.docx")

cat("\n✅ Table 1 已保存至 outputs/Table1_Baseline_Core.docx\n")
