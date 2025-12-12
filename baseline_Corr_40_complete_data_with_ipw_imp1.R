############################################################
## Table 1 生成脚本 - 核心变量版 (40个变量)
## 功能: 三组基线特征比较 (总人群/分析人群/剔除人群)
## 输入: complete_data_with_ipw_imp1.rds
## 输出: Table1_Baseline_Core. docx
############################################################

## ═══════════════════════════════════════════════════════
## 环境设置
## ═══════════════════════════════════════════════════════
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(tableone)
library(flextable)
library(officer)
library(here)
library(pacman)
library(tibble)

cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║     Table 1 生成 - 核心变量版                       ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 1: 加载数据和变量标签
## ═══════════════════════════════════════════════════════
cat("【步骤 1】加载数据\n")

# 加载第1个插补数据集 (完整1124人)
data_path <- here:: here("outputs", "complete_data_with_ipw_imp1.rds")
load("variable_config.RData")

data_complete <- readRDS(data_path)

cat("  数据加载成功\n")
cat("  总样本数: ", nrow(data_complete), "\n")
cat("  变量数:    ", ncol(data_complete), "\n\n")

# 验证关键变量
required_vars <- c("has_fu_echo", "LVEDV_baseline", "EF_baseline")

if (!all(required_vars %in% names(data_complete))) {
  missing <- setdiff(required_vars, names(data_complete))
  stop("❌ 缺少关键变量: ", paste(missing, collapse = ", "))
}

# 统计三组样本量
n_total <- nrow(data_complete)
n_analysis <- sum(data_complete$has_fu_echo == "Yes", na.rm = TRUE)
n_excluded <- sum(data_complete$has_fu_echo == "No", na.rm = TRUE)

cat("  队列划分:\n")
cat("    总人群:    N = ", n_total, "\n", sep = "")
cat("    分析人群:  N = ", n_analysis, " (", round(n_analysis/n_total*100, 1), "%)\n", sep = "")
cat("    剔除人群: N = ", n_excluded, " (", round(n_excluded/n_total*100, 1), "%)\n\n", sep = "")

## ═══════════════════════════════════════════════════════
## 步骤 2: 定义Table 1包含的40个核心变量
## ═══════════════════════════════════════════════════════
cat("【步骤 2】定义核心变量\n")

## 变量分组定义
table1_vars <- list(
  
  ## 【1】人口学特征 (5个)
  demographics = c(
    "age",           # 年龄
    "gender",        # 性别
    "BMI",           # 体质指数
    "resident",      # 居住地
    "Career"         # 职业
  ),
  
  ## 【2】危险因素 (3个)
  risk_factors = c(
    "DM",            # 糖尿病
    "hypertension",  # 高血压
    "smoking"        # 吸烟
  ),
  
  ## 【3】心梗特征 (5个)
  mi_characteristics = c(
    "STEMI",         # ST段抬高型心梗
    "pPCI",          # 急诊PCI
    "ST_dev",        # ST段抬高
    "IN_killip",     # 入院Killip分级
    "GRACE_in"       # 入院GRACE评分
  ),
  
  ## 【4】基线心功能 (6个)
  cardiac_function = c(
    "EF_baseline",      # 基线射血分数
    "LVEDV_baseline",   # 基线左室舒张末容积
    "LVESV_baseline",   # 基线左室收缩末容积
    "EF_fu",            # 随访射血分数
    "LVEDV_fu",         # 随访左室舒张末容积
    "LVESV_fu"          # 随访左室收缩末容积
  ),
  
  ## 【5】心肌标志物 (3个)
  cardiac_markers = c(
    "cTnIpeak",         # 峰值肌钙蛋白I
    "NTproBNP_peak",    # 峰值NT-proBNP
    "CKMB"              # 肌酸激酶同工酶
  ),
  
  ## 【6】血常规 (3个)
  blood_routine = c(
    "WBC",           # 白细胞
    "HGB",           # 血红蛋白
    "PLT"            # 血小板
  ),
  
  ## 【7】肾功能 (2个)
  renal_function = c(
    "SCR",           # 血肌酐
    "EGFR"           # 肾小球滤过率
  ),
  
  ## 【8】炎症指标 (1个)
  inflammation = c(
    "CRP"            # C反应蛋白
  ),
  
  ## 【9】主要金属 (5个) - 使用原始值展示
  metals_primary = c(
    "Cu",            # 铜
    "Zn",            # 锌
    "Fe",            # 铁
    "Se",            # 硒
    "Pb"             # 铅
  ),
  
  ## 【10】用药 (5个)
  medications = c(
    "Aspirin",       # 阿司匹林
    "Statin",        # 他汀
    "ACEIorARB",     # ACEI/ARB
    "β_block",       # β受体阻滞剂
    "Clopidogrel"    # 氯吡格雷
  )
)

## 展平为向量
all_table1_vars <- unlist(table1_vars, use.names = FALSE)

## 检查变量存在性
vars_exist <- intersect(all_table1_vars, names(data_complete))
vars_missing <- setdiff(all_table1_vars, names(data_complete))

cat("  计划包含变量: ", length(all_table1_vars), "个\n", sep = "")
cat("  实际存在变量: ", length(vars_exist), "个\n", sep = "")

if (length(vars_missing) > 0) {
  cat("  ⚠️ 缺失变量: ", paste(vars_missing, collapse = ", "), "\n")
}

cat("\n")

## ═══════════════════════════════════════════════════════
## 步骤 3: 定义变量标签 (中英文对照)
## ═══════════════════════════════════════════════════════
cat("【步骤 3】定义变量标签\n")

## 变量标签映射表
label_mapping <- data.frame(
  variable = character(),
  label_cn = character(),
  stringsAsFactors = FALSE
)

## 【1】人口学
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("age", "gender", "BMI", "resident", "Career"),
  label_cn = c("年龄 (岁)", "男性", "体质指数 (kg/m²)", "城镇居民", "在职")
))

## 【2】危险因素
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("DM", "hypertension", "smoking"),
  label_cn = c("糖尿病", "高血压", "吸烟")
))

## 【3】心梗特征
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("STEMI", "pPCI", "ST_dev", "IN_killip", "GRACE_in"),
  label_cn = c("ST段抬高型心肌梗死", "急诊PCI", "ST段抬高", "入院Killip分级", "入院GRACE评分")
))

## 【4】基线心功能
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("EF_baseline", "LVEDV_baseline", "LVESV_baseline",
               "EF_fu", "LVEDV_fu", "LVESV_fu"),
  label_cn = c("基线左室射血分数 (%)", "基线左室舒张末容积 (mL)", "基线左室收缩末容积 (mL)",
               "随访左室射血分数 (%)", "随访左室舒张末容积 (mL)", "随访左室收缩末容积 (mL)")
))

## 【5】心肌标志物
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("cTnIpeak", "NTproBNP_peak", "CKMB"),
  label_cn = c("峰值肌钙蛋白I (ng/mL)", "峰值NT-proBNP (pg/mL)", "肌酸激酶同工酶 (U/L)")
))

## 【6】血常规
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("WBC", "HGB", "PLT"),
  label_cn = c("白细胞计数 (×10⁹/L)", "血红蛋白 (g/L)", "血小板计数 (×10⁹/L)")
))

## 【7】肾功能
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("SCR", "EGFR"),
  label_cn = c("血肌酐 (μmol/L)", "肾小球滤过率 (mL/min/1.73m²)")
))

## 【8】炎症
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("CRP"),
  label_cn = c("C反应蛋白 (mg/L)")
))

## 【9】金属
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("Cu", "Zn", "Fe", "Se", "Pb"),
  label_cn = c("血清铜 (μg/L)", "血清锌 (μg/L)", "血清铁 (μg/L)", 
               "血清硒 (μg/L)", "血清铅 (μg/L)")
))

## 【10】用药
label_mapping <- rbind(label_mapping, data.frame(
  variable = c("Aspirin", "Statin", "ACEIorARB", "β_block", "Clopidogrel"),
  label_cn = c("阿司匹林", "他汀类药物", "ACEI/ARB", "β受体阻滞剂", "氯吡格雷")
))

cat("  ✓ 标签定义完成:  ", nrow(label_mapping), "个变量\n\n", sep = "")

## ═══════════════════════════════════════════════════════
## 步骤 4: 识别分类变量
## ═══════════════════════════════════════════════════════
cat("【步骤 4】识别分类变量\n")

## 定义分类变量
categorical_vars <- c(
  # 人口学
  "gender", "resident", "Career",
  # 危险因素
  "DM", "hypertension", "smoking",
  # 心梗特征
  "STEMI", "pPCI", "ST_dev",
  # 用药
  "Aspirin", "Statin", "ACEIorARB", "β_block", "Clopidogrel"
)

# 有序分类变量
ordinal_vars <- c("IN_killip")

categorical_vars_exist <- intersect(categorical_vars, vars_exist)
ordinal_vars_exist <- intersect(ordinal_vars, vars_exist)

cat("  分类变量: ", length(categorical_vars_exist), "个\n", sep = "")
cat("  有序变量: ", length(ordinal_vars_exist), "个\n\n", sep = "")

## ═══════════════════════════════════════════════════════
## 步骤 5: 数据预处理
## ═══════════════════════════════════════════════════════
cat("【步骤 5】数据预处理\n")

## 创建分组变量
data_complete <- data_complete %>%
  mutate(
    cohort = case_when(
      has_fu_echo == "Yes" ~ "Analysis",
      has_fu_echo == "No" ~ "Excluded",
      TRUE ~ NA_character_
    ),
    cohort = factor(cohort, levels = c("Analysis", "Excluded"))
  )

## 确保分类变量是因子
for (var in categorical_vars_exist) {
  if (! is.factor(data_complete[[var]])) {
    data_complete[[var]] <- as.factor(data_complete[[var]])
  }
}

## 处理有序变量
for (var in ordinal_vars_exist) {
  if (!is.ordered(data_complete[[var]])) {
    data_complete[[var]] <- factor(data_complete[[var]], ordered = TRUE)
  }
}

cat("  ✓ 数据预处理完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 6: 生成TableOne对象
## ═══════════════════════════════════════════════════════
cat("【步骤 6】生成TableOne\n")

## 创建TableOne
tab1 <- CreateTableOne(
  vars = vars_exist,
  strata = "cohort",
  data = data_complete,
  factorVars = categorical_vars_exist,
  test = TRUE,           # 进行统计检验
  smd = TRUE,            # 计算SMD
  addOverall = TRUE      # 添加总体列
)

cat("  ✓ TableOne创建完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 7: 提取并格式化表格
## ═══════════════════════════════════════════════════════
cat("【步骤 7】格式化表格\n")

## 提取表格矩阵
tab1_matrix <- print(
  tab1,
  showAllLevels = FALSE,     # 二分类变量只显示一个水平
  smd = TRUE,                # 包含SMD
  test = TRUE,               # 包含p值
  exact = "stage",           # 精确检验
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE
)

## 转换为数据框
tab1_df <- as.data.frame(tab1_matrix)

## 添加变量名列
tab1_df$Variable <- rownames(tab1_df)
rownames(tab1_df) <- NULL

## 重新排列列顺序
tab1_df <- tab1_df %>%
  select(Variable, Overall, Analysis, Excluded, p, SMD)

## 重命名列
colnames(tab1_df) <- c(
  "特征",
  paste0("总人群\n(N=", n_total, ")"),
  paste0("分析人群\n(N=", n_analysis, ")"),
  paste0("剔除人群\n(N=", n_excluded, ")"),
  "P值",
  "SMD"
)

cat("  ✓ 表格格式化完成\n")
cat("    行数: ", nrow(tab1_df), "\n", sep = "")
cat("    列数: ", ncol(tab1_df), "\n\n", sep = "")

## ═══════════════════════════════════════════════════════
## 步骤 8: 替换变量标签
## ═══════════════════════════════════════════════════════
cat("【步骤 8】替换中文标签\n")

## 创建标签映射
for (i in 1:nrow(label_mapping)) {
  var_name <- label_mapping$variable[i]
  var_label <- label_mapping$label_cn[i]
  
  # 替换完全匹配的变量名
  tab1_df$特征[tab1_df$特征 == var_name] <- var_label
  
  # 替换包含变量名的行 (如 "gender = Male")
  tab1_df$特征 <- gsub(
    pattern = paste0("^", var_name, " "),
    replacement = paste0(var_label, " "),
    x = tab1_df$特征
  )
}

cat("  ✓ 标签替换完成\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 9: 添加分组小标题
## ═══════════════════════════════════════════════════════
cat("【步骤 9】添加分组小标题\n")

## 定义小标题
group_headers <- list(
  list(title = "人口学特征", vars = table1_vars$demographics),
  list(title = "危险因素", vars = table1_vars$risk_factors),
  list(title = "心肌梗死特征", vars = table1_vars$mi_characteristics),
  list(title = "基线心功能", vars = c(table1_vars$cardiac_function[1:3])),
  list(title = "随访心功能", vars = c(table1_vars$cardiac_function[4:6])),
  list(title = "心肌标志物", vars = table1_vars$cardiac_markers),
  list(title = "血常规", vars = table1_vars$blood_routine),
  list(title = "肾功能", vars = table1_vars$renal_function),
  list(title = "炎症指标", vars = table1_vars$inflammation),
  list(title = "血清金属元素", vars = table1_vars$metals_primary),
  list(title = "用药情况", vars = table1_vars$medications)
)

## 在每组第一个变量前插入小标题行
tab1_final <- data.frame()

for (group in group_headers) {
  # 获取该组的中文标签
  group_labels <- label_mapping$label_cn[label_mapping$variable %in% group$vars]
  
  if (length(group_labels) == 0) next
  
  # 找到该组第一个变量在表格中的位置
  first_var_row <- which(grepl(group_labels[1], tab1_df$特征, fixed = TRUE))[1]
  
  if (is.na(first_var_row)) next
  
  # 如果还没有添加数据,先添加之前的行
  if (nrow(tab1_final) == 0) {
    if (first_var_row > 1) {
      tab1_final <- tab1_df[1:(first_var_row - 1), ]
    }
  } else {
    # 添加从上次结束到当前组开始的行
    last_row <- nrow(tab1_final)
    prev_last_var <- tab1_df$特征[last_row]
    prev_last_row_in_orig <- which(tab1_df$特征 == prev_last_var)[1]
    
    if (first_var_row > prev_last_row_in_orig + 1) {
      tab1_final <- rbind(
        tab1_final,
        tab1_df[(prev_last_row_in_orig + 1):(first_var_row - 1), ]
      )
    }
  }
  
  # 插入小标题行
  header_row <- tab1_df[first_var_row, ]
  header_row$特征 <- paste0("**", group$title, "**")
  header_row[, 2: 6] <- ""
  
  tab1_final <- rbind(tab1_final, header_row)
  
  # 添加该组的变量
  group_rows <- which(sapply(group_labels, function(label) {
    any(grepl(label, tab1_df$特征, fixed = TRUE))
  }))
  
  for (label in group_labels) {
    matching_rows <- which(grepl(label, tab1_df$特征, fixed = TRUE))
    if (length(matching_rows) > 0) {
      tab1_final <- rbind(tab1_final, tab1_df[matching_rows, ])
    }
  }
}

## 添加剩余的行
if (nrow(tab1_final) < nrow(tab1_df)) {
  remaining_rows <- setdiff(1:nrow(tab1_df), 
                            which(tab1_df$特征 %in% tab1_final$特征))
  if (length(remaining_rows) > 0) {
    tab1_final <- rbind(tab1_final, tab1_df[remaining_rows, ])
  }
}

cat("  ✓ 分组小标题已添加\n\n")

## ═══════════════════════════════════════════════════════
## 步骤 10: 生成Word文档
## ═══════════════════════════════════════════════════════
cat("【步骤 10】生成Word文档\n")

## 创建flextable
ft <- flextable(tab1_final)

## 设置表头
ft <- set_header_labels(
  ft,
  特征 = "特征",
  .   = names(tab1_final)[2:6]
)

## 设置列宽
ft <- width(ft, j = 1, width = 3)    # 特征列
ft <- width(ft, j = 2:4, width = 1.5) # 数据列
ft <- width(ft, j = 5, width = 0.8)  # P值
ft <- width(ft, j = 6, width = 0.8)  # SMD

## 设置字体
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- font(ft, fontname = "SimSun", part = "body", 
           i = ~ grepl("[\u4e00-\u9fa5]", 特征))  # 中文用宋体

## 设置对齐
ft <- align(ft, align = "left", part = "body", j = 1)
ft <- align(ft, align = "center", part = "body", j = 2:6)
ft <- align(ft, align = "center", part = "header")

## 三线表样式
ft <- border_remove(ft)
ft <- hline_top(ft, border = fp_border(color = "black", width = 2), part = "header")
ft <- hline_bottom(ft, border = fp_border(color = "black", width = 1), part = "header")
ft <- hline_bottom(ft, border = fp_border(color = "black", width = 2), part = "body")

## 小标题行加粗
ft <- bold(ft, i = ~ grepl("^\\*\\*", 特征), j = 1)

## 移除小标题行的星号
tab1_final$特征 <- gsub("\\*\\*", "", tab1_final$特征)
ft <- flextable(tab1_final)
ft <- compose(ft, i = ~ grepl("^人口学|^危险|^心肌|^基线|^随访|^血|^肾|^炎|^用药", 特征),
              j = "特征",
              value = as_paragraph(as_b(特征)))

## 添加脚注
ft <- add_footer_lines(ft, values = c(
  "注: 连续变量表示为均值 ± 标准差或中位数 [四分位距]; 分类变量表示为例数 (百分比)。",
  "P值基于t检验、Mann-Whitney U检验或χ²检验。",
  "SMD = 标准化均数差 (Standardized Mean Difference); SMD < 0.1 表示组间差异可忽略。",
  "缺失数据已通过多重插补处理 (m=20)。",
  paste0("分析人群: 有完整随访超声数据的患者 (N=", n_analysis, ")。"),
  paste0("剔除人群: 缺失随访超声数据的患者 (N=", n_excluded, ")。")
))

ft <- fontsize(ft, size = 9, part = "footer")
ft <- align(ft, align = "left", part = "footer")

## 自动调整
ft <- autofit(ft)

## 保存为Word
output_file <- here::here("outputs", "Table1_Baseline_Core.docx")

doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)

print(doc, target = output_file)

cat("  ✓ Word文档已保存\n")
cat("    文件路径: ", output_file, "\n\n", sep = "")

## ═══════════════════════════════════════════════════════
## 步骤 11: 生成CSV备份
## ═══════════════════════════════════════════════════════
cat("【步骤 11】保存CSV备份\n")

csv_file <- here::here("outputs", "Table1_Baseline_Core.csv")

write.csv(tab1_final, 
           file = csv_file,
           row.names = FALSE,
           fileEncoding = "UTF-8")

cat("  ✓ CSV文件已保存\n")
cat("    文件路径:  ", csv_file, "\n\n", sep = "")

## ═══════════════════════════════════════════════════════
## 步骤 12: 生成统计摘要
## ═══════════════════════════════════════════════════════
cat("【步骤 12】生成统计摘要\n")

## 计算显著差异的变量数
sig_vars <- sum(as.numeric(gsub("[<>]", "", tab1_final$P值)) < 0.05, na.rm = TRUE)

## 计算高SMD的变量数
high_smd <- sum(as.numeric(tab1_final$SMD) > 0.1, na.rm = TRUE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("📊 Table 1 统计摘要\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("样本量:\n")
cat("  总人群:    ", n_total, "\n", sep = "")
cat("  分析人群: ", n_analysis, " (", round(n_analysis/n_total*100, 1), "%)\n", sep = "")
cat("  剔除人群:  ", n_excluded, " (", round(n_excluded/n_total*100, 1), "%)\n\n", sep = "")

cat("变量统计:\n")
cat("  总变量数:          ", length(vars_exist), "\n", sep = "")
cat("  连续变量:          ", length(vars_exist) - length(categorical_vars_exist), "\n", sep = "")
cat("  分类变量:         ", length(categorical_vars_exist), "\n", sep = "")
cat("  显著差异变量 (p<0.05): ", sig_vars, "\n", sep = "")
cat("  高SMD变量 (>0.1):      ", high_smd, "\n\n", sep = "")

cat("关键发现:\n")

## 识别SMD最大的5个变量
tab1_smd <- tab1_final %>%
  filter(! is.na(SMD), SMD != "") %>%
  mutate(SMD_num = as.numeric(SMD)) %>%
  arrange(desc(SMD_num)) %>%
  head(5)

if (nrow(tab1_smd) > 0) {
  cat("  组间差异最大的5个变量:\n")
  for (i in 1:min(5, nrow(tab1_smd))) {
    cat("    ", i, ".  ", tab1_smd$特征[i], " (SMD = ", 
        sprintf("%.3f", tab1_smd$SMD_num[i]), ")\n", sep = "")
  }
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n\n")

## ═══════════════════════════════════════════════════════
## 完成
## ═══════════════════════════════════════════════════════
cat("╔═══════════════════════════════════════════════════════╗\n")
cat("║            Table 1 生成完成！                        ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

cat("输出文件:\n")
cat("  1. Table1_Baseline_Core. docx  - Word格式三线表\n")
cat("  2. Table1_Baseline_Core.csv   - CSV备份\n\n")

cat("下一步:\n")
cat("  1. 检查Word文档格式\n")
cat("  2. 生成Figure S1 (IPW平衡性Love Plot)\n")
cat("  3. 进行主要结局分析\n\n")

cat("完成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")