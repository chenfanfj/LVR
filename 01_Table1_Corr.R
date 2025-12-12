# -----------------------------------------------------------------------------
# 脚本名称: 01_Table1_Corr.R
# 功能: 生成基线特征三线表 (Table 1) 和 金属相关性分析
# 输入: mice_imputation_final.rds (插补数据), imputed_data_with_ipw_weights.rds (权重数据)
# 输出: Table1_Baseline.docx, Correlation_Plot.png
# -----------------------------------------------------------------------------

# 1. 加载必要的包
# -----------------------------------------------------------------------------
if(!require(pacman)) install.packages("pacman")
pacman::p_load(
  mice,       # 处理插补数据
  tableone,   # 生成基线表
  survey,     # 处理加权数据
  flextable,  # 生成Word三线表
  officer,    # Word导出辅助
  dplyr,      # 数据清洗
  corrplot,   # 相关性绘图
  Hmisc       # 相关性计算
)

# 2. 数据加载与准备
# -----------------------------------------------------------------------------
# 读取插补后的权重数据 (包含 sw_trunc 权重列)
# 注意：这里假设您有一个文件包含插补后的完整数据和权重。
# 如果没有，请使用: data_imp <- complete(readRDS("outputs/mice_imputation_final.rds"), 1)
# 然后 merge 权重。
# 下面使用第1个插补数据集作为展示代表 (Table 1 通常展示单次插补或原始数据，但加权需要完整数据)

# 假设路径，请根据实际情况修改
data_path <- "outputs/imputed_data_with_ipw_weights.rds" 

if(file.exists(data_path)){
  # 如果是RDS文件
  data_complete <- readRDS(data_path)
} else if(file.exists("outputs/imputed_data_with_ipw_weights.csv")){
  # 如果是CSV文件
  data_complete <- read.csv("outputs/imputed_data_with_ipw_weights.csv")
} else {
  # 回退方案：从MICE对象提取第1个插补集 (仅作演示)
  mi_object <- readRDS("outputs/mice_imputation_final.rds")
  data_complete <- complete(mi_object, 1)
  # 模拟权重 (实际请加载您的权重文件)
  data_complete$sw_trunc <- runif(nrow(data_complete), 0.5, 2.0) 
  warning("警告：使用了模拟权重，请确保加载了真实的 sw_trunc 权重！")
}

# 3. 定义变量标签与顺序 (模拟 variable lable and order.txt)
# -----------------------------------------------------------------------------
# 【用户请注意】：请在此处对照您的 txt 文件调整顺序和中文名称
# 格式： "变量名" = "中文标签"

var_map <- c(
  # --- 人口学特征 ---
  "age"             = "年龄 (岁)",
  "gender"          = "性别 (男性)",
  "BMI"             = "体重指数 (kg/m²)",
  "smoking"         = "吸烟史",
  "resident"        = "居住地 (城市)",
  
  # --- 既往史 ---
  "hypertension"    = "高血压",
  "DM"              = "糖尿病",
  "dyslipidemia"    = "血脂异常",
  
  # --- 入院与治疗 ---
  "Killip_class"    = "Killip 分级", # 假设变量名，请修正
  "STEMI"           = "STEMI 诊断",
  "pPCI"            = "行急诊 PCI",
  
  # --- 实验室检查 ---
  "cTnIpeak"        = "肌钙蛋白 I 峰值 (ng/mL)",
  "NTproBNP_peak"   = "NT-proBNP 峰值 (pg/mL)",
  "WBC"             = "白细胞计数 (×10⁹/L)",
  "HGB"             = "血红蛋白 (g/L)",
  "PLT"             = "血小板计数 (×10⁹/L)",
  "CRP"             = "C反应蛋白 (mg/L)",
  
  # --- 超声心动图 ---
  "EF_baseline"     = "基线 LVEF (%)",
  "LVEDV_baseline"  = "基线 LVEDV (mL)",
  
  # --- 环境金属暴露 (对数转换值或原始值，视分析而定) ---
  "log_Cd"          = "血镉 (ln ug/L)",
  "log_Pb"          = "血铅 (ln ug/L)",
  "log_As"          = "血砷 (ln ug/L)",
  "log_Cu"          = "血铜 (ln ug/L)",
  "log_Zn"          = "血锌 (ln ug/L)"
)

# 提取变量名列表
myVars <- names(var_map)
# 确保数据集中存在这些变量 (取交集)
myVars <- intersect(myVars, names(data_complete))
# 识别分类变量
catVars <- c("gender", "smoking", "resident", "hypertension", "DM", "dyslipidemia", "STEMI", "pPCI", "Killip_class")
catVars <- intersect(catVars, names(data_complete))

# 4. 构建三个队列数据
# -----------------------------------------------------------------------------
# 队列 1: 总人群 (Total Cohort)
cohort_total <- data_complete

# 队列 2: 分析人群 - 未加权 (Analysis Cohort - Unweighted)
# 假设有随访超声的标志是 has_fu_echo == 1 (参考 Source 390)
cohort_analysis <- subset(data_complete, has_fu_echo == 1)

# 队列 3: 分析人群 - 加权 (Analysis Cohort - Weighted IPW)
# 需要建立 survey design 对象
if("sw_trunc" %in% names(cohort_analysis)){
  svy_design <- svydesign(ids = ~1, data = cohort_analysis, weights = ~sw_trunc)
} else {
  stop("错误：分析数据集中未找到 'sw_trunc' 权重列。")
}

# 5. 使用 TableOne 生成数据
# -----------------------------------------------------------------------------
message("正在生成 Table 1...")

# (1) 总人群
tab1 <- CreateTableOne(vars = myVars, data = cohort_total, factorVars = catVars, test = FALSE)

# (2) 分析人群 (未加权)
tab2 <- CreateTableOne(vars = myVars, data = cohort_analysis, factorVars = catVars, test = FALSE)

# (3) 分析人群 (加权) - 计算 SMD 需要一点技巧，svyCreateTableOne 可以做
tab3 <- svyCreateTableOne(vars = myVars, data = svy_design, factorVars = catVars, test = FALSE)

# 6. 格式化与合并
# -----------------------------------------------------------------------------
# 自定义渲染函数：Mean ± SD
my_renderer <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  sprintf("%.1f ± %.1f", mean_val, sd_val)
}

# 提取打印矩阵 (TableOne 默认输出 Mean (SD), 这里我们后续在 dataframe 处理)
mat1 <- print(tab1, showAllLevels = FALSE, printToggle = FALSE, smd = FALSE)
mat2 <- print(tab2, showAllLevels = FALSE, printToggle = FALSE, smd = TRUE) # 这里的SMD是相对于谁？通常不显示
mat3 <- print(tab3, showAllLevels = FALSE, printToggle = FALSE, smd = TRUE) # 这里的SMD需要单独计算

# 提取 SMD (Standardized Mean Difference)
# 比较 "未加权分析组" vs "剔除组" (通常做法) 或者 "加权后" 是否平衡
# 这里为了简化，我们直接展示三列数值。SMD通常单独计算展示。

# 将矩阵转为数据框
df_final <- cbind(
  mat1[, 1, drop=FALSE], # Total
  mat2[, 1, drop=FALSE], # Analysis Unweighted
  mat3[, 1, drop=FALSE]  # Analysis Weighted
)
colnames(df_final) <- c("总队列 (N=1344)", "分析队列 (未加权, N=661)", "分析队列 (IPW加权, N=661)")

# 转为 data.frame 并处理行名
df_out <- as.data.frame(df_final)
df_out$Variable <- rownames(df_final)
rownames(df_out) <- NULL
df_out <- df_out[, c("Variable", "总队列 (N=1344)", "分析队列 (未加权, N=661)", "分析队列 (IPW加权, N=661)")]

# 替换中文标签
# 这是一个简单的替换循环
for(eng_name in names(var_map)){
  # 替换变量名行
  df_out$Variable <- gsub(eng_name, var_map[eng_name], df_out$Variable)
}
# 替换 levels (如 0/1) 为是/否，如果需要，可以在这里手动替换
df_out$Variable <- gsub(" = 1 \\(100\\.0\\%\\)", "", df_out$Variable) # 去除多余的标签

# 7. 导出为 Word 三线表 (使用 flextable)
# -----------------------------------------------------------------------------
ft <- flextable(df_out)

# 设置表头和格式
ft <- set_header_labels(ft, Variable = "特征")
ft <- autofit(ft)
ft <- font(ft, fontname = "Times New Roman", part = "all") # 学术常用字体
ft <- font(ft, fontname = "SimSun", part = "body", i = ~ grepl("[\u4e00-\u9fa5]", Variable)) # 中文宋体

# 三线表样式 (APA 风格)
ft <- border_remove(ft)
ft <- hline_top(ft, border = fp_border(color = "black", width = 1.5))
ft <- hline(ft, i = 1, part = "header", border = fp_border(color = "black", width = 1.5)) 
ft <- hline_bottom(ft, part = "body", border = fp_border(color = "black", width = 1.5))

# 添加脚注
ft <- add_footer_lines(ft, values = "注：连续变量表示为 Mean ± SD；分类变量表示为 n (%)。IPW：逆概率加权。")

# 保存 Word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "outputs/Table1_Baseline.docx")

message("Table 1 已保存至 outputs/Table1_Baseline.docx")

# 8. 金属相关性分析 (Spearman)
# -----------------------------------------------------------------------------
message("正在进行金属相关性分析...")

# 定义金属变量列表 (请确保这些列名在数据中存在)
metal_vars <- c("log_Cd", "log_Pb", "log_As", "log_Cu", "log_Zn")
metal_vars <- intersect(metal_vars, names(data_complete))

if(length(metal_vars) > 1){
  # 提取金属数据
  metal_data <- data_complete[, metal_vars]
  
  # 计算 Spearman 相关系数
  M <- cor(metal_data, method = "spearman", use = "pairwise.complete.obs")
  
  # 保存图片
  png("outputs/Correlation_Plot.png", width = 2000, height = 2000, res = 300)
  
  # 绘制热图
  corrplot.mixed(M, 
                 lower = "number", 
                 upper = "circle",
                 tl.col = "black",
                 tl.pos = "lt",
                 number.cex = 0.8,
                 title = "金属暴露 Spearman 相关性热图",
                 mar = c(0,0,2,0))
  dev.off()
  
  message("相关性热图已保存至 outputs/Correlation_Plot.png")
} else {
  warning("未找到足够的金属变量进行相关性分析，请检查变量名。")
}

# 结束
message("所有任务完成。")