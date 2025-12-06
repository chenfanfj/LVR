############################################################
## 文件名: missing_diagnostics_LVR.R
## 目的：
##   对 AMI_PCI_metal_clean3.xlsx（纳排后 n = 1124）进行
##   缺失数据的初步诊断，重点关注 LVEDV_fu（随访 LVEDV，
##   主要结局，缺失较多）以及其他变量，在“未进行任何插补
##   (MI)”之前完成。
##
##   输出内容：
##     - 缺失模式热图 & 各变量缺失比例柱状图（保存到 plots/）
##     - 按 has_fu_echo（是否有随访超声/LVEDV）分层的基线对比表
##       （保存到 outputs/）
##     - Little 的 MCAR 检验结果
##     - 以 has_fu_echo 为因变量的逻辑回归模型结果
##
## 使用的核心 R 包：
##   - naniar, visdat：缺失数据可视化（GitHub/CRAN 上维护活跃、
##     被广泛使用）
##   - mice 或 MissMech：Little 的 MCAR 检验（缺失数据方法学中的
##     标准工具）
##   - gtsummary：分层描述性表格（带 p 值，适合发表）
############################################################

## ---------------------------
## 1. 环境设置 & 数据读取
## ---------------------------

## 安装并加载需要的 R 包（若尚未安装则先 install.packages）
pkgs <- c(
  "here",       # 项目根目录路径
  "readxl",     # 读取 Excel
  "dplyr",      # 数据整理
  "naniar",     # 缺失可视化和辅助函数
  "visdat",     # 数据类型与缺失可视化
  "gtsummary",  # 分层描述性表格与 p 值
  "broom",      # 模型结果整洁输出
  "mice",       # Little 检验辅助（md.pattern 等）
  "VIM",        # 另一套缺失模式可视化工具
  "glmnet"      # 非必须，这里保持简洁，可视情况去掉
)
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

## 若想用 MissMech 的 Little MCAR 检验：
if (!requireNamespace("MissMech", quietly = TRUE)) {
  install.packages("MissMech")
}
library(MissMech)

## 设置项目根目录（假定 R project 在 D:/R_try/2026_APMM）
## 若在该 Rproj 内运行，here::here() 会自动识别根目录
root_dir <- here::here()
cat("检测到的项目根目录为:", root_dir, "\n")

## 定义数据路径
data_path <- here::here("data", "AMI_PCI_metal_clean3.xlsx")

## 读取 Excel 数据
dat_raw <- readxl::read_excel(
  path      = data_path,
  col_names = TRUE
)

## 可选：检查样本量
cat("数据集行数(样本数):", nrow(dat_raw), "\n")
## 期望为 1124

## 工作副本
dat <- dat_raw

## 第一列为 ID，记录变量名
id_var <- names(dat)[1]
cat("假定 ID 变量为:", id_var, "\n")

## 创建 has_fu_echo: 是否有随访 LVEDV（即 LVEDV_fu 是否非缺失）
if (!"LVEDV_fu" %in% names(dat)) {
  stop("数据集中未找到 LVEDV_fu 变量，请检查变量名。")
}

dat <- dat %>%
  mutate(
    has_fu_echo = if_else(!is.na(LVEDV_fu), "Yes", "No")
  )

## 将 has_fu_echo 转为因子，便于后续分层分析
dat <- dat %>%
  mutate(
    has_fu_echo = factor(has_fu_echo, levels = c("Yes", "No"))
  )

## ---------------------------
## 2. 缺失情况可视化
## ---------------------------
## 目标：
##   - 用热图展示各变量缺失模式
##   - 用柱状图展示每个变量缺失比例
##   - 特别关注 LVEDV_fu 缺失与其他变量缺失的共现

## 确保 plots/ 与 outputs/ 目录存在
plots_dir   <- here::here("plots")
outputs_dir <- here::here("outputs")
if (!dir.exists(plots_dir))   dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE)

## 2.1 使用 naniar 绘制缺失热图（vis_miss）
## 去掉 ID 列，只看其他变量
dat_no_id <- dat %>% select(-all_of(id_var))

## 缺失热图：白色 = 有值, 红色 = 缺失；根据缺失量排序变量
p_heat <- naniar::vis_miss(dat_no_id, sort_miss = TRUE) +
  ggplot2::labs(
    title    = "变量缺失模式热图",
    subtitle = "白色 = 非缺失, 红色 = 缺失"
  )

ggplot2::ggsave(
  filename = file.path(plots_dir, "missingness_heatmap_all_vars.png"),
  plot     = p_heat,
  width    = 12, height = 8, dpi = 300
)

## 2.2 每个变量缺失数/比例的柱状图：naniar::gg_miss_var
p_bar <- naniar::gg_miss_var(dat_no_id) +
  ggplot2::labs(
    title = "各变量缺失数量/比例",
    x     = "变量",
    y     = "缺失观测数"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

ggplot2::ggsave(
  filename = file.path(plots_dir, "missingness_bar_all_vars.png"),
  plot     = p_bar,
  width    = 12, height = 8, dpi = 300
)

cat("缺失热图和柱状图已保存至:", plots_dir, "\n")

## ---------------------------
## 3. 按 has_fu_echo 分层的基线特征对比
## ---------------------------
## 目标：
##   - 比较有随访超声 (has_fu_echo = Yes) 与 无随访超声 (No)
##     两组在基线特征上的差异。
##   - 使用 gtsummary::tbl_summary(by = has_fu_echo)，并通过 add_p()
##     计算组间差异的统计检验（t 检验/Wilcoxon、卡方）。

## 给定的变量列表：
bin_vars <- c(
  "resident", "Career", "gender", "smoking", "DM", "hypertension",
  "Cancer", "his_stroke", "AF", "pPCI", "VF", "Cardio_shock",
  "Temporary_pacemaker", "STEMI", "Inferior_MI", "Anterior_MI",
  "EA_MI", "AT_MI", "Posterior_MI", "Lateral_MI", "HL_MI", "RV_MI",
  "LM", "TL_LM", "LAD", "TL_LAD", "LCX", "TL_LCX", "RCA", "TL_RCA",
  "Ticagrelor", "Heparin", "Clopidogrel", "Aspirin", "Statin",
  "ACEIorARB", "β_block", "CCB", "diuretics", "insulin", "metformin",
  "tirofiban", "Cardiac_arrest_in", "ST_dev", "ST_dep", "Thrombo_Burden"
)
ord_vars <- c(
  "IN_killip", "OUT_killip", "GRACE_in_str", "Grace_out_str"
)
con_vars <- c(
  "age", "height", "weight", "BMI", "SBP", "DBP", "HR", "Stent_no",
  "Lesion_no", "GRACE_in", "GRACE_out", "WBC", "NE_", "LY_", "MO_",
  "EO_", "BA_", "NE", "LY", "EO", "BA", "RBC", "HGB", "HCT", "MCV",
  "MCH", "MCHC", "RDW", "PLT", "MPV", "PCT", "PDW", "MO",
  "NTproBNP_baseline", "NTproBNP_peak", "cTnI_baseline", "cTnIpeak",
  "cTnI_out", "FT3", "FT4", "S_TSH", "h_CT", "PAB_SH", "AST_ALT",
  "BUN_CREA", "CO2", "GGT", "ALT", "LDH", "LDL", "AST", "BUN", "URIC",
  "CHOL", "TBIL", "TP", "PHOS", "CL", "OSM", "GLB", "TG", "A_G", "Alb",
  "DBIL", "ALP", "SCR", "EGFR", "CK", "CKMB", "GLU", "APROA",
  "APROB_APROA", "APROB", "Ca", "Na", "K", "Mg", "IBIL", "AG", "HDL",
  "HbAlc", "DD", "PT_sec", "TT", "INR", "ATⅢ", "APTT", "FIB", "IL_6",
  "CRP", "EF_baseline", "LVEDV_baseline", "LVESV_baseline",
  "EF_fu", "LVEDV_fu", "LVESV_fu", "ΔLVEDV", "Echo_fu_day",
  "Echo_fu_month", "Al", "As", "B", "Ba", "Cr", "Cu", "Fe", "Li", "Mn",
  "Mo", "Ni", "Pb", "Rb", "Sb", "Se", "Sn", "Sr", "V", "Zn"
)

## 对“基线对比”而言，我们希望排除明显的随访变量
fu_vars <- c("EF_fu", "LVEDV_fu", "LVESV_fu", "ΔLVEDV",
             "Echo_fu_day", "Echo_fu_month")

baseline_vars <- setdiff(c(bin_vars, ord_vars, con_vars), fu_vars)
baseline_vars <- intersect(baseline_vars, names(dat))

## 把二元与有序变量转换为 factor，以便 gtsummary 正确识别为分类变量
dat_comp <- dat
for (v in bin_vars) {
  if (v %in% names(dat_comp)) dat_comp[[v]] <- factor(dat_comp[[v]])
}
for (v in ord_vars) {
  if (v %in% names(dat_comp)) dat_comp[[v]] <- factor(dat_comp[[v]], ordered = TRUE)
}

## 使用 gtsummary::tbl_summary，按 has_fu_echo 分层
tbl_comp <- dat_comp %>%
  select(all_of(c("has_fu_echo", baseline_vars))) %>%
  gtsummary::tbl_summary(
    by        = has_fu_echo,
    missing   = "ifany",
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",      # 连续变量：均值 ± 标准差
      all_categorical() ~ "{n} ({p}%)"         # 分类变量：n (百分比)
    )
  ) %>%
  gtsummary::add_p(
    test = list(
      all_continuous() ~ "t.test",   # 或根据分布改为 "wilcox.test"
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)
  ) %>%
  gtsummary::modify_header(label = "**变量**") %>%
  gtsummary::bold_labels()

## 控制台查看（gt 对象）
tbl_comp

## 保存 Word 与 CSV
comp_docx_path <- file.path(outputs_dir, "baseline_by_has_fu_echo.docx")
comp_csv_path  <- file.path(outputs_dir, "baseline_by_has_fu_echo.csv")

gtsummary::tbl_save(
  x        = tbl_comp,
  filename = comp_docx_path
)

tbl_comp_df <- as.data.frame(tbl_comp, stringsAsFactors = FALSE)
write.csv(tbl_comp_df, file = comp_csv_path, row.names = FALSE)

cat("按 has_fu_echo 分层的基线对比表已保存至:\n",
    comp_docx_path, "\n", comp_csv_path, "\n")

## ---------------------------
## 4. 高级诊断
## ---------------------------

## 4.1 Little 的 MCAR 检验
## 目标：
##   - 统计学上检验“缺失模式是否符合 MCAR（完全随机缺失）”。
##   - 若 p < 0.05，则拒绝 MCAR，说明缺失模式显著偏离 MCAR。

## Little 检验对维度较敏感，通常选取重要的连续/近似正态变量子集
mcar_vars <- c("age", "gender", "DM", "hypertension",
               "EF_baseline", "LVEDV_baseline", "LVESV_baseline")

mcar_vars <- intersect(mcar_vars, names(dat_comp))

dat_mcar <- dat_comp %>%
  select(all_of(mcar_vars))

## MissMech::TestMCARNormality 要求矩阵为数值型
dat_mcar_num <- dat_mcar
for (v in mcar_vars) {
  if (is.factor(dat_mcar_num[[v]]) || is.character(dat_mcar_num[[v]])) {
    dat_mcar_num[[v]] <- as.numeric(dat_mcar_num[[v]])
  }
}

## 使用 MissMech::TestMCARNormality 做 Little MCAR 检验
## 注意：该检验假定多元正态；在真实数据中仅是近似判断。
mcar_test <- MissMech::TestMCARNormality(
  as.matrix(dat_mcar_num),
  imputed = FALSE
)

print(mcar_test)
## 解释（在注释中）：
##   - 检查输出中的 p.value。
##   - 若 p < 0.05，则拒绝 MCAR，说明整体缺失模式与完全随机缺失不符。
##   - 若 p >= 0.05，则不能拒绝 MCAR（不代表一定是 MCAR，也可能是 MAR/MNAR）。

## 4.2 逻辑回归：has_fu_echo ~ 基线协变量
## 目标：
##   - 判断 LVEDV_fu 缺失（has_fu_echo 为 No）是否与观测到的
##     基线变量相关，从而初步判断是否为 MAR。
##   - 以 has_fu_echo 为因变量，拟合简单的 logistic 回归。

## 选取一些重要、临床相关的基线预测变量
logit_vars <- c("age", "gender", "DM", "hypertension",
                "EF_baseline", "LVEDV_baseline", "LVESV_baseline")

logit_vars <- intersect(logit_vars, names(dat_comp))

## 构建分析数据集，去除因变量或自变量缺失的记录
logit_dat <- dat_comp %>%
  select(all_of(c("has_fu_echo", logit_vars))) %>%
  filter(!is.na(has_fu_echo))

## 拟合 logistic 回归模型：has_fu_echo ~ 预测变量
## 注意：has_fu_echo 因子水平为 Yes/No，glm(family=binomial) 默认
## 以第一个水平为参照。
logit_formula <- as.formula(
  paste("has_fu_echo ~", paste(logit_vars, collapse = " + "))
)

fit_logit <- glm(logit_formula, data = logit_dat, family = binomial)

summary(fit_logit)

## 计算伪 R²（如 McFadden / Nagelkerke），看模型解释度
if (!requireNamespace("pscl", quietly = TRUE)) {
  install.packages("pscl")
}
library(pscl)
pseudo_r2 <- pscl::pR2(fit_logit)
cat("伪 R² (McFadden / Nagelkerke 等):\n")
print(pseudo_r2)

## 使用 broom::tidy 以 OR 形式输出系数
fit_logit_tidy <- broom::tidy(fit_logit, exponentiate = TRUE, conf.int = TRUE)
print(fit_logit_tidy)
## 解释：
##   - OR（优势比）> 1 且 p < 0.05，说明该协变量与“有随访超声
##     (has_fu_echo = Yes)”显著相关。
##   - 若有多个协变量显著，且伪 R² 明显 > 0，则说明 LVEDV_fu 的
##     缺失情形很大程度上可由这些已观测的基线变量预测，
##     更加符合 MAR 机制。
##   - 若协变量均不显著、伪 R² 接近 0，则缺失可能更接近 MCAR，
##     或由未观测因素驱动（即 MNAR）。

## ---------------------------
## 5. 结果解读要点（摘要性注释）
## ---------------------------

# 1. 缺失热图与柱状图：
#    - 可直观查看哪些变量缺失比例较高；
#    - 可观察 LVEDV_fu 缺失是否与其他重要变量的缺失共现。
#
# 2. 按 has_fu_echo 分层的基线对比表：
#    - 若在年龄、合并症、病情严重程度（如 EF_baseline）等基线变量上
#      双组有显著差异（p < 0.05），说明“谁有随访超声”并非随机，
#      而是与这些已观测协变量相关。
#    - 这更符合条件 MAR（给定协变量后缺失可视为随机），提示后续
#      分析中应采用能够处理“结局缺失非随机”的方法：如对 LVEDV_fu
#      使用 IPW（逆概率加权），并对基线协变量使用 MI。
#
# 3. Little 的 MCAR 检验：
#    - 若 p < 0.05，则整体缺失模式显著偏离 MCAR，进一步支持“不应
#      假定完全随机缺失”的结论；
#    - 若 p >= 0.05，则不能拒绝 MCAR，但并不排除 MAR/MNAR 的可能，
#      仍需结合其他诊断信息。
#
# 4. 逻辑回归 has_fu_echo ~ 基线协变量：
#    - 若年龄、EF_baseline 等基线变量的 OR 显著 > 1 或 < 1，
#      且伪 R² 非常 > 0，则说明“是否有随访超声”在统计上可被
#      这些基线因素预测，缺失机制更接近 MAR。
#    - 若所有系数基本不显著，且伪 R² 很小，缺失可能更接近 MCAR，
#      或是由未观测变量驱动的 MNAR。
#
# 综上：
#   - 这些诊断是“插补和 IPW 之前”的探索性步骤；
#   - 最终如何建模缺失（例如 IPW 模型中加入哪些协变量，
#     MI 模型中纳入哪些预测变量）应综合这些结果来确定：
#       * 哪些基线变量显著预测 has_fu_echo；
#       * 哪些变量缺失比例较高；
#       * Little 检验是否否定 MCAR。
#   - 在本 LVR 研究情境下：
#       * 若 MCAR 被否定，且逻辑回归表明缺失与基线变量强相关，
#         则应在后续分析中充分利用 IPW 处理 LVEDV_fu 缺失，并对
#         基线协变量、金属和代谢组学变量进行 MI（插补模型中
#         应包括预测缺失的变量），以减少偏倚并提高精度。
############################################################