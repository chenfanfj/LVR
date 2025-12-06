############################################################
## LVR_baseline_gtsummary.R – 使用 gtsummary 直接导出 Word
## 目的:
##   1) 读取 AMI_PCI_metal_clean3. xlsx (n=1124)
##   2) 使用 gtsummary 构建基线表
##   3) 应用中文标签和三线表格式
##   4) 直接导出到 outputs/ 目录的 Word (. docx) 文件
############################################################
rm(list = ls(all.names = TRUE))
gc()
## 加载所需包
pkgs <- c(
  "here",       # 项目相对路径
  "readxl",     # 读取 Excel
  "dplyr",      # 数据整理
  "gtsummary",  # 集成汇总表
  "flextable",  # 微调表格格式
  "officer"     # Word 文档导出
)

for (p in pkgs) {
  if (! requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

## ---------------------------------
## 1. 项目路径与数据导入
## ---------------------------------

data_path <- here::here("data", "AMI_PCI_metal_clean3.xlsx")
dat <- readxl::read_excel(path = data_path, col_names = TRUE)

## 可选:验证样本量
# stopifnot(nrow(dat) == 1124)

## ------------------------------------------------
## 2. 定义变量列表(使用您的确切变量名)
## ------------------------------------------------

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
  "age", "height", "weight", "BMI", "SBP", "DBP", "HR",
  "Stent_no", "Lesion_no", "GRACE_in", "GRACE_out",
  "WBC", "NE_", "LY_", "MO_", "EO_", "BA_", "NE", "LY", "EO", "BA",
  "RBC", "HGB", "HCT", "MCV", "MCH", "MCHC", "RDW", "PLT", "MPV",
  "PCT", "PDW", "MO",
  "NTproBNP_baseline", "NTproBNP_peak",
  "cTnI_baseline", "cTnIpeak", "cTnI_out",
  "FT3", "FT4", "S_TSH", "h_CT", "PAB_SH", "AST_ALT", "BUN_CREA",
  "CO2", "GGT", "ALT", "LDH", "LDL", "AST", "BUN", "URIC", "CHOL",
  "TBIL", "TP", "PHOS", "CL", "OSM", "GLB", "TG", "A_G", "Alb",
  "DBIL", "ALP", "SCR", "EGFR", "CK", "CKMB", "GLU",
  "APROA", "APROB_APROA", "APROB",
  "Ca", "Na", "K", "Mg", "IBIL", "AG", "HDL", "HbAlc", "DD",
  "PT_sec", "TT", "INR", "ATⅢ", "APTT", "FIB", "IL_6", "CRP",
  "EF_baseline", "LVEDV_baseline", "LVESV_baseline",
  "EF_fu", "LVEDV_fu", "LVESV_fu", "ΔLVEDV", "Echo_fu_day", "Echo_fu_month",
  "Al", "As", "B", "Ba", "Cr", "Cu", "Fe", "Li", "Mn", "Mo", "Ni", "Pb",
  "Rb", "Sb", "Se", "Sn", "Sr", "V", "Zn"
)

all_vars <- c(bin_vars, ord_vars, con_vars)

## 检查缺失变量
missing_in_data <- setdiff(all_vars, names(dat))
if (length(missing_in_data) > 0) {
  warning("以下变量不在数据集中:\n",
          paste(missing_in_data, collapse = ", "))
}

## ----------------------------------------
## 3. 数据准备和中文标签
## ----------------------------------------

dat_tbl <- dat

## 将二元变量转为因子
for (v in bin_vars) {
  if (v %in% names(dat_tbl)) {
    dat_tbl[[v]] <- factor(dat_tbl[[v]])
  }
}

## 将有序变量转为有序因子
for (v in ord_vars) {
  if (v %in% names(dat_tbl)) {
    dat_tbl[[v]] <- factor(dat_tbl[[v]], ordered = TRUE)
  }
}

## 示例:将选定的二元变量重新编码为中文标签
if ("gender" %in% names(dat_tbl)) {
  dat_tbl$gender <- factor(dat_tbl$gender,
                           levels = c(0, 1),
                           labels = c("女性", "男性"))
}
if ("smoking" %in% names(dat_tbl)) {
  dat_tbl$smoking <- factor(dat_tbl$smoking,
                            levels = c(0, 1),
                            labels = c("从不/已戒烟", "目前吸烟"))
}
if ("DM" %in% names(dat_tbl)) {
  dat_tbl$DM <- factor(dat_tbl$DM,
                       levels = c(0, 1),
                       labels = c("无", "有"))
}
if ("hypertension" %in% names(dat_tbl)) {
  dat_tbl$hypertension <- factor(dat_tbl$hypertension,
                                 levels = c(0, 1),
                                 labels = c("无", "有"))
}

## 定义中文变量标签(根据需要扩展)
var_labels_zh <- c(
  resident        = "是否本地居民",
  Career          = "职业",
  gender          = "性别",
  smoking         = "吸烟状态",
  DM              = "糖尿病",
  hypertension    = "高血压",
  Cancer          = "恶性肿瘤病史",
  his_stroke      = "既往卒中史",
  AF              = "房颤",
  age             = "年龄(岁)",
  height          = "身高(cm)",
  weight          = "体重(kg)",
  BMI             = "体质指数(kg/m²)",
  SBP             = "收缩压(mmHg)",
  DBP             = "舒张压(mmHg)",
  HR              = "心率(次/分)",
  EF_baseline     = "基线LVEF(%)",
  LVEDV_baseline  = "基线LVEDV(mL)",
  LVESV_baseline  = "基线LVESV(mL)",
  EF_fu           = "随访LVEF(%)",
  LVEDV_fu        = "随访LVEDV(mL)",
  LVESV_fu        = "随访LVESV(mL)",
  ΔLVEDV          = "LVEDV变化(mL)",
  Echo_fu_day     = "随访超声间隔(天)",
  Echo_fu_month   = "随访超声间隔(月)"
  ## 根据需要添加更多变量
)

## 将标签应用于数据集(gtsummary 将使用这些标签)
for (v in names(var_labels_zh)) {
  if (v %in% names(dat_tbl)) {
    attr(dat_tbl[[v]], "label") <- var_labels_zh[[v]]
  }
}

## ----------------------------------------
## 4. 使用 gtsummary 构建基线表
## ----------------------------------------

vars_for_table <- intersect(all_vars, names(dat_tbl))

tbl_summary_gt <- dat_tbl %>%
  select(all_of(vars_for_table)) %>%
  gtsummary::tbl_summary(
    missing = "always",           # 显示缺失值
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",  # 连续变量:均值(标准差)
      all_categorical() ~ "{n} ({p}%)"     # 分类变量:n (%)
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    )
  ) %>%
  gtsummary::modify_header(label = "**变量**", stat_0 = "**整体队列 (N={N})**") %>%
  gtsummary::bold_labels()

## ----------------------------------------
## 5. 转换为 flextable 并应用三线表格式
## ----------------------------------------

ft <- gtsummary::as_flex_table(tbl_summary_gt)

## 定义三线表边框
top_bottom_border    <- officer::fp_border(color = "black", width = 2)
header_bottom_border <- officer::fp_border(color = "black", width = 1)

## 首先移除所有边框
ft <- flextable::border_remove(ft)

## 顶部边框:第一行
ft <- flextable::hline_top(ft, part = "header", border = top_bottom_border)

## 表头底部(细线)
ft <- flextable::hline_bottom(ft, part = "header", border = header_bottom_border)

## 底部边框:主体的最后一行
ft <- flextable::hline_bottom(ft, part = "body", border = top_bottom_border)

## 可选:设置中文字体(例如 SimSun/宋体)
ft <- flextable::font(ft, fontname = "SimSun", part = "all")

## 自动调整列宽
ft <- flextable::autofit(ft)

## ----------------------------------------
## 6. 直接导出到 Word
## ----------------------------------------

output_dir <- here::here("outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_docx <- file.path(output_dir, "table1_baseline_gtsummary_cn.docx")

## 创建 Word 文档
doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "表1 基线特征", style = "heading 1")
doc <- flextable::body_add_flextable(doc, value = ft)

## 保存文档
print(doc, target = output_docx)

cat("基线表 Table 1 已保存至:", output_docx, "\n")

############################################################