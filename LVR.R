############################################################
## LVR_baseline_gtsummary.R – 使用 gtsummary 直接导出 Word
## 目的:
##   1) 读取 AMI_PCI_metal_clean3. xlsx (n=1124)
##   2) 使用 gtsummary 构建基线表
##   3) 应用三线表格式
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
## 3. 数据准备
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
