############################################################
## LVR_baseline_分组版.  R – 带分组小标题和自定义顺序
############################################################

rm(list = ls(all.names= TRUE))
gc()

library(here)
library(readxl)
library(dplyr)
library(gtsummary)
library(flextable)
library(officer)

## ---------------------------------
## 1. 读取数据
## ---------------------------------
data_path <- "D:/R_try/2026_APMM/data/AMI_PCI_metal_clean3.xlsx"
dat <- readxl::read_excel(path = data_path, col_names = TRUE)

cat("数据读取成功:", nrow(dat), "行,", ncol(dat), "列\n")

## ---------------------------------
## 2. 定义变量类型(保持原样)
## ---------------------------------

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

## ---------------------------------
## 3. 数据预处理(保持原样)
## ---------------------------------

if ("ID" %in% names(dat)) {
  dat_work <- dat %>% select(-ID)
} else {
  dat_work <- dat
}

bin_vars_exist <- intersect(bin_vars, names(dat_work))
ord_vars_exist <- intersect(ord_vars, names(dat_work))
con_vars_exist <- intersect(con_vars, names(dat_work))

## 二元变量重新编码
if ("gender" %in% names(dat_work)) {
  dat_work$gender <- factor(dat_work$gender,
                            levels = c(2, 1),
                            labels = c("女性", "男性"))
}

if ("smoking" %in% names(dat_work)) {
  dat_work$smoking <- factor(dat_work$smoking,
                             levels = c(0, 1),
                             labels = c("不吸烟", "吸烟"))
}

if ("DM" %in% names(dat_work)) {
  dat_work$DM <- factor(dat_work$DM,
                        levels = c(0, 1),
                        labels = c("无", "有"))
}

if ("hypertension" %in% names(dat_work)) {
  dat_work$hypertension <- factor(dat_work$hypertension,
                                  levels = c(0, 1),
                                  labels = c("无", "有"))
}

if ("resident" %in% names(dat_work)) {
  dat_work$resident <- factor(dat_work$resident,
                              levels = c(0, 1),
                              labels = c("农村", "城镇"))
}

if ("Career" %in% names(dat_work)) {
  dat_work$Career <- factor(dat_work$Career,
                            levels = c(0, 1),
                            labels = c("无业/退休", "在职"))
}

other_bin_vars <- setdiff(bin_vars_exist, 
                          c("gender", "smoking", "DM", "hypertension", 
                            "resident", "Career"))

for (v in other_bin_vars) {
  dat_work[[v]] <- factor(dat_work[[v]],
                          levels = c(0, 1),
                          labels = c("否", "是"))
}

## 有序变量
if ("GRACE_in_str" %in% names(dat_work)) {
  dat_work$GRACE_in_str <- factor(dat_work$GRACE_in_str,
                                  levels = c(1, 2, 3),
                                  labels = c("低危", "中危", "高危"),
                                  ordered = TRUE)
}

if ("Grace_out_str" %in% names(dat_work)) {
  dat_work$Grace_out_str <- factor(dat_work$Grace_out_str,
                                   levels = c(1, 2, 3),
                                   labels = c("低危", "中危", "高危"),
                                   ordered = TRUE)
}

if ("IN_killip" %in% names(dat_work)) {
  dat_work$IN_killip <- factor(dat_work$IN_killip,
                               levels = c(1, 2, 3, 4),
                               labels = c("I级", "II级", "III级", "IV级"),
                               ordered = TRUE)
}

if ("OUT_killip" %in% names(dat_work)) {
  dat_work$OUT_killip <- factor(dat_work$OUT_killip,
                                levels = c(1, 2, 3, 4),
                                labels = c("I级", "II级", "III级", "IV级"),
                                ordered = TRUE)
}

for (v in con_vars_exist) {
  if (!  is.numeric(dat_work[[v]])) {
    dat_work[[v]] <- as.numeric(dat_work[[v]])
  }
}

cat("数据预处理完成\n")

## ---------------------------------
## 4. ★★★ 定义分组和变量顺序 ★★★
## ---------------------------------

## 定义每个分组的变量和小标题
var_groups <- list(
  list(
    title = "人口学特征",
    vars = c("age", "gender", "height", "weight", "BMI", "resident", "Career")
  ),
  list(
    title = "生命体征",
    vars = c("SBP", "DBP", "HR")
  ),
  list(
    title = "危险因素",
    vars = c("smoking", "DM", "hypertension", "Cancer", "his_stroke", "AF")
  ),
  list(
    title = "入院评估",
    vars = c("IN_killip", "GRACE_in", "GRACE_in_str")
  ),
  list(
    title = "出院评估",
    vars = c("OUT_killip", "GRACE_out", "Grace_out_str")
  ),
  list(
    title = "心电图特征",
    vars = c("ST_dev", "ST_dep")
  ),
  list(
    title = "心肌梗死类型",
    vars = c("STEMI", "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI", 
             "EA_MI", "AT_MI", "HL_MI", "RV_MI")
  ),
  list(
    title = "冠脉病变情况",
    vars = c("Lesion_no", "LAD", "LCX", "RCA", "LM", 
             "TL_LAD", "TL_LCX", "TL_RCA", "TL_LM", "Thrombo_Burden")
  ),
  list(
    title = "介入治疗",
    vars = c("pPCI", "Stent_no", "Temporary_pacemaker")
  ),
  list(
    title = "并发症",
    vars = c("VF", "Cardio_shock", "Cardiac_arrest_in")
  ),
  list(
    title = "血常规",
    vars = c("WBC", "RBC", "HGB", "PLT", "HCT",
             "NE_", "LY_", "MO_", "EO_", "BA_",
             "NE", "LY", "MO", "EO", "BA",
             "MCV", "MCH", "MCHC", "RDW", "MPV", "PCT", "PDW")
  ),
  list(
    title = "心肌标志物",
    vars = c("cTnI_baseline", "cTnIpeak", "cTnI_out", 
             "NTproBNP_baseline", "NTproBNP_peak", "CK", "CKMB")
  ),
  list(
    title = "肝肾功能",
    vars = c("ALT", "AST", "AST_ALT", "TBIL", "DBIL", "IBIL",
             "SCR", "EGFR", "BUN", "BUN_CREA", "URIC",
             "PAB_SH", "TP", "GLB", "A_G", "Alb", "GGT", "LDH", "ALP")
  ),
  list(
    title = "血糖代谢",
    vars = c("GLU", "HbAlc")
  ),
  list(
    title = "血脂",
    vars = c("CHOL", "TG", "LDL", "HDL", "APROA", "APROB", "APROB_APROA")
  ),
  list(
    title = "炎症指标",
    vars = c("CRP", "IL_6", "h_CT")
  ),
  list(
    title = "心脏超声",
    vars = c("EF_baseline", "LVEDV_baseline", "LVESV_baseline",
             "EF_fu", "LVEDV_fu", "LVESV_fu", "ΔLVEDV",
             "Echo_fu_day", "Echo_fu_month")
  ),
  list(
    title = "甲状腺功能",
    vars = c("FT3", "FT4", "S_TSH")
  ),
  list(
    title = "电解质",
    vars = c("K", "Na", "Ca", "Mg", "PHOS", "CL", "CO2", "OSM", "AG")
  ),
  list(
    title = "凝血功能",
    vars = c("DD", "APTT", "PT_sec", "TT", "INR", "ATⅢ", "FIB")
  ),
  list(
    title = "血清金属元素",
    vars = c("Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
             "Li", "Mo", "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B")
  ),
  list(
    title = "用药情况",
    vars = c("Aspirin", "Clopidogrel", "Ticagrelor", "tirofiban",
             "Statin", "ACEIorARB", "β_block", "CCB", "Heparin",
             "diuretics", "insulin", "metformin")
  )
)

## 提取所有变量(按分组顺序)
ordered_vars <- unlist(lapply(var_groups, function(g) g$vars))

## 只保留数据中存在的变量
ordered_vars_exist <- intersect(ordered_vars, names(dat_work))

cat("按分组排序后共", length(ordered_vars_exist), "个变量\n")

## ---------------------------------
## 5.  构建基线表(使用排序后的变量)
## ---------------------------------

cat("构建基线表...\n")

tbl <- dat_work %>%
  select(all_of(ordered_vars_exist)) %>%
  gtsummary::tbl_summary(
    type = list(
      all_of(intersect(con_vars_exist, ordered_vars_exist)) ~ "continuous",
      all_of(intersect(c(bin_vars_exist, ord_vars_exist), ordered_vars_exist)) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany",
    missing_text = "缺失数据",
    digits = list(
      all_continuous() ~ 2,
      all_categorical() ~ c(0, 1)
    )
  ) %>%
  gtsummary::modify_header(
    label = "**特征**",
    stat_0 = "**整体队列 (N = {N})**"
  ) %>%
  gtsummary::bold_labels()

## ---------------------------------
## 6. ★★★ 替换中文标签并插入分组小标题 ★★★
## ---------------------------------

## 标签映射(与之前相同)
label_mapping <- tibble::tribble(
  ~variable,          ~new_label,
  "age",              "年龄(岁)",
  "gender",           "性别",
  "height",           "身高(cm)",
  "weight",           "体重(kg)",
  "BMI",              "体质指数(kg/m²)",
  "resident",         "居住地",
  "Career",           "职业状态",
  "SBP",              "收缩压(mmHg)",
  "DBP",              "舒张压(mmHg)",
  "HR",               "心率(次/分)",
  "smoking",          "吸烟状态",
  "DM",               "糖尿病",
  "hypertension",     "高血压",
  "Cancer",           "恶性肿瘤史",
  "his_stroke",       "既往脑卒中",
  "AF",               "心房颤动",
  "IN_killip",        "入院Killip分级",
  "GRACE_in",         "入院GRACE评分",
  "GRACE_in_str",     "入院GRACE危险分层",
  "OUT_killip",       "出院Killip分级",
  "GRACE_out",        "出院GRACE评分",
  "Grace_out_str",    "出院GRACE危险分层",
  "ST_dev",           "ST段抬高",
  "ST_dep",           "ST段压低",
  "STEMI",            "ST段抬高型心肌梗死",
  "Anterior_MI",      "前壁心梗",
  "Inferior_MI",      "下壁心梗",
  "Lateral_MI",       "侧壁心梗",
  "Posterior_MI",     "后壁心梗",
  "EA_MI",            "心尖部心梗",
  "AT_MI",            "前间壁心梗",
  "HL_MI",            "高侧壁心梗",
  "RV_MI",            "右室心梗",
  "Lesion_no",        "病变数量",
  "TL_LM",            "左主干完全闭塞",
  "TL_LAD",           "前降支完全闭塞",
  "TL_LCX",           "回旋支完全闭塞",
  "TL_RCA",           "右冠完全闭塞",
  "LAD",              "前降支病变",
  "LCX",              "回旋支病变",
  "RCA",              "右冠病变",
  "LM",               "左主干病变",
  "Thrombo_Burden",   "血栓负荷",
  "pPCI",             "急诊PCI",
  "Stent_no",         "支架数量",
  "Temporary_pacemaker", "临时起搏器",
  "VF",               "心室颤动",
  "Cardio_shock",     "心源性休克",
  "Cardiac_arrest_in", "院内心脏骤停",
  "WBC",              "白细胞计数(×10⁹/L)",
  "RBC",              "红细胞计数(×10¹²/L)",
  "HGB",              "血红蛋白(g/L)",
  "PLT",              "血小板计数(×10⁹/L)",
  "NE_",              "中性粒细胞百分比(%)",
  "LY_",              "淋巴细胞百分比(%)",
  "MO_",              "单核细胞百分比(%)",
  "EO_",              "嗜酸性粒细胞百分比(%)",
  "BA_",              "嗜碱性粒细胞百分比(%)",
  "NE",               "中性粒细胞计数(×10⁹/L)",
  "LY",               "淋巴细胞计数(×10⁹/L)",
  "MO",               "单核细胞计数(×10⁹/L)",
  "EO",               "嗜酸性粒细胞计数(×10⁹/L)",
  "BA",               "嗜碱性粒细胞计数(×10⁹/L)",
  "HCT",              "红细胞压积(%)",
  "MCV",              "平均红细胞体积(fL)",
  "MCH",              "平均血红蛋白量(pg)",
  "MCHC",             "平均血红蛋白浓度(g/L)",
  "RDW",              "红细胞分布宽度(%)",
  "MPV",              "平均血小板体积(fL)",
  "PCT",              "血小板压积(%)",
  "PDW",              "血小板分布宽度(%)",
  "cTnI_baseline",    "基线肌钙蛋白I(ng/mL)",
  "cTnIpeak",         "峰值肌钙蛋白I(ng/mL)",
  "cTnI_out",         "出院肌钙蛋白I(ng/mL)",
  "NTproBNP_baseline", "基线NT-proBNP(pg/mL)",
  "NTproBNP_peak",    "峰值NT-proBNP(pg/mL)",
  "CK",               "肌酸激酶(U/L)",
  "CKMB",             "肌酸激酶同工酶(U/L)",
  "ALT",              "丙氨酸转氨酶(U/L)",
  "AST",              "天冬氨酸转氨酶(U/L)",
  "AST_ALT",          "AST/ALT比值",
  "TBIL",             "总胆红素(μmol/L)",
  "DBIL",             "直接胆红素(μmol/L)",
  "IBIL",             "间接胆红素(μmol/L)",
  "SCR",              "血肌酐(μmol/L)",
  "EGFR",             "肾小球滤过率(mL/min/1.73m²)",
  "BUN",              "尿素氮(mmol/L)",
  "BUN_CREA",         "尿素氮/肌酐比值",
  "URIC",             "尿酸(μmol/L)",
  "PAB_SH",           "前白蛋白(g/L)",
  "TP",               "总蛋白(g/L)",
  "GLB",              "球蛋白(g/L)",
  "A_G",              "白球比",
  "Alb",              "白蛋白(g/L)",
  "GGT",              "γ-谷氨酰转移酶(U/L)",
  "LDH",              "乳酸脱氢酶(U/L)",
  "ALP",              "碱性磷酸酶(U/L)",
  "GLU",              "葡萄糖(mmol/L)",
  "HbAlc",            "糖化血红蛋白(%)",
  "CHOL",             "总胆固醇(mmol/L)",
  "TG",               "甘油三酯(mmol/L)",
  "LDL",              "低密度脂蛋白(mmol/L)",
  "HDL",              "高密度脂蛋白(mmol/L)",
  "APROA",            "载脂蛋白A(g/L)",
  "APROB",            "载脂蛋白B(g/L)",
  "APROB_APROA",      "载脂蛋白B/A比值",
  "CRP",              "C反应蛋白(mg/L)",
  "IL_6",             "白介素-6(pg/mL)",
  "h_CT",             "降钙素(ng/L)",
  "EF_baseline",      "基线左室射血分数(%)",
  "LVEDV_baseline",   "基线左室舒张末容积(mL)",
  "LVESV_baseline",   "基线左室收缩末容积(mL)",
  "EF_fu",            "随访左室射血分数(%)",
  "LVEDV_fu",         "随访左室舒张末容积(mL)",
  "LVESV_fu",         "随访左室收缩末容积(mL)",
  "ΔLVEDV",           "左室舒张末容积变化百分比(%)",
  "Echo_fu_day",      "超声随访时间(天)",
  "Echo_fu_month",    "超声随访时间(月)",
  "FT3",              "游离三碘甲状腺原氨酸(pmol/L)",
  "FT4",              "游离甲状腺素(pmol/L)",
  "S_TSH",            "促甲状腺激素(mIU/L)",
  "K",                "钾(mmol/L)",
  "Na",               "钠(mmol/L)",
  "Ca",               "钙(mmol/L)",
  "Mg",               "镁(mmol/L)",
  "PHOS",             "磷(mmol/L)",
  "CL",               "氯(mmol/L)",
  "CO2",              "二氧化碳结合力(mmol/L)",
  "OSM",              "渗透压(mOsm/L)",
  "AG",               "阴离子间隙(mmol/L)",
  "DD",               "D-二聚体(mg/L)",
  "APTT",             "活化部分凝血活酶时间(秒)",
  "PT_sec",           "凝血酶原时间(秒)",
  "TT",               "凝血酶时间(秒)",
  "INR",              "国际标准化比值",
  "ATⅢ",              "抗凝血酶Ⅲ(%)",
  "FIB",              "纤维蛋白原(g/L)",
  "Al",               "铝(μg/L)",
  "As",               "砷(μg/L)",
  "B",                "硼(μg/L)",
  "Ba",               "钡(μg/L)",
  "Cr",               "铬(μg/L)",
  "Cu",               "铜(μg/L)",
  "Fe",               "铁(μg/L)",
  "Li",               "锂(μg/L)",
  "Mn",               "锰(μg/L)",
  "Mo",               "钼(μg/L)",
  "Ni",               "镍(μg/L)",
  "Pb",               "铅(μg/L)",
  "Rb",               "铷(μg/L)",
  "Sb",               "锑(μg/L)",
  "Se",               "硒(μg/L)",
  "Sn",               "锡(μg/L)",
  "Sr",               "锶(μg/L)",
  "V",                "钒(μg/L)",
  "Zn",               "锌(μg/L)",
  "Aspirin",          "阿司匹林",
  "Clopidogrel",      "氯吡格雷",
  "Ticagrelor",       "替格瑞洛",
  "tirofiban",        "替罗非班",
  "Statin",           "他汀类药物",
  "ACEIorARB",        "ACEI/ARB",
  "β_block",          "β受体阻滞剂",
  "CCB",              "钙通道阻滞剂",
  "Heparin",          "肝素",
  "diuretics",        "利尿剂",
  "insulin",          "胰岛素",
  "metformin",        "二甲双胍"
)

## 替换标签
tbl_body <- tbl$table_body
for (i in 1:nrow(label_mapping)) {
  var_name <- label_mapping$variable[i]
  new_label <- label_mapping$new_label[i]
  
  rows_to_update <- which(tbl_body$variable == var_name & tbl_body$row_type == "label")
  
  if (length(rows_to_update) > 0) {
    tbl_body$label[rows_to_update] <- new_label
  }
}

## ★★★ 插入分组小标题 ★★★
## 找到每个分组第一个变量的位置,插入小标题行
insert_positions <- list()
for (g in var_groups) {
  first_var <- intersect(g$vars, ordered_vars_exist)[1]
  if (! is.na(first_var)) {
    first_row <- which(tbl_body$variable == first_var & tbl_body$row_type == "label")[1]
    if (length(first_row) > 0 && ! is.na(first_row)) {
      insert_positions[[length(insert_positions) + 1]] <- list(
        position = first_row,
        title = g$title
      )
    }
  }
}

## 从后往前插入(避免位置偏移)
insert_positions <- insert_positions[order(sapply(insert_positions, function(x) x$position), decreasing = TRUE)]

for (pos_info in insert_positions) {
  header_row <- tbl_body[pos_info$position, ]
  header_row$label <- paste0("**", pos_info$title, "**")
  header_row$row_type <- "label"
  header_row$stat_0 <- ""  # 清空统计值
  
  # 插入到指定位置之前
  tbl_body <- rbind(
    tbl_body[1:(pos_info$position - 1), ],
    header_row,
    tbl_body[pos_info$position:nrow(tbl_body), ]
  )
}

## 更新表格
tbl$table_body <- tbl_body

## 添加脚注
tbl <- tbl %>%
  gtsummary::modify_footnote(
    all_stat_cols() ~ "连续变量以均数 ± 标准差表示,分类变量以例数(百分比)表示"
  )

cat("✓ 基线表构建完成\n")

## ---------------------------------
## 7. 格式化为三线表
## ---------------------------------

ft <- gtsummary::as_flex_table(tbl)

border_top_bottom <- officer::fp_border(color = "black", width = 2)
border_header <- officer::fp_border(color = "black", width = 1)

ft <- flextable::border_remove(ft)
ft <- flextable::hline_top(ft, part = "header", border = border_top_bottom)
ft <- flextable::hline_bottom(ft, part = "header", border = border_header)
ft <- flextable::hline_bottom(ft, part = "body", border = border_top_bottom)

## ★★★ 修复:使用 label 而不是 特征 ★★★
ft <- flextable::bold(ft, 
                      i = ~ grepl("^\\*\\*.*\\*\\*$", label),  ## ← 改为 label
                      j = 1)

ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
ft <- flextable::fontsize(ft, size = 10, part = "all")
ft <- flextable::autofit(ft)

## ---------------------------------
## 8. 导出Word
## ---------------------------------

output_dir <- "D:/R_try/2026_APMM/outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "Table1_Baseline_Grouped.docx")

doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "表1 研究对象基线特征", style = "heading 1")
doc <- flextable::body_add_flextable(doc, value = ft)
doc <- officer::body_add_par(doc, "")
doc <- officer::body_add_par(doc, 
                             "缩写: BMI, 体质指数; SBP, 收缩压; DBP, 舒张压; HR, 心率; DM, 糖尿病; PCI, 经皮冠状动脉介入治疗; LVEF, 左室射血分数; LVEDV, 左室舒张末容积; LVESV, 左室收缩末容积; GRACE, 全球急性冠脉事件注册评分。",
                             style = "Normal"
)

print(doc, target = output_file)

cat("✓✓✓ 成功导出分组基线表至:", output_file, "\n")

############################################################