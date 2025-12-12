# =============================================================================
# 脚本名称: 00_Generate_Variable_Config.R
# 功能: 从原始规范文件提取变量定义,生成标准化配置文件
# 输入: variable_label_and_order.txt (原始文件)
# 输出: variable_config.R (标准配置文件)
# 作者: Statistical Consultant
# 日期: 2025-12-12
# =============================================================================

library(tibble)

# -----------------------------------------------------------------------------
# 第一部分: 定义变量标签映射
# -----------------------------------------------------------------------------

label_mapping <- tribble(
  ~variable,          ~new_label,
  # 人口学特征
  "age",              "年龄(岁)",
  "gender",           "性别",
  "height",           "身高(cm)",
  "weight",           "体重(kg)",
  "BMI",              "体质指数(kg/m²)",
  "resident",         "居住地",
  "Career",           "职业状态",
  
  # 生命体征
  "SBP",              "收缩压(mmHg)",
  "DBP",              "舒张压(mmHg)",
  "HR",               "心率(次/分)",
  
  # 危险因素
  "smoking",          "吸烟状态",
  "DM",               "糖尿病",
  "hypertension",     "高血压",
  "Cancer",           "恶性肿瘤史",
  "his_stroke",       "既往脑卒中",
  "AF",               "心房颤动",
  
  # 入院评估
  "IN_killip",        "入院Killip分级",
  "GRACE_in",         "入院GRACE评分",
  "GRACE_in_str",     "入院GRACE危险分层",
  
  # 出院评估
  "OUT_killip",       "出院Killip分级",
  "GRACE_out",        "出院GRACE评分",
  "Grace_out_str",    "出院GRACE危险分层",
  
  # 心电图特征
  "ST_dev",           "ST段抬高",
  "ST_dep",           "ST段压低",
  
  # 心肌梗死类型
  "STEMI",            "ST段抬高型心肌梗死",
  "Anterior_MI",      "前壁心梗",
  "Inferior_MI",      "下壁心梗",
  "Lateral_MI",       "侧壁心梗",
  "Posterior_MI",     "后壁心梗",
  "EA_MI",            "心尖部心梗",
  "AT_MI",            "前间壁心梗",
  "HL_MI",            "高侧壁心梗",
  "RV_MI",            "右室心梗",
  
  # 冠脉病变情况
  "Lesion_no",        "病变数量",
  "LAD",              "前降支病变",
  "LCX",              "回旋支病变",
  "RCA",              "右冠病变",
  "LM",               "左主干病变",
  "TL_LAD",           "前降支完全闭塞",
  "TL_LCX",           "回旋支完全闭塞",
  "TL_RCA",           "右冠完全闭塞",
  "TL_LM",            "左主干完全闭塞",
  "Thrombo_Burden",   "血栓负荷",
  
  # 介入治疗
  "pPCI",             "急诊PCI",
  "Stent_no",         "支架数量",
  "Temporary_pacemaker", "临时起搏器",
  
  # 并发症
  "VF",               "心室颤动",
  "Cardio_shock",     "心源性休克",
  "Cardiac_arrest_in", "院内心脏骤停",
  
  # 血常规
  "WBC",              "白细胞计数(×10⁹/L)",
  "RBC",              "红细胞计数(×10¹²/L)",
  "HGB",              "血红蛋白(g/L)",
  "PLT",              "血小板计数(×10⁹/L)",
  "HCT",              "红细胞压积(%)",
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
  "MCV",              "平均红细胞体积(fL)",
  "MCH",              "平均血红蛋白量(pg)",
  "MCHC",             "平均血红蛋白浓度(g/L)",
  "RDW",              "红细胞分布宽度(%)",
  "MPV",              "平均血小板体积(fL)",
  "PCT",              "血小板压积(%)",
  "PDW",              "血小板分布宽度(%)",
  
  # 心肌标志物
  "cTnI_baseline",    "基线肌钙蛋白I(ng/mL)",
  "cTnIpeak",         "峰值肌钙蛋白I(ng/mL)",
  "cTnI_out",         "出院肌钙蛋白I(ng/mL)",
  "NTproBNP_baseline", "基线NT-proBNP(pg/mL)",
  "NTproBNP_peak",    "峰值NT-proBNP(pg/mL)",
  "CK",               "肌酸激酶(U/L)",
  "CKMB",             "肌酸激酶同工酶(U/L)",
  
  # 肝肾功能
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
  
  # 血糖代谢
  "GLU",              "葡萄糖(mmol/L)",
  "HbAlc",            "糖化血红蛋白(%)",
  
  # 血脂
  "CHOL",             "总胆固醇(mmol/L)",
  "TG",               "甘油三酯(mmol/L)",
  "LDL",              "低密度脂蛋白(mmol/L)",
  "HDL",              "高密度脂蛋白(mmol/L)",
  "APROA",            "载脂蛋白A(g/L)",
  "APROB",            "载脂蛋白B(g/L)",
  "APROB_APROA",      "载脂蛋白B/A比值",
  
  # 炎症指标
  "CRP",              "C反应蛋白(mg/L)",
  "IL_6",             "白介素-6(pg/mL)",
  "h_CT",             "降钙素(ng/L)",
  
  # 心脏超声
  "EF_baseline",      "基线左室射血分数(%)",
  "LVEDV_baseline",   "基线左室舒张末容积(mL)",
  "LVESV_baseline",   "基线左室收缩末容积(mL)",
  "EF_fu",            "随访左室射血分数(%)",
  "LVEDV_fu",         "随访左室舒张末容积(mL)",
  "LVESV_fu",         "随访左室收缩末容积(mL)",
  "ΔLVEDV",           "左室舒张末容积变化(%)",
  "Echo_fu_day",      "超声随访时间(天)",
  "Echo_fu_month",    "超声随访时间(月)",
  
  # 甲状腺功能
  "FT3",              "游离三碘甲状腺原氨酸(pmol/L)",
  "FT4",              "游离甲状腺素(pmol/L)",
  "S_TSH",            "促甲状腺激素(mIU/L)",
  
  # 电解质
  "K",                "钾(mmol/L)",
  "Na",               "钠(mmol/L)",
  "Ca",               "钙(mmol/L)",
  "Mg",               "镁(mmol/L)",
  "PHOS",             "磷(mmol/L)",
  "CL",               "氯(mmol/L)",
  "CO2",              "二氧化碳结合力(mmol/L)",
  "OSM",              "渗透压(mOsm/L)",
  "AG",               "阴离子间隙(mmol/L)",
  
  # 凝血功能
  "DD",               "D-二聚体(mg/L)",
  "APTT",             "活化部分凝血活酶时间(秒)",
  "PT_sec",           "凝血酶原时间(秒)",
  "TT",               "凝血酶时间(秒)",
  "INR",              "国际标准化比值",
  "ATⅢ",              "抗凝血酶Ⅲ(%)",
  "FIB",              "纤维蛋白原(g/L)",
  
  # 血清金属元素(原始值)
  "Cu",               "铜(μg/L)",
  "Zn",               "锌(μg/L)",
  "Fe",               "铁(μg/L)",
  "Se",               "硒(μg/L)",
  "Pb",               "铅(μg/L)",
  "Al",               "铝(μg/L)",
  "As",               "砷(μg/L)",
  "Cr",               "铬(μg/L)",
  "Mn",               "锰(μg/L)",
  "Ni",               "镍(μg/L)",
  "Li",               "锂(μg/L)",
  "Mo",               "钼(μg/L)",
  "Rb",               "铷(μg/L)",
  "Sb",               "锑(μg/L)",
  "Sn",               "锡(μg/L)",
  "Sr",               "锶(μg/L)",
  "V",                "钒(μg/L)",
  "Ba",               "钡(μg/L)",
  "B",                "硼(μg/L)",
  
  # 血清金属元素(对数转换值)
  "log_Cu",           "铜-对数(ln μg/L)",
  "log_Zn",           "锌-对数(ln μg/L)",
  "log_Fe",           "铁-对数(ln μg/L)",
  "log_Se",           "硒-对数(ln μg/L)",
  "log_Pb",           "铅-对数(ln μg/L)",
  "log_Al",           "铝-对数(ln μg/L)",
  "log_As",           "砷-对数(ln μg/L)",
  "log_Cr",           "铬-对数(ln μg/L)",
  "log_Mn",           "锰-对数(ln μg/L)",
  "log_Ni",           "镍-对数(ln μg/L)",
  
  # 用药情况
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

# -----------------------------------------------------------------------------
# 第二部分: 定义变量分组结构
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# 第三部分: 定义分类变量列表
# -----------------------------------------------------------------------------

categorical_vars <- c(
  # 二元变量
  "gender", "smoking", "DM", "hypertension", "resident", "Career",
  "Cancer", "his_stroke", "AF", "ST_dev", "ST_dep", "STEMI",
  "Anterior_MI", "Inferior_MI", "Lateral_MI", "Posterior_MI",
  "EA_MI", "AT_MI", "HL_MI", "RV_MI", "LAD", "LCX", "RCA", "LM",
  "TL_LAD", "TL_LCX", "TL_RCA", "TL_LM", "pPCI", 
  "Temporary_pacemaker", "VF", "Cardio_shock", "Cardiac_arrest_in",
  "Aspirin", "Clopidogrel", "Ticagrelor", "tirofiban", "Statin",
  "ACEIorARB", "β_block", "CCB", "Heparin", "diuretics", 
  "insulin", "metformin",
  # 有序变量
  "GRACE_in_str", "Grace_out_str", "IN_killip", "OUT_killip"
)

# -----------------------------------------------------------------------------
# 第四部分: 保存配置文件
# -----------------------------------------------------------------------------

save(label_mapping, var_groups, categorical_vars, 
     file = "variable_config.RData")

message("✅ 配置文件已生成:  variable_config.RData")
message(sprintf("  - 变量标签数:  %d", nrow(label_mapping)))
message(sprintf("  - 分组数量: %d", length(var_groups)))
message(sprintf("  - 分类变量数: %d", length(categorical_vars)))

