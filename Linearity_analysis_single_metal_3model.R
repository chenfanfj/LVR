############################################################
## 单金属主效应分析 - 分层回归 (Model 1 & Model 2)
## 修改内容：
##   1. LVR 定义: 基线 EF >= 50 且 随访 EF_fu < 50
##   2. Model 1: 单因素 (Metal)
##   3. Model 2: 核心调整 (Metal + Age + Gender + EF_baseline + cTnIpeak)
##   4. Model 3: 全部协变量
############################################################
## ═══════════════════════════════════════════════════════
## 步骤 1: 三模型同步分析
## ═══════════════════════════════════════════════════════

library(dplyr)
library(survey)
library(mitools)
library(here)

# 加载数据并预处理
outcome_data <- readRDS(here::here("outputs", "4 IPW2_complete", "outcome_analysis_data_with_weights_extended.rds"))
outcome_data_filtered <- outcome_data %>%
  filter(EF_baseline >= 50) %>%
  mutate(LVR = ifelse(EF_fu < 50, 1, 0))

n_imps <- length(unique(outcome_data_filtered$.imp))
metals <- c("Cu", "Zn", "Fe", "Se", "Pb", "Al", "As", "Cr", "Mn", "Ni",
            "Mo", "Rb", "Sb", "Sn", "Sr", "V", "Ba", "B", "Li")

# 协变量定义
covariates_core <- c("age", "gender", "EF_baseline", "cTnIpeak")
covariates_all <- c("age", "gender", "EF_baseline", "LVEDV_baseline", "cTnIpeak", 
                    "pPCI", "STEMI", "smoking", "DM", "hypertension", 
                    "NTproBNP_peak", "GRACE_in", "WBC", "HGB", "CRP", 
                    "CHOL", "LDL", "AST")

results_list <- list()

for(metal in metals) {
  log_metal <- paste0("log_", metal)
  cat("正在分析金属:", metal, "\n")
  
  fits_m1 <- list(); fits_m2 <- list(); fits_m3 <- list()
  
  for(imp_i in 1:n_imps) {
    dat_i <- outcome_data_filtered %>% filter(.imp == imp_i)
    
    # 计算 Model 3 所需的 PS
    ps_form <- as.formula(paste0(log_metal, " ~ ", paste(covariates_all, collapse = " + ")))
    dat_i$ps_score <- predict(lm(ps_form, data = dat_i))
    
    design_i <- svydesign(ids = ~1, weights = ~sw_trunc, data = dat_i)
    
    # 拟合三个模型
    m1 <- tryCatch({svyglm(as.formula(paste0("LVR ~ ", log_metal)), design = design_i, family = quasibinomial())}, error = function(e) NULL)
    m2 <- tryCatch({svyglm(as.formula(paste0("LVR ~ ", log_metal, " + ", paste(covariates_core, collapse = " + "))), design = design_i, family = quasibinomial())}, error = function(e) NULL)
    m3 <- tryCatch({svyglm(as.formula(paste0("LVR ~ ", log_metal, " + ps_score")), design = design_i, family = quasibinomial())}, error = function(e) NULL)
    
    if(!is.null(m1)) fits_m1[[length(fits_m1)+1]] <- m1
    if(!is.null(m2)) fits_m2[[length(fits_m2)+1]] <- m2
    if(!is.null(m3)) fits_m3[[length(fits_m3)+1]] <- m3
  }
  
  # 结果提取函数
  get_res <- function(fit_list) {
    pooled <- MIcombine(fit_list)
    s <- summary(pooled)
    beta <- s[log_metal, "results"]
    se <- s[log_metal, "se"]
    p <- 2 * (1 - pnorm(abs(beta / se)))
    return(c(sprintf("%.2f (%.2f, %.2f)", exp(beta), exp(beta - 1.96*se), exp(beta + 1.96*se)), 
             format.pval(p, digits = 3, eps = 0.001)))
  }
  
  res_m1 <- get_res(fits_m1)
  res_m2 <- get_res(fits_m2)
  res_m3 <- get_res(fits_m3)
  
  results_list[[metal]] <- data.frame(
    Metal = metal,
    Model1_OR = res_m1[1], Model1_P = res_m1[2],
    Model2_OR = res_m2[1], Model2_P = res_m2[2],
    Model3_OR = res_m3[1], Model3_P = res_m3[2]
  )
}

final_df <- bind_rows(results_list)
write.csv(final_df, here::here("outputs", "EF", "Three_Models_LVR_Comparison.csv"), row.names = FALSE)

## ═══════════════════════════════════════════════════════
## 步骤 2: 绘制中文三线表
## ═══════════════════════════════════════════════════════

if(!require(flextable)) install.packages("flextable"); library(flextable)
if(!require(officer)) library(officer)
library(officer)
# 1. 金属中文名称映射
metal_map <- c("Cu"="铜 (Cu)", "Zn"="锌 (Zn)", "Fe"="铁 (Fe)", "Se"="硒 (Se)", "Pb"="铅 (Pb)", 
               "Al"="铝 (Al)", "As"="砷 (As)", "Cr"="铬 (Cr)", "Mn"="锰 (Mn)", "Ni"="镍 (Ni)", 
               "Mo"="钼 (Mo)", "Rb"="铷 (Rb)", "Sb"="锑 (Sb)", "Sn"="锡 (Sn)", "Sr"="锶 (Sr)", 
               "V"="钒 (V)", "Ba"="钡 (Ba)", "B"="硼 (B)", "Li"="锂 (Li)")

table_df <- final_df %>%
  mutate(Metal = metal_map[Metal])

# 2. 构建三线表
ft <- flextable(table_df) %>%
  # 合并表头 (双层表头)
  add_header_row(values = c("", "Model 1 (未调整)", "Model 2 (核心调整)", "Model 3 (PS调整)"), colwidths = c(1, 2, 2, 2)) %>%
  set_header_labels(Metal = "金属", 
                    Model1_OR = "OR (95% CI)", Model1_P = "P值",
                    Model2_OR = "OR (95% CI)", Model2_P = "P值",
                    Model3_OR = "OR (95% CI)", Model3_P = "P值") %>%
  
  # 设置三线表样式
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 2)) %>% # 顶线
  hline(i = 1, part = "header", border = fp_border(width = 1)) %>% # 标题下横线
  hline_bottom(part = "all", border = fp_border(width = 2)) %>% # 底线
  
  # 对齐与字体设置
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  autofit()

# 3. 导出到 Word (可选)
doc_path <- here::here("outputs", "EF", "LVR_Three_Models_Table.docx")
read_docx() %>%
  body_add_flextable(value = ft) %>%
  print(target = doc_path)
