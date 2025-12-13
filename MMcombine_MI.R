rm(list = ls())
gc()
# ==============================================================================
# 步骤 1: 准备工作 (同前)
# ==============================================================================
library(mice)
library(dplyr)
library(readxl)

# 读取数据 (如果已经读入可跳过)
old_mids <- readRDS("outputs/mice_imputation_final.rds") 
key_table <- read_excel("data/AMI_PCI_metal_clean3.xlsx") %>% 
  select(NO, ID) %>% mutate(ID=as.character(ID), NO=as.character(NO))
metabo_data <- read_excel("data/metabolism.xlsx", sheet = "metabo") %>% 
  mutate(NO=as.character(NO))

M_count <- old_mids$m 
new_imputed_list <- list()

# ==============================================================================
# 步骤 2: 极速循环插补 (含进度监控)
# ==============================================================================

cat("正在启动优化后的插补程序...\n")
cat("总共需处理数据集数量:", M_count, "\n")
start_time <- Sys.time() # 记录开始时间

for(i in 1:M_count) {
  
  # --- A. 数据准备 ---
  temp_data <- complete(old_mids, action = i)
  temp_data$ID <- as.character(temp_data$ID)
  
  # 合并
  temp_merged <- left_join(temp_data, key_table, by = "ID")
  temp_merged <- left_join(temp_merged, metabo_data, by = "NO")
  
  # --- B. 关键优化: 构建轻量级预测矩阵 ---
  # exclude: 排除 ID 和 NO 参与预测
  # mincor=0.1: 只有相关系数绝对值 > 0.1 的变量才会被用来做预测模型
  # 这步操作能极大减少计算量
  pred_matrix <- quickpred(temp_merged, 
                           mincor = 0.1, 
                           exclude = c("ID", "NO"))
  
  # --- C. 执行插补 ---
  # print=FALSE 保持安静
  temp_imp_obj <- mice(temp_merged, 
                       m = 1, 
                       method = 'pmm', 
                       maxit = 5, 
                       predictorMatrix = pred_matrix, # 使用优化后的矩阵
                       print = FALSE)
  
  new_imputed_list[[i]] <- complete(temp_imp_obj)
  
  # --- D. 进度报告 ---
  current_time <- Sys.time()
  time_spent <- round(difftime(current_time, start_time, units = "mins"), 2)
  avg_time_per_loop <- round(as.numeric(time_spent) / i, 2)
  est_time_left <- round((M_count - i) * avg_time_per_loop, 2)
  
  cat(sprintf("[%s] 进度: %d/%d | 已用时: %s分 | 预计剩余: %s分\n", 
              format(current_time, "%H:%M:%S"), 
              i, M_count, time_spent, est_time_left))
}

# ==============================================================================
# 步骤 3: 封装保存 (同前)
# ==============================================================================
cat("\n正在封装数据...\n")
long_data <- bind_rows(new_imputed_list, .id = ".imp")
long_data$.imp <- as.numeric(long_data$.imp)

original_data_w_na <- complete(old_mids, action = 0)
original_data_w_na$ID <- as.character(original_data_w_na$ID)
original_data_w_na <- left_join(original_data_w_na, key_table, by="ID")
original_data_w_na <- left_join(original_data_w_na, metabo_data, by="NO")
original_data_w_na$.imp <- 0
original_data_w_na$.id <- 1:nrow(original_data_w_na)

long_data <- long_data %>% group_by(.imp) %>% mutate(.id = 1:n()) %>% ungroup()
final_long_df <- bind_rows(original_data_w_na, long_data)

new_mids <- as.mids(final_long_df)
saveRDS(new_mids, "outputs/mice_imputation_with_metabo_fast.rds")

cat("\n全部完成！文件已保存。\n")