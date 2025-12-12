##八、金属之间的spearman分析
rm(list=ls())
# 安装并加载所需的包（如果未安装）
if (!require("readxl")) install.packages("readxl")
if (!require("corrplot")) install.packages("corrplot")
if (!require("psych")) install.packages("psych")
library(readxl)
library(corrplot)
library(psych)
data <- read_excel("D:/VR/EFLVR/EFLVR.xlsx")
metal_columns <- c("Al", "As", "B", "Ba", "Cr", "Cu", "Fe", "Li", "Mn", "Mo", "Ni", "Pb", "Rb", "Sb", "Se", "Sn", "Sr", "V", "Zn")
metal_data <- data[, colnames(data) %in% metal_columns]

# 计算Spearman相关系数矩阵及显著性检验结果
corr_result <- corr.test(metal_data, method = "spearman")
cor_matrix <- corr_result$r
p_matrix <- corr_result$p

output_file_path <- "D:/VR/EFLVR/金属变量相关性0.png"
png(output_file_path, width = 30 * 100, height = 30 * 100, res = 300)

# 定义颜色生成函数，使用colorRampPalette创建从红到白到蓝的颜色渐变，用于绘图
#addcol <- colorRampPalette(c("#4C0000", "white", "#040676"))

# 使用corrplot绘制相关性矩阵图，上三角用颜色表示并标记显著性，使用Spearman分析得到的p值矩阵
corrplot(cor_matrix, method = "circle", 
         #col = addcol(100),
         tl.col = "black", tl.cex = 0.8, tl.srt = 0, tl.pos = "lt",
         p.mat = p_matrix,  # 这里修改为使用Spearman分析得到的p_matrix
         diag = T, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

# 下三角用数字显示相关系数，同样使用Spearman分析得到的相关系数矩阵
corrplot(cor_matrix, method = "number", type = "lower", 
         #col = addcol(100),
         tl.col = "n", tl.cex = 0.8, tl.pos = "n", order = 'AOE',
         add = T)
dev.off()