# ==============================================================================
# 分析2.2: GSEA-Metabolite
# 方法:  基于代谢物与LVR的关联强度进行有序富集分析
# ==============================================================================

rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(fgsea)  # Fast Gene Set Enrichment Analysis

setwd("D:/R_try/2026_APMM/")

# ------------------------------------------------------------------------------
# 步骤1: 准备排序列表
# ------------------------------------------------------------------------------
cat("===== 步骤1: 构建代谢物排序列表 =====\n")

# 1.1 读取组间差异结果
metabo_diff <- read.csv("outputs/metabolite_group_differences_with_names.csv")

# 1.2 创建排序指标
# 使用 log2(OR) 作为排序指标(正值表示LVR组升高,负值表示降低)
metabolite_ranking <- metabo_diff %>%
  filter(! is.na(Log2FC), !is.na(Metabolite_Name)) %>%
  arrange(desc(Log2FC)) %>%
  select(Metabolite_Name, Log2FC, P_value, FDR)

# 转换为命名向量(fgsea要求格式)
ranked_metabolites <- setNames(metabolite_ranking$Log2FC, 
                               metabolite_ranking$Metabolite_Name)

cat("- 排序代谢物数:", length(ranked_metabolites), "\n")
cat("- Log2FC范围:", round(min(ranked_metabolites), 3), "~", 
    round(max(ranked_metabolites), 3), "\n\n")

# ------------------------------------------------------------------------------
# 步骤2: 定义代谢物集合(Gene Sets)
# ------------------------------------------------------------------------------
cat("===== 步骤2: 定义代谢物通路集合 =====\n")

# 使用与ORA相同的通路数据库
metabolite_sets <- list(
  
  "Oxidative_Stress" = c(
    "Glutathione", "Glutathione disulfide", "4-Hydroxynonenal",
    "Malondialdehyde", "8-iso-Prostaglandin F2alpha"
  ),
  
  "TCA_Cycle" = c(
    "Citric acid", "Isocitric acid", "Succinic acid", "Fumaric acid",
    "Malic acid", "2-Oxoglutaric acid"
  ),
  
  "Fatty_Acid_Metabolism" = c(
    "Palmitamide", "Stearamide", "Palmitic acid", "Stearic acid",
    "Acetylcarnitine", "Palmitoylcarnitine", "9-Oxooctadecanoic acid",
    "Hexadecanedioic acid", "Octadecanedioic acid"
  ),
  
  "Purine_Metabolism" = c(
    "Hypoxanthine", "Xanthine", "Uric acid", "Allopurinol"
  ),
  
  "BCAA_Metabolism" = c(
    "Valine", "Leucine", "Isoleucine", "Norleucine",
    "Ketoleucine", "2-Ketocaproic acid", "3-Methyl-2-oxovaleric acid"
  ),
  
  "Inflammatory_Lipids" = c(
    "Arachidonic acid", "Prostaglandin E2", "Leukotriene B4"
  ),
  
  "Tryptophan_Metabolism" = c(
    "Kynurenine", "Indoleacetic acid"
  ),
  
  "Nicotinate_Metabolism" = c(
    "Trigonelline", "3-Pyridylacetic acid", "4-Aminobenzoic acid"
  ),
  
  # 基于相关分析新增的集合
  "Metal_Associated_Positive" = c(
    "Trigonelline", "3-Pyridylacetic acid", "4-Aminobenzoic acid",
    "Pentaethylene glycol", "Tetraethylene glycol"
  ),
  
  "Metal_Associated_Negative" = c(
    "Stearamide", "Palmitamide", "2-Hydroxystearic acid"
  )
)

# 过滤:  仅保留在排序列表中的代谢物
metabolite_sets_filtered <- lapply(metabolite_sets, function(set) {
  intersect(set, names(ranked_metabolites))
})

# 移除空集合
metabolite_sets_filtered <- metabolite_sets_filtered[sapply(metabolite_sets_filtered, length) > 0]

cat("- 可用代谢物集合:", length(metabolite_sets_filtered), "个\n")
for(set_name in names(metabolite_sets_filtered)) {
  cat("  ", set_name, ":", length(metabolite_sets_filtered[[set_name]]), "个代谢物\n")
}

# ------------------------------------------------------------------------------
# 步骤3: 运行GSEA
# ------------------------------------------------------------------------------
cat("\n===== 步骤3: 运行GSEA =====\n")

set.seed(123)  # 可重复性

gsea_results <- fgsea(
  pathways = metabolite_sets_filtered,
  stats = ranked_metabolites,
  minSize = 2,        # 最小集合大小
  maxSize = 500,      # 最大集合大小
  nperm = 10000       # 排列次数
)

# 按FDR排序
gsea_results <- gsea_results %>%
  arrange(padj) %>%
  mutate(
    Direction = ifelse(ES > 0, "Enriched in LVR", "Enriched in Non-LVR")
  )

write.csv(gsea_results, 
          "outputs/15 pathway_enrichment/GSEA_metabolite_results.csv",
          row.names = FALSE)

cat("\n✅ GSEA完成\n")
cat("显著富集通路 (FDR < 0.05):", sum(gsea_results$padj < 0.05, na.rm = TRUE), "条\n")
cat("显著富集通路 (FDR < 0.25):", sum(gsea_results$padj < 0.25, na.rm = TRUE), "条\n\n")

# 打印结果
print(gsea_results %>% select(pathway, ES, NES, pval, padj, Direction))

# ------------------------------------------------------------------------------
# 步骤4: 可视化
# ------------------------------------------------------------------------------
cat("\n===== 步骤4: GSEA可视化 =====\n")

# 4.1 GSEA表格图
if(nrow(gsea_results %>% filter(padj < 0.25)) > 0) {
  
  plot_data_gsea <- gsea_results %>%
    filter(padj < 0.25) %>%
    arrange(NES)
  
  p_gsea_bar <- ggplot(plot_data_gsea, aes(x = NES, y = reorder(pathway, NES))) +
    geom_bar(stat = "identity", aes(fill = Direction), width = 0.7) +
    scale_fill_manual(values = c("Enriched in LVR" = "#E64B35FF", 
                                 "Enriched in Non-LVR" = "#4DBBD5FF")) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    labs(
      title = "GSEA:  Metabolite Set Enrichment in LVR",
      x = "Normalized Enrichment Score (NES)",
      y = "",
      fill = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  
  ggsave("outputs/15 pathway_enrichment/GSEA_barplot.png",
         p_gsea_bar, width = 10, height = max(6, nrow(plot_data_gsea)*0.5), dpi = 300)
  
  cat("✅ GSEA条形图已保存\n")
  
} else {
  cat("⚠️ 无显著富集(FDR<0.25),可能需要更多代谢物或调整参数\n")
}

# 4.2 GSEA经典富集图(选择最显著的通路)
if(nrow(gsea_results %>% filter(padj < 0.25)) > 0) {
  
  top_pathway <- gsea_results %>%
    filter(padj < 0.25) %>%
    arrange(padj) %>%
    slice(1) %>%
    pull(pathway)
  
  png("outputs/15 pathway_enrichment/GSEA_enrichment_plot_top.png",
      width = 10, height = 6, units = "in", res = 300)
  
  plotEnrichment(metabolite_sets_filtered[[top_pathway]], ranked_metabolites) +
    labs(title = paste("GSEA:", top_pathway),
         subtitle = paste0("NES = ", round(gsea_results$NES[gsea_results$pathway == top_pathway], 2),
                           ", FDR = ", formatC(gsea_results$padj[gsea_results$pathway == top_pathway], 
                                               format = "e", digits = 2))) +
    theme_minimal(base_size = 12)
  
  dev.off()
  
  cat("✅ GSEA富集图已保存(Top通路:", top_pathway, ")\n")
}

# ------------------------------------------------------------------------------
# 步骤5: Leading Edge分析
# ------------------------------------------------------------------------------
cat("\n===== 步骤5: Leading Edge分析 =====\n")

# 提取Leading Edge代谢物(对富集贡献最大的代谢物)
if(nrow(gsea_results %>% filter(padj < 0.05)) > 0) {
  
  leading_edge_summary <- gsea_results %>%
    filter(padj < 0.05) %>%
    select(pathway, leadingEdge) %>%
    tidyr::unnest(leadingEdge) %>%
    group_by(leadingEdge) %>%
    summarise(
      Pathways = paste(pathway, collapse = "; "),
      N_Pathways = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Pathways)) %>%
    rename(Metabolite = leadingEdge)
  
  write.csv(leading_edge_summary, 
            "outputs/15 pathway_enrichment/GSEA_leading_edge_metabolites.csv",
            row. names = FALSE)
  
  cat("✅ Leading Edge代谢物已保存\n")
  cat("   核心代谢物(出现在≥2条通路):", 
      sum(leading_edge_summary$N_Pathways >= 2), "个\n")
  
  print(head(leading_edge_summary, 10))
  
} else {
  cat("⊙ 无显著富集通路,无Leading Edge结果\n")
}

cat("\n========== 分析2.2完成 ==========\n")