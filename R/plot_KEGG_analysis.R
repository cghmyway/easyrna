




# 定义函数进行KEGG分析和可视化
#' @title KEGG analysis
#'
#' @param des_res The results of differentially expressed genes from desCompare function.
#' @param pvalueCutoff Threshold for KEGG analysis. Default:0.05
#' @param qvalueCutoff Threshold for KEGG analysis. Default:0.05
#' @param up_color KEGG term color of up-regulated genes
#' @param down_color KEGG term color of up-regulated genes
#' @param font_size Font size of KEGG items
#' @param point_size Point Size
#'
#' @return ggplot2
#' @export
#'
#' @examples run_KEGG_analysis(des_res, pvalueCutoff = 0.2, qvalueCutoff = 0.2, up_color = "blue", down_color = "red", font_size =2, point_size = 8)
run_KEGG_analysis <- function(des_res, pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                              up_color = 'purple', down_color = 'orange',
                              font_size = 10, point_size = 5) {

  library(dplyr)
  library(clusterProfiler)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(ggrepel)

  # 筛选上调和下调基因
  deg_up <- des_res %>% filter(log2FoldChange > 1) %>% mutate(Group = "Up")
  deg_down <- des_res %>% filter(log2FoldChange < -1) %>% mutate(Group = "Down")

  # 基因转换为 ENTREZID
  kegg_trance_up <- bitr(deg_up$symbol, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)
  kegg_trance_down <- bitr(deg_down$symbol, fromType = "SYMBOL", toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

  # 上调和下调基因进行 KEGG 富集分析
  kegg_up <- enrichKEGG(gene = kegg_trance_up$ENTREZID, organism = "hsa",
                        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)

  kegg_down <- enrichKEGG(gene = kegg_trance_down$ENTREZID, organism = "hsa",
                          pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)

  # 转换为数据框
  kegg_up <- as.data.frame(kegg_up)
  kegg_down <- as.data.frame(kegg_down)

  # 增加 Group 信息
  kegg_up$Group <- "Up"
  kegg_down$Group <- "Down"

  # 计算 -log10(pvalue)，并合并上下调数据框
  kegg_up <- kegg_up %>% mutate(Log = -log10(.$pvalue))
  kegg_down <- kegg_down %>% mutate(Log = -log10(.$pvalue))
  kegg_all <- bind_rows(kegg_up, kegg_down)

  # 绘图
  ggplot(kegg_all, aes(x = Log, y = Description)) +
    geom_point(aes(color = Group, size = Count)) +
    scale_color_manual(values = c(down_color,up_color)) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = point_size))) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = font_size),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text.y.left = element_blank()) +
    geom_text_repel(aes(label = Description), size = font_size, show.legend = F) +
    labs(x = "-Log10(P value)", y = "KEGG Enrichment")
}

# 使用示例
# run_KEGG_analysis(des_res, pvalueCutoff = 0.01, qvalueCutoff = 0.01, up_color = "blue", down_color = "red", font_size = 12, point_size = 6)













