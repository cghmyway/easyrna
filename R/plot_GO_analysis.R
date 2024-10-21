



# 定义函数进行GO分析和可视化
#' @title GO analysis and visulization
#'
#' @param des_res The results of differentially expressed genes from desCompare function.
#' @param pvalueCutoff Threshold for GO analysis. Default:0.05
#' @param qvalueCutoff Threshold for GO analysis. Default:0.05
#' @param topn How many GO items should be displayed?
#' @param up_color GO term color of up-regulated genes
#' @param down_color GO term color of down-regulated genes
#' @param font_size Font size of GO entries. Default:12
#' @param title_size Title size
#'
#' @return ggplot2
#' @export
#'
#' @examples run_GO_analysis(des_res = des_res, topn = 15, up_color = "red4", down_color = "green4", font_size = 12)
run_GO_analysis <- function(des_res, pvalueCutoff = 0.05, qvalueCutoff = 0.05, topn = 10,
                            up_color = 'purple', down_color = 'orange', font_size = 3, title_size = 14) {


  library(dplyr)
  library(clusterProfiler)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)

  # 筛选上调和下调基因
  deg_up <- des_res %>% dplyr::filter(log2FoldChange > 1)
  deg_down <- des_res %>% dplyr::filter(log2FoldChange < -1)

  # 上调基因的转换和GO富集分析
  deg_trance_up <- bitr(deg_up$symbol, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  deg_go_up <- enrichGO(gene = str_to_title(deg_trance_up$SYMBOL),
                        OrgDb = org.Mm.eg.db, ont = "BP",
                        keyType = "SYMBOL", pAdjustMethod = "BH",
                        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)

  # 下调基因的转换和GO富集分析
  deg_trance_down <- bitr(deg_down$symbol, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  deg_go_down <- enrichGO(gene = str_to_title(deg_trance_down$SYMBOL),
                          OrgDb = org.Mm.eg.db, ont = "BP",
                          keyType = "SYMBOL", pAdjustMethod = "BH",
                          pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)

  # 转化为数据框
  deg_go_up <- as.data.frame(deg_go_up)
  deg_go_down <- as.data.frame(deg_go_down)

  # 添加分组信息
  deg_go_up$Group <- "Up"
  deg_go_down$Group <- "Down"

  # 计算log，上调的为正数，下调的为负数
  deg_go_up$LogP <- -log10(x = deg_go_up$p.adjust)
  deg_go_down$LogP <- log10(x = deg_go_down$p.adjust)
  deg_go_all <- bind_rows(deg_go_up, deg_go_down)

  # 挑选上下调前topn条目
  deg_go_all <- deg_go_all %>% dplyr::group_by(Group) %>%
    top_n(n = topn, wt = abs(LogP))

  # 设置因子水平，调整绘图顺序
  deg_go_all$Description <- factor(deg_go_all$Description, levels = rev(deg_go_all$Description))

  # 提取上下调条目
  up.title <- deg_go_all %>% dplyr::filter(Group == "Up")
  down.title <- deg_go_all %>% dplyr::filter(Group == "Down")

  # 绘图
  ggplot(deg_go_all, aes(x = LogP, y = Description, fill = Group)) +
    geom_col() +
    theme_bw() +
    theme(legend.position = "right",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = "grey60", size = 1.2),
          axis.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = title_size)) +
    geom_text(data = up.title, aes(x = -0.2, y = Description, label = Description), size =font_size, hjust = 1) +
    geom_text(data = down.title, aes(x = 0.2, y = Description, label = Description), size =font_size, hjust = 0) +
    scale_x_continuous(breaks = seq(-15, 15, 3)) +
    labs(x = "Log(P value)", y = "", title = "Enriched GO Biological Processes") +
    scale_fill_manual(values = c(down_color,up_color)) +
    geom_text(data = up.title, aes(label = Count), hjust = -0.1, color = up_color, size = 4) +
    geom_text(data = down.title, aes(label = Count), hjust = 1.1, color = down_color, size = 4)
}































