

# 定义函数
#' @title process_rMATS_results
#' @param event_files A list showing the Alternative Splicing Event files
#' @param pvalue_cutoff p-value threshold for differential splicing events
#' @param diff_cutoff PSI threshold for differential splicing events
#' @param orgdb GO and KEGG enrichment analysis databases. Default: org.Hs.eg.db
#' @param go_analysis Whether to conduct GO analysis
#' @param kegg_analysis Whether to conduct  KEGG analysis
#' @param save_results Whether to save the result in csv file
#' @param output_prefix Output file prefix
#'
#' @return
#' @export
#'
#' @examples event_files <- list(SE = "SE.MATS.JC.txt",RI = "RI.MATS.JC.txt",MXE = "MXE.MATS.JC.txt",A5SS = "A5SS.MATS.JC.txt",A3SS = "A3SS.MATS.JC.txt");
#' process_rMATS_results(event_files)




# 定义函数
process_rMATS_results <- function(event_files,
                                  pvalue_cutoff = 0.05,
                                  diff_cutoff = 0.1,
                                  orgdb = org.Hs.eg.db,
                                  go_analysis = TRUE,
                                  kegg_analysis = TRUE,
                                  save_results = TRUE,
                                  output_prefix = "rMATS_results") {

  # 必要的库
  library(tidyverse)
  library(paletteer)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)


  # Palette setup
  pal <- paletteer::paletteer_d("Polychrome::palette36")[3:16]

  # Helper function to process each event file
  process_event_file <- function(file, type) {
    data <- read.table(file, sep = "\t", header = TRUE)
    sig_data <- data %>%
      filter(PValue < pvalue_cutoff & abs(IncLevelDifference) > diff_cutoff) %>%
      mutate(status = if_else(IncLevelDifference > 0, "Activated", "Depressed")) %>%
      dplyr::select(geneSymbol, PValue, IncLevelDifference, status) %>%
      mutate(type = type)
    return(sig_data)
  }

  # 读取并处理所有事件类型
  sigSE <- process_event_file(event_files$SE, "SE")
  sigRI <- process_event_file(event_files$RI, "RI")
  sigMXE <- process_event_file(event_files$MXE, "MXE")
  sigA5SS <- process_event_file(event_files$A5SS, "A5SS")
  sigA3SS <- process_event_file(event_files$A3SS, "A3SS")

  # 合并所有事件
  all_events <- bind_rows(sigSE, sigRI, sigMXE, sigA5SS, sigA3SS)

  # 保存汇总结果
  if (save_results) {
    write.csv(all_events, file = paste0(output_prefix, "_sig_all_summary.csv"))
  }

  # 统计各个状态的事件数量
  event_summary <- all_events %>%
    group_by(type, status) %>%
    summarise(count = n()) %>%
    spread(status, count, fill = 0) %>%
    mutate(total = Activated + Depressed,
           prop = round(total / sum(total) * 100, 2))

  if (save_results) {
    write.csv(event_summary, file = paste0(output_prefix, "_summary_data.csv"))
  }

  # 绘制事件比例饼图
  p1 <- ggplot(event_summary, aes(x = "", y = prop, fill = type)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = pal) +
    theme_void() +
    ggrepel::geom_text_repel(aes(label = paste0(prop, "%")),
                             position = position_stack(vjust = 0.5),
                             box.padding = 0.7, size = 6) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12))

  if (save_results) {
    ggsave(paste0(output_prefix, "_Event_prop.pdf"), p1, width = 8, height = 6)
  }

  # 如果进行 GO 富集分析
  if (go_analysis) {
    all_activated <- all_events %>% filter(status == "Activated")
    all_depressed <- all_events %>% filter(status == "Depressed")

    activated_genes <- all_activated$geneSymbol
    depressed_genes <- all_depressed$geneSymbol

    activated_GO <- enrichGO(gene = activated_genes, OrgDb = orgdb, ont = "BP", keyType = "SYMBOL",
                             pAdjustMethod = "BH", pvalueCutoff = 0.2)
    depressed_GO <- enrichGO(gene = depressed_genes, OrgDb = orgdb, ont = "BP", keyType = "SYMBOL",
                             pAdjustMethod = "BH", pvalueCutoff = 0.2)

    activated_GO <- as.data.frame(activated_GO)
    depressed_GO <- as.data.frame(depressed_GO)

    activated_GO$Group <- "Activated"
    depressed_GO$Group <- "Depressed"

    all_GO <- bind_rows(activated_GO, depressed_GO)
    all_GO$LogP <- ifelse(all_GO$Group == "Activated", -log10(all_GO$p.adjust), log10(all_GO$p.adjust))

    # 只取前12个显著的GO条目
    top_GO <- all_GO %>%
      group_by(Group) %>%
      top_n(n = 12, wt = abs(LogP)) %>%
      filter(!duplicated(Description)) %>%
      arrange(desc(LogP))

    p2 <- ggplot(top_GO, aes(x = LogP, y = fct_rev(Description), fill = Group)) +
      geom_col() +
      theme_bw() +
      theme(legend.position = "right", axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), panel.grid = element_blank(),
            panel.border = element_blank(), axis.line.x = element_line(color = "grey60", size = 1.2),
            axis.text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14)) +
      geom_text(data = top_GO %>% filter(Group == "Activated"),
                aes(x = -0.2, label = Description), size = 4, hjust = 1) +
      geom_text(data = top_GO %>% filter(Group == "Depressed"),
                aes(x = 0.2, label = Description), size = 4, hjust = 0) +
      labs(x = "Log(P value)", y = "", title = "Enriched GO Biological Processes")

    if (save_results) {
      ggsave(paste0(output_prefix, "_GO_enrichment.pdf"), p2, width = 12, height = 8)
    }
  }

  # 如果进行 KEGG 富集分析
  if (kegg_analysis) {
    activated_genes <- bitr(all_activated$geneSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)$ENTREZID
    depressed_genes <- bitr(all_depressed$geneSymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)$ENTREZID

    activated_KEGG <- enrichKEGG(gene = activated_genes, organism = "hsa", pvalueCutoff = 1)
    depressed_KEGG <- enrichKEGG(gene = depressed_genes, organism = "hsa", pvalueCutoff = 1)

    activated_KEGG <- as.data.frame(activated_KEGG)
    depressed_KEGG <- as.data.frame(depressed_KEGG)

    activated_KEGG$Group <- "Activated"
    depressed_KEGG$Group <- "Depressed"

    all_KEGG <- bind_rows(activated_KEGG, depressed_KEGG)

    p3 <- ggplot(all_KEGG, aes(x = -log10(pvalue), y = Description)) +
      geom_point(aes(color = Group, size = Count)) +
      scale_color_manual(values = c('red3', 'grey3')) +
      theme_minimal() +
      geom_text_repel(aes(label = Description), size = 3) +
      labs(x = "-Log10(P value)", y = "KEGG Enrichment", title = "Enriched KEGG Pathways")

    if (save_results) {
      ggsave(paste0(output_prefix, "_KEGG_enrichment.pdf"), p3, width = 12, height = 8)
    }
  }
}

# # # 使用示例
# event_files <- list(
#   SE = "SE.MATS.JC.txt",
#   RI = "RI.MATS.JC.txt",
#   MXE = "MXE.MATS.JC.txt",
#   A5SS = "A5SS.MATS.JC.txt",
#   A3SS = "A3SS.MATS.JC.txt"
# )
#
# process_rMATS_results(event_files)

