



#' @title DEGs betweent two groups
#' @description
#' Conducting differential analysis using DESeq2 with ReadCount.
#'
#' @param readcount A data frame containing ReadCount, with rows representing gene names and columns representing sample names.
#' @param group1.name Custom group names, such as "OV".
#' @param group1.number The number of samples in Group 1.
#' @param group2.name Custom group names, such as "Control".
#' @param group2.number The number of samples in Group 2.
#' @param log2FC.threshold Set the threshold for the absolute value of log2 Fold Change (|log2FC|) for differentially expressed genes, with a default value of 1.
#' @param pvalue.threshold The p-value for differentially expressed genes, with a default of 0.05.
#' @param padj.threshold The adjusted p-value for differentially expressed genes, with a default of 0.05.
#' @param top.n The number of gene names ranked at the top (with the smallest p-value and the largest log2FC).
#' @param dataset Results from the refer function.
#' @param type The type of gene names in the ReadCount file: ENSEMBL or SYMBOL.
#'
#' @return List
#' @export
#'
#' @examples res_list <- desCompare(readcount=data, group1.name="A", group1.number=3, group2.name="B", group2.number=3,dataset=dataset, type="ENSEMBL")
desCompare <- function(readcount = NULL,
                       group1.name = NULL,  #
                       group1.number = NULL,  #

                       group2.name = NULL,  #
                       group2.number = NULL,  #

                       log2FC.threshold = 1,  # log2FC
                       pvalue.threshold = 0.05,  # pvalue
                       padj.threshold = 0.05,  # padj
                       top.n = 10,  #
                       dataset = NULL,  # findname
                       type = "ENSEMBL"  # findname
){

  #
  library(tidyverse)
  library(DESeq2)
  library(clusterProfiler)
  library(ggrepel)
  library(AnnotationHub)
  library(org.Gg.eg.db)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(DOSE)

  ########## 1.readcount ###########

  ############ 2.Set data categories. ###################
  condition <- factor(c(rep(group1.name, group1.number), rep(group2.name, group2.number)),
                      levels = c(group1.name, group2.name))
  condition

  #
  colData <- data.frame(row.names = colnames(readcount), condition)

  #
  dds <- DESeqDataSetFromMatrix(readcount, colData, design = ~ condition)
  dds <- DESeq(dds)

  #
  res <- results(dds, contrast = c("condition", group1.name, group2.name))
  res <- res[order(res$pvalue), ]

  # Filter differentially expressed genes (using custom thresholds: padj, pvalue, and log2FC).
  diff_gene_deseq2 <- subset(res, padj < padj.threshold & abs(log2FoldChange) > log2FC.threshold)
  dim(diff_gene_deseq2)  #

  ####
  deg.data <- as.data.frame(res)
  deg.data$log10p <- -log10(deg.data$pvalue)
  deg.data$log10q <- -log10(deg.data$padj)
  deg.data$Group <- "Not-Sig"
  deg.data$Group[which((deg.data$pvalue < pvalue.threshold) & (deg.data$log2FoldChange > log2FC.threshold))] <- 'UP'
  deg.data$Group[which((deg.data$pvalue < pvalue.threshold) & (deg.data$log2FoldChange < -log2FC.threshold))] <- 'DOWN'
  table(deg.data$Group)  #

  deg.data$label <- ""  #
  deg.data <- deg.data[order(deg.data$pvalue), ]  #

  # Use the findname function to obtain gene symbols.
  deg.data$symbol <- findname(rownames(deg.data), dataset = dataset, type = type)

  # Extract the top n genes with the smallest p-adjusted (padj) value and the largest log2 Fold Change (log2FC) among the upregulated genes.
 # up.gene <- head(deg.data$symbol[deg.data$Group == "UP" & order(deg.data$padj, -deg.data$log2FoldChange)], top.n)
  up.gene <- head(deg.data$symbol[deg.data$Group == "UP"][order(deg.data$padj[deg.data$Group == "UP"], -deg.data$log2FoldChange[deg.data$Group == "UP"])], top.n)
  # Extract the top n genes with the smallest padj and the largest log2FC among the downregulated genes.
  down.gene <- head(deg.data$symbol[deg.data$Group == "DOWN" & order(deg.data$padj, deg.data$log2FoldChange)], top.n)

  # Fill in the gene names into the labels.
  deg.top10.gene <- c(as.character(up.gene), as.character(down.gene))
  deg.data$label[match(deg.top10.gene, deg.data$symbol)] <- deg.top10.gene

  return(list(diff_genes = diff_gene_deseq2, DEG_Dataframe = deg.data))  # Return the differential genes and the complete data table.
}







