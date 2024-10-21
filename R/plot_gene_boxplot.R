

# 定义函数
#' @title Plot gene expression with boxplot
#'
#' @param data A data frame of the TPM or FPKM values of genes, with row names as gene names and column names as sample names.
#' @param groups Group names, in the form of a vector.
#' @param group_positions A list for the distribution of samples of each group in the data frame.
#' @param gene_of_interest Gene of interest.
#' @param method Methods for Differential Analysis Between Groups.
#' @param colors The color of each group.
#'
#' @return ggplot2
#' @export
#'
#' @examples plot_gene_boxplot(data = data_tpm, groups = c("A","B"), group_positions = list(1:3,4:6), gene_of_interest = c("ENSGALG00000012172"), method = "t.test")
plot_gene_boxplot <- function(data, groups, group_positions, gene_of_interest, method = "t.test", colors = NULL) {

  library(ggplot2)
  library(reshape2)  # 可能需要用到 melt 函数
  library(ggpubr)    # 用于添加显著性标记
  library(RColorBrewer)  # 用于调色板

  # 提取该基因在各组的数据
  gene_data_list <- lapply(seq_along(groups), function(i) {
    data[gene_of_interest, group_positions[[i]], drop = FALSE]
  })

  # 为数据添加组别标签
  melted_data_list <- lapply(seq_along(gene_data_list), function(i) {
    melted_data <- melt(gene_data_list[[i]], varnames = "Sample", value.name = "TPM", id.vars = NULL)
    melted_data$Group <- groups[i]
    melted_data$Gene <- gene_of_interest  # 添加基因列
    return(melted_data)
  })

  # 合并所有组的数据
  combined_gene_data <- do.call(rbind, melted_data_list)

  # 添加 log10(TPM + 1) 列
  combined_gene_data$log_TPM <- log10(combined_gene_data$TPM + 1)

  # 如果未提供颜色，则使用默认调色板
  if (is.null(colors)) {
    colors <- brewer.pal(min(length(groups), 3), "Set1")  # 默认颜色
  }

  # 绘制箱线图，使用 log10(TPM + 1)
  p <- ggplot(combined_gene_data, aes(x = Group, y = log_TPM, fill = Group)) +
    geom_boxplot(color = "black") +  # 只设置箱线轮廓颜色
    theme_dose() +
    labs(x = NULL, y = "Gene expression \n log10(TPM + 1)", title = "") +
    scale_fill_manual(values = colors) +  # 使用用户提供的颜色
    geom_jitter(aes(fill = Group), shape = 21, position = position_jitter(0.2), size = 3, color = "black",
                alpha = 0.5) +  # 直接设置透明度
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ Gene)  # 按基因分面显示

  # 根据选择的检验方法进行显著性差异标记
  if (method == "t.test") {
    p <- p + stat_compare_means(method = "t.test", label = "p.signif",
                                comparisons = combn(groups, 2, simplify = FALSE),
                                aes(group = Group))
  } else if (method == "anova") {
    p <- p + stat_compare_means(method = "anova", label = "p.signif", aes(group = Group))
  } else {
    stop("Unsupported method. Please use 't.test' or 'anova'.")
  }

  return(p)
}




# 使用函数示例
# 假设 data_tpm 是你的原始 TPM 数据框
# groups <- c("A", "B")  # 分组名称
# group_positions <- list(1:3, 4:6)  # 每个组的数据列索引
# gene_of_interest <- "ENSGALG00000000055"
# plot_gene_boxplot(data_tpm, groups, group_positions, gene_of_interest, method = "t.test")




