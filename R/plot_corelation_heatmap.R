





# 定义函数
#' @title Plot correlation heatmap
#'
#' @param data A data frame of the TPM or FPKM values of genes, with row names as gene names and column names as sample names.
#' @param method Correlation calculation method: "pearson", "spearman", "kendall"
#' @param cluster_rows Whether to cluster rows
#' @param cluster_cols Whether to cluster the column
#' @param sample_annotation Sample annotation data frame. The row names are sample names, and the column names are different traits or groups.
#' @param annotation_colors Sample Annotation Color. Specify a list of colors for members in different groups, used for legend colors
#' @param heatmap_colors Heatmap color scheme
#' @param show_numbers Whether to display the correlation coefficient value
#' @param main_title Heatmap title
#'
#' @return ggplot2
#' @export
#'
#' @examples plot_correlation_heatmap(data = data_tpm, method = "spearman", sample_annotation = sample_annotation, annotation_colors = list(Group=c("A"="red3","B"="blue3"),Gender=c("Male"="purple", "Female"="yellow2")))

plot_correlation_heatmap <- function(data,
                                     method = "pearson",  # 相关性计算方法：pearson, spearman, kendall
                                     cluster_rows = TRUE,  # 是否聚类行
                                     cluster_cols = TRUE,  # 是否聚类列
                                     sample_annotation = NULL,  # 样本注释数据框
                                     annotation_colors = NULL,  # 样本注释颜色
                                     heatmap_colors = colorRampPalette(c("blue", "white", "red"))(100),  # 热图颜色
                                     show_numbers = TRUE,  # 是否显示相关系数值
                                     main_title = "Sample Correlation Heatmap") {  # 热图标题

  library(pheatmap)
  library(corrplot)

  # Step 1: 将TPM数据转为矩阵
  tpm_matrix <- as.matrix(data)

  # Step 2: 计算样本间的相关性矩阵
  cor_matrix <- cor(tpm_matrix, method = method)

  # Step 3: 绘制相关性热图
  pheatmap(cor_matrix,
           cluster_rows = cluster_rows,   # 对行进行聚类
           cluster_cols = cluster_cols,   # 对列进行聚类
           display_numbers = show_numbers, # 显示相关系数值
           annotation_col = sample_annotation,  # 样本注释信息
           annotation_colors = annotation_colors,  # 样本注释的颜色
           color = heatmap_colors,  # 热图颜色方案
           main = main_title)  # 热图标题
}
