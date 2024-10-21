
# 定义PCA分析和绘图函数
#' @title PCA analysis of bulk RNA-seq
#'
#' @param data A data frame of the TPM or FPKM values of genes, with row names as gene names and column names as sample names.
#' @param groups Group names, in the form of a vector.
#' @param group_positions A list for the distribution of samples of each group in the data frame.
#' @param point.size point size. Default:4
#' @param label.size label size. Default:4
#' @param colors point colors
#' @param ellipse Whether to draw a confidence ellipse. Default:FALSE
#'
#' @return
#' @export
#'
#' @examples pca_plot(data = data_tpm, groups = groups, group_positions = group_positions, colors = colors, ellipse = TRUE)
pca_plot <- function(data, groups, group_positions, point.size=4,label.size=4, colors, ellipse = FALSE) {

  # 1. 对TPM数据进行 log2(TPM + 1) 转换
  log_data_tpm <- log2(data + 1)

  # 2. 提取分组数据，根据group_positions确定哪些列属于哪个分组
  selected_data <- log_data_tpm[, unlist(group_positions)]

  # 3. 转置数据，使基因在列上，样本在行上
  log_data_tpm_t <- t(selected_data)

  # 4. 进行PCA分析，并且标准化数据
  pca_result <- prcomp(log_data_tpm_t, scale. = TRUE)

  # 5. 计算每个主成分解释的变异比例
  explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pc1_var <- round(explained_variance[1] * 100, 2)
  pc2_var <- round(explained_variance[2] * 100, 2)

  # 6. 提取主成分得分
  pca_data <- as.data.frame(pca_result$x)

  # 7. 为每个样本添加分组信息
  group_labels <- unlist(lapply(1:length(group_positions), function(i) {
    rep(groups[i], length(group_positions[[i]]))
  }))

  pca_data$Group <- factor(group_labels, levels = groups)

  # 8. 绘制PCA结果
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = point.size) +
    geom_text_repel(aes(label = rownames(pca_data)), size = label.size) +  # 显示样本名称
    labs(title = paste0("PCA analysis"),
         x = paste0("PC1, ", pc1_var, "%", " variation"),
         y = paste0("PC2, ", pc2_var, "%", " variation")) +
    scale_color_manual(values = colors) +  # 手动设置组别颜色
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1)))

  # 9. 如果启用置信椭圆，则添加置信椭圆，填充颜色为该组颜色的30%透明度
  if (ellipse) {
    fill_colors <- sapply(colors, function(color) adjustcolor(color, alpha.f = 0.3))

    p <- p + stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.3, color = NA) +
      scale_fill_manual(values = fill_colors)
  }

  # 返回绘图对象
  return(p)
}

# # 示例调用
# # 假设data_tpm为输入的TPM数据，group_positions定义了不同分组的列位置
# data_tpm <- data.frame(
#   a1_TPM = c(1.3045040, 16.7206632, 0.1591585, 0.3312856, 3.4896403, 103.8641133),
#   a2_TPM = c(0.4195933, 19.0681500, 1.535799e-11, 0.9590212, 3.212506, 87.633683),
#   a3_TPM = c(0.6752132, 17.095720, 0.1588771, 8.267492e-12, 4.964944, 97.219623),
#   b1_TPM = c(0.4274642, 20.920130, 1.564608e-11, 8.141756e-12, 5.480899, 90.287460),
#   b2_TPM = c(0.2405087, 20.480700, 1.584561e-11, 0.0824559, 3.753776, 105.14445),
#   b3_TPM = c(0.5566577, 20.918630, 1.528112e-11, 0.0795184, 4.274740, 96.86149)
# )
#
# groups <- c("Group1", "Group2")
# group_positions <- list(1:3, 4:6)
# colors <- c("Group1" = "blue", "Group2" = "red")
#
# # 调用函数
# pca_plot(data = data_tpm, groups = groups, group_positions = group_positions, colors = colors, ellipse = TRUE)














