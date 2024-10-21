#' @title plot volcano using output from desCompare
#'
#' @param deg.data desCompare's DEG_Dataframe
#' @param title volcano's title
#' @param xlab xlab's title. Default:log2FoldChange
#' @param ylab ylab's title. Default:-log10(P-value)
#' @param x_cutoff The threshold line on the X-axis.
#' @param y_cutoff The threshold line on the Y-axis.
#' @param label_column "label" column. Do not change.
#' @param group_column "Group" column. Do not change.
#' @param colors point colors. Default:colors = c("purple", "#BBBBBB", "orange")
#' @param log_scale Should points be plotted using log10p or log10q for Y-axis values? Default: log10p.
#' @param point_size Default:2
#' @param label_size Default:3
#'
#' @return volcano plot
#' @export
#'
#' @examples plot_volcano(deg.data = res[["DEG_Dataframe"]], title = "HC:AA_vs_XH", log_scale="log10p")

plot_volcano <- function(deg.data, title = "Volcano Plot", xlab = "log2FoldChange", ylab = "-log10(P-value)",
                         x_cutoff = c(-1, 1), y_cutoff = 1.3, label_column = "label", group_column = "Group",
                         log_scale = "log10p", point_size = 2, label_size = 3, colors = c("purple", "#BBBBBB", "orange")) {

  # 确认是使用log10p还是log10q
  if (!log_scale %in% c("log10p", "log10q")) {
    stop("log_scale must be either 'log10p' or 'log10q'")
  }

  ggplot(deg.data, aes(x = log2FoldChange, y = !!sym(log_scale), color = !!sym(group_column))) +
    geom_point(size = point_size) +  # 自定义点大小
    scale_color_manual(values = colors) +  # 自定义颜色
    xlab(xlab) +
    ylab(ylab) +
    geom_hline(yintercept = y_cutoff, linetype = "dashed", color = "black") +  # 添加水平线
    geom_vline(xintercept = x_cutoff, linetype = "dashed", color = "black") +  # 添加垂直线
    theme_classic() +
    geom_text_repel(aes(label = !!sym(label_column)), size = label_size, max.overlaps = Inf) +  # 自定义标签大小
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 16))
}




