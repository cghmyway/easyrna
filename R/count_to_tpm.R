


#' @title Generate TPM value using Gene ReadCount data
#'
#' @param data A dataframe generated from 'geneLength()', with 'est_len' column
#' @param effLen_col Specify the column representing gene lengths in the input data frame
#'
#' @return Dataframe
#' @export
#'
#' @examples
#'
countToTpm <- function(data, effLen_col = "est_len") {
  # 内部的 countToTpm 函数
  countToTpm <- function(counts, effLen) {
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }

  # 确保指定的基因长度列存在于数据框中
  if (!(effLen_col %in% names(data))) {
    stop(paste("Column", effLen_col, "not found in the data"))
  }

  # 获取所有 read count 列（排除基因长度列）
  read_count_cols <- setdiff(names(data), effLen_col)

  # 计算 TPM
  result <- data %>%
    mutate(across(all_of(read_count_cols),
                  ~ countToTpm(., .data[[effLen_col]]),
                  .names = "{.col}_TPM"))

  return(result)
}

