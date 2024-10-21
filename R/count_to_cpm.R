


#' @title Generate TPM value using Gene ReadCount data
#'
#' @param expmat A dataframe generated from 'geneLength()', with 'est_len' column
#'
#' @return Dataframe
#' @export
#'
#' @examples
countToCPM <- function(expmat) {
  # CPM
  cpm_rate <- function(counts, total_counts) {
    log_counts <- log(pmax(counts, 1e-10))  # 处理 counts 为 0 的情况
    log_total <- log(total_counts)
    log_counts + log(1e6) - log_total
  }

  # counts
  total_counts <- colSums(expmat)

  #  counts 计算
  expmat %>%
    mutate(
      across(everything(),  #
             .fns = ~ {
               colname <- cur_column()  #
               log_cpm <- cpm_rate(.x, total_counts[colname])  #  counts
               exp(log_cpm)  # CPM
             },
             .names = "{.col}_CPM")
    )
}



