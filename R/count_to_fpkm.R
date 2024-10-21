


#' @title Generate FPKM value using Gene ReadCount data
#'
#' @param expmat A dataframe generated from 'geneLength()', with 'est_len' column
#' @param effLen_col Specify the column representing gene lengths in the input data frame
#'
#' @return Dataframe
#' @export
#'
#' @examples
countToFpkm <- function(expmat, effLen_col) {
  #
  fpkm_rate <- function(counts, effLen, total_counts) {
    log_counts <- log(pmax(counts, 1e-10))  # counts = 0
    log_effLen <- log(pmax(effLen, 1e-10))  # effLen = 0
    log_counts - log_effLen - log(total_counts)
  }

  #
  effLen <- expmat[[effLen_col]]

  #  counts
  counts_cols <- names(expmat)[!names(expmat) %in% effLen_col]
  total_counts <- colSums(expmat[counts_cols])

  #
  expmat %>%
    mutate(
      across(all_of(counts_cols),  #  counts
             .fns = ~ {
               colname <- cur_column()  #
               log_fpkm <- fpkm_rate(.x, effLen, total_counts[colname])  # total counts
               exp(log_fpkm + log(1e9))  # FPKM
             },
             .names = "{.col}_FPKM")
    ) %>%
    select(-all_of(effLen_col))  #
}






