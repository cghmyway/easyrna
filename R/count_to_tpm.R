


#' @title Generate TPM value using Gene ReadCount data
#'
#' @param expmat A dataframe generated from 'geneLength()', with 'est_len' column
#' @param effLen_col Specify the column representing gene lengths in the input data frame
#'
#' @return Dataframe
#' @export
#'
#' @examples
#'
countToTpm <- function(expmat, effLen_col) {

  rate <- function(counts, effLen) {
    log_counts <- log(pmax(counts, 1e-10))
    log_effLen <- log(pmax(effLen, 1e-10))
    log_counts - log_effLen
  }


  effLen <- expmat[[effLen_col]]

  #
  expmat %>%
    mutate(
      across(!all_of(effLen_col),  #
             .fns = ~ {
               log_rate <- rate(.x, effLen)
               denom <- log(sum(exp(log_rate)))  #
               exp(log_rate - denom + log(1e6))  #
             },
             .names = "{.col}_TPM")
    ) %>%
    dplyr::select(-all_of(effLen_col))  #
}

