



#' @title Generate TPM value using FPKM
#'
#' @param expmat
#'
#' @return Dataframe
#' @export
#'
#' @examples
fpkmToTpm <- function(expmat) {
  #
  tpm_rate <- function(fpkm) {
    log_fpkm <- log(pmax(fpkm, 1e-10))
    log_total_fpkm <- log(sum(exp(log_fpkm)))
    log_fpkm - log_total_fpkm + log(1e6)
  }

  #
  expmat %>%
    mutate(
      across(everything(),
             .fns = ~ {
               log_tpm <- tpm_rate(.x)
               exp(log_tpm)
             },
             .names = "{.col}_TPM")
    )
}
















