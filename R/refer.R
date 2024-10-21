
#' @title Make a gtf database
#' @description This function is used to construct a database from gtf file.
#'
#'
#' @param gtf A gtf for chicken from ENSEMBL
#'
#' @return a global object df
#' @export
#'
#' @examples dataset <- refer("Gallus_gallus.GRCg6a.10.gtf")
#'
refer <- function(gtf){
  library(tidyverse)
  library(data.table)
  res <- read.table(gtf, sep = "\t", header = F)
  res <- res %>% dplyr::filter(V3=="gene")

  info_gene <- res %>% dplyr::filter(str_detect(V9,"gene_name")) %>%as.data.frame()
  info_ensembl <- res %>%
    dplyr::filter(str_detect(V9,"gene_name")=="FALSE") %>%
    as.data.frame()
  # make info gene file which has gene symbol
  info_gene %>% separate(col = "V9",
                         into = c("ENSEMBL","version","gene_name","source","type"),
                         sep = ";") -> info_gene
  info_gene %>% dplyr::select("V4","V5","V7","ENSEMBL","gene_name")-> info_gene
  info_gene %>% separate(col = "ENSEMBL", into = c("no1","ENSEMBL"),sep = " ") %>% separate(
    col = "gene_name", into = c("no2","no3","SYMBOL"),sep = " ")-> info_gene
  info_gene <- info_gene %>% dplyr::select("V4","V5","V7","ENSEMBL","SYMBOL")

  # make info ensemb file without symbol
  info_ensembl %>% separate(col = "V9",
                            into = c("ENSEMBL","no1","no2","no3"),sep = ";")%>%
    dplyr::select("V4","V5","V7","ENSEMBL")-> info_ensembl
  info_ensembl %>% separate(col = "ENSEMBL", into = c("no1","ENSEMBL"),sep =" " )-> info_ensembl
  info_ensembl %>% dplyr::select("V4","V5","V7","ENSEMBL")%>% mutate(SYMBOL=ENSEMBL)-> info_ensembl

  # conbine the row
  colnames(info_gene) <- c("gene_start","gene_end","strand","ENSEMBL","SYMBOL")
  colnames(info_ensembl) <- c("gene_start","gene_end","strand","ENSEMBL","SYMBOL")
  df <-as.data.frame(bind_rows(info_gene, info_ensembl))

}




