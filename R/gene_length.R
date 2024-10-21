


#' @title A GeneLength Used for Counting
#' @description Before transform the count into other, a gene length should be made.
#'
#'
#' @param gtf Input gtf.file for your species
#' @param df Here is the ReadCount data for each sample, with gene names as row names and sample names as column names.
#'
#' @return Return a data frame containing the corresponding gene lengths and the ReadCount for each sample.
#' @export
#' @import tidyverse
#' @examples length <- geneLength(gtf)
geneLength <- function(df,gtf){
  library(tidyverse)
  gtf <- read_tsv(gtf,comment = "#",
                  col_names = c("chr","source","type","start","end",
                                "score","strand","phase","attributes"))%>%
 filter(type=="exon")%>%mutate(len=end-start+1)%>% dplyr::select(start,end,attributes,len)
gtf$attributes %>% str_extract(., "gene_id \"[\\w|\\.]+") %>% str_remove(., "gene_id \"") -> gtf$gene_id
gtf <- gtf %>% dplyr::select(start, end, gene_id, len) %>%
  distinct(start,end,gene_id, .keep_all = T) %>%
  dplyr::select(gene_id,len) %>% group_by(gene_id) %>%
  summarise(est_len=sum(len)) %>% column_to_rownames(var = "gene_id")

expmat <- df %>% bind_cols(gtf) %>% drop_na()
}







