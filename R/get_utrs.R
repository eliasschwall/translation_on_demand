
get_utrs <- function(gene_df, mart.object = mart) {
  library(tidyverse)
  # fetching the biomaRt attributes we need separately for both UTRs
  utr_5 <- biomaRt::getBM(attributes = c("5utr", "mgi_symbol", "ensembl_gene_id"),
                          filters = "mgi_symbol",
                          values = gene_df$mgi_symbol,
                          mart = mart.object)
  
  utr_3 <- biomaRt::getBM(attributes = c("3utr", "mgi_symbol", "ensembl_gene_id"),
                          filters = "mgi_symbol",
                          values = gene_df$mgi_symbol,
                          mart = mart.object)
  
  # Change column names
  colnames(utr_5)[1] <- "five_prime_utr"
  colnames(utr_3)[1] <- "three_prime_utr"

  # Process 5' UTRs
  utr_5 <- utr_5 %>% 
    tibble() %>% 
    dplyr::filter(!stringr::str_detect(five_prime_utr, "Sequence unavailable")) %>%
    mutate(
      utr_5_length = nchar(five_prime_utr),
      A_count_absolute = str_count(five_prime_utr, "A"),
      A_count_percentage = str_count(five_prime_utr, "A") / utr_5_length * 100,
      T_count_absolute = str_count(five_prime_utr, "T"),
      T_count_percentage = str_count(five_prime_utr, "T") / utr_5_length * 100,
      G_count_absolute = str_count(five_prime_utr, "G"),
      G_count_percentage = str_count(five_prime_utr, "G") / utr_5_length * 100,
      C_count_absolute = str_count(five_prime_utr, "C"),
      C_count_percentage = str_count(five_prime_utr, "C") / utr_5_length * 100,
      A_to_T_ratio = str_count(five_prime_utr, "A") / str_count(five_prime_utr, "T"),
      G_to_C_ratio = str_count(five_prime_utr, "G") / str_count(five_prime_utr, "C"),
      GC_content_percentage = ((str_count(five_prime_utr, "G") + str_count(five_prime_utr, "C")) / 
                               (str_count(five_prime_utr, "A") + str_count(five_prime_utr, "T") + 
                                  str_count(five_prime_utr, "G") + str_count(five_prime_utr, "C"))*100)
    )
    
  # Process 3' UTRs  
  utr_3 <- utr_3 %>% 
    tibble() %>% 
    dplyr::filter(!stringr::str_detect(three_prime_utr, "Sequence unavailable")) %>%
    mutate(
      utr_3_length = nchar(three_prime_utr),
      A_count_absolute = str_count(three_prime_utr, "A"),
      A_count_percentage = str_count(three_prime_utr, "A") / utr_3_length * 100,
      T_count_absolute = str_count(three_prime_utr, "T"),
      T_count_percentage = str_count(three_prime_utr, "T") / utr_3_length * 100,
      G_count_absolute = str_count(three_prime_utr, "G"),
      G_count_percentage = str_count(three_prime_utr, "G") / utr_3_length * 100,
      C_count_absolute = str_count(three_prime_utr, "C"),
      C_count_percentage = str_count(three_prime_utr, "C") / utr_3_length * 100,
      A_to_T_ratio = str_count(three_prime_utr, "A") / str_count(three_prime_utr, "T"),
      G_to_C_ratio = str_count(three_prime_utr, "G") / str_count(three_prime_utr, "C"),
      GC_content_percentage = ((str_count(three_prime_utr, "G") + str_count(three_prime_utr, "C")) / 
                               (str_count(three_prime_utr, "A") + str_count(three_prime_utr, "T") + 
                                  str_count(three_prime_utr, "G") + str_count(three_prime_utr, "C"))*100)
    )           
  
  return(list(utr_5 = utr_5,
              utr_3 = utr_3))
}
