#' get_utrs
#'
#' @param gene_list \code{list} containing the genes you want the UTRs for
#' @param biomart_attribute_of_gene_list \code{character} of which biomart attribute your gene list based on eg. ensembl_gene_id, mgi_symbol etc.
#'
#' @return \code{list} comprised of two tibbles. One for all 5' UTRs and one for all 3' UTRs + ensembl_gene_id + mgi_symbol + ensembl_transcript_id
#' @export
#'
#' @examples
#' get_utrs(gene_list, 
#' biomart_attribute_of_gene_list = "mgi_symbol",
#' mart.object = mart)
get_utrs <- function(gene_list, 
                     biomart_attribute_of_gene_list = "mgi_symbol",
                     mart.object = mart) {
  
  
  # fetching the biomaRt attributes we need seperatly for both UTRs
  utr_5 <- biomaRt::getBM(attributes = c("5utr", "mgi_symbol", "ensembl_gene_id", "ensembl_transcript_id"),
                          filters = biomart_attribute_of_gene_list,
                          values = gene_list,
                          mart = mart.object)
  
  utr_3 <- biomaRt::getBM(attributes = c("3utr", "mgi_symbol", "ensembl_gene_id", "ensembl_transcript_id"),
                          filters = biomart_attribute_of_gene_list,
                          values = gene_list,
                          mart = mart.object)
  
  # we need to change the column name because a number as the first symbols causes problems 
  colnames(utr_5)[1] <- "five_prime_utr"
  colnames(utr_3)[1] <- "three_prime_utr"

  
  utr_5 <- utr_5 %>% 
    tibble() %>% 
    dplyr::filter(!stringr::str_detect(five_prime_utr, "Sequence unavailable")) %>% # filtering out rows where no sequence was found
    mutate(
      utr_5_length = nchar(five_prime_utr), # adding the length of the UTR to the tibble
      A_count_absolute = str_count(five_prime_utr, "A"), # adding nucleotide compositions
      A_count_percentage = str_count(five_prime_utr, "A") / utr_5_length * 100,
      T_count_absolute = str_count(five_prime_utr, "T"),
      T_count_percentage = str_count(five_prime_utr, "T") / utr_5_length * 100,
      G_count_absolute = str_count(five_prime_utr, "G"),
      G_count_percentage = str_count(five_prime_utr, "G") / utr_5_length * 100,
      C_count_absolute = str_count(five_prime_utr, "C"),
      C_count_percentage = str_count(five_prime_utr, "C") / utr_5_length * 100,
      A_to_T_ratio = str_count(five_prime_utr, "A") / str_count(five_prime_utr, "T"), # adding A/T ratio
      G_to_C_ratio= str_count(five_prime_utr, "G") / str_count(five_prime_utr, "C"), # adding G/c ratio
      GC_content_percentage = ((str_count(five_prime_utr, "G") + str_count(five_prime_utr, "C")) / # adding the GC content percentage 
                                 (str_count(five_prime_utr, "A") + str_count(five_prime_utr, "T") + 
                                    str_count(five_prime_utr, "G") + str_count(five_prime_utr, "C"))*100) 
    )
    
  utr_3 <- utr_3 %>% 
    tibble() %>% 
    dplyr::filter(!stringr::str_detect(three_prime_utr, "Sequence unavailable")) %>% # filtering out rows where no sequence was found
    mutate(
      utr_3_length = nchar(three_prime_utr), # adding the length of the UTR to the tibble
      A_count_absolute = str_count(three_prime_utr, "A"), # adding nucleotide compositions
      A_count_percentage = str_count(three_prime_utr, "A") / utr_3_length * 100,
      T_count_absolute = str_count(three_prime_utr, "T"),
      T_count_percentage = str_count(three_prime_utr, "T") / utr_3_length * 100,
      G_count_absolute = str_count(three_prime_utr, "G"),
      G_count_percentage = str_count(three_prime_utr, "G") / utr_3_length * 100,
      C_count_absolute = str_count(three_prime_utr, "C"),
      C_count_percentage = str_count(three_prime_utr, "C") / utr_3_length * 100,
      A_to_T_ratio = str_count(three_prime_utr, "A") / str_count(three_prime_utr, "T"), # adding A/T ratio
      G_to_C_ratio= str_count(three_prime_utr, "G") / str_count(three_prime_utr, "C"), # adding G/c ratio
      GC_content_percentage = ((str_count(three_prime_utr, "G") + str_count(three_prime_utr, "C")) / # adding the GC content percentage 
                                 (str_count(three_prime_utr, "A") + str_count(three_prime_utr, "T") + 
                                    str_count(three_prime_utr, "G") + str_count(three_prime_utr, "C"))*100) 
    )           
  
  return(list(utr_5 = utr_5,
              utr_3 = utr_3))
}
