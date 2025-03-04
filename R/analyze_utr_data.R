# Function to analyze UTR data
analyze_utr_data <- function(utr_data, input_genes, time_point, utr_type) {
  
  # Get the UTR column name
  utr_col <- if(utr_type == "five_prime") "five_prime_utr" else "three_prime_utr"
  
  # Define regulatory elements
  regulatory_elements <- list(
    ARE = c("ATTTA", "ATTTTA", "TTATTTATT", "TATTTAT"),
    polyA = c("AATAAA", "ATTAAA", "AGTAAA", "TATAAA"),
    TOP = c("CTTTCC", "CTTTTC", "CTTTCT"),
    IRES = c("GGATCC", "CCGCGG", "GCGCGC")
  )
  
  # Add regulatory element counts
  utr_results <- utr_data %>%
    mutate(
      ARE_count = sapply(get(utr_col), function(seq) {
        sum(sapply(regulatory_elements$ARE, function(x) count_non_overlapping_patterns(seq, x)))
      }),
      polyA_count = sapply(get(utr_col), function(seq) {
        sum(sapply(regulatory_elements$polyA, function(x) count_non_overlapping_patterns(seq, x)))
      }),
      TOP_count = sapply(get(utr_col), function(seq) {
        sum(sapply(regulatory_elements$TOP, function(x) count_non_overlapping_patterns(seq, x)))
      }),
      IRES_count = sapply(get(utr_col), function(seq) {
        sum(sapply(regulatory_elements$IRES, function(x) count_non_overlapping_patterns(seq, x)))
      }),
      # Add columns showing which specific patterns were found
      ARE_patterns = sapply(get(utr_col), function(seq) {
        found <- sapply(regulatory_elements$ARE, function(x) count_non_overlapping_patterns(seq, x) > 0)
        paste(regulatory_elements$ARE[found], collapse = ";")
      }),
      polyA_patterns = sapply(get(utr_col), function(seq) {
        found <- sapply(regulatory_elements$polyA, function(x) count_non_overlapping_patterns(seq, x) > 0)
        paste(regulatory_elements$polyA[found], collapse = ";")
      }),
      TOP_patterns = sapply(get(utr_col), function(seq) {
        found <- sapply(regulatory_elements$TOP, function(x) count_non_overlapping_patterns(seq, x) > 0)
        paste(regulatory_elements$TOP[found], collapse = ";")
      }),
      IRES_patterns = sapply(get(utr_col), function(seq) {
        found <- sapply(regulatory_elements$IRES, function(x) count_non_overlapping_patterns(seq, x) > 0)
        paste(regulatory_elements$IRES[found], collapse = ";")
      })
    )
  
  # Basic summary
  cat(sprintf("\n=== Analysis for %s %s UTRs ===\n", time_point, utr_type))
  cat(sprintf("\nInput genes: %d", length(input_genes)))
  cat(sprintf("\nGenes with UTRs: %d (%0.1f%%)", 
              length(unique(utr_results$mgi_symbol)),
              length(unique(utr_results$mgi_symbol))/length(input_genes)*100))
  cat(sprintf("\nTotal UTR sequences (isoforms): %d", nrow(utr_results)))
  cat(sprintf("\nMean isoforms per gene: %0.1f", 
              nrow(utr_results)/length(unique(utr_results$mgi_symbol))))
  
  # Regulatory element analysis
  cat("\n\nRegulatory Elements Summary:")
  elements <- c("ARE", "polyA", "TOP", "IRES")
  for(element in elements) {
    genes_with_element <- utr_results %>%
      filter(get(paste0(element, "_count")) > 0) %>%
      pull(mgi_symbol) %>%
      unique()
    
    cat(sprintf("\n\nGenes with %s: %d (%0.1f%% of genes with UTRs)", 
                element, 
                length(genes_with_element),
                length(genes_with_element)/length(unique(utr_results$mgi_symbol))*100))
    
    cat("\nTop 5 genes by count:")
    top_genes <- utr_results %>%
      group_by(mgi_symbol) %>%
      summarise(
        total_count = sum(get(paste0(element, "_count"))),
        isoform_count = n(),
        isoforms_with_element = sum(get(paste0(element, "_count")) > 0),
        patterns = paste(unique(unlist(strsplit(get(paste0(element, "_patterns")), ";"))), collapse = ";")
      ) %>%
      filter(total_count > 0) %>%
      arrange(desc(total_count)) %>%
      head(5)
    print(top_genes)
  }
  
  return(utr_results)
}