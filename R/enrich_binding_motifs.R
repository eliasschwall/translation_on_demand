enrich_binding_motifs <- function(utrs) {
  
  # needed for the mouse background motifs
  library(PWMEnrich.Mmusculus.background)
  data(PWMLogn.mm9.MotifDb.Mmus)
  
  
  # converting utrs strings into Biostrings to work with PWMEnrich
  utrs_dna <- Biostrings::DNAStringSet(utrs)
  
  # Extract the list of PWM objects and determine the maximum PWM length
  my_pwms <- PWMLogn.mm9.MotifDb.Mmus@pwms
  motif_lengths <- sapply(my_pwms, function(x) ncol(x@pfm))
  max_motif_len <- max(motif_lengths)
  
  # Filter out sequences shorter than the longest PWM
  utrs_dna_filtered <- utrs_dna[Biostrings::width(utrs_dna) >= max_motif_len]
  
  # doing the actuall motif enrichment
  res <- PWMEnrich::motifEnrichment(utrs_dna_filtered, PWMLogn.mm9.MotifDb.Mmus)
  
  # getting the report sorted for top motifs
  report <- PWMEnrich::groupReport(res, by.top.motifs=TRUE)
  
  return(report)
}






