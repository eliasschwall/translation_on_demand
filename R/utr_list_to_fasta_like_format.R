utr_list_to_fasta_like_format <- function(utr_tibble, path_for_saving_file) {
  
  # Open a connection to the file for appending
  con <- file(paste0(path_for_saving_file, ".FASTA"), open = "w")
  
  apply(utr_tibble, MARGIN = 1, function(x) {
    utr_transcript_id <- x[4]
    utr_sequence <- x[1]
    
    # Write the gene name and sequence in FASTA format
    if(Biostrings::width(utr_sequence) >=8 ){
      writeLines(paste0(">", utr_transcript_id), con)
      writeLines(utr_sequence, con)
    }

  })
  
  # Close the file connection
  close(con)
}
