tod_list_ame_wrapper <- function(todc_list,
                                 background_sequences,
                                 motif_database_path) {
  
  background_5_utr <- background_sequences$utr_5$five_prime_utr %>% 
    Biostrings::DNAStringSet() %>% 
    setNames(background_sequences$utr_5$mgi_symbol)
  
  background_3_utr <- background_sequences$utr_3$three_prime_utr %>% 
    Biostrings::DNAStringSet() %>% 
    setNames(background_sequences$utr_3$mgi_symbol)
  
  # Initialize an empty list to store the results
  ame_results <- list()
  
  ame_results <- lapply(todc_list, function(time_interval){
    
    utr_5_seq <- time_interval$utr_5$five_prime_utr %>% 
      Biostrings::DNAStringSet() %>% 
      setNames(time_interval$utr_5$mgi_symbol)
    
    utr_3_seq <- time_interval$utr_3$three_prime_utr %>% 
      Biostrings::DNAStringSet() %>% 
      setNames(time_interval$utr_3$mgi_symbol)
    
  
    # Initialize a list to store results for the current time interval
    interval_results <- list()
    
    
    # Run AME and save the output in the interval_results list
    interval_results[["utr_5"]] <- runAme(input = utr_5_seq,
                                          control = background_5_utr,
                                          database = motif_database_path,
                                          meme_path = "/meme_installation/meme-5.5.7/meme/bin/")
    
    interval_results[["utr_3"]] <- runAme(input = utr_3_seq,
                                          control = background_3_utr,
                                          database = motif_database_path,
                                          meme_path = "/meme_installation/meme-5.5.7/meme/bin/")
    
    
    return(interval_results)  # Return the results for the current time interval
  })
  
  return(ame_results)  # Return the combined results
}
