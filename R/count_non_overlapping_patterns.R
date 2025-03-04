# Function to count non-overlapping pattern occurrences
count_non_overlapping_patterns <- function(sequence, pattern) {
  matches <- gregexpr(pattern, sequence, fixed = TRUE)[[1]]
  if (matches[1] == -1) return(0)  # No matches found
  
  # Get positions and lengths of matches
  match_positions <- matches
  match_lengths <- attr(matches, "match.length")
  
  # Remove overlapping matches
  non_overlapping <- numeric(0)
  last_end <- 0
  
  for(i in seq_along(match_positions)) {
    start_pos <- match_positions[i]
    end_pos <- start_pos + match_lengths[i] - 1
    
    if(start_pos > last_end) {
      non_overlapping <- c(non_overlapping, start_pos)
      last_end <- end_pos
    }
  }
  
  return(length(non_overlapping))
}