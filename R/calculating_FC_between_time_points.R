calculating_FC_between_time_points <- function(count_matrix) {
  
  # Extract the time points from the column names
  time_points <- colnames(count_matrix)
  
  # Initialize an empty dataframe with the same number of rows as count_matrix
  FC_df <- data.frame(matrix(ncol = 0, nrow = nrow(count_matrix)))
  
  # Calculate fold changes between consecutive time points
  for(i in (2:ncol(count_matrix)-1)){
    
    # Calculate fold change between current column and previous column
    FC_column <- count_matrix[,i+1] - count_matrix[,i]
    
    # Add the fold-change column to FC_df with the appropriate name
    FC_df[[paste0("FC_", time_points[i], "_to_", time_points[i+1])]] <- FC_column
  }
  
  # Set row names of FC_df to match the original count_matrix
  rownames(FC_df) <- rownames(count_matrix)
  
  return(FC_df)
}