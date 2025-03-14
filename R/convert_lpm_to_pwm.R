convert_to_pwm <- function(lpm, background_freq = rep(0.25, nrow(lpm))) {
  # Check if the input is a matrix
  if (!is.matrix(lpm)) {
    stop("Input must be a matrix.")
  }
  
  # Get the number of positions and letters
  num_letters <- nrow(lpm)
  num_positions <- ncol(lpm)
  
  # Initialize the PWM matrix
  pwm <- matrix(0, nrow = num_letters, ncol = num_positions)
  
  # Calculate the PWM using log2 transformation
  for (j in 1:num_positions) {
    for (i in 1:num_letters) {
      # Calculate PWM value
      pwm[i, j] <- log2(lpm[i, j] / background_freq[i])
    }
  }
  
  # Return the PWM matrix
  return(pwm)
}

