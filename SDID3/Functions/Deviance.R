poisson_deviance <- function(y, mu) {
  # Ensure y and mu are numeric vectors of same length
  if (length(y) != length(mu)) stop("y and mu must be the same length")
  
  # Replace 0s in y with 1 for log term, then set term to 0 where y == 0
  log_term <- ifelse(y == 0, 0, y * log(y / mu))
  
  deviance <- 2 * (log_term - (y - mu))
  return(sum(deviance, na.rm = TRUE))
}