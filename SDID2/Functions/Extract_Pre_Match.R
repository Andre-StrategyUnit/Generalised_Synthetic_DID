# Plot Pre Period Match # 

Extract_Pre_Match <- function(treated_pre, controls_pre, 
                              optim_result){
  
  
  
  # extract weights from the optimal solution
  intercept_opt <- optim_result$solution[1]
  optimal_weights <- optim_result$solution[-1]
  optimal_weights
  
  # --- Predict synthetic control ----
  
  # Clamp and transform
  log_controls <- log(pmin(pmax(controls_pre, 1e-10), 1e+10))  # T x J
  linear_pred  <- intercept_opt + as.vector(log_controls %*% optimal_weights)
  
  mu <- exp(linear_pred)
  mu <- pmin(pmax(mu, 1e-8), 1e8)  # clamp predictions
  
  # --- Plot fit ----
  
  pre_match <- tibble(
    time = 1:length(treated_pre),
    synthetic_counts = mu,
    treated_counts = treated_pre)
  
  return(pre_match)
  
}
