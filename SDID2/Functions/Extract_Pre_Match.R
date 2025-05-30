# Plot Pre Period Match # 

Extract_Pre_Match <- function(treated_pre, controls_pre, 
                              optim_result){
  
  
  
  # extract weights from the optimal solution
  intercept_opt <- optim_result$solution[1]
  optimal_weights <- optim_result$solution[-1]
  optimal_weights
  
  # --- Predict synthetic control ----
  
  synthetic_linear_pred <- intercept_opt + as.vector(controls_pre %*% optimal_weights)
  synthetic_counts <- exp(synthetic_linear_pred)
  
  # --- Plot fit ----
  
  pre_match <- tibble(
    time = 1:length(treated_pre),
    synthetic_counts = synthetic_counts,
    treated_counts = treated_pre)
  
  return(pre_match)
  
}
