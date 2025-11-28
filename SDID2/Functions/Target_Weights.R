Target_Weights <- function(controls_pre){

  geom_means <- exp(colMeans(log(pmin(pmax(controls_pre, 1e-10), 1e10))))
  inverse_geom_means <- 1 / geom_means
  target_weights <- inverse_geom_means / sum(inverse_geom_means)
  
  return(target_weights)
  
}

Target_Intercept <- function(controls_pre, treated_pre){
  
  # matches using the DID target weights 
  target_weights <- Target_Weights(controls_pre)
  optim_result2 <- NULL
  optim_result2$solution <- c(0, target_weights)
  
  # matches using the did method, identify the control weights
  pre_match2 <- Extract_Pre_Match(treated_pre = treated_pre, 
                                  optim_result = optim_result2, 
                                  controls_pre = controls_pre) |> 
    rename(did_estimate = synthetic_counts)
  
  intercept <- sum(pre_match2$treated_counts)/sum(pre_match2$did_estimate)
  
}