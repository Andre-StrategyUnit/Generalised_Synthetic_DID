Target_Weights <- function(controls_pre){

  geom_means <- exp(colMeans(log(pmin(pmax(controls_pre, 1e-10), 1e10))))
  inverse_geom_means <- 1 / geom_means
  target_weights <- inverse_geom_means / sum(inverse_geom_means)
  
  return(target_weights)
  
}
