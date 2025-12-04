# Standard SDID Method # 

Athey_Method <- function(d){
  
  d2 <- d |>
    mutate(treated = as.logical(treated))
  
  setup <- synthdid::panel.matrices(
    panel = as.data.frame(d2), 
    unit = 1, 
    time = 2, 
    outcome = 9, 
    treatment = 5, 
    treated.last = FALSE)
  tau.hat = synthdid::synthdid_estimate(setup$Y, setup$N0, setup$T0)
  return(tau.hat[1])
  
}
