# Pre Wrangle Data #

Extract_Treated <- function(d, treated_unit_id, treated_time){
  
  # Raw treated outcomes (counts)
  treated_pre <- d %>%
    subset(unit_id == treated_unit_id & time < treated_time) %>%
    arrange(time) %>%
    pull(actual.poisson)
  
  return(treated_pre)
}

Extract_Treated_Post <- function(d, treated_unit_id, treated_time){
  
  # Raw treated outcomes (counts)
  treated_post <- d %>%
    subset(unit_id == treated_unit_id & time >= treated_time) %>%
    arrange(time) %>%
    pull(actual.poisson)
  
  return(treated_post)
}

Extract_Controls_Pre <- function(d, treated_unit_id, treated_time){
  
  controls_pre <- d %>%
    subset(unit_id != treated_unit_id & time < treated_time) %>%
    select(unit_id, time, actual.poisson) %>%
    pivot_wider(names_from = unit_id, values_from = actual.poisson) %>%
    arrange(time) %>%
    select(-time) %>%
    as.matrix()
  
  return(controls_pre)
}

Extract_Controls_Post <- function(d, treated_unit_id, treated_time){
  
  controls_pre <- d %>%
    subset(unit_id != treated_unit_id & time >= treated_time) %>%
    select(unit_id, time, actual.poisson) %>%
    pivot_wider(names_from = unit_id, values_from = actual.poisson) %>%
    arrange(time) %>%
    select(-time) %>%
    as.matrix()
  
  return(controls_pre)
}
Extract_VC_Pre_Poisson2 <- function(d, optim_result, treated_unit_id, treated_time) {
  
  
  controls_pre <- Extract_Controls_Pre(d, treated_unit_id, treated_time)
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  # clamp then log-transform
  y_mat        <- pmin(pmax(controls_pre, 1e-10), 1e+10)
  log_controls <- log(y_mat)
  
  # affine on log scale → log μ
  log_pred         <- intercept_opt + log_controls %*% weights_opt
  expected_success <- exp(log_pred)
  
  as.vector(expected_success)
}

Extract_VC_Post_Poisson2 <- function(d, optim_result, treated_unit_id, treated_time) {
  
  
  controls_post <- Extract_Controls_Post(d, treated_unit_id, treated_time)
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  y_mat        <- pmin(pmax(controls_post, 1e-10), 1e+10)
  log_controls <- log(y_mat)
  
  log_pred         <- intercept_opt + log_controls %*% weights_opt
  expected_success <- exp(log_pred)
  
  as.vector(expected_success)
}

