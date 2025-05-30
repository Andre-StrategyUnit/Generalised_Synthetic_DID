# Pre Wrangle Data #

Extract_Treated.bn <- function(d, 
                               treated_unit_id = treated_unit_id, 
                               treated_time = treated_time){
  
  # Raw treated outcomes (counts)
  treated_pre <- d %>%
    subset(unit_id == treated_unit_id & time < treated_time) %>%
    arrange(time) %>%
    pull(actual.binomial)
  
  return(treated_pre)
}

Extract_Treated.bn_trials <- function(d, 
                                      treated_unit_id = treated_unit_id, 
                                      treated_time = treated_time){
  
  # Raw treated outcomes (counts)
  treated_pre_trials <- d %>%
    subset(unit_id == treated_unit_id & time < treated_time) %>%
    arrange(time) %>%
    pull(trials)
  
  return(treated_pre_trials)
}


Extract_Controls_Pre.bn <- function(d, 
                                    treated_unit_id = treated_unit_id, 
                                    treated_time = treated_time){
  
  controls_pre <- d %>%
    subset(unit_id != treated_unit_id & time < treated_time) %>%
    select(unit_id, time, actual.binomial) %>%
    pivot_wider(names_from = unit_id, values_from = actual.binomial) %>%
    arrange(time) %>%
    select(-time) %>%
    as.matrix()
  
  return(controls_pre)
}

Extract_Controls_Pre.bn_trials <- function(d, 
                                           treated_unit_id = treated_unit_id, 
                                           treated_time = treated_time){
  
  controls_pre_trials <- d %>%
    subset(unit_id != treated_unit_id & time < treated_time) %>%
    select(unit_id, time, trials) %>%
    pivot_wider(names_from = unit_id, values_from = trials) %>%
    arrange(time) %>%
    select(-time) %>%
    as.matrix()
  
  return(controls_pre_trials)
}

Extract_Treated.bn_post_trials <- function(d, 
                                           treated_unit_id = treated_unit_id, 
                                           treated_time = treated_time){
  
  # Raw treated outcomes (counts)
  treated_pre_trials <- d %>%
    subset(unit_id == treated_unit_id & time >= treated_time) %>%
    arrange(time) %>%
    pull(trials)
  
  return(treated_pre_trials)
}


Extract_VC_Pre_Binomial <- function(d, optim_result, treated_unit_id, treated_time, treated_pre_trials) {
  # Extract pre-period control successes and trials
  controls_pre <- d %>%
    subset(unit_id != treated_unit_id & time < treated_time) %>%
    select(unit_id, time, actual.binomial, trials) %>%
    arrange(time) %>%
    pivot_wider(names_from = unit_id, values_from = c(actual.binomial, trials)) 
  
  # Reshape into matrices
  y_mat <- controls_pre %>%
    select(starts_with("actual.binomial_")) %>%
    as.matrix()
  
  n_mat <- controls_pre %>%
    select(starts_with("trials_")) %>%
    as.matrix()
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  logit_p <- intercept_opt + as.vector(qlogis((y_mat / n_mat) %*% weights_opt))
  p <- plogis(logit_p)
  
  # Return expected successes (treated trials * predicted p)
  expected_successes <- p * treated_pre_trials
  return(expected_successes)
}

Extract_VC_Post_Binomial <- function(d, optim_result, treated_unit_id, treated_time, treated_post_trials) {
  controls_post <- d %>%
    subset(unit_id != treated_unit_id & time >= treated_time) %>%
    select(unit_id, time, actual.binomial, trials) %>%
    arrange(time) %>%
    pivot_wider(names_from = unit_id, values_from = c(actual.binomial, trials)) 
  
  # Reshape into matrices
  y_mat <- controls_post %>%
    select(starts_with("actual.binomial_")) %>%
    as.matrix()
  
  n_mat <- controls_post %>%
    select(starts_with("trials_")) %>%
    as.matrix()
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  logit_p <- intercept_opt + as.vector(qlogis((y_mat / n_mat) %*% weights_opt))
  p <- plogis(logit_p)
  
  expected_successes <- p * treated_post_trials
  return(expected_successes)
}

Extract_VC_Pre_Binomial2 <- function(d, optim_result, treated_unit_id, treated_time, treated_pre_trials) {
  # Extract pre-period control successes and trials
  controls_pre <- d %>%
    subset(unit_id != treated_unit_id & time < treated_time) %>%
    select(unit_id, time, actual.binomial, trials) %>%
    arrange(time) %>%
    pivot_wider(names_from = unit_id, values_from = c(actual.binomial, trials))
  
  # Reshape into matrices
  y_mat <- controls_pre %>%
    select(starts_with("actual.binomial_")) %>%
    as.matrix()
  
  n_mat <- controls_pre %>%
    select(starts_with("trials_")) %>%
    as.matrix()
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  # Convert to logit scale
  p_mat <- pmin(pmax(y_mat / n_mat, 1e-8), 1 - 1e-8)
  logit_mat <- qlogis(p_mat)
  
  # Predict logits using affine model
  logit_pred <- intercept_opt + as.vector(logit_mat %*% weights_opt)
  p_pred     <- plogis(logit_pred)
  
  expected_successes <- p_pred * treated_pre_trials
  return(expected_successes)
}


Extract_VC_Post_Binomial2 <- function(d, optim_result, treated_unit_id, treated_time, treated_post_trials) {
  controls_post <- d %>%
    subset(unit_id != treated_unit_id & time >= treated_time) %>%
    select(unit_id, time, actual.binomial, trials) %>%
    arrange(time) %>%
    pivot_wider(names_from = unit_id, values_from = c(actual.binomial, trials)) 
  
  y_mat <- controls_post %>%
    select(starts_with("actual.binomial_")) %>%
    as.matrix()
  
  n_mat <- controls_post %>%
    select(starts_with("trials_")) %>%
    as.matrix()
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  p_mat <- pmin(pmax(y_mat / n_mat, 1e-8), 1 - 1e-8)
  logit_mat <- qlogis(p_mat)
  
  logit_pred <- intercept_opt + as.vector((logit_mat %*% weights_opt))
  p_pred     <- plogis(logit_pred)
  
  expected_successes <- p_pred * treated_post_trials
  return(expected_successes)
}
