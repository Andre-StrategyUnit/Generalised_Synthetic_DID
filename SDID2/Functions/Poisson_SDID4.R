Poisson_SDID4 <- function(
    d, 
    n_starts = 1, 
    silent = TRUE, 
    lambda_list = c(-Inf, 
                    -5, -3,-1,0, 1,3,5)){
  

  
  # minimal wrangling
  treated_unit_id <- 1
  treated_time <- min(subset(d$time, d$treated == 1))

  
  # Raw treated outcomes (counts)
  treated_pre <- Extract_Treated(d = d, 
                                 treated_unit_id = treated_unit_id, 
                                 treated_time = treated_time)
  
  # Raw controls outcomes (counts)
  controls_pre <- Extract_Controls_Pre(d = d, 
                                       treated_unit_id = treated_unit_id, 
                                       treated_time = treated_time)
  
  # Number of control units
  n_controls <- ncol(controls_pre)
  
  # from optimiser using no regularisation
  optim_result <- Optimiser_Wrapper(
    treated_pre = treated_pre, 
    controls_pre = controls_pre, 
    lambda_list = lambda_list, 
    n_starts = n_starts, 
    silent = silent)
  
  # estimate the post intervention period avg using time weights
  cf.est <- Post_Intervention_Mu3(
    optim_result = optim_result, 
    d = d, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time, 
    controls_pre = controls_pre,
    treated_pre = treated_pre, 
    controls_post = controls_post) |> 
    as.numeric()
  
  # define the obs post intervention mean
  obs <- mean(subset(d$actual.poisson, d$unit_id == 1 & d$treated == 1))
  
  # all the results
  results <- tibble(
    obs_outcome = obs, 
    est_cf_outcome = cf.est, 
    est_att = obs - cf.est, 
    lambda_used = optim_result$lambda_used, 
    post_timepoints = length(subset(d$actual.poisson, d$unit_id == 1 & d$treated == 1))) 

}

Poisson_SDID_solution <- function(
    d, 
    n_starts = 1, 
    silent = TRUE, 
    lambda_list = c(-Inf, 
                    -5, -3,-1,0, 1,3,5)){
  
  
  
  # minimal wrangling
  treated_unit_id <- 1
  treated_time <- min(subset(d$time, d$treated == 1))
  
  
  # Raw treated outcomes (counts)
  treated_pre <- Extract_Treated(d = d, 
                                 treated_unit_id = treated_unit_id, 
                                 treated_time = treated_time)
  
  # Raw controls outcomes (counts)
  controls_pre <- Extract_Controls_Pre(d = d, 
                                       treated_unit_id = treated_unit_id, 
                                       treated_time = treated_time)
  
  # Number of control units
  n_controls <- ncol(controls_pre)
  
  # from optimiser using no regularisation
  optim_result <- Optimiser_Wrapper(
    treated_pre = treated_pre, 
    controls_pre = controls_pre, 
    lambda_list = lambda_list, 
    n_starts = n_starts, 
    silent = silent)
  
  return(optim_result$solution)
}
