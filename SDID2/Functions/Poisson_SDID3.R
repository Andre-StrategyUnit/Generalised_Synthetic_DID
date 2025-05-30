# Poisson SDID # 

# packages
library(tidyverse)
library(nloptr)

Poisson_SDID3 <- function(
    d, 
    n_starts = 1, 
    silent = TRUE, 
    lower_lambda = -10, 
    upper_lambda = 20, 
    lambda_jump = 2, 
    sum_to_one_tolerance = 0.01){
  
  # import data
  d <- d
  
  # minimal wrangling
  treated_unit_id <- 1
  treated_time <- min(d$time[d$treated == 1])
  
  # Raw treated outcomes (counts)
  treated_pre <- Extract_Treated(d, treated_unit_id, treated_time)
  treated_post <- Extract_Treated_Post(d = d, 
                                       treated_unit_id = treated_unit_id, 
                                       treated_time = treated_time)
  
  # Raw controls outcomes (counts)
  controls_pre <- Extract_Controls_Pre(d, treated_unit_id, treated_time)
  controls_post <- Extract_Controls_Post(d, treated_unit_id, treated_time)
  
  # Number of control units
  n_controls <- ncol(controls_pre)
  
  # Extract optimal model -
  
  optim_result <- Optimiser_Wrapper(
    treated_pre = treated_pre, 
    controls_pre = controls_pre, 
    n_starts = n_starts, 
    silent = silent)
  
  # extract VC pre/post
  vc_pre  <- Extract_VC_Pre_Poisson2(d, 
                                     optim_result, 
                                     treated_unit_id, 
                                     treated_time)
  vc_post <- Extract_VC_Post_Poisson2(d, 
                                      optim_result, 
                                      treated_unit_id, 
                                      treated_time)
  
  # pre_match_measure <- Extract_Pre_Match(treated_pre = treated_pre,
  #                                        optim_result = optim_result,
  #                                        controls_pre = controls_pre)
  # ggplot(data = pre_match_measure,
  #        mapping = aes(x = time)) +
  #   geom_line(mapping = aes(y = synthetic_counts)) +
  #   geom_line(mapping = aes(y = treated_counts),
  #             col = "pink")
  
  # try a few values of lambda 
  time_lambdas <- Determine_Sum_To_One_Lambda(
    d = d, 
    optim_result = optim_result, 
    controls_pre = controls_pre, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time,
    list_of_lambdas = c(exp(seq(from = lower_lambda, 
                                to = upper_lambda, 
                                by = lambda_jump))), 
    tolerance = sum_to_one_tolerance) 
  
  # select the best one
  best_lambda_row <- time_lambdas |>
    subset(within_tolerance == TRUE) |>
    slice_min(order_by = unpenalised_objective, 
              n = 1, 
              with_ties = FALSE)
  
  # estimate the effects
  cf_est <- Post_Intervention_Mu2(
    optim_result = optim_result, 
    d = d, 
    treated_pre = treated_pre, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time,
    controls_pre = controls_pre, 
    controls_post = controls_post, 
    best_lambda_row = best_lambda_row)
  
  obs <-     (subset(d$actual.poisson, 
                     d$unit_id == 1 & d$treated == 1))
  results <- tibble(
    cf.est = cf_est,
    obs = mean(obs),
    vc_mean = mean(vc_post),
    poisson_sdid_att = obs - cf_est, 
    post_timepoints = length(vc_post))
  
  return(results)
}

# d <- Create_SC_Data(time_n = 100,
#                     treated_time = 81,
#                     units_n = 6,
#                     scaling_factor = 25,
#                     unit_fe_sd = 3,
#                     impact_factor = 0.5)
# Poisson_SDID3(d)




