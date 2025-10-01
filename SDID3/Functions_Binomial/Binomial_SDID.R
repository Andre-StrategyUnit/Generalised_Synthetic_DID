# Binomial SDID # 

library(nloptr)

# import functions
source("Functions_Binomial/Create_Data_Binomial.R") # create dataset
source("Functions_Binomial/Pre_Wrangle_Binomial.R")
source("Functions_Binomial/Gradient_Functions_Binomial.R") # gradient functions
source("Functions_Binomial/Multi_Start_Nloptr_Binomial.R") # multi start process
source("Functions_Binomial/Optimiser_Wrapper_Binomial.R") # wrapper for extracting the optimised weights with multiple starts
source("Functions_Binomial/Time_Gradient_Functions_Binomial.R") # sdid time gradient functions
source("Functions_Binomial/Optimise_Time_Lambdas_Binomial.R") # optimise time lambda
source("Functions_Binomial/Sum_To_One_Lambda_Binomial.R") # try few lambdas to ensure time weights sum to and give performance of them
source("Functions_Binomial/Post_Intervention_Mu_Binomial.R") # calculate the post intervention mean using the time weights


  
  # import data
  d <- Create_SC_Data_Binomial(time_n = 100, 
                               units_n = 15, 
                               trials_lambda = 1000)
  
  Binomial_SDID <- function(
    d, 
    n_starts = 10,
    silent = TRUE,
    lower_lambda = -10,
    upper_lambda = 20,
    lambda_jump = 2,
    sum_to_one_tolerance = 0.01){
  
  # minimal wrangling
  treated_unit_id <- 1
  treated_time <- min(d$time[d$treated == 1])
  
  # Raw treated outcomes (counts)
  treated_pre <- Extract_Treated.bn(d)
  treated_pre_trials <- Extract_Treated.bn_trials(d)
  treated_post_trials<- Extract_Treated.bn_post_trials(d)
  
  # Raw controls outcomes (counts)
  controls_pre <- Extract_Controls_Pre.bn(d)
  controls_pre_trials <- Extract_Controls_Pre.bn_trials(d)
  
  # Number of control units
  n_controls <- ncol(controls_pre)
  
  # Extract optimal model -----
  
  optim_result <- Optimiser_Wrapper_Binomial(
    controls_pre = controls_pre, 
    controls_trials_pre = controls_pre_trials, 
    treated_pre = treated_pre, 
    treated_trials_pre = treated_pre_trials,
    wt_lb = 0, wt_ub = +Inf, 
    n_starts = n_starts, 
    silent = silent)
  
  vc_pre <- Extract_VC_Pre_Binomial(
    d = d, 
    optim_result = optim_result, 
    treated_unit_id = 1, 
    treated_time = treated_time, 
    treated_pre_trials = treated_pre_trials)
  
  vc_post <- Extract_VC_Post_Binomial(
    d = d, 
    optim_result = optim_result, 
    treated_unit_id = 1, 
    treated_time = treated_time, 
    treated_post_trials = treated_post_trials)
  
  time_lambdas <- Determine_Sum_To_One_Lambda_Binomial(
    optim_result      = optim_result,
    y_pre_successes   = vc_pre,
    y_pre_trials      = treated_pre_trials,
    y_post_successes  = vc_post,
    y_post_trials     = treated_post_trials,
    list_of_lambdas   = exp(seq(from = lower_lambda, 
                                to   = upper_lambda, 
                                by   = lambda_jump)),
    tolerance         = sum_to_one_tolerance)

  best_lambda_row <- time_lambdas %>%
      filter(between(sum_par, 1 - sum_to_one_tolerance, 1 + sum_to_one_tolerance)) %>%
      slice_min(unpenalised_objective, n = 1)
  
  # estimate the effects
  cf.est <- Post_Intervention_Mu_Binomial(
      optim_result      = optim_result, 
      d                 = d, 
      treated_unit_id   = treated_unit_id, 
      treated_time      = treated_time,
      treated_post_trials = treated_post_trials,
      treated_pre_trials = treated_pre_trials, 
      best_lambda_row   = best_lambda_row)
    
  
  obs <- mean(subset(d$actual.binomial/d$trials, d$unit_id == 1 & d$treated == 1))
  cf.actual <- subset(d$cf.prob*d$trials, d$time >= treated_time & d$unit_id == 1)
    
  results <- tibble(
    cf.est = cf.est, 
    cf.actual = sum(cf.actual)/sum(treated_post_trials), 
    obs = obs, 
    vc_post_avg = sum(vc_post)/sum(treated_post_trials), 
    binomial_sdid_att = obs - cf.est,
    binomial_sdid_att.adj = obs - (cf.est / vc_post_avg), 
    true_att = obs - cf.actual) 
  
  # results |> view()
  
  return(results)

  
  # tibble(
  #   time = 1:max(d$time), 
  #   cf = subset(d$cf.prob, d$unit_id == 1), 
  #   actual =  subset(d$actual.prob, d$unit_id == 1), 
  #   vc = c(vc_pre, vc_post), 
  #   trials = c(treated_pre_trials, treated_post_trials), 
  #   vc_prob = vc/trials) |>
  #   ggplot(mapping = aes(x = time)) + 
  #   geom_line(mapping = aes(y = cf), 
  #             col = "blue") + 
  #   geom_line(mapping = aes(y = vc_prob), 
  #             col = "red") + 
  #   geom_vline(xintercept = treated_time)
  
  }

  
  # import data
  d <- Create_SC_Data_Binomial(time_n = 100, 
                               units_n = 15, 
                               trials_lambda = 1000)
  Binomial_SDID(d = d)
  