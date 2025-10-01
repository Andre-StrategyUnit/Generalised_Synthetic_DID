# Binomial SDID # 

library(nloptr)

# import functions
source("Functions_Binomial/Create_Data_Binomial.R") # create dataset
source("Functions_Binomial/Pre_Wrangle_Binomial.R")
source("Functions_Binomial/Gradient_Functions_Binomial.R") # gradient functions
source("Functions_Binomial/Multi_Start_Nloptr_Binomial.R") # multi start process
source("Functions_Binomial/Optimiser_Wrapper_Binomial.R") # wrapper for extracting the optimised weights with multiple starts
source("Functions_Binomial/Affine_Time_Weights_Binomial.R") # calculate the post intervention mean using the time weights



# -------------- top-level SDID2 --------------
Binomial_SDID2 <- function(
    d, 
    n_starts = 10,
    silent = TRUE,
    lower_lambda = -10,
    upper_lambda = 20,
    lambda_jump = 2,
    sum_to_one_tolerance = 0.01
) {
  # data wrangling as before
  treated_unit_id       <- 1
  treated_time          <- min(d$time[d$treated == 1])
  treated_pre           <- Extract_Treated.bn(d, 
                                              treated_unit_id = treated_unit_id, 
                                              treated_time = treated_time)
  treated_pre_trials    <- Extract_Treated.bn_trials(d, 
                                                     treated_unit_id = treated_unit_id, 
                                                     treated_time = treated_time)
  treated_post_trials   <- Extract_Treated.bn_post_trials(d, 
                                                          treated_unit_id = treated_unit_id, 
                                                          treated_time = treated_time)
  controls_pre          <- Extract_Controls_Pre.bn(d, 
                                                   treated_unit_id = treated_unit_id, 
                                                   treated_time = treated_time)
  controls_pre_trials   <- Extract_Controls_Pre.bn_trials(d, 
                                                          treated_unit_id = treated_unit_id, 
                                                          treated_time = treated_time)
  
  # unit-weight SC
  optim_result <- Optimiser_Wrapper_Binomial(
    controls_pre = controls_pre, 
    controls_trials_pre = controls_pre_trials, 
    treated_pre = treated_pre, 
    treated_trials_pre = treated_pre_trials,
    wt_lb = 0, wt_ub = +Inf, 
    n_starts = n_starts, 
    silent = silent
  )

  # extract VC pre/post
  vc_pre  <- Extract_VC_Pre_Binomial2(d, 
                                     optim_result, 
                                     treated_unit_id, 
                                     treated_time, 
                                     treated_pre_trials)
  vc_post <- Extract_VC_Post_Binomial2(d, 
                                      optim_result, 
                                      treated_unit_id, 
                                      treated_time, 
                                      treated_post_trials)
  
  # time-weight λ search
  time_lambdas <- Determine_Sum_To_One_Lambda_Binomial2(
    optim_result, vc_pre, treated_pre_trials, vc_post, treated_post_trials,
    exp(seq(lower_lambda, upper_lambda, by = lambda_jump)), sum_to_one_tolerance
  )
  best_lambda_row <- time_lambdas %>%
    subset(between(sum_par, 1 - sum_to_one_tolerance, 1 + sum_to_one_tolerance)) %>%
    slice_min(order_by = unpenalised_objective, 
              n = 1, 
              with_ties = FALSE)
  
  # post-intervention mean
  cf_est <- Post_Intervention_Mu_Binomial2(
    optim_result = optim_result, 
    d = d, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time,
    treated_post_trials = treated_post_trials, 
    treated_pre_trials = treated_pre_trials, 
    lambda_star = best_lambda_row$sum_to_one_lambda
  )
  
  # compile results as before
  obs <- mean(subset(d$actual.binomial/d$trials, 
                     d$unit_id == 1 & d$treated == 1))
  cf_actual <- subset(d$cf.prob*d$trials, 
                      d$time >= treated_time & d$unit_id == 1)
  results <- tibble(
    cf.est = cf_est,
    cf.actual = sum(cf_actual) / sum(treated_post_trials),
    obs = obs,
    vc_post_avg = sum(vc_post) / sum(treated_post_trials),
    binomial_sdid_att = obs - cf_est,
    true_att = obs - sum(d$cf.prob[d$unit_id==1 & d$time>=treated_time] * treated_post_trials) / sum(treated_post_trials)
  )
  return(results)
}


# import data
d <- Create_SC_Data_Binomial(time_n = 100, 
                             units_n = 15, 
                             trials_lambda = 1000)
Binomial_SDID2(d = d)


results <- tibble(
  cf.est = numeric(), 
  cf.actual = numeric(), 
  obs = numeric(), 
  vc_post_avg = numeric(), 
  binomial_sdid_att = numeric(),
  binomial_sdid_att_adj = numeric(), 
  true_att = numeric()
)

tictoc::tic()
for (i in 1:500){
  
  d <- Create_SC_Data_Binomial(time_n = 100, 
                               units_n = 15, 
                               trials_lambda = sample(x = c(100, 500, 1000), 
                                                      size = 1))
  
  results[i, ] = Binomial_SDID2(d = d)
  
}
tictoc::toc()
results |>
  view()

ggplot(data = results, 
       mapping = aes(x = true_att, y = binomial_sdid_att)) + 
  geom_point()
