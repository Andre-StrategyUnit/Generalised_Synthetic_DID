# Poisson SDID # 

# packages
library(tidyverse)
library(nloptr)

Poisson_SDID3_VC <- function(
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
  
  virtual_control <- c(vc_pre, vc_post)
  
  return(virtual_control)
}

# d <- Create_SC_Data(time_n = 100,
#                     treated_time = 81,
#                     units_n = 6,
#                     scaling_factor = 25,
#                     unit_fe_sd = 3,
#                     impact_factor = 0.5)
# Poisson_SDID3(d)




