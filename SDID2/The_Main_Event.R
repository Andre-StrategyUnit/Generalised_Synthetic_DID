
# --- Basic setup ----

options(scipen = 999999999)

# background packages
library(nloptr)
library(tidyverse)
library(ggplot2)

# import functions
source("Functions/Create_Data.R") # create dataset
source("Functions/Pre_Wrangle.R")
source("Functions/Gradient_Functions.R") # gradient functions
source("Functions/Time_Gradient_Functions.R") # sdid time gradient functions
source("Functions/Optimise_Unit_Lambdas.R") # optimise unit lambdas
source("Functions/Optimise_Time_Lambdas.R") # optimise time lambda
source("Functions/Multi_Start_Nloptr.R") # multi start process
source("Functions/Optimiser_Wrapper.R") # wrapper for extracting the optimised weights with multiple starts
source("Functions/Extract_Pre_Match.R") # extract data to assess pre period match
source("Functions/Sum_To_One_Lambda.R") # try few lambdas to ensure time weights sum to and give performance of them
source("Functions/Post_Intervention_Mu.R") # calculate the post intervention mean using the time weights
source("Functions/Athey_Method.R")
source("Functions/Target_Weights.R")
source("Functions/Improved_Placebo.R")
source("Functions/Poisson_SDID3.R")
source("Functions/Poisson_SDID4.R")

# --- Create dataset ----

good_ctrls.v <- 15
scale_fctr.v <- 1
impact_fctr.v <- 1

d <- Create_SC_Data(time_n = 100,
                    treated_time = 80, 
                    treated_units = 2,
                    units_n = 100, 
                    good_ctrls = good_ctrls.v, 
                    scaling_factor = scale_fctr.v, 
                    unit_fe_sd = 6, 
                    impact_factor = impact_fctr.v)

treated_unit_id <- 1
treated_time <- min(d$time[d$treated == 1])

# --- Extract pre-period data ----

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

# --- Extract solutions from different model -----

# from optimiser using no regularisation
optim_result <- Optimiser_Wrapper(
  treated_pre = treated_pre, 
  controls_pre = controls_pre, 
  lambda_list = c(-Inf), 
  n_starts = 1, 
  silent = TRUE)

# matches using the DID target weights 
target_weights <- Target_Weights(controls_pre)
optim_result2 <- NULL
optim_result2$solution <- c(0, target_weights)

# --- Predict the post intervention average ----

# estimate the effects
cf.est <- Post_Intervention_Mu3(
  optim_result = optim_result, 
  d = d, 
  treated_unit_id = treated_unit_id, 
  treated_time = treated_time, 
  controls_pre = controls_pre,
  treated_pre = treated_pre, 
  controls_post = controls_post) |> 
  as.numeric()

athey_att <- Athey_Method(d)

cf.true <- mean(subset(d$cf.poisson, d$unit_id == 1 & d$treated == 1))
obs <- mean(subset(d$actual.poisson, d$unit_id == 1 & d$treated == 1))

results <- tibble(
  obs = obs, 
  cf.estimate = cf.est, 
  cf.true = cf.true,
  true_att = obs - cf.true, 
  poisson_sdid_att = obs - cf.est,
  athey_sdid_att = athey_att) 


results |> view()




# --- Inference ----

placebos <- Improved_Placebo(d = d, 
                             starts_lambdas = 10)

k <- Poisson_SDID4(d = d, 
                   lambda_list = c(-Inf))
k
