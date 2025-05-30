
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

# --- Create dataset ----

d <- Create_SC_Data(time_n = 100,
                    treated_time = 80, 
                    units_n = 100, 
                    scaling_factor = 25, 
                    unit_fe_sd = 3, 
                    impact_factor = 1)

treated_unit_id <- 1
treated_time <- min(d$time[d$treated == 1])

# --- Extract pre-period data ----

# Raw treated outcomes (counts)
treated_pre <- Extract_Treated(d)

# Raw controls outcomes (counts)
controls_pre <- Extract_Controls_Pre(d)

# Number of control units
n_controls <- ncol(controls_pre)

# Extract optimal model -----

optim_result <- Optimiser_Wrapper(
  treated_pre = treated_pre, 
  controls_pre = controls_pre, 
  n_starts = 1, 
  silent = TRUE)

# Plot pre match ----

pre_match <- Extract_Pre_Match(treated_pre = treated_pre, 
                               optim_result = optim_result)
ggplot(data = pre_match,
       mapping = aes(x = time)) +
  geom_line(aes(y = synthetic_counts), color = "blue") +
  geom_line(aes(y = treated_counts), color = "red") +
  theme_minimal() +
  labs(y = "Counts",
       title = "Synthetic Control: Treated vs Synthetic (Raw Counts Model)",
       subtitle = "w/ intercept, optimized by MMA")


# --- Predict the post ----

# try a few values of lambda 
time_lambdas <- Determine_Sum_To_One_Lambda(
  d = d, 
  optim_result = optim_result, 
  controls_pre = controls_pre, 
  list_of_lambdas = c(exp(seq(from = -10, to = 20, 
                              by = 2))), 
  tolerance = 0.01) 

# select the best one
best_lambda_row <- time_lambdas |>
  subset(within_tolerance == TRUE) |>
  slice_min(order_by = unpenalised_objective, 
            n = 1)

# estimate the effects
cf.est <- Post_Intervention_Mu(
  optim_result = optim_result, 
  d = d, 
  treated_pre = treated_pre, 
  controls_post = controls_post, 
  best_lambda_row = best_lambda_row)

athey_att <- Athey_Method(d)

cf.true <- mean(subset(d$cf.poisson, d$unit_id == 1 & d$treated == 1))
obs <- mean(subset(d$actual.poisson, d$unit_id == 1 & d$treated == 1))

results <- tibble(
  obs = obs, 
  cf.true = cf.true,
  cf.est = cf.est,
  true_att = obs-cf.true, 
  poisson_sdid_att = obs - cf.est,
  poisson_sdid_att.adj = obs - (cf.est / best_lambda_row$predicted_post_mean), 
  athey_sdid_att = athey_att) 
  

results |> view()



