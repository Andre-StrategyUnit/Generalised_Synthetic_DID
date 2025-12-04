
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

# --- Create dataset ----

n = n+1
set.seed(12)
good_ctrls.v <- 40
scale_fctr.v <- 7.5
impact_fctr.v <- -1

d <- Create_SC_Data(time_n = 100,
                    treated_time = 80, 
                    treated_units = 1,
                    units_n = 40, 
                    good_ctrls = good_ctrls.v, 
                    scaling_factor = scale_fctr.v, 
                    unit_fe_sd = 7.5, 
                    impact_factor = impact_fctr.v)

treated_unit_id <- 1
treated_time <- min(d$time[d$treated == 1])

d |> 
  mutate(treatment_status = ifelse(test = treat_group == 0, 
                                   yes = "Control", 
                                   no = "Treated")) |> 
  ggplot() + 
  # geom_line(mapping = aes(x = time,
  #                         y = cf.poisson,
  #                         group = unit_id),
  #           col = "green") +
  geom_line(mapping = aes(x = time, 
                          y = (actual.poisson), 
                          group = unit_id, 
                          col = unit_id > 10
                          )) + 
  theme_minimal() + 
  facet_grid(treatment_status~.) + 
  geom_vline(xintercept = treated_time-1) + 
  ggtitle("Count of an example outcome over time") + 
  ylab("Count") + 
  xlab("Time") + 
  theme(legend.position="none")  + scale_color_manual(values = c("blue", "red"))
  

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
plot(target_weights, colMeans(controls_pre))


# --- Plot pre match ----

# matches using the optimiser with no regularisation
pre_match <- Extract_Pre_Match(treated_pre = treated_pre, 
                               optim_result = optim_result, 
                               controls_pre = controls_pre)

# matches using the did method, identify the control weights
pre_match2 <- Extract_Pre_Match(treated_pre = treated_pre, 
                               optim_result = optim_result2, 
                               controls_pre = controls_pre) |> 
  rename(did_estimate = synthetic_counts)

# matches using the did method, identify the intercept
pre_match2b <- pre_match2 |> 
  mutate(intercept = sum(pre_match2$treated_counts)/sum(pre_match2$did_estimate), 
         did_est2 = intercept * did_estimate)

# join
pre_match3 <- left_join(pre_match, pre_match2b)

ggplot(data = pre_match3,
       mapping = aes(x = time)) +
  geom_line(aes(y = (did_est2)), color = "darkgreen") +
  geom_line(aes(y = (synthetic_counts)), color = "blue") +
  geom_line(aes(y = (treated_counts)), color = "red") +
  theme_minimal() +
  labs(y = "Counts",
       title = "Pre Period Synthetic Control",
       subtitle = "Treated vs Unregularized Synthetic vs DID")

# --- Plot post match ----

# mathces using the optimiser with no regularisation
post_match <- Extract_Pre_Match(treated_pre = Extract_Treated_Post(d = d, 
                                                                   treated_unit_id = treated_unit_id, 
                                                                   treated_time = treated_time), 
                                optim_result = optim_result, 
                                controls_pre = Extract_Controls_Post(d = d, 
                                                                     treated_unit_id = treated_unit_id, 
                                                                     treated_time = treated_time))

# matches using the did method, identify the control weights
post_match2 <- Extract_Pre_Match(treated_pre = Extract_Treated_Post(d = d, 
                                                                   treated_unit_id = treated_unit_id, 
                                                                   treated_time = treated_time), 
                                optim_result = optim_result2, 
                                controls_pre = Extract_Controls_Post(d = d, 
                                                                     treated_unit_id = treated_unit_id, 
                                                                     treated_time = treated_time)) |> 
  rename(did_estimate = synthetic_counts)

# matches using the did method, identify the intercept
post_match2b <- post_match2 |> 
  mutate(intercept = sum(post_match2$treated_counts)/sum(post_match2$did_estimate), 
         did_est2 = intercept * did_estimate)

# join together
post_match3 <- left_join(post_match, post_match2b)

ggplot(data = post_match3,
       mapping = aes(x = time)) +
  geom_line(aes(y = (did_est2)), color = "darkgreen") +
  geom_line(aes(y = (synthetic_counts)), color = "blue") +
  geom_line(aes(y = (treated_counts)), color = "red") +
  theme_minimal() +
  labs(y = "Counts",
       title = "Post Period Synthetic Control",
       subtitle = "Treated vs Unregularized Synthetic vs DID")

# --- Predict the post intervention average ----

# # try a few values of lambda 
# time_lambdas <- Determine_Sum_To_One_Lambda(
#   d = d, 
#   treated_unit_id = treated_unit_id, 
#   treated_time = treated_time,
#   optim_result = optim_result, 
#   controls_pre = controls_pre, 
#   list_of_lambdas = c(0, 
#                       exp(seq(from = -10, to = 20, 
#                               by = 2))), 
#   tolerance = 0.01) 
# 
# # select the best one
# best_lambda_row <- time_lambdas |>
#   subset(within_tolerance == TRUE) |>
#   slice_min(order_by = unpenalised_objective, 
#             n = 1)

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



