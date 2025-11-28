library(tidyverse)


rm(list = ls())
a <- Sys.time()

# Simulate data ----
source("Functions/Poisson_Functions.R")
treated_units = c(1:2)
d <- Create_SC_Data(units_n = 22,
                    good_ctrls = 12,
                    time_n = 25, 
                    treated_time = 25,
                    treated_units = max(treated_units), 
                    scaling_factor = 5, 
                    impact_factor = 0)

# Run a single analysis ----
tictoc::tic()
output <- Multi_SDID_Poisson(d = d,
                             n_starts = 1,
                             treated_list = treated_units, 
                             lambda_list = c(-Inf, 
                                             -3, -1, 0, 1, 3, 5))
tictoc::toc()
# view(output)

central_overview <- Aggregate_Central(output)
# view(central_overview)

# Run all the individual placebo analyses ----

tictoc::tic()
placebos <- Improved_Placebo2(
  d = d, 
  main_result = output)
tictoc::toc()

# view(placebos)

# Aggregate all the placebo analyses ----

tictoc::tic()
output2 <- Aggregate_Placebos(
  placebos = placebos, 
  iter = 1000)
tictoc::toc()

# ggplot(data = output2, 
#        mapping = aes(x = (impact))) + 
#   geom_histogram(bins = 20)

# Calculate inferential statistics ----

inferential_results <- SDID_Inferential(
  central_overview = central_overview, 
  aggregate_placebos = output2,
  p_direction = "two_sided", 
  quantiles = c(0.25, 0.975)
)

inferential_results




subset(d, treated == 1) |> 
  summarise(cf = sum(cf.poisson), 
            obs = sum(actual.poisson), 
            impact = log(obs/cf))

b <- Sys.time()
b-a
