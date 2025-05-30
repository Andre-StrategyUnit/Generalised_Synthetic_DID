library(tidyverse)

# Simulate once ----
Create_SC_Data_Binomial <- function(
    units_n = 80, 
    scaling_factor = 5, 
    time_n = 50, 
    treated_units = 1,
    impact_factor = 1,
    unit_fe_sd = 1, 
    intercept2 = TRUE, 
    alpha2 = 0.5,
    treated_time = round(0.9*time_n, 0),
    trials_lambda = 250  # trials drawn from pois distr w this lambda
){
  
  # Latent factor simulation
  latent_data <- tibble(
    time = 1:time_n, 
    factor1 = arima.sim(model = list(order = c(2, 1, 2), 
                                     ar = c(+0.9, -0.05), 
                                     ma = c(-0.5, +0.5)), 
                        n = time_n)[1:time_n], 
    factor2 = arima.sim(model = list(order = c(2, 1, 2), 
                                     ar = c(-0.6, +0.2), 
                                     ma = c(-0.5, +0.5)), 
                        n = time_n)[1:time_n])
  
  factor_weightings <- tibble(
    unit_id = 1:units_n, 
    factor1.wts = c(runif(4, 0.7, 1), runif(units_n - 4, 0, 0.3)),
    factor2.wts = c(0, runif(3, -0.1, 0.1), runif(units_n - 4, 0.8, 1.2)))
  
  sample_data1 <- tibble(
    unit_id = rep(1:units_n, each = time_n), 
    time = rep(1:time_n, times = units_n), 
    treat_group = ifelse(unit_id <= treated_units, 1, 0), 
    start_date = ifelse(unit_id <= treated_units, treated_time, NA), 
    treated = ifelse(unit_id <= treated_units & time >= treated_time, 1, 0))
  
  sample_data2 <- left_join(sample_data1, factor_weightings, by = "unit_id")
  sample_data3 <- left_join(sample_data2, latent_data, by = "time")
  
  sample_data4 <- sample_data3 |>
    mutate(
      unit_fe = rep(rnorm(units_n, mean = 0, sd = unit_fe_sd), each = time_n),
      cf.gauss = unit_fe +
        (factor1.wts * factor1) +
        (factor2.wts * factor2) +
        rnorm(units_n * time_n, mean = 0, sd = 0.5),
      cf.gauss1 = ifelse(unit_id == 1,
                         unit_fe + factor1 + rnorm(units_n * time_n, mean = 0, sd = 0.5),
                         cf.gauss),
      impact = ifelse(treated == 1,
                      (time/time_n) * impact_factor * sd(factor1),
                      0),
      actual.gauss = cf.gauss1 + impact,
      cf.gauss2 = (cf.gauss1 / sd(cf.gauss1)),
      actual.gauss2 = (actual.gauss / sd(cf.gauss1)),
      cf.gauss3 = scales::rescale(x = cf.gauss2, 
                                  to = c(-2, +2), 
                                  from = c(min(actual.gauss2), max(actual.gauss2))),
      actual.gauss3 = scales::rescale(x = actual.gauss2, 
                                      to = c(-2, +2), 
                                      from = c(min(actual.gauss2), max(actual.gauss2))),
      cf.prob = plogis(cf.gauss3),  # logit link
      actual.prob = plogis(actual.gauss3),
      trials = rpois(n = units_n * time_n, 
                     lambda = trials_lambda),
      cf.binomial = rbinom(n = n(), size = trials, prob = cf.prob),
      actual.binomial = rbinom(n = n(), size = trials, prob = actual.prob)
    ) |>
    select(unit_id, time, 
           treat_group, start_date, treated, 
           cf.gauss3, actual.gauss3,
           cf.prob, actual.prob, 
           trials, 
           cf.binomial, actual.binomial)
  
  return(sample_data4)
}

# obj <- Create_SC_Data_Binomial(impact_factor = 1) 
# 
# obj |>
#   ggplot(mapping = aes(x = time,
#                        y = actual.binomial/trials,
#                        group = unit_id,
#                        col = time < start_date)) +
#   geom_line() +
#   facet_grid(unit_id==1~.)


