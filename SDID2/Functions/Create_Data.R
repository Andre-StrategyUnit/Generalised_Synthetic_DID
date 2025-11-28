
# Simulate once ----

Create_SC_Data <- function(
    units_n = 80, 
    good_ctrls = 40, 
    scaling_factor = 5, 
    time_n = 50, 
    treated_units = 1,
    impact_factor = 1,
    unit_fe_sd = 1, 
    intercept2 = TRUE, 
    alpha2 = 0.5, # this default is the elastic net
    treated_time = round(0.9*time_n, 0)){
  
  # create latent dataset
  latent_data <- tibble(
    time = 1:time_n, 
    factor1 = arima.sim(model = list(order = c(2, 1, 2), 
                                     ar = c(+0.9, -0.1), 
                                     ma = c(-0.5, +0.5)), 
                        n = time_n)[1:time_n], 
    factor2 = arima.sim(model = list(order = c(2, 1, 2), 
                                     ar = c(+0.9, -0.1), 
                                     ma = c(-0.5, +0.5)), 
                        n = time_n)[1:time_n])
  
  # create factor weightings
  factor_weightings <- tibble(
    unit_id = 1:units_n, 
    factor1.wts = c(runif(n = good_ctrls, 
                          min = 0.7, max = 1), 
                    runif(n = (units_n-good_ctrls), 
                          min = 0, max = 0.3)),
    factor2.wts = c(0,
                    runif(n = good_ctrls-1, 
                          min = -0.1, max = 0.1),
                    runif(n = (units_n-good_ctrls), 
                          min = 0.8, max = 1.2)))
  
  # create original sample outcome dataset
  sample_data1 <- tibble(
    unit_id = rep(x = c(1:units_n), 
                  each = time_n), 
    time = rep(x = c(1:time_n),
               times = units_n), 
    treat_group = ifelse(test = unit_id <= treated_units, 
                         yes = 1, 
                         no = 0), 
    start_date = ifelse(test = treat_group == 1, 
                        yes = treated_time, 
                        no = NA), 
    treated = ifelse(test = treat_group == 1 & time >= start_date, 
                     yes = 1, no = 0 )) 
  
  # join together with the factor weightings
  sample_data2 <- left_join(sample_data1, factor_weightings, by = "unit_id")
  sample_data3 <- left_join(sample_data2, latent_data, by = "time")
  
  # wrangle into final version of the sample dataset
  sample_data4 <- sample_data3 |>
    mutate(
      unit_fe = rep(x = rnorm(n = units_n, 
                              mean = 0, 
                              sd = unit_fe_sd), 
                    each = time_n), 
      cf.gauss = unit_fe + 
        (factor_weightings$factor1.wts[unit_id] * latent_data$factor1[time]) + 
        (factor_weightings$factor2.wts[unit_id] * latent_data$factor2[time]) + 
        rnorm(n = units_n * time_n, 
              mean = 0, 
              sd = 0.5), 
      cf.gauss1 = ifelse(test = unit_id == 1, 
                         yes = unit_fe + latent_data$factor1[time] + rnorm(n = units_n * time_n, mean = 0, sd = 0.5), 
                         no = cf.gauss), 
      impact = ifelse(test = treated == 1, 
                      yes = (time/time_n) * impact_factor * sd(latent_data$factor1), 
                      no = 0), 
      actual.gauss = cf.gauss1 + impact, 
      cf.gauss2 = cf.gauss1/sd(cf.gauss1), 
      actual.gauss2 = actual.gauss/sd(cf.gauss1)) |>
    select(unit_id, time, 
           treat_group, start_date, treated, 
           cf.gauss2, actual.gauss2) |>
    mutate(cf.poisson = round(x = scaling_factor * ((exp(1))^cf.gauss2), digits = 0), 
           actual.poisson = round(x = scaling_factor * ((exp(1))^actual.gauss2), digits = 0))
  
  return(sample_data4)
  
}


