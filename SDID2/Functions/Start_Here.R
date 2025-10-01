library(tidyverse)

source("Functions/Poisson_Functions.R")
treated_units = c(1:30)
d <- Create_SC_Data(units_n = 130,
                    good_ctrls = 80,
                    time_n = 50, 
                    treated_units = max(treated_units), 
                    scaling_factor = 10)
# view(d)

tictoc::tic()
output <- Multi_SDID_Poisson(d = d,
                             n_starts = 1,
                             treated_list = treated_units, 
                             lambda_list = -Inf)
tictoc::toc()
# view(output)

starts_lambdas <- tibble(
  equiv_unit = treated_units,
  starts = 1 + max(d$time) - output$post_timepoints, 
  lambdas = log10(output$units_lambda_used))

tictoc::tic()
placebos <- Improved_Placebo(
  d = d, 
  starts_lambdas = starts_lambdas)
tictoc::toc()

# view(placebos)

output2 <- tibble(
  cf = numeric(), 
  obs = numeric(), 
  impact = numeric())

for (i in 1:10000){
  
  output2[i,] <- placebos |> 
    group_by(equiv_unit) |> 
    sample_n(size = 1) |> 
    ungroup() |> 
    summarise(cf = sum(cf.est * post_timepoints), 
              obs = sum(obs * post_timepoints), 
              impact = log(obs/cf))
}


ggplot(data = output2, 
       mapping = aes(x = (impact))) + 
  geom_histogram()

mean(output2$impact)

log(sum(d$actual.poisson)/sum(d$cf.poisson))
