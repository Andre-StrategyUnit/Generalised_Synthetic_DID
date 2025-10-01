

library(nloptr)
library(tidyverse)

source("Functions/Poisson_Functions.R") # create dataset

tictoc::tic()

# simulate data ----
d <- Create_SC_Data(treated_units = 2,
                    time_n = 50,
                    units_n = 12, 
                    impact_factor = 0)

d <- d |> 
  select(unit_id, time,
         treated, treat_group,
         actual.poisson)

# central estimate based on full sample ----
full_results <- Multi_SDID_Poisson(d = d, 
                                    treated_list = c(1,2), 
                                    n_starts = 1)

# jack knife uncertainty ----
jack_results <- JackKnife_Poisson_SDID(d = d, 
                                        treated_list = c(1,2))


# inferential ----
inferential_results <- Poisson_SDID_Inferential(
  full_results = full_results, 
  jack_results = jack_results, 
  alpha = 0.05)
inferential_results
tictoc::toc()

# plot ----
plot_data <- Multi_SDID_Poisson_VC(d = d,
                                   treated_list = c(1,2),
                                   n_starts = 1) 

plot_data2 <- plot_data|>
  mutate(time = 1:nrow(plot_data))

plot_data3 <- plot_data2 |>
  pivot_longer(cols = c(predicted_sum, actual_sum), 
               names_to = "Estimand", 
               values_to = "Value") |>
  mutate(Estimand = case_when(
    Estimand == "predicted_sum" ~ "Counterfactual Count", 
    Estimand == "actual_sum" ~ "Observed Count"))


plot_data3 |>
  ggplot(mapping = aes(x = time, 
                       y = Value, 
                       col = Estimand)) +
  geom_line() +
  geom_vline(xintercept = min(subset(d$time,
                                     d$treated == 1))-0.5) +
  theme_minimal() + 
  ggtitle(label = "Title")
