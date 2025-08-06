

library(nloptr)
library(tidyverse)

source("Functions_Binomial/Binomial_Functions.R") # create dataset

# simulate data ----
d <- Create_SC_Data_Binomial(treated_units = 2,
                             time_n = 50,
                             units_n = 15,
                             trials_lambda = 5, 
                             impact_factor = 0)

d <- d |>
  select(unit_id, time,
         treated, treat_group,
         actual.binomial, trials)

# central estimate based on full sample ----
full_results <- Multi_SDID_Binomial(d = d, 
                                    treated_list = c(1,2), 
                                    n_starts = 10)


# jack knife uncertainty ----
jack_results <- JackKnife_Binomial_SDID(d = d, 
                                        treated_list = c(1,2))


# outputs
inferential_results <- Binomial_SDID_Inferential(
  full_results = full_results, 
  jack_results = jack_results, 
  alpha = 0.05
)

# plot ----
plot_data <- Multi_SDID_Binomial_VC(d = d,
                                   treated_list = c(1,2),
                                   n_starts = 1) 

plot_data2 <- plot_data|>
  mutate(time = 1:nrow(plot_data), 
         predicted_proportion = (predicted_sum/trials_sum), 
         actual_proportion = actual_sum/trials_sum)

plot_data3 <- plot_data2 |>
  pivot_longer(cols = c(predicted_sum, actual_sum, 
                        predicted_proportion, actual_proportion), 
               names_to = "Estimand", 
               values_to = "Value") |>
  mutate(Estimand = case_when(
    Estimand == "predicted_sum" ~ "Counterfactual Count", 
    Estimand == "actual_sum" ~ "Observed Count", 
    Estimand == "predicted_proportion" ~ "Counterfactual Rate", 
    Estimand == "actual_proportion" ~ "Observed Rate"))


plot_data3 |>
  subset(Estimand %in% c("Observed Rate", "Counterfactual Rate")) |>
  ggplot(mapping = aes(x = time, 
                       y = Value, 
                       col = Estimand)) +
  geom_line() +
  geom_vline(xintercept = min(subset(d$time,
                                     d$treated == 1))-0.5) +
  theme_minimal() + 
  scale_y_continuous(limits = c(0, 1)) + 
  ggtitle(label = "Title")

inferential_results 
