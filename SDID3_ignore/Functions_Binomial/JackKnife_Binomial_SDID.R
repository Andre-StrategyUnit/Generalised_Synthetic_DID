# 
# d <- Create_SC_Data_Binomial(treated_units = 2,
#                              time_n = 50,
#                              units_n = 10,
#                              trials_lambda = 50)
# 
# d <- d |>
#   select(unit_id, time, 
#          treated, treat_group, 
#          actual.binomial, trials)


JackKnife_Binomial_SDID <- function(d, treated_list){
  
control_units <- unique(subset(d$unit_id, d$treat_group == 0))
treated_list = treated_list

results2 <- tibble(
  cf.est = numeric(), 
  obs = numeric(),  
  vc_post_avg = numeric(),  
  treated_post_trials.sum = numeric(), 
  binomial_sdid_att = numeric(), 
  post_timepoint = numeric(), 
  excluded = numeric())

for (i in 1:length(control_units)){

d2 <- d |>
  subset(treat_group == 1 | treat_group == 0 & unit_id != control_units[i])

jack_result <- Multi_SDID_Binomial(d2, 
                                treated_list = treated_list, 
                                n_starts = 1)
jack_result$excluded <- control_units[i]

results2 <- bind_rows(results2, jack_result)

}

return(results2)

}
# 
# jack_results <- JackKnife_Binomial_SDID(d = d, treated_list = c(1,2))
# 
# jack_results2 <- jack_results |>
#   mutate(impact = binomial_sdid_att * treated_post_trials.sum) |>
#   group_by(excluded) |>
#   summarise(impact = sum(impact))
# 
# N <- length(jack_results2$impact)
# mean_impact <- mean(jack_results2$impact)
# var_jack <- (N - 1) / N * sum((jack_results2$impact- mean_impact)^2)
# se_jack <- sqrt(var_jack)
