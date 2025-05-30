Binomial_SDID_Inferential <- function(full_results, 
                                     jack_results, 
                                     alpha){
  
  central_estimate <- full_results |>
    mutate(impact = binomial_sdid_att * treated_post_trials.sum) |>
    summarise(impact = sum(impact)) |> 
    pull(impact)
  
  jack_results2 <- jack_results |>
    mutate(impact = binomial_sdid_att * treated_post_trials.sum) |>
    group_by(excluded) |>
    summarise(impact = sum(impact))
  
  # inferential statistics ----
  N <- length(jack_results2$impact)
  mean_impact <- mean(jack_results2$impact)
  
  # this is the formula for the jack knife variance, 
  # it is different to the standard formula 
  # that is a feature not a bug
  var_jack <- ((N - 1) / N) * sum((jack_results2$impact- mean_impact)^2)
  se_jack <- sqrt(var_jack)
  
  
  
  t_stat <- central_estimate/se_jack
  df <- N - 1
  
  # p value
  p_value <- 2 * pt(-abs(t_stat), df = df)
  
  # CIs
  alpha <- 0.1
  t_crit <- qt(1 - alpha / 2, df = df)
  
  ci_lower <- central_estimate - t_crit * se_jack
  ci_upper <- central_estimate + t_crit * se_jack
  
  results3 <- tibble(
    aggregate_impact = central_estimate, 
    ci_lower = ci_lower, 
    ci_upper = ci_upper, 
    p_value = p_value, 
    sum_trials = sum(full_results$treated_post_trials.sum),
    sum_timepoints = sum(full_results$post_timepoint),
    sum_units = length(unique(full_results$treated_unit)))
  
  return(results3)
  
}
