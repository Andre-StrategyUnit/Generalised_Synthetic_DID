Poisson_SDID_Inferential <- function(full_results, 
                                     jack_results, 
                                     alpha){
  
  central_estimate <- full_results |>
    mutate(impact = poisson_sdid_att * post_timepoints) |>
    summarise(full_impact = sum(impact)) |> 
    pull(full_impact)
  
  jack_results2 <- jack_results |>
    mutate(impact = poisson_sdid_att * post_timepoints) |>
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
  alpha <- alpha
  t_crit <- qt(1 - alpha / 2, df = df)
  
  ci_lower <- central_estimate - t_crit * se_jack
  ci_upper <- central_estimate + t_crit * se_jack
  
  results3 <- tibble(
    aggregate_impact = central_estimate, 
    ci_lower = ci_lower, 
    ci_upper = ci_upper, 
    p_value = p_value, 
    sum_timepoints = sum(full_results$post_timepoints),
    sum_units = length(unique(full_results$treated_unit)))
  
  return(results3)
  
}
