Poisson_SDID_Inferential.jack <- function(full_results, 
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

Poisson_SDID_Inferential.permutation <- function(full_results, 
                                                 perm_results, 
                                                 alpha){
  
  central_estimate <- full_results |>
    mutate(obs_outcome_aggr = obs * post_timepoints, 
           cf_outcome_aggr = cf.est * post_timepoints) |>
    summarise(obs_outcome_aggr = sum(obs_outcome_aggr), 
              cf_outcome_aggr = sum(cf_outcome_aggr), 
              scaled_effect = log(obs_outcome_aggr) - log(cf_outcome_aggr)) |> 
    ungroup()
  
  true_effect <- central_estimate$scaled_effect
  perm_effects <- as.vector(perm_results$scaled_effect)
  
  exact_p_value <- (1 + sum(abs(perm_effects) >= abs(true_effect)))/(length(perm_effects)+1)
  se <- sqrt(mean(perm_effects^2))
  n_permutations <- length(perm_effects)
  t_crit <- qt(1- (alpha/2), df = n_permutations - 1)
  ci.lwr <- true_effect - (t_crit * se)
  ci.upr <- true_effect + (t_crit * se)
  
  
  tbl <- tibble(
    obs_outcome_aggr = central_estimate$obs_outcome_aggr, 
    cf_outcome_aggr = central_estimate$cf_outcome_aggr, 
    meaningful_effect = central_estimate$obs_outcome_aggr - central_estimate$cf_outcome_aggr, 
    scaled_effect = log(central_estimate$obs_outcome_aggr) - log(central_estimate$cf_outcome_aggr), 
    log_ci.lwr = ci.lwr, 
    log_ci.upr = ci.upr, 
    cf_outcome_ci.lwr = (exp(log_ci.lwr) * obs_outcome_aggr), 
    cf_outcome_ci.upr = (exp(log_ci.upr) * obs_outcome_aggr), 
    exact_p =  (sum(abs(perm_effects) >= abs(true_effect)))/(length(perm_effects)), 
    adj_exact_p = (1 + sum(abs(perm_effects) >= abs(true_effect)))/(length(perm_effects)+1)
  )
    
    return(tbl)
    
}
