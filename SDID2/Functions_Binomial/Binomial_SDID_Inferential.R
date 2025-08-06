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
  alpha <- alpha
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


Binomial_SDID_Inferential.permutation <- function(full_results, 
                                                 perm_results, 
                                                 alpha){
  
  central_estimate <- full_results |>
    mutate(obs_success = round(obs* treated_post_trials.sum, 0) , 
           cf_success = round(cf.est * treated_post_trials.sum, 0)) |>
    summarise(obs_success = sum(obs_success), 
              cf_success = sum(cf_success), 
              trials = sum(treated_post_trials.sum)) |>
    ungroup() |>
    mutate(obs_odds = obs_success/(trials-obs_success), 
           cf_odds = cf_success/(trials-cf_success), 
           scaled_effect = log(obs_odds) - log(cf_odds)) 
  
  true_effect <- central_estimate$scaled_effect
  perm_effects <- as.vector(perm_results$scaled_effect)
  
  se <- sqrt(mean(perm_effects^2))
  n_permutations <- length(perm_effects)
  t_crit <- qt(1- (alpha/2), df = n_permutations - 1)
  ci.lwr <- true_effect - (t_crit * se)
  ci.upr <- true_effect + (t_crit * se)
  
  tbl <- central_estimate |>
    mutate(
      scaled_ci.lwr = ci.lwr, 
      scaled_ci.upr = ci.upr,
      cf_odds_ci.lwr = obs_odds/exp(scaled_ci.lwr),
      cf_odds_ci.upr = obs_odds/exp(scaled_ci.upr),
      cf_outcome_ci.lwr = (trials/(1+cf_odds_ci.lwr))*cf_odds_ci.lwr, 
      cf_outcome_ci.upr = (trials/(1+cf_odds_ci.upr))*cf_odds_ci.upr,
      exact_p =  (sum(abs(perm_effects) >= abs(true_effect)))/(length(perm_effects)), 
      adj_exact_p = (1 + sum(abs(perm_effects) >= abs(true_effect)))/(length(perm_effects)+1)
  )
  
  return(tbl)
  
}
