Post_Intervention_Mu_Binomial <- function(
    optim_result, 
    d, 
    treated_unit_id, 
    treated_time,
    treated_post_trials,
    treated_pre_trials, 
    best_lambda_row
) {
  # Extract post-period control successes and trials
  controls_post <- d %>%
    filter(unit_id != treated_unit_id, time >= treated_time) %>%
    select(unit_id, time, actual.binomial) %>%
    arrange(time) %>%
    pivot_wider(names_from = unit_id, values_from = actual.binomial) %>%
    select(-time) %>%
    as.matrix()
  
  controls_post_trials <- d %>%
    filter(unit_id != treated_unit_id, time >= treated_time) %>%
    select(unit_id, time, trials) %>%
    arrange(time) %>%
    pivot_wider(names_from = unit_id, values_from = trials) %>%
    select(-time) %>%
    as.matrix()
  
  # Extract control pre-period data from earlier pipeline
  controls_pre <- Extract_Controls_Pre.bn(d)
  controls_pre_trials <- Extract_Controls_Pre.bn_trials(d)
  
  # Get intercept and weights
  intercept_opt <- optim_result$solution[1]
  optimal_weights <- optim_result$solution[-1]
  
  # Compute synthetic logit-scale predictions (pre and post)
  synthetic_logit_post <- intercept_opt + as.vector((controls_post / controls_post_trials) %*% optimal_weights)
  synthetic_logit_pre <- intercept_opt + as.vector((controls_pre / controls_pre_trials) %*% optimal_weights)
  
  # Convert to predicted proportions
  synthetic_p_post <- plogis(synthetic_logit_post)
  synthetic_p_pre  <- plogis(synthetic_logit_pre)
  
  # Multiply by trials to get expected successes
  expected_successes_post <- synthetic_p_post * treated_post_trials
  expected_successes_pre  <- synthetic_p_pre * treated_pre_trials
  
  # Run time-weight optimization to weight pre-period estimates
  y_pre_successes <- expected_successes_pre
  y_pre_trials    <- treated_pre_trials 
  y_post_successes_mean <- mean(expected_successes_post)
  y_post_trials_mean    <- mean(treated_post_trials)
  
  result2 <- nloptr::nloptr(
    x0 = rep(1 / length(y_pre_successes), length(y_pre_successes)),
    eval_f = function(par) neg_ll_sdid_time_binomial(
      par,
      y_pre_successes = y_pre_successes,
      y_pre_trials = y_pre_trials,
      y_post_successes_mean = y_post_successes_mean,
      y_post_trials_mean = y_post_trials_mean,
      sum_to_one_lambda = best_lambda_row$sum_to_one_lambda
    ),
    eval_grad_f = function(par) neg_ll_sdid_time_grad_binomial(
      par,
      y_pre_successes = y_pre_successes,
      y_pre_trials = y_pre_trials,
      y_post_successes_mean = y_post_successes_mean,
      y_post_trials_mean = y_post_trials_mean,
      sum_to_one_lambda = best_lambda_row$sum_to_one_lambda
    ),
    lb = rep(0, length(y_pre_successes)),
    ub = rep(1, length(y_pre_successes)),
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel = 1e-6,
      maxeval = 2000
    )
  )
  
  # Return weighted average of post-period synthetic p
  time_weights <- result2$solution
  estimated_post_success_rate <- sum(time_weights * synthetic_p_pre)
  
  return(estimated_post_success_rate)
}
