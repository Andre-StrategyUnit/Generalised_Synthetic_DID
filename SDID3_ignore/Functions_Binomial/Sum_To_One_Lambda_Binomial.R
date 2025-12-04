Determine_Sum_To_One_Lambda_Binomial <- function(
    optim_result,
    
    # Pre‑period treated outcomes:
    y_pre_successes,    # vector length T_pre
    y_pre_trials,       # vector length T_pre
    
    # Post‑period treated outcomes (as vectors):
    y_post_successes,   # vector length T_post
    y_post_trials,      # vector length T_post
    
    list_of_lambdas     = exp(seq(-10, 10, by = 1)),
    tolerance           = 0.01
) {
  # 1) Extract the weights/intercept from your first‑stage SC fit
  intercept_opt   <- optim_result$solution[1]
  control_weights <- optim_result$solution[-1]
  
  # 2) Collapse post‑period into means
  y_post_mean_succ <- mean(y_post_successes)
  y_post_mean_trials <- mean(y_post_trials)
  
  # 3) Set up the time‑weight problem
  T_pre <- length(y_pre_successes)
  
  initial_par <- rep(1 / T_pre, T_pre)
  lb          <- rep(0, T_pre)
  ub          <- rep(1, T_pre)
  
  # 4) Grid‑search your λ penalty
  result_tbl <- optimise_time_lambdas_Binomial(
    initial_par               = initial_par,
    y_pre_successes           = y_pre_successes,
    y_pre_trials              = y_pre_trials,
    y_post_successes_mean     = y_post_mean_succ,
    y_post_trials_mean        = y_post_mean_trials,
    sum_to_one_lambda         = list_of_lambdas,
    lb                         = lb,
    ub                         = ub
  )
  
  # 5) Flag those within tolerance of sum‑to‑one
  result_tbl %>%
    mutate(
      within_tolerance = between(sum_par, 1 - tolerance, 1 + tolerance)
    )
}

