multi_start_nloptr_binomial <- function(
    n_starts = 10,
    controls_pre,         # T × J matrix of control successes
    controls_trials_pre,  # T × J matrix of control trials
    treated_pre,          # vector length T (successes)
    treated_trials_pre,   # vector length T (trials)
    sum_to_one_lambda = 0,
    ridge_lambda = 0, 
    lb, ub,
    silent = FALSE) {
  
  results <- list()
  
  for (i in 1:n_starts) {
    if (!silent) cat("Multi-start iteration", i, "\n")
    
    initial_intercept <- runif(1, min = -1, max = 1)
    initial_weights <- runif(ncol(controls_pre))
    initial_weights <- initial_weights / sum(initial_weights)  # normalize
    
    initial_par <- c(initial_intercept, initial_weights)
    
    optim_result <- nloptr::nloptr(
      x0 = initial_par,
      eval_f = function(par) neg_loglik_binomial_local(
        par = par,
        controls_pre = controls_pre, 
        controls_trials_pre = controls_trials_pre, 
        treated_pre = treated_pre, 
        treated_trials_pre = treated_trials_pre,
        ridge_lambda = ridge_lambda,
        sum_to_one_lambda = sum_to_one_lambda),
      eval_grad_f = function(par) neg_loglik_binomial_grad(
        par = par,
        controls_pre = controls_pre, 
        controls_trials_pre = controls_trials_pre, 
        treated_pre = treated_pre, 
        ridge_lambda = ridge_lambda,
        treated_trials_pre = treated_trials_pre,
        sum_to_one_lambda = sum_to_one_lambda),
      lb = lb,
      ub = ub,
      opts = list(
        algorithm = "NLOPT_LD_MMA",
        xtol_rel = 1e-6,
        maxeval = 2500,
        print_level = 0
      )
    )
    
    results[[i]] <- optim_result
  }
  
  # Pick the best solution (smallest objective)
  results[[n_starts + 1]] <- results[[which.min(sapply(results, function(x) x$objective))]]
  
  return(results)
}
