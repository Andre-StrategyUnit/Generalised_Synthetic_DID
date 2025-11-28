multi_start_nloptr <- function(
    n_starts = 10,
    controls_pre, treated_pre,
    sum_to_one_lambda = sum_to_one_lambda,
    ridge_lambda = ridge_lambda, 
    lb, ub,
    silent = FALSE) {

  
  results <- list()
  
  for (i in 1:n_starts) {
    if (!silent) cat("Multi-start iteration", i, "\n")
    
    # identify the starting intercept
    initial_intercept <- Target_Intercept(controls_pre = controls_pre, 
                                          treated_pre = treated_pre)
    
    # identify the starting weights
    # target weights
    initial_weights <- Target_Weights(controls_pre)
    
    initial_par <- c(initial_intercept, initial_weights)
    
    optim_result <- nloptr::nloptr(
      x0 = initial_par,
      eval_f = function(par) neg_loglik_poisson_local(
        par = par,
        controls_pre = controls_pre,
        treated_pre = treated_pre,
        sum_to_one_lambda = sum_to_one_lambda, 
        ridge_lambda = ridge_lambda),
      eval_grad_f = function(par) neg_loglik_poisson_grad(
        par = par,
        controls_pre = controls_pre,
        treated_pre = treated_pre,
        sum_to_one_lambda = sum_to_one_lambda, 
        ridge_lambda = ridge_lambda),
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
