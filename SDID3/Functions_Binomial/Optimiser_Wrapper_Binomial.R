# Run Multi Start Optimisation Wrapper # 

Optimiser_Wrapper_Binomial <- function(controls_pre = controls_pre, 
                                       controls_trials_pre = controls_pre_trials, 
                                       treated_pre = treated_pre, 
                                       treated_trials_pre = treated_pre_trials,
                                       n_starts, 
                                       silent, 
                                       wt_lb = 0, wt_ub = 1, 
                                       ridge_lambda = 0, 
                                       sum_to_one_lambda = 0) {
  
  
  nll_values <- numeric()
  
  initial_par <- c(0, rep(x = 1/ncol(controls_pre), 
                          each = ncol(controls_pre)))
  
  lb <- c(-Inf, rep(x = wt_lb, 
                    each = ncol(controls_pre)))
  ub <- c(+Inf, rep(x = wt_ub, 
                    each = ncol(controls_pre)))

  
  # run multi start optimisations
  multi_optim_result <- multi_start_nloptr_binomial(
    n_starts = n_starts,
    controls_pre = controls_pre, 
    controls_trials_pre = controls_trials_pre, 
    treated_pre = treated_pre, 
    treated_trials_pre = treated_trials_pre,
    sum_to_one_lambda = sum_to_one_lambda,
    ridge_lambda = ridge_lambda, 
    lb = lb,
    ub = ub, 
    silent = silent
  )
  
  # extract the best solution
  optim_result <- multi_optim_result[[n_starts+1]] # no here is just n_starts + 1
  
  return(optim_result)
}
