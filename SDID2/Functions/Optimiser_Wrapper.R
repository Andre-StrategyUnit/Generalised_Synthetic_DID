# Run Multi Start Optimisation Wrapper # 

Optimiser_Wrapper <- function(treated_pre, 
                              controls_pre, 
                              n_starts, silent) {
  
  
  nll_values <- numeric()
  
  initial_par <- c(0, rep(x = 1/ncol(controls_pre), 
                          each = ncol(controls_pre)))
  
  lb <- c(-Inf, rep(x = 0, 
                    each = ncol(controls_pre)))
  ub <- c(+Inf, rep(x = 1, 
                    each = ncol(controls_pre)))
  
  # run multi start optimisations
  multi_optim_result <- multi_start_nloptr(
    n_starts = n_starts,
    controls_pre = controls_pre,
    treated_pre = treated_pre,
    sum_to_one_lambda = 0,
    ridge_lambda = 0, 
    lb = lb,
    ub = ub, 
    silent = silent
  )
  
  # extract the best solution
  optim_result <- multi_optim_result[[n_starts+1]] # no here is just n_starts + 1
  
  return(optim_result)
}
