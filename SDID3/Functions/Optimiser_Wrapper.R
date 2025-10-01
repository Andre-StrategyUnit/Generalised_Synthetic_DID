# Run Multi Start Optimisation Wrapper # 

Optimiser_Wrapper <- function(treated_pre, 
                              controls_pre, 
                              n_starts, silent, 
                              lambda_list = c(-3,-1,0, 1,3,5)) {
  
  
  nll_values <- numeric()
  
  initial_par <- c(0, rep(x = 1/ncol(controls_pre), 
                          each = ncol(controls_pre)))
  
  lb <- c(-Inf, rep(x = 0, 
                    each = ncol(controls_pre)))
  ub <- c(+Inf, rep(x = 1, 
                    each = ncol(controls_pre)))
  
  controls_pre.odd <- controls_pre[c(seq(1, nrow(controls_pre), by = 2)),]
  controls_pre.even <- controls_pre[c(seq(2, nrow(controls_pre), by = 2)),]
  treated_pre.odd <- treated_pre[seq(1, length(treated_pre), by = 2)]
  treated_pre.even <- treated_pre[seq(2, length(treated_pre), by = 2)]
  
  results78 <- tibble(
    ridge_lambda = numeric(), 
    scaled_effect = numeric())
  
  iter = 0 
  
  for (l in 10^lambda_list){
  # run multi start optimisations
  multi_optim_result <- multi_start_nloptr(
    n_starts = n_starts,
    controls_pre = controls_pre.odd,
    treated_pre = treated_pre.odd,
    sum_to_one_lambda = 0,
    ridge_lambda = l, 
    lb = lb,
    ub = ub, 
    silent = silent
  )
  
  iter = iter + 1
  # extract the best solution
  optim_result <- multi_optim_result[[n_starts+1]] # no here is just n_starts + 1
  
  intercept_opt <- optim_result$solution[1]
  weights_opt   <- optim_result$solution[-1]
  
  # clamp then log-transform
  y_mat        <- pmin(pmax(controls_pre.even, 1e-10), 1e+10)
  log_controls <- log(y_mat)
  
  # affine on log scale → log μ
  log_pred         <- intercept_opt + log_controls %*% weights_opt
  expected_success <- exp(log_pred)
  
  results78b <- tibble(
    ridge_lambda = l, 
    scaled_effect = abs(log(sum(treated_pre.even)/sum(expected_success[,1]))))
  
  results78[iter,] <- results78b
  
  }
  
  best_lambda <- subset(results78, scaled_effect == min(results78$scaled_effect))$ridge_lambda
  
  multi_optim_result <- multi_start_nloptr(
    n_starts = n_starts,
    controls_pre = controls_pre,
    treated_pre = treated_pre,
    sum_to_one_lambda = 0,
    ridge_lambda = best_lambda, 
    lb = lb,
    ub = ub, 
    silent = silent
  )
  
  # extract the best solution
  optim_result <- multi_optim_result[[n_starts+1]] # no here is just n_starts + 1
  
  
  return(optim_result)
}
