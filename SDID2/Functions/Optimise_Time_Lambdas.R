# Function to optimize for a single lambda
optimise_single_lambda2 <- function(lambda, y_pre, y_post_mean){
  
  initial_par <- c(0, rep(1 / length(y_pre), length(y_pre)))
  lb <- c(-Inf, rep(0, length(y_pre)))
  ub <- c(+Inf, rep(1, length(y_pre)))
  
  optim_result <- nloptr::nloptr(
    x0 = initial_par,
    eval_f = function(par) neg_ll_sdid_time2(
      par = par,
      y_pre = y_pre, 
      y_post_mean = y_post_mean, 
      sum_to_one_lambda = lambda),
    eval_grad_f = function(par) neg_ll_sdid_time_grad2(
      par = par,
      y_pre = y_pre, 
      y_post_mean = y_post_mean, 
      sum_to_one_lambda = lambda),
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel = 1e-6,
      maxeval = 2000,
      print_level = 0
    )
  )
  
  # extract results
  intercept_hat <- optim_result$solution[1]
  weights_hat   <- optim_result$solution[-1]
  eta_hat       <- intercept_hat + sum(weights_hat * y_pre)
  mu_hat        <- exp(eta_hat)
  
  penalty <- lambda * (sum(weights_hat) - 1)^2
  
  tibble(
    sum_to_one_lambda     = lambda, 
    raw_objective         = optim_result$objective, 
    sum_par_hat           = (sum(weights_hat) - 1)^2,
    penalty               = penalty, 
    sum_par               = sum(weights_hat),
    unpenalised_objective = optim_result$objective - penalty,
    predicted_post_mean   = mu_hat / y_post_mean
  )
}


optimise_time_lambdas2 <- function(
    initial_par,
    controls_pre,
    treated_pre,
    y_pre, 
    y_post_mean, 
    sum_to_one_lambda,
    lb,
    ub
) {
    # Grid over lambdas
    purrr::map_dfr(sum_to_one_lambda, optimise_single_lambda2)
  }





