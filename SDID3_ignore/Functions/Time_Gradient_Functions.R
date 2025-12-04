# loss function for sdid on time weights

neg_ll_sdid_time <- function(par, y_pre, y_post_mean, sum_to_one_lambda = 0) {
  mu_pred <- sum(par * y_pre)
  mu_pred <- pmin(pmax(mu_pred, 1e-2), 1e8)
  
  nll <- -(y_post_mean * log(mu_pred) - mu_pred - lgamma(y_post_mean + 1))
  
  penalty <- sum_to_one_lambda * (sum(par) - 1)^2
  total <- nll + penalty
  
  if (!is.finite(total)) {
    print("---- Non-finite total detected ----")
    print(paste("mu_pred:", mu_pred))
    print(paste("y_post_mean:", y_post_mean))
    print(paste("log(mu_pred):", log(mu_pred)))
    print(paste("penalty:", penalty))
    print(paste("par (weights):", paste(round(par, 4), collapse = ", ")))
  }
  
  return(total)
}


# gradient of loss function for sdid on time weights ----

neg_ll_sdid_time_grad <- function(par, y_pre, y_post_mean, sum_to_one_lambda = 0) {
  mu_pred <- sum(par * y_pre)
  mu_pred <- pmin(pmax(mu_pred, 1e-8), 1e8)
  
  # Derivative of negative log-likelihood w.r.t. pi
  grad <- (1 - y_post_mean / mu_pred) * y_pre
  
  # Penalty: sum-to-one soft constraint
  penalty_grad <- 2 * sum_to_one_lambda * (sum(par) - 1)
  
  grad_final <- grad + penalty_grad
  grad_final[!is.finite(grad_final)] <- 0
  
  return(grad_final)
}
neg_ll_sdid_time2 <- function(par, y_pre, y_post_mean, sum_to_one_lambda = 0) {
  intercept <- par[1]
  weights   <- par[-1]
  
  # Log-transform y_pre (with clamping)
  log_y_pre <- log(pmin(pmax(y_pre, 1e-10), 1e10))
  
  # Linear predictor and predicted mean
  eta     <- intercept + sum(weights * log_y_pre)
  mu_pred <- exp(eta)
  mu_pred <- pmin(pmax(mu_pred, 1e-2), 1e8)
  
  # Poisson negative log-likelihood
  nll <- -(y_post_mean * log(mu_pred) - mu_pred - lgamma(y_post_mean + 1))
  
  # Penalty for sum(weights) ≈ 1
  penalty <- sum_to_one_lambda * (sum(weights) - 1)^2
  total   <- nll + penalty
  
  return(total)
}

neg_ll_sdid_time_grad2 <- function(par, y_pre, y_post_mean, sum_to_one_lambda = 0) {
  intercept <- par[1]
  weights   <- par[-1]
  
  # Log-transform y_pre (with clamping)
  log_y_pre <- log(pmin(pmax(y_pre, 1e-10), 1e10))
  
  # Linear predictor and predicted mean
  eta     <- intercept + sum(weights * log_y_pre)
  mu_pred <- exp(eta)
  mu_pred <- pmin(pmax(mu_pred, 1e-8), 1e8)
  
  # Derivative of Poisson NLL wrt eta
  d_nll_d_eta <- mu_pred - y_post_mean
  
  # Gradients
  grad_intercept <- d_nll_d_eta
  grad_weights   <- d_nll_d_eta * log_y_pre
  
  # Penalty gradient
  grad_weights <- grad_weights + 2 * sum_to_one_lambda * (sum(weights) - 1)
  
  return(c(grad_intercept, grad_weights))
}
