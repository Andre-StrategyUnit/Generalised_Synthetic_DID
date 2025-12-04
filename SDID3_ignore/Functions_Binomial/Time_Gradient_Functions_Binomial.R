# -------------- likelihood & gradient --------------
neg_ll_sdid_time_binomial <- function(par, 
                                      y_pre_successes, 
                                      y_pre_trials, 
                                      y_post_successes_mean, 
                                      y_post_trials_mean, 
                                      sum_to_one_lambda = 0) {
  # Pre‑period proportions
  mean_prop_pre <- y_pre_successes / y_pre_trials
  
  # Linear predictor
  eta_pred <- sum(par * mean_prop_pre)
  eta_pred <- pmin(pmax(eta_pred, -10), 10)
  
  # Predicted probability
  p_pred <- plogis(eta_pred)
  p_pred <- pmin(pmax(p_pred, 1e-8), 1 - 1e-8)
  
  # Binomial negative log‑likelihood
  nll <- -dbinom(
    x    = round(y_post_successes_mean), 
    size = round(y_post_trials_mean), 
    prob = p_pred, 
    log  = TRUE
  )
  
  # Sum‑to‑one penalty
  penalty <- sum_to_one_lambda * (sum(par) - 1)^2
  total   <- nll + penalty
  
  if (!is.finite(total)) {
    warning("Non‑finite total in neg_ll_sdid_time_binomial()")
  }
  total
}

neg_ll_sdid_time_grad_binomial <- function(par, 
                                           y_pre_successes, 
                                           y_pre_trials, 
                                           y_post_successes_mean, 
                                           y_post_trials_mean, 
                                           sum_to_one_lambda = 0) {
  mean_prop_pre <- y_pre_successes / y_pre_trials
  eta_pred <- sum(par * mean_prop_pre)
  eta_pred <- pmin(pmax(eta_pred, -10), 10)
  
  p_pred <- plogis(eta_pred)
  p_pred <- pmin(pmax(p_pred, 1e-8), 1 - 1e-8)
  
  # derivative of NLL wrt eta
  d_nll_d_eta <- (p_pred - (y_post_successes_mean / y_post_trials_mean)) * y_post_trials_mean
  
  # chain rule through eta = sum(w * mean_prop_pre)
  grad <- d_nll_d_eta * mean_prop_pre
  
  # penalty gradient
  grad <- grad + 2 * sum_to_one_lambda * (sum(par) - 1)
  grad[!is.finite(grad)] <- 0
  grad
}