# --- Define loss function (no log transform, raw counts) ----
neg_loglik_poisson_local <- function(par, 
                                     controls_pre, treated_pre, 
                                     sum_to_one_lambda = 0,
                                     ridge_lambda = 0.1) {
  intercept <- par[1]
  weights   <- par[-1]
  
  # Clamp and transform
  log_controls <- log(pmin(pmax(controls_pre, 1e-10), 1e+10))  # T x J
  linear_pred  <- intercept + as.vector(log_controls %*% weights)
  
  mu <- exp(linear_pred)
  mu <- pmin(pmax(mu, 1e-8), 1e8)  # clamp predictions
  
  # Poisson negative log-likelihood
  nll <- -sum(dpois(x = treated_pre, lambda = mu, log = TRUE))
  
  # target weights
  target_weights <- Target_Weights(controls_pre)
  
  # Penalties
  sum1_penalty  <- sum_to_one_lambda * (sum(weights) - 1)^2
  ridge_penalty <- ridge_lambda * sum((weights-target_weights)^2)
  
  return(nll + sum1_penalty + ridge_penalty)
}


# --- Define gradient of loss function ----
neg_loglik_poisson_grad <- function(par, 
                                    controls_pre, 
                                    treated_pre, 
                                    sum_to_one_lambda = 0, 
                                    ridge_lambda = 0) {
  intercept <- par[1]
  weights   <- par[-1]
  
  log_controls <- log(pmin(pmax(controls_pre, 1e-10), 1e+10))  # T x J
  eta          <- intercept + as.vector(log_controls %*% weights)
  mu           <- exp(eta)
  mu           <- pmin(pmax(mu, 1e-8), 1e8)
  
  residual <- mu - treated_pre  # derivative of -loglik wrt eta for Poisson
  
  # Gradients
  grad_intercept <- sum(residual)
  grad_weights   <- t(log_controls) %*% residual
  
  # target weights
  target_weights <- Target_Weights(controls_pre)
  
  # Penalty gradients
  grad_weights <- grad_weights +
    2 * sum_to_one_lambda * (sum(weights) - 1) +
    2 * ridge_lambda * (weights-target_weights)
  
  return(c(grad_intercept, as.vector(grad_weights)))
}
