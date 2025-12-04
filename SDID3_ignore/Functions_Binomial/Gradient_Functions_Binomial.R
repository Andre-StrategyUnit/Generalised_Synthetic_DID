neg_loglik_binomial_local <- function(par,
                                      controls_pre,         # T × J matrix of control successes
                                      controls_trials_pre,  # T × J matrix of control trials
                                      treated_pre,          # vector length T (successes)
                                      treated_trials_pre,   # vector length T (trials)
                                      sum_to_one_lambda = 0,
                                      ridge_lambda = 0) {
  intercept <- par[1]
  weights   <- par[-1]
  
  time_n <- nrow(controls_pre)
  nll <- 0
  
  for (t in seq_len(time_n)) {
    # Logit-scale synthetic controls
    prop_t <- controls_pre[t, ] / controls_trials_pre[t, ]
    prop_t <- pmin(pmax(prop_t, 1e-8), 1 - 1e-8)  # clamp
    logit_t <- qlogis(prop_t)
    
    eta_t <- intercept + sum(logit_t * weights)
    p_t <- plogis(eta_t)
    p_t <- pmin(pmax(p_t, 1e-8), 1 - 1e-8)  # clamp again for binomial
    
    y_t <- treated_pre[t]
    n_t <- treated_trials_pre[t]
    
    nll_t <- -dbinom(y_t, size = n_t, prob = p_t, log = TRUE)
    nll <- nll + nll_t
  }
  
  # Penalties
  sum1_penalty  <- sum_to_one_lambda * (sum(weights) - 1)^2
  ridge_penalty <- ridge_lambda * sum(weights^2)
  
  return(nll + sum1_penalty + ridge_penalty)
}

neg_loglik_binomial_grad <- function(par,
                                     controls_pre,         # T × J (successes)
                                     controls_trials_pre,  # T × J (trials)
                                     treated_pre,          # length T (successes)
                                     treated_trials_pre,   # length T (trials)
                                     sum_to_one_lambda = 0,
                                     ridge_lambda = 0.1) {
  intercept <- par[1]
  weights   <- par[-1]
  
  # Compute control logits (logit-scale)
  props <- controls_pre / controls_trials_pre
  props <- pmin(pmax(props, 1e-8), 1 - 1e-8)
  control_logits <- qlogis(props)  # T × J
  
  # Linear predictor
  eta <- intercept + as.vector(control_logits %*% weights)  # length T
  
  # Predicted probabilities
  p <- plogis(eta)
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)  # clamp
  
  y <- treated_pre
  n <- treated_trials_pre
  
  # Residuals
  resid <- y - n * p
  
  # Derivative of NLL w.r.t. eta
  dL_deta <- -resid * p * (1 - p)
  
  # Gradient
  grad_intercept <- sum(dL_deta)
  grad_weights   <- t(control_logits) %*% dL_deta  # J × 1
  
  # Penalty gradients
  grad_weights <- grad_weights +
    2 * sum_to_one_lambda * (sum(weights) - 1) +
    2 * ridge_lambda * weights
  
  return(c(grad_intercept, as.vector(grad_weights)))
}
