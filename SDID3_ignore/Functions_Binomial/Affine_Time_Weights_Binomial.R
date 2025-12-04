### Updated SDID Binomial functions with joint intercept estimation

# -------------- likelihood & gradient (with intercept) --------------
neg_ll_sdid_time_binomial2 <- function(par, 
                                       y_pre_successes, 
                                       y_pre_trials, 
                                       y_post_successes_mean, 
                                       y_post_trials_mean, 
                                       sum_to_one_lambda = 0) {
  # Parameters: intercept + time weights
  intercept <- par[1]
  weights   <- par[-1]
  
  # Convert pre-period mean proportions to logit scale
  prop_pre  <- y_pre_successes / y_pre_trials
  prop_pre  <- pmin(pmax(prop_pre, 1e-8), 1 - 1e-8)
  logit_pre <- log(prop_pre / (1 - prop_pre))
  
  # Affine predictor on logit scale
  eta_pred <- intercept + sum(weights * logit_pre)
  eta_pred <- pmin(pmax(eta_pred, -10), 10)  # stabilize
  
  # Map to probability
  p_pred <- plogis(eta_pred)
  p_pred <- pmin(pmax(p_pred, 1e-8), 1 - 1e-8)
  
  # Binomial likelihood
  nll <- -dbinom(
    x    = round(y_post_successes_mean),
    size = round(y_post_trials_mean),
    prob = p_pred,
    log  = TRUE
  )
  
  # Sum-to-one penalty
  penalty <- sum_to_one_lambda * (sum(weights) - 1)^2
  total <- nll + penalty
  
  if (!is.finite(total)) warning("Non-finite total in neg_ll_sdid_time_binomial2()")
  return(total)
}


neg_ll_sdid_time_grad_binomial2 <- function(par, 
                                            y_pre_successes, 
                                            y_pre_trials, 
                                            y_post_successes_mean, 
                                            y_post_trials_mean, 
                                            sum_to_one_lambda = 0) {
  intercept <- par[1]
  weights   <- par[-1]
  
  # Pre-period logits
  prop_pre  <- y_pre_successes / y_pre_trials
  prop_pre  <- pmin(pmax(prop_pre, 1e-8), 1 - 1e-8)
  logit_pre <- log(prop_pre / (1 - prop_pre))
  
  eta_pred <- intercept + sum(weights * logit_pre)
  eta_pred <- pmin(pmax(eta_pred, -10), 10)
  p_pred   <- plogis(eta_pred)
  p_pred   <- pmin(pmax(p_pred, 1e-8), 1 - 1e-8)
  
  # Gradient of NLL w.r.t. eta
  d_nll_d_eta <- (p_pred - (y_post_successes_mean / y_post_trials_mean)) * y_post_trials_mean
  
  # Gradient w.r.t. intercept and weights
  grad_intercept <- d_nll_d_eta
  grad_weights   <- d_nll_d_eta * logit_pre
  
  # Penalty gradient on weights
  grad_weights <- grad_weights + 2 * sum_to_one_lambda * (sum(weights) - 1)
  
  return(c(grad_intercept, grad_weights))
}


# -------------- optimiser over lambda --------------
optimise_time_lambdas_Binomial2 <- function(
    initial_par,
    y_pre_successes, 
    y_pre_trials, 
    y_post_successes_mean, 
    y_post_trials_mean, 
    sum_to_one_lambda,  # scalar or vector
    lb,
    ub
) {
  optimise_single2 <- function(lambda) {
    res <- nloptr::nloptr(
      x0         = initial_par,
      eval_f     = function(p) neg_ll_sdid_time_binomial2(
        par                       = p,
        y_pre_successes           = y_pre_successes,
        y_pre_trials              = y_pre_trials,
        y_post_successes_mean     = y_post_successes_mean,
        y_post_trials_mean        = y_post_trials_mean,
        sum_to_one_lambda         = lambda),
      eval_grad_f= function(p) neg_ll_sdid_time_grad_binomial2(
        par                       = p,
        y_pre_successes           = y_pre_successes,
        y_pre_trials              = y_pre_trials,
        y_post_successes_mean     = y_post_successes_mean,
        y_post_trials_mean        = y_post_trials_mean,
        sum_to_one_lambda         = lambda),
      lb         = lb,
      ub         = ub,
      opts       = list(
        algorithm  = "NLOPT_LD_MMA",
        xtol_rel   = 1e-6,
        maxeval    = 2000,
        print_level= 0
      )
    )
    # extract solution
    sol <- res$solution
    intercept <- sol[1]
    weights   <- sol[-1]
    # compute predicted p
    eta_pred <- intercept + sum(weights * (y_pre_successes / y_pre_trials))
    p_pred   <- plogis(eta_pred)
    
    tibble(
      sum_to_one_lambda     = lambda,
      raw_objective         = res$objective,
      sum_par_hat           = (sum(weights) - 1)^2,
      penalty               = lambda * sum_par_hat,
      unpenalised_objective = res$objective - penalty,
      sum_par               = sum(weights),
      predicted_post_p      = p_pred
    )
  }
  purrr::map_dfr(sum_to_one_lambda, optimise_single2)
}

# -------------- wrapper for λ grid search --------------
Determine_Sum_To_One_Lambda_Binomial2 <- function(
    optim_result,
    y_pre_successes,    # vector length T_pre
    y_pre_trials,       # vector length T_pre
    y_post_successes,   # vector length T_post
    y_post_trials,      # vector length T_post
    list_of_lambdas     = exp(seq(-10, 10, by = 1)),
    tolerance           = 0.01
) {
  # Extract VC from optim_result
  y_post_mean_succ   <- mean(y_post_successes)
  y_post_mean_trials <- mean(y_post_trials)
  
  T_pre <- length(y_pre_successes)
  # initial: include intercept + equal weights
  initial_par <- c(0, rep(1/T_pre, T_pre))
  # lower/upper: intercept unbounded, weights ≥0
  lb <- c(-Inf, rep(0, T_pre))
  ub <- c( Inf, rep(1, T_pre))
  
  result_tbl <- optimise_time_lambdas_Binomial2(
    initial_par,
    y_pre_successes,
    y_pre_trials,
    y_post_mean_succ,
    y_post_mean_trials,
    list_of_lambdas,
    lb,
    ub
  )
  result_tbl %>%
    mutate(within_tolerance = between(sum_par, 1 - tolerance, 1 + tolerance))
}

# -------------- post-intervention estimator --------------
Post_Intervention_Mu_Binomial2 <- function(
    optim_result, 
    d, 
    treated_unit_id, 
    treated_time,
    treated_post_trials,
    treated_pre_trials, 
    wt_lb = 0, wt_ub = +Inf, 
    lambda_star
) {
  # --- Extract synthetic control (unit-weighted) ---
  
  # Synthetic logit per time point (pre/post)
  exp_succ_pre <- Extract_VC_Pre_Binomial2(
    d = d, 
    optim_result = optim_result, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time, 
    treated_pre_trials = treated_pre_trials)
  exp_succ_post <- Extract_VC_Post_Binomial2(
    d = d, 
    optim_result = optim_result, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time, 
    treated_post_trials = treated_post_trials)
  
  # --- Fit time weights on pre-period logits ---
  T_pre <- length(exp_succ_pre)
  init2 <- c(0, rep(1 / T_pre, T_pre))   # intercept + equal weights
  lb2   <- c(-Inf, rep(wt_lb, T_pre))
  ub2   <- c(+Inf, rep(wt_ub, T_pre))
  
  res2 <- nloptr::nloptr(
    x0 = init2,
    eval_f = function(par) neg_ll_sdid_time_binomial2(
      par,
      y_pre_successes       = exp_succ_pre,
      y_pre_trials          = treated_pre_trials,
      y_post_successes_mean = mean(exp_succ_post),
      y_post_trials_mean    = mean(treated_post_trials),
      sum_to_one_lambda     = lambda_star
    ),
    eval_grad_f = function(par) neg_ll_sdid_time_grad_binomial2(
      par,
      y_pre_successes       = exp_succ_pre,
      y_pre_trials          = treated_pre_trials,
      y_post_successes_mean = mean(exp_succ_post),
      y_post_trials_mean    = mean(treated_post_trials),
      sum_to_one_lambda     = lambda_star
    ),
    lb   = lb2,
    ub   = ub2,
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel  = 1e-6,
      maxeval   = 2000
    )
  )
  
  # Final time-weighted logit prediction
  par2 <- res2$solution
  intercept_time <- par2[1]
  weights_time   <- par2[-1]
  
  treated_pre <- Extract_Treated.bn(
    d = d, 
    treated_unit_id = treated_unit_id, 
    treated_time = treated_time)
  
  # Calculate logits from fractional expected successes and trials
  prop_pre <- treated_pre / treated_pre_trials
  prop_pre <- pmin(pmax(prop_pre, 1e-8), 1 - 1e-8)
  logit_pre <- log(prop_pre / (1 - prop_pre))
  
  # Then apply the time weights:
  mean_logit_pred <- sum(weights_time * logit_pre) + intercept_time
  mean_pre_predict <- plogis(mean_logit_pred)
  
  return(mean_pre_predict)
}
