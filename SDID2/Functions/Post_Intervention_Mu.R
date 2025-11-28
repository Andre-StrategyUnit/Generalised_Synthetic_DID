# Estimate Post Treatment Mean # 

Post_Intervention_Mu <- function(
    optim_result, 
    d, 
    treated_pre, 
    treated_unit_id, treated_time, 
    controls_pre, controls_post, 
    best_lambda_row){
  
  controls_post <- d %>%
    filter(unit_id != treated_unit_id, time >= treated_time) %>%
    select(unit_id, time, cf.poisson) %>%
    pivot_wider(names_from = unit_id, values_from = cf.poisson) %>%
    arrange(time) %>%
    select(-time) %>%
    as.matrix()
  
  initial_par <- c(0, rep(x = 1/(nrow(controls_pre)-1), 
                          each = nrow(controls_pre)-1))
  
  lb <- rep(x = 0, each = nrow(controls_pre))
  ub <- rep(x = 1, each = nrow(controls_pre))
  
  # optimise sdid time weights
  
  # extract weights from the optimal solution
  intercept_opt <- optim_result$solution[1]
  optimal_weights <- optim_result$solution[-1]
  
  # extract VC pre/post
  vc_pre  <- Extract_VC_Pre_Poisson2(d, 
                                     optim_result, 
                                     treated_unit_id, 
                                     treated_time)
  vc_post <- Extract_VC_Post_Poisson2(d, 
                                      optim_result, 
                                      treated_unit_id, 
                                      treated_time)
  
  # Inputs
  y_pre <- vc_pre
  y_post_mean <- mean(vc_post)
  
  nll_values <- numeric()
  sum_to_ones <- numeric()
  predicteds <- numeric()
  
  result2 <- nloptr::nloptr(
    x0 = rep(1 / length(y_pre), length(y_pre)),
    eval_f = function(par) 
      neg_ll_sdid_time(par, y_pre, y_post_mean, 
                       best_lambda_row$sum_to_one_lambda),
    eval_grad_f = function(par) 
      neg_ll_sdid_time_grad(par, y_pre, y_post_mean, 
                            best_lambda_row$sum_to_one_lambda),
    lb = rep(0, length(y_pre)),
    ub = rep(1, length(y_pre)),
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel = 1e-6,
      maxeval = 2000
    )
  )
  
  obj = (result2$solution %*% treated_pre)[1,1]
  return(obj)
  
}
Post_Intervention_Mu2 <- function(
    optim_result, 
    d, 
    treated_pre, 
    treated_unit_id, 
    treated_time, 
    controls_pre, 
    controls_post, 
    best_lambda_row
) {
  # --- Get synthetic pre/post predictions ---
  vc_pre  <- Extract_VC_Pre_Poisson2(d, 
                                     optim_result, 
                                     treated_unit_id, 
                                     treated_time)
  vc_post <- Extract_VC_Post_Poisson2(d, 
                                      optim_result, 
                                      treated_unit_id, 
                                      treated_time)
  
  # --- Prepare inputs for time weighting ---
  y_pre       <- vc_pre                   # on count scale
  y_post_mean <- mean(vc_post)
  T_pre       <- length(y_pre)
  
  # Initial par and bounds (intercept + weights)
  x0 <- c(0, rep(1 / T_pre, T_pre))
  lb <- c(-Inf, rep(0, T_pre))
  ub <- c( Inf, rep(1, T_pre))
  
  # Optimize time weights (with intercept)
  res2 <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(par) 
      neg_ll_sdid_time2(par, y_pre, y_post_mean, 
                        sum_to_one_lambda = best_lambda_row$sum_to_one_lambda),
    eval_grad_f = function(par) 
      neg_ll_sdid_time_grad2(par, y_pre, y_post_mean, 
                             sum_to_one_lambda = best_lambda_row$sum_to_one_lambda),
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel  = 1e-6,
      maxeval   = 2000
    )
  )
  
  # Extract solution
  intercept_time <- res2$solution[1]
  weights_time   <- res2$solution[-1]
  
  # --- Apply to treated_pre (convert to log scale) ---
  log_treated_pre <- log(pmin(pmax(treated_pre, 1e-10), 1e10))  # clamp to avoid -Inf
  mu_pred <- exp(intercept_time + sum(weights_time * log_treated_pre))
  
  return(mu_pred)
}

Post_Intervention_Mu3 <- function(
    optim_result, 
    d, 
    treated_pre, 
    treated_unit_id, 
    treated_time, 
    controls_pre, 
    controls_post, 
    best_lambda_row
) {
  # --- Get synthetic pre/post predictions ---
  vc_pre  <- Extract_VC_Pre_Poisson2(d, 
                                     optim_result, 
                                     treated_unit_id, 
                                     treated_time)
  vc_post <- Extract_VC_Post_Poisson2(d, 
                                      optim_result, 
                                      treated_unit_id, 
                                      treated_time)
  
  # --- Prepare inputs for time weighting ---
  y_pre       <- vc_pre                   # on count scale
  y_post_mean <- mean(vc_post)
  T_pre       <- length(y_pre)
  
  # Initial par and bounds (intercept + weights)
  x0 <- c(0, rep(1 / T_pre, T_pre))
  lb <- c(-Inf, rep(0, T_pre))
  ub <- c( Inf, rep(1, T_pre))
  
  # Optimize time weights (with intercept)
  res2 <- nloptr::nloptr(
    x0 = x0,
    eval_f = function(par) 
      neg_ll_sdid_time3(par, 
                        y_pre, 
                        y_post_mean),
    eval_grad_f = function(par) 
      neg_ll_sdid_time_grad3(par, 
                             y_pre, 
                             y_post_mean),
    lb = lb,
    ub = ub,
    opts = list(
      algorithm = "NLOPT_LD_MMA",
      xtol_rel  = 1e-6,
      maxeval   = 2000
    )
  )
  
  # Extract solution
  intercept_time <- res2$solution[1]
  weights_time   <- res2$solution[-1]
  
  mu_pred <- intercept_time + sum(weights_time * y_pre)
  mu_pred <- ifelse(test = mu_pred > 0, 
                    yes = mu_pred, 
                    no = 0)
  
  return(mu_pred)
}
