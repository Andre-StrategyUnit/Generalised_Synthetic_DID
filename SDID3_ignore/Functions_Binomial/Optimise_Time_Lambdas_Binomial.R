# -------------- optimiser over lambda --------------
optimise_time_lambdas_Binomial <- function(
    initial_par,
    y_pre_successes, 
    y_pre_trials, 
    y_post_successes_mean, 
    y_post_trials_mean, 
    sum_to_one_lambda,  # scalar or vector
    lb,
    ub
) {
  optimise_single <- function(lambda) {
    res <- nloptr(
      x0         = initial_par,
      eval_f     = function(p) neg_ll_sdid_time_binomial(
        par                       = p,
        y_pre_successes           = y_pre_successes,
        y_pre_trials              = y_pre_trials,
        y_post_successes_mean     = y_post_successes_mean,
        y_post_trials_mean        = y_post_trials_mean,
        sum_to_one_lambda         = lambda),
      eval_grad_f= function(p) neg_ll_sdid_time_grad_binomial(
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
    
    # compute predicted p
    eta_pred <- sum(res$solution * (y_pre_successes / y_pre_trials))
    p_pred   <- plogis(eta_pred)
    
    tibble(
      sum_to_one_lambda    = lambda,
      raw_objective        = res$objective,
      sum_par_hat          = (sum(res$solution) - 1)^2,
      penalty              = lambda * sum_par_hat,
      unpenalised_objective= res$objective - (lambda * sum_par_hat),
      sum_par              = sum(res$solution),
      predicted_post_p     = p_pred
    )
  }
  
  # vectorize over lambda grid
  map_dfr(sum_to_one_lambda, optimise_single)
}
