optimise_unit_lambdas <- function(
    lambdas,
    initial_par,
    controls_pre,
    treated_pre,
    sum_to_one_lambda,
    lb,
    ub
) {
  # Function to optimize for a single lambda
  optimise_single_lambda <- function(ridge_lambda, sum_to_one_lambda) {
    optim_result <- nloptr::nloptr(
      x0 = initial_par,
      eval_f = function(par) neg_loglik_poisson_local(
        par = par,
        controls_pre = controls_pre,
        treated_pre = treated_pre,
        sum_to_one_lambda = sum_to_one_lambda,
        ridge_lambda = ridge_lambda),
      eval_grad_f = function(par) neg_loglik_poisson_grad(
        par = par,
        controls_pre = controls_pre,
        treated_pre = treated_pre,
        sum_to_one_lambda = sum_to_one_lambda,
        ridge_lambda = ridge_lambda),
      lb = lb,
      ub = ub,
      opts = list(
        algorithm = "NLOPT_LD_MMA",
        xtol_rel = 1e-6,
        maxeval = 2000,
        print_level = 0
      )
    )
    
    results <- tibble(
      ridge_lambda = ridge_lambda,
      sum_to_one_lambda = sum_to_one_lambda, 
      raw_objective = optim_result$objective, 
      sum_par_hat = sum(optim_result$solution[-1]^2),
      sum_par2_hat = (sum(optim_result$solution[-1]) - 1)^2,
      penalty = (ridge_lambda * sum_par_hat) + (sum_to_one_lambda * sum_par2_hat), 
      sum_par = sum(optim_result$solution[-1]),
      unpenalised_objective = raw_objective - penalty
    )
    
    return(results)
  }
  
  inputs <- tibble(
    ridge_lambda = rep(x = lambdas, 
                       each = length(sum_to_one_lambda)), 
    sum_to_one_lambda = rep(x = sum_to_one_lambda, 
                            times = length(lambdas))
  )
  
  # Apply optimization across all lambdas
  nll_values <- purrr::pmap(inputs, 
                            optimise_single_lambda) |>
    bind_rows()
  
  # Return results 
  return(nll_values)
  
}
