Determine_Sum_To_One_Lambda <- function(
    d, 
    optim_result, 
    controls_pre, 
    treated_unit_id, 
    treated_time,
    list_of_lambdas = exp(seq(from = -10, to = 10, by = 1)), 
    tolerance = 0.01
) {
  # Extract optimal weights
  intercept_opt   <- optim_result$solution[1]
  optimal_weights <- optim_result$solution[-1]
  
  # Create synthetic control predictions (pre/post)
  controls_post <- d %>%
    filter(unit_id != treated_unit_id, time >= treated_time) %>%
    select(unit_id, time, actual.poisson) %>%
    pivot_wider(names_from = unit_id, values_from = actual.poisson) %>%
    arrange(time) %>%
    select(-time) %>%
    as.matrix()
  
  # extract VC pre/post
  vc_pre  <- Extract_VC_Pre_Poisson2(d, 
                                     optim_result, 
                                     treated_unit_id, 
                                     treated_time)
  vc_post <- Extract_VC_Post_Poisson2(d, 
                                      optim_result, 
                                      treated_unit_id, 
                                      treated_time)
  
  # Setup optimization
  T_pre <- length(vc_pre)
  initial_par <- c(0, rep(1 / T_pre, T_pre))  # intercept + weights
  lb <- c(-Inf, rep(0, T_pre))
  ub <- c( Inf, rep(1, T_pre))
  
  # Grid search
  result_tbl <- optimise_time_lambdas2(
    initial_par       = initial_par,
    y_pre             = vc_pre, 
    y_post_mean       = mean(vc_post), 
    controls_pre      = controls_pre,
    treated_pre       = NA,  # no longer used inside
    sum_to_one_lambda = list_of_lambdas,
    lb                = lb,
    ub                = ub
  )
  
  result_tbl %>%
    mutate(within_tolerance = between(sum_par, 1 - tolerance, 1 + tolerance))
}
