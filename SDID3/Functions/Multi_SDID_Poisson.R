# d <- Create_SC_Data(treated_units = 0,
#                     scaling_factor = 25,
#                     time_n = 50,
#                     units_n = 12)
# d2 <- d |>
#   mutate(treated = case_when(
#     unit_id == 1 & time > 45 ~ 1,
#     unit_id == 1 & time <= 45 ~ 0,
#     unit_id == 2 & time > 25 ~ 1,
#     unit_id == 2 & time <= 25 ~ 0,
#     unit_id == 3 & time > 37 ~ 1,
#     unit_id == 3 & time <= 37 ~ 0,
#     .default = 0
#   ))

Multi_SDID_Poisson <- function(d, treated_list, 
                               n_starts, 
                               lambda_list = c(-3,-1,0, 1,3,5)){

  # tictoc::tic()
  
  results <- tibble(
    cf.est = numeric(), 
    obs = numeric(), 
    vc_post_avg = numeric(), 
    poisson_sdid_att = numeric(), 
    post_timepoints = numeric(), 
    pre_obs = numeric(), 
    pre_vc = numeric())

  treated_list <- treated_list

  for (i in treated_list){
    
    d2 <- subset(d, 
                 unit_id %notin% treated_list | unit_id == treated_list[i])
    d3 <- RelevelPoisson(d = d2, 
                          treated_unit_id = treated_list[i])
    results[i, ] <- Poisson_SDID3(d = d3, 
                                  n_starts = n_starts, 
                                  lambda_list = lambda_list)
    
  }
  # tictoc::toc()
  
  results <- results |>
    mutate(treated_unit = treated_list)
  
  return(results)
  
}

# Multi_SDID_Poisson(d = d2,
#                    treated_list = c(1,2), 
#                    n_starts = 1)

