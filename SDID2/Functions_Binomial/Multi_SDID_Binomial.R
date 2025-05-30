# d <- Create_SC_Data_Binomial(treated_units = 2,
#                              time_n = 50,
#                              units_n = 12,
#                              trials_lambda = 50)

Multi_SDID_Binomial <- function(d, treated_list, n_starts){

  # tictoc::tic()
  
  results <- tibble(
  cf.est = numeric(), 
  obs = numeric(), 
  vc_post_avg = numeric(), 
  treated_post_trials.sum = numeric(),
  binomial_sdid_att = numeric(), 
  post_timepoint = numeric())

  treated_list <- treated_list

  for (i in treated_list){
    
    d2 <- subset(d, 
                 unit_id %notin% treated_list | unit_id == treated_list[i])
    d3 <- RelevelBinomial(d = d2, 
                          treated_unit_id = treated_list[i])
    results[i, ] <- Binomial_SDID3(d = d3, 
                                   n_starts = n_starts)
    
  }
  # tictoc::toc()
  
  results <- results |>
    mutate(treated_unit = treated_list)
  
  return(results)
  
}

# Multi_SDID_Binomial(d, c(1,2), n_starts = 1)

