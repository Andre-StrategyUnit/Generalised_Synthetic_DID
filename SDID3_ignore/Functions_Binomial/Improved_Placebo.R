

Improved_Placebo.bn <- function(d){

# select staret dates
start_dates <- as.integer(unique(d$start_date))
start_dates <- subset(start_dates, is.na(start_dates) == FALSE)
  
  # create lookup
  d <- d
  # treated_list <- treated_list
  treated_list = unique(subset(d, treat_group == 1)$unit_id)
  #n_starts <- n_starts 
  n_starts = 1
  
  all_units <- unique(d$unit_id)
  ctrl_units <- setdiff(all_units, treated_list)
  
  
  list_results <- tibble(
    cf.est = numeric(), 
    obs = numeric(), 
    vc_post_avg = numeric(), 
    treated_post_trials.sum = numeric(), 
    binomial_sdid_att = numeric(), 
    post_timepoints = numeric(), 
    unit = numeric(), 
    start_date = numeric())
  
  
  for (j in start_dates){
  
    start_date_here = j
  
  for (i in ctrl_units){
    
    d2 <- subset(d, unit_id %in% ctrl_units) |>
      mutate(treat_group = ifelse(test = unit_id == i, 
                                  yes = 1, 
                                  no = unit_id), 
             unit_id = ifelse(test = unit_id == i,
                              yes = 1, 
                              no = unit_id)) 
    d3 <- d2|> 
      mutate(start_date = ifelse(test = unit_id == 1,
                                 yes = start_date_here, 
                                 no = NA),
             treated =  ifelse(test = unit_id == 1 & time >= start_date_here, 
                               yes = 1, 
                               no = 0))
    
    # run analysis
    output <- Multi_SDID_Binomial(d = d3, 
                                  treated_list = 1, 
                                  n_starts = n_starts) |> 
      mutate(unit = i, 
             start_date = start_date_here)

    
    list_results <-  rbind(list_results, output)
  }
  
  }
  
  return(list_results)
  
}
  

  

