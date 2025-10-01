
Binomial_Permutations <- function(d, treated_list, 
                                  n_starts = 1){
  
  # create lookup
  d <- d
  treated_list <- treated_list
  n_starts <- n_starts 
  
  all_units <- d$unit_id
  ctrl_units <- setdiff(all_units, treated_list)
  
  placebo_groups <- split(
    sample(ctrl_units), 
    ceiling(seq_along(ctrl_units)/length(treated_list))
  )
  
  if(length(placebo_groups[[length(placebo_groups)]]) < length(treated_list)){
    
    placebo_groups[[length(placebo_groups)]] <- NULL
  }
  
  
  start_date_tbl <- d |> 
    subset(unit_id %in% treated_list) |>
    group_by(unit_id) |>
    summarise(start_date2 = mean(start_date)) |>
    ungroup()
  
  list_results <- tibble(
    obs_success = numeric(), 
    cf_success = numeric(), 
    trials = numeric(), 
    obs_odds = numeric(), 
    cf_odds = numeric(), 
    scaled_effect = numeric())
  
  for (i in 1:length(placebo_groups)){
    
    lookups <- tibble(
      fake_treated = placebo_groups[[i]], 
      unit_id = treated_list) 
    lookups2 <- full_join(start_date_tbl, lookups, by = "unit_id") |> 
      rename(acting_as = unit_id) |>
      rename(unit_id = fake_treated)
    
    
    # re wrangle the dataset
    
    ctrls_d <- d |> 
      subset(unit_id %notin% treated_list)
    d2 <- left_join(ctrls_d, lookups2, by = "unit_id") |> 
      select(-start_date) |>
      rename(start_date = start_date2) |>
      mutate(unit_id = ifelse(test = is.na(start_date), 
                              yes = unit_id, 
                              no = acting_as), 
             treat_group = ifelse(test = is.na(start_date),
                                  yes = 0, 
                                  no = 1),
             treated = case_when(
               is.na(start_date) ~ 0, 
               time < start_date ~ 0, 
               time >= start_date ~ 1)) |> 
      arrange(unit_id, time)
    
    # run analysis
    
    output <- Multi_SDID_Binomial(d = d2, 
                                 treated_list = treated_list, 
                                 n_starts = n_starts)
    
    central_estimate <- output |>
      mutate(obs_success = round(obs* treated_post_trials.sum, 0) , 
             cf_success = round(cf.est * treated_post_trials.sum, 0)) |>
      summarise(obs_success = sum(obs_success), 
                cf_success = sum(cf_success), 
                trials = sum(treated_post_trials.sum)) |>
      ungroup() |>
      mutate(cf_success = case_when(
        cf_success == 0 ~ 1, 
        cf_success == trials ~ trials-1, 
        .default = cf_success)) |>
      mutate(obs_odds = obs_success/(trials-obs_success), 
             cf_odds = cf_success/(trials-cf_success), 
             scaled_effect = log(obs_odds/cf_odds)) 
    
    list_results[i, ] <- central_estimate
  }
  
  return(list_results)
  
}

