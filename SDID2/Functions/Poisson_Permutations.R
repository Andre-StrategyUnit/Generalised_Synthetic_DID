
Poisson_Permutations <- function(d, treated_list, 
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
    obs_outcome_aggr = numeric(), 
    cf_outcome_aggr = numeric(), 
    scaled_effect = numeric(), 
    vc_outcome_aggr = numeric()
  )
  
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
    
    output <- Multi_SDID_Poisson(d = d2, 
                                 treated_list = treated_list, 
                                 n_starts = n_starts)
    
    central_estimate <- output |>
      mutate(obs_outcome_aggr = obs * post_timepoints, 
             cf_outcome_aggr = cf.est * post_timepoints, 
             vc_outcome_aggr = vc_post_avg * post_timepoints) |>
      summarise(obs_outcome_aggr = sum(obs_outcome_aggr), 
                cf_outcome_aggr = sum(cf_outcome_aggr), 
                scaled_effect = log(obs_outcome_aggr) - log(cf_outcome_aggr), 
                vc_outcome_aggr = sum(vc_outcome_aggr)) |> 
      ungroup()
    
    list_results[i, ] <- central_estimate
  }
  
  return(list_results)
  
}


Poisson_Permutations.single <- function(d, treated_list, 
                                        n_starts = 1){
  
  # create lookup
  d <- d
  treated_list <- treated_list
  n_starts <- n_starts 
  
  all_units <- d$unit_id
  ctrl_units <- setdiff(all_units, treated_list)
  
  placebo_groups <- split(
    sample(ctrl_units), 
    ceiling(seq_along(ctrl_units)/1)
  )
  
  
  start_date_tbl <- d |> 
    subset(unit_id %in% treated_list) |>
    group_by(unit_id) |>
    summarise(start_date2 = mean(start_date)) |>
    ungroup()
  
  list_results <- tibble(
    obs_outcome_aggr = numeric(), 
    cf_outcome_aggr = numeric(), 
    scaled_effect = numeric(), 
    synth_control_aggr = numeric()
  )
  
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
    
    output <- Multi_SDID_Poisson(d = d2, 
                                 treated_list = treated_list, 
                                 n_starts = n_starts)
    
    central_estimate <- output |>
      mutate(obs_outcome_aggr = obs * post_timepoints, 
             cf_outcome_aggr = cf.est * post_timepoints, 
             synth_control_aggr = vc_post_avg * post_timepoints) |>
      summarise(obs_outcome_aggr = sum(obs_outcome_aggr), 
                cf_outcome_aggr = sum(cf_outcome_aggr), 
                scaled_effect = log(obs_outcome_aggr) - log(cf_outcome_aggr), 
                synth_control_aggr = sum(synth_control_aggr)) |> 
      ungroup()
    
    list_results[i, ] <- central_estimate
  }
  
  return(list_results)
  
}