# source("Functions/Poisson_Functions.R")

Improved_Placebo <- function(d, main_result){
  
  # identify treated unit 
  treated_units <- unique((subset(d, treated == 1)$unit_id))

  # Placebo Inference
  starts_lambdas <- tibble(
    equiv_unit = treated_units,
    starts = 1 + max(d$time) - main_result$post_timepoints, 
    lambdas = log10(main_result$lambda_used))

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
    vc_mean = numeric(), 
    poisson_sdid_att = numeric(), 
    post_timepoints = numeric(), 
    pre_obs = numeric(), 
    pre_vc = numeric(),
    unit = numeric(), 
    start_date = numeric(),
    equiv_unit = numeric() )
  
  
  for (j in 1:nrow(starts_lambdas)){
  
    equiv_unit = starts_lambdas$equiv_unit[j]
    start_date_here = starts_lambdas$starts[j]
    lambda = starts_lambdas$lambdas[j]
  
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
    output <- Poisson_SDID4(d = d3, 
                            n_starts = n_starts, 
                            lambda_list = c(lambda)) |> 
      mutate(unit = i, 
             start_date = start_date_here, 
             equiv_unit = equiv_unit)

    
    list_results <-  rbind(list_results, output)
  }
  
  }
  
  return(list_results)
  
}

Improved_Placebo2 <- function(d, main_result, 
                              cum_contr = 0.9){
  
  # identify treated unit 
  treated_units <- unique((subset(d, treated == 1)$unit_id))
  
  # Placebo Inference
  starts_lambdas <- tibble(
    equiv_unit = treated_units,
    starts = 1 + max(d$time) - main_result$post_timepoints, 
    lambdas = log10(main_result$lambda_used))
  
  # select start dates
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
    vc_mean = numeric(), 
    poisson_sdid_att = numeric(), 
    post_timepoints = numeric(), 
    pre_obs = numeric(), 
    pre_vc = numeric(),
    unit = numeric(), 
    start_date = numeric(),
    equiv_unit = numeric() )
  
  
  for (j in 1:nrow(starts_lambdas)){
    
    equiv_unit = starts_lambdas$equiv_unit[j]
    start_date_here = starts_lambdas$starts[j]
    lambda = starts_lambdas$lambdas[j]
    
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
      output <- Poisson_SDID4(d = d3, 
                              n_starts = n_starts, 
                              lambda_list = c(lambda)) |> 
        mutate(unit = i, 
               start_date = start_date_here, 
               equiv_unit = equiv_unit)
      
      
      list_results <-  rbind(list_results, output)
    }
    
  }
  
  # next task is to exclude all the control units which didnt contribute to the vc 
  # so first identify all the controls which didnt contribute to each treated units vc
  
  treated_list <- treated_units
  
  placebos_keep <- tibble(
    treated_unit = numeric(), 
    ctrl_unit = numeric(), 
    contribution_pcnt = numeric()
  )

  for (i in treated_list){
    
    d2 <- subset(d, 
                 unit_id %notin% treated_list | unit_id == treated_list[i])
    d3 <- RelevelPoisson(d = d2, 
                         treated_unit_id = treated_list[i])
    
    treated_time <- min(subset(d$time, d$treated == 1))
    # Raw controls outcomes (counts)
    controls_pre <- Extract_Controls_Pre(d = d3, 
                                         treated_unit_id = treated_list[i], 
                                         treated_time = treated_time) |> 
      colSums()
    
    solution <- Poisson_SDID_solution(d = d3, 
                          n_starts = 1, 
                          lambda_list = starts_lambdas$lambdas[i])
    
    s_total = sum(controls_pre * solution[-1]) |> 
      as.numeric()
    
    
    a <- tibble(
      treated_unit = i, 
      ctrl_unit = max(treated_list) + (1:length(controls_pre)), 
      unit_pre_sum = controls_pre, 
      solution_weight = solution[-1], 
      predicted_contribution = unit_pre_sum * solution_weight, 
      contribution_pcnt = predicted_contribution/s_total) |> 
      arrange(desc(contribution_pcnt)) |> 
      select(treated_unit, ctrl_unit, contribution_pcnt) |> 
      mutate(cumulative = cumsum(contribution_pcnt), 
             keep = cumulative < cum_contr) |> 
      subset(keep == TRUE) |>
      select(treated_unit, ctrl_unit, contribution_pcnt)
    
    placebos_keep <- rbind(placebos_keep, a)
  
    }
  
  list_results2 <- left_join(
    x = list_results, 
    y = placebos_keep, 
    by = c("equiv_unit" = "treated_unit", "unit" = "ctrl_unit")) |> 
    subset(is.na(contribution_pcnt) == FALSE)
  
  return(list_results2)
  
}
  

  

