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

Multi_SDID_Poisson_VC <- function(d, treated_list, n_starts, 
                                  lambda_list = c(-3,-1,0, 1,3,5)) {
  
  all_preds <- list()
  
  for (unit_id in treated_list) {
    d2 <- subset(d, unit_id %notin% treated_list | unit_id == unit_id)
    
    d3 <- RelevelPoisson(d = d2, treated_unit_id = unit_id)
    
    preds <- Poisson_SDID3_VC(d = d3, n_starts = n_starts, 
                              lambda_list = lambda_list)
    
    all_preds[[as.character(unit_id)]] <- preds
  }
  
  # Combine into a matrix (each column is a unit's prediction vector)
  pred_mat <- do.call(cbind, all_preds)
  
  # Sum across treated units (column-wise sum)
  summed_preds <- rowSums(pred_mat)
  
  # Optionally: get actual summed outcome for comparison
  treated_data <- subset(d, unit_id %in% treated_list)
  actuals <- treated_data |>
    group_by(time) |>
    summarise(sum_poisson = sum(actual.poisson))|>
    ungroup()
  
  tbl <- tibble(
    #time = 1:length(summed_preds),
    time = unique(d$time), 
    predicted_sum = summed_preds,
    actual_sum = actuals$sum_poisson, 
    residuals = actual_sum - predicted_sum
  )
  
  # Return both predictions and actuals
  return(tbl)
}
# 
# d <- Create_SC_Data(units_n = 50, 
#                     treated_units = 2, 
#                     impact_factor = 0)
# 
# Multi_SDID_Poisson_VC(d = d,
#                     treated_list = c(1,2),
#                     n_starts = 1) |>
#    ggplot(mapping = aes(x = 1:50)) +
#    geom_line(mapping = aes(y = predicted_sum),
#              col = "blue") +
#    geom_line(mapping = aes(y = actual_sum))
