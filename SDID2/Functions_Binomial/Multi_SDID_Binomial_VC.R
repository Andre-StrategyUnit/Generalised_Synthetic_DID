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

# d <- Create_SC_Data_Binomial(treated_units = 2)
Multi_SDID_Binomial_VC <- function(d, treated_list, n_starts) {
  
  all_preds <- list()
  
  # single imputation for all of the cases with zero trials
  d <- impute_dataset(d)
  
  for (i in treated_list) {
    d2 <- subset(d, unit_id %notin% treated_list | unit_id == i)
    
    d3 <- RelevelBinomial(d = d2, treated_unit_id = i)
    
    preds <- Binomial_SDID3_VC(d = d3, n_starts = n_starts)
    
    all_preds[[as.character(i)]] <- preds
  }
  
  # Combine into a matrix (each column is a unit's prediction vector)
  pred_mat <- do.call(cbind, all_preds)
  
  # Sum across treated units (column-wise sum)
  summed_preds <- rowSums(pred_mat)
  
  # Optionally: get actual summed outcome for comparison
  treated_data <- subset(d, unit_id %in% treated_list)
  actuals <- treated_data |>
    group_by(time) |>
    summarise(sum_binomial = sum(actual.binomial), 
              sum_trials = sum(trials))|>
    ungroup()
  
  tbl <- tibble(
    time = 1:length(summed_preds),
    predicted_sum = summed_preds,
    actual_sum = actuals$sum_binomial,
    trials_sum = actuals$sum_trials
  )
  
  # Return both predictions and actuals
  return(tbl)
}


# Multi_SDID_Binomial_VC(d = d,
#                    treated_list = c(1,2),
#                    n_starts = 1) |>
#   ggplot(mapping = aes(x = 1:max(d$time))) +
#   geom_line(mapping = aes(y = predicted_sum/trials_sum),
#             col = "blue") +
#   geom_line(mapping = aes(y = actual_sum/trials_sum))

