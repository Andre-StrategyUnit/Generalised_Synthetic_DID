# Impute Zero Trials 

impute_unit <- function(df_unit) {
  times <- df_unit$time
  rates <- df_unit$rate
  successes <- df_unit$actual.binomial
  trials <- df_unit$trials
  
  for (i in seq_along(times)) {
    if (trials[i] == 0) {
      # Find indices of time points with non-zero trials
      nonzero <- which(trials > 0 & !is.na(rates))
      
      if (length(nonzero) == 0) {
        stop(paste0("Unit ", unique(df_unit$unit_id), " has no non-zero trials to impute from at time ", times[i]))
      }
      
      # Calculate absolute distance to current time point
      distances <- abs(times[nonzero] - times[i])
      
      # Find index of the closest time point(s)
      min_dist <- min(distances)
      closest_indices <- nonzero[distances == min_dist]
      
      # Handle ties by choosing the first closest (could be improved if needed)
      nearest_time <- closest_indices[1]
      nearest_rate <- rates[nearest_time]
      nearest_rate <- pmin(pmax(nearest_rate, 1e-8), 1 - 1e-8)  # clamp
      
      # Impute successes and trials
      trials[i] <- 1
      successes[i] <- rbinom(1, size = 1, prob = nearest_rate)
    }
  }
  
  df_unit$actual.binomial <- successes
  df_unit$trials <- trials
  return(df_unit)
}

impute_dataset <- function(d){
  
  d2 <- d|>
    mutate(rate = actual.binomial/trials)
  
  d_imputed <- d2 |>
    group_by(unit_id) |>
    group_split() |>
    map_df(impute_unit)

  d_imputed |>
    mutate(rate = actual.binomial/trials)
  
  return(d_imputed)
}
