RelevelBinomial <- function(d, treated_unit_id){
  
  # Create a new vector of unique units, with the treated unit first
  unit_levels <- d$unit_id |> unique() |> setdiff(treated_unit_id) |> sort()
  new_levels <- c(treated_unit_id, unit_levels)
  
  # Create named lookup table: unit name → new ID
  lookup <- setNames(seq_along(new_levels), new_levels)
  
  # Apply the recoding
  d <- d |>
    mutate(unit_id = lookup[as.character(unit_id)])
  
  return(d)
}

# d <- Create_SC_Data_Binomial(units_n = 20,
#                              time_n = 100)
# d2 <- d|>
#   subset(unit_id != 1) |>
#   mutate(treated = ifelse(test = unit_id == 10 & time > 50,
#                           yes = 1,
#                           no = 0))
# 
# 
# d3 <- RelevelBinomial(d = d2,
#                      treated_unit_id = 10)
# Binomial_SDID3(d3)
