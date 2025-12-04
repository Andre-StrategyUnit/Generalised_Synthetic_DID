# Example poisson run through

rm(list = ls())
a <- Sys.time()

# Simulate data ----
source("Functions/Poisson_Functions.R")
treated_units = c(1:5)
d <- Create_SC_Data(units_n = 50,
                    good_ctrls = 40,
                    time_n = 100, 
                    treated_time = 80,
                    treated_units = max(treated_units), 
                    scaling_factor = 3, 
                    impact_factor = 0.2)

# Run Analysis

# estimate impact for each treated unit
output <- Multi_SDID_Poisson(d = d,
                             treated_list = treated_units)

# aggregate estimates of impact
aggr_output <- Aggregate_Central_Poisson(output)

# placebo tests 
placebos <- Placebo_Poisson(d = d, 
                            main_result = output)
placebos_aggregated <- Aggregate_Placebos_Poisson(placebos = placebos)

# construct inferential statistics
inferential_results <- SDID_Inferential_Poisson(
  central_overview = aggr_output, 
  aggregate_placebos = placebos_aggregated)

