Aggregate_Placebos <- function(placebos, 
                               iter=1000){
  
  output2 <- tibble(
    cf = numeric(), 
    obs = numeric(), 
    impact = numeric())
  
  for (i in 1:iter){
    
    output2[i,] <- placebos |> 
      group_by(equiv_unit) |> 
      sample_n(size = 1) |> 
      ungroup() |> 
      summarise(cf = sum(est_cf_outcome * post_timepoints), 
                obs = sum(obs_outcome * post_timepoints), 
                impact = log(obs/cf))
  }
  
  return(output2)
  
}


Aggregate_Central <- function(output){
  
  k <- output |> 
    group_by(treated_unit) |> 
    summarise(obs_outcome = obs_outcome * post_timepoints, 
              cf_outcome = est_cf_outcome * post_timepoints, 
              timepoints = post_timepoints) |> 
    ungroup() |> 
    summarise(obs_outcome = sum(obs_outcome), 
              cf_outcome = sum(cf_outcome), 
              unit_timepoints = sum(timepoints), 
              impact_effect_size = log(obs_outcome/cf_outcome))
  
  return(k)
  
}

SDID_Inferential <- function(central_overview, 
                             aggregate_placebos, 
                             p_direction,
                             quantiles = c(0.95)){
  
  inferential <- NULL 
  inferential$central <- central_overview |> 
    mutate(bias = mean(output2$impact), 
           Unbiased_Impact = impact_effect_size - bias, 
           Unbiased_CF_Outcome = exp(log(obs_outcome) - Unbiased_Impact))
  
  confidence_intervals_work <- tibble(
    Quantiles = quantiles,
    Orig_Impact_Effect = as.numeric(quantile(x = output2$impact, 
                                             probs = quantiles)),
    Unbiased_Impact_Effect = as.numeric(quantile(x = output2$impact - mean(output2$impact), 
                                                 probs = quantiles)), 
    Orig_CF_Outcome = exp(log(central_overview$obs_outcome) - Orig_Impact_Effect), 
    Unbiased_CF_Outcome = exp(log(central_overview$obs_outcome) - Unbiased_Impact_Effect)
  )
  
  inferential$confidence_intervals <- tibble(
    Quantiles = quantiles,
    Orig_Impact_Effect = inferential$central$impact_effect_size  + 
                         as.numeric(quantile(x = output2$impact, 
                                             probs = quantiles)),
    Unbiased_Impact_Effect = inferential$central$Unbiased_Impact + 
                             as.numeric(quantile(x = output2$impact - mean(output2$impact), 
                                                 probs = quantiles)), 
    Orig_CF_Outcome = exp(log(central_overview$obs_outcome) - Orig_Impact_Effect), 
    Unbiased_CF_Outcome = exp(log(central_overview$obs_outcome) - Unbiased_Impact_Effect)
  )
  
  adj_placebos <- c(output2$impact - mean(output2$impact))
  
  p_table <- tibble(
    p_direction = p_direction, 
    value = case_when(
      p_direction == "higher" ~ paste0("the observed impact effect size was more positive than ", 
                                       100 * mean(inferential$central$Unbiased_Impact > adj_placebos), 
                                       "% of the placebo draws"), 
      p_direction == "lower" ~ paste0("the observed impact effect size was more negative than ", 
                                      100 * mean(inferential$central$Unbiased_Impact < adj_placebos), 
                                      "% of the placebo draws"),  
      p_direction == "two_sided" ~ paste0("the observed impact effect size was larger in absolute value than ", 
                                          100 * mean(abs(inferential$central$Unbiased_Impact) > abs(adj_placebos)), 
                                          "% of the placebo draws"), 
      .default = "no p value direction stated"
      
    )
  )
  
  inferential$p_value = p_table$value
  
  return(inferential)
  
}
