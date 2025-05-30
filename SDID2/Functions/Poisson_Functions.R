# Import functions

source("Functions/Create_Data.R") # create dataset
source("Functions/Pre_Wrangle.R")
source("Functions/Gradient_Functions.R") # gradient functions
source("Functions/Time_Gradient_Functions.R") # sdid time gradient functions
source("Functions/Optimise_Unit_Lambdas.R") # optimise unit lambdas
source("Functions/Optimise_Time_Lambdas.R") # optimise time lambda
source("Functions/Multi_Start_Nloptr.R") # multi start process
source("Functions/Optimiser_Wrapper.R") # wrapper for extracting the optimised weights with multiple starts
source("Functions/Sum_To_One_Lambda.R") # try few lambdas to ensure time weights sum to and give performance of them
source("Functions/Post_Intervention_Mu.R") # calculate the post intervention mean using the time weights
source("Functions/Extract_Pre_Match.R")
source("Functions/Relevel_Poisson.R")
source("Functions/Poisson_SDID3.R")
`%notin%` <- Negate(`%in%`)
source("Functions/Multi_SDID_Poisson.R")
source("Functions/JackKnife_Poisson_SDID.R") 
source("Functions/Poisson_SDID_Inferential.R") 
source("Functions/Poisson_SDID3_VC.R") 
source("Functions/Multi_SDID_Poisson_VC.R") 
