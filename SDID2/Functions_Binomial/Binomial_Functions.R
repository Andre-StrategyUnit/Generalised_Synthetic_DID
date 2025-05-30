# import functions
source("Functions_Binomial/Create_Data_Binomial.R") # create dataset
source("Functions_Binomial/Pre_Wrangle_Binomial.R")
source("Functions_Binomial/Gradient_Functions_Binomial.R") # gradient functions
source("Functions_Binomial/Multi_Start_Nloptr_Binomial.R") # multi start process
source("Functions_Binomial/Optimiser_Wrapper_Binomial.R") # wrapper for extracting the optimised weights with multiple starts
source("Functions_Binomial/Affine_Time_Weights_Binomial.R") # calculate the post intervention mean using the time weights
source("Functions_Binomial/Relevel_Binomial.R") 
source("Functions_Binomial/Binomial_SDID3.R") 
`%notin%` <- Negate(`%in%`)
source("Functions_Binomial/Multi_SDID_Binomial.R") 
source("Functions_Binomial/JackKnife_Binomial_SDID.R") 
source("Functions_Binomial/Binomial_SDID_Inferential.R") 
source("Functions_Binomial/Binomial_SDID3_VC.R") 
source("Functions_Binomial/Multi_SDID_Binomial_VC.R") 


