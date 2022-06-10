### Still in progress - adapting 
## variPart() to be more generic and convertable to other JSDM approaches


# ## Fixed options: 
# 
# Type = 'III'
# R2adjust = TRUE
# indSite = FALSE
# verbose = TRUE
# HMSCprior = NULL
# 
# 
# ## need to be provided
# 
# model= hmsc= model ## fitted full hmsc model object
# groupX =  c(rep("env",6),rep("spa",10))  # specification of categories of the predictors into envrionmental and spatial
# # ...  = other options passed to the HMSC function
# 

# I make use of purrr package functions to simplify the list operations

### I've also taken out the option include autocorrelated random effects. 

## I have direct specification of MCMC parametetrs 

## results are returned in a simpler dataframe format, not a nested list 

variPart_simple <- function (model,
                             nEnv, ## number of environmental variables (including intercept!!)
                             nSpa, niter=4000, nburn=2000, thin=10) {
  
#  nsite <- nrow(model$data$Y)
#  nsp <- ncol(model$data$Y)
  nX <- ncol(model$data$X)  ## number of predictors, including intercept term 
  family <- attributes(model)$class[2]  # model family (probit)
  
  ### build models
  env_data = model$data$X[,1:nEnv  ]
  spa_data = model$data$X[, (nEnv+1):nX]
  ran_data = model$data$Random
  
  ## need to add intercept only where not already part of environmental data
  HMSC_data_list<-  list( 'env' = as.HMSCdata(Y = model$data$Y, X = env_data, Random = NULL, interceptX = FALSE),  # 
                          'spa' = as.HMSCdata(Y = model$data$Y, X = spa_data, Random = NULL, interceptX = TRUE),   # 
                          'ran' = as.HMSCdata(Y = model$data$Y, X = NULL, Random = ran_data, interceptX = TRUE),   # 
                          'env_spa' = as.HMSCdata(Y = model$data$Y, X = model$data$X, Random = NULL, interceptX = FALSE),  # 
                          'env_ran' = as.HMSCdata(Y = model$data$Y, X = env_data, Random = ran_data, interceptX = FALSE),  # 
                          'spa_ran' = as.HMSCdata(Y = model$data$Y, X = spa_data, Random = ran_data, interceptX = TRUE)  # 
  )
  
  submodels_list <- map(.x = HMSC_data_list, 
                        .f=hmsc,
                        family = family, niter = niter, nburn = nburn, # ..., 
                        thin = thin, verbose = FALSE)
  
  ## find r2 for each sub model
  r2_dataframe <- map_dfc(submodels_list, R2_simplified,type = "nakagawa",  adjust = TRUE)
  
  ### add r_2 of the full model
  r2_dataframe$all <-R2_simplified(model, type = "nakagawa", adjust = TRUE)
  
  VennFrac <- data.frame(species = colnames(model$data$Y))                                        
  VennFrac$'random' = r2_dataframe$all - r2_dataframe$env_spa  #c                                
  VennFrac$'env'    = r2_dataframe$all - r2_dataframe$spa_ran     #a                             
  VennFrac$'spa'    = r2_dataframe$all - r2_dataframe$env_ran       #b                           
  VennFrac$'spa-cod'= r2_dataframe$all - r2_dataframe$env - (VennFrac$'random'+VennFrac$'spa')  #e
  VennFrac$'env-cod'= r2_dataframe$all - r2_dataframe$spa - (VennFrac$'random'+VennFrac$'env')  # f
  VennFrac$'env-spa'= r2_dataframe$all - r2_dataframe$ran - (VennFrac$'spa'+VennFrac$'env'  )   # d
  VennFrac$'env-spa-cod' = r2_dataframe$all - rowSums(VennFrac[,-1])          #g                 
  
  ## 
  Partition = data.frame(species = colnames(model$data$Y) )
  Partition$env =   VennFrac$'env' + VennFrac$'env-cod' + 0.5 * (VennFrac$'env-spa' +  VennFrac$'env-spa-cod' )
  Partition$spa =  VennFrac$'spa'  +  VennFrac$'spa-cod' + 0.5 * ( VennFrac$'env-spa' +  VennFrac$'env-spa-cod')
  Partition$codist = VennFrac$'random'
  
  Partition$env    = ifelse(Partition$env < 0,    0, Partition$env)
  Partition$spa    = ifelse(Partition$spa < 0,    0, Partition$spa)
  Partition$codist = ifelse(Partition$codist < 0, 0, Partition$codist)
  Partition$R2_sum = Partition$env +Partition$spa   +Partition$codist 
  Partition$R2_full =  r2_dataframe$all
  
  
  return(list(VennFrac=as_tibble(VennFrac), 
              r2_dataframe=as_tibble(r2_dataframe),
              Partition = as_tibble(Partition)))
}







