## Script to source for building MCs [version]

## this one just returns the raw data, and does not do the processing which is not handled seperately 

BuildMC6 <- function(RunID, InputParameters){
  
  ## Hardcoded parameters 
  num_blocks = 1000 ## total time
  Ext_Thresh = 0.00001 
  num_steps_per_block = 10 # number of time intervals to advance model per block
  ArrivalPop = 0.0001  ## initial density of newly arrived species (made a lower to account for wider initial appearance. )
  ImmFrac = 1/250 ## proportion of global species pool that arrive each step (nb ceilinged)
  K = 10 # carrying capacity parameter 
  Length_CC = 25
  Region_Length = 14
  
  Disp     <- InputParameters$Disp[RunID]
  comp_str <- InputParameters$comp_str[RunID]
  R_max    <- InputParameters$R_max[RunID]
  comp_conn    <- InputParameters$comp_conn[RunID]
  S_global <- InputParameters$S_global[RunID]
  EnvNoiseSigma <- InputParameters$EnvNoiseSigma[RunID]
  
  
  ######################################  V-  NB switch to random intiation here
  Before_During_Tester<- InitiateSystem_Rand( Pars = list( seed =  RunID ,
                                                      Region_Length = Region_Length, # Size of Square of nodes
                                                      S_global = S_global, # number of species in global pool
                                                      comp_conn = comp_conn,  #  competitive connectance
                                                      mutu_conn = NULL,  ## positive interaction 
                                                      mutu_relstr = 0))  
  
  B_MAT   <- Before_During_Tester$B_MAT 
  CompMat <- Before_During_Tester$CompMat
  DispMat <- Before_During_Tester$DispMat
  E_nodes <- Before_During_Tester$E_nodes
  EnvOpt  <- Before_During_Tester$EnvOpt
  E_nodes_2 <- Before_During_Tester$E_nodes_2
  EnvOpt_2  <- Before_During_Tester$EnvOpt_2
  
  R_Mat_2 <- Before_During_Tester$R_Mat_2
  R_Mat_3 <- Before_During_Tester$R_Mat_3
  
  
  ### Derived Parameters
  timescale = 1/ num_steps_per_block  ## ham fisted adjustor
  N_sp <- nrow( B_MAT)
  t_max = num_steps_per_block +1
  N_Nodes <-  ncol(B_MAT)
  
  ### Determining climate change
  E_change_noise =c( rep(0, num_blocks-Length_CC),  ## initial period of no noise to build community in 
                     seq(0, 0.25,length.out = Length_CC)) + 
    c( rep(0, 500), 
       rnorm(num_blocks-500, sd = EnvNoiseSigma))
  
  #### run model
  #####
  ## Output results container   array time / species / nodes 
  Total_species = rep(NA, num_blocks)
  
  for(t in 1:num_blocks){
    
    #########
    ### Update Environment + growth rates
    # R_Mat <- Calc_RMat(E_nodes, EnvOpt, E_change=E_change_noise[t], R_max = R_max)
    # R_Mat <- Calc_RMat_2(E_nodes, EnvOpt, E_nodes_2, EnvOpt_2,
    #                     RelEnvInfl = RelEnvInfl, E_change_noise[t], R_max = R_max)
    
    R_Mat <- Calc_RMat_3(E_nodes, EnvOpt, 
                         E_change=E_change_noise[t],
                         R_max = R_max,
                         R_Mat_2, R_Mat_3)
    
    ############
    ## Dispersal
    B_MAT <- B_MAT +  TotalDispersal(Disp, B_MAT, DispMat,ImmFrac, N_sp, N_Nodes , ArrivalPop)
    ############
    ## Run Dynamics
    B_MAT <- Dynamics_Discrete_Sparse(B_MAT , t_max, timescale,R_Mat,comp_str, CompMat, K)
    B_MAT <- as.matrix(B_MAT)
    #######
    ## Remove locally extinct species 
    ## this option does it at a node level (prevents 'tunnelling') : ## B_MAT[B_MAT< Parameters$Ext_Thresh] <- 0  
    B_MAT[rowSums(B_MAT>Ext_Thresh)==0, ] <- 0  ## This option removes species below threshold in all sites
    #################
    ## Tracking into containers , and just before climate change]
    ### save jsut before CC
    if(t == num_blocks){  During_Matrix<- B_MAT} # save last
    ## Save at end.
    if(t == num_blocks-Length_CC -1  ){   Before_Matrix<- B_MAT} # pre CC save 
    Total_species[t]<- sum(rowSums(B_MAT>Ext_Thresh)>0)
    
    if(t%%100==0){print(t)}
  }

  
  Trial <- list(  Name = paste0('Run',InputParameters$RunName[1],'_',RunID),
                  Before_Matrix=Before_Matrix,
                  During_Matrix= During_Matrix ,
                  ModelSetUp = Before_During_Tester,
                  E_change_noise = E_change_noise,
                  Total_species= Total_species,
                  Region_Length = Region_Length
  )
  
  if(!dir.exists(paste0('../SavedObjects/', InputParameters$RunName[1]))){
    dir.create(paste0('../SavedObjects/', InputParameters$RunName[1]) )
  }
  
  save(Trial, file = paste0('../SavedObjects/', InputParameters$RunName[1], '/Run_',RunID,'_MCobject_PreProcess'))
  cat(paste0('\nDone Building: Run_',RunID, 'Time', Sys.time()))
  
  return( TRUE)
}
