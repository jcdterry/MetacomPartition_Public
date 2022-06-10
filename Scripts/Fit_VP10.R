Fit_VP10 <- function(RunID,TrialName,InputParameters, niter, nburn, thin, OccuMin =5, FalseAbsence=0 ){
  

  load(file = paste0('../SavedObjects/', InputParameters$RunName[1], '/Run_',RunID,'_MCobject_PreProcess'))
  
  #########
  # Processing to get focal species, binareised from focal nodes
  Side_length = Trial$Region_Length 
  
  Mat_to_keep <- matrix(TRUE, nrow = Side_length, ncol =Side_length)
  
  Mat_to_keep[c(1,2,Side_length-1, Side_length), ]<-FALSE
  Mat_to_keep[,c(1,2,Side_length-1, Side_length) ]<-FALSE
  
  sites_to_keep <-which(Mat_to_keep )
  
  Before_focalsites <- Trial$Before_Matrix[,sites_to_keep]  # originally rows = species, cols = sites 
  During_focalsites <- Trial$During_Matrix[,sites_to_keep]  # originally rows = species, cols = sites 
  
  ## Binarising
  Occ_thresh = InputParameters$Occ_thresh[RunID]
  Before_focalsites <- Before_focalsites> Occ_thresh
  During_focalsites <- During_focalsites> Occ_thresh

  ## Introducing False Absences [False Absence is a number between 0 and 1, of the number of real observations that are missed]
  if(FalseAbsence!= 0 ){
    BeforeCorrectlyObserved <- rbernoulli(n = nrow(Before_focalsites)*ncol(Before_focalsites),  p = 1-FalseAbsence)
    DuringCorrectlyObserved <- rbernoulli(n = nrow(During_focalsites)*ncol(During_focalsites),  p = 1-FalseAbsence)
    
    
    cat('True n_obs (before:)')
    cat(sum(Before_focalsites))
    
    Before_focalsites <- Before_focalsites*BeforeCorrectlyObserved
    
    
    cat('Observed n_obs (before:)')
    cat(sum(Before_focalsites))
    
    During_focalsites <- During_focalsites*DuringCorrectlyObserved
  }
  
  
  
  ### Then find species that can be straight dropped if they too rarely occur 
  Occupancy_before <- rowSums(Before_focalsites)
  Occupancy_during <- rowSums(During_focalsites)
  
  ## Not too rare or too common
  OccuMax = 100-OccuMin
  
  KeepSpecies <- Occupancy_before>=OccuMin & Occupancy_during>=OccuMin &  Occupancy_before <= OccuMax & Occupancy_during<=OccuMax 
  
  ## also do a final transpose so have sites as rows, and species as columns
  Before_focalsites_focalspecies_t <- t(Before_focalsites[KeepSpecies, ])*1
  During_focalsites_focalspecies_t <- t(During_focalsites[KeepSpecies, ])*1
  
  ##  Adding species names for back tracking if necessary
  speciesnames <- paste0('Sp_', which(KeepSpecies))
  
  colnames(Before_focalsites_focalspecies_t) <- speciesnames
  colnames(During_focalsites_focalspecies_t) <- speciesnames
  
  ## Adding node names to match original set up 
  nodenames <- paste0('Node_', sites_to_keep)
  rownames(Before_focalsites_focalspecies_t ) <- nodenames
  rownames(During_focalsites_focalspecies_t ) <- nodenames
  
  E_toUse= Trial$ModelSetUp$E_nodes_2[sites_to_keep]

  ChangingEVar = Trial$ModelSetUp$E_nodes[sites_to_keep]
  ### because changing var changes equally for all species and all sites, essentially just offset,
  ## and effect of CC would be lost in the scaling
  
  
  # End Processing
  #########
  
  ## Prepping envrionment and spatial data for fitting:
  
  MEM <- read.csv(paste0('../SavedObjects/MEMs/MEMof',Side_length,'side_focal',Side_length-4,'.csv'))
  X_formatted <- cbind(scale(E_toUse),
                       scale(E_toUse)^2,
                       scale(ChangingEVar),  # now with the 'climate' variavle (not including noise terms)
                       scale(ChangingEVar)^2,
                       MEM[,-1])  # description of sites . direct and squared environment + a moran eigenvector term (dropping node id numbers)
  
  X_groups <- c(rep("env",5),   # fitting five direct environmental terms (i.e. including intercept)
                rep("spa",ncol(MEM)-1))  ## and as many spatial terms as the eigen vector map does (NB name column dropped also)
  
  N_Sites<- (Side_length-4)^2
  
  
  ### Function to prepare all the data for HMSC function 
  formData_before <- as.HMSCdata(Y = Before_focalsites_focalspecies_t,                   
                                 X = X_formatted,  Random = as.factor(1:N_Sites),         
                                 scaleX = TRUE,  interceptX = TRUE)
  
  formData_during <- as.HMSCdata(Y = During_focalsites_focalspecies_t,                           
                                 X =X_formatted,  Random = as.factor(1:N_Sites),           
                                 scaleX = TRUE,  interceptX = TRUE) 
  
  ## fitting models
  model_before<-  hmsc(formData_before, family = "probit",
                       niter = niter, nburn = nburn, thin = thin )  # takes a few minutes
  
  model_during<-  hmsc(formData_during, family = "probit",
                       niter = niter, nburn = nburn, thin = thin  )  # takes a few minutes
  
  ## Variance partitioning
  vpSpp_before <-   variPart(model_before,groupX = X_groups, type = "III",  R2adjust = TRUE)  # takes 6x as long... 
  
  vpSpp_during <-   variPart(model_during, groupX = X_groups, type = "III",  R2adjust = TRUE)  # takes 6x as long... 
  
  ## Data reformatting and partitioning 
  vpSpp_before %>% 
    map(as_tibble) %>%
    bind_cols() %>% 
    mutate(Sample = 'Before_WithTempVar',
           species =colnames(Before_focalsites_focalspecies_t) ) %>% 
    left_join(data.frame(prevalence = colSums(Before_focalsites_focalspecies_t), 
                         species  = colnames(Before_focalsites_focalspecies_t)   ) ,
              by = "species") -> SPECIES_VP_DATA_before
  
  vpSpp_during %>% 
    map(as_tibble) %>%
    bind_cols() %>% 
    mutate(Sample = 'During_WithTempVar',
           species =colnames(During_focalsites_focalspecies_t) ) %>% 
    left_join(data.frame(prevalence = colSums(During_focalsites_focalspecies_t), 
                         species  = colnames(During_focalsites_focalspecies_t)   )    ,
              by = "species") -> SPECIES_VP_DATA_during
  
  # bind_rows(SPECIES_VP_DATA_before,SPECIES_VP_DATA_during ) %>%
  #   set_names(c( "c", "b", "a", "e", "f", "d", "g", "Sample", "species", "prevalance")) %>% 
  #   transmute(species = species,
  #             Sample = Sample,
  #             env = a + f + 0.5 * d + 0.5 * g,
  #             env = ifelse(env < 0, 0, env),
  #             spa = b + e + 0.5 * d + 0.5 * g,
  #             spa = ifelse(spa < 0, 0, spa),
  #             codist = c,
  #             codist = ifelse(codist < 0, 0, codist),
  #             r2 = env + spa + codist,
  #             iteration = paste0('Run_', RunID)) -> VP_Both
  VP_Both <-  bind_rows(SPECIES_VP_DATA_before,SPECIES_VP_DATA_during ) 
  VP_Both$Occ_thresh <- Occ_thresh
  VP_Both$RunID = RunID
  
  if(!dir.exists(paste0('../SavedObjects/', TrialName , '_VP'))){
    dir.create(paste0('../SavedObjects/', TrialName, '_VP') )
  }

  write_csv(VP_Both, file = paste0('../SavedObjects/', TrialName, '_VP/Run_',RunID,'.csv'))
  cat(paste0('\nDone VPing:Run_',RunID, 'Time:', Sys.time()))
  return( TRUE)
}
