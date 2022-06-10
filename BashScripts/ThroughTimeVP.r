
TrialName="TrackThroughTime"

.libPaths("~/LibraryCentos/")

sessionInfo()

Start<- Sys.time()
print('Starting at:')
Start

## Load packages
library(tidyverse)
library(igraph)
library(HMSC)
select <- dplyr::select

cat(getwd())
setwd('~/GitHub/MC_coherence/BashScripts/')

### Load all functions
walk( list.files('../Scripts/', full.names = TRUE), source)

RunID =  1

t_test = as.numeric(Sys.getenv("SGE_TASK_ID"))

InputParameters<- read_csv(paste0('../Parameters/',TrialName, '.csv'),
                           show_col_types = FALSE)

load(file = paste0('../SavedObjects/', InputParameters$RunName[1], '/Run_',RunID,'_ArrayOutput'))

for(Occ_thresh in c(0.01, 0.1, 1)  ){
  
  #########
  # Processing to get focal species, binareised from focal nodes
  Side_length = Trial$Region_Length 
  
  Mat_to_keep <- matrix(TRUE, nrow = Side_length, ncol =Side_length)
  
  Mat_to_keep[c(1,2,Side_length-1, Side_length), ]<-FALSE
  Mat_to_keep[,c(1,2,Side_length-1, Side_length) ]<-FALSE
  
  sites_to_keep <-which(Mat_to_keep )
  
  
  slice <- Trial$ArrayContainer[,,t_test]
  slice_focalsites <- slice[,sites_to_keep]  # originally rows = species, cols = sites 
  
  ## Binarising
  slice_focalsites <- slice_focalsites> Occ_thresh
  
  ### Then find species that can be straight dropped if they too rarely occur 
  
  focal_sites_allslices <- Trial$ArrayContainer[,sites_to_keep,]
  
  
  ### need to do some more calculations about which species to keep - need to occur at at least 20 sites accross all time samples
  ## will be a bit lower as more times
  
  Occupancy_allslice <- apply(focal_sites_allslices>Occ_thresh, MARGIN = c(1,3), sum)  # make a species x time slice occupancy 
  KeepSpecies<- apply(Occupancy_allslice>=20, MARGIN = 1, all)
  
  ## also do a final transpose so have sites as rows, and species as columns
  slice_focalsites_focalspecies_t <- t(slice_focalsites[KeepSpecies, ])*1
  
  ##  Adding species names for back tracking if necessary
  speciesnames <- paste0('Sp_', which(KeepSpecies))
  
  colnames(slice_focalsites_focalspecies_t) <- speciesnames
  
  ## Adding node names to match original set up 
  nodenames <- paste0('Node_', sites_to_keep)
  rownames(slice_focalsites_focalspecies_t ) <- nodenames
  
  E_toUse= Trial$ModelSetUp$E_nodes_2[sites_to_keep]
  
  # End Processing
  #########
  
  ## Prepping envrionment and spatial data for fitting:
  
  MEM <- read.csv(paste0('../SavedObjects/MEMs/MEMof',Side_length,'side_focal',Side_length-4,'.csv'))
  X_formatted <- cbind(scale(E_toUse),
                       scale(E_toUse)^2,
                       MEM[,-1])  # description of sites . direct and squared environment + a moran eigenvector term (dropping node id numbers)
  
  X_groups <- c(rep("env",3),   # fitting three direct environmental terms (i.e. including intercept)
                rep("spa",ncol(MEM)-1))  ## and as many spatial terms as the eigen vector map does (NB name column dropped also)
  
  N_Sites<- (Side_length-4)^2
  
  ### Function to prepare all the data for HMSC function 
  formData_slice <- as.HMSCdata(Y = slice_focalsites_focalspecies_t,                   
                                X = X_formatted,  Random = as.factor(1:N_Sites),         
                                scaleX = TRUE,  interceptX = TRUE)
  
  
  ## fitting models
  model_slice<-  hmsc(formData_slice, family = "probit",
                      niter = 20000, nburn = 5000, thin = 10 )  # takes a few minutes
  
  ## Variance partitioning
  vpSpp_slice <-   variPart(model_slice,groupX = X_groups, type = "III",  R2adjust = TRUE)  # takes 6x as long... 
  
  ## Data reformatting and partitioning 
  vpSpp_slice %>% 
    map(as_tibble) %>%
    bind_cols() %>% 
    mutate(Sample = t_test,
           species =colnames(slice_focalsites_focalspecies_t) ) %>% 
    left_join(data.frame(prevalence = colSums(slice_focalsites_focalspecies_t), 
                         species  = colnames(slice_focalsites_focalspecies_t)   ) ,
              by = "species")  %>%
    set_names(c( "c", "b", "a", "e", "f", "d", "g", "Sample", "species", "prevalance")) %>% 
    transmute(species = species,
              Sample = Sample,
              env = a + f + 0.5 * d + 0.5 * g,
              env = ifelse(env < 0, 0, env),
              spa = b + e + 0.5 * d + 0.5 * g,
              spa = ifelse(spa < 0, 0, spa),
              codist = c,
              codist = ifelse(codist < 0, 0, codist),
              r2 = env + spa + codist,
              iteration = paste0('Run_', RunID)) -> VP_Both
  
  VP_Both$Occ_thresh <- Occ_thresh
  
  write_csv(VP_Both, file = paste0('../SavedObjects/', TrialName, '_VP/TimeSlice_',t_test,'_VP_Occ_', Occ_thresh))
  cat(paste0('\nDone VPing:TimeSlice_',t_test,
             'OccupancyThresh:', Occ_thresh,
             'Time:', Sys.time()))
}



print('Finishing at:')
Sys.time()

print('Total time:')
Sys.time() - Start


