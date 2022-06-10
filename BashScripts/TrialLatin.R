
TrialName="TrialLatin"

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

RunID =  as.numeric(Sys.getenv("SGE_TASK_ID"))

InputParameters<- read_csv(paste0('../Parameters/',TrialName, '.csv'),
                           show_col_types = FALSE)

print(t(InputParameters[RunID,]))

BuildMC6(RunID, InputParameters)  ##   # can comment out if already done

Sys.time()

Fit_VP10(RunID,TrialName,InputParameters,niter = 20000, nburn = 10000, thin = 10,OccuMin = 5 )  # full 
#Fit_VP7(RunID,TrialName,InputParameters,niter = 100, nburn = 50, thin = 1) # quick for testing


print('Finishing at:')
Sys.time()

print('Total time:')
Sys.time() - Start