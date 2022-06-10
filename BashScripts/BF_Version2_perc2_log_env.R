
.libPaths("~/LibraryCentos/")

### repeated using same code with more years 
## 4 May
## using different species and variable cuts after talk with Will
## removes burnet campion, log+1's Woodland, Calc Grass, Urban (not farming or temperature)
## Saves full models for fit examination
## Saving raw R2s, not the partitions

RUN_NAME <- 'HMSC_Bf2_perc2_log_env'

sessionInfo()

Start<- Sys.time()
print('Starting at:')
Start

## Load packages
library(tidyverse)
library(HMSC)
select <- dplyr::select

cat(getwd())
setwd('~/GitHub/MC_coherence/BashScripts/')

RunID =  as.numeric(Sys.getenv("SGE_TASK_ID"))

YearSetToTest <- 2000:2020
YearToTest <- YearSetToTest[RunID]

load('../../Butterflydata/YearSlices2000x2020_v2_2perc_nomoths')   ### using the more cut down data

SiteSelection <- read_csv('../../Butterflydata/Site_selection2_CHALK.csv')
MEM   <- read_csv('../../Butterflydata/Data forButterflies/MEM_justfocalsites_first10.csv')

colnames( Year_Slices[[1]])


## filter by year

Year_slice<-  arrange(Year_Slices[[RunID]],  OSGR)
Site_by_Sp_Mat <- as.matrix(Year_slice[,-1])

if(nrow(SiteSelection) != nrow(Year_slice)){stop('Number of Sites has got misaligned')}

## Environment Data
SiteSelection %>%  
  arrange(OSGR)-> EnvData

## Improving the distribution of three of the envrionmental variables
EnvData$Woodland <- log(EnvData$Woodland+1)
EnvData$Calc_Grass <- log(EnvData$Calc_Grass+1)
EnvData$Urban <- log(EnvData$Urban+1)


### Spatial Data
MEM%>%
  filter( OSGR %in% SiteSelection$OSGR)  %>%  # just in case
  arrange(OSGR)  -> MEMBySite

## Combining Predictors
X_formatted <- left_join(EnvData,
                         MEMBySite, by = "OSGR" ) %>%
  select(-OSGR) %>%
  as.matrix() 

### prepare all the data for HMSC function 
formData <- as.HMSCdata(Y = Site_by_Sp_Mat ,                   
                        X = X_formatted,
                        Random = as.factor(1:nrow(SiteSelection)),         
                        scaleX = TRUE, 
                        interceptX = TRUE)

set.seed(1)
## HMSC  Fitting:
model<-  hmsc(formData, family = "probit",
              niter = 100000, nburn = 50000, thin = 50 )  # takes a few minutes

if(!dir.exists(paste0('../SavedObjects/', RUN_NAME,'_fits'))){
  dir.create(paste0('../SavedObjects/', RUN_NAME,'_fits') )
}
save(model, file = paste0('../SavedObjects/',RUN_NAME,'_fits/model_',YearToTest))

## Second model for convergence checking 
set.seed(2)
model2<-  hmsc(formData, family = "probit",
               niter = 100000, nburn = 50000, thin = 50 )  # takes a few minutes

save(model2, file = paste0('../SavedObjects/',RUN_NAME,'_fits/model_',YearToTest, '_2'))

## Variance partitioning
X_groups <- c(rep("env",7),   # fitting six direct environmental terms + intercept
              rep("spa",10))  ## and num terms used from MEM

source( '../Scripts/variPart_returnboth.R')  ## variPart_returnboth

both_fracandr2<-variPart_returnboth(model,groupX = X_groups,
                                    type = "III",  R2adjust = TRUE) 

vpSpp <- both_fracandr2$res
vpSpp %>% 
  map(as_tibble) %>%
  bind_cols() %>%
  mutate(Sample = RUN_NAME,
         species  = colnames(Year_slice)[-1] ,
         Year= YearToTest) -> Venn

Venn %>%
  set_names(c( "c", "b", "a", "e", "f", "d", "g", "Sample", "species" , 'Year')) %>% 
  transmute(env = a + f + 0.5 * d + 0.5 * g,
            env = ifelse(env < 0, 0, env),
            spa = b + e + 0.5 * d + 0.5 * g,
            spa = ifelse(spa < 0, 0, spa),
            codist = c,
            codist = ifelse(codist < 0, 0, codist),
            r2 = env + spa + codist) -> vpDF

both_fracandr2$R2model %>% 
  map(as_tibble) %>%
  bind_cols() -> R2s


if(!dir.exists(paste0('../SavedObjects/', RUN_NAME))){
  dir.create(paste0('../SavedObjects/', RUN_NAME) )
}
bind_cols( Venn,vpDF ,R2s)%>%
  select(  Species = species,
           Year=Year, 
           seg_c= 1,
           seg_b = 2,
           seg_a = 3,
           seg_e=  4,
           seg_f= 5, 
           seg_d= 6, 
           seg_g =7 ,
           Part_env = 11,
           Part_spa = 12, 
           Part_codist = 13, 
           Part_R2_sum = r2, 
           Model_env = 15,
           Model_spa = 16,
           Model_random = 17,
           Model_EnvSpa = 18,
           Model_EnvRan = 19,
           Model_SpaRan= 20, 
           Model_All = 21) %>%
  write_csv(paste0('../SavedObjects/', RUN_NAME,'/vpDF_Y_', YearToTest, '.csv'))

print('Finishing at:')
Sys.time()

print('Total time:')
Sys.time() - Start
