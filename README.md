# Code and data to support 'Codistribution as an indicator of whole metacommunity response to environmental change.' (Terry et al 2023, Ecography)
---

This public repository contains all the scripts behind the paper *Codistribution as an indicator of whole metacommunity response to environmental change* (Terry et al 2023, DOI: 10.1111/ecog.06605). Do contact me (Chris Terry, currently christopher.terry[at]biology.ox.ac.uk) if anything is missing, breaks or is unclear. 

You are free to make use of the code in any way you wish, although it would be good to cite the paper. Most of the underlying environmental data is fully free (Open Government License), but it is probably better to start from scratch for most uses anyway. See the original sources for the butterfly and bird data citations - it is open data but needs proper citation.  

For the simulation study, all code, model inputs and summary results as required to fully replicate the study are included, but the actual JSDM model fits are not included because they are very large (multiple Gb). They should be totally reconstituitable from the code however.

For the two test cases (UK birds and butterflies), complete original raw data are not included as the files are very large and openly available elsewhere as described in the paper. Processed species x site data and environmental variables are included to allow full replication of the model fitting part of the paper.

All results were obtained using R Version 4.1.1, mostly using CRAN packages. The main exception was `HMSC`, which is a github package (https://github.com/guiblanchet/HMSC). Not that this is not the same as the `Hmsc` package on CRAN!

The simulations made heavy use of high-performance computing and are not suitable for rerunning on a normal computer. 

Core scripts are presented as `.rmd` files, and code is organised into simulation, bird and butterfly folders, with hopefully informative names. Most probably won't just knit straight out of ht box because I have not included the certain very large files. (They could be regenerated with the included code if needed).

The output figures are collated into a single SI document with `SI_Figures_and_Captions.rmd`. Helper functions to be sourced in are stored in the `Scripts` folder. The `BashScripts` folder includes both the bash scripts used to set up the HPC and the R scripts that are called to run the simulations and or fitting on HPC. 


## Description of the data and file structure

### Empirical Data

#### Butterflies

Processed data for the butterfly study is in the `Data/` folder. `ButterflyMarkdowns/ButterflyDataPrep2.rmd` describes the original source of the data and the transformation and filterring steps carried out. Transects are identified by their Ordnance Survey 1km grid number (OSGR). 

Contents:

- `ButterflyNames.csv` - Reference table of scientific and common names
- `MEM_justfocalsites_first10.csv` - table of spatial predictors for the focal sites
- `Site_selection2_CHALK.csv` - table of environmental predictors for the focal sites
- `SiteDetails.csv` table of location and land cover for all transect squares
- `YearSlices2000x2020_v2_2perc_nomoths` - R list object. each list entry is a species x site matrix passed to the the JSDM scripts
- `UKBMS_selected_sites_FilledGaps_v2_2perc_nomoths.csv`   - as above, but in a long form .csv format 

#### Birds

Processed data for the bird study is in `BirdData_clean/` folder (Apologies for the naming clash). 

Contents:

- `Bird_X_formatted.csv` Environmental and spatial predictor variables
- `BOU_British_List.csv`  List used to determine native status of species
- `focalBirds_AllC_1970.csv` species by site presence absence matrices 1970 all data quality
- `focalBirds_AllC_2010.csv` species by site presence absence matrices 2010 all data quality
- `focalBirds_High_1970.csv` species by site presence absence matrices 1970 excluding uncertain observations
- `focalBirds_High_2010.csv` species by site presence absence matrices 2010 excluding uncertain observations
- `finalSP_HQ` R object with list of species to include in the final analysis
- `finalSP_all` R object with list of species to include in the final analysis
- `species_lookup.csv` Species names table (from BTO)

### Model fits and outputs

#### Simulations


`SavedObjects/` folder contains .csv files including the R2 of each of the sub models in the main trial (`LatinSimulation,csv`) and the false- absence trial (`TrialObsQual.csv`). These are analysed with the .rmd files in `SimulationMarkdowns`.

The TrackThroughTime folders contain results from the trial that was used to make the example figure in the main text. See `SimulationExample.rmd`

#### Butterflies

The summary variance partitioning results for each year are saved as .csvs in `SavedObjects/HMSC_Bf2_perc2_log_env`.
The full model fits are too large. 

#### Birds

The summary variance partitioning results for each year are saved as R objects in `SavedObjects/BirdVPs`.
The full model fits are too large. Code to make use of the R objects is in `Birds Analysis.Rmd`

## Sharing/Access information


You are free to make use of the code in any way you wish, although it would be good to cite the paper. Most of the underlying environmental data is fully free (Open Government License), but it is probably better to start from scratch for most uses anyway. See the original sources for the butterfly and bird data citations - it is open data but needs proper citation.  


