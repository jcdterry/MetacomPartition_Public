#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N  prelimOQ_run 
#$ -t 1-1620

module -s load R/4.1.1
module -s load gcc/10.2.0  

# Replace the following line with a program or command

R CMD BATCH   GitHub/MC_coherence/BashScripts/TrialObsQual.R    GitHub/MC_coherence/ROuts/TrialObsQual.R-${SGE_TASK_ID}.Rout