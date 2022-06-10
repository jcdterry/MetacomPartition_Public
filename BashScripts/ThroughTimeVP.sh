#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N  ThroughTimeVP1 
#$ -t 1-75

module -s load R/4.1.1

# Replace the following line with a program or command

R CMD BATCH   GitHub/MC_coherence/BashScripts/ThroughTimeVP.r    GitHub/MC_coherence/ROuts/ThroughTimeVP.r-${SGE_TASK_ID}.Rout