#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -N  BF_v2_1
#$ -t 1-21

module -s load R/4.1.1
module -s load gcc/10.2.0  

# Replace the following line with a program or command

R CMD BATCH   GitHub/MC_coherence/BashScripts/BF_Version2_perc2_log_env.R    GitHub/MC_coherence/ROuts/BF_Version2_perc2_log_env-${SGE_TASK_ID}.Rout