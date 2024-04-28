#!/bin/bash

#SBATCH -A naiss2023-22-1378
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J makeADiHS
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

cd /crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/
module load  R/4.2.1
module load  R_packages/4.2.1
Rscript --vanilla /crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/7-makefileAD_v3.R
