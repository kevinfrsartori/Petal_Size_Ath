#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J makeADiHS
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

cd /$PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/
module load  R/4.2.1
module load  R_packages/4.2.1
Rscript --vanilla $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/makefileAD_v3.R
