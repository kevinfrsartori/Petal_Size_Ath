#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J est_sfs
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

genelist=$1


cd $PATH-TO-SERVER/software/est_sfs/

#./est-sfs config-kimura.txt data-file.$SOC_x.txt seed-file.txt output-file-sfs.$SOC_x.txt output-file-pvalues.$SOC_x.txt
#./est-sfs config-JC.txt data-file.$SOC_x.txt seed-file.txt output-file-sfs.$SOC_x.txt output-file-pvalues.$SOC_x.txt
./est-sfs config-rate6.txt data-file.$genelist.txt seed-file.txt output-file-sfs.$genelist.txt output-file-pvalues.$genelist.txt
