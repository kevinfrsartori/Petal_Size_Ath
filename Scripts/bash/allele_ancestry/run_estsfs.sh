#!/bin/bash

#SBATCH -A naiss2023-22-1378
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J est_sfs
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

genelist=$1


cd /crex/proj/snic2020-16-182/software/est_sfs/

#./est-sfs config-kimura.txt data-file.$SOC_x.txt seed-file.txt output-file-sfs.$SOC_x.txt output-file-pvalues.$SOC_x.txt
#./est-sfs config-JC.txt data-file.$SOC_x.txt seed-file.txt output-file-sfs.$SOC_x.txt output-file-pvalues.$SOC_x.txt
./est-sfs config-rate6.txt data-file.$genelist.txt seed-file.txt output-file-sfs.$genelist.txt output-file-pvalues.$genelist.txt
