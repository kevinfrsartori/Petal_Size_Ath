#!/bin/bash
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 5-00:00:00
#SBATCH -J EHH
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

# How to run:
# sbatch scan_EHH.sh [chromosome number] ["piece" ID -or- SNP pos]
# How to parallelize...
# ... for chrm pieces:
#for i in 1 2 3 4 5; do while read name piece; do sbatch scan_EHH.sh $i $piece; done < $PATH-TO-SERVER/A_thaliana/fastphase/splitlistletters.$i.txt; done
# ... for interpieces:
# while read name chrm snp ref alt swin send; do sbatch scan_EHH.sh $chrm $snp; done < interpiece_snps.txt

chrm=$1
piece=$2

module load R/4.2.1
module load R_packages/4.2.1
module load RStudio/2022.07.1-554

Rscript --vanilla scan_rehh.R $chrm $piece



