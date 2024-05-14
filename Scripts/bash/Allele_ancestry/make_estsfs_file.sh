#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J GFFVCFFA
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

genelist=$1
cd $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/
source progress_bar.sh
cd $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/$genelist.files/

## - identify allele of the outgroups, make est-sfs file 

touch data-file.$genelist.txt
touch data-file.snps.$genelist.txt
echo -e "focalspecies\toutgroup1\toutgroup2" > data-file.$genelist.txt
echo -e "snpID\tgeneID" > data-file.snps.$genelist.txt

module load R/4.2.1
module load R_packages/4.2.1
module load RStudio/2022.07.1-554

while read HOG OG N ahal alyr atha
do

Rscript --vanilla $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/spot_alleles.R $atha $genelist

current=$(grep -n $HOG $genelist | cut -f1 -d:)
total=$(wc -l $genelist)
show_progress $current $total

done < $genelist


