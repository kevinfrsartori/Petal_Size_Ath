#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J makeAD
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

# Add iHs and ADstate info to snp lists
#--------------------------

# how to run:
#ls snp_lists/ | grep "LD" | sed -e 's/\.txt$//' > snpfiles.txt
# 
#while read filename; do echo $filename; bash 4.1-snp_iHs_ADstate.sh $filename; done < snpfiles.txt

filename=$1

head -n 1 $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/Genome_wide_iHS.txt > snp_lists/$filename.iHS.AD.txt

while read rs chr ps
do
grep $rs $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/Genome_wide_iHS.txt >> snp_lists/$filename.iHS.AD.txt
done < snp_lists/$filename.txt

# How to make final dataset for further analyses:
#ml R/4.2.1 R_packages/4.2.1 RStudio/2022.07.1-554
#dir=$PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/annotations/
#while read filename; do echo $filename; Rscript --vanilla 4.2-snp_iHs_ADstate.R $filename $dir; done < snpfiles.txt

# How to add Gene Ontology annotation
#ml R/4.2.1 R_packages/4.2.1 RStudio/2022.07.1-554
#dir=$PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/annotations/
#while read filename; do echo $filename; Rscript --vanilla 4.3-func.annot.R $filename $dir; done < snpfiles.txt

