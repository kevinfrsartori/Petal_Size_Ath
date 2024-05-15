#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J PiNPiS_Ath
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

# this script is made to keep track of the steps that allows the whole EHH analysis
# don't run it all, follow the steps

## Compute EHH scan genome wide

## 1 - convert bed files to fastphase

# we split the genome by the chromosomes and the chromosomes will be split in pieces for faster computation
for i in 1 2 3 4 5
do
sbatch rehh/bed2fastphase_chrm.sh $i
done
# This can take long time still

### IMPORTANT: the metrics computed at the begining and end of the chromosome pieces will be wrong because the haplotype length is amputated.
# need to run fastphase for the snps located at the end of each chromosome piece (but the last piece)
# for this, run interpieces.sh
bash rehh/interpieces.sh

## 2  - compute the metrics by using the r package rehh
# chrm pieces:
for i in 1 2 3 4 5; do while read name piece; do sbatch scan_EHH.sh $i $piece; done < $PATH-TO-SERVER/A_thaliana/fastphase/splitlistletters.$i.txt; done
# interpieces:
while read name chrm snp ref alt swin send; do sbatch scan_EHH.sh $chrm $snp; done < interpiece_snps.txt

# It takes a LOT of time for R to open the fastphase, this jobs might be killed before rehh managed till the end.
# Run the following if some pieces were not entirely read 
# Check if everything run
cd $PATH-TO-SERVER/A_thaliana/fastphase/
for i in 1 2 3 4 5
do
more splitlistletters.$i.txt
done
# See all the piece names and check in your data whats missing
cd $PATH-TO-SERVER/private/Athaliana_flowersize/rehh
# make file with missing pieces
echo -e "1 ab\n2 ae\n2 af\n3 aa\n4 aa\n4 ab\n4 af\n5 aa\n5 ab\n" > undo.txt
# Run again...
while read chrm piece
do
sbatch scan_EHH.sh $chrm $piece
done < undo.txt

# 3 - combine the chromosome pieces
# 3.1 split interpiece snp lists
for i in 1 2 3 4 5
do
awk -F' ' -v "var=$i" '$2==var' interpiece_snps.txt > interpiece_snps.$i.txt
done
# 3.2 run auto_combine
for i in 1 2 3 4 5
do
Rscript --vanilla auto_combine_pieces_scan_EHH.R $i
done

## 3 - make final table
for trait in Ovule_Number Long_Stamens Short_Stamens Petal_Area Petal_Length Petal_Width Sepal_Area Sepal_Length Sepal_Width Leaf_Area Leaf_Length Leaf_Width flowering_time
do
echo $trait
Rscript --vanilla rehh/final-table_v2.R $trait
done


# 4 - SNP EHH
# it is possible to run rehh and produce plots for specific SNPs by running the following
while read SNP_pos chrm swin ewin 
do
sbatch rehh/bed2fastphase2ehhplot.sh $SNP_pos $chrm $swin $ewin
done < PA_relevant_snps.txt

