#!/bin/bash -l
 
#SBATCH -A ...
#SBATCH ...

# how to run :
#for trait in Ovule_Number Long_Stamens Short_Stamens Petal_Area Petal_Length Petal_Width Sepal_Area Sepal_Length Sepal_Width Leaf_Area Leaf_Length Leaf_Width flowering_time; do sbatch 1-GWAs.sh $trait; done

module load bioinfo-tools
module load GEMMA/0.98.1

trait=$1

# cd refers to the file organization in the cluster
cd A_thaliana/GWAs/Ushape/

#mkdir $trait
cd $trait

# 1 - GWAs

#cp A_thaliana/bed/SNP_1001g_filtered.bed SNP_1001g_filtered.bed
#cp A_thaliana/bed/SNP_1001g_filtered.bim SNP_1001g_filtered.bim
#cp A_thaliana/phenotypes/SNP_1001g_filtered_$trait.fam SNP_1001g_filtered.fam

# relatedness matrix
#gemma -bfile SNP_1001g_filtered -gk 2 -o SNP_1001g_filtered_$trait

# gwas
gemma -bfile SNP_1001g_filtered -k output/SNP_1001g_filtered_$trait.sXX.txt -lmm 4 -o SNP_1001g_filtered_$trait
