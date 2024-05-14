#!/bin/bash -l
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J mvgwas
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

#how to run:
# while read trait1 trait2; do sbatch 1-2-mvGWAs.sh $trait1 $trait2; done < $PATH-TO-SERVER/private/Athaliana_flowersize/gwas/multi/Pairwise_traits_comb.txt

trait1=$1
trait2=$2

module load bioinfo-tools
module load GEMMA/0.98.1

#mkdir $PATH-TO-SERVER/kevinfrs/private/Athaliana_flowersize/gwas/multi/
cd $PATH-TO-SERVER/kevinfrs/private/Athaliana_flowersize/gwas/multi/

# 1 - GWAs

#cp $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered.bed SNP_1001g_filtered.bed
#cp $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered.bim SNP_1001g_filtered.bim
#cp $PATH-TO-SERVER/A_thaliana/phenotypes/SNP_1001g_filtered_multi.fam SNP_1001g_filtered.fam

# relatedness matrix
#recovered from previous computation

# mvgwas
gemma -bfile SNP_1001g_filtered -k $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/flowering_time/output/SNP_1001g_filtered_flowering_time.sXX.txt -lmm 4 -n $trait1 $trait2 -o SNP_1001g_filtered_multi_$trait1.$trait2


