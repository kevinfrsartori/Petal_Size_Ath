#!/bin/bash -l
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -J gwas
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

# how to run :
#for trait in Ovule_Number Long_Stamens Short_Stamens Petal_Area Petal_Length Petal_Width Sepal_Area Sepal_Length Sepal_Width Leaf_Area Leaf_Length Leaf_Width flowering_time; do sbatch 1-GWAs.sh $trait; done

ml bioinfo-tools GEMMA/0.98.1 plink/1.90b4.9 

trait=$1
echo $trait

cd $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/

#mkdir $trait
cd $trait

# 1 - GWAs

#cp $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered.bed SNP_1001g_filtered.bed
#cp $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered.bim SNP_1001g_filtered.bim
#cp $PATH-TO-SERVER/A_thaliana/phenotypes/SNP_1001g_filtered_$trait.fam SNP_1001g_filtered.fam

# relatedness matrix
#gemma -bfile SNP_1001g_filtered -gk 2 -o SNP_1001g_filtered_$trait

#eigen decomposition
#gemma -bfile SNP_1001g_filtered -k output/SNP_1001g_filtered_$trait.sXX.txt -eigen -o eigen.$trait

#LD pruning
for r in {1..10}
do
plink --bfile SNP_1001g_filtered --indep 20 5 $r
done

# gwas
#gemma -bfile SNP_1001g_filtered -k output/SNP_1001g_filtered_$trait.sXX.txt -lmm 4 -o SNP_1001g_filtered_$trait

# to run when gwas are done:
# ls $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/*/output/SNP_1001g_filtered_*.assoc.txt > assoc_dir.txt
# while read dir; do cp $dir $PATH-TO-SERVER/kevinfrs/private/Athaliana_flowersize/gwas/ ; done < assoc_dir.txt

# SNPs effect
#gemma -bfile SNP_1001g_filtered -k output/SNP_1001g_filtered_$trait.sXX.txt -bslmm 1 -o bslmm_SNP_1001g_filtered_$trait

# 2 - extract effects
# keep snp lists

#head -n 1 output/bslmm_SNP_1001g_filtered_$trait.param.txt > output/bslmm_flct_$trait.param.txt
#while read snp chr pos
#do
#grep $snp output/bslmm_SNP_1001g_filtered_$trait.param.txt >> output/bslmm_flct_$trait.param.txt
#done < $PATH-TO-SERVER/kevinfrs/private/Athaliana_flowersize/snp_lists/flct_Hits_$trait.txt

#head -n 1 output/bslmm_SNP_1001g_filtered_$trait.param.txt > output/bslmm_pvl4_$trait.param.txt
#while read snp chr pos
#do
#grep $snp output/bslmm_SNP_1001g_filtered_$trait.param.txt >> output/bslmm_pvl4_$trait.param.txt
#done < $PATH-TO-SERVER//kevinfrs/private/Athaliana_flowersize/snp_lists/pvl4_Hits_$trait.txt

# to run after:
# ls $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/*/output/bslmm_flct_*.param.txt > param_dir.txt
# ls $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/*/output/bslmm_pvl4_*.param.txt >> param_dir.txt
# while read dir; do cp $dir $PATH-TO-SERVER/kevinfrs/private/Athaliana_flowersize/gwas/ ; done < param_dir.txt

