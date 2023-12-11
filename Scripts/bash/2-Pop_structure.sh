#!/bin/bash
 
#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J Cr_adm
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL


# Script started 2023-05-08, Kevin Sartori
# Analysis of A th sequences
# individuals, variants

# 1 - Prepare files and run admixture for K in {1..15}
#-----------------------------------------------------

cd /crex/proj/snic2020-16-182/

mkdir A_thaliana/admixture/Petal_Size_Ath

cp /domus/h1/kevinfrs/private/Athaliana_flowersize/studied_acc.txt A_thaliana/admixture/Petal_Size_Ath/studied_acc.txt

ml bioinfo-tools plink/1.90b4.9 ADMIXTURE/1.3.0

# - Prunne dataset
plink --bfile A_thaliana/bed/SNP_1001g_filtered --allow-extra-chr --snps-only --maf 0.05 --geno 0.05 --keep A_thaliana/admixture/Petal_Size_Ath/studied_acc.txt --make-bed --out A_thaliana/admixture/Petal_Size_Ath/SNP_1001g_filtered_Petal_Size_Ath

# - ADMIXTURE
for K in {1..15};
do 
sbatch admixture.sh $K /crex/proj/snic2020-16-182/A_thaliana/admixture/Petal_Size_Ath SNP_1001g_filtered_Petal_Size_Ath.bed
done

# 2 - When 1 is done, combine results and choose best K
#------------------------------------------------------

grep -h CV A_thaliana/admixture/Petal_Size_Ath/log*.out
echo "admixture analysis done."

# Plot result in command line :
grep -h CV A_thaliana/admixture/Petal_Size_Ath/log*.out | awk -F'[ ]' '{print $4}' > CV.txt
grep -h CV A_thaliana/admixture/Petal_Size_Ath/log*.out | awk -F'[=]' '{print $2}'| awk -F'[)]' '{print $1}' > K.txt
paste K.txt CV.txt | sort -n -k 1 | head -n 1 > adm.txt
paste K.txt CV.txt | sort -n -k 1 > adm2.txt
cat adm2.txt >> adm.txt 
awk '{print $2}' adm.txt | gnuplot -e "set terminal dumb; set xrange [0:10]; set title 'JAG ADMIXTURE'; plot '-' with impulses ls -1"



# - Best solution is 6 groups, 6 haplotypes
mv SNP_1001g_filtered_Petal_Size_Ath.6.Q A_thaliana/admixture/Petal_Size_Ath/SNP_1001g_filtered_Petal_Size_Ath.6.Q
mv SNP_1001g_filtered_Petal_Size_Ath.6.P A_thaliana/admixture/Petal_Size_Ath/SNP_1001g_filtered_Petal_Size_Ath.6.P
rm SNP_1001g_filtered_Petal_Size_Ath.*.Q
rm SNP_1001g_filtered_Petal_Size_Ath.*.P

# - Plot
#module load  R/4.2.1
#module load RStudio/2022.07.1-554
#rstudio
# Open plot_admixture.R 

