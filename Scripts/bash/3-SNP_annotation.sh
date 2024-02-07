#!/bin/bash
 
#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J Cr_adm
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

ml bioinfo-tools plink/1.90b4.9


# Annotate SNPs with snpEff
#--------------------------

# rename variants if needed 
#head /crex/proj/snic2020-16-182/A_thaliana/bed/SNP_1001g_filtered.bim
#awk -F '\t' 'BEGIN{OFS="\t";} {$2="snp_"$1"_"$4}1' /crex/proj/snic2020-16-182/A_thaliana/bed/SNP_1001g_filtered.bim > /crex/proj/snic2020-16-182/A_thaliana/bed/SNP_1001g_filtered_.bim
#rm /crex/proj/snic2020-16-182/A_thaliana/bed/SNP_1001g_filtered.bim
#mv /crex/proj/snic2020-16-182/A_thaliana/bed/SNP_1001g_filtered_.bim /crex/proj/snic2020-16-182/A_thaliana/bed/SNP_1001g_filtered.bim

#cd /domus/h1/kevinfrs/private/Athaliana_flowersize
#ls /domus/h1/kevinfrs/private/Athaliana_flowersize/snp_lists | awk 'BEGIN { FS = "." } ; { print $1 }' > snps_lists.txt
#mkdir /crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/

while read name
do
cd /crex/proj/snic2020-16-182/
# 4.1 make VCF from snp list
plink --bfile A_thaliana/bed/SNP_1001g_filtered --extract /domus/h1/kevinfrs/private/Athaliana_flowersize/snp_lists/$name.txt --make-bed --out A_thaliana/GWAs/Petal_Size_Ath/annotations/$name
plink --bfile A_thaliana/GWAs/Petal_Size_Ath/annotations/$name --recode vcf bgz --out A_thaliana/GWAs/Petal_Size_Ath/annotations/$name
# 4.2 run snpEff
cd /home/kevinfrs/snpEff/
java -Xmx8g -jar snpEff.jar Arabidopsis_thaliana /crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/$name.vcf > /crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/$name.ann.vcf
done < snps_lists.txt

# make and filter datasets with R
#--------------------------------

cd /domus/h1/kevinfrs/private/Athaliana_flowersize
module load  R/4.3.1
module load   R_packages/4.3.1
#module load RStudio/2023.06.2-561
dir=/crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/

while read name
do
Rscript --vanilla 3-snp_annotation_to_genes_v3.R $name $dir
done < snps_lists.txt

# Extract gff around hits for figures
#------------------------------------
grep "Hit" snps_lists.txt > hits_lists.txt
while read name
do
chrm=Chr$(awk 'BEGIN { FS = "_" } ; { print $2 }' snp_lists/$name.txt)
min=$(awk 'BEGIN { FS = "_" } ; { print $3-20000 }' snp_lists/$name.txt)
max=$(awk 'BEGIN { FS = "_" } ; { print $3+20000 }' snp_lists/$name.txt)
grep $chrm /crex/proj/snic2020-16-182/A_thaliana/gff/TAIR10_GFF3_genes.gff | grep -w "gene" | awk -v var=$min '$5>var' | awk -v var=$max '$4<var' > snp_lists/$name.gff
done < hits_lists.txt


