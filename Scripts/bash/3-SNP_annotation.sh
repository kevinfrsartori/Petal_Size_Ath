#!/bin/bash
 
#SBATCH -A naiss2023-22-1378
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:20:00
#SBATCH -J prepare
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

ml bioinfo-tools plink/1.90b4.9

# Annotate SNPs with snpEff
#--------------------------
# how to run:
#ls snp_lists/ | sed -e 's/\.txt$//' > snpfiles.txt
# while read filename; do bash 3-SNP_annotation.sh $filename; done < snpfiles.txt

filename=$1
here=/crex/proj/snic2020-16-182/A_thaliana/temp/temp
there=/crex/proj/snic2020-16-182

# 1 make VCF from snp list
plink --bfile $there/A_thaliana/bed/SNP_1001g_filtered --extract $here/snp_lists/$filename.txt --make-bed --out $there/A_thaliana/GWAs/Petal_Size_Ath/annotations/$filename
plink --bfile $there/A_thaliana/GWAs/Petal_Size_Ath/annotations/$filename --recode vcf bgz --out $there/A_thaliana/GWAs/Petal_Size_Ath/annotations/$filename
# 2 run snpEff
cd $there/software/snpEff
java -Xmx8g -jar snpEff.jar Arabidopsis_thaliana $there/A_thaliana/GWAs/Petal_Size_Ath/annotations/$filename.vcf > $there/A_thaliana/GWAs/Petal_Size_Ath/annotations/$filename.ann.vcf

# How to make dataset for further analyses:
#ml R/4.2.1 R_packages/4.2.1 RStudio/2022.07.1-554
#dir=/crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/
#while read filename; do echo $filename; Rscript --vanilla 3-snp_annotation_to_genes_v4.R $filename $dir; done < snpfiles.txt
