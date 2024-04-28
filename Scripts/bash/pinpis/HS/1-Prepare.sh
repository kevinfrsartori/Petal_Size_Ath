#!/bin/bash

#SBATCH -M snowy
#SBATCH -A naiss2023-5-501
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J prepare
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

HS=$1
cd /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/$HS/

# 1- Prepare dataset
#################
# One vcf per chromosome. VCF has to be made from the reference dataset used later (fasta, gff)
#ml bioinfo-tools vcftools/0.1.16
#sed -i 's/ /_/g' accessions_$HS.txt
#vcftools --vcf ../../A_thaliana/1001genomes_snp_2all_mac1_rel.95_Chr.vcf --keep accessions_$HS.txt --recode --out 1001genomes_snp_2all_mac1_rel.95_Chr_$HS

for i in 1 2 3 4 5
do
sbatch Split_Chrm.sh $i $HS
done
