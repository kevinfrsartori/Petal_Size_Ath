#!/bin/bash

#SBATCH -M XXXXX
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J prepare
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

HS=$1
cd $PATH-TO-SERVER

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
