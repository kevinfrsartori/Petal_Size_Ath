#!/bin/bash

#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J freq_VCF_HS
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL


# Update on data used 11/14/2023
# - Full SNP dataset
# - Samples filtered by relatedness
# - Habitat suitability ranges redefined 
# see 0-make_files.sh


ml bioinfo-tools plink/1.90b4.9 vcftools/0.1.16


# Compute freq per HS classe
## How to launch this script ...
##cd private/Athaliana_flowersize/sfs/2-handDAF/
##ls /domus/h1/kevinfrs/private/Athaliana_flowersize/nichemodelling/accessions/ | grep txt > HSlist.txt
##while read HS ; do sbatch SFS_HS.sh $HS ; done < HSlist.txt
#...from here...
HS=$1
#mkdir /crex/proj/snic2020-16-182/A_thaliana/VCF/HS_vcf
cd /crex/proj/snic2020-16-182/A_thaliana/VCF
awk -F" " 'BEGIN {OFS=" "} {print $1"_"$1}' /domus/h1/kevinfrs/private/Athaliana_flowersize/nichemodelling/accessions/$HS > HS_vcf/$HS
vcftools --vcf 1001genomes_snp_2all_mac1_rel.95.vcf --keep HS_vcf/$HS --freq --out HS_vcf/$HS
# ... till here.
