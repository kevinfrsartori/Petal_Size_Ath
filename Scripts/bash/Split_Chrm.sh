#!/bin/bash

#SBATCH -M snowy
#SBATCH -A naiss2023-5-501
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J chrm_vcf
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0


chrm=$1
HS=$2
cd /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/$HS/

#make vcf per chromosome, filter out non variant and indels
vcftools --vcf 1001genomes_snp_2all_mac1_rel.95_Chr_$HS.recode.vcf --out 1001genomes_snp_2all_mac1_rel.95_Chr_$HS_chrm_$chrm --chr Chr$chrm --recode
#make index
gatk --java-options "-Xmx8G" IndexFeatureFile -I 1001genomes_snp_2all_mac1_rel.95_Chr_$HS_chrm_$chrm.recode.vcf

