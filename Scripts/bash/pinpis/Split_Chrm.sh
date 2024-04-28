#!/bin/bash

#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J Make_filter_vcf
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0


chrm=$1
spp=$2

#make vcf per chromosome, filter out non variant and indels
vcftools --vcf /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_rel099.vcf --out /crex/proj/snic2020-16-182/$spp/VCF/chrm_$chrm --chr $chrm --recode --remove-indels --mac 1
#make index
gatk --java-options "-Xmx8G" IndexFeatureFile -I /crex/proj/snic2020-16-182/$spp/VCF/chrm_$chrm.recode.vcf

