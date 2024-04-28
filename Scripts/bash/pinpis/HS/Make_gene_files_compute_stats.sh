#!/bin/bash

#SBATCH -M snowy
#SBATCH -A naiss2023-5-501
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J Makefiles_ComputePiNPiS
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=FAIL

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0
module load SeqKit/0.15.0
module load R/4.3.1
module load R_packages/4.3.1


piece=$1
HS=$2

# make fasta and vcf per gene, compute stats with R
echo "make fasta and vcf for:" 
while read chrm gene start end strand
do
echo $gene
echo $chrm

seqkit grep -p $gene /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/A_thaliana/Athaliana_167_TAIR10.cds.fa > fasta.fasta

vcftools --gzvcf /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/$HS/1001genomes_snp_2all_mac1_rel.95_Chr_$chrm.recode.vcf --chr Chr$chrm --from-bp $start --to-bp $end --recode --out vcf

# Extract gff 
grep $gene /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/A_thaliana/Athaliana_167_TAIR10.gene.gff > $gene.tmp.gff
grep 'CDS' $gene.tmp.gff > gff.gff
rm $gene.tmp.gff


Rscript --vanilla ../StdFiles_to_PiNPiS.R $piece
rm fasta.fasta vcf.recode.vcf gff.gff full_fasta.fasta vcf.log

done < ../$piece