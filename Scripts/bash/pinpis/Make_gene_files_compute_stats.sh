#!/bin/bash

#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J Makefiles_ComputePiNPiS
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0
module load SeqKit/0.15.0
module load R/4.2.1
module load RStudio/2022.07.1-554


piece=$1

# make fasta and vcf per gene, compute stats with R
echo "make fasta and vcf for:" 
while read chrm gene start end strand
do
echo $gene
echo $chrm

seqkit grep -p $gene /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.cds.fa > fasta.fasta

vcftools --gzvcf /crex/proj/snic2020-16-182/A_thaliana/VCF/chrm_$chrm.recode.vcf --chr $chrm --from-bp $start --to-bp $end --recode --out vcf

# Extract gff 
grep $gene /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.gene.gff > $gene.tmp.gff
grep 'CDS' $gene.tmp.gff > gff.gff
rm $gene.tmp.gff


Rscript --vanilla ../StdFiles_to_PiNPiS.R $piece
rm fasta.fasta vcf.recode.vcf gff.gff full_fasta.fasta vcf.log

done < ../$piece