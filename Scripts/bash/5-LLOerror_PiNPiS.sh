#!/bin/bash

#SBATCH -M snowy
#SBATCH -A naiss2023-5-501
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J PiNPiS_Ath
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL

# the lists of genotypes that belongs to the habitat suitability quantiles are there
#ll /crex/proj/snic2020-16-182/A_thaliana/PINPIS/PINPIS_HS/ 
cd /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/

module load bioinfo-tools
module load samtools/1.12
module load vcftools/0.1.16
module load GATK/4.2.0.0
module load SeqKit/0.15.0
module load R/4.3.1
module load R_packages/4.3.1
module load RStudio/2022.07.1-554

gene=$1
# make new folder for gene
mkdir PINPIS_HS/$gene
cd PINPIS_HS/$gene
cp /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/StdFiles_to_PiNPiS_leave20pout.R StdFiles_to_PiNPiS_leave20pout.R

grep $gene /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/HS02/5.gene.list.txt | head -n 1 > infogene.txt

while read chrm gene start end strand
do

# prepare result file
touch result_$gene.txt
echo "HS Replicate Pis Pin PinPis" > result_$gene.txt

for HS in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
# make files
echo $gene
echo $chrm
# Make fasta
seqkit grep -p $gene /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/A_thaliana/Athaliana_167_TAIR10.cds.fa > fasta.fasta
# Make vcf
vcftools --vcf /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/$HS/1001genomes_snp_2all_mac1_rel.95_Chr_$chrm.recode.vcf --chr Chr$chrm --from-bp $start --to-bp $end --recode --out vcf
# make gff 
grep $gene /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/A_thaliana/Athaliana_167_TAIR10.gene.gff > $gene.tmp.gff
grep 'CDS' $gene.tmp.gff > gff.gff
rm $gene.tmp.gff
# Run R
Rscript --vanilla StdFiles_to_PiNPiS_leave20pout.R $gene $HS
rm fasta.fasta vcf.recode.vcf gff.gff full_fasta.fasta
done

done < infogene.txt