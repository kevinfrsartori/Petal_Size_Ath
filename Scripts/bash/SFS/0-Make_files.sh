#!/bin/bash

#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J freq_VCF_HS
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL


ml bioinfo-tools plink/1.90b4.9 vcftools/0.1.16

cd /crex/proj/snic2020-16-182/A_thaliana/VCF/
#vcftools --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf --remove-indels --max-alleles 2 --mac 1 --recode --out 1001genomes_snp_2all_mac1
#After filtering, kept 1135 out of 1135 Individuals
#After filtering, kept 10700811 out of a possible 12883854 Sites

#plink --vcf 1001genomes_snp_2all_mac1.recode.vcf --make-bed --out ../bed/1001genomes_snp_2all_mac1
## Relationship-based prunning 
#for i in 0.1 0.25 0.50 0.75 0.9 0.95 0.99 0.999 0.9999
#do
#plink --bfile ../bed/1001genomes_snp_2all_mac1 --rel-cutoff $i
#echo $i >> cutoff.txt
#grep 'excluded' plink.log >> excluded.txt
#rm plink.log
#done

#paste cutoff.txt excluded.txt > result.txt
#rm cutoff.txt excluded.txt
#rm plink.nosex

# Plot result in command line :
#awk '{print $2}' result.txt | gnuplot -e "set terminal dumb 50 50; set xrange [-1:10]; set title 'Nb excluded accessions'; plot '-' with impulses ls -1"
# plateau reached for --rel-cutoff 0.95

# Make file
#plink --bfile ../bed/1001genomes_snp_2all_mac1 --rel-cutoff 0.95 --make-bed --out ../bed/1001genomes_snp_2all_mac1_rel.95
#plink --bfile ../bed/1001genomes_snp_2all_mac1 --rel-cutoff 0.95 --recode vcf --out 1001genomes_snp_2all_mac1_rel.95


# VCF modifs
grep '#' /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95.vcf > /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.vcf
sed -i 's/contig=<ID=/contig=<ID=Chr/g' /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.vcf
grep -v '#' /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95.vcf > /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_table.vcf
cat /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_table.vcf | awk '$1="Chr"$1' OFS='\t' >> /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.vcf
ml bioinfo-tools GATK/4.3.0.0
gatk IndexFeatureFile -I /crex/proj/snic2020-16-182/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.vcf


# check correspondance between fasta and 1001g vcf
# fasta is organized in lines of 79 nucleotides
# check begining of chr1
grep -v "##" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | grep "Chr1" | head -n 10 | awk '{print $1,$2,$3,$4,$5,$6}' OFS="\t" > vcftop10.txt
awk '{printf"%.10f\n",$2/79;}' vcftop10.txt | awk -F"." '{print $1}' > time79top.txt
paste vcftop10.txt time79top.txt > vcftop10X.txt
awk '$8=79*$7' vcftop10X.txt | awk '$9=$7+2' | awk '$10=$2-$8' | awk '{print $4, $9,$10}' OFS="\t" > vcftop10XY.txt
while read vcf X Y
do
echo "vcf has $vcf"
echo "fasta has "
head /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa -n $X | tail -n 1 | cut -c $Y
echo
done < vcftop10XY.txt


#check end
grep -v "##" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | grep "Chr1" | tail -n 20 | awk '{print $1,$2,$3,$4,$5,$6}' OFS="\t" > vcftail10.txt
awk '{printf"%.10f\n",$2/79;}' vcftail10.txt | awk -F"." '{print $1}' > time79tail.txt
paste vcftail10.txt time79tail.txt > vcftail10X.txt
awk '$8=79*$7' vcftail10X.txt | awk '$9=$7+2' | awk '$10=$2-$8' | awk '{print $4,$5, $9,$10}' OFS="\t" > vcftail10XY.txt
while read REF ALT X Y
do
echo "vcf has $REF $ALT"
echo "fasta has "
head /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa -n $X | tail -n 1 | cut -c $Y
echo
done < vcftail10XY.txt

while read vcf X Y
do
echo "vcf has $vcf"
echo "fasta has "
head /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa -n $X | tail -n 1 | cut -c $Y
echo
done < vcftail10XY.txt


#check all
grep -v "##" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | grep "Chr1" | awk '{print $1,$2,$3,$4,$5,$6}' OFS="\t" > vcf.txt
awk '{printf"%.10f\n",$2/79;}' vcf.txt | awk -F"." '{print $1}' > time79.txt
paste vcf.txt time79.txt > vcfX.txt
awk '$8=79*$7' vcfX.txt | awk '$9=$7+2' | awk '$10=$2-$8' | awk '{print $4, $9,$10}' OFS="\t" > vcfXY.txt
touch corresp.txt
while read vcf X Y
do
fasta=$(head /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa -n $X | tail -n 1 | cut -c $Y)
if [ $vcf != $fasta ]
then
echo "$X $Y" >> corresp.txt
fi
done < vcfXY.txt


# test other ref from ensembl
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/
# check correspondance between fasta and 1001g vcf
# fasta is organized in lines of 60 nucleotides
# check begining of chr1
grep -v "##" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | grep "Chr1" | head -n 100 | awk '{print $1,$2,$3,$4,$5,$6}' OFS="\t" > vcftop10.txt
awk '{printf"%.10f\n",$2/60;}' vcftop10.txt | awk -F"." '{print $1}' > time60top.txt
paste vcftop10.txt time60top.txt > vcftop10X.txt
awk '$8=60*$7' vcftop10X.txt | awk '$9=$7+2' | awk '$10=$2-$8' | awk '{print $4,$5, $9,$10}' OFS="\t" > vcftop10XY.txt
while read REF ALT X Y
do
echo "vcf has $REF $ALT"
echo "fasta has "
head /crex/proj/snic2020-16-182/A_thaliana/REF/ensembl/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa -n $X | tail -n 1 | cut -c $Y
echo
done < vcftop10XY.txt

#check end
grep -v "##" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | grep "Chr1" | tail -n 100 | awk '{print $1,$2,$3,$4,$5,$6}' OFS="\t" > vcftail10.txt
awk '{printf"%.10f\n",$2/60;}' vcftail10.txt | awk -F"." '{print $1}' > time60tail.txt
paste vcftail10.txt time60tail.txt > vcftail10X.txt
awk '$8=60*$7' vcftail10X.txt | awk '$9=$7+2' | awk '$10=$2-$8' | awk '{print $4,$5, $9,$10}' OFS="\t" > vcftail10XY.txt

while read REF ALT X Y
do
echo "vcf has $REF $ALT"
echo "fasta has "
head /crex/proj/snic2020-16-182/A_thaliana/REF/ensembl/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa -n $X | tail -n 1 | cut -c $Y
echo
done < vcftail10XY.txt

#########################
# The REF allele in the 1001g VCF is not always the one in the TAIR10 reference fasta
# but, allele in TAIT10 ref fasta is always one of the two in the vcf, the REF or the ALT
# So it should be ok...
# Can we just shorten Chr1 of the fasta or pretend CHR1 is longer in the vcf ?
# what is at the end of the fasta?
# below is error from GATK
#Found contigs with the same name but different lengths:
#contig reference = Chr1 / 30427671 <- the fasta
#contig features = Chr1 / 30427614 <- the vcf from 1001g
# fasta is 57 nuc shorter
#grep ">" /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa
#grep "Chr2" -B 5 /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa
#grep "Chr1" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | tail -n 1
# vcf last SNP is 30427613, it explains why the vcf head states that Chr1 is 30427614 nuc long ?
# I will just modify the vcf head by hand...
#grep "Chr1" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf -n | head -n 1
#sed -i  's/##contig=<ID=Chr1,length=30427614>/##contig=<ID=Chr1,length=30427671>/g' /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf

# same shit for Chr2
#grep -v "##" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf | grep "Chr2" | head -n 10 | awk '{print $1,$2,$3,$4,$5,$6}' OFS="\t" > vcftop10.txt
#awk '{printf"%.10f\n",$2/79;}' vcftop10.txt | awk -F"." '{print $1}' > time79top.txt
#paste vcftop10.txt time79top.txt > vcftop10X.txt
#grep ">" -n /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa
# chr2 start line 385163
#awk '$8=79*$7' vcftop10X.txt | awk '$9=$7+385164' | awk '$10=$2-$8' | awk '{print $4,$5, $9,$10}' OFS="\t" > vcftop10XY.txt
#while read REF ALT X Y
#do
#echo "vcf has $REF $ALT"
#echo "fasta has "
#head /crex/proj/snic2020-16-182/A_thaliana/REF/Athaliana_167_TAIR10.fa -n $X | tail -n 1 | cut -c $Y
#echo
#done < vcftop10XY.txt
#grep "Chr2" /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf -n | head -n 1
#  contig reference = Chr2 / 19698289
#  contig features = Chr2 / 19697998.
#sed -i  's/##contig=<ID=Chr2,length=19697998>/##contig=<ID=Chr2,length=19698289>/g' /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf

#Chr3
#sed -i  's/##contig=<ID=Chr3,length=23459135>/##contig=<ID=Chr3,length=23459830>/g' /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf

#Chr4
#sed -i  's/##contig=<ID=Chr4,length=18585005>/##contig=<ID=Chr4,length=18585056>/g' /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf

#Chr5
#sed -i  's/##contig=<ID=Chr5,length=26975451>/##contig=<ID=Chr5,length=26975502>/g' /crex/proj/snic2020-16-182/A_thaliana/VCF/SNP_1001g_filtered_Chr.vcf

