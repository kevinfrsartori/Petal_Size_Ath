#!/bin/bash
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -J EHH
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL


# Extended Haplotype Homozygosity (EHH)

# How to run:
# sbatch bed2fastphase_chrm.sh [chromosome number]

chrm=$1

cd $PATH-TO-SERVER/A_thaliana/fastphase/

# 1 - from bed file to fastphase files
# (suppose that you made and filtered your .bed files before)
# 1.1 - make list of snp for rehh R package 
module load bioinfo-tools
module load plink/1.90b4.9
plink --bfile $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered --allow-extra-chr --chr $chrm --make-just-bim --out SNP_1001g_filtered.$chrm
#The maximum number of characters per line which can be read by fastPHASE is currently set to 500,000
# wc -l SNP_1001g_filtered.1.bim
# it should work but nope. Works with smaller pieces
awk -F"\t" '{ print $2 }' < SNP_1001g_filtered.$chrm.bim > fullsnplist.chrm-$chrm.txt
split -l 50000 fullsnplist.chrm-$chrm.txt snplist.chrm-$chrm.
ls snplist.chrm-$chrm.* > splitlist.$chrm.txt
awk -F"." '{print $3 }' < splitlist.$chrm.txt > letters.$chrm.txt
paste splitlist.$chrm.txt letters.$chrm.txt | column -s $'\t' -t > splitlistletters.$chrm.txt

while read dir let
do
plink --bfile $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered --chr $chrm --extract $dir --recode 12 fastphase --out $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$let
plink --bfile $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered --chr $chrm --extract $dir --make-just-bim --out $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$let
cat $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$let.bim | awk -F" " '{ print $2, $1, $4, $5, $6 }' > $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$let.bimlike
# 2 - run fastphase
sbatch $PATH-TO-SERVER/private/Athaliana_flowersize/rehh/fastphase.sh SNP_1001g_filtered.$chrm.$let.chr-$chrm.recode.phase.inp Chrm.$chrm.$let
done < splitlistletters.$chrm.txt

### IMPORTANT: the metrics computed at the begining and end of the chromosome pieces will be wrong because the haplotype length is amputated.
# need to run fastphase for the snps located at the end of each chromosome piece (but the last piece)
# bash interpieces.sh





