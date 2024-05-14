#!/bin/bash
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J EHH
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL


# Extended Haplotype Homozygosity (EHH)

# How to run:
# sbatch bed2fastphase2ehhplot.sh [snp position in bp] [chromosome number] [start window] [end window]

snp=$1
chrm=$2
frombp=$3
tobp=$4

# 1 - from bed file to fastphase files
# (suppose that you made and filtered your .bed files before)
# 1.1 - make list of snp for rehh R package 
module load bioinfo-tools
module load plink/1.90b4.9
plink --bfile $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered --allow-extra-chr --chr $chrm --from-bp $frombp --to-bp $tobp --make-just-bim --out $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$snp
# 1.2 - make fastphase format file
plink --bfile $PATH-TO-SERVER/A_thaliana/bed/SNP_1001g_filtered --allow-extra-chr --chr $chrm --from-bp $frombp --to-bp $tobp --recode 12 fastphase --out $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$snp
cat $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$snp.bim | awk -F" " '{ print $2, $1, $4, $5, $6 }' > $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$chrm.$snp.bimlike
# 2 - run fastphase
module load bioinfo-tools
module load fastPHASE/1.4.8
cd $PATH-TO-SERVER/A_thaliana/fastphase/
fastPHASE -T10 -oChrm.$chrm.$snp SNP_1001g_filtered.$chrm.$snp.chr-$chrm.recode.phase.inp

# 3 - Compute indices (ehh, ihs etc) and plot figure
module load R/4.2.1
module load R_packages/4.2.1
module load RStudio/2022.07.1-554
cd $PATH-TO-SERVER/private/Athaliana_flowersize/rehh/
Rscript --vanilla rehh.R $snp $chrm

