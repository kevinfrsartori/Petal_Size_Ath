#!/bin/bash
 
#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J Cr_adm
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL


K=$1
dir=$2
file=$3

ml bioinfo-tools ADMIXTURE/1.3.0

admixture --cv $dir/$file $K | tee $dir/log${K}.out
