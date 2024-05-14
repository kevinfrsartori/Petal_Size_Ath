#!/bin/bash
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 4-00:00:00
#SBATCH -J EHH
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

# How to run:
# sbatch fastphase.sh [FILENAME.phase.inp] [OUTPUTNAME]

file=$1
out=$2

module load bioinfo-tools
module load fastPHASE/1.4.8

fastPHASE -T10 -o$out $file



