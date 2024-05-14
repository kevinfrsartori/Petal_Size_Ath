#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J GFFVCFFA
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

genelist=$1
cd $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447
source progress_bar.sh
cd $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/$genelist.files/



# MODULES
#--------
echo "loading modules..."

module load ABINIT/8.10.3
module load Amber/20
module load Armadillo/9.700.2
module load CDO/1.9.5
module load GOTM/5.3-221-gac7ec88d
module load MUMPS/5.5.0
module load MUMPS/5.5.0-hybrid
module load NCL-graphics/6.6.2
module load NCO/4.8.1
module load NCO/4.9.2
module load NCO/4.9.3
module load OpenFOAM/6
module load OpenFOAM/7
module load OpenFOAM/v1912
module load PRISMS-PF/2.1.1
module load PROJ/8.1.0
module load QGIS/3.4.12
module load Rosetta/3.7
module load Siesta/4.1-MaX-1.0
module load Siesta/4.1-b4
module load Singular/4.1.2
module load WPS/4.1
module load WRF/4.1.3
module load WRF/4.1.3-dm+sm
module load WRF/4.1.3-dmpar
module load XCrySDen/1.5.60
module load XCrySDen/1.6.2
module load bioinfo-tools  augustus/3.3.3-CGP
module load deal.II/9.1.1-gcc
module load deal.II/9.1.1-intel
module load ncview/2.1.7
module load ncview/2.1.7-intel-2019b
module load wrf-python/1.3.1
module load MAFFT/7.429-GCC-8.2.0-2.31.1-with-extensions

echo "# modules loaded"
echo
echo "# Multiple sequence alignement"
echo
while read HOG OG N ahal alyr atha
do
mafft $atha.fasta > mafft_$atha
current=$(grep -n $HOG $genelist | cut -f1 -d:)
total=$(wc -l $genelist)
show_progress $current $total

done < $genelist


