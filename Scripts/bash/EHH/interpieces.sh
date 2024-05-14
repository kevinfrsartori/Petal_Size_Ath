#!/bin/bash
 
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:01:00
#SBATCH -J interpiece
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL


cd $PATH-TO-SERVER/private/Athaliana_flowersize/rehh/

touch interpiece_snps.txt  
for i in 1 2 3 4 5
do
while read piece
do
if [ $piece != "aa" ]
then
awk 'NR==1{print $1, $2, $3, $4, $5, $3-100000, $3+100000}' $PATH-TO-SERVER/A_thaliana/fastphase/SNP_1001g_filtered.$i.$piece.bimlike >> interpiece_snps.txt
fi
done < $PATH-TO-SERVER/A_thaliana/fastphase/letters.$i.txt
done

while read rs chrm pos ref alt swin ewin
do
sbatch bed2fastphase2ehhplot.sh $pos $chrm $swin $ewin
done < interpiece_snps.txt

