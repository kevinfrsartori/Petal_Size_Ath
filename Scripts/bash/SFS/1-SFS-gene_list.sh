#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -J SFS_HS
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

#--------------------------------
# Compute SFS for a list of genes
#--------------------------------

#for trait in Ovule_Number Long_Stamens Short_Stamens Petal_Area Petal_Length Petal_Width Sepal_Area Sepal_Length Sepal_Width Leaf_Area Leaf_Length Leaf_Width flowering_time; do sbatch SFS-gene_list.sh $trait; done


trait=$1

echo
echo $trait
echo

# genes list
#-----------
ls $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/annotations/ | grep "annotated_simple_flct_Hits" | grep $trait | grep -v "AD" > $trait.ls.txt
while read list
do
sort -k 4 $PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/annotations/$list | awk '{print $4}' | uniq >> $trait.genelist.t.txt
done < $trait.ls.txt

sort -k 1 $trait.genelist.t.txt | uniq > $trait.genelist.txt
rm $trait.genelist.t.txt

grep "#" $PATH-TO-SERVER/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.Chr1.recode.vcf > $trait.genes.vcf

# combine vcfs of each gene
#--------------------------

#make gene gff list
while read gene
do
echo
echo $gene
echo 
grep $gene $PATH-TO-SERVER/A_thaliana/REF/447/Athaliana_447_Araport11.gene.gff3 | grep -w gene >> $trait.gene.list.gff
done < $trait.genelist.txt

ml bioinfo-tools vcftools/0.1.16
while read chrm source type start end rs strand ps info
do
echo
echo $chrm.$start
echo 
vcftools --vcf $PATH-TO-SERVER/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.$chrm.recode.vcf --chr $chrm --from-bp $start --to-bp $end --recode --out $trait.gene
grep -v "#" $trait.gene.recode.vcf >> $trait.genes.vcf
done < $trait.gene.list.gff


# Split the vcf per HS ranges
#----------------------------

# prepare HS lists

ls $PATH-TO-SERVER/private/Athaliana_flowersize/hs/ | grep txt > HSlist.txt
#while read HS
#do
#awk -F" " 'BEGIN {OFS=" "} {print $1"_"$1}' $PATH-TO-SERVER/private/Athaliana_flowersize/hs/$HS > updated.$HS
#done < HSlist.txt

# extract with vcftools

ml bioinfo-tools vcftools/0.1.16
while read HS
do
echo
echo writing $trait.$HS ...
echo
vcftools --vcf $trait.genes.vcf --keep updated.$HS --freq --out $trait.$HS
done < HSlist.txt



# Then use R script
#------------------

ml R/4.2.1 R_packages/4.2.1 
Rscript --vanilla handSFS_HS_genelist_v2.R $trait

