#!/bin/bash

#########
#
# Running PINPIS-HS pipeline from scratch
#
#########

cd /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/

# the lists of genotypes that belongs to the habitat suitability quantiles are there
ll HS/

# 1 - Make new folder
mkdir PINPIS_HS/
for i in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
mkdir PINPIS_HS/$i
cp HS/accessions_$i.txt PINPIS_HS/$i/accessions_$i.txt 
done

# 2 - Prepare files, folders and genetic data
for i in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
cp 1-Prepare.sh PINPIS_HS/$i/
cp Split_Chrm.sh PINPIS_HS/$i/
done

for i in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
cd PINPIS_HS/$i/
sbatch 1-Prepare.sh $i
cd ../..
done

# 3 - Compute PINPIS
for i in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
cp 2-Make.sh PINPIS_HS/$i/
cp gene.list.txt PINPIS_HS/$i/
cp Make_gene_files_compute_stats.sh PINPIS_HS/$i/
cp StdFiles_to_PiNPiS.R PINPIS_HS/$i/ 
done

for i in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
echo $i
cd PINPIS_HS/$i/
bash 2-Make.sh $i
cd ../..
done

# 4 - ICI FAIRE UN SCRIPT POUR COMBINER LES PIECES
for i in HS01 HS02 HS03 HS04 HS05 HS06 HS07 HS08 HS09 HS10
do
echo $i
cd PINPIS_HS/$i/
ls /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/$i/result* > reslist.txt
echo "Gene_ID Nb_sites Nb_s_sites Nb_ns_sites Pis Pin PinPis" > wholegenome_pinpis_$i.txt
while read dir
do
cat $dir | sed '1d' >> wholegenome_pinpis_$i.txt
done < reslist.txt
cd ../..
done


# 5 - leave 20% error calculation for gene lists
while read gene name snp chrm rs signif allele effect type pos
do
echo $gene
sbatch 5-LLOerror_PiNPiS.sh $gene
done < annotated_simple_top100_Petal_Area.txt

# 6 - leave 20% error calculation for GWAs significant genes
for gene in AT1G77080 AT2G32370 AT3G03580 
do
echo $gene
sbatch 5-LLOerror_PiNPiS.sh $gene
done 

# extract final datasets
cp /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/*/wholegenome* outputs
cp /crex/proj/snic2020-6-185/nobackup/private/Violette/R/temp/PINPIS_HS/AT*/result_AT* outputs



