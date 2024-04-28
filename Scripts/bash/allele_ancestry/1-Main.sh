#!/bin/bash

#SBATCH -A snic2022-22-1195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -J WG_AD
#SBATCH --mail-user kevin.sartori@slu.se
#SBATCH --mail-type=ALL


## 1 - Run Orthofinder
ml bioinfo-tools BEDTools/2.29.2 SeqKit/0.15.0
#Alyr:
grep '>' /crex/proj/snic2020-16-182/proteomes/Alyr.protein.fa | head
cp /crex/proj/snic2020-16-182/proteomes/Alyr.protein.fa /crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/Alyr.protein.fa
#Ahal
grep '2.p' /crex/proj/snic2020-16-182/proteomes/Ahal.protein.fa | head
grep '1.p' /crex/proj/snic2020-16-182/proteomes/Ahal.protein.fa | awk '{print $1}' | awk '{gsub(/^.{1}/,"");}1' > Ahalgenelist.txt
seqkit grep -f /crex/proj/snic2020-16-182/proteomes/Ahalgenelist.txt /crex/proj/snic2020-16-182/proteomes/Ahal.protein.fa -o /crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/Ahal.protein.fa
#Atha 447:
grep '>' /crex/proj/snic2020-16-182/A_thaliana/REF/447/Athaliana_447_Araport11.protein.fa | awk '{print $1}' | awk '{gsub(/^.{1}/,"");}1' > /crex/proj/snic2020-16-182/A_thaliana/REF/447/Athagenelist.txt
# filtered for 1 row per gene with R	ml R/4.2.1 R_packages/4.2.1 RStudio/1.4.1106
seqkit grep -f /crex/proj/snic2020-16-182/A_thaliana/REF/447/Atha_gene_name_list.txt /crex/proj/snic2020-16-182/A_thaliana/REF/447/Athaliana_447_Araport11.protein.fa -o /crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/Atha.protein.fa

# update 2024-01: Re run orthofinder with option -y to split paralogs 
# sbatch 2-orthofinder.sh

## 2 - Make list of single copy ortholog with one gene per species (only and at least)  ml R/4.2.1 R_packages/4.2.1 RStudio/1.4.1106
grep AT /crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/OrthoFinder/Results_Jan30/Phylogenetic_Hierarchical_Orthogroups/N0.tsv > /crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/OrthoFinder/Results_Jan30/Phylogenetic_Hierarchical_Orthogroups/AT_N0.txt
#see 3-make_full_list.R

## 2.2 - Check that they all exist
bash gffexists.sh
# missings.txt stayed empty, so no gene missing in theory

## 3 - Make files needed fasta, gff and vcf per gene
split -l 100 one_ortho_per_species.txt gene_list.
ls gene_list.* > gene_lists.txt
while read list
do
sbatch 4-make_gff_vcf_fasta_v3.sh $list
done < gene_lists.txt

## 4 - run multialignement mafft
while read list
do
sbatch multi_seq_align.sh $list
done < gene_lists.txt

## 5 - extract Ath polymorphic sites for all species
while read list
do
sbatch make_estsfs_file.sh $list
done < gene_lists.txt


## 6 - run est-sfs
## 6.1 transfer files to est-sfs folder
while read list
do
cp /crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/$list.files/data-file.$list.txt /crex/proj/snic2020-16-182/software/est_sfs/data-file.$list.txt
sed -i '1d' /crex/proj/snic2020-16-182/software/est_sfs/data-file.$list.txt
done < gene_lists.txt
## 6.2 - run est-sfs
while read list
do
sbatch run_estsfs.sh $list
done < gene_lists.txt

## 7 - make file and merge with EHH
sbatch makefileAD.sh
