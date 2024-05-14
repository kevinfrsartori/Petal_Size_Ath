#!/bin/bash

#SBATCH -M XXXXX
#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 03:00:00
#SBATCH -J Make_Compute
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=FAIL

HS=$1
cd $PATH-TO-SERVER

# 2- make gene list
###################
# important to optimize computation time
# one list per chromosome, chromosomes will be treated in parallel

# 2.1
#echo "make list of genes available in gff file"
#grep -E '>' $PATH-TO-SERVER/A_thaliana/REF/Athaliana_167_TAIR10.cds.fa | sed -r 's/^.{1}//' > gff.gene.list.txt
#echo "done"

# 2.2
#echo "make list of genes + infos"
#touch chrom.txt
#touch gene.txt
#touch start.txt
#touch end.txt
#touch strand.txt
#while read ATXG pacid poly loc id ann
#do
#grep $ATXG $PATH-TO-SERVER/A_thaliana/REF/Athaliana_167_TAIR10.gene.gff | cut -f1 | head -1 >> chrom.txt
#grep $ATXG $PATH-TO-SERVER/A_thaliana/REF/Athaliana_167_TAIR10.gene.gff | cut -f7 | head -1 >> strand.txt
#grep $ATXG $PATH-TO-SERVER/A_thaliana/REF/Athaliana_167_TAIR10.gene.gff > $ATXG.gff3
#cut -f4 -d$'\t' $ATXG.gff3 > $ATXG.bornes.txt
#cut -f5 -d$'\t' $ATXG.gff3 >> $ATXG.bornes.txt
#sort -n $ATXG.bornes.txt | head -1 >> start.txt
#sort -n $ATXG.bornes.txt | tail -1 >> end.txt
#echo $ATXG >> gene.txt
#rm $ATXG.gff3 $ATXG.bornes.txt
#done < gff.gene.list.txt
#paste chrom.txt gene.txt start.txt end.txt strand.txt | column -s $'\t' -t > gene.list.txt
#rm start* end* chrom* strand* gene.txt

#2.3
#echo "split list with 1000 genes per list"
#grep -v "ChrC" gene.list.txt | grep -v "ChrM" > 5Chrm.gene.list.txt
#cat 5Chrm.gene.list.txt | cut -c 4- > 5.gene.list.txt
#split -l 1000 5.gene.list.txt split.gene.list.
#ls split.gene.list.* > listlist.txt
#while read piece
#do
#mkdir temp_$piece
#touch result_$piece.txt
#echo "Gene_ID Nb_sites Nb_s_sites Nb_ns_sites Pis Pin PinPis" > result_$piece.txt
#echo $piece
#done < listlist.txt


# make fasta and vcf per gene, compute stats with R
while read piece
do
cd temp_$piece/
sbatch ../Make_gene_files_compute_stats.sh $piece $HS
cd ..
done < listlist.txt

# 3 combine pieces
#-----------------

ll $PATH-TO-SERVER/A_thaliana/PINPIS/wholegenome/pieces/*/result*
touch pinpis_snp_2all_mac1_rel.95.txt
echo "Gene_ID Nb_sites Nb_s_sites Nb_ns_sites Pis Pin PinPis" > pinpis_snp_2all_mac1_rel.95.txt
while read piece
do
sed '1d' $PATH-TO-SERVER/A_thaliana/PINPIS/wholegenome/pieces/$piece/result_$piece.txt >> pinpis_snp_2all_mac1_rel.95.txt
done < splitlist.txt

