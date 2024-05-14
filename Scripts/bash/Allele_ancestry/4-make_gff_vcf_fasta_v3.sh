#!/bin/bash

#SBATCH -A XXXXX
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J GFFVCFFA
#SBATCH --mail-user xxxxx
#SBATCH --mail-type=ALL

genelist=$1

module load bioinfo-tools
module load BEDTools/2.29.2
module load SeqKit/0.15.0
module load vcftools/0.1.16

cd $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447

source progress_bar.sh

#mkdir $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/$genelist.files/
#mv $genelist $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/$genelist.files/
cd $PATH-TO-SERVER/A_thaliana/temp/temp/ancestry447/$genelist.files/

while read HOG OG N ahal alyr atha
do

echo $HOG

# atha
echo "atha"
touch $atha.fasta
grep $atha $PATH-TO-SERVER/A_thaliana/REF/447/Athaliana_447_Araport11.gene.gff3 | grep -w "gene" > $atha.t.gff
if [ ! -s $atha.t.gff ]; then
  grep $atha $PATH-TO-SERVER/A_lyrata/REF/Alyrata_107_v1.0.gene.gff3 | grep -w mRNA > $atha.t.gff
fi
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $4-5000, $5+5000}' $atha.t.gff > $atha.bed
awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4-5000, $5+5000, $6, $7, $8, $9}' $atha.t.gff > $atha.gff
rm $atha.t.gff
scaf=$(grep $atha $PATH-TO-SERVER/A_thaliana/REF/447/Athaliana_447_Araport11.gene.gff3 | head -n 1 | awk '{print $1}')
max=$(seqkit grep -p $scaf $PATH-TO-SERVER/A_thaliana/REF/447/Athaliana_447_TAIR10.fa | sed '1d' | tr -d '\n' | wc -c)
echo -e "$scaf\t0\t$max" >> $atha.bed
start=$(sort -nrk2 $atha.bed | sed -n 1p | awk '{print $2}')
stop=$(sort -nk3 $atha.bed | sed -n 1p | awk '{print $3}')
echo -e "$scaf\t$start\t$stop" > $atha.bed
bedtools getfasta -fi $PATH-TO-SERVER/A_thaliana/REF/447/Athaliana_447_TAIR10.fa -bed $atha.bed > s_$atha.fasta
cat s_$atha.fasta >> $atha.fasta

if [ -f $atha.recode.vcf ];
 then
 vcflines=$(grep -v "#" $atha.recode.vcf | wc -l)
  if [ $vcflines -lt 1 ];
   then
   rm $atha.recode.vcf
  fi
fi

if [ ! -f $atha.recode.vcf ];
then
vcftools --vcf $PATH-TO-SERVER/A_thaliana/VCF/1001genomes_snp_2all_mac1_rel.95_Chr.vcf --chr $scaf --from-bp $start --to-bp $stop --recode --out $atha
else
echo "vcf exists"
fi

# alyr
echo "alyr"
if [ $alyr == "NA" ]
  then
echo -e ">$alyr\nATG" > $alyr.fasta
cat $alyr.fasta >> $atha.fasta
  else
    alyrid=$(awk -F"\t" -v "var=$alyr" '$4 == var' $PATH-TO-SERVER/A_lyrata/REF/Alyrata_107_v1.0.annotation_info.txt | awk '{print $3}')
    grep ID=$alyrid. $PATH-TO-SERVER/A_lyrata/REF/Alyrata_107_v1.0.gene.gff3 | grep gene > $alyr.t.gff
    if [ ! -s $alyr.t.gff ]; then
    grep $alyrid $PATH-TO-SERVER/A_lyrata/REF/Alyrata_107_v1.0.gene.gff3 | grep mRNA > $alyr.t.gff
    fi
    awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $4-5000, $5+5000}' $alyr.t.gff > $alyr.bed
    awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4-5000, $5+5000, $6, $7, $8, $9}' $alyr.t.gff > $alyr.gff
    rm $alyr.t.gff
    scaf=$(grep ID=$alyrid. $PATH-TO-SERVER/A_lyrata/REF/Alyrata_107_v1.0.gene.gff3 | grep gene | awk '{print $1}')
    max=$(seqkit grep -p $scaf $PATH-TO-SERVER/A_lyrata/REF/Alyrata_384_v1.fa | sed '1d' | tr -d '\n' | wc -c)
    echo -e "$scaf\t0\t$max" >> $alyr.bed
    start=$(sort -nrk2 $alyr.bed | sed -n 1p | awk '{print $2}')
    stop=$(sort -nk3 $alyr.bed | sed -n 1p | awk '{print $3}')
    echo -e "$scaf\t$start\t$stop" > $alyr.bed
    bedtools getfasta -fi $PATH-TO-SERVER/A_lyrata/REF/Alyrata_384_v1.fa -bed $alyr.bed > s_$alyr.fasta
    cat s_$alyr.fasta >> $atha.fasta
fi 

# ahal
echo "ahal"
if [ $ahal == "NA" ]
  then
echo -e ">$ahal\nATG" > $ahal.fasta
cat $ahal.fasta >> $atha.fasta
  else
    ahalid=$(awk -F"\t" -v "var=$ahal" '$4 == var' $PATH-TO-SERVER/A_halleri/REF/Ahalleri_264_v1.1.annotation_info.txt | awk '{print $3}')
    grep $ahalid $PATH-TO-SERVER/A_halleri/REF/Ahalleri_264_v1.1.gene.gff3 | grep gene > $ahal.t.gff
    if [ ! -s $ahal.t.gff ]; then
    grep $ahalid $PATH-TO-SERVER/A_halleri/REF/Ahalleri_264_v1.1.gene.gff3 | grep mRNA > $ahal.t.gff
    fi
    awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $4-5000, $5+5000}' $ahal.t.gff > $ahal.bed
    awk -F"\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4-5000, $5+5000, $6, $7, $8, $9}' $ahal.t.gff > $ahal.gff
    rm $ahal.t.gff
    scaf=$(cat $ahal.gff | awk '{print $1}')
    max=$(seqkit grep -p $scaf $PATH-TO-SERVER/A_halleri/REF/Ahalleri_264_v1.fa | sed '1d' | tr -d '\n' | wc -c)
    echo -e "$scaf\t0\t$max" >> $ahal.bed
    start=$(sort -nrk2 $ahal.bed | sed -n 1p | awk '{print $2}')
    stop=$(sort -nk3 $ahal.bed | sed -n 1p | awk '{print $3}')
    echo -e "$scaf\t$start\t$stop" > $ahal.bed
    bedtools getfasta -fi $PATH-TO-SERVER/A_halleri/REF/Ahalleri_264_v1.fa -bed $ahal.bed > s_$ahal.fasta
    cat s_$ahal.fasta >> $atha.fasta
fi 


rm *$alyr* *$ahal* *.bed *.log s_* NA.fasta


current=$(grep -n $HOG $genelist | cut -f1 -d:)
total=$(wc -l $genelist)
show_progress $current $total

done < $genelist