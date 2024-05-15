#!/usr/bin/env Rscript
argum=commandArgs(trailingOnly=T)
#argum=c("Petal_Area")
##########
#
# Extended haplotype homozygosity
# Make final dataset
# 2023-09-28
#
##########
library(vcfR)

# Starting the file from scratch
# this script might replace Athaliana_flowersize/7-Make_final_ADstate_table.R
# The script is not optimal because it uses the scanEHH files first to recover ref and alt alleles
# but then it uses WDiHS that is made from scanEHH but lost ref and alt in the script 7-makefileAD.R
# This has to be simplyfied.

# 1 - load EHH scan file (concatenate chromosomes)
for (i in 1:5) {
 if (i==1) {
   scanehh<-read.table(paste0("$PATH-TO-SERVER/A_thaliana/fastphase/scanEHH.chrm-",i,".txt"),h=T,sep=" ",dec=".")
 }else{
   scanehh<-rbind(scanehh,read.table(paste0("$PATH-TO-SERVER/A_thaliana/fastphase/scanEHH.chrm-",i,".txt"),h=T,sep=" ",dec="."))
 }
}
scanehh$snpID<-paste0("snp_",scanehh$CHR,"_",scanehh$POSITION)

# Open ancestry state file of the studied SNPs
ADstate<-read.table(paste0("$PATH-TO-SERVER/private/Athaliana_flowersize/ADstate_",argum[1],".csv"),sep=";",h=T)[,c(1:8)] 
# merge iHH from scan
ADstate<-merge(x = ADstate,y = scanehh[,c(3,4,5,7)],by = "snpID",all.x = T,sort = F)
# Create new column that we fill later on
ADstate$FREQ_REF<-NA
ADstate$FREQ_ALT<-NA
ADstate$a_is<-NA # the "a" in "FREQ_a" refer to REF or ALt allele ? we answer this based on reported frequencies
ADstate$iHH_REF<-NA
ADstate$iHH_ALT<-NA

# Report the REF and ALT frequencies from vcf files
# make loop
for (i in 1:dim(ADstate)[1]) {
vcf<-read.vcfR(paste0("$PATH-TO-SERVER/A_thaliana/Ushape_genomics/VCF/",ADstate$snpID[i],".vcf.gz"))
# all alleles
GT<-matrix(data = vcf@gt[,-1],nrow = 1, ncol = 1135,byrow = F  )
alleles<-c(substr(x = GT,start = 1,stop = 1),substr(x = GT,start = 3,stop = 3))
# allele freqs
freqs<-as.vector(table(alleles))
ADstate$FREQ_REF[i]<-freqs[1]/2270 
ADstate$FREQ_ALT[i]<-freqs[2]/2270
# for some unknown reasons, the allele frequencies are not exact compare with the frequency reported by est_sfs (removing NA doesnt make it better)
# so we use the closest frequency
ADstate$a_is[i]<-c("REF","ALT")[which.min(c(abs(ADstate$FREQ_a[i]-ADstate$FREQ_REF[i] ),abs(ADstate$FREQ_a[i]-ADstate$FREQ_ALT[i])))]
}

# if "a" is REF, then iHH_REF is IHH_a, and so on
if (any(ADstate$a_is =="ALT")){ ADstate$iHH_REF[which(ADstate$a_is =="ALT")]<-ADstate$IHHd[which(ADstate$a_is =="ALT")] } 
if (any(ADstate$a_is =="REF")){ ADstate$iHH_REF[which(ADstate$a_is =="REF")]<-ADstate$IHHa[which(ADstate$a_is =="REF")] } 
if (any(ADstate$a_is =="ALT")){ ADstate$iHH_ALT[which(ADstate$a_is =="ALT")]<-ADstate$IHHa[which(ADstate$a_is =="ALT")] } 
if (any(ADstate$a_is =="REF")){ ADstate$iHH_ALT[which(ADstate$a_is =="REF")]<-ADstate$IHHd[which(ADstate$a_is =="REF")] } 

# Load genome wide iHS (with standardized the statistic)
# iHS data
WDiHS<-read.table("$PATH-TO-SERVER/private/Athaliana_flowersize/rehh/Genome_wide_iHS.txt",h=T)
str(WDiHS)

# Merge
ADstate<-merge(ADstate,WDiHS[,c("snpID","FREQ_ancestral","FREQ_derived","IHH_ancestral","IHH_derived","iHS","Pvalue")] ,by="snpID",all.x = T,sort = F)

# for the missing ones
# Somesnps can be missing if there was no orthologs to compute the ancestry state
# In that case, we keep ref as the ancestral (often true), compute the iHS, and compute the Pvalue by reusing the scrit ihh2ihs

library(rehh)
# I made a custom version of rehh::ihh2ihs to compute the statistic and Pvalue for extra snps based on whole genome computation
source("$PATH-TO-SERVER/private/Athaliana_flowersize/rehh/whole_genome_ADstate/2-Blast_SingleCopyOrtho/ihh2ihs_ks.R")
# first make dataset with missing ancestry in the format required by rehh
res_ihh<-WDiHS[,c("chrm","pos","FREQ_ancestral","IHH_ancestral","IHH_derived","uniHS","iHS")]
mis_ihh<-ADstate[which(is.na(ADstate$Pvalue)),c("CHROM","POS","FREQ_REF","iHH_REF","iHH_ALT")]
# run modified function ihh2ihs_ks
mis_ihs<-ihh2ihs_ks(res_ihh,mis_ihh)
# then merge info
rownames(ADstate)<-ADstate$snpID 
rownames(mis_ihs)<-paste0("snp_",mis_ihs$CHROM,"_",mis_ihs$POS)

ADstate[rownames(mis_ihs),"iHS"]<-mis_ihs$iHS
ADstate[rownames(mis_ihs),"Pvalue"]<-mis_ihs$Pvalue
ADstate[rownames(mis_ihs),"FREQ_ancestral"]<-mis_ihs$FREQ_REF
ADstate[rownames(mis_ihs),"FREQ_derived"]<-1-mis_ihs$FREQ_REF
ADstate[rownames(mis_ihs),"IHH_ancestral"]<-mis_ihs$iHH_REF
ADstate[rownames(mis_ihs),"IHH_derived"]<-mis_ihs$iHH_ALT

# reorganize dataset
ADstate<-ADstate[,c("snpID","CHROM","POS","geneID","REF","ALT","FREQ_a","a_is","FREQ_ancestral","FREQ_derived","IHH_ancestral","IHH_derived","iHS","Pvalue")] 

# Now it is a bit complex
# We need to attribute the iHHa value to the correct allele, to determine, according to the sign of iHHd/iHHa which allele has the highest iHH
# btw, all of this would not be necessary if we had a clear understanding on how rehh decide who is the "a" who is the "d"
# if "a" is "ALT" (based on freqs), and "a" is derived, then derived is ALT and iHHd stand for iHH of iHH_alternative, and if iHHancestral > iHHderived, then ref allele is selected.
# There are 7 more cases, those are the 8 paths of a the decision tree:
P1<-which(ADstate$a_is == "ALT" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)<abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral > ADstate$IHH_derived )
P2<-which(ADstate$a_is == "ALT" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)<abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral < ADstate$IHH_derived )
P3<-which(ADstate$a_is == "ALT" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)>abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral > ADstate$IHH_derived )
P4<-which(ADstate$a_is == "ALT" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)>abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral < ADstate$IHH_derived )
P5<-which(ADstate$a_is == "REF" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)<abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral > ADstate$IHH_derived )
P6<-which(ADstate$a_is == "REF" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)<abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral < ADstate$IHH_derived )
P7<-which(ADstate$a_is == "REF" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)>abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral > ADstate$IHH_derived )
P8<-which(ADstate$a_is == "REF" & abs(ADstate$FREQ_a-ADstate$FREQ_derived)>abs(ADstate$FREQ_a-ADstate$FREQ_ancestral) & ADstate$IHH_ancestral < ADstate$IHH_derived )

ref<-c(P1,P4,P6,P7)
alt<-c(P2,P3,P5,P8)
# Now print the allele with the highest iHH ...
ADstate$selected_allele<-NA 
if(length(ref)>0){ ADstate$selected_allele[ref]<-ADstate$REF[ref] } 
if(length(alt)>0){ ADstate$selected_allele[alt]<-ADstate$ALT[alt] }   
#  ... and filter by the Pvalue. Note that the Pvalue reported by ihh2ihs is -log10-transformed
ADstate$selected_allele[which(ADstate$Pvalue<(-log10(0.1)))]<-NA

# Take advantage of that to print ancestral alleles
ADstate$ANC<-NA 
ADstate$ANC[c(P1,P2,P7,P8)]<-ADstate$REF[c(P1,P2,P7,P8)]
ADstate$ANC[c(P3,P4,P5,P6)]<-ADstate$ALT[c(P3,P4,P5,P6)]

# reorganize dataset
ADstate<-ADstate[,c(1:6,16,9:15)]

# Write
write.csv2(x = ADstate , file = paste0("$PATH-TO-SERVER/private/Athaliana_flowersize/ADstate_iHS_",argum[1],".csv"),quote = F,row.names = F)
