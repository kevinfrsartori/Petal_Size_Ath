#!/usr/bin/env Rscript
rm(list = ls())
args = commandArgs(trailingOnly = T)
#args<-c("flct_Hits_flowering_time","$PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/annotations/")#test
trait<-paste0(strsplit(args[1],split = "_")[[1]][3],"_",strsplit(args[1],split = "_")[[1]][4])

####################
# SNPEFF A. thaliana 
# Petal Area project
# 2022-05-08 updated 2023-01
####################
# extract info regarding hit SNPs detected by GWAs

library(vcfR)
vcf<-read.vcfR(paste0(args[2],args[1] ,".ann.vcf"))
#vcf<-read.vcfR(paste0("$PATH-TO-SERVER/A_thaliana/admixture/Ushape/annotated_top100_",args[1] ,".ann.vcf"))
#vcf<-read.vcfR(paste0("$PATH-TO-SERVER/A_thaliana/admixture/Ushape/annotated_full_snplist.ann.vcf"))

anno<-data.frame(vcf@fix)
anno$INFO<-as.character(anno$INFO)
full_anno<-data.frame(ALT=NA,effect=NA,impact=NA,gene_name=NA,gene_ID=NA,feature_type=NA,feature_ID=NA,transcript_biotype=NA,
                      rank=NA,HGVS_c=NA,HGVS_p=NA,cDNA_position=NA,CDS_position=NA,protein_postion=NA,distance_to_feature=NA,warnings=NA,snp=NA)
anno[,8]<-unlist(lapply(strsplit(anno[,8],"="), `[[`, 2))
anno2<-strsplit(anno[,8], ",")

for (i in 1:length(anno2)) {
  anno3<-strsplit(anno2[[i]],"\\|")
  
  anno_temp<-data.frame(matrix(data = NA,nrow = length(anno3),ncol = 16,byrow = T))
  colnames(anno_temp)<-colnames(full_anno)[1:16]
  for (j in 1:length(anno3)) {
    anno_temp[j,1:length(anno3[[j]])]<-  anno3[[j]]
  }
  anno_temp$snp<-rep(anno$ID[i],length(anno3) )
  anno_temp[anno_temp == ""] <-NA
  full_anno<-rbind(full_anno,anno_temp)
}
rm(anno_temp,anno2,anno3)
full_anno<-full_anno[-1,]

#remove intergenic
if(length(grep("intergenic",full_anno$feature_type))>0){ full_anno<-full_anno[-grep("intergenic",full_anno$feature_type ),] } 

#remove problematic genes
if(length(grep("WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS",full_anno$warnings))>0){ full_anno<-full_anno[-grep("WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS",full_anno$warnings ),] }
if(length(grep("WARNING_TRANSCRIPT_NO_START_CODON",full_anno$warnings))>0){ full_anno<-full_anno[-grep("WARNING_TRANSCRIPT_NO_START_CODON",full_anno$warnings ),] }

# group snps per LD blocks
full_anno$chrm<-as.factor(unlist(lapply(strsplit(full_anno$snp,split = "_"),'[[',2)))
full_anno$pos<-as.numeric(unlist(lapply(strsplit(full_anno$snp,split = "_"),'[[',3)))
full_anno$block<-NA 
for (chrm in levels(full_anno$chrm )) {
  anno_temp<-full_anno[which(full_anno$chrm==chrm),]
  anno_temp<-anno_temp[order(anno_temp$pos),]
  for (i in 1:length(anno_temp$pos)) {
    if (i==1) {
      j<-1
      block<-paste0(chrm,".",j)
      anno_temp$block[i]<-block 
    }else if( (anno_temp$pos[i]-anno_temp$pos[i-1])<20000){
      anno_temp$block[i]<-block
    }else{
        j<-j+1
        block<-paste0(chrm,".",j)
        anno_temp$block[i]<-block
      }  
  }
  full_anno[which(full_anno$chrm==chrm),]<-anno_temp
}

full_anno$block<-as.factor(full_anno$block)

# remove splice variants
for (block in levels(full_anno$block)) {
  # proceed per block
  anno_temp<-full_anno[which(full_anno$block==block),]
  # if block more than one line
  if (dim(anno_temp)[1]>1){
  # keep only geneID.1
  if(length(grep("\\.1",anno_temp$feature_ID))>0){anno_temp<-anno_temp[grep("\\.1",anno_temp$feature_ID),]}
    #remove pseudogene too
    if("pseudogene" %in% anno_temp$transcript_biotype ){anno_temp<-anno_temp[-which(anno_temp$transcript_biotype=="pseudogene"),] }
    if("lincRNA" %in% anno_temp$transcript_biotype ){anno_temp<-anno_temp[-which(anno_temp$transcript_biotype=="lincRNA"),] }
    
  }
  # make new file
  if (!exists("full_anno.1")) { full_anno.1<-anno_temp }else{full_anno.1<-rbind(full_anno.1,anno_temp)} 
}

full_anno.1$gene_name<-gsub(pattern = "'",replacement = ".",x = full_anno.1$gene_name)
full_anno.1$gene_name<-gsub(pattern = "%",replacement = ".",x = full_anno.1$gene_name)

trait
print(paste0(length(unique(full_anno.1$gene_ID)), " genes"))                                

write.table(x = full_anno.1[,c(17,18,19,5,4,3,2,1)],file = paste0(args[2],"annotated_simple_",args[1],".txt"),quote = F,row.names = F,col.names = F,sep = "\t")

print("")
print("")

