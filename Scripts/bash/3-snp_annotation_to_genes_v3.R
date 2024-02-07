#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)
#args<-c("top100_Petal_Area","/crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/")
trait<-paste0(strsplit(args[1],split = "_")[[1]][2],"_",strsplit(args[1],split = "_")[[1]][3])

####################
# SNPEFF A. thaliana 
# ushaped project
# 2022-05-08
####################
# extract info regarding top 100 SNPs detected by GWAs

# automatize the selection of the best candidate gene in a 20Kb LD window 

library(vcfR)
vcf<-read.vcfR(paste0(args[2],args[1] ,".ann.vcf"))
#vcf<-read.vcfR(paste0("/crex/proj/snic2020-16-182/A_thaliana/admixture/Ushape/annotated_top100_",args[1] ,".ann.vcf"))
#vcf<-read.vcfR(paste0("/crex/proj/snic2020-16-182/A_thaliana/admixture/Ushape/annotated_full_snplist.ann.vcf"))

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

#gwas pval keep highest pval only
snps<-data.frame(snp=unique(full_anno$snp),Pval=NA)
assoc<-read.table(paste0(args[2],"/../",trait,"/output/SNP_1001g_filtered_",trait,".assoc.txt"),h=T,sep="\t",dec=".")
for (i in 1:length(snps$snp)) {
  snps$Pval[i]<-assoc$p_lrt[which(assoc$rs == snps$snp[i])]      
}
full_anno.1<-merge(full_anno.1,snps,by="snp",all.x = T,sort = F)

for (block in levels(full_anno.1$block)) {
  # proceed per block
  anno_temp<-full_anno.1[which(full_anno.1$block==block),]
  # if block more than one line
  if (dim(anno_temp)[1]>1){
    # keep snp with best pval
    snp<-anno_temp$snp[which.min(anno_temp$Pval)]
    anno_temp<-anno_temp[which(anno_temp$snp == snp),] 
  }
  # make new file
  if (!exists("full_anno.hit")) { full_anno.hit<-anno_temp }else{full_anno.hit<-rbind(full_anno.hit,anno_temp)} 
}

#Keep only the closest gene
full_anno.hit$distance_to_feature[which(is.na(full_anno.hit$distance_to_feature))]<-0  
for (block in levels(full_anno.hit$block)) {
  # proceed per block
  anno_temp<-full_anno.hit[which(full_anno.hit$block==block),]
  # if block more than one line
  if (dim(anno_temp)[1]>1){
    # keep snp with best pval
    gene_ID<-anno_temp$gene_ID[which.min(anno_temp$distance_to_feature)]
    anno_temp<-anno_temp[which(anno_temp$gene_ID == gene_ID),] 
  }
  # make new file
  if (!exists("full_anno.closest")) { full_anno.closest<-anno_temp }else{full_anno.closest<-rbind(full_anno.closest,anno_temp)} 
}

full_anno.closest<-full_anno.closest[,c(6,5,1,18,19,21,2,3,4,16)]

write.table(x = full_anno.closest,file = paste0(args[2],"annotated_simple_",args[1],".txt"),quote = F,row.names = F,col.names = F,sep = "\t")


