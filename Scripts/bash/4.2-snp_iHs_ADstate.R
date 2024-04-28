#!/usr/bin/env Rscript
rm(list = ls())
args = commandArgs(trailingOnly = T)
#args<-c("flct_Hits_flowering_time","/crex/proj/snic2020-16-182/A_thaliana/GWAs/Petal_Size_Ath/annotations/")
trait<-paste0(strsplit(args[1],split = "_")[[1]][3],"_",strsplit(args[1],split = "_")[[1]][4])

####################
# iHS and Allele ancestry 
# Petal Area project
# 2022-05-08 updated 2023-01
####################

# This holds for snps located in a gene only, meaning that we could compute ancestry state

ihs<-read.table(file = paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/snp_lists/",args[1],".iHS.AD.txt"),h=T )
genelist<-read.table(paste0(args[2],"annotated_simple_",args[1],".txt"))
                
genelist$label<-paste0(genelist$V1,"_",genelist$V4)
ihs$label<-paste0(ihs$snpID,"_",ihs$geneID)

print("Gene number:")
length(unique(genelist$V4))

aprx_match<-which(genelist$V1 %in% ihs$snpID)
print("Gene with ancestry number:")
length(unique(genelist[aprx_match,]$V5))

perf_match<-which(genelist$label %in% ihs$label)
print("Gene with exact ancestry number:")
length(unique(genelist[perf_match,]$V5))

aprx_match<-aprx_match[-which(aprx_match %in% perf_match)]
# j<-1
# i<-perf_match[j]
# ihs[which(ihs$label==genelist$label[i]),]
# j<-j+1
# 
# j<-1
# i<-aprx_match[j]
# ihs[which(ihs$snpID==genelist$V1[i]),]
# j<-j+1

# OK now make dataset...

genelist_pm<-genelist[perf_match,]
genelist_pm<-merge(genelist_pm,ihs[,5:13],by="label",all.x=T)[,-1]

genelist_am<-genelist[aprx_match,]
for (i in 1:length(aprx_match)) {
  snp<-genelist_am$V1[i]
  
  gl_dat<-genelist_am[i,1:8]
  
  ihs.t<-ihs[which(ihs$snpID==snp),]
if(mean(ihs.t$iHS)==ihs.t$iHS[1]){   ihs_dat<-ihs.t[1,5:12] }else{ihs_dat<-rep(NA,8)}
  
  genelist_pm[dim(genelist_pm)[1]+1,]<-c(gl_dat,ihs_dat)
}

colnames(genelist_pm)[1:8]<-c("snpID","chr","ps","geneID","genename","effect","type","ALT")
colnames(genelist_pm)[16]<-"iHS_pval"

genelist_r<-genelist[-c(aprx_match,perf_match),]

for (i in 1:dim(genelist_r)[1]) {
  snp<-genelist_r$V1[i]
  gl_dat<-genelist_r[i,1:8]
  ihs_dat<-rep(NA,8)
  genelist_pm[dim(genelist_pm)[1]+1,]<-c(gl_dat,ihs_dat)
}

write.table(x = genelist_pm, file = paste0(args[2],"annotated_simple_",args[1],".iHS.AD.txt"),quote = F,sep = "\t",row.names = F,col.names = T)
