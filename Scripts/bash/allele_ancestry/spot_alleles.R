#!/usr/bin/env Rscript
argum=commandArgs(trailingOnly=T)
#argum=c("AT1G03670","gene_list.aa")

##########
# A. thaliana Petal size variation
# Whole genome REHH computation
#     identify alleles for the outgroups
# KS - 2023-09-22
##########

# Allele counts for A.tha
library(vcfR)
if(file.exists(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2],".files/",argum[1],".recode.vcf"))) {
  
  vcf<-read.vcfR(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2],".files/",argum[1],".recode.vcf"))
  
  alt<-data.frame(vcf@fix)[,c(2,4,5)]
  alt$POS<-as.numeric(alt$POS)
  #table of genotypes
  gt<-matrix(vcfR::extract.gt(vcf),nrow = dim(alt)[1])
  
  Y<-strsplit(alt$ALT,split = ",")
  # initiate new columns in ALT
  X<-rep(NA,length(alt$ALT))
  alt<-data.frame(alt,ALT1=X,ALT2=X,ALT3=X,ALT4=X,ALT5=X)
  # input new columns
  for (i in 1:length(alt$ALT)) {
    for (j in 1:5) {
      alt[i,3+j]<-Y[[i]][j]
    }
  }
  # remove the former "ALT" column 
  alt<-alt[,-3]
  # split maternal paternal alleles
  allm<-substr(gt,1,1)
  allp<-substr(gt,3,3)
  #replace number by nucleic acid
  for (i in 1:dim(gt)[1] ) {
    if("0" %in% allm[i,]){allm[i,which(allm[i,]=="0")]<-alt[i,2]} 
    if("1" %in% allm[i,]){allm[i,which(allm[i,]=="1")]<-alt[i,3]} 
    if("2" %in% allm[i,]){allm[i,which(allm[i,]=="2")]<-alt[i,4]} 
    if("3" %in% allm[i,]){allm[i,which(allm[i,]=="3")]<-alt[i,5]} 
    if("4" %in% allm[i,]){allm[i,which(allm[i,]=="4")]<-alt[i,6]} 
    if("0" %in% allp[i,]){allp[i,which(allp[i,]=="0")]<-alt[i,2]} 
    if("1" %in% allp[i,]){allp[i,which(allp[i,]=="1")]<-alt[i,3]} 
    if("2" %in% allp[i,]){allp[i,which(allp[i,]=="2")]<-alt[i,4]} 
    if("3" %in% allp[i,]){allp[i,which(allp[i,]=="3")]<-alt[i,5]} 
    if("4" %in% allp[i,]){allp[i,which(allp[i,]=="4")]<-alt[i,6]} 
  }
  alld<-rbind(data.frame(t(allm)),data.frame(t(allp)))
  
  for (i in 1:dim(gt)[1] ) {
    alld[,i]<-factor(x = alld[,i],levels = c("A","C","G","T"))
    if (i==1) {
      focal<-table(alld[,i])
    }else{
        focal<-rbind(focal,table(alld[,i]))
                     } 
  }

}

resfocal<-data.frame(focal=paste0(focal[,1],",",focal[,2],",",focal[,3],",",focal[,4]))

# outgroup alleles
library(pegas,warn.conflicts = F)
if(file.exists(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2],".files/",argum[1],".gff"))) {
  
  gff<-read.table(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2],".files/",argum[1],".gff"))
  # note: +1 because bedtools behavior
  for (i in 1:dim(gff)[1] ) {   
    if (i==1) {
      pos<-seq(gff$V4[i]+1 ,gff$V5[i],1)
    }else{
      pos<-c(pos,seq(gff$V4[i]+1 ,gff$V5[i],1))
    }  
  }
  
  if (gff$V7[1]=="-") { pos<-rev(pos) }
  
  fasta<-read.FASTA(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2] ,".files/","mafft_",argum[1]))
  names(fasta)<-c("atha","alyr","ahal")
  fasta<-data.frame(atha=fasta$atha,alyr=fasta$alyr,ahal=fasta$ahal)
  
  fasta_atcg<-as.matrix(fasta)
  fasta_atcg[which(fasta_atcg=="88")]<-"A"
  fasta_atcg[which(fasta_atcg=="18")]<-"T"
  fasta_atcg[which(fasta_atcg=="48")]<-"G"
  fasta_atcg[which(fasta_atcg=="28")]<-"C"
  fasta_atcg[which(fasta_atcg=="04")]<-"-"
  
  fasta<-data.frame(fasta_atcg)
  fasta$pos<-"NA" 
  nucaci<-length(which(fasta$atha=="-"))
  if(nucaci>0){
    fasta$pos[-which(fasta$atha=="-")]<-pos
  }else{fasta$pos<-pos}  
  
  snppos<-vcf@fix[,2]
  
  res<-fasta[which(fasta$pos %in% snppos),2:4]
  
for (i in 1:dim(res)[1] ) {
  resi<-res[i,] 
  resout<-data.frame(matrix(data = 0,nrow =2,ncol = 4,byrow = T,dimnames = list(c("alyr","ahal"),c("A","C","G","T"))))
  for (j in 1:2) {
    resout[rownames(resout)[j],resi[,rownames(resout)[j]]]<-1
  }
  
  outgroup1<-paste0(resout[1,1] ,",",resout[1,2],",",resout[1,3],",",resout[1,4] )
  outgroup2<-paste0(resout[2,1] ,",",resout[2,2],",",resout[2,3],",",resout[2,4] )

  if (i==1) {
    outgroups<-data.frame(outgroup1=outgroup1,outgroup2=outgroup2)
  }else{
    outgroups<-rbind(outgroups,data.frame(outgroup1=outgroup1,outgroup2=outgroup2))
  } 
}
  
  result<-cbind(resfocal,outgroups)
  mapped<-grep("1",paste(result[,2],result[,3],sep = ","))
  snpslist<-data.frame(snpID=paste0("snp_",substr(vcf@fix[,1],4,4),"_",vcf@fix[,2]),geneID=argum[1])
 
# remove snp if none of the outgroup has an allele
  if (length(mapped)>0) {
    result<-result[mapped,]
    snpslist<-snpslist[mapped,]
    
  resultall<-read.table(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2] ,".files/","data-file.",argum[2] ,".txt"),h=T,sep="\t")
  resultall<-rbind(resultall,result)
  
  snpslistall<-read.table(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2] ,".files/","data-file.snps.",argum[2] ,".txt"),h=T,sep="\t")
  snpslistall<-rbind(snpslistall,snpslist)
  
  write.table(x = resultall,file = paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2] ,".files/","data-file.",argum[2] ,".txt"),sep = "\t",quote = F,row.names = F,col.names = T)
  write.table(x = snpslistall,file = paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",argum[2] ,".files/","data-file.snps.",argum[2] ,".txt"),sep = "\t",quote = F,row.names = F,col.names = T)
  
  }  
}


rm(list = ls())


