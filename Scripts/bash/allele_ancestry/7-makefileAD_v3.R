##########
# A. thaliana Petal size variation
# Whole genome REHH computation
#     identify alleles for the outgroups
# KS - 2023-09-25
##########

# load iHH data
print("load iHH data")
for (i in 1:5) {
print(paste0(i,"/5"))
  if (i==1) {
    chrm<-read.table(paste0("/crex/proj/snic2020-16-182/A_thaliana/fastphase/scanEHH.chrm-",i,".txt"),h=T,sep=" ",dec=".")
  }else{
    chrm<-rbind(chrm,read.table(paste0("/crex/proj/snic2020-16-182/A_thaliana/fastphase/scanEHH.chrm-",i,".txt"),h=T,sep=" ",dec="."))
  }
}


head(chrm)
chrm$snpID<-paste0("snp_",chrm$CHR,"_",chrm$POSITION)
length(unique(chrm$snpID))

# the dataset mentions IHHa for ancestral allele and IHHd for derived allele but we don't have computed the ancestry state yet.
# the idea is to make new dataset with IHH_minor for minor allele and IHH_major for major allele
# and define the ancestry afterwards
hist(chrm$FREQ_a)

chrm$IHHminor[which(chrm$FREQ_a<=.5)]<-chrm$IHHa[which(chrm$FREQ_a<=.5)]
chrm$IHHminor[which(chrm$FREQ_a>.5)]<-chrm$IHHd[which(chrm$FREQ_a>.5)]

chrm$IHHmajor[which(chrm$FREQ_a<=.5)]<-chrm$IHHd[which(chrm$FREQ_a<=.5)]
chrm$IHHmajor[which(chrm$FREQ_a>.5)]<-chrm$IHHa[which(chrm$FREQ_a>.5)]

chrm$FREQ_minor[which(chrm$FREQ_a<=0.5)]<-chrm$FREQ_a[which(chrm$FREQ_a<=0.5)]
chrm$FREQ_minor[which(chrm$FREQ_a>0.5)]<-(1-chrm$FREQ_a[which(chrm$FREQ_a>0.5)])

#load AD data
print("load AD data")
#genelist.f<-read.table("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/gene_lists.txt")[,1] 
genelist<-substr(grep("output-file-pvalues",list.files("/crex/proj/snic2020-16-182/software/est_sfs/"),value = T),21,32)
#rerun<-genelist.f[-which(genelist.f %in% genelist)]
#write.table(x = rerun, file = "rerun.txt",quote = F,col.names = F,row.names = F)
for (i in 1:length(genelist)) {
 print(paste0(i,"/",length(genelist)))
  estsfs<-read.table(paste0("/crex/proj/snic2020-16-182/software/est_sfs/output-file-pvalues.",genelist[i],".txt"),skip=8)[,3]
  estprob<-read.table(paste0("/crex/proj/snic2020-16-182/software/est_sfs/data-file.",genelist[i],".txt"),h=F)
  estsnp<-read.table(paste0("/crex/proj/snic2020-16-182/A_thaliana/temp/temp/ancestry447/",genelist[i],".files/data-file.snps.",genelist[i],".txt"),h=T)
  est.t<-data.frame(chrm=unlist(lapply(strsplit(estsnp$snpID,split = "_" ),"[[",2)),
                  pos=unlist(lapply(strsplit(estsnp$snpID,split = "_" ),"[[",3)),
                  snpID=estsnp$snpID, geneID=estsnp$geneID,
                  P_major_ancestral=estsfs)
  est.t<-cbind(est.t,estprob)
  if (i==1) {
    est<-est.t
  }else{
    est<-rbind(est,est.t)
  }
  rm(est.t,estprob,estsnp,estsfs)
}
head(est)
length(unique(est$snpID))
#! duplicated snps
x<-est$snpID[which(duplicated(est$snpID))][3]
est[which(est$snpID==x),]
# they belong to different genes
colnames(est)[6:8]<-c("focal","outgroup1","outgroup2")
# split nucleic acid info per species
est$AT_A<-as.numeric(unlist(lapply(strsplit(est$focal,split = ","),"[[",1)))
est$AT_C<-as.numeric(unlist(lapply(strsplit(est$focal,split = ","),"[[",2)))
est$AT_G<-as.numeric(unlist(lapply(strsplit(est$focal,split = ","),"[[",3)))
est$AT_T<-as.numeric(unlist(lapply(strsplit(est$focal,split = ","),"[[",4)))

# est$AL_A<-as.numeric(unlist(lapply(strsplit(est$outgroup1,split = ","),"[[",1)))
# est$AL_C<-as.numeric(unlist(lapply(strsplit(est$outgroup1,split = ","),"[[",2)))
# est$AL_G<-as.numeric(unlist(lapply(strsplit(est$outgroup1,split = ","),"[[",3)))
# est$AL_T<-as.numeric(unlist(lapply(strsplit(est$outgroup1,split = ","),"[[",4)))
# est$AH_A<-as.numeric(unlist(lapply(strsplit(est$outgroup2,split = ","),"[[",1)))
# est$AH_C<-as.numeric(unlist(lapply(strsplit(est$outgroup2,split = ","),"[[",2)))
# est$AH_G<-as.numeric(unlist(lapply(strsplit(est$outgroup2,split = ","),"[[",3)))
# est$AH_T<-as.numeric(unlist(lapply(strsplit(est$outgroup2,split = ","),"[[",4)))

# remove alleles with 0 everywhere
rsum<-apply(est[,c("AT_A","AT_C","AT_G","AT_T")],1,sum,na.rm=T)
torm<-which(rsum==0)
if (length(torm>0)) { est<-est[-torm,]}

# Replace 0 by NA in order to spot the minor allele frequence
if(any(est$AT_A==0)){est$AT_A[which(est$AT_A==0)]<-NA }
if(any(est$AT_C==0)){est$AT_C[which(est$AT_C==0)]<-NA }
if(any(est$AT_G==0)){est$AT_G[which(est$AT_G==0)]<-NA }
if(any(est$AT_T==0)){est$AT_T[which(est$AT_T==0)]<-NA }

# vector of nucleic acids in order
acid<-c("A","C","G","T")
# Spot the min and maj alleles
print("Spot the min and maj alleles")
est$MAJ<-NA
est$MIN<-NA
steps<-c(seq(1,dim(est)[1],1000),dim(est)[1])
for(i in 2:length(steps)){
est$MAJ[steps[i-1]:steps[i]]<-acid[apply(est[steps[i-1]:steps[i],c("AT_A","AT_C","AT_G","AT_T")],1,which.max)]
est$MIN[steps[i-1]:steps[i]]<-acid[apply(est[steps[i-1]:steps[i],c("AT_A","AT_C","AT_G","AT_T")],1,which.min)]
print(i)
}



# Spot the ancestral and derived allele
print("pot the ancestral and derived allele")
est$ANC[which(est$P_major_ancestral>.5)]<-est$MAJ[which(est$P_major_ancestral>.5)]
est$DER[which(est$P_major_ancestral>.5)]<-est$MIN[which(est$P_major_ancestral>.5)]
est$ANC[which(est$P_major_ancestral<=.5)]<-est$MIN[which(est$P_major_ancestral<=.5)]
est$DER[which(est$P_major_ancestral<=.5)]<-est$MAJ[which(est$P_major_ancestral<=.5)]

head(est)
l_est<-length(est$chrm)
estihh<-merge(est[,-c(6,7,8)] ,chrm[,7:10],by = "snpID",all = F)
head(estihh)
# estihh[which(estihh$snpID==x),]

# chrm$IHHminor[which(chrm$FREQ_a<=.5)]<-chrm$IHHa[which(chrm$FREQ_a<=.5)]

estihh$IHH_ancestral[which(estihh$P_major_ancestral>.5)]<-estihh$IHHmajor[which(estihh$P_major_ancestral>.5)]
estihh$IHH_ancestral[which(estihh$P_major_ancestral<=.5)]<-estihh$IHHminor[which(estihh$P_major_ancestral<=.5)]
estihh$IHH_derived[which(estihh$P_major_ancestral>.5)]<-estihh$IHHminor[which(estihh$P_major_ancestral>.5)]
estihh$IHH_derived[which(estihh$P_major_ancestral<=.5)]<-estihh$IHHmajor[which(estihh$P_major_ancestral<=.5)]

head(estihh)

estihh$FREQ_ancestral[which(estihh$P_major_ancestral>.5)]<-1-estihh$FREQ_minor[which(estihh$P_major_ancestral>.5)]
estihh$FREQ_ancestral[which(estihh$P_major_ancestral<=.5)]<-estihh$FREQ_minor[which(estihh$P_major_ancestral<=.5)]
estihh$FREQ_derived[which(estihh$P_major_ancestral>.5)]<-estihh$FREQ_minor[which(estihh$P_major_ancestral>.5)]
estihh$FREQ_derived[which(estihh$P_major_ancestral<=.5)]<-1-estihh$FREQ_minor[which(estihh$P_major_ancestral<=.5)]

# estihh[which(estihh$snpID==x),]

# remove singleton
f<-which(estihh$FREQ_derived<=1/1135)
estihh<-estihh[-f,] 
# relation frequency ~ iHH
# par(mfrow=c(1,2))
# plot(estihh$IHH_ancestral ~ estihh$FREQ_ancestral,pch=16,col=rgb(0,.5,.4,.3),log="y",ylim=c(10,1000000))
# plot(estihh$IHH_derived ~ estihh$FREQ_derived,pch=16,col=rgb(0,.4,.5,.3),log="y",ylim=c(10,1000000))

# compute the statistic ultimately
print("compute the statistic")
estihh$iHS<-log(estihh$IHH_ancestral/estihh$IHH_derived)
# par(mfrow=c(1,1))
# plot(estihh$iHS ~ estihh$FREQ_derived,pch=16,col=rgb(0,.4,.5,.3),ylab="unstandardized iHS",xlab="Derived allele frequency")

# Write table
write.table(x = estihh[,-c(6:9)],file = "Genome_wide_uniHS.txt",quote = F,sep = "\t",row.names = F,col.names = T)

# Compute standardized iHS
library(rehh)
scanhh<-estihh[,c("chrm","pos","FREQ_ancestral","IHH_ancestral","IHH_derived")]
minmaf=0
scanhh = scanhh[scanhh[, 3] > minmaf & scanhh[, 3] < (1 - minmaf), ]
estihs<-ihh2ihs(res_ihh = scanhh,minmaf = minmaf)[[1]]
 
# merge datasets
colnames(estihh)[dim(estihh)[2]]<-"uniHS"
# estihs$snpID<-paste0("snp_",estihs$chrm,"_",estihs$pos) 
# estihh<-merge(estihh,estihs[,-c(1,2)],by = "snpID",all = F)
estihh$iHS<-estihs$iHS
estihh$Pvalue<-estihs$Pvalue
# plot(estihh$iHS ~ estihh$FREQ_derived,pch=16,col=rgb(0,.4,.5,.3),ylab="standardized iHS",xlab="Derived allele frequency")

# estihh[which(estihh$snpID==x),]

# Write table
estihh<-estihh[,c("snpID","chrm","pos","geneID","ANC","DER","IHH_ancestral","IHH_derived","FREQ_ancestral","FREQ_derived","iHS","Pvalue")]
write.table(x = estihh,file = "Genome_wide_iHS.txt",quote = F,sep = "\t",row.names = F,col.names = T)


# Manhattan plot
##estihh<-read.table("Genome_wide_iHS.txt",h=T)
#range(estihh$Pvalue)
#length(unique(estihh$snpID))
#-log10(0.05/length(unique(estihh$snpID)))
# estihh$chrm<-as.numeric(estihh$chrm)
# estihh$pos<-as.numeric(estihh$pos)
# estihh$rpos<-estihh$pos  
# for (j in 2:5) {
#   estihh$rpos[estihh$chrm==j]<-max(estihh$rpos[estihh$chrm==j-1])+estihh$rpos[estihh$chrm==j]
# }
# estihh$rpos<-estihh$rpos/1000000   
# palette(c("grey35","grey65","grey25","grey80","grey50"))
# estihh$Pvalue[which(estihh$Pvalue<2)]<-NA  
# plot(estihh$rpos,estihh$Pvalue, cex=.5,pch=16,col=estihh$chrm,ylab="-log10(P-value iHS)",
#      xlab="",las=1,main="Manhattan plot of standardized iHS P-values")
# abline(h=-log10(0.05/dim(estihh)[1]),lty=1,col="black",lwd=1)
# axis(side = 1,at = 60,labels = "Position (Mb)",line = 1.5,tick = F)
# estihh$Pvalue[order(estihh$Pvalue,decreasing = T)][100]  
# 
# estihh$snpID[order(estihh$Pvalue,decreasing = T)][1]  
