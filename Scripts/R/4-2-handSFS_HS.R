##########
# A. thaliana Petal size variation
# Whole genome REHH computation
#     identify alleles for the outgroups
# KS - 2023-09-25
##########

DAF1135<-read.table(file = "../large_files/Ath_Petal_size/sfs/Genome_wide_DAF.txt",h=T)

hsqt<-c(paste0("HS0",1:9),"HS10")

freq_class<-c(seq(.1,.5,.05),1)
#freq_class<-c(seq(.1,.95,.1),1)

handtable<-matrix(data = NA,nrow = length(freq_class),ncol = length(hsqt),byrow = T,dimnames = list(hsqt,freq_class))


for (i in 1:10) {
  
freq<-read.table(paste0("../large_files/Ath_Petal_size/sfs/accessions_",hsqt[i],".txt.frq"),h=F,fill=T,skip = 1)
head(freq)
colnames(freq)<-c("chrm","pos","N_alleles","N_chr","all1","all2")
freq$nuca1<-substr(freq$all1,1,1)
freq$freq1<-as.numeric(sub('..','',freq$all1))
freq$nuca2<-substr(freq$all2,1,1)
freq$freq2<-as.numeric(sub('..','',freq$all2))
freq$snpID<-paste0("snp_",freq$chrm,"_",freq$pos)

DAF<-merge(x = DAF1135,y = freq,by = "snpID",all.x = T)
DAF$HSfreq<-NA 
DAF$HSfreq[which(DAF$DER == DAF$nuca1)]<-DAF$freq1[which(DAF$DER == DAF$nuca1)]
DAF$HSfreq[which(DAF$DER == DAF$nuca2)]<-DAF$freq2[which(DAF$DER == DAF$nuca2)]

DAF$HSDAF<-freq_class[10]
for (j in 9:1) {
DAF$HSDAF[which(DAF$HSfreq<freq_class[j])]<-freq_class[j]
}

ggtable.t<-data.frame(counts=as.vector(table(DAF$HSDAF)),freq=names(table(DAF$HSDAF)),HS=rep(hsqt[i],length(table(DAF$HSDAF))))
if(i==1){ggtable<-ggtable.t}else{ggtable<-rbind(ggtable,ggtable.t)}
handtable[i,]<-as.vector(table(DAF$HSDAF))
}

library(ggplot2)
ggplot(ggtable, aes(fill=HS,y=counts,x=freq))+
  geom_bar(position="dodge",stat="identity")

shade<-gray.colors(10,start = .8,end = .2,gamma = 1)
barplot(height = handtable[,1:9],beside = T,col=shade,las=1,main = "Allele counts per frequency")
legend(x = 80,y = 4e+05,legend = hsqt,fill = shade)

write.table(x = handtable,file = "Genetics/sfs/DAF_table_rel.95_.1to.5.txt",quote = F)
write.table(x = ggtable,file = "Genetics/sfs/DAF_table_rel.95_.1to.5.gg.txt",quote = F)
