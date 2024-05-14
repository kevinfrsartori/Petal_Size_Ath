#!/usr/bin/env Rscript

##########
# A. thaliana Petal size variation
# Whole genome REHH computation
#     identify alleles for the outgroups
# KS - 2023-09-25
##########
args = commandArgs(trailingOnly = T)
trait=args[1] 
#trait="Petal"

# load ancestral state data from rehh analysis
AD<-read.table(file = "$PATH-TO-SERVER/Allele_ancestry/Genome_wide_iHS.txt",h=T)

# habitat suitability range names
hsqt<-c(paste0("HS0",1:9),"HS10")

# frequency bins to use for the plot:
freq_class<-c(seq(.1,.5,.05),1)

# Prepare result table
handtable<-matrix(data = NA,nrow = length(freq_class),ncol = length(hsqt),byrow = T,dimnames = list(hsqt,freq_class))

# for each HS range
for (i in 1:10){
  
freq<-read.table(paste0(trait,".accessions_",hsqt[i],".txt.frq"),h=F,fill=T,skip = 1)
colnames(freq)<-c("chrm","pos","N_alleles","N_chr","all1","all2")
freq$chrm<-substr(freq$chrm,4,4)
freq$snpID<-paste0("snp_",freq$chrm,"_",freq$pos)
if(length(which(duplicated(freq$snpID)))>0){freq<-freq[-which(duplicated(freq$snpID)),]}

freq$nuca1<-substr(freq$all1,1,1)
freq$freq1<-as.numeric(sub('..','',freq$all1))
freq$nuca2<-substr(freq$all2,1,1)
freq$freq2<-as.numeric(sub('..','',freq$all2))

ADi<-merge(x = AD,y = freq,by = "snpID",all.y = T)
ADi<-na.omit(ADi)
ADi<-ADi[-which(duplicated(ADi$snpID)),]
ADi$HSfreq<-NA 
ADi$HSfreq[which(ADi$DER == ADi$nuca1)]<-ADi$freq1[which(ADi$DER == ADi$nuca1)]
ADi$HSfreq[which(ADi$DER == ADi$nuca2)]<-ADi$freq2[which(ADi$DER == ADi$nuca2)]

ADi$HSDAF<-freq_class[10]
for (j in 9:1) {
  ADi$HSDAF[which(ADi$HSfreq<freq_class[j])]<-freq_class[j]
}

ggtable.t<-data.frame(counts=as.vector(table(ADi$HSDAF)),freq=names(table(ADi$HSDAF)),HS=rep(hsqt[i],length(table(ADi$HSDAF))))
if(i==1){ggtable<-ggtable.t}else{ggtable<-rbind(ggtable,ggtable.t)}
handtable[i,]<-as.vector(table(ADi$HSDAF))
}


svglite::svglite(filename = paste0("unfolded_SFS_",trait,".svg"),width = 4.5, height = 5.5)
hscol<-colorRampPalette(rev(viridis::inferno(4)))
barplot(height = handtable[,1:4],beside = T,col=hscol(20)[6:15],las=1,
        main = paste0(trait," genes (", dim(ADi)[1], " SNPs)" ),ylim=c(0,1000),
        names.arg = paste0("<",colnames(handtable[,1:4])),
        xlab = "Allele frequency",ylab="Allele count")
dev.off()

write.table(x = handtable, file = paste0("unfolded_SFS_",trait,".txt"), quote = F, row.names = T, col.names = T)
