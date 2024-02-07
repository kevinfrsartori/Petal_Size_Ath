#------------------
# Petal Size Ath
# Supplementary figures
# 2024-01-18
#------------------

# Pi
#------------
#----------------------------------
hsranges<-c(paste0("HS0",1:9),"HS10")
#pinpis<-read.table("../large_files/Ath_Petal_size/pinpis/outputs/wholegenome_pinpis_HS01.txt",h=T,na.strings = c("NA","Inf","-Inf"))
#pinpis$log10_PinPis<-log10(pinpis$PinPis)
#pinpis$log10_PinPis[which(pinpis$log10_PinPis==-Inf)]<-NA
#hist(pinpis$log10_PinPis)


library(vioplot)
plot(0, 0, xlim=c(.5,10.5), ylim=c(-4,1), las=1,xaxt="n",yaxt="n",
     xlab="Habitat suitability ranges",ylab="PiS (log scale)")
axis(side = 1, at = 1:10, labels = hsranges, las=2)
axis(side = 2, at = c(-3,-2,-1,0,1,2), labels = c(0.001,0.01,0.1,1,10,100),las=1)

for (i in 1:10) {
  pinpis<-read.table(paste0("../large_files/Ath_Petal_size/pinpis/outputs/wholegenome_pinpis_",hsranges[i],".txt"),
                     h=T,na.strings = c("NA","Inf","-Inf"))
  pinpis$log10_Pis<-log10(pinpis$Pis)
  pinpis$log10_Pis[which(pinpis$log10_Pis==-Inf)]<-NA
  
  vioplot(pinpis$log10_Pis,add = T,at = i, 
          col = "grey85", border = "grey50", rectCol = "grey85", lineCol = "grey30")
}
title(main = "Synonymous diversity")

library(vioplot)
plot(0, 0, xlim=c(.5,10.5), ylim=c(-5,0), las=1,xaxt="n",yaxt="n",
     xlab="Habitat suitability ranges",ylab="PiN (log scale)")
axis(side = 1, at = 1:10, labels = hsranges, las=2)
axis(side = 2, at = c(-3,-2,-1,0,1,2), labels = c(0.001,0.01,0.1,1,10,100),las=1)

for (i in 1:10) {
  pinpis<-read.table(paste0("../large_files/Ath_Petal_size/pinpis/outputs/wholegenome_pinpis_",hsranges[i],".txt"),
                     h=T,na.strings = c("NA","Inf","-Inf"))
  pinpis$log10_Pin<-log10(pinpis$Pin)
  pinpis$log10_Pin[which(pinpis$log10_Pin==-Inf)]<-NA
  
  vioplot(pinpis$log10_Pin,add = T,at = i, 
          col = "grey85", border = "grey50", rectCol = "grey85", lineCol = "grey30")
}
title(main = "Non-Synonymous diversity")




annot<-read.table("Genetics/functionnal_annotation_Petal_Area.csv",h=T,sep=",",dec=".")[,c("gene_ID","gene_name","candidate")]
annot<-annot[which(annot$candidate == "yes"),]
# Select informative genes
annot<-annot[c(2,3,10,1,4,6),]
# vioplot PiS
par(mfrow=c(2,3),oma=c(2,2,0,0),mar=c(1,2,1,1))
mycol=c(rep(rgb(0,.6,.5),5),rgb(.7,.4,.3))
for(i in 1:dim(annot)[1]){
  piehsgene<-read.table(paste0("../large_files/Ath_Petal_size/pinpis/outputs/result_",annot$gene_ID[i],".1.txt"),h=T)
  piehsgene$log10_Pis<-log10(piehsgene$Pis)
  piehsgene$log10_Pis[which(piehsgene$log10_Pis==-Inf)]<-NA
  piehsgene$log10_Pis[which(piehsgene$log10_Pis==Inf)]<-NA
  vioplot(piehsgene$log10_Pis ~ piehsgene$HS,xaxt="n",las=1,col=mycol[i])
  axis(side = 1,at = 1:10,labels = rep("",10))
  title(main = annot$gene_name[i],outer = )
}
title(xlab = "Habitat suitability ranges",ylab="PiS (log10)",outer = T,line = 1,cex.lab=1)

# vioplot PiN
par(mfrow=c(2,3),oma=c(2,2,0,0),mar=c(1,2,1,1))
mycol=c(rep(rgb(0,.6,.5),5),rgb(.7,.4,.3))
for(i in 1:dim(annot)[1]){
  piehsgene<-read.table(paste0("../large_files/Ath_Petal_size/pinpis/outputs/result_",annot$gene_ID[i],".1.txt"),h=T)
  piehsgene$log10_Pin<-log10(piehsgene$Pin)
  piehsgene$log10_Pin[which(piehsgene$log10_Pin==-Inf)]<-NA
  piehsgene$log10_Pin[which(piehsgene$log10_Pin==Inf)]<-NA
  vioplot(piehsgene$log10_Pin ~ piehsgene$HS,xaxt="n",las=1,col=mycol[i])
  axis(side = 1,at = 1:10,labels = rep("",10))
  title(main = annot$gene_name[i],outer = )
}
title(xlab = "Habitat suitability ranges",ylab="PiN (log10)",outer = T,line = 1,cex.lab=1)


# PiNPiS-HS for all genes
traits<-c("Ovule_Number",colnames(read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T))[6:16])[c(4:12,1:3)]
respinpie<-data.frame(trait=traits, negative=NA, NS=NA, positive=NA)
respispie<-data.frame(trait=traits, negative=NA, NS=NA, positive=NA)
respinpispie<-data.frame(trait=traits, negative=NA, NS=NA, positive=NA)


for (j in 1:length(traits)) {

genelist<-unique(read.table(paste0("Genetics/annotated_simple_top100_",traits[j],".txt"),h=F,sep="\t",dec=".")[,1])
for (i in 1:10) {
  pinpis.t<-read.table(paste0("../large_files/Ath_Petal_size/pinpis/outputs/wholegenome_pinpis_",hsranges[i],".txt"),
                     h=T,na.strings = c("NA","Inf","-Inf"))
  pinpis.t<-pinpis.t[which(pinpis.t$Gene_ID %in% genelist),]
  pinpis.t$Gene_ID<-as.factor(pinpis.t$Gene_ID)
  pinpis.t2<-data.frame(Gene_ID=levels(pinpis.t$Gene_ID),
                        Pis=tapply(pinpis.t$Pis,pinpis.t$Gene_ID,mean,na.rm=T),
                        Pin=tapply(pinpis.t$Pin,pinpis.t$Gene_ID,mean,na.rm=T),
                        PinPis=tapply(pinpis.t$PinPis,pinpis.t$Gene_ID,mean,na.rm=T))
  
  
  pinpis.t2$log10_PinPis<-log10(pinpis.t2$PinPis)
  pinpis.t2$log10_PinPis[which(pinpis.t2$log10_PinPis %in% c(Inf,-Inf,NaN))]<-NA
  pinpis.t2$log10_Pis<-log10(pinpis.t2$Pis)
  pinpis.t2$log10_Pis[which(pinpis.t2$log10_Pis %in% c(Inf,-Inf,NaN))]<-NA
  pinpis.t2$log10_Pin<-log10(pinpis.t2$Pin)
  pinpis.t2$log10_Pin[which(pinpis.t2$log10_Pin %in% c(Inf,-Inf,NaN))]<-NA
  pinpis.t2$hsrange<-i
  
  if(i==1){ pinpis<-pinpis.t2 }else{ pinpis<-rbind(pinpis,pinpis.t2) }
  
}

N<-length(levels(as.factor(pinpis$Gene_ID)))
result<-data.frame(gene=rep(NA,N),pin_slope=rep(NA,N),pis_slope=rep(NA,N),pinpis_slope=rep(NA,N))
for (i in 1:N) {
  pinpis.t<-pinpis[which(pinpis$Gene_ID == levels(as.factor(pinpis$Gene_ID))[i]),]
  result$gene[i]<-pinpis.t$Gene_ID[1]
  
  #Pin
  if(length(na.omit(pinpis.t$log10_Pin)) > 3 ) {
  mod<-summary(lm(log10_Pin ~ hsrange, data = pinpis.t))
  if(mod$coefficients[2,4]<(0.1)){result$pin_slope[i]<-mod$coefficients[2,1]}else{result$pin_slope[i]<-0}
  }
  
  #Pis
  if(length(na.omit(pinpis.t$log10_Pis)) > 3 ) {
  mod<-summary(lm(log10_Pis ~ hsrange, data = pinpis.t))
  if(mod$coefficients[2,4]<(0.1)){result$pis_slope[i]<-mod$coefficients[2,1]}else{result$pis_slope[i]<-0}
  }
  
  if(length(na.omit(pinpis.t$log10_PinPis)) > 3 ) {
  mod<-summary(lm(log10_PinPis ~ hsrange, data = pinpis.t))
  if(mod$coefficients[2,4]<(0.1)){result$pinpis_slope[i]<-mod$coefficients[2,1]}else{result$pinpis_slope[i]<-0}
  }
  
}        

respinpie[j,]<-c(traits[j],as.vector(table(factor(x = sign(result$pin_slope),levels = c(-1,0,1)))) )
respispie[j,]<-c(traits[j],as.vector(table(factor(x = sign(result$pis_slope),levels = c(-1,0,1)))) )
respinpispie[j,]<-c(traits[j],as.vector(table(factor(x = sign(result$pinpis_slope),levels = c(-1,0,1)))) )

}


#pies
par(mar=c(0,0,0,0),oma=c(0,0,2,0))
layout(mat = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = F))
for (i in 1:dim(respie)[1]) {
  pie(as.numeric(respinpie[i,2:4]),col = c(rgb(.5,.5,1),rgb(.7,.7,.7),rgb(1,.5,.5)),labels = respinpie[i,2:4] )
  axis(side = 1,at = 0,labels = gsub(pattern = "_",replacement = " ",traits[i]),tick = F,line = -2)
}
title(main = "Pin/Pis - HS slopes",outer = T)


# Bonferroni 
hsranges<-c(paste0("HS0",1:9),"HS10")
genelist<-c("AT1G77080","AT2G32370","AT3G03580")

for (i in 1:10) {
  pinpis.t<-read.table(paste0("../large_files/Ath_Petal_size/pinpis/outputs/wholegenome_pinpis_",hsranges[i],".txt"),
                       h=T,na.strings = c("NA","Inf","-Inf"))
  pinpis.t<-pinpis.t[which(pinpis.t$Gene_ID %in% genelist),]
  pinpis.t$Gene_ID<-as.factor(pinpis.t$Gene_ID)
  pinpis.t2<-data.frame(Gene_ID=levels(pinpis.t$Gene_ID),
                        Pis=tapply(pinpis.t$Pis,pinpis.t$Gene_ID,mean,na.rm=T),
                        Pin=tapply(pinpis.t$Pin,pinpis.t$Gene_ID,mean,na.rm=T),
                        PinPis=tapply(pinpis.t$PinPis,pinpis.t$Gene_ID,mean,na.rm=T))
  
  
  pinpis.t2$log10_PinPis<-log10(pinpis.t2$PinPis)
  pinpis.t2$log10_PinPis[which(pinpis.t2$log10_PinPis %in% c(Inf,-Inf,NaN))]<-NA
  pinpis.t2$log10_Pis<-log10(pinpis.t2$Pis)
  pinpis.t2$log10_Pis[which(pinpis.t2$log10_Pis %in% c(Inf,-Inf,NaN))]<-NA
  pinpis.t2$log10_Pin<-log10(pinpis.t2$Pin)
  pinpis.t2$log10_Pin[which(pinpis.t2$log10_Pin %in% c(Inf,-Inf,NaN))]<-NA
  pinpis.t2$hsrange<-i
  
  if(i==1){ pinpis<-pinpis.t2 }else{ pinpis<-rbind(pinpis,pinpis.t2) }
  
}
par(mfrow=c(1,3),mar=c(3,3,1,1))
plot(pinpis$PinPis~pinpis$hsrange,pch=16,col=as.factor(pinpis$Gene_ID))
legend("topleft",legend = levels(as.factor(pinpis$Gene_ID)),col = palette(),pch = 16)
title(main = "PinPis")

plot(pinpis$log10_Pis~pinpis$hsrange,pch=16,col=as.factor(pinpis$Gene_ID))
legend("topleft",legend = levels(as.factor(pinpis$Gene_ID)),col = palette(),pch = 16)
title(main = "Pis")

plot(pinpis$log10_Pin~pinpis$hsrange,pch=16,col=as.factor(pinpis$Gene_ID))
legend("topleft",legend = levels(as.factor(pinpis$Gene_ID)),col = palette(),pch = 16)
title(main = "Pin")


# run bootstrap for the 3 genes

# run GWAs for all 107 phenotypes

phenotypes<-read.table("phenotypes/rawfiles/phenotypes107_gmeans.csv",h=T,sep=",",as.is = 1)
fam<-read.table("Genetics/SNP_1001g_filtered.fam")

hist(apply(X = phenotypes, MARGIN = 2, FUN = function(x){length(na.omit(x))}),xlab="N")

traits<-colnames(phenotypes)[-1]
for (j in 1:length(traits)) {
  fam<-read.table("Genetics/SNP_1001g_filtered.fam")
  for (i in 1:length(fam$V1)) { if(fam$V1[i] %in% phenotypes$accession_id){fam$V6[i]<-phenotypes[,traits[j]][which(phenotypes$accession_id==fam$V1[i])]}}
  fam$V6[which(fam$V6==-9)]<-NA
  write.table(fam,file = paste0("../large_files/Ath_Petal_size/gwas/pheno107/SNP_1001g_filtered_",traits[j],".fam"),sep = "\t",quote = F,row.names = F,col.names = F)
}
write.table(x = traits,file = "../large_files/Ath_Petal_size/gwas/pheno107/list_traits_pheno107.txt",quote = F,row.names = F)

# EHH along HS

# Run snpeff for the 3 genes with raw genetic data

# 