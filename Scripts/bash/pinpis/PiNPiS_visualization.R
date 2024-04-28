###################
#
# Plooting the results of PiNPiS analysis
# Kevin Sartori 2022
#
###################

results<-list.files()[grep("result",list.files())] 
pinpis<-read.table(results[1],h = T)[0,] 
for (i in results) {
  pinpisbis<-read.table(i,h=T,na.strings = c("NaN","Inf","-Inf"))
  pinpis<-rbind(pinpis,pinpisbis)
}
rm(pinpisbis)


length(which(pinpis$PinPis>2 ))
pinpis$PinPis[which(pinpis$PinPis>2)]<-2


hist(pinpis$PinPis,breaks = 50,col="white",border = "white",main = "PiN/PiS distribution",ylab="",xlab = "PiN/PiS",las=1,cex.lab=1.3,freq=F,xaxt="n")
title(ylab = "Density", cex.lab=1.3,line = 2.5)
d <- density(na.omit(pinpis$PinPis),n = 1000)
polygon(d$x,d$y,col="grey75",border = "black")
axis(1,at=c(0,0.5,1,1.5,2),labels=c(0,0.5,1,1.5,">=2"))
abline(v=median(pinpis$PinPis,na.rm=T))

pc1<-read.table("/domus/h1/kevinfrs/private/ReversedGWAs_1001g/Pre_analysis/1001g_pc1genelist_20221206.txt")
pc1pinpis<-pinpis[which(pinpis$Gene_ID %in% pc1$V1),] 
d <- density(na.omit(pc1pinpis$PinPis),n = 1000)
polygon(d$x,d$y,col=rgb(0,0.5,0,0.2),border = "black")
abline(v=median(pc1pinpis$PinPis,na.rm=T))


pc2<-read.table("/domus/h1/kevinfrs/private/ReversedGWAs_1001g/Pre_analysis/1001g_pc2genelist_20221206.txt")
pc2pinpis<-pinpis[which(pinpis$Gene_ID %in% pc2$V1),] 
d <- density(na.omit(pc2pinpis$PinPis),n = 1000)
polygon(d$x,d$y,col=rgb(0,0,0.5,0.2),border = "black")
abline(v=median(pc2pinpis$PinPis,na.rm=T))

write.csv2(pinpis,file = "Ath_wholegenome_pinpis.txt",quote = F,row.names = F)

#spot genes

genes<-c("AT3G63300","AT5G43870","AT3G22810","AT4G14740","AT4G32780","AT5G57770","AT4G17350","AT5G47440","AT4G16670")
pinpis[which(pinpis$Gene_ID %in% genes),]
