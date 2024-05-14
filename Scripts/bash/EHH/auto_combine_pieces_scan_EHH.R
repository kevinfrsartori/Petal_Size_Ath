#!/usr/bin/env Rscript
argum=commandArgs(trailingOnly=T)
#argum=c("4")

##########
#
# Extended haplotype homozygosity
# plot-scan 
# 2023-07-21
#
##########

print(paste0("Combining pieces of ",argum[1] ))

setwd("$PATH-TO-SERVER/A_thaliana/fastphase/")
let<-read.table(paste0("letters.",argum[1] ,".txt"))
for (i in 1:dim(let)[1] ) {
  if(file.exists(paste0("$PATH-TO-SERVER/A_thaliana/fastphase/scanEHH.chrm-",argum[1] ,".",let[i,1],".txt"))){
  scanhh_i<-read.table(paste0("$PATH-TO-SERVER/A_thaliana/fastphase/scanEHH.chrm-",argum[1] ,".",let[i,1],".txt"),h=T)
  if (!"scanhh" %in% ls()) { scanhh<-scanhh_i }else{ scanhh<-rbind(scanhh,scanhh_i) 
  }
  }else{
    print(paste0("! missing piece: ",as.character(let[i,1])))}   
}
rm(scanhh_i)  
#plot(scanhh$IES~scanhh$POSITION, type="l",ylim=c(0,10000))

ip<-read.table(paste0("$PATH-TO-SERVER/private/Athaliana_flowersize/rehh/interpiece_snps.",argum[1],".txt"))
ip<-ip[which(ip$V2==substr(argum[1],1,1)),]
for (i in 1:dim(ip)[1] ) {
  inter_i<-read.table(paste0("$PATH-TO-SERVER/A_thaliana/fastphase/scanEHH.chrm-",ip$V2[i],".",ip$V3[i],".txt"),h=T)
  colnames(inter_i)[6]<-"IES_ip"
  colnames(inter_i)[4]<-"IHHa_ip"
  colnames(inter_i)[5]<-"IHHd_ip"
  scanhh<-merge(scanhh,inter_i[,c(2,6,4,5)],by="POSITION",all.x = T)
  scanhh$IHHa<-apply(scanhh[,c(4,8)],1,max,na.rm=T)
  scanhh$IHHd<-apply(scanhh[,c(5,9)],1,max,na.rm=T)
  scanhh$IES<-apply(scanhh[,c(6,7)],1,max,na.rm=T)
  scanhh<-scanhh[,-c(7,8,9)] 
}

#lines(scanhh$IES ~  scanhh$POSITION,type="l",lwd=2,col="red")

write.table(x = scanhh,file = paste0("$PATH-TO-SERVER/A_thaliana/fastphase/scanEHH.chrm-",argum[1],".txt"),quote = F,row.names = F,sep=" ")

print("done")