#------------------
# Petal Size Ath
# Manuscript figures
# Figure 4 - Climate relationship and fitness tradeoff
# 2024-01-02
#------------------

# Panel A - Habitat suitability map

library(raster)
library(viridis)
HS<-raster("../large_files/Ath_Petal_size/habitat_suitability/Niche_modelling/Habitat_suitability_Ath_2023-04-27.grd")
ath_dist<-as(extent(-10, 50, 35, 69), 'SpatialPolygons')
crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
HSath<-crop(HS, ath_dist)
# visualize the map
plot(HS$layer,col=hscol(255),las=1)
hscol<-colorRampPalette(rev(viridis::inferno(4)))
plot(HSath$layer,col=hscol(255),las=1)
# studied accessions
acc<-na.omit(read.table("Genetics/studied_acc.txt"))
g1001<-na.omit(read.table("Genetics/accessions.txt",h=F,sep=",",as.is = 4)[,c(1,6,7,11)])
colnames(g1001)<-c("accession_name","latitude","longitude","group")
g1001<-g1001[which(g1001$accession_name %in% acc$V1),]
# plot the collecting sites on the map
points(g1001$longitude,g1001$latitude,pch=21,col=rgb(0,.6,.6),cex=1,lwd=3)

# Panel B - Trait-HS relationship
library(quantreg)
# data from this study
phenotypes<-read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
phenotypes<-phenotypes[-which(is.na(phenotypes$Petal_Area)),]
habsuit<-read.table("Niche_Modelling/accessions_1001g_habitatsuitability.csv",h=T,sep=",")
# few accessions were removed due to poor climatic resolution in site surounded by sea, leading to a potential poor estimation of HS
habsuit<-habsuit[-which(habsuit$HS < 0.2),]
phenotypes<-merge(phenotypes,habsuit,by.x="Genotype",by.y="accession_name",all.x=T)
par(mar=c(4,4,0,0))
# Petal
layout(matrix(c(1,1,1,1,2,3,4,5),byrow = T,nrow=2),heights = c(2,1))
plot(phenotypes$Petal_Area~phenotypes$HS,pch=21,lwd=1,bg=rgb(0,.6,.6),cex=1.5,
     xlab="Habitat suitability",ylab="Petal Area")
#LM
summod<-summary(lm(phenotypes$Petal_Area~phenotypes$HS))
summod$coefficients
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*summod$coefficients[2,1]+summod$coefficients[1,1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*summod$coefficients[2,1]+summod$coefficients[1,1],
         lwd=2,lty=1)
#QR
rqfit95 <- rq(phenotypes$Petal_Area~phenotypes$HS,tau=0.95)
rqfit05 <- rq(phenotypes$Petal_Area~phenotypes$HS,tau=0.05)
anova(rqfit95,rqfit05)
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*rqfit95$coefficients[2]+rqfit95$coefficients[1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*rqfit95$coefficients[2]+rqfit95$coefficients[1],
         lwd=1,lty=2)
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*rqfit05$coefficients[2]+rqfit05$coefficients[1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*rqfit05$coefficients[2]+rqfit05$coefficients[1],
         lwd=1,lty=2)

# Camemberts
# function:
datapie<-function(traits,dataset){
for (i in 1:(length(traits)-1)) {
quantiprob<-length(na.omit(unique(dataset[,i])))/length(na.omit(dataset[,i]))
if(quantiprob<.1){
  if(i==1){
  res<-data.frame(trait=traits[i],slope="untested",triangle="untested",type="untested")
  }else{res<-rbind(res,c(traits[i],"untested","untested","untested"))}
}else{
# lm
summod<-summary(lm(dataset[,i]~dataset$HS))
if(summod$coefficients[2,4]<0.01){
slope<-c("negative","neutral","positive")[sign(summod$coefficients[2,1])+2]
}else{slope<-"flat"}
# qr
rqfit95 <- rq(dataset[,i]~dataset$HS,tau=0.95)
rqfit05 <- rq(dataset[,i]~dataset$HS,tau=0.05)
if(0 %in% c(rqfit05$coefficients[2],rqfit95$coefficients[2])){
  triangle<-"untested"
  type<-slope
}else{
tritest<-anova(rqfit95,rqfit05)
if(tritest$table[1,4]<0.01){
  triangle<-"triangular"
  type<-paste0(slope," ",triangle)
}else{
  triangle<-NA
  type<-slope}
}
# make dataset
if(i==1){
  res<-data.frame(trait=traits[i],slope=slope,triangle=triangle,type=type)
}else{res<-rbind(res,c(traits[i],slope,triangle,type))}
}
}
res$type<-factor(x = res$type,levels = c("negative triangular","negative","flat triangular","flat","positive triangular","positive","untested"))
return(res)
}

# plots
par(mar=c(0,0,0,0))
# All data from this study
traits<-colnames(phenotypes)[c(5:16,22,27)]
dataset<-phenotypes[,traits]
(res<-datapie(traits,dataset))
palette(c(rgb(0,0,1),rgb(.6,.6,1),rgb(0,0,0),rgb(.6,.6,.6),rgb(1,0,0),rgb(1,.6,.6),rgb(1,1,1)))
pie(table(res$type),col = palette(),labels = sub("0"," ",table(res$type)))
axis(side = 1,at = 0,labels = "This study",tick = F,line = -2)
# 107 phenotypes
pheno107<-read.table("Phenotypes/rawfiles/phenotypes107_gmeans.csv",h=T,sep=",")
pheno107<-merge(pheno107,habsuit[,c(1,5)],by.x="accession_id",by.y="accession_name",all.x=T)
traits<-colnames(pheno107)[-1]
dataset<-pheno107[,traits]
(res<-datapie(traits = traits,dataset = dataset))
palette(c(rgb(0,0,1),rgb(.6,.6,1),rgb(0,0,0),rgb(.6,.6,.6),rgb(1,0,0),rgb(1,.6,.6),rgb(1,1,1)))
pie(table(res$type),col = palette(),labels = sub("0"," ",table(res$type)))
axis(side = 1,at = 0,labels = "Atwell et al. 2010",tick = F,line = -2)
# Przybylska data
Przybylska<-read.table("Phenotypes/rawfiles/phenotypic_datarecord.txt",h=T,sep="\t",dec=",")
Przybylska$X1001g_ID<-as.factor(Przybylska$X1001g_ID)
traits<-levels(as.factor(Przybylska$traitName))
accessions<-levels(as.factor(Przybylska$X1001g_ID))
dataset<-data.frame(matrix(data = NA,nrow = length(accessions),ncol = length(traits)+1,dimnames = list(accessions,c(traits,"accession_name"))))
dataset$accession_name<-accessions
for (i in 1:length(traits)) {
dataset.t<-na.omit(Przybylska[which(Przybylska$traitName==traits[i]),c(2,8,9)])
dataset[,i]<-as.numeric(tapply(X = dataset.t$traitValue, INDEX = dataset.t$X1001g_ID,FUN = mean))
}
dataset<-merge(dataset,habsuit[,c(1,5)],by="accession_name",all.x=T)
traits<-colnames(dataset)[-1]
dataset<-dataset[,traits]
(res<-datapie(traits = traits,dataset = dataset))
palette(c(rgb(0,0,1),rgb(.6,.6,1),rgb(0,0,0),rgb(.6,.6,.6),rgb(1,0,0),rgb(1,.6,.6),rgb(1,1,1)))
pie(table(res$type),col = palette(),labels = sub("0"," ",table(res$type)))
axis(side = 1,at = 0,labels = "Przybylska et al. 2023",tick = F,line = -2)
# Legend
par(mar=c(0,0,0,0))
plot(rep(1,7),1:7,xlim=c(0,10),ylim=c(0,9),pch=22,bg=rev(palette()),cex=4,bty="n",xaxt="n",yaxt="n")
text(x = c(rep(2,7),1),y = c(1:7,8.5), pos=4, labels=c(rev(levels(res$type)),"Legend"),font=c(rep(1,7),2))

# Panel C - PiN/PiS-HS relationship
# !!!!! need to rerun pinpis-HS avec les bonnes donnees

# Panel D - SFS-HS
handtable<-as.matrix(read.table("Genetics/sfs/DAF_table_rel.95_.1to.5.txt"))
colnames(handtable)<-substr(colnames(handtable),2,5)
shade<-gray.colors(10,start = .8,end = .2,gamma = 1)
barplot(height = handtable[,1:9],beside = T,col=shade,las=1,main = "Allele counts per frequency")
legend(x = 75,y = 4e+05,legend = rownames(handtable),fill = shade,cex = .75,ncol = 2)

# Panel E - Large petal allele freq-HS
snpeffect <- read.table("Genetics/bslmm_top100_Leaf_Area.param.txt",h=T,sep="\t",dec=".")
snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma



