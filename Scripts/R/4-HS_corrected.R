  #------------------
# Petal Size Ath
# Manuscript figures
# Figure 3 - Climate relationship and fitness tradeoff
# 2024-01-02
#------------------

# Panel a - Habitat suitability map 650 x 500
#----------------------------------
par(mar=c(3,3,2,0),oma=c(0,0,0,0))
library(raster)
library(rasterVis)
library(viridis)
HS<-raster("../large_files/Ath_Petal_size/habitat_suitability/Niche_modelling/Habitat_suitability_Ath_2023-04-27.grd")
# crop
ath_dist<-as(extent(-10, 50, 35, 69), 'SpatialPolygons')
crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
HSath<-crop(HS, ath_dist)
# visualize the map
hscol<-colorRampPalette(rev(viridis::inferno(4)))
#plot(HS$layer,col=hscol(255),las=1) # whole
plot(HSath$layer,col=hscol(255),las=1, legend=F) # cropped
# studied accessions
acc<-na.omit(read.table("Genetics/studied_acc.txt"))
g1001<-na.omit(read.table("Genetics/accessions.txt",h=F,sep=",",as.is = 4)[,c(1,6,7,11)])
colnames(g1001)<-c("accession_name","latitude","longitude","group")
g1001<-g1001[which(g1001$accession_name %in% acc$V1),]
# plot the studied sites on the map
points(g1001$longitude,g1001$latitude,pch=21,col=rgb(0,.6,.6),cex=.8,lwd=2)

# Crop seas out
nosea<-rnaturalearthdata::map_units50
HSathEUnosea<-terra::mask(HSath,mask = nosea)
hscol<-colorRampPalette( rev( c( viridis::inferno(4)[c(1,2,3,4,4)], "#FFFFFF" ) ) )
plot(HSathEU,col=hscol(255),las=1, legend=T) # cropped,

#points(g1001$longitude,g1001$latitude,pch=21,bg="white",cex=1)
#points(g1001$longitude,g1001$latitude,pch=21,bg=rgb(0,0,0,0),col=rgb(0,.6,.6),cex=1,lwd=2)
points(g1001$longitude,g1001$latitude,pch=21,bg=rgb(0,.6,.6),cex=1)

# For Supp, full map 1300 300
# crop ylim=c(15,70),xlim=c(-150,150)
ath_dist<-as(extent(-150, 150, 15, 70), 'SpatialPolygons')
crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
HSath<-crop(HS, ath_dist)
# visualize the map
hscol<-colorRampPalette( rev( c( viridis::inferno(4)[c(1,2,3,4)], "#FFFFFF" ) ) )
#plot(HS$layer,col=hscol(255),las=1) # whole
plot(HSath$layer,col=hscol(255),las=1, legend=F) # cropped
#crop seas out
HSathnosea<-terra::mask(HSath,mask = nosea)
plot(HSathnosea,col=hscol(255),las=1)

#For supp: Limiting factors
lim<-raster("../large_files/Ath_Petal_size/habitat_suitability/Niche_modelling/Limiting_facgtors_Ath_2023-07-07.grd")
# visualize the map
longname<-c("isothermality","Temp_coldest_month","Annual_temp_range","Temp_wettest_quarter","Temp_warmest_quarter","Prec_seasonality","Prec_wettest_quarter","Prec_driest_quarter","Altitude")
mycol=c("purple","coral","wheat4","coral3","coral4","slategray2","skyblue2","skyblue4","black")
ath_dist<-as(extent(-10, 50, 35, 69), 'SpatialPolygons')
crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
limath<-crop(lim, ath_dist)
limath<-ratify(limath)
levelplot(limath,col.regions=mycol,att="ID",colorkey=F)

lim<-ratify(lim)
levelplot(lim,col.regions=mycol,att="ID",colorkey=F)

lim<-ratify(limath)
levelplot(limath,col.regions=mycol,att="ID",colorkey=F)

lim<-terra::mask(lim,mask = europe)
levelplot(lim,col.regions=mycol,att="ID",colorkey=F)


# Legend
plot(rep(0,9),1:9,xlim=c(0,10),ylim=c(9,-1),cex=3,pch=22,col="black",bty="n",xaxt="n",yaxt="n",ylab="",xlab="",
     bg=c("purple","coral","wheat4","coral3","coral4","slategray2","skyblue2","skyblue4","black"))
text(1,-.5,"Legend",font=2)
text(rep(.5,9),1:9,pos=4,labels=longname)

# Petal area as a function of limiting factor
phenotypes<-read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
limiting<-read.table("Niche_Modelling/accessions_1001g_habitatsuitability.csv",h=T,sep=",")
phenotypes<-merge(phenotypes,limiting,by.x="Genotype",by.y="accession_name",all.x=T)
levels(as.factor(phenotypes$LF))


longname<-rev(c("isothermality","Temp_coldest_month","A,nnual_temp_range","Temp_wettest_quarter","Temp_warmest_quarter","Prec_seasonality","Prec_wettest_quarter","Prec_driest_quarter","Altitude"))
mycol=rev(c("purple","coral","wheat4","coral3","coral4","slategray2","skyblue2","skyblue4","black"))
phenotypes$LF<-factor(x = phenotypes$LF,levels = longname)
par(mar=c(4,10,1,1))
boxplot(phenotypes$Petal_Area ~ phenotypes$LF,horizontal=T,las=1,col=mycol,ylab="",ylim=c(0,4),xlim=c(.5,10),xlab="Petal area (mm2)")
axis(side = 2,at = 10,labels = "Limiting factor : ",las=1,tick = F)
#anova
anova<-aov(Petal_Area ~ LF, data = phenotypes)
(tukey<-TukeyHSD(anova))
library(multcomp)
glht_mod<-glht(anova, linfct=mcp(LF="Tukey"))
letters<-cld(glht_mod)$mcletters$Letters
text(0,c(1,2,3,4,7,8,9),letters,pos=4)

# Panel b - Trait-HS relationship 500 x 500
#--------------------------------

library(quantreg)
# 1 - data from this study
phenotypes<-read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
phenotypes<-phenotypes[-which(is.na(phenotypes$Petal_Area)),]
# habitat suitability data
habsuit<-read.table("Niche_Modelling/accessions_1001g_habitatsuitability.csv",h=T,sep=",")
# few accessions were removed due to poor climatic resolution in site surrounded by sea, leading to a potential poor estimation of HS
habsuit<-habsuit[-which(habsuit$HS < 0.2),]
phenotypes<-merge(phenotypes,habsuit,by.x="Genotype",by.y="accession_name",all.x=T)
# Biplot Petal Area ~ habitat suitability
par(mar=c(4,5.5,.1,.1))
layout(matrix(c(1,1,1,1,2,3,4,5),byrow = T,nrow=2),heights = c(2,1))

plot(phenotypes$Petal_Area~phenotypes$HS,pch=21,lwd=1,bg=rgb(0,.6,.6),cex.axis=2,las=1,cex=1.5,
     xlab="Habitat suitability",ylab=expression(Petal ~ Area ~ (mm^2)),cex.lab=2)
# linear model
summod<-summary(lm(phenotypes$Petal_Area~phenotypes$HS))
summod$coefficients
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*summod$coefficients[2,1]+summod$coefficients[1,1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*summod$coefficients[2,1]+summod$coefficients[1,1],
         lwd=2,lty=1)
# quantile regressions
rqfit95 <- rq(phenotypes$Petal_Area~phenotypes$HS,tau=0.95)
rqfit05 <- rq(phenotypes$Petal_Area~phenotypes$HS,tau=0.05)
anova(rqfit95,rqfit05)
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*rqfit95$coefficients[2]+rqfit95$coefficients[1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*rqfit95$coefficients[2]+rqfit95$coefficients[1],
         lwd=1,lty=2)
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*rqfit05$coefficients[2]+rqfit05$coefficients[1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*rqfit05$coefficients[2]+rqfit05$coefficients[1],
         lwd=1,lty=2)

# 2 - Pie chart for other datasets
# make function for reproducibility
datapie<-function(traits,dataset){
for (i in 1:(length(traits)-1)) {
# estimate probability that a variable is qualitative (if there are too few categories compare to the number of data, we do not process it)
quantiprob<-length(na.omit(unique(dataset[,i])))/length(na.omit(dataset[,i]))
if(quantiprob<.1){
  if(i==1){
  res<-data.frame(trait=traits[i],slope="untested",triangle="untested",type="untested")
  }else{res<-rbind(res,c(traits[i],"untested","untested","untested"))}
}else{
# linear model
summod<-summary(lm(dataset[,i]~dataset$HS))
if(summod$coefficients[2,4]<0.01){
slope<-c("negative","neutral","positive")[sign(summod$coefficients[2,1])+2]
}else{slope<-"flat"}
# quantile regression
rqfit95 <- rq(dataset[,i]~dataset$HS,tau=0.95)
rqfit05 <- rq(dataset[,i]~dataset$HS,tau=0.05)
# if one slope is 0, the anova crash. We do not run the anova in that case.
if(0 %in% c(rqfit05$coefficients[2],rqfit95$coefficients[2])){
  triangle<-"untested"
  type<-slope
}else{
# anova
tritest<-anova(rqfit95,rqfit05)
# If slopes are significantly different, the relationship is triangular
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
# define res$type as factor with predefined levels
res$type<-factor(x = res$type,levels = c("negative triangular","negative","flat triangular","flat","positive triangular","positive","untested"))
return(res)
}

# panel c - Plot the pies
par(mar=c(0,0,0,0))
# 2.1 - All data from this study
traits<-colnames(phenotypes)[c(5:16,22,27)]
dataset<-phenotypes[,traits]
(res<-datapie(traits,dataset))
palette(c(rgb(0,0,1),rgb(.6,.6,1),rgb(0,0,0),rgb(.6,.6,.6),rgb(1,0,0),rgb(1,.6,.6),rgb(1,1,1)))
pie(table(res$type),col = palette(),labels = sub("0"," ",table(res$type)))
axis(side = 1,at = 0,labels = "This study",tick = F,line = -2,cex.axis=1.5)
# 2.2 - 107 phenotypes from Atwell et al. 2010
pheno107<-read.table("Phenotypes/rawfiles/phenotypes107_gmeans.csv",h=T,sep=",")
pheno107<-merge(pheno107,habsuit[,c(1,5)],by.x="accession_id",by.y="accession_name",all.x=T)
traits<-colnames(pheno107)[-1]
dataset<-pheno107[,traits]
(res<-datapie(traits = traits,dataset = dataset))
palette(c(rgb(0,0,1),rgb(.6,.6,1),rgb(0,0,0),rgb(.6,.6,.6),rgb(1,0,0),rgb(1,.6,.6),rgb(1,1,1)))
pie(table(res$type),col = palette(),labels = sub("0"," ",table(res$type)))
axis(side = 1,at = 0,labels = "Atwell et al. 2010",tick = F,line = -2,cex.axis=1.5)
# 2.3 - Przybylska et al. 2023
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
axis(side = 1,at = 0,labels = "Przybylska et al. 2023",tick = F,line = -2,cex.axis=1.5)
# Legend
par(mar=c(0,0,0,0))
plot(rep(1,7),1:7,xlim=c(0,10),ylim=c(0,9),pch=22,bg=rev(palette()),cex=4,bty="n",xaxt="n",yaxt="n")
type<-c("untested","positive","pos. triangular","flat","flat triangular","negative","neg. triangular")
text(x = c(rep(2,7),1),y = c(1:7,8.5), pos=4, labels=c(type,"Legend"),font=c(rep(1,7),2),cex=1.5)


dev.off()

# Panel d - Large petal allele freq-HS 350 x 200
#-------------------------------------
# allele 0 versus 1 recovered from glm data
assoc<-read.table("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_Petal_Area.assoc.txt",h=T,sep="\t",dec=".")
# effect recovered from bslmm
snpeffect <- read.table("Genetics/bslmm_flct_Petal_Area.param.txt",h=T,sep="\t",dec=".")
snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma
snpeffect<-merge(snpeffect,assoc[,c(2,5,6)],by="rs",all.x=T)
# define "large petal allele"
snpeffect$large_petal_allele<-snpeffect$allele1
snpeffect$large_petal_allele[which(snpeffect$snpeffect<0)]<-snpeffect$allele0[which(snpeffect$snpeffect<0)]
rm(assoc)
# frq data per HS
hsranges<-c(paste0("HS0",1:9),"HS10")
for (i in hsranges) {
hs<-read.table(paste0("Genetics/frq/Petal_Area.",i,".frq"))
hs<-hs[which(hs$V7 %in% snpeffect$rs),]
hs$allele_a<-unlist(lapply(X = strsplit(x = hs$V5, split = ":"),"[[",1))
hs$frq_allele_a<-as.numeric(unlist(lapply(X = strsplit(x = hs$V5, split = ":"),"[[",2)))
hs$allele_b<-unlist(lapply(X = strsplit(x = hs$V6, split = ":"),"[[",1))
hs$frq_allele_b<-as.numeric(unlist(lapply(X = strsplit(x = hs$V6, split = ":"),"[[",2)))
# merge large petal allele
hs<-merge(hs,snpeffect[,c("rs","large_petal_allele")],by.x="V7",by.y="rs",all.x=T)
# spot large petal allele freq
hs$large_petal_allele_frq<-NA
hs$large_petal_allele_frq[which(hs$allele_a==hs$large_petal_allele)]<-hs$frq_allele_a[which(hs$allele_a==hs$large_petal_allele)]
hs$large_petal_allele_frq[which(hs$allele_b==hs$large_petal_allele)]<-hs$frq_allele_b[which(hs$allele_b==hs$large_petal_allele)]
# merge with snpeffect
snpeffect<-merge(snpeffect,hs[,c("V7","large_petal_allele_frq")],by.x="rs",by.y="V7",all.x = T)
# rename
colnames(snpeffect)[dim(snpeffect)[2]]<-paste0("large_petal_allele_frq_",i)
}
# Barplot all large petal alleles 500 x 250 
par(mar=c(1,1,1,1), oma=c(2,3,0,0))
large_petal_allele_frq_per_hsrange<-apply(snpeffect[,grep("HS",colnames(snpeffect))],MARGIN = 2,FUN = sum,na.rm=T)/length(na.omit(snpeffect$large_petal_allele_frq_HS01))
hscol<-colorRampPalette(rev(viridis::inferno(4)))
barplot(large_petal_allele_frq_per_hsrange,ylim=c(0,.25),xpd=F,las=1,col=hscol(20)[6:15],space=0,xaxt="n")
axis(side = 1,at = 1:10-.5,labels = rep("",10),las=2)
axis(side = 1,at = 5,labels = "less          -          Habitat suitability          -          more",tick = F,line = 0)
axis(side = 2,at = .125,labels = "Large Petal Allele frequency",tick = F,line = 2)


# Panel e - SFS-HS 450 x 600
#-----------------
dev.off()
par(mfrow=c(1,1))

# PETAL
trait="Petal"
handtable<-as.matrix(read.table(paste0("Genetics/sfs/unfolded_SFS_",trait,".txt")))
colnames(handtable)<-substr(colnames(handtable),2,5)
#shade<-gray.colors(10,start = .8,end = .2,gamma = 1)
hscol<-colorRampPalette(rev(viridis::inferno(4)))
par(mar=c(4,4,4,0),oma=c(.1,.1,.1,.1))
barplot(height = handtable[,1:4],beside = T,col=hscol(20)[6:15],las=1,
        main = paste0(trait," genes (", sum(handtable[1,]), " SNPs)" ),ylim=c(0,1000),
        names.arg = paste0("<",colnames(handtable[,1:4])),
        xlab = "Allele frequency",ylab="Allele count")

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v <- c( (v[1]+v[2])/2.3, v[2], (v[3]+v[4])/2.3, v[4] )
par( fig=v, new=TRUE, mar=c(3,3,0,0) )
library(e1071)
handtable<-handtable[,1:9]
plot(x = 1:10, y = apply(X = handtable, MARGIN = 1, FUN = skewness),ylim=c(1,2),
     las=1,ylab="skewness",xlab="HS",pch=21,bg=hscol(20)[6:15],cex=2,yaxt="n",xaxt="n")
axis(side = 2,at = c(1,2), labels = c(1,2), las=1)
axis(side = 2,at = 1.5, labels = "SFS Skewness", tick = F)
axis(side = 1,at = c(1,10), labels = c(1,10),tick = F)
axis(side = 1,at = 1:10, labels = c("","","","","","","","","",""))
axis(side = 1,at = 5.5, labels = "HS", tick = F)


y<-apply(X = handtable, MARGIN = 1, FUN = skewness)
x<-1:10
reg<-summary(lm(y~x))

segments(x0 = 1,y0 = reg$coefficients[1,1]+reg$coefficients[2,1]*1,x1 = 10,y1 = reg$coefficients[1,1]+reg$coefficients[2,1]*10, lwd = 2)

#LEAF
dev.off()
trait="Leaf"
handtable<-as.matrix(read.table(paste0("Genetics/sfs/unfolded_SFS_",trait,".txt")))
colnames(handtable)<-substr(colnames(handtable),2,5)
#shade<-gray.colors(10,start = .8,end = .2,gamma = 1)
hscol<-colorRampPalette(rev(viridis::inferno(4)))
par(mar=c(4,4,4,0),oma=c(.1,.1,.1,.1))
barplot(height = handtable[,1:4],beside = T,col=hscol(20)[6:15],las=1,
        main = paste0(trait," genes (", sum(handtable[1,]), " SNPs)" ),ylim=c(0,1000),
        names.arg = paste0("<",colnames(handtable[,1:4])),
        xlab = "Allele frequency",ylab="Allele count")

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v <- c( (v[1]+v[2])/2.3, v[2], (v[3]+v[4])/2.3, v[4] )
par( fig=v, new=TRUE, mar=c(3,3,0,0) )
library(e1071)
handtable<-handtable[,1:9]
plot(x = 1:10, y = apply(X = handtable, MARGIN = 1, FUN = skewness),ylim=c(1,2),
     las=1,ylab="skewness",xlab="HS",pch=21,bg=hscol(20)[6:15],cex=2,yaxt="n",xaxt="n")
axis(side = 2,at = c(1,2), labels = c(1,2), las=1)
axis(side = 2,at = 1.5, labels = "SFS Skewness", tick = F)
axis(side = 1,at = c(1,10), labels = c(1,10),tick = F)
axis(side = 1,at = 1:10, labels = c("","","","","","","","","",""))
axis(side = 1,at = 5.5, labels = "HS", tick = F)
y<-apply(X = handtable, MARGIN = 1, FUN = skewness)
x<-1:10
reg<-summary(lm(y~x))

segments(x0 = 1,y0 = mean(y),x1 = 10,y1 = mean(y), lwd = 2, lty = 2)
