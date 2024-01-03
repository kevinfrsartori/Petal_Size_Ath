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
habsuit<-read.table("Niche_Modelling/accessions_1001g_habitatsuitability.csv",h=T,sep=",")
# few accessions were removed due to poor climatic resolution in site surounded by sea, leading to a potential poor estimation of HS
habsuit<-habsuit[-which(habsuit$HS < 0.2),]
phenotypes<-merge(phenotypes,habsuit,by.x="Genotype",by.y="accession_name",all.x=T)
par(mar=c(4,4,0,0))
# Sepal
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
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*rqfit95$coefficients[2]+rqfit95$coefficients[1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*rqfit95$coefficients[2]+rqfit95$coefficients[1],
         lwd=1,lty=2)
rqfit05 <- rq(phenotypes$Petal_Area~phenotypes$HS,tau=0.05)
segments(min(phenotypes$HS,na.rm = T),min(phenotypes$HS,na.rm = T)*rqfit05$coefficients[2]+rqfit05$coefficients[1],
         max(phenotypes$HS,na.rm = T),max(phenotypes$HS,na.rm = T)*rqfit05$coefficients[2]+rqfit05$coefficients[1],
         lwd=1,lty=2)
anova(rqfit95,rqfit05)


# data from other studies



# Panel C - PiN/PiS-HS relationship

# Panel D - SFS-HS
handtable<-as.matrix(read.table("Genetics/sfs/DAF_table_rel.95_.1to.5.txt"))
colnames(handtable)<-substr(colnames(handtable),2,5)
shade<-gray.colors(10,start = .8,end = .2,gamma = 1)
barplot(height = handtable[,1:9],beside = T,col=shade,las=1,main = "Allele counts per frequency")
legend(x = 75,y = 4e+05,legend = rownames(handtable),fill = shade,cex = .75,ncol = 2)

# Panel E - Large petal allele freq-HS
