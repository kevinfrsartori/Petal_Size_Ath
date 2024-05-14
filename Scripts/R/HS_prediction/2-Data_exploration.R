####################
#
# Hbitat suitability data exploration
# KS - 2023-04-27
#
####################

# import raster image with HS score for A. thaliana distribution range 
HS<-stack("Habitat_suitability_Ath_2023-04-27.grd")

# visualize the map
plot(HS$layer)

# import dataset with A. thaliana collecting sites
g1001<-na.omit(read.table("accessions.txt",h=F,sep=",",as.is = 4)[,c(1,6,7,11)])
colnames(g1001)<-c("accession_name","latitude","longitude","group")

# plot the collecting sites on the map
#points(g1001$longitude,g1001$latitude,pch=".",cex=2)

# make the locations as a spatial object coordinates
coor<- SpatialPoints(g1001[,c("longitude","latitude")],proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# extract HS values for the coordinates and convert to dataset
#values<-data.frame(extract(HS$layer,coor,method="bilinear"))
values<-data.frame(extract(HS$layer,coor,buffer=1000,fun=max))
names(values)<-"HS"

# merge with accession dataset
g1001<-cbind(g1001,values)
write.csv(x = g1001,file = "accessions_1001g_habitatsuitability.csv",quote = F,row.names = F)
# make the boxplot
par(mar=c(6,4,1,1),oma=c(6,1,1,1))
plot(g1001$HS~g1001$group,las=2,xlab="",ylab="Habitat Suitability")

# Zoom
plot(HS$layer,xlim=c(-10,70),ylim=c(30,70))
points(g1001$longitude,g1001$latitude,pch=".",cex=2)



phenotype<-read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",sep=",",h=T)
g1001<-merge(g1001,phenotype,by.x="accession_name",by.y="Genotype",all=T)
par(oma=c(2,2,1,1))
plot(g1001$Petal_Area~g1001$HS,xlim=c(0,0.8),pch=16,col=g1001$group,cex=1,ylab="Petal Area (mm2)",xlab="Habitat Suitability")
plot(g1001$Seed~g1001$HS,xlim=c(0,0.8),pch=16,col=g1001$group,cex=1,ylab="Ovule Number",xlab="Habitat Suitability")
plot(g1001$Sepal_Area~g1001$HS,xlim=c(0,0.8),pch=16,col=g1001$group,cex=1,ylab="Sepal Area (mm2)",xlab="Habitat Suitability")
plot(g1001$Leaf_Area~g1001$HS,xlim=c(0,0.8),pch=16,col=g1001$group,cex=1,ylab="Leaf Area (mm2)",xlab="Habitat Suitability")
plot(g1001$flowering_time~g1001$HS,xlim=c(0,0.8),pch=16,col=g1001$group,cex=1,ylab="Leaf Area (mm2)",xlab="Habitat Suitability")



hist(g1001$HS)
quantile(g1001$HS,c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
#0%        10%        20%        30%        40%        50%        60%        70%        80%        90%       100% 
#0.01421496 0.45431006 0.52085227 0.56103659 0.58834815 0.60821909 0.61850846 0.64275384 0.66511321 0.68610811 0.77811801 


g01<-g1001$accession_name[which(g1001$HS < quantile(g1001$HS,.1))]
g02<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.1) & g1001$HS < quantile(g1001$HS,.2))]
g03<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.2) & g1001$HS < quantile(g1001$HS,.3))]
g04<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.3) & g1001$HS < quantile(g1001$HS,.4))]
g05<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.4) & g1001$HS < quantile(g1001$HS,.5))]
g06<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.5) & g1001$HS < quantile(g1001$HS,.6))]
g07<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.6) & g1001$HS < quantile(g1001$HS,.7))]
g08<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.7) & g1001$HS < quantile(g1001$HS,.8))]
g09<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.8) & g1001$HS < quantile(g1001$HS,.9))]
g10<-g1001$accession_name[which(g1001$HS > quantile(g1001$HS,.9) & g1001$HS < quantile(g1001$HS,1))]

write.table(x = data.frame(V1=g01,v2=g01),file = "accessions_HS01.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g02,v2=g02),file = "accessions_HS02.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g03,v2=g03),file = "accessions_HS03.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g04,v2=g04),file = "accessions_HS04.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g05,v2=g05),file = "accessions_HS05.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g06,v2=g06),file = "accessions_HS06.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g07,v2=g07),file = "accessions_HS07.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g08,v2=g08),file = "accessions_HS08.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g09,v2=g09),file = "accessions_HS09.txt",quote = F,col.names = F,row.names = F)
write.table(x = data.frame(V1=g10,v2=g10),file = "accessions_HS10.txt",quote = F,col.names = F,row.names = F)


# I am making lists of accessions 
# I want to know if accession from low HS have less genetic diversity than accessions from high HS


#### Limiting factor
# import raster image with HS score for A. thaliana distribution range 
lim<-stack("Limiting_facgtors_Ath_2023-07-07.grd")
# visualize the map
plot(lim)
# import dataset with A. thaliana collecting sites
g1001<-na.omit(read.table("accessions_1001g_habitatsuitability.csv",h=T,sep=","))

# make the locations as a spatial object coordinates
coor<- SpatialPoints(g1001[,c("longitude","latitude")],proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# extract limiting factor for the coordinates and convert to dataset
values<-data.frame(extract(lim$layer,coor,method="simple"))
names(values)<-"limiting_factor"

# convert number by bioclim names
ath_rast<-stack("Bioclim_Ath_1001g.grd")
shortname<-names(ath_rast)
longname<-c("isothermality","Temp_coldest_month","Annual_temp_range","Temp_wettest_quarter","Temp_warmest_quarter","Prec_seasonality","Prec_wettest_quarter","Prec_driest_quarter","Altitude")
shortlong<-data.frame(shortname=shortname,longname=longname)

g1001$LF<-NA
for (i in 1:dim(g1001)[1]) {
g1001$LF[i]<-shortlong$longname[values[i,]]  
}

write.csv(x = g1001,file = "accessions_1001g_habitatsuitability.csv",quote = F,row.names = F)

par(oma=c(5,5,1,1))
boxplot(g1001$HS~g1001$LF,las=2)
