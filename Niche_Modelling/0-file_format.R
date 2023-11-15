####################
#
# Convert .tiff to .asc for use in MaxEnt
# KS - 2023-04-21
#
####################

getwd()
locs<-read.table("../../../Kevin/Science_project/2-Projects/3-Arabidopsis_flower_size/data/newloc_jesse_lasky.txt",h=T)
lomin<-min(locs$NewLongitude)-1
lomax<-max(locs$NewLongitude)+1
lamin<-min(locs$NewLatitude)-1
lamax<-max(locs$NewLatitude)+1

# Make location dataset
locs$species<-rep("A_thaliana",dim(locs)[2])
locs<-locs[,c(3,1,2)]
colnames(locs)<-c("species","longitude","latitude")
write.csv(locs,file = "C:/Users/kesi0002/Downloads/lasky_locs.csv",quote = F,row.names = F)

# Which is to use : proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")
library(raster)
# extract variables used in https://www.biorxiv.org/content/10.1101/2022.03.06.483202v1.full.pdf
variables<-c("bio3","bio6","bio7","bio8","bio10","bio15","bio16","bio17")

# test. to be removed
#lomin<-(-20)
#lomax<-20
#lamin<-30
#lamax<-60
#variables<-c("ai","bio1","bio12")
#locs_test<-locs[-which(locs$longitude<=lomin),]
#locs_test<-locs_test[-which(locs_test$longitude>=lomax),]
#locs_test<-locs_test[-which(locs_test$latitude<=lamin),]
#locs_test<-locs_test[-which(locs$latitude>=lamax),]
#write.csv(locs_test,file = "lasky_locs_test.csv",quote = F,row.names = F)

#### Make loop for each variable
for(i in 1:length(variables)){
  #1-download tiff image
  tiff<- raster(paste0("/vsicurl/https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_",variables[i],"_1981-2010_V.2.1.tif"))
  #2-crop
  extt<- as(extent(lomin,lomax,lamin,lamax), 'SpatialPolygons')
  crs(extt) <- "+proj=longlat +datum=WGS84 +no_defs"
  crpd<-crop(tiff, extt)
  #plot(crpd)
  #2-make raster
  rast<- raster(crpd)
  #4-make asc
  writeRaster(crpd, paste0("layers/for_ath_software/",variables[i],".asc"), format="ascii")
  #5-make .grd for use in R 
  writeRaster(crpd, paste0("layers/for_ath_R/",variables[i],".grd"), format="raster")
}

# fast visualization of a grid
plot(aggr)
points(locs$NewLongitude,locs$NewLatitude,pch=16,col=rgb(.5,.5,.5,.2))

# Make location dataset
locs$species<-rep("A_thaliana",dim(locs)[2])
locs<-locs[,c(3,1,2)]
colnames(locs)<-c("species","longitude","latitude")
write.csv(locs,file = "lasky_locs.csv",quote = F,row.names = F)


# check where are the HS outliers
clim<-raster("layers/for_ath_R/bio10.grd")
phenotypes<-read.table("C:/Users/kesi0002/GitProjects/U_shape/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
adm<-read.table("C:/Users/kesi0002/GitProjects/U_shape/GWAs/results_admixture.txt",h=T,sep = ";",dec = ".")
refs<-which(!duplicated(adm$group))
adm$group<-factor(x = adm$group,levels = adm$group[refs])
adm$couleur<-factor(x = adm$couleur,levels = adm$couleur[refs])
phenotypes<-merge(phenotypes,adm,by.x="Genotype",by.y="AthID",all=T)
pal <- colorRampPalette(c("white","black"))

names(clim)<-"layer"
plot(clim$layer,xlim=c(47,51),ylim=c(37,40),col=pal(10))
points(phenotypes$longitude,phenotypes$latitude,cex=1,lwd=2,pch=21,col="red")

plot(clim$layer,xlim=c(-30,-10),ylim=c(10,30),col=pal(10))
points(phenotypes$longitude,phenotypes$latitude,cex=1,lwd=2,pch=21,col="red")
