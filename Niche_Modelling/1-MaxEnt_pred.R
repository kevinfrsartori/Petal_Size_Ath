####################
#
# Estimate Habitat suitability
# with MaxEnt
# KS - 2023-04-25
#
####################

###### 1- Building the model

#library(devtools)
#install_github('johnbaums/rmaxent')
library(rmaxent)
# import localisation file (the ones for modelling)
occ <- read.table("samples/lasky_locs.csv",header = T,sep = ",")[-1]
colnames(occ)<-c("lon","lat")
#import rasters with climatic variables (.grd files made in script file_format.R)
library(raster)
pred_files<-list.files("layers/for_ath_R/", '\\.grd$',full.names = T)
#stack the layers 
predictors <- stack(pred_files)
#check names
names(predictors)
# Run MaxEnt
library(dismo)
library(rJava)
me <- maxent(predictors, occ)
#save the model
save(me,file="MaxEnt_Model_Ath.Rdata")
rm(list=ls())


###### 2- Building the raster file containing layers with climatic variables for A. thaliana range

# which range covers A th in 1001g?
g1001<-na.omit(read.table("accessions.txt",h=F,sep=",")[,c(1,6,7)])
colnames(g1001)<-c("accession","lon","lat")
range(g1001$lon)
range(g1001$lat)
# make an extent around this range
ath_dist<-as(extent(-140, 140, 15, 69), 'SpatialPolygons')
crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
# define variables to download
variables<-c("bio3","bio6","bio7","bio8","bio10","bio15","bio16","bio17")
# make loop
for (i in 1:length(variables)) {
 #download
varras<- raster(paste0("/vsicurl/https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_",variables[i],"_1981-2010_V.2.1.tif"))
 #crop with extent
varras<-crop(varras, ath_dist)
if (i==1) {
  # first variable, store in new raster
  ath_rast<-varras
}else{
  # next variables, stack on top of previous layers
  ath_rast<-stack(ath_rast,varras)
  }
}
# Extract altitude too
elev<-raster("layers/GMTED2010.grd")
varras<-crop(elev, ath_dist)
names(varras)<-"GMTED2010_altitude"
# stack
ath_rast<-stack(ath_rast,varras)
# Visualize the layers
plot(ath_rast)
writeRaster(ath_rast, "Bioclim_Ath_1001g.grd", options="INTERLEAVE=BAND", overwrite=TRUE)
rm(list=ls())

###### 3- Prediction

# load clim data
ath_rast<-stack("Bioclim_Ath_1001g.grd")
plot(ath_rast)
# load model
load("MaxEnt_Model_Ath.Rdata")

# The predictive model cannot run on the whole map in once (at least on my computer)
# Need to crop the raster in frames of 20 degrees of longitude 
crops<-seq(-140,120,20)
# Make loop for each frame
for (i in 1:length(crops)) {
  print(paste0(i,"/",length(crops)),quote=F)
  #define new extent
  ath_dist<-as(extent( crops[i],crops[i]+20,15, 69), 'SpatialPolygons')
  crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
  #crop the raster "ath_rast" with new extent
  cropi<-crop(ath_rast,ath_dist)
  
  #predict habitat suitability on the cropped raster
  prediction <- project(me, cropi)
  
   if (i==1) {
     # first frame, store in new raster
     HS<-prediction$prediction_logistic
   }else{
     # next frames, merge aside the previous frames
     HS<-merge(HS,prediction$prediction_logistic)}
}
plot(HS)
# save the raster
writeRaster(HS, "Habitat_suitability_Ath_2023-04-27.grd", options="INTERLEAVE=BAND", overwrite=TRUE)




###### 4- Limiting factor

# load predictors
ath_rast<-stack("Bioclim_Ath_1001g.grd")
plot(ath_rast)
# load model
load("MaxEnt_Model_Ath.Rdata")

# The predictive model cannot run on the whole map in once (at least on my computer)
# Need to crop the raster in frames of 20 degrees of longitude 
crops<-seq(-140,120,20)
# Make loop for each frame
for (i in 1:length(crops)) {
  print(paste0(i,"/",length(crops)),quote=F)
  #define new extent
  ath_dist<-as(extent( crops[i],crops[i]+20,15, 69), 'SpatialPolygons')
  crs(ath_dist) <- "+proj=longlat +datum=WGS84 +no_defs"
  #crop the raster "ath_rast" with new extent
  cropi<-crop(ath_rast,ath_dist)
  
  #identify which variable is most responsible for decreasing suitability
  lim_i <- limiting(cropi,me)
  
  if (i==1) {
    # first frame, store in new raster
    lim<-lim_i
  }else{
    # next frames, merge aside the previous frames
    lim<-merge(lim,lim_i)}
}
plot(lim)
# save the raster
writeRaster(lim, "Limiting_facgtors_Ath_2023-07-07.grd", options="INTERLEAVE=BAND", overwrite=TRUE)


lim<-stack("Limiting_facgtors_Ath_2023-07-07.grd")
rasterVis::levelplot(lim,col.regions=rainbow)
names(lim)


# habitat suitability of a certain location


HS<-raster("Habitat_suitability_Ath_2023-04-27.grd")
loc<-t(data.frame(tub=c(9.057645,48.521637), mad=c(-3.703790,40.416775), mtp=c(3.862038,43.62505), spn=c(2.957075,41.72091), swe=c(13.207352,55.71226)))
colnames(loc)<-c("longitude","latitude")

coor<- SpatialPoints(loc,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# extract HS values for the coordinates and convert to dataset
#values<-data.frame(extract(HS$layer,coor,method="bilinear"))
values<-data.frame(extract(HS$layer,coor,buffer=1000,fun=max))
names(values)<-"HS"

loc
values
