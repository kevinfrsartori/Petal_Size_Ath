####################
#
# Obtain altitude
# KS - 2023-05-03
#
####################
getwd()

# download altitude
# Obtained with manual download from https://earthexplorer.usgs.gov/ 
# require log in with account (kevinfrs Ks.0202031312)
# go to "Digital Elevation"; then "GMTED2010"; (see https://www.usgs.gov/coastal-changes-and-impacts/gmted2010)
# resolution 30arc second, similarly to the CHELSA data sets
# N10
list_N10<-grep("N10",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_N10)) {
  print(paste0(i,"/",length(list_N10)),quote=F)
  #load tiff
  tiff<- raster(list_N10[i])
  if (i==1) {
    # first frame, store in new raster
    N10<-tiff
  }else{
    # next frames, merge aside the previous frames
    N10<-merge(tiff,N10)}
}
#N30
list_N30<-grep("N30",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_N30)) {
  print(paste0(i,"/",length(list_N30)),quote=F)
  #load tiff
  tiff<- raster(list_N30[i])
  if (i==1) {
    # first frame, store in new raster
    N30<-tiff
  }else{
    # next frames, merge aside the previous frames
    N30<-merge(tiff,N30)}
}
#N50
list_N50<-grep("N50",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_N50)) {
  print(paste0(i,"/",length(list_N50)),quote=F)
  #load tiff
  tiff<- raster(list_N50[i])
  if (i==1) {
    # first frame, store in new raster
    N50<-tiff
  }else{
    # next frames, merge aside the previous frames
    N50<-merge(tiff,N50)}
}
#N70
list_N70<-grep("N70",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_N70)) {
  print(paste0(i,"/",length(list_N70)),quote=F)
  #load tiff
  tiff<- raster(list_N70[i])
  if (i==1) {
    # first frame, store in new raster
    N70<-tiff
  }else{
    # next frames, merge aside the previous frames
    N70<-merge(tiff,N70)}
}
#S10
list_S10<-grep("S10",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_S10)) {
  print(paste0(i,"/",length(list_S10)),quote=F)
  #load tiff
  tiff<- raster(list_S10[i])
  if (i==1) {
    # first frame, store in new raster
    S10<-tiff
  }else{
    # next frames, merge aside the previous frames
    S10<-merge(tiff,S10)}
}
#S30
list_S30<-grep("S30",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_S30)) {
  print(paste0(i,"/",length(list_S30)),quote=F)
  #load tiff
  tiff<- raster(list_S30[i])
  if (i==1) {
    # first frame, store in new raster
    S30<-tiff
  }else{
    # next frames, merge aside the previous frames
    S30<-merge(tiff,S30)}
}
#S50
list_S50<-grep("S50",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_S50)) {
  print(paste0(i,"/",length(list_S50)),quote=F)
  #load tiff
  tiff<- raster(list_S50[i])
  if (i==1) {
    # first frame, store in new raster
    S50<-tiff
  }else{
    # next frames, merge aside the previous frames
    S50<-merge(tiff,S50)}
}
#S70
list_S70<-grep("S70",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_S70)) {
  print(paste0(i,"/",length(list_S70)),quote=F)
  #load tiff
  tiff<- raster(list_S70[i])
  if (i==1) {
    # first frame, store in new raster
    S70<-tiff
  }else{
    # next frames, merge aside the previous frames
    S70<-merge(tiff,S70)}
}
#S90
list_S90<-grep("S90",value = T,list.files(path = paste0("layers/altitude_GMTED2010/",(list.files("layers/altitude_GMTED2010/"))),pattern = "mea",full.names = T))
# Make loop for each frame
for (i in 1:length(list_S90)) {
  print(paste0(i,"/",length(list_S90)),quote=F)
  #load tiff
  tiff<- raster(list_S90[i])
  if (i==1) {
    # first frame, store in new raster
    S90<-tiff
  }else{
    # next frames, merge aside the previous frames
    S90<-merge(tiff,S90)}
}
rm(tiff)
#merge by pairs
N1030<-merge(N10,N30)
rm(N10,N30)
N5070<-merge(N50,N70)
rm(N50,N70)
S1030<-merge(S10,S30)
rm(S10,S30)
S5070<-merge(S50,S70)
rm(S50,S70)
S507090<-merge(S5070,S90)
rm(S90,S5070)
N<-merge(N1030,N5070)
rm(N1030,N5070)
S<-merge(S1030,S507090)
rm(S1030,S507090)

elev<-merge(N,S)
plot(elev)
rm(N,S)
writeRaster(elev, "layers/GMTED2010.grd", format="raster")

elev<-raster("layers/GMTED2010.grd")

locs<-read.table("newloc_jesse_lasky.txt",h=T)
lomin<-min(locs$NewLongitude)-1
lomax<-max(locs$NewLongitude)+1
lamin<-min(locs$NewLatitude)-1
lamax<-max(locs$NewLatitude)+1

extt<- as(extent(lomin,lomax,lamin,lamax), 'SpatialPolygons')
crs(extt) <- "+proj=longlat +datum=WGS84 +no_defs"
crpd<-crop(elev, extt)
plot(crpd)
names(crpd)<-"GMTED2010_altitude"
writeRaster(crpd, "layers/for_ath_R/altitude.grd", format="raster")



# check where are the HS outliers
elev<-raster("layers/GMTED2010.grd")
phenotypes<-read.table("C:/Users/kesi0002/GitProjects/U_shape/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
adm<-read.table("C:/Users/kesi0002/GitProjects/U_shape/GWAs/results_admixture.txt",h=T,sep = ";",dec = ".")
refs<-which(!duplicated(adm$group))
adm$group<-factor(x = adm$group,levels = adm$group[refs])
adm$couleur<-factor(x = adm$couleur,levels = adm$couleur[refs])
phenotypes<-merge(phenotypes,adm,by.x="Genotype",by.y="AthID",all=T)
pal <- colorRampPalette(c("white","black"))

plot(elev$layer,xlim=c(47,51),ylim=c(37,40),col=pal(10))
points(phenotypes$longitude,phenotypes$latitude,cex=1,lwd=2,pch=21,col="black")

plot(elev$layer,xlim=c(-30,-10),ylim=c(10,30),col=pal(10))
points(phenotypes$longitude,phenotypes$latitude,cex=1,lwd=2,pch=21,col="black")

# level zero, as if in the sea. Low value come from bad elevation data ?