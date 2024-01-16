#------------------
# Petal Size Ath
# Manuscript figures
# Figure 1 - Trait variation and covariation
# 2023-12-06
#------------------

# Fig1.A - principal component analysis 600 x 500
#--------------------------------------

phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
# remove row when missing petal area values
phenotypes<-phenotypes[-which(is.na(phenotypes$Petal_Area)),]

library(FactoMineR)
colnames(phenotypes)
datapca<-phenotypes[,c("Leaf_Area","Leaf_Length","Leaf_Width",
                       "Petal_Area","Petal_Length","Petal_Width",
                       "Sepal_Area","Sepal_Length","Sepal_Width",
                       "Long_Stamens","Short_Stamens","Ovule_Number",
                       "flowering_time")]
row.names(datapca)<-phenotypes$Genotype
library(missMDA)
nb<-estim_ncpPCA(datapca,scale=T)
plot(nb$criterion)
comp<-imputePCA(datapca,ncp=3,scale=T)
res.pca<-PCA(comp$completeObs,quanti.sup = 13)

### Coordinates
layout(mat = matrix(data = c(1,1,2,3),nrow = 2,ncol = 2,byrow = T),heights = c(.4,.6),widths = c(.4,.6))

colors<-c("green4","green4","green4","white","white","white","greenyellow","greenyellow","greenyellow","lightblue","lightblue","lightpink","black")
lett<-c("Area","Length","Width","Area","Length","Width","Area","Length","Width","Long","Short","Ovule number","Flowering time")

par(mar=c(.1,.1,2,.1),oma=c(0,0,0,0))
barplot(c(res.pca$var$cor[,1],res.pca$quanti.sup$cor[,1]),
        xaxt="n",yaxt="n",horiz = T,xlim = c(-.35,1.1),space = 0,col=colors)

axis(side = 3,at = c(0,.5,1), labels = c("","",""),line = 0)
axis(side = 3,at = c(0,.5,1), labels = c("0",".5","1"),line = -.5,cex.axis=.75,tick = F)
axis(side = 3,at = c(.25), labels = c(paste0("PC1 - ",round(res.pca$eig[1,2],1),"%")),line = 0,cex.axis=1,tick = F)

X<-c(res.pca$var$cor[,1],res.pca$quanti.sup$cor[,1])
X[which(X<0)]<-0
text(x = X,y= 1:13-.5,labels = lett,pos = 4,cex = .75)

img<-png::readPNG("phenotypes/rawfiles/images/legend_manhattan.png")
rasterImage(img,xleft = -.3, ybottom = 4, xright = -.01, ytop = 12)

par(mar=c(.1,2,.1,.1))

barplot(c(res.pca$var$cor[,2],res.pca$quanti.sup$cor[,2]),
        xaxt="n",yaxt="n",ylim=c(-.2,1.2),space=0,col=colors)
axis(side = 2,at = c(0,.5,1), labels = c("","",""),line = 0)
axis(side = 2,at = c(0,.5,1), labels = c("0",".5","1"),line = -.5,cex.axis=.75,tick = F,las=1)
axis(side = 2,at = c(.25), labels = c(paste0("PC2 - ",round(res.pca$eig[2,2],1),"%")),line = 0,cex.axis=1,tick = F)
Y<-c(res.pca$var$cor[,2],res.pca$quanti.sup$cor[,2])
Y[which(Y<0)]<-0
text(y = Y,x= (0:12),labels = lett,pos = 4,cex = .75,srt=90)

### Custom PCA - save at size 800x800
library(plotrix)
par(mar=c(0,0,0,0))
plot(0,0,xlim=c(-5.2,5.2),ylim=c(-5.2,5.2),pch=3,frame.plot = F,xlab = paste0("PC1 ",round(res.pca$eig[1,2],2),"%"),
     ylab=paste0("PC2 ",round(res.pca$eig[2,2],2),"%"),asp=1,xaxt="n",yaxt="n")
# images
img<-png::readPNG("phenotypes/rawfiles/images/5165_contour.png")
rasterImage(img,-7.5,-2.5,-2.5,2.5)
img<-png::readPNG("phenotypes/rawfiles/images/6133_contour.png")
rasterImage(img,7,-2.5,2.5,2.5)
img<-png::readPNG("phenotypes/rawfiles/images/8419_contour.png")
rasterImage(img,-2.5,6.8,2.5,2.3)
img<-png::readPNG("phenotypes/rawfiles/images/9762_contour.png")
rasterImage(img,-2.5,-7.3,2.5,-2.3)
# circle
draw.circle(0,0,5)
arrows(-6,0,6,0,length = 0.1,angle = 20)
arrows(0,-6,0,5.5,length = 0.1,angle = 20)

points(res.pca$ind$coord[,1], res.pca$ind$coord[,2],col=rgb(0,.4,.4,.5),pch=16,cex=1.5)
text(x = -2.5,y = -.3,labels = "PC1",pos=3)
text(x = 0,y = -2.5,labels = "PC2",pos=2,srt=90)

# pairwise correlation with flowering time
summary(lm(datapca$Ovule_Number ~ datapca$flowering_time))
# TO DO make table with H2 etc


# Fig1.B - phenotypes allometry 650 x 600
#----------------------------------------

# log10 and normalize traits
smadt<-data.frame(log10(comp$completeObs))

library(smatr)
layout(matrix(c(1,4,2,3),2,2,T))

# Sepal area ~ Petal area
par(mar=c(0,5,5,0))
plot(smadt$Sepal_Area ~ smadt$Petal_Area,
     ylab=expression(Sepal ~ Area ~ mm^2 ~ "(log10)"),
     pch=16,col=rgb(0,.4,.4,.5),xaxt="n",las=1)
smamod<-smatr::sma(smadt$Sepal_Area ~ smadt$Petal_Area,slope.test = 1)
p<-smamod$pval[[1]]
p<-"P-value < 0.01"
r<-round(smamod$r2[[1]],2)
s<-round(smamod$coef[[1]][2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])*0
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=4)

# Leaf Area ~ Petal area
par(mar=c(5,5,0,0))
plot(smadt$Leaf_Area ~ smadt$Petal_Area,
     ylab=expression(Leaf ~ Area ~ mm^2 ~ "(log10)"),
     xlab=expression(Petal ~ Area ~ mm^2 ~ "(log10)"),
     pch=16,col=rgb(0,.4,.4,.5),las=1)
smamod<-smatr::sma(smadt$Leaf_Area ~ smadt$Petal_Area,slope.test = 1)
p<-round(smamod$pval[[1]],2)

x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])*0
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.1
text(x,y,paste0("P-value = ",p),pos=4)

# Leaf Area ~ Setal area
par(mar=c(5,0,0,5))
plot(smadt$Leaf_Area ~ smadt$Sepal_Area,
     ylab="",yaxt="n",
     xlab=expression(Sepal ~ Area ~ mm^2 ~ "(log10)"),
     pch=16,col=rgb(0,.4,.4,.5),las=1)
smamod<-smatr::sma(smadt$Leaf_Area ~ smadt$Sepal_Area,slope.test = 1)
p<-round(smamod$pval[[1]],2)

x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])*0
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.1
text(x,y,paste0("P-value = ",p),pos=4)

# other slopes

par(mar=c(2,2,4,4))
plot(0:1,0:1,xlim=c(0,2),type='l',lty=2,col=gray(.25),
     xaxt="n",yaxt="n",ylab="",xlab="",lwd=1.5)
text(1,1,"y=x",pos=4,col=gray(.25))

smamod<-smatr::sma(smadt$Long_Stamens ~ smadt$Petal_Length,slope.test = 1)
segments(x0 = 0,y0 = smamod$coef[[1]][1,1],x1 = 1,
         y1 = smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1],lwd = 2)
text(1,smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1],"Sta ~ Pet (L)",pos=4)

smamod<-smatr::sma(smadt$Sepal_Length ~ smadt$Long_Stamens,slope.test = 1)
segments(x0 = 0,y0 = smamod$coef[[1]][1,1],x1 = 1,
         y1 = smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1])
text(1,smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1],"Sep ~ Sta (L)",pos=4)

smamod<-smatr::sma(smadt$Sepal_Area ~ smadt$Petal_Area,slope.test = 1)
segments(x0 = 0,y0 = smamod$coef[[1]][1,1],x1 = 1,
         y1 = smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1])
text(1,smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1],"Sep ~ Pet (A)",pos=4)

smamod<-smatr::sma(smadt$Sepal_Width ~ smadt$Petal_Width,slope.test = 1)
segments(x0 = 0,y0 = smamod$coef[[1]][1,1],x1 = 1,
         y1 = smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1])
text(1,smamod$coef[[1]][1,1]+smamod$coef[[1]][2,1],"Sep ~ Pet (W)",pos=4)



