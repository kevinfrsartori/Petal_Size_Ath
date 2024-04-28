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

par(mar=c(.1,.1,3,.1),oma=c(0,0,0,0))
barplot(c(res.pca$var$cor[,1],res.pca$quanti.sup$cor[,1]),
        xaxt="n",yaxt="n",horiz = T,xlim = c(-.35,1.1),space = 0,col=colors)

axis(side = 3,at = c(0,.5,1), labels = c("","",""),line = 0)
axis(side = 3,at = c(0,.5,1), labels = c("0",".5","1"),line = -.5,cex.axis=1,tick = F)
axis(side = 3,at = .5, labels = c(paste0("Trait correlation with PC1 (",round(res.pca$eig[1,2],1),"% var.)")),line = .5,cex.axis=1.25,tick = F)

X<-c(res.pca$var$cor[,1],res.pca$quanti.sup$cor[,1])
X[which(X<0)]<-0
text(x = X,y= 1:13-.5,labels = lett,pos = 4,cex = 1.25)

img<-png::readPNG("phenotypes/rawfiles/images/leaf_flower_dissect.png")
rasterImage(img,xleft = -.3, ybottom = 4.5, xright = -.01, ytop = 11.5)

par(mar=c(.1,3,.1,.1))

barplot(c(res.pca$var$cor[,2],res.pca$quanti.sup$cor[,2]),
        xaxt="n",yaxt="n",ylim=c(-.2,1.2),space=0,col=colors)
axis(side = 2,at = c(0,.5,1), labels = c("","",""),line = 0)
axis(side = 2,at = c(0,.5,1), labels = c("0",".5","1"),line = -.5,cex.axis=1,tick = F,las=1)
axis(side = 2,at = .5, labels = c(paste0("Trait correlation with PC2 (",round(res.pca$eig[2,2],1),"% var.)")),line = .5,cex.axis=1.25,tick = F)
Y<-c(res.pca$var$cor[,2],res.pca$quanti.sup$cor[,2])
Y[which(Y<0)]<-0
text(y = Y,x= (0:12),labels = lett,pos = 4,cex = 1.25,srt=90)

### Custom PCA 
library(plotrix)
par(mar=c(0,0,0,0))
plot(0,0,xlim=c(-5.2,5.2),ylim=c(-5.2,5.2),pch=3,frame.plot = F,xlab = paste0("PC1 ",round(res.pca$eig[1,2],2),"%"),
     ylab=paste0("PC2 ",round(res.pca$eig[2,2],2),"%"),asp=1,xaxt="n",yaxt="n")
# images
img<-png::readPNG("phenotypes/rawfiles/images/5165_contour_col.png")
rasterImage(img,-7.5,-2.5,-2.5,2.5)
img<-png::readPNG("phenotypes/rawfiles/images/6133_contour_col.png")
rasterImage(img,7.5,-2.5,2.5,2.5)
img<-png::readPNG("phenotypes/rawfiles/images/8419_contour_col.png")
rasterImage(img,-2.5,6.8,2.5,1.8)
img<-png::readPNG("phenotypes/rawfiles/images/9762_contour_col.png")
rasterImage(img,-2.5,-7.3,2.5,-2.3)
# circle
draw.circle(0,0,5)
arrows(-6,0,6,0,length = 0.1,angle = 20)
arrows(0,-6,0,5.5,length = 0.1,angle = 20)

points(res.pca$ind$coord[,1], res.pca$ind$coord[,2],pch=21,bg=rgb(0,0,0,.1),cex=1.5)
#text(x = -2.5,y = -.3,labels = "PC1",pos=3)
#text(x = 0,y = -2.5,labels = "PC2",pos=2,srt=90)


# Fig1.B - phenotypes allometry 650 x 600
#----------------------------------------

# log10 and normalize traits
smadt<-data.frame(log10(comp$completeObs))

library(smatr)
layout(matrix(c(1,4,2,3),2,2,T))

# Sepal area ~ Petal area
par(mar=c(0,5,5,0))
plot(smadt$Sepal_Area ~ smadt$Petal_Area,
     ylab=expression(Sepal ~ Area ~ mm^2 ~ "(log10)"),cex.lab=1.25,
     pch=21,bg=rgb(0,0,0,.1),xaxt="n",las=1,cex=1.5)
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
     ylab=expression(Leaf ~ Area ~ mm^2 ~ "(log10)"),cex.lab=1.25,
     xlab=expression(Petal ~ Area ~ mm^2 ~ "(log10)"),
     pch=21,bg=rgb(0,0,0,.1),las=1,cex=1.5)
smamod<-smatr::sma(smadt$Leaf_Area ~ smadt$Petal_Area,slope.test = 1)
p<-round(smamod$pval[[1]],2)

x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])*0
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.1
text(x,y,paste0("P-value = ",p),pos=4)

# Leaf Area ~ Setal area
par(mar=c(5,0,0,5))
plot(smadt$Leaf_Area ~ smadt$Sepal_Area,
     ylab="",yaxt="n",
     xlab=expression(Sepal ~ Area ~ mm^2 ~ "(log10)"),cex.lab=1.25,
     pch=21,bg=rgb(0,0,0,.1),las=1, cex=1.5)
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



# Correlation with other traits

#this study
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)[,c(1,5:16,22)]
names(phenotypes)[2]<-"Ovule_Number"

# 107 pheno
pheno107<-read.table("Phenotypes/rawfiles/phenotypes107_gmeans.csv",h=T,sep=",")

# Przybylska
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
rm(Przybylska,dataset.t)

# merge
btw1<-merge(phenotypes,pheno107,by.x="Genotype",by.y="accession_id",all.x=T)
btw2<-merge(phenotypes,dataset,by.x="Genotype",by.y="accession_name",all.x=T)

# Correlation within
library(corrgram)
corrgram(btw1[,2:14], order=NULL, lower.panel=panel.shade, upper.panel=NULL, text.panel=panel.txt, main="Within this study")
cormat<-cor(btw1[,2:14],use = "pairwise.complete.obs")
range(cormat[-which(cormat==1)],na.rm = T)
# significance within
#install.packages("Hmisc")
library("Hmisc")
wt_t<-0.05/length(2:14)
sigmat<-rcorr(as.matrix(btw1[,2:14]))
sigwt<-sigmat$P
cormat[which(sigwt<wt_t)]
cormat[-which(sigwt<wt_t)]<-"NS"
range(cormat[-which(cormat=="NS")],na.rm = T)

write.csv2(cormat,file = "Phenotypes/Within_this_study_correlation.csv",quote = F,row.names = T)

# Correlations btw1
library(corrgram)
corrgram(btw1[,-1], order=NULL, lower.panel=panel.shade, upper.panel=NULL, text.panel=panel.txt, main="This study versus Przybylska")
cormat<-cor(btw1[,-1],use = "pairwise.complete.obs")
corbtw1<-cormat[14:120,1:13]
range(corbtw1,na.rm = T)
corbtw1_wo_ft<-cormat[14:120,1:12]
range(corbtw1_wo_ft,na.rm = T)
# significance btw1
btw1_t<-0.05/(length(14:120)*length(1:13))

sigmat<-rcorr(as.matrix(btw1[,-1]))
sigbtw1<-sigmat$P[14:120,1:13]
corbtw1[which(sigbtw1<btw1_t)]
corbtw1[-which(sigbtw1<btw1_t)]<-"NS"
range(corbtw1[-which(corbtw1=="NS")],na.rm = T)

write.csv2(corbtw1,file = "Phenotypes/With_107pheno_correlation.csv",quote = F,row.names = T)


# Correlations btw2
library(corrgram)
corrgram(btw2[,-1], order=NULL, lower.panel=panel.shade, upper.panel=NULL, text.panel=panel.txt, main="This study versus 107 phenotypes")
cormat<-cor(btw2[,-1],use = "pairwise.complete.obs")
corbtw2<-cormat[14:29,1:13]
range(corbtw2,na.rm = T)
corbtw2_wo_ft<-cormat[14:29,1:12]
range(corbtw2_wo_ft,na.rm = T)

# significance btw2
btw2_t<-0.05/(length(14:29)*length(1:13))
sigmat<-rcorr(as.matrix(btw2[,-1]))
sigbtw2<-sigmat$P[14:29,1:13]
corbtw2[which(sigbtw2<btw2_t)]
corbtw2[-which(sigbtw2<btw2_t)]<-"NS"
range(corbtw2[-which(corbtw2=="NS")],na.rm = T)
range(corbtw2[,-13][-which(corbtw2=="NS")],na.rm = T)

write.csv2(corbtw2,file = "Phenotypes/With_Przybylska_correlation.csv",quote = F,row.names = T)
