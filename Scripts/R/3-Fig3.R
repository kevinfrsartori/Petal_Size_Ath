#------------------
# Petal Size Ath
# Manuscript figures
# Figure 3 - PINPIS BETA and FITNESS
# 2023-12-06
#------------------

pinpis<-read.table("../large_files/Ath_Petal_size/pinpis/pinpis_snp_2all_mac1_rel.95.txt",h=T,dec=".",na.strings = c("NaN","NA","-Inf","Inf"))
head(pinpis)
pinpis$logpinpis<-log10(pinpis$PinPis)
pinpis$logpinpis[which(pinpis$logpinpis=="-Inf")]<-NA
hist(pinpis$logpinpis)

# 1 - PINPIS plot
par(mar=c(3.1,1,1,1))
d<-density(x = na.omit(pinpis$logpinpis),bw = .1)
plot(d,main = "PiN/PiS density distribution",ylim=c(-.2,.75),xlim=c(-4,2),
     yaxt="n",xaxt="n",xlab="")
axis(side = 1,at = c(-3,-2,-1,0,1,2),labels = c(0.001,0.01,0.1,1,10,100))
axis(side = 1,at = -1,labels = "PiN/PiS (log scale)",line = 1,tick = F)

polygon(d,col = rgb(0,.6,.6))

petal_list<-read.table("Genetics/functionnal_annotation_Petal_Area.csv",h=T,sep=",")
candidates<-petal_list$gene_ID[which(petal_list$candidate=="yes")]
petal_diff<-petal_list$gene_ID[which(petal_list$expr_petal_diff_expansion=="yes")]

d2<-density(log10(pinpis$PinPis[which(pinpis$Gene_ID %in% petal_list$gene_ID)]),bw=.25)
polygon(d2$x,d2$y*.9-.1,col = rgb(0,.5,.5))

d3<-density(log10(pinpis$PinPis[which(pinpis$Gene_ID %in% candidates)]),bw=.25)
polygon(d3$x,d3$y*.6-.2,col = rgb(0,.4,.4))

abline(v=mean(pinpis$logpinpis,na.rm = T),lwd=2,lty=2,col="white")

text(-4.2,0.02,"Whole genome",pos=4)
text(-4.2,-.1+0.02,"Top 39 genes",pos=4)
text(-4.2,-.2+0.02,"10 relevant genes",pos=4)

# 2 - SNP effect
snpeffect <- read.table("Genetics/bslmm_top100_Leaf_Area.param.txt",h=T,sep="\t",dec=".")
snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma
hist(snpeffect$snpeffect,breaks = 20)
  
# 3 - Petal and fitness
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)

moises<-read.table("Phenotypes/rawfiles/Moises_etal_2019.csv",h=T,sep=";",dec=",")
moises<-moises[which(moises$indpop=="i"),]
moises<-moises[which(moises$water=="h"),]
moises.m<-moises[which(moises$site=="madrid"),]
# zero value for fitness seems to be a mistake
moises.m$Fitness[which(moises.m$Fitness==0)]<-NA
moises.t<-moises[which(moises$site=="tuebingen"),]
which(duplicated(moises.m$id))
which(duplicated(moises.t$id))

lietal<-read.table("Phenotypes/rawfiles/Li.et.al.2010.csv",h=T,sep=",",dec=".")

phenotypes<-merge(phenotypes,moises.m[,c("id","Fitness")],by.x = "Genotype", by.y = "id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"seeds_madrid"
phenotypes<-merge(phenotypes,moises.t[,c("id","Fitness")],by.x = "Genotype", by.y = "id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"seeds_tuebingen"
phenotypes<-merge(phenotypes,lietal[,c("accession_id","Yield.spain.2009..1st.experiment.")],by.x = "Genotype", by.y = "accession_id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"yield_spain"
phenotypes<-merge(phenotypes,lietal[,c("accession_id","Yield.sweden.2009..1st.experiment.")],by.x = "Genotype", by.y = "accession_id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"yield_sweden"
phenotypes$seeds_madrid<-phenotypes$seeds_madrid/1000
phenotypes$seeds_tuebingen<-phenotypes$seeds_tuebingen/1000

par(oma=c(3,1,1,1),mfrow=c(2,2))
# tuebingen
par(mar=c(1,4,1,1))

plot(phenotypes$seeds_tuebingen ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",xaxt="n",pch=16,col=rgb(0,.4,.4))
axis(side = 1, at = 1:4,labels = c("","","",""))
axis(side = 2, at = 10,labels = "Seed set (K)",tick = F,line = 1.25)
mod<-summary(lm(phenotypes$seeds_tuebingen ~ phenotypes$Petal_Area))

p<-mod$coefficients[2,4]
if(p<0.05){
segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coefficients[1,1] + mod$coefficients[2,1] * min(phenotypes$Petal_Area,na.rm = T),
         x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coefficients[1,1] + mod$coefficients[2,1] * max(phenotypes$Petal_Area,na.rm = T))
p<-"P-value < 0.05" }else{
p<-paste0("P-value = ",round(p,3))
}
r<-round(mod$r.squared,2)
s<-round(mod$coefficients[2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Tuebingen - HS=0.62",)

#madrid
par(mar=c(1,2,1,3))
plot(phenotypes$seeds_madrid ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",xaxt="n",pch=16,col=rgb(.8,.6,0))
axis(side = 1, at = 1:4,labels = c("","","",""))
mod<-summary(lm(phenotypes$seeds_madrid ~ phenotypes$Petal_Area))
p<-mod$coefficients[2,4]
if(p<0.05){
  segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coefficients[1,1] + mod$coefficients[2,1] * min(phenotypes$Petal_Area,na.rm = T),
           x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coefficients[1,1] + mod$coefficients[2,1] * max(phenotypes$Petal_Area,na.rm = T))
  p<-"P-value < 0.05" }else{
    p<-paste0("P-value = ",round(p,3))
  }
r<-round(mod$r.squared,2)
s<-round(mod$coefficients[2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Madrid - HS=0.52")

#sweden
par(mar=c(1,4,1,1))

plot(phenotypes$yield_sweden ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",pch=16,col=rgb(0,.4,.4))
axis(side = 2, at = .3,labels = "dry seed weight (g)",tick = F,line = 1.25)
axis(side = 1, at = 2.5,labels = expression(Petal ~ Area ~ "(" ~ mm^2 ~ ")"),tick = F,line = 1.25)

mod<-summary(lm(phenotypes$yield_sweden ~ phenotypes$Petal_Area))
p<-mod$coefficients[2,4]
if(p<0.05){
  segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coefficients[1,1] + mod$coefficients[2,1] * min(phenotypes$Petal_Area,na.rm = T),
           x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coefficients[1,1] + mod$coefficients[2,1] * max(phenotypes$Petal_Area,na.rm = T))
  p<-"P-value < 0.05" }else{
    p<-paste0("P-value = ",round(p,3))
  }
r<-round(mod$r.squared,2)
s<-round(mod$coefficients[2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Sweden - HS=0.67")

#spain
par(mar=c(1,2,1,3))

plot(phenotypes$yield_spain ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",pch=16,col=rgb(.8,.6,0))
axis(side = 1, at = 2.5,labels = expression(Petal ~ Area ~ "(" ~ mm^2 ~ ")"),tick = F,line = 1.25)
mod<-summary(lm(phenotypes$yield_spain ~ phenotypes$Petal_Area))
p<-mod$coefficients[2,4]
if(p<0.05){
  segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coefficients[1,1] + mod$coefficients[2,1] * min(phenotypes$Petal_Area,na.rm = T),
           x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coefficients[1,1] + mod$coefficients[2,1] * max(phenotypes$Petal_Area,na.rm = T))
  p<-"P-value < 0.05" }else{
    p<-paste0("P-value = ",round(p,3))
  }
r<-round(mod$r.squared,2)
s<-round(mod$coefficients[2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Spain - HS=0.43")


