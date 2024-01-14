#------------------
# Petal Size Ath
# Manuscript figures
# Figure 3 - PINPIS BETA and FITNESS
# 2023-12-06
#------------------

# 3a - SNP effect
#----------------

snpeffect <- read.table("Genetics/bslmm_top100_Petal_Area.param.txt",h=T,sep="\t",dec=".")
snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma
# Spot allele 1 in GWAs, change sign of effect if allele 1 is not derived
Assoc <- read.table("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_Petal_Area.assoc.txt",h=T,sep="\t",dec=".")
snpeffect<-merge(snpeffect,Assoc[,c("rs","allele0","allele1","p_lrt")],by="rs",all.x = T,sort = F)
ancestry <- read.table("Genetics/ADstate_Petal_Area.csv",h=T,sep=";")
snpeffect<-merge(snpeffect,ancestry,by.x="rs",by.y="snpID",all.x=T,sort = F)
snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]<-(-1)*snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]
#candidates
petal_list<-read.table("Genetics/functionnal_annotation_Petal_Area.csv",h=T,sep=",")
candidates<-petal_list$SNP[which(petal_list$candidate=="yes")]
# plot
par(mar=c(3.1,1,1,1),oma=c(0,0,0,0))
a<-hist(snpeffect$snpeffect,breaks = 20,main = "",ylim=c(0,8),las=1,col = rgb(0,.7,.7),
        yaxt="n",ylab="",xlab="")
b<-hist(snpeffect$snpeffect[which(snpeffect$rs %in% candidates)],breaks = a$breaks,add=T,col = rgb(0,.35,.35))
c<-hist(snpeffect$snpeffect[which(snpeffect$snpeffect > 0)],breaks = a$breaks,add=T,col = rgb(.9,.7,.2))
d<-hist(snpeffect$snpeffect[which(snpeffect$snpeffect > 0 & snpeffect$rs %in% candidates)],breaks = a$breaks,add=T,col = rgb(.5,.4,0))

for (i in 1:max(a$counts)) {
  segments(x0 = a$breaks[-length(a$breaks)],y0 = a$counts-i,x1 = a$breaks[-1],y1 = a$counts-i,col="black")
}
axis(side = 1,at = 0.0015,labels = "SNP total effect size",line = 1,tick = F)

# Pies for other traits
traits<-c("Ovule_Number",colnames(read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T))[6:16])[c(4:12,1:3)]
respie<-data.frame(trait=traits, negative=NA, positive=NA)

for (i in 1:length(traits)) {
snpeffect <- read.table(paste0("Genetics/bslmm_top100_",traits[i],".param.txt"),h=T,sep="\t",dec=".")
snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma
# Spot allele 1 in GWA
Assoc <- read.table(paste0("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_",traits[i],".assoc.txt"),h=T,sep="\t",dec=".")
snpeffect<-merge(snpeffect,Assoc[,c("rs","allele0","allele1","p_lrt")],by="rs",all.x = T,sort = F)
# spot ancestry
ancestry <- read.table(paste0("Genetics/ADstate_",traits[i],".csv"),h=T,sep=";")
snpeffect<-merge(snpeffect,ancestry,by.x="rs",by.y="snpID",all.x=T,sort = F)
# change sign of effect if allele 1 is not derived
snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]<-(-1)*snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]
# make a dataset
respie[i,2:3]<-table(sign(snpeffect$snpeffect))
}
#plot
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
layout(mat = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = F))
for (i in 1:dim(respie)[1]) {
pie(as.numeric(respie[i,2:3]),col = c(rgb(0,.7,.7),rgb(.9,.7,.2)),labels = "")
axis(side = 1,at = 0,labels = gsub(pattern = "_",replacement = " ",traits[i]),tick = F,line = -2)
}

# 3b - PINPIS
#------------

pinpis<-read.table("../large_files/Ath_Petal_size/pinpis/pinpis_snp_2all_mac1_rel.95.txt",h=T,dec=".",na.strings = c("NaN","NA","-Inf","Inf"))
head(pinpis)
pinpis$logpinpis<-log10(pinpis$PinPis)
pinpis$logpinpis[which(pinpis$logpinpis=="-Inf")]<-NA
hist(pinpis$logpinpis)

#  plot 500 x 400
par(mar=c(3.1,1,1,1))
d<-density(x = na.omit(pinpis$logpinpis),bw = .1)
plot(d,main = "PiN/PiS density distribution",ylim=c(-.2,.75),xlim=c(-4,2),
     yaxt="n",xaxt="n",xlab="")
axis(side = 1,at = c(-3,-2,-1,0,1),labels = c(0.001,0.01,0.1,1,10))
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


# 3c - Petal and fitness 700 x 500
#-----------------------

library(smatr)

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

hscol<-colorRampPalette(rev(viridis::inferno(4)))

par(oma=c(3,1,1,1),mfrow=c(2,2))
# tuebingen
par(mar=c(1,4,1,1))

plot(phenotypes$seeds_tuebingen ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",xaxt="n",pch=16,col=hscol(100)[62])
axis(side = 1, at = 1:4,labels = c("","","",""))
axis(side = 2, at = 10,labels = "Seed set (K)",tick = F,line = 1.25)
mod<-sma(phenotypes$seeds_tuebingen ~ phenotypes$Petal_Area)

p<-mod$pval[[1]]
if(p<0.05){
segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coef[[1]][1,1] + mod$coef[[1]][2,1] * min(phenotypes$Petal_Area,na.rm = T),
         x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coef[[1]][1,1] + mod$coef[[1]][2,1] * max(phenotypes$Petal_Area,na.rm = T))
p<-"P-value < 0.05" }else{
p<-paste0("P-value = ",round(p,3))
}
r<-round(mod$r2[[1]],2)
s<-round(mod$coef[[1]][2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Tuebingen - HS=0.62",)

#madrid
par(mar=c(1,2,1,3))
plot(phenotypes$seeds_madrid ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",xaxt="n",pch=16,col=hscol(100)[52])
axis(side = 1, at = 1:4,labels = c("","","",""))
mod<-sma(phenotypes$seeds_madrid ~ phenotypes$Petal_Area)

p<-mod$pval[[1]]
if(p<0.05){
  segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coef[[1]][1,1] + mod$coef[[1]][2,1] * min(phenotypes$Petal_Area,na.rm = T),
           x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coef[[1]][1,1] + mod$coef[[1]][2,1] * max(phenotypes$Petal_Area,na.rm = T))
  p<-"P-value < 0.05" }else{
    p<-paste0("P-value = ",round(p,3))
  }
r<-round(mod$r2[[1]],2)
s<-round(mod$coef[[1]][2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Madrid - HS=0.52")

#sweden
par(mar=c(1,4,1,1))

plot(phenotypes$yield_sweden ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",pch=16,col=hscol(100)[67])
axis(side = 2, at = .3,labels = "dry seed weight (g)",tick = F,line = 1.25)
axis(side = 1, at = 2.5,labels = expression(Petal ~ Area ~ "(" ~ mm^2 ~ ")"),tick = F,line = 1.25)

mod<-sma(phenotypes$yield_sweden ~ phenotypes$Petal_Area)
p<-mod$pval[[1]]
if(p<0.05){
  segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coefficients[1,1] + mod$coefficients[2,1] * min(phenotypes$Petal_Area,na.rm = T),
           x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coefficients[1,1] + mod$coefficients[2,1] * max(phenotypes$Petal_Area,na.rm = T))
  p<-"P-value < 0.05" }else{
    p<-paste0("P-value = ",round(p,3))
  }
r<-round(mod$r2[[1]],2)
s<-round(mod$coef[[1]][2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Sweden - HS=0.67")

#spain
par(mar=c(1,2,1,3))

plot(phenotypes$yield_spain ~ phenotypes$Petal_Area,las=1,
     ylab="",xlab="",pch=16,col=hscol(100)[43])
axis(side = 1, at = 2.5,labels = expression(Petal ~ Area ~ "(" ~ mm^2 ~ ")"),tick = F,line = 1.25)
mod<-sma(phenotypes$yield_spain ~ phenotypes$Petal_Area)
p<-mod$pval[[1]]
if(p<0.05){
  segments(x0 = min(phenotypes$Petal_Area,na.rm = T),y0 = mod$coef[[1]][1,1] + mod$coef[[1]][2,1] * min(phenotypes$Petal_Area,na.rm = T),
           x1 = max(phenotypes$Petal_Area,na.rm = T),y1 = mod$coef[[1]][1,1] + mod$coef[[1]][2,1] * max(phenotypes$Petal_Area,na.rm = T))
  p<-"P-value < 0.05" }else{
    p<-paste0("P-value = ",round(p,3))
  }
r<-round(mod$r2[[1]],2)
s<-round(mod$coef[[1]][2,1],2)
x<-par("usr")[1]+(par("usr")[2]-par("usr")[1])
y<-par("usr")[3]+(par("usr")[4]-par("usr")[3])*.8
text(x,y,paste0(" Slope = ",s,"\n R-squared = ",r,"\n ",p),pos=2)

title(main = "Spain - HS=0.43")


#Wilczek
Wilczek<-read.table("Phenotypes/rawfiles/Wilczek.et.al.2014.csv",h=T,sep=",")
Wilczek$Stock.Number<-gsub(pattern = "cs", replacement = "CS", x = Wilczek$Stock.Number)
g1001acc<-read.table("Phenotypes/rawfiles/FullList1001g.csv",h=T,sep=";")
Wilczek<-merge(Wilczek, g1001acc[,c(1,2)], by.x="Accession.Abbreviation", by.y="name",all.x=T)
Wilczek<-Wilczek[,c("idAccession",grep("Fitness",colnames(Wilczek),value = T))]
phenotypes<-merge(phenotypes,Wilczek,by.x = "Genotype", by.y = "idAccession",all.x = T)

# Norwich
plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Norwich.Autumn,pch=16,col="orange")
points(phenotypes$Petal_Area, phenotypes$Fitness.in.Norwich.Spring,pch=16,col="green")
points(phenotypes$Petal_Area, phenotypes$Fitness.in.Norwich.Summer,pch=16,col="red")

summary(lm(phenotypes$Fitness.in.Norwich.Autumn ~ phenotypes$Petal_Area))
sma(phenotypes$Fitness.in.Norwich.Spring ~ phenotypes$Petal_Area)
sma(phenotypes$Fitness.in.Norwich.Summer ~ phenotypes$Petal_Area)

which(phenotypes$Fitness.in.Halle.Autumn > 70000)
which(phenotypes$Petal_Area > 4)
phenotypes$Fitness.in.Halle.Autumn[238]<-NA
plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Halle.Autumn)
summary(lm(phenotypes$Fitness.in.Halle.Autumn ~ phenotypes$Petal_Area))

plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Valencia.Autumn)
summary(lm(phenotypes$Fitness.in.Valencia.Autumn ~ phenotypes$Petal_Area))


# not enough data for the rest of the dataset
#plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Oulu.Autumn)
#plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Halle.2007)
#plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Valencia.2007)


# try with modelisation
phenotypes$Genotype<-as.factor(phenotypes$Genotype)

mad<-na.omit(phenotypes[,c("Genotype","seeds_madrid","Petal_Area")])
colnames(mad)[2]<-"fitness"
mad$site<-"madrid"
mad$exp<-1

tub<-na.omit(phenotypes[,c("Genotype","seeds_tuebingen","Petal_Area")])
colnames(tub)[2]<-"fitness"
tub$site<-"tuebingen"
tub$exp<-1

spa<-na.omit(phenotypes[,c("Genotype","yield_spain","Petal_Area")])
colnames(spa)[2]<-"fitness"
spa$site<-"spain"
spa$exp<-2

swe<-na.omit(phenotypes[,c("Genotype","yield_sweden","Petal_Area")])
colnames(swe)[2]<-"fitness"
swe$site<-"sweden"
swe$exp<-2

newdata<-droplevels(rbind(mad,tub))
newdata$exp<-as.factor(newdata$exp)

summary(lm(fitness ~ Petal_Area + exp + site, data = newdata))

lme4::lmer(fitness ~ Petal_Area * site + (1|Genotype), data = newdata)
summary(lme4::lmer(fitness ~ Petal_Area * site + (1|Genotype), data = newdata))
