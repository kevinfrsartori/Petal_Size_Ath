# 3c - Petal and fitness 700 x 500
#-----------------------

phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)

# Moises
#--------

moises<-read.table("Phenotypes/rawfiles/Moises_etal_2019.csv",h=T,sep=";",dec=",")
# keep individual fitness only
moises<-moises[which(moises$indpop=="i"),]
# split dataset for merge with pheno
moises.h<-moises[which(moises$water=="h"),]
moises.l<-moises[which(moises$water=="l"),]
moises.mh<-moises.h[which(moises.h$site=="madrid"),]
moises.ml<-moises.l[which(moises.l$site=="madrid"),]
moises.th<-moises.h[which(moises.h$site=="tuebingen"),]
moises.tl<-moises.l[which(moises.l$site=="tuebingen"),]
# zero value for fitness seems to be a mistake
moises.mh$Fitness[which(moises.mh$Fitness==0)]<-NA
# merge
phenotypes<-merge(phenotypes,moises.mh[,c("id","Seeds")],by.x = "Genotype", by.y = "id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"seeds_madrid_high"
phenotypes<-merge(phenotypes,moises.th[,c("id","Seeds")],by.x = "Genotype", by.y = "id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"seeds_tuebingen_high"
phenotypes<-merge(phenotypes,moises.ml[,c("id","Seeds")],by.x = "Genotype", by.y = "id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"seeds_madrid_low"
phenotypes<-merge(phenotypes,moises.tl[,c("id","Seeds")],by.x = "Genotype", by.y = "id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"seeds_tuebingen_low"

par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(4,4,0,0))
plot(phenotypes$Petal_Area, phenotypes$seeds_madrid_high,ylim=c(1000,30000),log="y",pch=16,col="red")
summary(lm(log10(phenotypes$seeds_madrid_high) ~ phenotypes$Petal_Area))


points(phenotypes$Petal_Area, phenotypes$seeds_tuebingen_high,pch=16,col="blue")
summary(lm(log10(phenotypes$seeds_tuebingen_high) ~ phenotypes$Petal_Area))

points(phenotypes$Petal_Area, phenotypes$seeds_tuebingen_low,pch=17,col="blue")
summary(lm(log10(phenotypes$seeds_tuebingen_low) ~ phenotypes$Petal_Area))

points(phenotypes$Petal_Area, phenotypes$seeds_madrid_low,pch=17,col="red")
summary(lm(log10(phenotypes$seeds_madrid_low) ~ phenotypes$Petal_Area))

# Model
modmoises<-merge(moises[,c("id","site","water","Seeds")],phenotypes[,c("Genotype","Petal_Area")],by.x="id",by.y="Genotype",all.x=T)
modmoises$Seeds_log<-log10(modmoises$Seeds)
which(modmoises$Seeds_log %in% c(Inf,-Inf))
hist(modmoises$Seeds)
summary(lm(Seeds_log ~ site * water + Petal_Area, data = modmoises))

modmoises$treatment<-as.factor(paste0(modmoises$site,"_",modmoises$water,"w"))
mod<-summary(lm(Seeds_log ~ treatment * Petal_Area, data = modmoises))
mod
opt_gro<-c(mod$coefficients[1,1],mod$coefficients[1,1]+mod$coefficients[2:4,1])
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
plot(slope ~ opt_gro)

#opt_gro<-tapply(modmoises$Seeds_log, INDEX = modmoises$treatment, FUN = median, na.rm=T)
#plot(slope ~ opt_gro)

#nlm
mod<-summary(lme4::lmer(Seeds_log ~ treatment * Petal_Area + (1|id), data = modmoises))
mod$coefficients
opt_gro<-c(mod$coefficients[1,1],mod$coefficients[1,1]+mod$coefficients[2:4,1])
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
plot(slope ~ opt_gro)

# Li et al
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)

lietal<-read.table("Phenotypes/rawfiles/Li.et.al.2010.csv",h=T,sep=",",dec=".")
grep("ield",colnames(lietal),value = T)

lietal.sp.1st<-na.omit(lietal[,c("accession_id","Yield.spain.2009..1st.experiment.")])
colnames(lietal.sp.1st)[2]<-"yield"
lietal.sp.1st$country<-"spain"
lietal.sp.1st$exp<-"1st"
lietal.sp.2nd<-na.omit(lietal[,c("accession_id","Yield.spain.2009..2nd.experiment.")])
colnames(lietal.sp.2nd)[2]<-"yield"
lietal.sp.2nd$country<-"spain"
lietal.sp.2nd$exp<-"2nd"
lietal.sw.1st<-na.omit(lietal[,c("accession_id","Yield.sweden.2009..1st.experiment.")])
colnames(lietal.sw.1st)[2]<-"yield"
lietal.sw.1st$country<-"sweden"
lietal.sw.1st$exp<-"1st"
lietal.sw.2nd<-na.omit(lietal[,c("accession_id","Yield.sweden.2009..2nd.experiment.")])
colnames(lietal.sw.2nd)[2]<-"yield"
lietal.sw.2nd$country<-"sweden"
lietal.sw.2nd$exp<-"2nd"
modlietal<-rbind(lietal.sp.1st,lietal.sp.2nd,lietal.sw.1st,lietal.sw.2nd)

modlietal<-na.omit(merge(modlietal,phenotypes[,c("Genotype","Petal_Area")],by.x = "accession_id", by.y = "Genotype",all.x = T))
#keepacc<-table(modlietal$accession_id)
#modlietal<-modlietal[which(modlietal$accession_id %in% names(keepacc)[which(keepacc>2)]),]


modlietal$yield_log<-log10(modlietal$yield*1000)
hist(modlietal$yield_log)
which(modlietal$yield_log %in% c(Inf,-Inf))
summary(lm(yield_log ~ country * exp + Petal_Area, data = modlietal))
modlietal$treatment<-as.factor(paste0(modlietal$country,"_",modlietal$exp))
(mod<-summary(lm(yield_log ~ 0 + treatment * Petal_Area, data = modlietal)))

opt_gro<-mod$coefficients[1:4,1]
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
plot(slope ~ opt_gro)

opt_gro<-tapply(modlietal$yield_log, INDEX = modlietal$treatments)
plot(slope ~ opt_gro)

#nlm
mod<-(summary(lme4::lmer(yield ~ treatment * Petal_Area + (1|accession_id), data = modlietal)))
mod$coefficients
opt_gro<-c(mod$coefficients[1,1],mod$coefficients[1,1]+mod$coefficients[2:4,1])
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
plot(slope ~ opt_gro)

#Wilczek
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)

Wilczek<-read.table("Phenotypes/rawfiles/Wilczek.et.al.2014.csv",h=T,sep=",")
g1001acc<-read.table("Phenotypes/rawfiles/FullList1001g.csv",h=T,sep=";")
Wilczek<-merge(Wilczek, g1001acc[,c(1,2)], by.x="Accession.Abbreviation", by.y="name",all.x=T)
Wilczek<-Wilczek[,c("idAccession",grep("Fitness",colnames(Wilczek),value = T)[1:5])]

Wilczek_fitness<-log10(as.matrix(Wilczek[,grep("Fitness",colnames(Wilczek),value = T)]))
hist(Wilczek_fitness)
Wilczek_fitness[which(Wilczek_fitness %in% c(Inf,-Inf))]<-NA
Wilczek[,-1]<-Wilczek_fitness

Wilczek_NorAut<-na.omit(Wilczek[,c("idAccession","Fitness.in.Norwich.Autumn")])
colnames(Wilczek_NorAut)[2]<-"Fitness"
Wilczek_NorAut$country<-"Norwich.Autumn"

Wilczek_HalAut<-na.omit(Wilczek[,c("idAccession","Fitness.in.Halle.Autumn")])
colnames(Wilczek_HalAut)[2]<-"Fitness"
Wilczek_HalAut$country<-"Halle.Autumn"

Wilczek_ValAut<-na.omit(Wilczek[,c("idAccession","Fitness.in.Valencia.Autumn")])
colnames(Wilczek_ValAut)[2]<-"Fitness"
Wilczek_ValAut$country<-"Valencia.Autumn"

Wilczek_NorSpr<-na.omit(Wilczek[,c("idAccession","Fitness.in.Norwich.Spring")])
colnames(Wilczek_NorSpr)[2]<-"Fitness"
Wilczek_NorSpr$country<-"Norwich.Spring"

Wilczek_NorSum<-na.omit(Wilczek[,c("idAccession","Fitness.in.Norwich.Summer")])
colnames(Wilczek_NorSum)[2]<-"Fitness"
Wilczek_NorSum$country<-"Norwich.Summer"

modWilczek<-rbind(Wilczek_NorAut,Wilczek_HalAut,Wilczek_ValAut,Wilczek_NorSpr,Wilczek_NorSum)
modWilczek$country<-as.factor(modWilczek$country)
modWilczek<-na.omit(merge(modWilczek,phenotypes[,c("Genotype","Petal_Area")],by.x = "idAccession", by.y = "Genotype",all.x = T))

mod<-summary(lm(Fitness ~ 0 + country * Petal_Area, data = modWilczek))
mod
opt_gro<-mod$coefficients[1:5,1]
slope<-c(mod$coefficients[6,1],mod$coefficients[6,1]+mod$coefficients[7:10,1])
plot(slope ~ opt_gro)
plot(modWilczek$Fitness ~ modWilczek$Petal_Area, col=modWilczek$country,pch=16)

opt_gro<-tapply(modWilczek$Fitness, INDEX = modWilczek$country, FUN = mean)
plot(slope ~ opt_gro)

#nlm
(mod<-(summary(lme4::lmer(Fitness ~ country * Petal_Area + (1|idAccession), data = modWilczek))))
mod$coefficients
opt_gro<-c(mod$coefficients[1,1],mod$coefficients[1,1]+mod$coefficients[2:5,1])
slope<-c(mod$coefficients[6,1],mod$coefficients[6,1]+mod$coefficients[7:10,1])
plot(slope ~ opt_gro)


colnames(phenotypes)[dim(phenotypes)[2]]<-"yield_spain_1st"
phenotypes<-merge(phenotypes,lietal[,c("accession_id","Yield.spain.2009..2nd.experiment.")],by.x = "Genotype", by.y = "accession_id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"yield_spain_2nd"

plot(phenotypes$Petal_Area, phenotypes$yield_spain_1st,pch=16)
summary(lm(phenotypes$yield_spain_1st ~ phenotypes$Petal_Area))
points(phenotypes$Petal_Area, phenotypes$yield_spain_2nd)
summary(lm(phenotypes$yield_spain_2nd ~ phenotypes$Petal_Area))


phenotypes<-merge(phenotypes,lietal[,c("accession_id","Yield.sweden.2009..1st.experiment.")],by.x = "Genotype", by.y = "accession_id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"yield_sweden_1st"
phenotypes<-merge(phenotypes,lietal[,c("accession_id","Yield.sweden.2009..2nd.experiment.")],by.x = "Genotype", by.y = "accession_id",all.x = T)
colnames(phenotypes)[dim(phenotypes)[2]]<-"yield_sweden_2nd"

plot(phenotypes$Petal_Area, phenotypes$yield_sweden_1st,pch=16)
summary(lm(phenotypes$yield_sweden_1st ~ phenotypes$Petal_Area))
points(phenotypes$Petal_Area, phenotypes$yield_sweden_2nd)
summary(lm(phenotypes$yield_sweden_2nd ~ phenotypes$Petal_Area))


#Wilczek
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)

Wilczek<-read.table("Phenotypes/rawfiles/Wilczek.et.al.2014.csv",h=T,sep=",")
g1001acc<-read.table("Phenotypes/rawfiles/FullList1001g.csv",h=T,sep=";")
Wilczek<-merge(Wilczek, g1001acc[,c(1,2)], by.x="Accession.Abbreviation", by.y="name",all.x=T)
Wilczek<-Wilczek[,c("idAccession",grep("Fitness",colnames(Wilczek),value = T)[1:5])]

#Wilczek_fitness<-log10(as.matrix(Wilczek[,grep("Fitness",colnames(Wilczek),value = T)]))
#Wilczek_fitness[which(Wilczek_fitness==-Inf)]<-NA
#Wilczek[,-1]<-Wilczek_fitness


phenotypes<-merge(phenotypes,Wilczek,by.x = "Genotype", by.y = "idAccession",all.x = T)



# Norwich
plot(phenotypes$Petal_Area, phenotypes$Fitness.in.Norwich.Autumn,pch=16,col="orange",ylim=c(1,60000))
#points(phenotypes$Petal_Area, phenotypes$Fitness.in.Norwich.Spring,pch=16,col="green")
#points(phenotypes$Petal_Area, phenotypes$Fitness.in.Norwich.Summer,pch=16,col="red")

summary(lm(phenotypes$Fitness.in.Norwich.Autumn ~ phenotypes$Petal_Area))
summary(lm(phenotypes$Fitness.in.Norwich.Spring ~ phenotypes$Petal_Area))
summary(lm(phenotypes$Fitness.in.Norwich.Summer ~ phenotypes$Petal_Area))

#which(phenotypes$Fitness.in.Halle.Autumn > 70000)
#which(phenotypes$Petal_Area > 4)
#phenotypes$Fitness.in.Halle.Autumn[238]<-NA
points(phenotypes$Petal_Area, phenotypes$Fitness.in.Halle.Autumn,pch=16,col="blue")
summary(lm(log10(phenotypes$Fitness.in.Halle.Autumn) ~ phenotypes$Petal_Area))

points(phenotypes$Petal_Area, phenotypes$Fitness.in.Valencia.Autumn,pch=16,col="black")
summary(lm(log10(phenotypes$Fitness.in.Valencia.Autumn) ~ phenotypes$Petal_Area))







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




library(missMDA)
nb = estim_ncpPCA(Wilczek[,-1],ncp.max=5)
resMI = MIPCA(Wilczek[,-1],ncp=1)
Wilczek[,-1]<-resMI$res.imputePCA

library(chemometrics)
mout<-Moutlier(Wilczek[,-1],quantile =0.95, plot=T)

outliers<-as.vector(which(mout$md>mout$cutoff))
outliers<-as.vector(which(mout$rd>mout$cutoff))
Wilczek<-Wilczek[-outliers,]


trAsh
# keep only if 4 replicates
modmoisesna<-na.omit(modmoises)
keepacc<-table(modmoisesna$id)


modmoisesna<-modmoisesna[which(modmoisesna$id %in% names(keepacc)[which(keepacc>2)]),]
summary(lme4::lmer(Seeds_log ~ treatment * Petal_Area + (1|id), data = modmoisesna))
