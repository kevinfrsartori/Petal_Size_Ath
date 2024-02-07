# 3c - Petal and fitness 700 x 500
#---------------------------------


layout(mat = matrix(data=c(1,2,3,4,5,6),nrow = 2,ncol = 3,byrow = F),heights = c(5,2,5,2,5,2))
optcol<-colorRampPalette(rev(viridis::inferno(4)))

# 3.1 Li et al 2010 data
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
lietal<-read.table("Phenotypes/rawfiles/Li.et.al.2010.csv",h=T,sep=",",dec=".")
#grep("ield",colnames(lietal),value = T)

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
# convert to mg for better log transfo
modlietal$yield_log<-log10(modlietal$yield*1000)
#hist(modlietal$yield_log)
#which(modlietal$yield_log %in% c(Inf,-Inf))
modlietal$treatment<-as.factor(paste0(modlietal$country,"_",modlietal$exp))

# Fitness ~ Petal Area for significant relationship only
levels(modlietal$treatment)
l<-which(modlietal$treatment=="spain_1st")
lietal_plot<-na.omit(modlietal[l,])
par(mar=c(4,4,0,0))
plot(lietal_plot$yield ~ lietal_plot$Petal_Area, col=rgb(0,.5,.5), pch=16, las=1,
     ylab="Yield",xlab="")
axis(side = 1, at = 1.75, labels = "Petal Area (mm2)", tick = F, line = 1)
#lme
mod<-summary(nlme::lme(yield ~ Petal_Area, random=~1|accession_id, data = lietal_plot))
anova(mod)
mod$coefficients$fixed
segments(x0 = min(lietal_plot$Petal_Area),
         y0 = mod$coefficients$fixed[1]+min(lietal_plot$Petal_Area)*mod$coefficients$fixed[2],
         x1 = max(lietal_plot$Petal_Area),
         y1 = mod$coefficients$fixed[1]+max(lietal_plot$Petal_Area)*mod$coefficients$fixed[2],
         col = ,lwd = 2)

#lmer
mod<-(summary(lme4::lmer(yield ~ 0 + treatment * Petal_Area + (1|accession_id), data = modlietal)))
mod$coefficients
# reorder levels
ord<-order(mod$coefficients[1:4,1])
modlietal$treatment<-factor(x = modlietal$treatment, levels = levels(modlietal$treatment)[ord])
# plot slopes
par(mar=c(3,5,0,0))
opt_gro<-mod$coefficients[1:4,1]
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
plot(slope[ord] ~ opt_gro[ord],type="l",lwd=2,ylab="",xlab="",las=1)
axis(side = 1,at = 0.17,labels = "growth conditions' optimality",line=1,tick = F)
axis(side = 2,at = -0.020,labels = "Slope",line=2.5,tick = F)


# 3.2 Moises et al 2019 data
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
moises<-read.table("Phenotypes/rawfiles/Moises_etal_2019.csv",h=T,sep=";",dec=",")
modmoises<-merge(moises[,c("id","site","water","Seeds")],phenotypes[,c("Genotype","Petal_Area")],by.x="id",by.y="Genotype",all.x=T)
modmoises$Seeds_log<-log10(modmoises$Seeds)
which(modmoises$Seeds_log %in% c(Inf,-Inf))
#hist(modmoises$Seeds)
modmoises$treatment<-as.factor(paste0(modmoises$site,"_",modmoises$water,"w"))

# Fitness ~ Petal Area for significant relationship only
l<-which(modmoises$treatment=="tuebingen_hw")
moises_plot<-na.omit(modmoises[l,])
par(mar=c(4,4,0,0))
plot(moises_plot$Seeds_log ~ moises_plot$Petal_Area, col=rgb(0,.5,.5), pch=16, las=1,
     ylab="Seed number (log10)",xlab="",ylim=c(3.7,4.5))
axis(side = 1, at = 2.5, labels = "Petal Area (mm2)", tick = F, line = 1)
# lmer 
mod<-summary(nlme::lme(Seeds_log ~ Petal_Area, random=~1|id, data = moises_plot))
anova(mod)
mod$coefficients$fixed
segments(x0 = min(moises_plot$Petal_Area),
         y0 = mod$coefficients$fixed[1]+min(moises_plot$Petal_Area)*mod$coefficients$fixed[2],
         x1 = max(moises_plot$Petal_Area),
         y1 = mod$coefficients$fixed[1]+max(moises_plot$Petal_Area)*mod$coefficients$fixed[2],
         col = ,lwd = 2)

#lmer
mod<-summary(lme4::lmer(Seeds_log ~ 0 + treatment * Petal_Area + (1|id), data = modmoises))
mod$coefficients
# reorder levels
ord<-order(mod$coefficients[1:4,1])
modmoises$treatment<-factor(x = modmoises$treatment, levels = levels(modmoises$treatment)[ord])
# plot fit  
par(mar=c(3,5,0,0))
opt_gro<-mod$coefficients[1:4,1]
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
plot(slope[ord] ~ opt_gro[ord],type="l",lwd=2,ylab="",xlab="",las=1)
axis(side = 1,at = 3.9,labels = "growth condition optimality",line=1,tick = F)
axis(side = 2,at = -0.020,labels = "Slope",line=2.5,tick = F)


# 3.3 Wilczek et al 2014
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
Wilczek<-read.table("Phenotypes/rawfiles/Wilczek.et.al.2014.csv",h=T,sep=",")
g1001acc<-read.table("Phenotypes/rawfiles/FullList1001g.csv",h=T,sep=";")
Wilczek<-merge(Wilczek, g1001acc[,c(1,2)], by.x="Accession.Abbreviation", by.y="name",all.x=T)
Wilczek<-Wilczek[,c("idAccession",grep("Fitness",colnames(Wilczek),value = T)[1:5])]

Wilczek_fitness<-log10(as.matrix(Wilczek[,grep("Fitness",colnames(Wilczek),value = T)]))
#hist(Wilczek_fitness)
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

# Fitness ~ Petal Area for one significant relationship only
levels(modWilczek$country)
l<-which(modWilczek$country=="Norwich.Spring")
Wilczek_plot<-na.omit(modWilczek[l,])
par(mar=c(4,4,0,0))
plot(Wilczek_plot$Fitness ~ Wilczek_plot$Petal_Area, col=rgb(0,.5,.5), pch=16, las=1,
     ylab="Fitness (log10)",xlab="")
axis(side = 1, at = 2.5, labels = "Petal Area (mm2)", tick = F, line = 1)
# lmer 
mod<-summary(nlme::lme(Fitness ~ Petal_Area, random=~1|idAccession, data = Wilczek_plot))
anova(mod)
mod$coefficients$fixed
segments(x0 = min(Wilczek_plot$Petal_Area),
         y0 = mod$coefficients$fixed[1]+min(Wilczek_plot$Petal_Area)*mod$coefficients$fixed[2],
         x1 = max(Wilczek_plot$Petal_Area),
         y1 = mod$coefficients$fixed[1]+max(Wilczek_plot$Petal_Area)*mod$coefficients$fixed[2],
         col = ,lwd = 2)

#lmer
(mod<-(summary(lme4::lmer(Fitness ~ 0 + country * Petal_Area + (1|idAccession), data = modWilczek))))
mod$coefficients

ord<-order(mod$coefficients[1:5,1])
# reorder levels
modWilczek$country<-factor(x = modWilczek$country, levels = levels(modWilczek$country)[ord])

# plot slopes
par(mar=c(3,5,0,0))
opt_gro<-mod$coefficients[1:5,1]
slope<-c(mod$coefficients[6,1],mod$coefficients[6,1]+mod$coefficients[7:10,1])
plot(slope[ord] ~ opt_gro[ord],type="l",lwd=2,xlab="",ylab="",las=1)
axis(side = 1,at = 4,labels = "growth condition optimality",line=1,tick = F)
axis(side = 2,at = -0.1,labels = "Slope",line=2.5,tick = F)

#test
#moi
mod<-summary(lme4::lmer(Seeds_log ~ 0 + treatment * Petal_Area + (1|id), data = modmoises))
mod$coefficients
opt_gro<-mod$coefficients[1:4,1]
slope<-c(mod$coefficients[5,1],mod$coefficients[5,1]+mod$coefficients[6:8,1])
# plot fit  
plot(slope ~ opt_gro,ylim=c(-.2,0),xlim=c(3.25,4.6),pch=16,col="blue")

#wil
(mod<-(summary(lme4::lmer(Fitness ~ 0 + country * Petal_Area + (1|idAccession), data = modWilczek))))
mod$coefficients
opt_gro<-mod$coefficients[1:5,1]
slope<-c(mod$coefficients[6,1],mod$coefficients[6,1]+mod$coefficients[7:10,1])
points(opt_gro,slope,pch=16,col="red")


# 
# 2.3 - Przybylska et al. 2023
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
hist(dataset$fruit.number*dataset$fruit.length)
seedperfruit<-(dataset$fruit.length-mean(dataset$fruit.length,na.rm=T))/sd(dataset$fruit.length,na.rm=T)
hist(seedperfruit)
seedperfruit<-(seedperfruit*11.2)+28.3
hist(seedperfruit)
mean(seedperfruit,na.rm = T)
sd(seedperfruit,na.rm = T)
dataset$fitness<-dataset$fruit.number+seedperfruit
hist(dataset$fitness)
plot(dataset$fruit.number~dataset$fruit.length)

dataset<-merge(phenotypes, dataset, by.x="Genotype", by.y="accession_name",all.x=T)
plot(dataset$fruit.number~dataset$Petal_Area)
cor.test(dataset$fruit.number,dataset$Petal_Area)

# herbivory 2016
h2016<-read.csv2("Phenotypes/rawfiles/Manip_Arabidopsis_Herbivory_data_BRUT_060916_final_1-ligne-par-pot.csv",
                 h=T,dec=",",sep=";",stringsAsFactors=FALSE, fileEncoding="latin1",na.strings = ".")
h2016<-na.omit(h2016[,c("Bloc","Treatment","Accession_ID","NbSiliquePlante")])
h2016<-h2016[-which(h2016$Bloc=="A"),]
modh2016<-na.omit(merge(h2016,phenotypes[,c("Genotype","Petal_Area")],by.x = "Accession_ID", by.y = "Genotype",all.x = T))
modh2016$fitness<-log10(modh2016$NbSiliquePlante*28.3)
#lme
mod<-nlme::lme(fitness ~ Petal_Area * Treatment + Bloc, random=~1|Accession_ID, data = modh2016)
summod<-summary(mod)
mod<-summary(lme4::lmer(fitness ~ 0 + Treatment * Petal_Area + Bloc + (1|Accession_ID), data = modh2016))
mod$coefficients

mod$coefficients

plot((log10(modh2016$NbSiliquePlante*28.3))~modh2016$Petal_Area,pch=16,col=as.factor(modh2016$Treatment))

#vasseur 2018
# herbivory 2016
V2018<-read.table("Phenotypes/rawfiles/Vasseur2018.csv",sep=",",h=T)
V2018<-na.omit(V2018[,c("idAccession","FruitNumber")])
modV2018<-na.omit(merge(V2018,phenotypes[,c("Genotype","Petal_Area")],by.x = "idAccession", by.y = "Genotype",all.x = T))
modV2018$fitness<-log10(modV2018$FruitNumber*28.3)
#lme
mod<-nlme::lme(fitness ~ Petal_Area, random=~1|idAccession, data = modV2018)
summod<-summary(mod)
plot(modV2018$fitness ~ modV2018$Petal_Area)
mod<-summary(lme4::lmer(fitness ~ Petal_Area + (1|idAccession), data = modV2018))
mod$coefficients
