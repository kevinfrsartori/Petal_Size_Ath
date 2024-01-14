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
#---------------------------------

# Moises et al 2019 data
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
moises<-read.table("Phenotypes/rawfiles/Moises_etal_2019.csv",h=T,sep=";",dec=",")
modmoises<-na.omit(merge(moises[,c("id","site","water","Seeds")],phenotypes[,c("Genotype","Petal_Area")],by.x="id",by.y="Genotype",all.x=T))
modmoises$Seeds_log<-log10(modmoises$Seeds)
which(modmoises$Seeds_log %in% c(Inf,-Inf))
modmoises$treatment<-as.factor(paste0(modmoises$site,"_",modmoises$water,"w"))
rm(moises)
# LME moises
mod<-summary(nlme::lme(Seeds_log ~ treatment + Petal_Area, random=~1|id, data = modmoises))
anova(mod)

#LMER
mod<-lme4::lmer(Seeds_log ~ 0 + treatment * Petal_Area + (1|id), data = modmoises)
summod<-summary(mod)
summod$coefficients

# plot slopes
par(mar=c(5,5,0,0),oma=c(.5,.5,.5,.5))
opt_gro<-summod$coefficients[1:4,1]
slope<-c(summod$coefficients[5,1],summod$coefficients[5,1]+summod$coefficients[6:8,1])
sdslope<-c(summod$coefficients[5:8,2])
plot(slope,opt_gro,las=1,ylim=c(3.4,4.5),xlim=c(-.3,.05),pch=".",
     ylab="Growing conditions optimality\n(experiment's average fitness)",
     xlab="Fitness ~ Petal area slope")
abline(v=0,lwd=2,col="grey",lty=2)
segments(x0 = slope-sdslope, y0 = opt_gro, x1 = slope+sdslope, y1 = opt_gro, col = "black")
points(slope, opt_gro, pch=16,col="blue")

# Wilczek et al 2014
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
rm(g1001acc,Wilczek,Wilczek_fitness,Wilczek_NorSum,Wilczek_NorSpr,Wilczek_NorAut,Wilczek_HalAut,Wilczek_ValAut)
# LME
mod<-summary(nlme::lme(Fitness ~ Petal_Area, random=~1|idAccession, data = modWilczek))
anova(mod)

# LMER
mod<-lme4::lmer(Fitness ~ 0 + country * Petal_Area + (1|idAccession), data = modWilczek)
summod<-summary(mod)
summod$coefficients

# plot slopes
opt_gro<-summod$coefficients[1:5,1]
slope<-c(summod$coefficients[6,1],summod$coefficients[6,1]+summod$coefficients[7:10,1])
sdslope<-c(summod$coefficients[6,2],summod$coefficients[7:10,2])
segments(x0 = slope-sdslope, y0 = opt_gro, x1 = slope+sdslope, y1 = opt_gro, col = "black")
points(slope, opt_gro, pch=16,col="red")

# Przybylska et al. 2023
Przybylska<-read.table("Phenotypes/rawfiles/phenotypic_datarecord.txt",h=T,sep="\t",dec=",")
Przybylska$X1001g_ID<-as.factor(Przybylska$X1001g_ID)
traits<-levels(as.factor(Przybylska$traitName))
accessions<-levels(as.factor(Przybylska$X1001g_ID))
modprzy<-na.omit(Przybylska[which(Przybylska$traitName=="fruit number"),c("X1001g_ID","traitValue","HerbivoryIndex")])
modprzy<-na.omit(merge(modprzy,phenotypes[,c("Genotype","Petal_Area")],by.x="X1001g_ID",by.y="Genotype",all.x=T))
modprzy$HerbivoryIndex<-as.factor(modprzy$HerbivoryIndex)
modprzy$fitness<-log10(modprzy$traitValue*28.3)
rm(Przybylska)
#LME
mod<-nlme::lme(fitness ~ Petal_Area + HerbivoryIndex, random=~1|X1001g_ID, data = modprzy)
summod<-summary(mod)
anova(mod)
#plot(modprzy$fitness ~ modprzy$HerbivoryIndex)
#table(modprzy$HerbivoryIndex)

# LMER
mod<-lme4::lmer(fitness ~ Petal_Area + HerbivoryIndex + (1|X1001g_ID), data = modprzy)
summod<-summary(mod)
summod$coefficients
opt_gro<-summod$coefficients[1,1]
slope<-summod$coefficients[2,1]
points(slope, opt_gro, pch=16,col="grey")

# herbivory 2016
h2016<-read.csv2("Phenotypes/rawfiles/Manip_Arabidopsis_Herbivory_data_BRUT_060916_final_1-ligne-par-pot.csv",
                 h=T,dec=",",sep=";",stringsAsFactors=FALSE, fileEncoding="latin1",na.strings = ".")
h2016<-na.omit(h2016[,c("Bloc","Treatment","Accession_ID","NbSiliquePlante")])
h2016<-h2016[-which(h2016$Bloc=="A"),]
modh2016<-na.omit(merge(h2016,phenotypes[,c("Genotype","Petal_Area")],by.x = "Accession_ID", by.y = "Genotype",all.x = T))
modh2016$fitness<-log10(modh2016$NbSiliquePlante*28.3)
rm(h2016)
# LME
mod<-nlme::lme(fitness ~ Petal_Area * Treatment + Bloc, random=~1|Accession_ID, data = modh2016)
summod<-summary(mod)
anova(mod)
# LMER
mod<-lme4::lmer(fitness ~ 0 + Treatment * Petal_Area + Bloc + (1|Accession_ID), data = modh2016)
summod<-summary(mod)
summod$coefficients
opt_gro<-summod$coefficients[1:2,1]
slope<-c(summod$coefficients[3,1],summod$coefficients[3,1]+summod$coefficients[6,1])
sdslope<-c(summod$coefficients[3,2],summod$coefficients[6,2])
segments(x0 = slope-sdslope, y0 = opt_gro, x1 = slope+sdslope, y1 = opt_gro, col = "black")
points(slope, opt_gro, pch=16,col="green")

#vasseur 2018
V2018<-read.table("Phenotypes/rawfiles/Vasseur2018.csv",sep=",",h=T)
V2018<-na.omit(V2018[,c("idAccession","FruitNumber")])
modV2018<-na.omit(merge(V2018,phenotypes[,c("Genotype","Petal_Area")],by.x = "idAccession", by.y = "Genotype",all.x = T))
modV2018$fitness<-log10(modV2018$FruitNumber*28.3)
# LME
mod<-nlme::lme(fitness ~ Petal_Area, random=~1|idAccession, data = modV2018)
summod<-summary(mod)
anova(mod)
# LMER
mod<-lme4::lmer(fitness ~ Petal_Area + (1|idAccession), data = modV2018)
summod<-summary(mod)
summod$coefficients
opt_gro<-summod$coefficients[1,1]
slope<-summod$coefficients[2,1]
sdslope<-c(summod$coefficients[2,2])
segments(x0 = slope-sdslope, y0 = opt_gro, x1 = slope+sdslope, y1 = opt_gro, col = "black")
points(slope, opt_gro, pch=16,col="orange")





