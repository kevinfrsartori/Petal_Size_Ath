#------------------
# Petal Size Ath
# Manuscript figures
# Figure 3 - PINPIS BETA and FITNESS
# 2023-12-06
#------------------

# 2a - SNP effect
#----------------

snpeffect <- read.table("Genetics/bslmm_flct_Petal_Area.param.txt",h=T,sep="\t",dec=".")
snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma
# Spot allele 1 in GWAs, change sign of effect if allele 1 is not derived
Assoc <- read.table("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_Petal_Area.assoc.txt",h=T,sep="\t",dec=".")
snpeffect<-merge(snpeffect,Assoc[,c("rs","allele0","allele1","p_lrt")],by="rs",all.x = T,sort = F)
#ancestry and annot
ancestry <- na.omit(read.table("../large_files/Ath_Petal_size/tables/annotated_simple_flct_Hits_Petal_Area.iHS.AD.GO.txt",h=T,sep="\t"))

snpeffect<-merge(snpeffect,ancestry,by.x="rs",by.y="snpID",all.x=T,sort = F)
snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]<-(-1)*snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]
#candidates
candidates<-ancestry$snpID[which(ancestry$growth_dev==1)]
snpeffect<-snpeffect[order(snpeffect$growth_dev,decreasing = T),]

# Filter SNPs
# 1-Keep one snp per gene, the one with longest haplotype
snpeffect_nodup<-snpeffect[order(snpeffect$iHS_pval,decreasing = T),]
snpeffect_nodup<-snpeffect_nodup[-which(duplicated(snpeffect_nodup$geneID)),]
# 2-Keep one representative per SNP if ancestry is conserved
snpeffect_nodup<-snpeffect_nodup[-which(duplicated(snpeffect_nodup[,c("rs","iHS")])),]
# hist
hist(snpeffect_nodup$snpeffect)

# plot
par(mar=c(3.1,2,1,1),oma=c(0,0,0,0))
a<-hist(snpeffect_nodup$snpeffect,breaks = seq(-0.01,0.01,0.0005),main = "",
        ylim=c(0,6),xlim=range(-0.005,0.01),las=1,col = rgb(0,.7,.7),
      yaxt="n",ylab="",xlab="")
b<-hist(snpeffect_nodup$snpeffect[which(snpeffect_nodup$snpeffect > 0)],breaks = seq(-0.01,0.01,0.0005),
        add=T,col = rgb(.9,.7,.2))

for (i in 1:max(a$counts)) {
  segments(x0 = a$breaks[-length(a$breaks)],y0 = a$counts-i,x1 = a$breaks[-1],y1 = a$counts-i,col="black")
}
axis(side = 1,at = 0.0025,labels = "Derived allele total effect size",line = 1,tick = F)
title(main="Petal Area alleles counts")

# which allele is selected
# snpeffect.s<-snpeffect_nodup[which(snpeffect_nodup$iHS_pval>1),]
# snpeffect.s$selected<-"derived"
# x<-which(snpeffect.s$IHH_ancestral>snpeffect.s$IHH_derived)
# if(length(x)>0){snpeffect.s$selected[x]<-"ancestral"}
# # discretize
# discretize<-function(X,breaks){
#   interval<-(breaks[2]-breaks[1])/2
#   Y<-breaks[which(breaks > X)[1]]-interval
# }
# 
# breaks<-seq(-0.01,0.01,0.0005)
# x<-unlist(lapply(X = snpeffect.s$snpeffect,FUN = discretize,breaks=breaks))
# 
# y<-rep(1,length(x))
# xy<-data.frame(x=x,y=y)
# z<-which(duplicated(xy))
# while (length(z)>0) {
# xy$y[z]<-xy$y[z]+1
# z<-which(duplicated(xy))
# }
# xy$type<-factor(x = snpeffect.s$selected,levels = c("derived","ancestral"))
# points(xy$x,xy$y-.5,pch=c(24,25)[xy$type],bg=c("white","black")[xy$type],cex=2)
# 

# How to do for the other traits?

# 2b - Other traits
traits<-c("Ovule_Number",colnames(read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T))[6:16])[c(4:12,1:3)]
respie<-data.frame(trait=traits, negative=NA, positive=NA)

for (i in 1:length(traits)) {
  cat(rep("_",length(traits)),fill = T)
  cat(rep("#",i),fill=T)
  snpeffect <- read.table(paste0("Genetics/bslmm_flct_",traits[i],".param.txt"),h=T,sep="\t",dec=".")
  snpeffect$snpeffect<-snpeffect$alpha+snpeffect$beta*snpeffect$gamma
  # Spot allele 1 in GWA
  Assoc <- read.table(paste0("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_",traits[i],".assoc.txt"),h=T,sep="\t",dec=".")
  snpeffect<-merge(snpeffect,Assoc[,c("rs","allele0","allele1","p_lrt")],by="rs",all.x = T,sort = F)
  # spot ancestry
  ancestry <- na.omit(read.table(paste0("../large_files/Ath_Petal_size/tables/annotated_simple_flct_Hits_",traits[i],".iHS.AD.GO.txt"),h=T,sep="\t"))
  snpeffect<-merge(snpeffect,ancestry,by.x="rs",by.y="snpID",all.x=T,sort = F)
  snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]<-(-1)*snpeffect$snpeffect[which(!snpeffect$DER==snpeffect$allele1)]
  #candidates
  candidates<-ancestry$snpID[which(ancestry$growth_dev==1)]
  snpeffect<-snpeffect[order(snpeffect$growth_dev,decreasing = T),]
  
  # Filter SNPs
  # 1-Keep one snp per gene, the one with longest haplotype
  snpeffect_nodup<-snpeffect[order(snpeffect$iHS_pval,decreasing = T),]
  snpeffect_nodup<-snpeffect_nodup[-which(duplicated(snpeffect_nodup$geneID)),]
  # 2-Keep one representative per SNP if ancestry is conserved
  snpeffect_nodup<-snpeffect_nodup[-which(duplicated(snpeffect_nodup[,c("rs","iHS")])),]
  
  # make a dataset
  respie[i,2:3]<-table(sign(snpeffect_nodup$snpeffect))
}
# barplots
par(mar=c(.1,4,1,1),oma=c(0,0,0,0))

respie$neg_per<-(respie$negative/(respie$negative+respie$positive))*100
respie$pos_per<-(respie$positive/(respie$negative+respie$positive))*100
barplot(respie$pos_per,space = c(0,0,0,.2,0,0,.2,0,0,.2,.2,0),ylim = c(-100,100),
        yaxt="n",col = rgb(.9,.7,.2))
barplot(-respie$neg_per,space = c(0,0,0,.2,0,0,.2,0,0,.2,.2,0),
        yaxt="n",col = rgb(0,.7,.7),add=T)
axis(side = 2,at = c(-50,0,50), labels = c("50%","0%","50%"), las=1)
abline(h=c(50,-50),lty=2)
text(x=c(.5,1.5,2.5, 3.7,4.7,5.7, 6.9,7.9,8.9, 10.1, 11.3,12.3)-.3, y=(-85),srt=90,pos=4,
     labels = c("Area","Length","Width","Area","Length","Width","Area","Length","Width","","Long","Short") )

text(x = c(1.5,4.7,7.9,10.1,11.8),y=c(-85,-85,-85,-75,-85),c("Petal","Sepal","Leaf","Ovule\nNumber","Stamen"),pos=1)
segments(x0 = c(.5,3.7,6.9,11.3),y0 = -87,x1 = c(2.5,5.7,8.9,12.3),y1 = -87)


# 2c - PINPIS
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
polygon(d,col = rgb(0,.6,.6))

axis(side = 1,at = c(-3,-2,-1,0,1),labels = c(0.001,0.01,0.1,1,10))
axis(side = 1,at = -1,labels = "PiN/PiS (log scale)",line = 1,tick = F)

petal_list<-read.table("../large_files/Ath_Petal_size/tables/annotated_simple_flct_Hits_Petal_Area.iHS.AD.GO.txt",h=T,sep="\t")
candidates<-petal_list$geneID[which(petal_list$growth_dev==1)]
#petal_diff<-petal_list$gene_ID[which(petal_list$expr_petal_diff_expansion=="yes")]
allgenes<-tapply(X = log10(pinpis$PinPis[which(pinpis$Gene_ID %in% petal_list$geneID)]),
                 INDEX = pinpis$Gene_ID[which(pinpis$Gene_ID %in% petal_list$geneID)],
                 FUN =  mean)
  
d2<-density(allgenes,bw=.25)
polygon(d2$x,d2$y*.6-.1,col = rgb(0,.5,.5))

d3<-density(allgenes[which(names(allgenes) %in% candidates)],bw=.25)
polygon(d3$x,d3$y*.3-.2,col = rgb(0,.4,.4))

abline(v=mean(pinpis$logpinpis,na.rm = T),lwd=2,lty=2,col="white")

text(-4.2,0.02,"Whole genome",pos=4)
text(-4.2,-.1+0.02,"Petal area genes",pos=4)
text(-4.2,-.2+0.02,"Petal Growth Dev. genes",pos=4)


# Use genome wide distrib to compute P-value
1-pnorm(sort(allgenes), mean = mean(pinpis$logpinpis,na.rm=T), sd(pinpis$logpinpis,na.rm=T) )
gd_genes<-allgenes[which(names(allgenes) %in% candidates)]
1-pnorm(sort(gd_genes), mean = mean(pinpis$logpinpis,na.rm=T), sd(pinpis$logpinpis,na.rm=T) )

# 2d - iHS
#---------
# EHH plot from cluster

# Plot SNPs effects 275x300
library(vcfR)
# snp_1_28960616
vcf<-read.vcfR("../large_files/Ath_Petal_size/gwas/flct_Hits_Petal_Area.ann.vcf")
snp_1_28960616<-data.frame(genotypeID=colnames(vcf@gt),snp_1_28960616=vcf@gt[which(vcf@fix[,3] == "snp_1_28960616"),])[-1,]
snp_1_28960616$genotypeID<-unlist(strsplit(snp_1_28960616$genotypeID,"_"))[(1:length(snp_1_28960616$genotypeID))*2]
vcffix<-as.data.frame(vcf@fix[,1:5])
snp_1_28960616$snp_1_28960616[which(snp_1_28960616$snp_1_28960616 == "0/0")]<-vcffix[which(vcffix$ID=="snp_1_28960616"),4]
snp_1_28960616$snp_1_28960616[which(snp_1_28960616$snp_1_28960616 == "1/1")]<-vcffix[which(vcffix$ID=="snp_1_28960616"),5]
snp_1_28960616$snp_1_28960616<-factor(x = snp_1_28960616$snp_1_28960616,levels = vcffix[which(vcffix$ID=="snp_1_28960616"),4:5])

phenotypes<-read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T)
phenotypes<-merge(phenotypes,snp_1_28960616,by.x="Genotype",by.y="genotypeID")

par(mar=c(4,4,1,1))
boxplot(phenotypes$Petal_Area~phenotypes$snp_1_28960616,las=1,xlab="",
        ylab="Petal Area (mm2)",col=c(rgb(0,0,1,.5),rgb(1,0,0,.5)), )
axis(1,at = 1:2,labels = c("Ancestral","Derived"),tick = F,line = 1)


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


