#------------------
# Petal Size Ath
# Manuscript figures
# Figure 3 - PINPIS BETA and FITNESS
# 2023-12-06
#------------------

# 2a - SNP effect 600 x 200
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
par(mar=c(3.1,3.1,1,1),oma=c(0,0,0,0))
par(bg = "#EAE7D7")
a<-hist(snpeffect_nodup$snpeffect,breaks = seq(-0.01,0.01,0.0005),main = "",
        ylim=c(0,6),xlim=range(-0.005,0.01),las=1,col = rgb(.85,.85,.85),
      yaxt="n",ylab="",xlab="")
b<-hist(snpeffect_nodup$snpeffect[which(snpeffect_nodup$snpeffect > 0)],breaks = seq(-0.01,0.01,0.0005),
        add=T,col = rgb(1,1,1))

for (i in 1:max(a$counts)) {
  segments(x0 = a$breaks[-length(a$breaks)],y0 = a$counts-i,x1 = a$breaks[-1],y1 = a$counts-i,col="black")
}
axis(side = 1,at = 0.0025,labels = "Derived allele total effect size",line = 1,tick = F,cex.axis=1.25)
axis(side = 2,at = 3,labels = "Frequency",line = 1,tick = F,cex.axis=1.25)
axis(side = 2,at = 1:6-.5,labels = 1:6,las=1)

title(main="Petal Area's derived alleles effect")

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
hist(snpeffect_nodup$iHS[which(snpeffect_nodup$snpeffect>0)])
hist(snpeffect_nodup$iHS[which(snpeffect_nodup$snpeffect<0)])

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
  
  hist(snpeffect_nodup$iHS[which(snpeffect_nodup$snpeffect>0)],main = traits[i])
  
  
  # make a dataset
  respie[i,2:3]<-table(sign(snpeffect_nodup$snpeffect))
}
respie<-respie[c(1:9,11,12,10),]

# barplots 600 x 400

par(mar=c(4,1,1,1),oma=c(0,0,0,0))
par(bg = "#EAE7D7")
colors<-c("white","white","white","greenyellow","greenyellow","greenyellow","green4","green4","green4","lightblue","lightblue","lightpink")

respie$neg_per<-(respie$negative/(respie$negative+respie$positive))*100
respie$pos_per<-(respie$positive/(respie$negative+respie$positive))*100

barplot(respie$pos_per,space = c(0,0,0,.2,0,0,.2,0,0,.2,0,.2),xlim = c(-70,120),ylim=c(13,0),
        xaxt="n",col = colors,horiz = T)
axis(side = 1,at = 0,labels = "Sign of derived alleles' effect (relative abundance)",line = 1.5,tick = F,cex.axis=1.25)

barplot(-respie$neg_per,space = c(0,0,0,.2,0,0,.2,0,0,.2,0,.2),
        xaxt="n",col = colors,add=T,horiz = T)
barplot(-respie$neg_per,space = c(0,0,0,.2,0,0,.2,0,0,.2,0,.2),
        xaxt="n",col = rgb(0,0,0,.2),add=T,horiz = T)
axis(side = 1,at = c(-50,-25,0,25,50), labels = c("50%","negative","0%","positive","50%"), las=1,tick = F)
axis(side = 1,at = c(-50,0,50), labels = c("","",""))
abline(v=c(50,-50),lty=2)

text(y=c(.5,1.5,2.5, 3.7,4.7,5.7, 6.9,7.9,8.9, 10.1, 10.3,11.3), x=(80),srt=0,pos=4,
     labels = c("Area","Length","Width","Area","Length","Width","Area","Length","Width","","Long","Short") )

text(y = c(1.5,4.7,7.9,10.8,12.3)-.6,x=c(95,95,95,95,80)+15,c("Petal","Sepal","Leaf","Stamen","Ovule Number"),pos=1)


# 2c - PINPIS
#------------

pinpis<-read.table("../large_files/Ath_Petal_size/pinpis/pinpis_snp_2all_mac1_rel.95.txt",h=T,dec=".",na.strings = c("NaN","NA","-Inf","Inf"))
head(pinpis)
pinpis$logpinpis<-log10(pinpis$PinPis)
pinpis$logpinpis[which(pinpis$logpinpis=="-Inf")]<-NA
hist(pinpis$logpinpis)

#  plot 500 x 400
par(mar=c(3.1,3.1,1,1))
d<-density(x = na.omit(pinpis$logpinpis),bw = .1)
plot(d,ylim=c(-.2,.75),xlim=c(-4,2),
     main="",
     #     main = "PiN/PiS density distribution",
     yaxt="n",xaxt="n",xlab="",ylab = "")
polygon(d,col = rgb(.25,.25,.25))

axis(side = 1,at = c(-3,-2,-1,0,1),labels = c(0.001,0.01,0.1,1,10))
axis(side = 1,at = -1,labels = "PiN/PiS (log scale)",line = 1,tick = F,cex.axis=1.25)
axis(side = 2,at = .3,labels = "Density",line = 0,tick = F,cex.axis=1.25)

petal_list<-read.table("../large_files/Ath_Petal_size/tables/annotated_simple_flct_Hits_Petal_Area.iHS.AD.GO.txt",h=T,sep="\t")
candidates<-petal_list$geneID[which(petal_list$growth_dev==1)]
#petal_diff<-petal_list$gene_ID[which(petal_list$expr_petal_diff_expansion=="yes")]
allgenes<-tapply(X = log10(pinpis$PinPis[which(pinpis$Gene_ID %in% petal_list$geneID)]),
                 INDEX = pinpis$Gene_ID[which(pinpis$Gene_ID %in% petal_list$geneID)],
                 FUN =  mean)

d2<-density(allgenes,bw=.25)
polygon(d2$x,d2$y*.6-.1,col = rgb(.75,.75,.75))

d3<-density(allgenes[which(names(allgenes) %in% candidates)],bw=.25)
polygon(d3$x,d3$y*.3-.2,col = rgb(1,1,1))

abline(v=mean(pinpis$logpinpis,na.rm = T),lwd=2,lty=1,col=rgb(.25,.25,.25))
abline(v=mean(pinpis$logpinpis,na.rm = T),lwd=2,lty=2,col="white")

text(-4.2,0.02,"Whole genome",pos=4)
text(-4.2,-.1+0.02,"All petal area genes",pos=4)
text(-4.2,-.2+0.02,"Petal area candidates",pos=4)

sort(allgenes)  
sort(allgenes[which(names(allgenes) %in% candidates)])

genelist<-names(which(allgenes<0))

# Use genome wide distrib to compute P-value
1-pnorm(sort(allgenes), mean = mean(pinpis$logpinpis,na.rm=T), sd(pinpis$logpinpis,na.rm=T) )
gd_genes<-allgenes[which(names(allgenes) %in% candidates)]
1-pnorm(sort(gd_genes), mean = mean(pinpis$logpinpis,na.rm=T), sd(pinpis$logpinpis,na.rm=T) )


# all traits
traits<-c("Ovule_Number",colnames(read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T))[c(6:16,22)])[c(4:12,1:3,13)]
colors<-c("white","white","white","greenyellow","greenyellow","greenyellow","green4","green4","green4","lightpink","lightblue","lightblue","black")

d<-density(x = na.omit(pinpis$logpinpis),bw = .1)
plot(d,main = "PiN/PiS density distribution",ylim=c(-1.3,.75),xlim=c(-4,2),
     yaxt="n",xaxt="n",xlab="")
polygon(d,col = rgb(0,.6,.6))

for(i in 1:length(traits)){
list<-read.table(paste0("../large_files/Ath_Petal_size/tables/annotated_simple_flct_Hits_",traits[i],".iHS.AD.GO.txt"),h=T,sep="\t")
allgenes<-tapply(X = log10(pinpis$PinPis[which(pinpis$Gene_ID %in% list$geneID)]),
                 INDEX = pinpis$Gene_ID[which(pinpis$Gene_ID %in% list$geneID)],
                 FUN =  mean,na.rm=T)
if(any(allgenes %in% c(Inf,-Inf))){allgenes<-allgenes[-which(allgenes %in% c(Inf,-Inf))]}

if(any(is.na(allgenes))){allgenes<-allgenes[-which(is.na(allgenes))]}

d2<-density(allgenes,bw=.25)
polygon(d2$x,d2$y*.8-(.1*i),col = colors[i])

Ttest<-t.test(allgenes,na.omit(pinpis$logpinpis))
if(Ttest$p.value<0.1){ text(-4.2,-(.1*i)+0.02,paste0(traits[i]," (P < 0.1)"),pos=4) }else{ text(-4.2,-(.1*i)+0.02,traits[i],pos=4) }
}

abline(v=mean(pinpis$logpinpis,na.rm = T),lwd=2,lty=2,col="white")
library(statpsych)
test.skew(allgenes)

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


# 3c - Petal and fitness 780 x 550
#---------------------------------

layout(mat = matrix(data = c(1,1,1,1,2,1,1,1,1,3,1,1,1,1,4,1,1,1,1,5),nrow = 4,ncol = 5,byrow = T))

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
par(mar=c(9,7,0,0),oma=c(0,.5,.5,.5))
opt_gro<-summod$coefficients[1:4,1]
slope<-c(summod$coefficients[5,1],summod$coefficients[5,1]+summod$coefficients[6:8,1])
sdslope<-c(summod$coefficients[5:8,2])

plot(slope~opt_gro,las=1,xlim=c(3.4,4.5),ylim=c(-.3,.3),pch=".",xaxt="n",yaxt="n",xlab="",ylab="")

axis(side = 1,at = seq(3.4,4.4,.2), labels = seq(3.4,4.4,.2), cex.axis=1.5)
axis(side = 1,at = 4,labels = "Growing conditions optimality\n(experiment's average fitness)",tick = F,line = 4,cex.axis=2)

axis(side = 2,at = 0,labels = "Seed set ~ Petal area slope",tick = F,line = 3,cex.axis=2)
axis(side = 2,at = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3), labels = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3),las=1,cex.axis=1.5)
abline(h=0,lwd=2,col="grey",lty=2)
segments(y0 = slope-sdslope, x0 = opt_gro, y1 = slope+sdslope, x1 = opt_gro, col = "black",lwd=3)
points( opt_gro, slope, pch=21,bg="#17BEBB",cex=2.5,lwd=2)

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
segments(y0 = slope-sdslope, x0 = opt_gro, y1 = slope+sdslope, x1 = opt_gro, col = "black",lwd=3)
points( opt_gro,slope, pch=21,bg="#E4572E",cex=2.5,lwd=2)

# Przybylska et al. 2023
# OUT OF THE GRAPH
#Przybylska<-read.table("Phenotypes/rawfiles/phenotypic_datarecord.txt",h=T,sep="\t",dec=",")
#Przybylska$X1001g_ID<-as.factor(Przybylska$X1001g_ID)
#traits<-levels(as.factor(Przybylska$traitName))
#accessions<-levels(as.factor(Przybylska$X1001g_ID))
#modprzy<-na.omit(Przybylska[which(Przybylska$traitName=="fruit number"),c("X1001g_ID","traitValue","HerbivoryIndex")])
#modprzy<-na.omit(merge(modprzy,phenotypes[,c("Genotype","Petal_Area")],by.x="X1001g_ID",by.y="Genotype",all.x=T))
#modprzy$HerbivoryIndex<-as.factor(modprzy$HerbivoryIndex)
#modprzy$fitness<-log10(modprzy$traitValue*28.3)
#rm(Przybylska)
#LME
#mod<-nlme::lme(fitness ~ Petal_Area + HerbivoryIndex, random=~1|X1001g_ID, data = modprzy)
#summod<-summary(mod)
#anova(mod)
#plot(modprzy$fitness ~ modprzy$HerbivoryIndex)
#table(modprzy$HerbivoryIndex)

# LMER
#mod<-lme4::lmer(fitness ~ Petal_Area + HerbivoryIndex + (1|X1001g_ID), data = modprzy)
#summod<-summary(mod)
#summod$coefficients
#opt_gro<-summod$coefficients[1,1]
#slope<-summod$coefficients[2,1]
#points(slope, opt_gro, pch=16,col="grey")

# herbivory 2016 # NOT PUBLISHED YET
#h2016<-read.csv2("Phenotypes/rawfiles/Manip_Arabidopsis_Herbivory_data_BRUT_060916_final_1-ligne-par-pot.csv",
#                 h=T,dec=",",sep=";",stringsAsFactors=FALSE, fileEncoding="latin1",na.strings = ".")
#h2016<-na.omit(h2016[,c("Bloc","Treatment","Accession_ID","NbSiliquePlante")])
#h2016<-h2016[-which(h2016$Bloc=="A"),]
#modh2016<-na.omit(merge(h2016,phenotypes[,c("Genotype","Petal_Area")],by.x = "Accession_ID", by.y = "Genotype",all.x = T))
#modh2016$fitness<-log10(modh2016$NbSiliquePlante*28.3)
#rm(h2016)
# LME
#mod<-nlme::lme(fitness ~ Petal_Area * Treatment + Bloc, random=~1|Accession_ID, data = modh2016)
#summod<-summary(mod)
#anova(mod)
# LMER
#mod<-lme4::lmer(fitness ~ 0 + Treatment * Petal_Area + Bloc + (1|Accession_ID), data = modh2016)
#summod<-summary(mod)
#summod$coefficients
#opt_gro<-summod$coefficients[1:2,1]
#slope<-c(summod$coefficients[3,1],summod$coefficients[3,1]+summod$coefficients[6,1])
#sdslope<-c(summod$coefficients[3,2],summod$coefficients[6,2])
#segments(x0 = slope-sdslope, y0 = opt_gro, x1 = slope+sdslope, y1 = opt_gro, col = "black")
#points(slope, opt_gro, pch=16,col="green")

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
segments(y0 = slope-sdslope, x0 = opt_gro, y1 = slope+sdslope, x1 = opt_gro, col = "black",lwd=3)
points(opt_gro, slope, pch=21,bg="#FFC914",cex=2.5,lwd=2)

# mini plots for explanation
par(mar=c(1,1,1,1))

plot(0:1,0:1,type="l",xlim=c(0,1),ylim=c(-.5,1.5),lwd=2,xaxt="n",yaxt="n",xlab = "",ylab="")
plot(0:1,c(.5,.5),type="l",xlim=c(0,1),ylim=c(-.5,1.5),lwd=2,xaxt="n",yaxt="n",xlab = "",ylab="")
plot(0:1,1:0,type="l",xlim=c(0,1),ylim=c(-.5,1.5),lwd=2,xaxt="n",yaxt="n",xlab = "",ylab="")

# Bi plot for SUpp figure
#------------------------
par(mfrow=c(3,2),mar=c(4,4,1,1))

#H2016
plot(modh2016$Petal_Area,modh2016$fitness,pch=21,col="black",bg=c("green","darkgreen")[as.factor(modh2016$Treatment)],
     ylim=c(3,4.25),las=1,ylab="Fitness",xlab="Petal Area (mm2)")
mod<-nlme::lme(fitness ~ Petal_Area * Treatment + Bloc, random=~1|Accession_ID, data = modh2016)
summod<-summary(mod)
segments(min(modh2016$Petal_Area,na.rm = T),min(modh2016$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         max(modh2016$Petal_Area,na.rm = T),max(modh2016$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         lwd=2,lty=1,col="darkgreen")
segments(min(modh2016$Petal_Area,na.rm = T),min(modh2016$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[6])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[3]),
         max(modh2016$Petal_Area,na.rm = T),max(modh2016$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[6])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[3]),
         lwd=2,lty=1,col="green")
#Moises
plot(modmoises$Petal_Area,modmoises$Seeds_log,pch=21,col="black",bg=c("blue","blue4","cyan","cyan4")[as.factor(modmoises$treatment)],
     las=1,ylab="Fitness",xlab="Petal Area (mm2)")
mod<-nlme::lme(Seeds_log ~ Petal_Area * treatment , random=~1|id, data = modmoises)
summod<-summary(mod)
segments(min(modmoises$Petal_Area,na.rm = T),min(modmoises$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         max(modmoises$Petal_Area,na.rm = T),max(modmoises$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         lwd=2,lty=1,col="blue")
segments(min(modmoises$Petal_Area,na.rm = T),min(modmoises$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[6])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[3]),
         max(modmoises$Petal_Area,na.rm = T),max(modmoises$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[6])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[3]),
         lwd=2,lty=1,col="blue4")
segments(min(modmoises$Petal_Area,na.rm = T),min(modmoises$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[7])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[4]),
         max(modmoises$Petal_Area,na.rm = T),max(modmoises$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[7])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[4]),
         lwd=2,lty=1,col="cyan")
segments(min(modmoises$Petal_Area,na.rm = T),min(modmoises$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[8])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[5]),
         max(modmoises$Petal_Area,na.rm = T),max(modmoises$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[8])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[5]),
         lwd=2,lty=1,col="cyan4")
#Wilczek
plot(modWilczek$Petal_Area,modWilczek$Fitness,pch=21,col="black",bg=c("chocolate4","firebrick3","firebrick1","rosybrown1","chocolate1")[as.factor(modWilczek$country)]
     ,ylim=c(2,5),las=1,ylab="Fitness",xlab="Petal Area (mm2)")
mod<-nlme::lme(Fitness ~ Petal_Area * country, random=~1|idAccession, data = modWilczek)
summod<-summary(mod)
segments(min(modWilczek$Petal_Area,na.rm = T),min(modWilczek$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         max(modWilczek$Petal_Area,na.rm = T),max(modWilczek$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         lwd=2,lty=1,col="chocolate4")
segments(min(modWilczek$Petal_Area,na.rm = T),min(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[7])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[3]),
         max(modWilczek$Petal_Area,na.rm = T),max(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[7])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[3]),
         lwd=2,lty=1,col="firebrick3")
segments(min(modWilczek$Petal_Area,na.rm = T),min(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[8])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[4]),
         max(modWilczek$Petal_Area,na.rm = T),max(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[8])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[4]),
         lwd=2,lty=1,col="firebrick1")
segments(min(modWilczek$Petal_Area,na.rm = T),min(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[9])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[5]),
         max(modWilczek$Petal_Area,na.rm = T),max(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[9])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[5]),
         lwd=2,lty=1,col="rosybrown1")
segments(min(modWilczek$Petal_Area,na.rm = T),min(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[10])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[6]),
         max(modWilczek$Petal_Area,na.rm = T),max(modWilczek$Petal_Area,na.rm = T)*(summod$coefficients$fixed[2]+summod$coefficients$fixed[10])+(summod$coefficients$fixed[1]+summod$coefficients$fixed[6]),
         lwd=2,lty=1,col="chocolate1")
# Przy
plot(modprzy$Petal_Area,modprzy$fitness,pch=21,col="black",bg="grey",las=1,ylab="Fitness",xlab="Petal Area (mm2)")
mod<-nlme::lme(fitness ~ Petal_Area + HerbivoryIndex, random=~1|X1001g_ID, data = modprzy)
summod<-summary(mod)
# Vasseur
plot(modV2018$Petal_Area,modV2018$fitness,pch=21,col="black",bg="orange",las=1,ylab="Fitness",xlab="Petal Area (mm2)")
mod<-nlme::lme(fitness ~ Petal_Area, random=~1|idAccession, data = modV2018)
summod<-summary(mod)
segments(min(modV2018$Petal_Area,na.rm = T),min(modV2018$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         max(modV2018$Petal_Area,na.rm = T),max(modV2018$Petal_Area,na.rm = T)*summod$coefficients$fixed[2]+summod$coefficients$fixed[1],
         lwd=2,lty=1,col="orange")
# Legend
plot(rep(0,13),1:13,xlim=c(0,10),ylim=c(13,-1),cex=2,pch=21,col="black",bty="n",xaxt="n",yaxt="n",ylab="",xlab="",
     bg=c("green","darkgreen","blue","blue4","cyan","cyan4","chocolate4","firebrick3","firebrick1","rosybrown1","chocolate1","grey","orange"))
text(2,-.5,"Legend",font=2)
text(rep(0,13),1:13,pos=4,labels=c("control","herbivory","high water","low water","high water","low water","Halle Autumn","Norwich Autumn","Norwich Spring","Norwich Summer","Valencia Autumn","Przybylska et al. 2023","Vasseur et al. 2018"))
text(c(2.5,2.7,2.7),c(1.5,3.5,5.5),pos=4,labels=c("Unpublished","Madrid","Tuebingen"))
segments(x0 = 2.5,y0 = 1,x1 = 2.5,y1 = 2)
segments(x0 = 2.7,y0 = 3,x1 = 2.7,y1 = 4)
segments(x0 = 2.7,y0 = 5,x1 = 2.7,y1 = 6)
text(4.2,9,pos=4,"Wilczek et al. 2014")
segments(x0 = 4.2,y0 = 7,x1 = 4.2,y1 = 11)
text(4.5,4.5,"Exposito-Alonso et al. 2019",pos=4)
abline(h=c(.5,2.5,6.5,11.5,12.5,13.5),col="grey")
