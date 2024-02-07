###################
# 
# Final report figures
# CFM - KS - 2023-05-05
#
###################

phenotypes<-read.table("U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"

# Figure 1 - Made on QGIS

#make list of studied accessions
# I keep thes one that have flower traits, since most the analysis is flower focused
# 642 accessions in total, 407 accessions with flower organ traits

ushapelist<-data.frame(famID=phenotypes$Genotype[-which(is.na(phenotypes$Petal_Area))],wifamID=phenotypes$Genotype[-which(is.na(phenotypes$Petal_Area))])
#write.table(x = ushapelist,file = "ushapelist.txt",quote = F,row.names = F,col.names = F)
# Figure 2 - table and image


## CV ####

library(dplyr)
#install.packages("formattable")
require(formattable)

# Calculate the coefficient of variation (CV) of Arabidopsis thaliana traits:
# CV = sd/mean
# Data: "phenotypes"

#View(phenotypes)
#pheno <- phenotypes[,c(2:ncol(phenotypes))]
pheno <- phenotypes[,c(5:16)]
pheno.cv <- as.data.frame(apply(pheno,2,sd,na.rm=T))
colnames(pheno.cv) <- "sd"
pheno.cv$mean <- apply(pheno,2,mean,na.rm=T)
pheno.cv$cv <- (pheno.cv$sd / pheno.cv$mean)*100 
pheno.cv <- pheno.cv %>% arrange(desc(cv))
pheno.cv <- round(pheno.cv, 3)


# pheno.cv is a new table with sd, mean and cv for every trait.
# Get the table:
colnames(pheno.cv) <- c("SD", "Mean", "CV")
formattable(pheno.cv, list(CV=color_bar(color="#71CA97", fun="proportion")),
            align =c("c","c"))

# Take away leaf parts 
remove <- c("Leaf_Area", "Leaf_Length", "Leaf_Width")
pheno.cv.flower <- pheno.cv[!rownames(pheno.cv) %in% remove,]
#View(pheno.cv.flower)

formattable(pheno.cv.flower, list(CV=color_bar(color="#71CA97", fun="proportion")),
            align =c("c","c"))

# Adding the PVE (phenotypic variation explained) from GEMMA GWAs
pheno.cv$PVE<-c(0.72, 0.1, 0.12, 0.55, 0.04, 0.85, 0.06, 0.73, 0.43, 0.27, 0.10, 0.58)
formattable(pheno.cv, list(CV=color_bar(color="#71CA97", fun="proportion")),
            align =c("c","c"))
pheno.cv <- pheno.cv %>% arrange(desc(CV))

formattable(pheno.cv[c(1,2,6),], list(CV=color_bar(color="#71CA97", fun="proportion")),
            align =c("c","c"))

## trait correlations ####
library(corrplot)

str(phenotypes)
phenotypes$Genotype <- as.character(phenotypes$Genotype)
pheno <- phenotypes[,c(14,15,16,8,9,10,11,12,13,7,6,5)]
colnames(pheno) <- c("LA", "LL", "LW", "PA", "PL", "PW", "SA", "SL", "SW", "SS", "LS", "O")

corrplot.mixed(cor(pheno, method = "pearson",use = "pairwise.complete.obs"))
corrplot(cor(pheno, method = "pearson",use = "pairwise.complete.obs"))
corrplot(cor(pheno, method = "pearson",use = "pairwise.complete.obs"), method="number", type = "lower")

## trait correlation including flowering time
pheno <- phenotypes[,c(14,8,5,22)]
colnames(pheno) <- c("LA", "PA","O","FT")
corrplot.mixed(cor(pheno, method = "pearson",use = "pairwise.complete.obs"))

## Histogram ####

# petal area:
# col=rgb(1,0.8,0.7,0.9)
# leaf area:
# col=rgb(0.2,0.8,0.5,0.5)

# Layout to split the screen
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))

# Draw the boxplot and the histogram 
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(phenotypes$Leaf_Area , horizontal=TRUE , ylim=c(0,120), xaxt="n" , col=rgb(0.2,0.8,0.5,0.5) , frame=F)
par(mar=c(4, 3.1, 1.1, 2.1))
hist(phenotypes$Leaf_Area , breaks=40 , col=rgb(0.2,0.8,0.5,0.5) , border=F , main="" , xlab=bquote(Leaf~Area~(mm^2)), xlim=c(0,120), ylim=c(0,40))
dev.off()

## pa ~ la
fit <- lm(Petal_Area ~ Leaf_Area, data= phenotypes)
plot(Petal_Area ~ Leaf_Area, data= phenotypes)
abline(fit, col="red")

summary(lm(Petal_Area ~ Leaf_Area, data= phenotypes))

## GWAS control: flowering time ####
# genetic data in use: raw 1001g with following plink filtering
#--allow-extra-chr --bfile SNP_1001g --maf 0.05 --make-bed --out SNP_1001g_filtered --snps-only just-acgt

fam<-read.table("SNP_1001g_filtered.fam")
for (i in 1:length(fam$V1)) { if(fam$V1[i] %in% phenotypes$Genotype){fam$V6[i]<-phenotypes$flowering_time[which(phenotypes$Genotype==fam$V1[i])]}}
fam$V6[which(fam$V6==-9)]<-NA
write.table(fam,file = "SNP_1001g_filtered_FT.fam",sep = "\t",quote = F,row.names = F,col.names = F)

# Run GWAs on cluster
# pve estimate 0.94

# Import data 
Assoc <- read.table("GWAs/SNP_1001g_filtered_flowering_time.assoc.txt",h=T,sep="\t",dec=".")

Assoc$Manhattan<-(-log10(Assoc$p_lrt))
Assoc$pos<-Assoc$ps
for(i in 2:5){
  Assoc$pos[Assoc$chr==i]<-max(Assoc$pos[Assoc$chr==i-1])+Assoc$pos[Assoc$chr==i]
}
Assoc$pos<-Assoc$pos/1000000
Assoc$ps[which(Assoc$Manhattan>-log10(0.05/length(Assoc$Manhattan)))]
Assoc$Manhattan[which(Assoc$Manhattan<2)]<-NA
Assoc$rs<-paste0("snp",Assoc$chr,"_",Assoc$ps)
# Hits
SNP<-which(Assoc$Manhattan>-log10(0.05/length(Assoc$Manhattan)))
FT_snp_list<-Assoc$rs[SNP]
# top 100
FT_top100<-paste0("snp_",
Assoc$chr[order(Assoc$Manhattan,decreasing = T)][1:100],"_",
Assoc$ps[order(Assoc$Manhattan,decreasing = T)][1:100])
write.table(FT_top100,file = "FT_top100.txt",quote = F,row.names = F,col.names = F)

palette(c("grey35","grey65","grey25","grey80","grey50"))
plot(Assoc$pos,Assoc$Manhattan,cex=0.5,pch=16,col=Assoc$chr,ylab="-log10(p-value)",xlab="",las=1,
     main="GWAs",ylim=c(2,-log10(0.05/length(Assoc$Manhattan))+1))
abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="red")

# FLC Chr5:3173080..3179752 (6.67 Kb)
Assoc.chr <- subset(Assoc, chr == "5")
plot(Assoc.chr$ps, Assoc.chr$Manhattan, cex=0.5, pch=16, ylab="-log10(p-value)", xlab="", las=1,
     main="GWAS flowering time. chr5",
     ylim=c(2,-log10(0.05/length(Assoc$Manhattan))+1),
     xlim=c(3172080,3180752))
# no hit in FLC

# some of the best hits:
(SNP<-Assoc$ps[which(Assoc$Manhattan>6)])

# Chr5:18589247..18591247
Assoc.chr <- subset(Assoc, chr == "5")
plot(Assoc.chr$ps, Assoc.chr$Manhattan, cex=0.5, pch=16, ylab="-log10(p-value)", xlab="", las=1,
     main="GWAS flowering time. chr5",
     ylim=c(2,-log10(0.05/length(Assoc$Manhattan))+1),
     xlim=c(18589247,18591247))
# AT5G45830 delay of germination 1

# run the local score on Flowering time
Assoc <- read.table("SNP_1001g_filtered_FT.assoc.txt",h=T,sep="\t",dec=".")
source("run_localscore.R")
sigZones05


# GWAs for other traits
traits<-colnames(phenotypes)[c(5:16,22)]
for (j in 5:16) {
  fam<-read.table("SNP_1001g_filtered.fam")
  for (i in 1:length(fam$V1)) { if(fam$V1[i] %in% phenotypes$Genotype){fam$V6[i]<-phenotypes[,j][which(phenotypes$Genotype==fam$V1[i])]}}
  fam$V6[which(fam$V6==-9)]<-NA
  write.table(fam,file = paste0("SNP_1001g_filtered_",colnames(phenotypes)[j],".fam"),sep = "\t",quote = F,row.names = F,col.names = F)
}
# Make files for all
for (i in 1:length(traits)) {
Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[i],".assoc.txt"),h=T,sep="\t",dec=".")
Assoc$Manhattan<-(-log10(Assoc$p_lrt))
SNP<-which(Assoc$Manhattan>-log10(0.05/length(Assoc$Manhattan)))
#Hits if any
snp_list<-Assoc$rs[SNP]
if (length(snp_list)>0) {
write.table(snp_list,file = paste0("Hits_",traits[i],".txt"),quote = F,row.names = F,col.names = F)}
# top 100
top100<-paste0("snp_",
                  Assoc$chr[order(Assoc$Manhattan,decreasing = T)][1:100],"_",
                  Assoc$ps[order(Assoc$Manhattan,decreasing = T)][1:100])
write.table(top100,file = paste0("top100_",traits[i],".txt"),quote = F,row.names = F,col.names = F)
source("run_localscore.R", local=TRUE)
# Sigzones
sigZones05$chr<-paste0("Chr",sigZones05$chr)
write.table(sigZones05[-which(sigZones05$peak==0),],file = paste0("Sig_zones_",traits[i],".txt"),quote = F,row.names = F,col.names = F)
}

# plot Manhattan for
traits
i<-7
dev.off()
Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[i],".assoc.txt"),h=T,sep="\t",dec=".")
Assoc$Manhattan<-(-log10(Assoc$p_lrt))
Assoc$pos<-Assoc$ps
for(i in 2:5){
  Assoc$pos[Assoc$chr==i]<-max(Assoc$pos[Assoc$chr==i-1])+Assoc$pos[Assoc$chr==i]
}
Assoc$pos<-Assoc$pos/1000000
Assoc$ps[which(Assoc$Manhattan>-log10(0.05/length(Assoc$Manhattan)))]
Assoc$Manhattan[which(Assoc$Manhattan<2)]<-NA
Assoc$rs<-paste0("snp",Assoc$chr,"_",Assoc$ps)
palette(c("grey35","grey65","grey25","grey80","grey50"))
plot(Assoc$pos,Assoc$Manhattan,cex=0.5,pch=16,col=Assoc$chr,ylab="-log10(p-value)",xlab="",las=1,
     main=paste0("GWAs ",traits[i]),ylim=c(2,-log10(0.05/length(Assoc$Manhattan))+1))
abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="red")
i<-7
Assoc$Manhattan<-(-log10(Assoc$p_lrt))
plot(Assoc$ps,Assoc$Manhattan,cex=0.5,pch=16,ylab="-log10(p-value)",xlab="",las=1,
     main=paste0("GWAs ",traits[i]),xlim=c(15627234-10000,15627234+10000),ylim=c(0,-log10(0.05/length(Assoc$Manhattan))+1))
abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="red")

plot(Assoc$ps,Assoc$Manhattan,cex=0.5,pch=16,ylab="-log10(p-value)",xlab="",las=1,
    main=paste0("GWAs ",traits[i]),xlim=c(1601862-10000,1601862+10000),ylim=c(0,-log10(0.05/length(Assoc$Manhattan))+1))
abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="red")


# Detected genes
# top100 snps were converted into genes with snpEff (linux)
phenotypes<-read.table("U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
traits<-colnames(phenotypes)[c(5:16,22)]

for (i in 1:length(traits)) {
  gene_temp<-read.table(paste0("GWAs/gene_list_",traits[i],".txt"))
  gene_temp$trait<-traits[i]
  if (i==1) { gene<-gene_temp } else { gene<-rbind(gene,gene_temp) }
}
pleiotrope<-data.frame(gene_ID=names(table(gene$V1)),times=as.vector(table(gene$V1))) 
pleiotrope<-pleiotrope[-which(pleiotrope$times==1),]
pleiotrope<-pleiotrope[order(pleiotrope$times,decreasing = T),]

# 5 times genes
# AT1G36942 - expressed specifically in flowers at early stages
gene[which(gene$V1=="AT1G36942"),]
# AT1G36950 - expressed specifically in flowers at early stages
gene[which(gene$V1=="AT1G36950"),]
# AT5G10940 - ALTERED SEED GERMINATION 2
gene[which(gene$V1=="AT5G10940"),]
# AT5G10945 - regulates vegetative phase change - microRNA that target SPL genes - SPL genes are transcription factor that is required for the initiation of both micro- and megagametogenesis and is expressed in the sporogenous tissue of the anther and the ovule
gene[which(gene$V1=="AT5G10945"),]
# AT5G10946 - expressed in flower late stage
gene[which(gene$V1=="AT5G10946"),]

# 4 times genes
# AT5G03455 - CDC25 cell cycle regulator capable of reducing the mitotic cell length, expressed in leaves and sepals
gene[which(gene$V1=="AT5G03455"),]
# AT5G03480 -  	RNA-binding specifically expressed in flowers early stages
gene[which(gene$V1=="AT5G03480"),]

# make lists of 2 times + genes
write.table(x = pleiotrope$gene_ID,file = "pleiotrope_2plus.txt",row.names = F,col.names = F,quote = F)
write.table(x = pleiotrope$gene_ID[which(pleiotrope$times>2)],file = "pleiotrope_3plus.txt",row.names = F,col.names = F,quote = F)

# GWAs hits
# petal_Width snp_1_28960616
upstream<-"AT1G77080" #AGAMOUS-LIKE 27
locatedin<-"AT1G77090" # PSII PsbP family protein
downstream<-"AT1G77093" # defensin-like (DEFL) expressed in seeds
petal_Width<-data.frame(gene_ID=c(upstream,locatedin,downstream),times="hit")
# sepal area snp_4_15627234
upstream<-"AT4G32360" # ADRENODOXIN REDUCTASE maternal gametophytic control of embryogenesis
locatedin<-"AT4G32370" # Pectin lyase-like expressed in seeds
downstream<-"AT4G32375" #  	Pectin lyase  expressed in seeds
sepal_area<-data.frame(gene_ID=c(upstream,locatedin,downstream),times="hit")
#sepal length snp_4_1601862
upstream<-"AT4G03580"# hypothetical protein expressed in flower only before opening
locatedin<-"AT4G03590" # Cystatin/monellin 
downstream<-"AT4G03600" # pyrroline-5-carboxylate reductase
sepal_length<-data.frame(gene_ID=c(upstream,locatedin,downstream),times="hit")

pleiotrope_3plus_hits<-rbind(petal_Width,sepal_area,sepal_length,pleiotrope[which(pleiotrope$times>2),])
write.table(x = pleiotrope_3plus_hits,file = "GWAs/pleiotrope_3plus_hits.txt",row.names = F,col.names = F,quote = F)

petal_Width$times<-"Petal_Area"
sepal_area$times<-"Sepal_Area"
sepal_length$times<-"Sepal_Length"
colnames(gene)<-c("gene_ID","trait")
colnames(petal_Width)<-c("gene_ID","trait")
colnames(sepal_area)<-c("gene_ID","trait")
colnames(sepal_length)<-c("gene_ID","trait")
gene<-rbind(petal_Width,sepal_area,sepal_length,gene)

write.table(x = gene,file = "GWAs/gene_trait_list.txt",row.names = F,col.names = F,quote = F)

#### See further GWAs plot in a dedicated Rscript
#################################################

## PCA clima ####

library(factoextra)
library(FactoMineR)
library(ggplot2)
library(corrplot)

ultimate <- read.csv ("Ultimate_2023-05-05.csv")
#View(ultimate)

clima <- ultimate[, c(6:76)]
res.pca.clima <- PCA(clima[, c(3:length(clima))], graph = FALSE)
print(res.pca.clima)

# See the amount of variance explained by PC (eigenvalues)
get_eigenvalue(res.pca.clima)
fviz_eig(res.pca.clima, addlabels = TRUE, ylim = c(0, 50))

# See the eigenvalues explained. Representation of the variables.
library("corrplot")
var <- get_pca_var(res.pca.clima)
corrplot(var$cos2, is.corr=FALSE)

# See PCA of the variables
fviz_pca_var(res.pca.clima, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE)

# See PCA of the genotypes
fviz_pca_ind(res.pca.clima, label="none", habillage=as.factor(ultimate$Location), addEllipses=TRUE, ellipse.level=0.95)
# We see no grouping because of phenotypical differences.



##
### color by leaf area 
##

# We want to color individuals by leaf area. Therefore, we need to assign numbers to quartiles.
# Code rows (individuals) based on whether they are below, between or above the 1st and 3rd quartile
# Get the quartiles
summary(ultimate$Leaf_Area)
# Min = 11.78
Q0<-summary(ultimate$Leaf_Area)[1]
# 1stQ = 34.06
Q1<-summary(ultimate$Leaf_Area)[2]
# 2ndQ = 42.10
Q2<-summary(ultimate$Leaf_Area)[3]
# 3rdQ = 50.15
Q3<-summary(ultimate$Leaf_Area)[5]
# Max = 116.42
Q4<-summary(ultimate$Leaf_Area)[6]
# I want 4 groups
ultimate$Leaf.Quartile <- with(ultimate,
                                ifelse(Leaf_Area <= Q1, 1,
                                       ifelse(Leaf_Area > Q1 & Leaf_Area <= Q2, 2,
                                              ifelse(Leaf_Area > Q2 & Leaf_Area <= Q3, 3,
                                                     ifelse(Leaf_Area > Q3, 4, NA)))))
ultimate$Leaf.Quartile <- factor(ultimate$Leaf.Quartile)
Leaf.Quartile <- ultimate$Leaf.Quartile

# PCA clima (dim1 & dim2) + color by quartile:
palette(viridis::inferno(6)[2:5])
fviz_pca_ind(res.pca.clima, axes = c(1,2), label="none", habillage=ultimate$Leaf.Quartile, addEllipses=TRUE, ellipse.level=0.80,palette = palette(),select.ind = list(name = rownames(res.pca.clima$ind$coord)[-which(is.na(ultimate$Leaf.Quartile))]))
# Legende:
palette(viridis::inferno(6,alpha = .5)[2:5])
boxplot(ultimate$Leaf_Area,las=2,frame=F,border = F,ylab="Leaf Area")
polygon(x = c(.65,1.35,1.35,.65),y = c(Q0,Q0,Q1,Q1),col = 1,border = F)
polygon(x = c(.65,1.35,1.35,.65),y = c(Q1,Q1,Q2,Q2),col = 2,border = F)
polygon(x = c(.65,1.35,1.35,.65),y = c(Q2,Q2,Q3,Q3),col = 3,border = F)
polygon(x = c(.65,1.35,1.35,.65),y = c(Q3,Q3,Q4,Q4),col = 4,border = F)
boxplot(ultimate$Leaf_Area,las=2,frame=F,add=T,col = rgb(1,1,1,.2),pch=16,lwd=2)

##
### color by petal area 
##

# We want to color individuals by petal area. Therefore, we need to assign numbers to quartiles just like we did for the leaf area
# Code rows (individuals) based on whether they are below, between or above the 1st and 3rd quatil
# Get the quartiles
summary(ultimate$Petal_Area)
# Min = 0.5595
Q0<-summary(ultimate$Petal_Area)[1]
# 1stQ = 1.3740
Q1<-summary(ultimate$Petal_Area)[2]
# 2ndQ = 1.7315
Q2<-summary(ultimate$Petal_Area)[3]
# 3rdQ = 2.1705
Q3<-summary(ultimate$Petal_Area)[5]
# Max = 4.3533
Q4<-summary(ultimate$Petal_Area)[6]
# I want 4 groups
ultimate$Petal.Quartile <- with(ultimate,
                                 ifelse(Petal_Area <= Q1, 1,
                                        ifelse(Petal_Area > Q1 & Petal_Area <= Q2, 2,
                                               ifelse(Petal_Area > Q2 & Petal_Area <= Q3, 3,
                                                      ifelse(Petal_Area >Q3, 4, NA)))))
ultimate$Petal.Quartile <- factor(ultimate$Petal.Quartile)
Petal.Quartile <- ultimate$Petal.Quartile

# PCA clima (dim1 & dim2) + color by quartile:
palette(viridis::inferno(6)[2:5])
fviz_pca_ind(res.pca.clima, axes = c(1, 2), label="none", habillage=as.factor(ultimate$Petal.Quartile), addEllipses=TRUE, ellipse.level=0.80,
             palette = palette(),select.ind = list(name = rownames(res.pca.clima$ind$coord)[-which(is.na(ultimate$Petal.Quartile))]))
# Legende:
palette(viridis::inferno(6,alpha = .5)[2:5])
boxplot(ultimate$Petal_Area,las=2,frame=F,border = F,ylab="Petal Area")
polygon(x = c(.65,1.35,1.35,.65),y = c(Q0,Q0,Q1,Q1),col = 1,border = F)
polygon(x = c(.65,1.35,1.35,.65),y = c(Q1,Q1,Q2,Q2),col = 2,border = F)
polygon(x = c(.65,1.35,1.35,.65),y = c(Q2,Q2,Q3,Q3),col = 3,border = F)
polygon(x = c(.65,1.35,1.35,.65),y = c(Q3,Q3,Q4,Q4),col = 4,border = F)
boxplot(ultimate$Petal_Area,las=2,frame=F,add=T,col = rgb(1,1,1,.2),pch=16,lwd=2)


##
### color by flowering time 
##

# We want to color individuals by petal area. Therefore, we need to assign numbers to quartiles just like we did for the petal area
# Code rows (individuals) based on whether they are below, between or above the 1st and 3rd quatil
# Get the quartiles
summary(ultimate$flowering_time)
Q0<-summary(ultimate$flowering_time)[1]
Q1<-summary(ultimate$flowering_time)[2]
Q2<-summary(ultimate$flowering_time)[3]
Q3<-summary(ultimate$flowering_time)[5]
Q4<-summary(ultimate$flowering_time)[6]
# I want 4 groups
ultimate$Flower.Quartile <- with(ultimate,
                                    ifelse(flowering_time <= Q1, 1,
                                           ifelse(flowering_time > Q1 & flowering_time <= Q2, 2,
                                                  ifelse(flowering_time > Q2 & flowering_time < Q3, 3,
                                                         ifelse(flowering_time >=Q3, 4, NA)))))

ultimate$Flower.Quartile <- factor(ultimate$Flower.Quartile)
Flower.Quartile <- ultimate$Flower.Quartile

# PCA clima (dim1 & dim2) + color by quartile:
fviz_pca_ind(res.pca.clima, axes = c(1, 2), label="none", habillage=as.factor(ultimate$Flower.Quartile), addEllipses=TRUE, ellipse.level=0.95)


## Colors Ad_Gr####

# I want to associate a color to all admixture group. And the color will be the one in the 
# admixture geographic map 

Colors <- data.frame(c("western_europe", "sweden", "spain", "relict", "italy_balkan_caucasus",
                       "germany", "central_europe", "asia", "admixed"))
colnames(Colors) <- "admixture_group"
col <- c("deeppink4", "blue", "orange", "yellow", "deeppink2", "turquoise", "deepskyblue", "red", "darkgreen")
Colors <- cbind(Colors, col)
write.csv(Colors, "colors.csv", row.names = F)

col_palette <- c("western_europe" = "deeppink4", "sweden" = "blue", "spain" = "orange", "relict" =  "yellow",
                 "italy_balkan_caucasus" = "deeppink2", "germany" = "turquoise", "central_europe" = "deepskyblue",
                 "asia" = "red", "admixed" = "darkgreen")


## Boxplot PA ~ admixture_group ####


##
#### Boxplot PA ~ admixture_group
##

ratios <- read.csv("/Users/claudia/Desktop/M_thesis/Data/ratios.csv")
ultimate <- read.csv("/Users/claudia/Desktop/M_thesis/Data/ultimate.csv")

# Take away the admixed group:
box_pa <- subset(ratios, Admixture_Group != "admixed")
# Add latitude information
a <- ultimate[,c("Genotype", "Latitude")]
box_pa <- left_join(box_pa, a, by="Genotype")
View(box_pa)
box_pa <- box_pa[,c(1:3, 22, 4:21)]
# Order by latitude
box_pa$Admixture_Group <- with(box_pa, reorder(Admixture_Group, Latitude, median))


# Plot
install.packages("EnvStats")
require(EnvStats)

give.n <- function(x){
  return(c(y = 1, label = paste0("n = ",length(x))))
}

ggplot(box_pa, aes(x = Admixture_Group, y = Petal_Area)) + geom_boxplot() +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, color = "black"), 
                     axis.text.y = element_text(color = "black"),  ) +
  stat_n_text(size = 3)


# Color by admixture group

col_palette <- c("western_europe" = "deeppink4", "sweden" = "blue", "spain" = "orange", "relict" =  "yellow",
                 "italy_balkan_caucasus" = "deeppink2", "germany" = "turquoise", "central_europe" = "deepskyblue",
                 "asia" = "red", "admixed" = "darkgreen")

ggplot(box_pa, aes(x = Admixture_Group, y = Petal_Area, color=Admixture_Group)) + geom_boxplot() +
  scale_color_manual(values = col_palette) +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, color = "black"), 
                     axis.text.y = element_text(color = "black"),  ) +
  stat_n_text(size = 3)


# Are they different?
# 1-way ANOVA type III (corrected by different sample sizes)
require(car)
anova <-Anova(lm(Petal_Area ~ Admixture_Group, data=box_pa), type = "III")
anova


## Do permutation to see if there are significant different petal area sizes per ad-group
require(wPerm)
perm.oneway.anova(x = box_pa$Petal_Area, y = box_pa$Admixture_Group, trim=0, ford = NULL, R = 100)


# Group Ad_Gr by Tukey
# Get the group:

library(multcompView)

model <- lm(box_pa$Petal_Area ~ box_pa$Admixture_Group)
ANOVA <- aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'box_pa$Admixture_Group', conf.level = 0.95)

generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$Admixture_Group=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Admixture_Group) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , 'box_pa$Admixture_Group')
View(LABELS)


# Plot the group

ggplot(box_pa, aes(x = Admixture_Group, y = Petal_Area, color=Admixture_Group)) + geom_boxplot() +
  scale_color_manual(values = col_palette) +
  theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1, color = "black"), 
                     axis.text.y = element_text(color = "black"),  ) +
  stat_n_text(size = 3) +
  geom_text(data = LABELS, aes(x = Admixture_Group, y = 0.4, label = Letters),
            color = "black", size = 3.2, vjust = -0.5)



## Density plots ####

ratios <- read.csv("/Users/claudia/Desktop/M_thesis/Data/ratios.csv")
View(ratios)

# to make a density plot
library(ggplot2)

col_palette <- c("western_europe" = "deeppink4", "sweden" = "blue", "spain" = "orange", "relict" =  "yellow",
                 "italy_balkan_caucasus" = "deeppink2", "germany" = "turquoise", "central_europe" = "deepskyblue",
                 "asia" = "red", "admixed" = "darkgreen")


# all separated
ggplot(data=ratios, aes(x=PaSa, group = Admixture_Group, fill=Admixture_Group)) +
  scale_fill_manual(values = col_palette) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black")) +
  xlim (0, 3.5) +
  geom_density(adjust=1.5, alpha=0.3) +
  facet_wrap(~Admixture_Group)

# remove low sample size groups
a <- subset(ratios, Admixture_Group != "relict")
a <- subset(a, Admixture_Group != "asia")


# Plot together LaPa and PaSa
install.packages("patchwork")
require(patchwork)

LaPa <- ggplot(data=a, aes(x=LaPa, group = Admixture_Group, fill=Admixture_Group)) +
  scale_fill_manual(values = col_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  geom_density(adjust=1.5, alpha=0.3) + 
  scale_x_continuous(limits = c(-20,105))


PaSa <- ggplot(data=a, aes(x=PaSa, group = Admixture_Group, fill=Admixture_Group)) +
  scale_fill_manual(values = col_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black")) +
  geom_density(adjust=1.5, alpha=0.3) + 
  scale_x_continuous(limits = c(0,3.5))

combined <- LaPa + PaSa + plot_layout(nrow=1, guides="collect") 
combined


# Difference Sw to the rest
sw <- subset(box_pa, Admixture_Group == "sweden")
rest <- subset(box_pa, Admixture_Group != "sweden")

# mean and sd LaPa sw (43 accessions)
mean(sw$LaPa)
sd(sw$LaPa)
# mean LaPa rest (286 accessions)
mean(rest$LaPa)
sd(rest$LaPa)


### T-tests
# See if the variances are equal
var.test(sw$LaPa, rest$LaPa) # Nej
# Welch two sample t-test
t.test(sw$LaPa, rest$LaPa)

# See if the variances are equal
var.test(sw$PaSa, rest$PaSa) # Yes
# Welch two sample t-test
t.test(sw$LaPa, rest$LaPa, var.equal = TRUE)

