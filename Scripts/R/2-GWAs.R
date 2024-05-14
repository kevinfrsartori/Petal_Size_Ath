#------------------
# Petal Size Ath
# Manuscript figures
# Figure 2 - Genome wide association studies 
# 2024-01-22
#------------------


# 0- Prepare files
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
fam<-read.table("Genetics/SNP_1001g_filtered.fam")
traits<-colnames(phenotypes)[c(5:16,22)]
for (j in 1:length(traits)) {
  fam<-read.table("Genetics/SNP_1001g_filtered.fam")
  for (i in 1:length(fam$V1)) { if(fam$V1[i] %in% phenotypes$Genotype){fam$V6[i]<-phenotypes[,traits[j]][which(phenotypes$Genotype==fam$V1[i])]}}
  fam$V6[which(fam$V6==-9)]<-NA
  write.table(fam,file = paste0("Genetics/SNP_1001g_filtered_",traits[j],".fam"),sep = "\t",quote = F,row.names = F,col.names = F)
}

phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
traits<-colnames(phenotypes)[c(5:16,22)]
library(GWASTools)

# 1- GWAs were run with GEMMA (linux) see bash scripts

# 2- Outlier detection

# pre set threshold with FLC and FT

i<-13
Assoc <- read.table(paste0("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_",traits[i],".assoc.txt"),h=T,sep="\t",dec=".")
Assoc$Manhattan<-(-log10(Assoc$p_lrt))
# Bonferroni
SNP_b<-which(Assoc$Manhattan>-log10(0.05/length(Assoc$Manhattan)))
# FLC p-value
FLCassoc<-Assoc[which(Assoc$chr==5),]
FLCassoc<-FLCassoc[which(FLCassoc$ps < (3179448+10000)),]
FLCassoc<-FLCassoc[which(FLCassoc$ps > (3173497-10000)),]
plot(FLCassoc$Manhattan ~ FLCassoc$ps,las=1, xlab="Position along chromosome 5",
     ylab="-log10(P-value)",pch=16,col="black")
abline(v=c(3179448,3173497))
flc_threshold<-max(FLCassoc$Manhattan) # 4.285
FLCassoc$p_lrt[which.max(FLCassoc$Manhattan)]
SNP_flc<-which(Assoc$Manhattan>flc_threshold)


# Run snp detection for all

for (j in 1:length(traits)) {
  cat(traits[j])
  Assoc <- read.table(paste0("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  B_threshold<-(-log10(0.05/length(Assoc$Manhattan)))
  # Bonferroni
  SNP<-which(Assoc$Manhattan>=-log10(0.05/length(Assoc$Manhattan)))
  # FLC p-value
  SNP_flc<-which(Assoc$Manhattan>=flc_threshold)
  
  #Hits if any
  snp_list<-Assoc$rs[SNP]
  print(length(snp_list))
  if (length(snp_list)>0) {
    write.table(snp_list,file = paste0("Genetics/Hits_",traits[j],".txt"),quote = F,row.names = F,col.names = F,sep = "\t")}
  
  # FLC p-value hits if any
  snp_list<-Assoc$rs[SNP_flc]
  print(length(snp_list))
  snp_list<-Assoc[SNP_flc,c("rs","chr","ps")]
  if (dim(snp_list)[1]>0) {
    write.table(snp_list,file = paste0("Genetics/flct_Hits_",traits[j],".txt"),quote = F,row.names = F,col.names = F,sep = "\t")}
  
  # Filter by LD...
  # Update 16-02-2024
  #------------------
  snp_list<-Assoc$rs[SNP_4]
  print(length(snp_list))
  snp_list<-Assoc[SNP_4,c("rs","chr","ps","p_lrt")]
  snp_list$chr<-as.factor(snp_list$chr)
  snp_list$block<-NA 
  # name the blocks
  for (chr in levels(snp_list$chr )) {
    temp<-snp_list[which(snp_list$chr==chr),]
    temp<-temp[order(temp$ps),]
    for (i in 1:length(temp$ps)) {
      if (i==1) {
        b<-1
        block<-paste0(chr,".",b)
        temp$block[i]<-block 
      }else if( (temp$ps[i]-temp$ps[i-1])<10000){
        temp$block[i]<-block
      }else{
        b<-b+1
        block<-paste0(chr,".",b)
        temp$block[i]<-block
      }  
    }
    snp_list<-snp_list[-which(snp_list$chr==chr),]
    snp_list<-rbind(temp,snp_list)
  }
  #filter by plrt
  snp_list$block<-as.factor(snp_list$block)
  for (block in levels(snp_list$block )) {
    temp<-snp_list[which(snp_list$block==block),]
    temp<-temp[which.min(temp$p_lrt),]
    snp_list<-snp_list[-which(snp_list$block==block),]
    snp_list<-rbind(temp,snp_list)
  }
  if (dim(snp_list)[1]>0) {
    write.table(snp_list,file = paste0("Genetics/pvl4_LD_",traits[j],".txt"),quote = F,row.names = F,col.names = F,sep = "\t")}
  
}


# 3 - Plot Hits with all traits GWAs
#-----------------------------------


# 3.1 prepare files

hitsdir<-grep(pattern="txt",grep(pattern = "Genetics/Hits",x = list.files("Genetics/",full.names = T),value = T),value = T)
for (i in 1:length(hitsdir)) {if(i==1){hits<-read.table(hitsdir[i],sep = "_")}else{hits<-rbind(hits,read.table(hitsdir[i],sep = "_"))}}
hitsdir<-grep(pattern="gff",grep(pattern = "Genetics/Hits",x = list.files("Genetics/",full.names = T),value = T),value = T)
hitsnames<-unlist(lapply(strsplit(hitsdir,split = "\\."),"[[",1))
hitsnames<-paste0(unlist(lapply(strsplit(hitsnames,split = "_"),"[[",2)),"_",unlist(lapply(strsplit(hitsnames,split = "_"),"[[",3)))

# 3.2 make list of gff files for each hit,
# and list of list of association file for each trait in each hit 

gffs<-list()
gwas<-list()

for (i in 1:dim(hits)[1]) {
  
  #gff
  gff<-read.table(hitsdir[i])
  gff$V9<-substr(unlist(lapply(strsplit(gff$V9,split = "\\;"),"[[",1)),8,12)
  gff$V4<-gff$V4/1000000
  gff$V5<-gff$V5/1000000
  gff$V10<-gff$V4+((gff$V5-gff$V4)/2)
  gff<-gff[order(gff$V4),]
  gff$V11<-rep(c(-1.8,-2),10)[1:dim(gff)[1]]
  gffs[[i]]<-gff
  gwas[[i]]<-list()
} 
# make list of gwas output datasets around regions of interest
# This takes some time, better use "load" below

for (j in 1:length(traits)) {
  Assoc <- read.table(paste0("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$chr<-paste0("Chr",Assoc$chr)
  bonferonni=-log10(0.05/dim(Assoc)[1])
  ymax<-(bonferonni+1)
  for (i in 1:dim(hits)[1]) {
    gff<-gffs[[i]]
    Assoc.t<-Assoc[which(Assoc$chr==unique(gff$V1)),]
    Assoc.t<-Assoc.t[which(Assoc.t$ps > (hits$V3[i]-20000)),]
    Assoc.t<-Assoc.t[which(Assoc.t$ps < (hits$V3[i]+20000)),]
    Assoc.t$ps<-Assoc.t$ps/1000000
    Assoc.t$Manhattan<-(-log10(Assoc.t$p_lrt))
    gwas[[i]][[j]]<-Assoc.t
  }
}

save(gwas,file = "Genetics/Manhattan_data.Rdata")

# 3.3 Manhattan plots

load("Genetics/Manhattan_data.Rdata")
## number of analyzed SNPs/var = 540092
bonferonni<-(-log10(.05/540092))
flct<-4.285
ymax=8
traits
colors<-c("lightpink","lightblue","lightblue","white","white","white","greenyellow","greenyellow","greenyellow","green4","green4","green4","black")
pch<-c(23,22,22,21,22,24,21,22,24,21,22,24,23)
mars<-list(c(5,5,2,0),c(5,2.5,2,2.5),c(5,0,2,5))
par(mfrow=c(1,3),oma=c(0,0,0,0))

for (i in 1:dim(hits)[1]) {
  par(mar=mars[[i]])
  
  for (j in 1:length(traits)) {
    
    Assoc <- gwas[[i]][[j]]
    if (j==1) {
      plot(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=pch[j],bg=colors[j],ylab="-log10(p-value)",
           xlab=paste0("Position along Chromosome ",hits[i,2]," (Mbp)"),las=1,
           main=paste0("Manhattan plot at locus ",i),ylim=c(-2,ymax),yaxt="n",
           xlim=c(min(hits[i,3]/1000000)-.02,max(hits[i,3]/1000000)+.02))
    }else{points(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=pch[j],bg=colors[j])}
  }
  axis(side = 2,at = c(0,2,4,6,8),labels = c(0,2,4,6,8),las=1)
  abline(h=bonferonni,lty=2,col="darkred")
  abline(h=flct,lty=2,col="darkblue")
  abline(h=-1,lty=1,col="black")
  gff<-gffs[[i]]
  for (k in 1:dim(gff)[1]) {
    polygon(x=c(gff[k,4],gff[k,4],gff[k,5],gff[k,5]),y=c(-.5,-1.5,-1.5,-.5),col = rgb(0,.5,.5,1))
    text(gff[k,10],gff[k,11],gff[k,9],cex = .75)
  }
}

# 4- multivariate models
# 4-1- load data
# this study
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)[,c(1,5:16,22)]
names(phenotypes)[2]<-"Ovule_Number"
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
btw<-merge(phenotypes,dataset,by.x="Genotype",by.y="accession_name",all.x=T)
# Keep only highly pairwise correlated ones (without FT related traits)
btw<-btw[,c("Genotype","Ovule_Number", "Long_Stamens" , "Petal_Area", "Sepal_Area", "Leaf_Area", "leaf.area.per.leaf.dry.mass", "leaf.nitrogen.content.per.leaf.dry.mass", "whole.plant.dry.mass")]
# make fam file
traits<-colnames(btw)[-1]
fam<-read.table("Genetics/SNP_1001g_filtered.fam",na.strings = "-9")

for (j in 1:(dim(btw)[2]-1)) {
  fam[,5+j]<-NA
  for (i in 1:length(fam$V1)) { 
    if(fam$V1[i] %in% btw$Genotype){
      fam[i,5+j]<-btw[,traits[j]][which(btw$Genotype==fam$V1[i])]
    }
  }
}

write.table(fam,file = "Genetics/SNP_1001g_filtered_multi.fam",sep = "\t",quote = F,row.names = F,col.names = F)
comb<-t(combn(1:8,2))
write.table(comb,file = "Genetics/Pairwise_traits_comb.txt",quote = F,row.names = F,col.names = F,sep = "\t")

# 4-2- See bash scripts for running multivariate models

# 4-3- Compute the genetic correlation matrices

traits<-c("Ovule_Number", "Long_Stamens" , "Petal_Area", "Sepal_Area", "Leaf_Area", "leaf.area.per.leaf.dry.mass", "leaf.nitrogen.content.per.leaf.dry.mass", "whole.plant.dry.mass")
# correlation from covariance matrix function
# data
covmat<-as.matrix(read.table("Genetics/REML_multi.txt"))
colnames(covmat)<-traits
rownames(covmat)<-traits
covmat[upper.tri(covmat)]<-t(covmat)[upper.tri(t(covmat))]
# Convert
cormat<-cov2cor(covmat)
cormat.r<-round(cormat,2)
# Se
semat<-as.matrix(read.table("Genetics/SE_REML_multi.txt"))
diag(semat)<-diag(covmat)
semat<-cov2cor(semat)
#z-score
Z<-cormat/semat
#P-value
P<-2*pnorm(abs(Z),mean = 0,sd = 1,lower.tail = F)
diag(P)<-NA
Pmat<-matrix(data = NA,nrow = 8,ncol = 8,byrow = T,dimnames = dimnames(cormat))
Pmat[which(P>0.05)]<-"NS"
Pmat[which(P<0.01)]<-"<0.01***"
which(P>0.01 & P<0.05)
diag(Pmat)<-""
Pmat[upper.tri(Pmat)]<-cormat.r[upper.tri(cormat.r)]
write.csv2(Pmat,file = "Genetics/Genetic_correlation_sig.csv",quote = F,row.names = T)



# 5 - VENN DIAGRAMS
#------------------
library(VennDiagram)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# Petal sepal ovule stamen leaf
# Venn diagram - SNPs

threshold<-c("flct","pvl4")[1]

petal<-grep(pattern = "Petal",grep(pattern = threshold,list.files("Genetics/"),value = T),value = T)
petal<-petal[grep(pattern = "Hits",petal)]
sepal<-grep(pattern = "Sepal",grep(pattern = threshold,list.files("Genetics/"),value = T),value = T)
sepal<-sepal[grep(pattern = "Hits",sepal)]
ovule<-grep(pattern = "Ovule",grep(pattern = threshold,list.files("Genetics/"),value = T),value = T)
ovule<-ovule[grep(pattern = "Hits",ovule)]
stamen<-grep(pattern = "Stamen",grep(pattern = threshold,list.files("Genetics/"),value = T),value = T)
stamen<-stamen[grep(pattern = "Hits",stamen)]
leaf<-grep(pattern = "Leaf",grep(pattern = threshold,list.files("Genetics/"),value = T),value = T)
leaf<-leaf[grep(pattern = "Hits",leaf)]

petal_list<-unique(rbind(read.table(paste0("Genetics/",petal[1])),read.table(paste0("Genetics/",petal[2])),read.table(paste0("Genetics/",petal[3]))))[,1]
sepal_list<-unique(rbind(read.table(paste0("Genetics/",sepal[1])),read.table(paste0("Genetics/",sepal[2])),read.table(paste0("Genetics/",sepal[3]))))[,1]
ovule_list<-read.table(paste0("Genetics/",ovule[1]))[,1]
stamen_list<-unique(rbind(read.table(paste0("Genetics/",stamen[1])),read.table(paste0("Genetics/",stamen[2])) ))[,1]
leaf_list<-unique(rbind(read.table(paste0("Genetics/",leaf[1])),read.table(paste0("Genetics/",leaf[2])),read.table(paste0("Genetics/",leaf[3]))))[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = paste0("../large_files/Ath_Petal_size/figures/Venn_diagramm_",threshold,"_SNPs.png"),
  output=TRUE,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 100,
  compression = "lzw",
  lwd = 3,
  col=c("green4","greenyellow","lightblue","grey50"),
  fill = c(alpha("green4",0.3), alpha("greenyellow",0.3), alpha("lightblue",0.3), alpha("grey50",0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(0, 0,0,0),
  #cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("green4","greenyellow","lightblue","grey50"),
)


# Venn diagram - GENES

threshold<-c("flct","pvl4")[1]

petal<-grep(pattern = "Petal",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
petal<-petal[grep(pattern = "GO",petal)]

sepal<-grep(pattern = "Sepal",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
sepal<-sepal[grep(pattern = "GO",sepal)]

ovule<-grep(pattern = "Ovule",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
ovule<-ovule[grep(pattern = "GO",ovule)]

stamen<-grep(pattern = "Stamen",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
stamen<-stamen[grep(pattern = "GO",stamen)]

leaf<-grep(pattern = "Leaf",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
leaf<-leaf[grep(pattern = "GO",leaf)]

ft<-grep(pattern = "flower",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
ft<-ft[grep(pattern = "GO",ft)]

petal_list<-unique(rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",petal[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",petal[2]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",petal[3]),h=T))[,5])
sepal_list<-unique(rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",sepal[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",sepal[2]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",sepal[3]),h=T))[,5])
ovule_list<-unique(read.table(paste0("../large_files/Ath_Petal_size/tables/",ovule[1]),h=T))[,5]
stamen_list<-unique(rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",stamen[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",stamen[2]),h=T) )[,5])
leaf_list<-unique(rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",leaf[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",leaf[2]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",leaf[3]),h=T))[,5])
ft_list<-unique(read.table(paste0("../large_files/Ath_Petal_size/tables/",ft[1]),h=T))[,5]

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = paste0("../large_files/Ath_Petal_size/figures/Venn_diagramm_",threshold,"_GENES.png"),
  output=TRUE,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 100,
  compression = "lzw",
  lwd = 3,
  col=c("green4","greenyellow","lightblue","grey50"),
  fill = c(alpha("green4",0.3), alpha("greenyellow",0.3), alpha("lightblue",0.3), alpha("grey50",0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(0, 0,0,0),
  #cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("green4","greenyellow","lightblue","grey50"),
)

intersect(leaf_list,sepal_list)

# Venn diagram - a priori candidates

threshold<-c("flct","pvl4")[1]

petal<-grep(pattern = "Petal",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
petal<-petal[grep(pattern = "GO",petal)]

sepal<-grep(pattern = "Sepal",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
sepal<-sepal[grep(pattern = "GO",sepal)]

ovule<-grep(pattern = "Ovule",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
ovule<-ovule[grep(pattern = "GO",ovule)]

stamen<-grep(pattern = "Stamen",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
stamen<-stamen[grep(pattern = "GO",stamen)]

leaf<-grep(pattern = "Leaf",grep(pattern = threshold,list.files("../large_files/Ath_Petal_size/tables/"),value = T),value = T)
leaf<-leaf[grep(pattern = "GO",leaf)]


petal_list<-rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",petal[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",petal[2]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",petal[3]),h=T))
petal_list<-unique(petal_list[which(petal_list$growth_dev==1),"geneID"])

sepal_list<-rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",sepal[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",sepal[2]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",sepal[3]),h=T))
sepal_list<-unique(sepal_list[which(sepal_list$growth_dev==1),"geneID"])

ovule_list<-read.table(paste0("../large_files/Ath_Petal_size/tables/",ovule[1]),h=T)
ovule_list<-unique(ovule_list[which(ovule_list$growth_dev==1),"geneID"])

stamen_list<-rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",stamen[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",stamen[2]),h=T))
stamen_list<-unique(stamen_list[which(stamen_list$growth_dev==1),"geneID"])

leaf_list<-rbind(read.table(paste0("../large_files/Ath_Petal_size/tables/",leaf[1]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",leaf[2]),h=T),read.table(paste0("../large_files/Ath_Petal_size/tables/",leaf[3]),h=T) )
leaf_list<-unique(leaf_list[which(leaf_list$growth_dev==1),"geneID"])

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = paste0("../large_files/Ath_Petal_size/figures/Venn_diagramm_",threshold,"_growth_dev_GENES.png"),
  output=TRUE,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 100,
  compression = "lzw",
  lwd = 3,
  col=c("green4","greenyellow","lightblue","grey50"),
  fill = c(alpha("green4",0.3), alpha("greenyellow",0.3), alpha("lightblue",0.3), alpha("grey50",0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(0, 0,0,0),
  #cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("green4","greenyellow","lightblue","grey50"),
)

intersect(leaf_list,sepal_list) #"AT1G15340" "AT1G15350"
intersect(petal_list,stamen_list) #"AT1G68610" "AT1G68640" "AT1G68620"


