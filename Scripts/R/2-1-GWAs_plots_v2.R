#------------------
# Petal Size Ath
# Manuscript figures
# Figure 2 - Genome wide association studies
# 2023-12-06
#------------------


# 1 - Screen all GWAs for Hits

phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
traits<-colnames(phenotypes)[c(5:16,22)]

for (i in 1:length(traits)) {
  Assoc <- read.table(paste0("Genetics/SNP_1001g_filtered_",traits[i],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  SNP<-which(Assoc$Manhattan>-log10(0.05/length(Assoc$Manhattan)))
  #Hits if any
  snp_list<-Assoc$rs[SNP]
  if (length(snp_list)>0) {
    write.table(snp_list,file = paste0("Genetics/Hits_",traits[i],".txt"),quote = F,row.names = F,col.names = F)}
  # top 100
  top100<-paste0("snp_",
                 Assoc$chr[order(Assoc$Manhattan,decreasing = T)][1:100],"_",
                 Assoc$ps[order(Assoc$Manhattan,decreasing = T)][1:100])
  write.table(top100,file = paste0("Genetics/top100_",traits[i],".txt"),quote = F,row.names = F,col.names = F)
}

# Above data are processed under linux to get gff files around the hit (among others)

# 2 - Plot Hits with all traits GWAs

# 2.1 prepare files

hitsdir<-grep(pattern="txt",grep(pattern = "Hits",x = list.files("Genetics/",full.names = T),value = T),value = T)
for (i in 1:length(hitsdir)) {if(i==1){hits<-read.table(hitsdir[i],sep = "_")}else{hits<-rbind(hits,read.table(hitsdir[i],sep = "_"))}}
hitsdir<-grep(pattern="gff",grep(pattern = "Hits",x = list.files("Genetics/",full.names = T),value = T),value = T)
hitsnames<-grep(pattern="gff",grep(pattern = "Hits",x = list.files("Genetics/",full.names = F),value = T),value = T)
hitsnames<-unlist(lapply(strsplit(hitsnames,split = "\\."),"[[",1))
hitsnames<-paste0(unlist(lapply(strsplit(hitsnames,split = "_"),"[[",2)),"_",unlist(lapply(strsplit(hitsnames,split = "_"),"[[",3)))

# 2.2 make list of gff files for each hit,
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
  #gwas
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

# 2.3 Manhattan plots
load("Genetics/Manhattan_data.Rdata")
## number of analyzed SNPs/var = 540092
bonferonni<-(-log10(.05/540092))
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
abline(h=-1,lty=1,col="black")
  gff<-gffs[[i]]
for (k in 1:dim(gff)[1]) {
polygon(x=c(gff[k,4],gff[k,4],gff[k,5],gff[k,5]),y=c(-.5,-1.5,-1.5,-.5),col = rgb(0,.5,.5,1))
text(gff[k,10],gff[k,11],gff[k,9],cex = .75)
}
}


####################### HERE
# Venn diagram
petal<-grep(pattern = "Petal",grep(pattern = "gene_list",list.files("GWAs/"),value = T),value = T)
sepal<-grep(pattern = "Sepal",grep(pattern = "gene_list",list.files("GWAs/"),value = T),value = T)
ovule<-grep(pattern = "Ovule",grep(pattern = "gene_list",list.files("GWAs/"),value = T),value = T)
stamen<-grep(pattern = "Sta",grep(pattern = "gene_list",list.files("GWAs/"),value = T),value = T)
leaf<-grep(pattern = "Leaf",grep(pattern = "gene_list",list.files("GWAs/"),value = T),value = T)

petal_list<-unique(rbind(read.table(paste0("GWAs/",petal[1])),read.table(paste0("GWAs/",petal[2])),read.table(paste0("GWAs/",petal[3]))))[,1]
sepal_list<-unique(rbind(read.table(paste0("GWAs/",sepal[1])),read.table(paste0("GWAs/",sepal[2])),read.table(paste0("GWAs/",sepal[3]))))[,1]
ovule_list<-read.table(paste0("GWAs/",ovule[1]))[,1]
stamen_list<-unique(rbind(read.table(paste0("GWAs/",stamen[1])),read.table(paste0("GWAs/",stamen[2])) ))[,1]
leaf_list<-unique(rbind(read.table(paste0("GWAs/",leaf[1])),read.table(paste0("GWAs/",leaf[2])),read.table(paste0("GWAs/",leaf[3]))))[,1]

library(VennDiagram)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

venn.diagram(
  x = list(petal_list,sepal_list,ovule_list,stamen_list,leaf_list),
  category.names = c("Leaf","Ovule","Petal","Sepal","Stamen"),
  filename = "Figures/Venn_diagramm_top100.png",
  output=TRUE,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 100,
  compression = "lzw",
  lwd = 3,
  col=c("green4","lightpink","grey50","greenyellow","lightblue"),
  fill = c(alpha("green4",0.3), alpha("lightpink",0.3), alpha("grey50",0.3), alpha("greenyellow",0.3), alpha("lightblue",0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-50, 10, 120,90,20),
  #cat.dist = c(0.055, 0.055, 0.085,0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("green4","lightpink","grey50","greenyellow","lightblue"),
)

# Plot SNPs effects
library(vcfR)

# snp_1_28960616
vcf<-read.vcfR("GWAs/annotated_top100_Petal_Width.ann.vcf")
snp_1_28960616<-data.frame(genotypeID=colnames(vcf@gt),snp_1_28960616=vcf@gt[which(vcf@fix[,3] == "snp_1_28960616"),])[-1,]
snp_1_28960616$genotypeID<-unlist(strsplit(snp_1_28960616$genotypeID,"_"))[(1:length(snp_1_28960616$genotypeID))*2]
#  snp_4_15627234
vcf<-read.vcfR("GWAs/annotated_top100_Sepal_Area.ann.vcf")
snp_4_15627234<-data.frame(genotypeID=colnames(vcf@gt), snp_4_15627234=vcf@gt[which(vcf@fix[,3] == "snp_4_15627234"),])[-1,]
snp_4_15627234$genotypeID<-unlist(strsplit( snp_4_15627234$genotypeID,"_"))[(1:length( snp_4_15627234$genotypeID))*2]
# snp_4_1601862
vcf<-read.vcfR("GWAs/annotated_top100_Sepal_Length.ann.vcf")
snp_4_1601862<-data.frame(genotypeID=colnames(vcf@gt),snp_4_1601862=vcf@gt[which(vcf@fix[,3] == "snp_4_1601862"),])[-1,]
snp_4_1601862$genotypeID<-unlist(strsplit(snp_4_1601862$genotypeID,"_"))[(1:length(snp_4_1601862$genotypeID))*2]

genotypes<-merge(snp_1_28960616,snp_4_15627234)
genotypes<-merge(genotypes,snp_4_1601862)

phenotypes<-read.table("U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
phenotypes<-merge(phenotypes,genotypes,by.x="Genotype",by.y="genotypeID")
boxplot(phenotypes$Petal_Width~phenotypes$snp_1_28960616)
boxplot(phenotypes$Sepal_Area~phenotypes$snp_4_15627234)
boxplot(phenotypes$Sepal_Length~phenotypes$snp_4_1601862)





# Not used
# Manhattan for each trait separately
i<-1
(geneid<-gff$gene[i])
(traits<-unique(gene[which(gene$gene==gff$gene[i]),"trait"]))
(info<-gff[which(gff$gene==geneid),][1,])

par(mfrow=c(length(traits),1))
for (j in 1:length(traits)) {
  
  Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$chr<-paste0("Chr",Assoc$chr)
  Assoc<-Assoc[which(Assoc$chr==info$V1[1]),]
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  ymax<-(-log10(0.05/length(Assoc$Manhattan))+1)
  
  plot(Assoc$ps,Assoc$Manhattan,cex=0.5,pch=16,col=1,ylab="-log10(p-value)",xlab="",las=1,
       main=paste0("GWAs ",traits[j]),ylim=c(0,ymax),
       xlim=c(min(info[,4:5])-20000,max(info[,4:5])+20000))
  abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="red")
  
  polygon(x=c(info[,4],info[,4],info[,5],info[,5]),y=c(0,.2,.2,0),col = 1)
}

dev.off()

# 
