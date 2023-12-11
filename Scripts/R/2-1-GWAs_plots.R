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


# 2 - Plot Hits with all traits GWAs
hitsdir<-grep(pattern = "Hits",x = list.files("Genetics/",full.names = T),value = T)
for (i in 1:length(hitslist)) {if(i==1){hits<-read.table(hitsdir[i],sep = "_")}else{hits<-rbind(hits,read.table(hitsdir[i],sep = "_"))}}



#first hit
i<-1
(hit<-p3ph[(1:3)+(3*i)-3,])
(gffh<-gff[which(gff$gene %in% hit$V1),])
#(gffh<-gffh[-which(duplicated(gffh$gene)),])
traits<-levels(as.factor(gene$trait))
colors<-c("black","green4","green4","green4","lightblue","lightpink","white","white","white","greenyellow","greenyellow","greenyellow","lightblue")

for (j in 1:length(traits)) {
  
  Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$chr<-paste0("Chr",Assoc$chr)
  Assoc<-Assoc[which(Assoc$chr==unique(gffh$V1)),]
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  ymax<-(-log10(0.05/length(Assoc$Manhattan))+1)
  
  if (j==1) {
  plot(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j],ylab="-log10(p-value)",xlab="",las=1,
       main=paste0("Manhattan plot at locus ",i),ylim=c(-.2,ymax),
       xlim=c(min(gffh[2,4:5])-20000,max(gffh[2,4:5])+20000))
  }else{points(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j])}
}

abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="darkred")
abline(h=-.1,lty=1,col="black")

polygon(x=c(gffh[1,4],gffh[1,4],gffh[1,5],gffh[1,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[1,4],-.1,gffh[1,"gene"],pos=1)
polygon(x=c(gffh[2,4],gffh[2,4],gffh[2,5],gffh[2,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[2,4],-.1,gffh[2,"gene"],pos=1)
polygon(x=c(gffh[3,4],gffh[3,4],gffh[3,5],gffh[3,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[3,4]+2000,-.1,gffh[3,"gene"],pos=1)


#second locus
i<-2
(hit<-p3ph[(1:3)+(3*i)-3,])
(gffh<-gff[which(gff$gene %in% hit$V1),])
#(gffh<-gffh[-which(duplicated(gffh$gene)),])
traits<-levels(as.factor(gene$trait))
colors<-c("black","green4","green4","green4","lightblue","lightpink","white","white","white","greenyellow","greenyellow","greenyellow","lightblue")

for (j in 1:length(traits)) {
  
  Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$chr<-paste0("Chr",Assoc$chr)
  Assoc<-Assoc[which(Assoc$chr==unique(gffh$V1)),]
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  ymax<-(-log10(0.05/length(Assoc$Manhattan))+1)
  
  if (j==1) {
    plot(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j],ylab="-log10(p-value)",xlab="",las=1,
         main=paste0("Manhattan plot at locus ",i),ylim=c(-.2,ymax),
         xlim=c(min(gffh[2,4:5])-20000,max(gffh[2,4:5])+30000))
  }else{points(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j])}
}

abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="darkred")
abline(h=-.1,lty=1,col="black")


polygon(x=c(gffh[2,4],gffh[2,4],gffh[2,5],gffh[2,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[2,4],-.1,gffh[2,"gene"],pos=1)
polygon(x=c(gffh[3,4],gffh[3,4],gffh[3,5],gffh[3,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[3,4]+2000,-.1,gffh[3,"gene"],pos=1)
polygon(x=c(15632712,15632712,15635494,15635494),y=c(0,-.2,-.2,0),col = 1)
text(15632712+3000,-.1,"AT4G32380",pos=1)


# third locus
i<-3
(hit<-p3ph[(1:3)+(3*i)-3,])
(gffh<-gff[which(gff$gene %in% hit$V1),])
#(gffh<-gffh[-which(duplicated(gffh$gene)),])
traits<-levels(as.factor(gene$trait))
colors<-c("black","green4","green4","green4","lightblue","lightpink","white","white","white","greenyellow","greenyellow","greenyellow","lightblue")

for (j in 1:length(traits)) {
  
  Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$chr<-paste0("Chr",Assoc$chr)
  Assoc<-Assoc[which(Assoc$chr==unique(gffh$V1)),]
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  ymax<-(-log10(0.05/length(Assoc$Manhattan))+1)
  
  if (j==1) {
    plot(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j],ylab="-log10(p-value)",xlab="",las=1,
         main=paste0("Manhattan plot at locus ",i),ylim=c(-.2,ymax),
         xlim=c(min(gffh[2,4:5])-20000,max(gffh[2,4:5])+20000))
  }else{points(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j])}
}

abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="darkred")
abline(h=-.1,lty=1,col="black")

polygon(x=c(gffh[1,4],gffh[1,4],gffh[1,5],gffh[1,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[1,4]-1000,-.1,gffh[1,"gene"],pos=1)
polygon(x=c(gffh[2,4],gffh[2,4],gffh[2,5],gffh[2,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[2,4],-.1,gffh[2,"gene"],pos=1)
polygon(x=c(gffh[3,4],gffh[3,4],gffh[3,5],gffh[3,5]),y=c(0,-.2,-.2,0),col = 1)
text(gffh[3,4]+2000,-.1,gffh[3,"gene"],pos=1)

# FT hit
# Chr5:18589247..18591247

for (j in 1:length(traits)) {
  
  Assoc <- read.table(paste0("GWAs/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  Assoc$chr<-paste0("Chr",Assoc$chr)
  Assoc<-Assoc[which(Assoc$chr=="Chr5"),]
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  ymax<-(-log10(0.05/length(Assoc$Manhattan))+1)
  
  if (j==1) {
    plot(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j],ylab="-log10(p-value)",xlab="",las=1,
         main=paste0("Manhattan plot at DOG1"),ylim=c(-.2,ymax),
         xlim=c(18589247-20000,18591247+20000))
  }else{points(Assoc$ps,Assoc$Manhattan,cex=1.5,pch=21,bg=colors[j])}
}

abline(h=-log10(0.05/length(Assoc$Manhattan)),lty=2,col="darkred")
abline(h=-.1,lty=1,col="black")

polygon(x=c(18589247,18589247,18591247,18591247),y=c(0,-.2,-.2,0),col = 1)
text(18589247-1000,-.1,"AT5G45830 (DOG1)",pos=1)

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
# manhattan for each trait separetely
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
