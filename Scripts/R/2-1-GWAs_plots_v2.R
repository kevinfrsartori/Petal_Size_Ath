#------------------
# Petal Size Ath
# Manuscript figures
# Figure 2 - Genome wide association studies
# 2023-12-06
#------------------


# 1 - Screen all GWAs for Hits
#-----------------------------

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
#-----------------------------------

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


# 3 - VENN DIAGRAMS
#------------------
library(VennDiagram)
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# All traits
# Venn diagram - SNPs
petal<-grep(pattern = "Petal",grep(pattern = "top100",list.files("Genetics/"),value = T),value = T)
petal<-petal[-grep(pattern = "annotated",petal)]
sepal<-grep(pattern = "Sepal",grep(pattern = "top100",list.files("Genetics/"),value = T),value = T)
sepal<-sepal[-grep(pattern = "annotated",sepal)]
ovule<-grep(pattern = "Ovule",grep(pattern = "top100",list.files("Genetics/"),value = T),value = T)
ovule<-ovule[-grep(pattern = "annotated",ovule)]
stamen<-grep(pattern = "Stamen",grep(pattern = "top100",list.files("Genetics/"),value = T),value = T)
stamen<-stamen[-grep(pattern = "annotated",stamen)]
leaf<-grep(pattern = "Leaf",grep(pattern = "top100",list.files("Genetics/"),value = T),value = T)
leaf<-leaf[-grep(pattern = "annotated",leaf)]

petal_list<-unique(rbind(read.table(paste0("Genetics/",petal[1])),read.table(paste0("Genetics/",petal[2])),read.table(paste0("Genetics/",petal[3]))))[,1]
sepal_list<-unique(rbind(read.table(paste0("Genetics/",sepal[1])),read.table(paste0("Genetics/",sepal[2])),read.table(paste0("Genetics/",sepal[3]))))[,1]
ovule_list<-read.table(paste0("Genetics/",ovule[1]))[,1]
stamen_list<-unique(rbind(read.table(paste0("Genetics/",stamen[1])),read.table(paste0("Genetics/",stamen[2])) ))[,1]
leaf_list<-unique(rbind(read.table(paste0("Genetics/",leaf[1])),read.table(paste0("Genetics/",leaf[2])),read.table(paste0("Genetics/",leaf[3]))))[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = "Figures/Venn_diagramm_All_top_100_SNPs.png",
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



# Venn diagram - Genes
petal<-grep(pattern = "Petal",grep(pattern = "annotated",list.files("Genetics/"),value = T),value = T)
sepal<-grep(pattern = "Sepal",grep(pattern = "annotated",list.files("Genetics/"),value = T),value = T)
ovule<-grep(pattern = "Ovule",grep(pattern = "annotated",list.files("Genetics/"),value = T),value = T)
stamen<-grep(pattern = "Stamen",grep(pattern = "annotated",list.files("Genetics/"),value = T),value = T)
leaf<-grep(pattern = "Leaf",grep(pattern = "annotated",list.files("Genetics/"),value = T),value = T)

petal_list<-unique(rbind(read.table(paste0("Genetics/",petal[1])),read.table(paste0("Genetics/",petal[2])),read.table(paste0("Genetics/",petal[3]))))[,1]
sepal_list<-unique(rbind(read.table(paste0("Genetics/",sepal[1])),read.table(paste0("Genetics/",sepal[2])),read.table(paste0("Genetics/",sepal[3]))))[,1]
ovule_list<-read.table(paste0("Genetics/",ovule[1]))[,1]
stamen_list<-unique(rbind(read.table(paste0("Genetics/",stamen[1])),read.table(paste0("Genetics/",stamen[2])) ))[,1]
leaf_list<-unique(rbind(read.table(paste0("Genetics/",leaf[1])),read.table(paste0("Genetics/",leaf[2])),read.table(paste0("Genetics/",leaf[3]))))[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = "Figures/Venn_diagramm_All_genes.png",
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

# Area
# Venn diagram - SNPs

petal_list<-read.table("Genetics/top100_Petal_Area.txt")[,1]
sepal_list<-read.table("Genetics/top100_Sepal_Area.txt")[,1]
leaf_list<-read.table("Genetics/top100_Leaf_Area.txt")[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,petal_list),
  category.names = c("Leaf","Sepal","Petal"),
  filename = "Figures/Venn_diagramm_Area_top_100_SNPs.png",
  output=TRUE,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 100,
  compression = "lzw",
  lwd = 3,
  col=c("green4","greenyellow","grey50"),
  fill = c(alpha("green4",0.3), alpha("greenyellow",0.3),  alpha("grey50",0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(25,35,45),
  #cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("green4","greenyellow","grey50"),
)

# Venn diagram - Genes
petal_list<-read.table("Genetics/annotated_simple_top100_Petal_Area.txt")[,1]
sepal_list<-read.table("Genetics/annotated_simple_top100_Sepal_Area.txt")[,1]
leaf_list<-read.table("Genetics/annotated_simple_top100_Leaf_Area.txt")[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,petal_list),
  category.names = c("Leaf","Sepal","Petal"),
  filename = "Figures/Venn_diagramm_Area_genes.png",
  output=TRUE,
  imagetype="png" ,
  height = 800 , 
  width = 800 , 
  resolution = 100,
  compression = "lzw",
  lwd = 3,
  col=c("green4","greenyellow","grey50"),
  fill = c(alpha("green4",0.3), alpha("greenyellow",0.3),  alpha("grey50",0.3)),
  cex = 2,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(25,35,45),
  #cat.dist = c(0.055, 0.055, 0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("green4","greenyellow","grey50"),
)


# Length
# Venn diagram - SNPs
petal_list<-read.table("Genetics/top100_Petal_Length.txt")[,1]
sepal_list<-read.table("Genetics/top100_Sepal_Length.txt")[,1]
stamen_list<-read.table("Genetics/top100_Long_Stamens.txt")[,1]
leaf_list<-read.table("Genetics/top100_Leaf_Length.txt")[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = "Figures/Venn_diagramm_Length_top_100_SNPs.png",
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



# Venn diagram - Genes
petal_list<-read.table("Genetics/annotated_simple_top100_Petal_Length.txt")[,1]
sepal_list<-read.table("Genetics/annotated_simple_top100_Sepal_Length.txt")[,1]
stamen_list<-read.table("Genetics/annotated_simple_top100_Long_Stamens.txt")[,1]
leaf_list<-read.table("Genetics/annotated_simple_top100_Leaf_Length.txt")[,1]

venn.diagram(
  x = list(leaf_list,sepal_list,stamen_list,petal_list),
  category.names = c("Leaf","Sepal","Stamen","Petal"),
  filename = "Figures/Venn_diagramm_Length_genes.png",
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



# Petal Area Manhattan plot
par(mfrow=c(1,1),mar=c(4,3,2,1),oma=c(0,0,0,0))

Assoc <- read.table("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_Petal_Area.assoc.txt",h=T,sep="\t",dec=".")
bonferonni<-(-log10(0.05/dim(Assoc)[1]))
Assoc$pos<-Assoc$ps
for (i in 2:5) { Assoc$pos[which(Assoc$chr==i)]<-Assoc$pos[which(Assoc$chr==i)]+max(Assoc$pos[which(Assoc$chr==i-1)]) }
#Assoc$chr<-paste0("Chr",Assoc$chr)
Assoc$Manhattan<-(-log10(Assoc$p_lrt))
ymax<-(bonferonni+1)
Assoc$pos<-Assoc$pos/1000000
Assoc<-Assoc[which(Assoc$Manhattan>1),]

palette(c("gray30","gray60","gray40","gray70","gray50"))
plot(Assoc$pos,Assoc$Manhattan,cex=0.5,pch=16,col=Assoc$chr,ylab="",xlab="",las=1,
     main="Petal Area GWAs",ylim=c(1,ymax))
axis(side = 2,at = 4,labels = "-log10(p-value)",line = 1,tick = F)
axis(side = 1,at = 60,labels = "SNPs relative position (Mbp)",line = 1,tick = F)

abline(h=bonferonni,lty=2,col="darkred")

top100<-order(Assoc$Manhattan,decreasing = T)[1:100]
palette(c(rgb(0,.5,.5,1),rgb(0,.8,.8,1),rgb(0,.6,.6,1),rgb(0,.9,.9,1),rgb(0,.7,.7,1)))
points(Assoc$pos[top100],Assoc$Manhattan[top100],cex=0.75,pch=16,col=Assoc$chr[top100])

# Plot top 100 for fig 2C
par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(Assoc$pos[top100],Assoc$Manhattan[top100],
     cex=2,pch=16,col=Assoc$chr[top100],xaxt="n", yaxt="n",bty="n",ylab="",xlab="")


# All traits Manhattan plots

par(mfrow=c(6,1),mar=c(1,3,0,.2),oma=c(4,0,0,0))
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
traits<-colnames(phenotypes)[c(5:16,22)]
traits<-traits[-which(traits=="Petal_Area")]
for (j in 1:length(traits)) {
  
  Assoc <- read.table(paste0("../large_files/Ath_Petal_size/gwas/SNP_1001g_filtered_",traits[j],".assoc.txt"),h=T,sep="\t",dec=".")
  bonferonni<-(-log10(0.05/dim(Assoc)[1]))
  Assoc$pos<-Assoc$ps
  for (i in 2:5) { Assoc$pos[which(Assoc$chr==i)]<-Assoc$pos[which(Assoc$chr==i)]+max(Assoc$pos[which(Assoc$chr==i-1)]) }
  Assoc$Manhattan<-(-log10(Assoc$p_lrt))
  ymax<-(bonferonni+1)
  Assoc$pos<-Assoc$pos/1000000
  Assoc<-Assoc[which(Assoc$Manhattan>1),]
  
  palette(c("gray30","gray60","gray40","gray70","gray50"))
  plot(Assoc$pos,Assoc$Manhattan,cex=0.5,pch=16,col=Assoc$chr,
       ylab="",xlab="",las=1,xaxt="n",
       ylim=c(1,ymax))
  axis(side = 2,at = 4,labels = "-log10(p-value)",line = 1,tick = F)
  axis(1,at = c(0,20,40,60,80,100,120),labels = c("","","","","","",""))
  
  text(0,7.5,gsub(pattern = "_",replacement = " ",traits[j]),pos=4)
  abline(h=bonferonni,lty=2,col="darkred")
  if (j %in% c(6,12)) {   
    axis(1,at = c(0,20,40,60,80,100,120),labels = c(0,20,40,60,80,100,120))
    axis(side = 1,at = 60,labels = "SNPs relative position (Mbp)",line = 1,tick = F)
}
  
}


# 4 - Filter genes by function
petal_list<-read.table("Genetics/annotated_simple_top100_Petal_Area.txt")[,1:3]
colnames(petal_list)<-c("gene_ID","gene_name","SNP")
petal_list$link<-paste0("https://www.arabidopsis.org/servlets/TairObject?type=locus&name=",petal_list$gene_ID)
write.csv2(petal_list,file = "Genetics/functionnal_annotation_Petal_Area.csv",quote = F,row.names = F)

