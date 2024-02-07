#------------------
# Petal Size Ath
# Manuscript figures
# Figure 2 - Genome wide association studies
# 2023-12-06
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

# multivariate model
# data
#this study
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
#merge
btw<-merge(phenotypes,dataset,by.x="Genotype",by.y="accession_name",all.x=T)
# Keep only highly pairwise correlated ones (without FT related traits)
btw<-btw[,c("Genotype","Ovule_Number", "Long_Stamens" , "Petal_Area", "Sepal_Area", "Leaf_Area", "leaf.area.per.leaf.dry.mass", "leaf.nitrogen.content.per.leaf.dry.mass", "whole.plant.dry.mass")]
#make fam file
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


# Compute the genetic correlations
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

