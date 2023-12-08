#------------------
# Petal Size Ath
# Manuscript figures
# Figure 2 - Genome wide association studies
# 2023-12-06
#------------------

# 0- Prepare files
phenotypes<-read.table("phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
names(phenotypes)[5]<-"Ovule_Number"
fam<-read.table("SNP_1001g_filtered.fam")
traits<-colnames(phenotypes)[c(5:16,22)]
for (j in 5:16) {
  fam<-read.table("SNP_1001g_filtered.fam")
  for (i in 1:length(fam$V1)) { if(fam$V1[i] %in% phenotypes$Genotype){fam$V6[i]<-phenotypes[,j][which(phenotypes$Genotype==fam$V1[i])]}}
  fam$V6[which(fam$V6==-9)]<-NA
  write.table(fam,file = paste0("SNP_1001g_filtered_",colnames(phenotypes)[j],".fam"),sep = "\t",quote = F,row.names = F,col.names = F)
}