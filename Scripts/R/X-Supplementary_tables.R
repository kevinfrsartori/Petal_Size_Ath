# Combine tables
#install.packages("xlsx")
library(xlsx)

#annotated SNPs
traits<-c(colnames(read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1))[c(22,8,9,10,11,12,13,14,15,16,6,7)],"Ovule_Number")
for (i in 1:length(traits)) {
ann<-read.table(paste0("../large_files/Ath_Petal_size/tables/annotated_simple_flct_Hits_",traits[i],".iHS.AD.GO.txt"),h=T)
ann<-ann[order(ann$iHS_pval,decreasing = T),]
if(length(which(duplicated(ann)))>0){ ann<-ann[-which(duplicated(ann)),] }
write.xlsx(ann, file="../large_files/Ath_Petal_size/tables/Table_S4_SNPs_annotated.xlsx", sheetName=traits[i], row.names=FALSE,append = T)
}
# gene ontologies
for (i in 1:length(traits)) {
  go<-read.table(paste0("../large_files/Ath_Petal_size/tables/gene_ontology_flct_Hits_",traits[i],".txt"),sep = "\t")
  write.xlsx(go, file="../large_files/Ath_Petal_size/tables/Table_S5_Gene_ontologies.xlsx", sheetName=traits[i], row.names=FALSE,append = T)
}


# Accession list and traits
phenotypes<-read.table("Phenotypes/U_Shaped_Data_corrected_2023-05-05.csv",h=T,as.is = 1)
colnames(phenotypes)[5]<-"Ovule_Number"
colnames(phenotypes)[1]<-"Accession_ID"
phenotypes<-phenotypes[-which(is.na(phenotypes$Petal_Area)),]
write.xlsx(phenotypes[,c(1:16,22)], file="../large_files/Ath_Petal_size/tables/Table_S1_Accession_list_phenotypes.xlsx", sheetName="Accession_list", row.names=FALSE)
