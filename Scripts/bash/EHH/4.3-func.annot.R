#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
#args=c("flct_Hits_Petal_Area","$PATH-TO-SERVER/A_thaliana/GWAs/Petal_Size_Ath/annotations/")

##########
#
# Make table with info about genes PO and GO
# KS - 2023-08-24
#
##########
gene<-read.table(paste0(args[2],"annotated_simple_",args[1],".iHS.AD.txt"),
                 h=T,sep="\t")

# Gene ontology
library(topGO)
library("biomaRt")
#collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'plants.ensembl.org')
# Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))

all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
int.genes <- unique(gene$geneID)
int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP'
              , allGenes = int.genes
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO
              ,nodeSize = 10
)
results <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
results.tab <- GenTable(object = go.obj, elimFisher = results,topNodes=10,numChar=100)
 # Keep gene ontology table for supplement
write.table(x = results.tab, file = paste0(args[2],"gene_ontology_",args[1],".txt"),
            quote=F,row.names = F,col.names = T,sep="\t")

# Use GOslim term to filter gene lists
GO<-read.table("$PATH-TO-SERVER/A_thaliana/REF/ATH_GO_GOSLIM.txt",skip = 4,sep="\t",fill = T,quote = "")
#keep only biological process
GO<-GO[which(GO$V8=="P"),]
#make list for each categories
growth_dev_GO<-levels(as.factor(GO$V9))[c(11,17,29,33,2,10,20,22,8)]
flower_dev_GO<-levels(as.factor(GO$V9))[c(17)]

growth_dev_genes<-unique(GO$V1[which(GO$V9 %in% growth_dev_GO)])
flower_dev_genes<-unique(GO$V1[which(GO$V9 %in% flower_dev_GO)])

# growth dev. genes:
print("growth dev. genes:")
unique(gene$geneID[which(gene$geneID %in% growth_dev_genes)])
# flower dev. genes:
print("flower dev. genes:")
unique(gene$geneID[which(gene$geneID %in% flower_dev_genes)])

#make filtering columns
gene$flower_dev<-0
gene$flower_dev[which(gene$geneID %in% flower_dev_genes)]<-1
gene$growth_dev<-0
gene$growth_dev[which(gene$geneID %in% growth_dev_genes)]<-1

#write table
write.table(gene,paste0(args[2],"annotated_simple_",args[1],".iHS.AD.GO.txt"),sep="\t",quote = F,row.names = F,col.names = T)
