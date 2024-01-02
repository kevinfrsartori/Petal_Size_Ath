#------------------
# Petal Size Ath
# Manuscript figures
# Figure 4 - Climate relationship and fitness tradeoff
# 2024-01-02
#------------------

# Panel A - Habitat suitability map

# Panel B - Trait-HS relationship

# Panel C - PiN/PiS-HS relationship

# Panel D - SFS-HS
handtable<-as.matrix(read.table("Genetics/sfs/DAF_table_rel.95_.1to.5.txt"))
colnames(handtable)<-substr(colnames(handtable),2,5)
shade<-gray.colors(10,start = .8,end = .2,gamma = 1)
barplot(height = handtable[,1:9],beside = T,col=shade,las=1,main = "Allele counts per frequency")
legend(x = 75,y = 4e+05,legend = rownames(handtable),fill = shade,cex = .75,ncol = 2)

# Panel E - Large petal allele-HS
