# 
#  
# Make list of orthologs
# Maximizing the number of genes for Athaliana
# 

# AT_N0.tsv was made with grep AT from N0.tsv

ortholist<-read.table(file = "/crex/proj/snic2020-16-182/proteomes/A.tha447.hal.lyr.013024/OrthoFinder/Results_Jan30/Phylogenetic_Hierarchical_Orthogroups/AT_N0.txt",h=F,sep = "\t",na.strings = "")
colnames(ortholist)<-c("HOG","OG","N","Ahal","Alyr","Atha")

# Keep only one gene Ahal and Alyr
ortholist$Ahal<-unlist(lapply(strsplit(ortholist$Ahal,split = ", "),"[[",1))
ortholist$Alyr<-unlist(lapply(strsplit(ortholist$Alyr,split = ", "),"[[",1))

# Remove rows if no ortholog at all
naAhal<-which(is.na(ortholist$Ahal))
naAlyr<-which(is.na(ortholist$Alyr))
ortholist<-ortholist[-intersect(naAhal,naAlyr),]

# Fill the NA with gene from the same OG if possible
ortholist$OG<-as.factor(ortholist$OG)
for (i in 1:length(levels(ortholist$OG))) {
# for (i in 1:20) {
  ortholist.t<-ortholist[which(ortholist$OG==levels(ortholist$OG)[i]),]
if(dim(ortholist.t)[1]>1){
# Fill Ahal
  if(length(na.omit(ortholist.t$Ahal))>0){
    ortholist.t$Ahal[1]<-na.omit(ortholist.t$Ahal)[1]

    for (j in 2:dim(ortholist.t)[1]) {
      if(is.na(ortholist.t$Ahal[j])) { ortholist.t$Ahal[j]<-ortholist.t$Ahal[j-1] }
    }
  }
# Fill Alyr
  if(length(na.omit(ortholist.t$Alyr))>0){
    ortholist.t$Alyr[1]<-na.omit(ortholist.t$Alyr)[1]
    for (j in 2:dim(ortholist.t)[1]) {
      if(is.na(ortholist.t$Alyr[j])) { ortholist.t$Alyr[j]<-ortholist.t$Alyr[j-1] }
    }
  }
}
# Make new dataset
if (i==1) { ortholist.cor<-ortholist.t }else{ ortholist.cor<-rbind(ortholist.cor,ortholist.t) }
  
}


# treat apart one ATxG versus several ATxG
ortholist_multi<-ortholist.cor[grep(",",ortholist.cor$Atha),]
ortholist_single<-ortholist.cor[-grep(",",ortholist.cor$Atha),]

# Split Ath multi gene OG
for (i in 1:dim(ortholist_multi)[1]) {
  ortholist_multi.t<-ortholist_multi[i,]
  Atha<-unlist(strsplit(ortholist_multi.t$Atha,split = ", "))
  # copy the orthologs for all
  ortholist_multi.t2<-data.frame(HOG=ortholist_multi.t$HOG, OG=ortholist_multi.t$OG,
                                 N=ortholist_multi.t$N, Ahal=ortholist_multi.t$Ahal,
                                 Alyr=ortholist_multi.t$Alyr, Atha=Atha)
# Make new dataset
  if (i==1) { ortholist_multi.cor<-ortholist_multi.t2 }else{ ortholist_multi.cor<-rbind(ortholist_multi.cor,ortholist_multi.t2) }
  
}

ortholist.final<-rbind(ortholist_single,ortholist_multi.cor)

dim(na.omit(ortholist.final))

ortholist.final$Atha<-substr(ortholist.final$Atha,1,9)

write.table(x = ortholist.final, file = "one_ortho_per_species.txt",quote = F,col.names = F,row.names = F)
